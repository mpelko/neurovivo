'''
Created on May 29, 2012

@author: mpelko
'''
import neurovivo.common as cmn
import numpy as np

def recursive_analysis(parameters):
    PATH_RESULTS = parameters["PATH_RESULTS"]
    chosen_parameters=parameters["chosen_sweep_params"]
    if "PATH_RESULTS_REMOTE" in parameters.keys():
        PATH_RESULTS = parameters["PATH_RESULTS_REMOTE"]
    else:
        PATH_RESULTS = PATH_RESULTS + "/results/" 
    results = cmn.load_results(PATH_RESULTS, chosen_parameters)
    #print chosen_parameters
    for i, function_spec in enumerate(parameters["analysis_functions"]):
        if not "kwargs" in function_spec.keys():
            function_spec["kwargs"] = {}
        result=function_spec["func"](results, parameters, **function_spec["kwargs"])
        parameters["analysis_results"][i].append(result)

# the mother of all analysis functions.
def batch_analysis(PATH_RESULTS, sweep_parameter_keys, fixed_parameters, functions, save_analysis=False):
    """
    TODO: document this
    
    If the original simulation sweep was more extensive then the analysis query, simply fix some of the
    originally sweeped parameters by specifying their values in a dictionary input under fixed_parameters.
    
    function will be the function which goes into the sweep recursion and produce the results.
    """

    # checking if all required parameters are supplied
    simulation_sweep_parameters = cmn.load_yaml(PATH_RESULTS+"/sweep_params/sweep_params.yml")
    simulation_sweep_parameters_keys = simulation_sweep_parameters.keys()
    all_incoming_parameters=sweep_parameter_keys[:]
    all_incoming_parameters.extend(fixed_parameters.keys())
    simulation_sweep_parameters_keys.sort()
    all_incoming_parameters.sort()
    #print fixed_parameters,all_incoming_parameters, simulation_sweep_parameters_keys
    assert simulation_sweep_parameters_keys==all_incoming_parameters, "{}, {}".format(simulation_sweep_parameters_keys,all_incoming_parameters)
    
    try:
        metadata = cmn.load_metadata(PATH_RESULTS, method="yaml")
    except:
        metadata = cmn.load_metadata(PATH_RESULTS)
    
    parameters_additional = {
                  "PATH_RESULTS": PATH_RESULTS,
                  "sweep_parameter_keys": sweep_parameter_keys,
                  "chosen_parameters": fixed_parameters,
                  "analysis_results": [[] for _ in functions],
                  "analysis_functions": functions,
                  }

    parameters = dict(metadata.items()+parameters_additional.items())
    
    # fills in all the possible values of sweep_parameters
    sweep_parameters_list = []
    sweep_parameters = {}
    sweep_parameter_keys.sort()
    for key in sweep_parameter_keys:
        sweep_parameters_list.append({key:simulation_sweep_parameters[key]})
        sweep_parameters[key]=simulation_sweep_parameters[key]
        
    PATH = PATH_RESULTS+"/analysis/"
    cmn.mkdir_p(PATH)
    cmn.parameter_sweep(recursive_analysis, parameters, sweep_parameters, fixed_parameters)

    result = None
    for i,function_spec in enumerate(functions):
        result = parameters["analysis_results"][i]
        kwargs = cmn.get_default_function_parameters(function_spec["func"])
        if not "kwargs" in function_spec.keys():
            function_spec["kwargs"]={}
        for kwarg_key in function_spec["kwargs"].keys():
            kwargs[kwarg_key]=function_spec["kwargs"][kwarg_key]
        metadata={
                  "sweep_parameter_keys":sweep_parameter_keys,
                  "fixed_parameters":fixed_parameters,
                  "function_name":function_spec["func"].__name__,
                  "func_kwargs":kwargs,
                  }
        if save_analysis:
            PATH_SAVE = PATH + function_spec["func"].__name__
            PATH_SAVE = cmn.mkdir_safe(PATH_SAVE+"/analys_results")
            cmn.save_metadata(PATH_SAVE, metadata)
            cmn.save_pickle(PATH_SAVE+"/analysis.pkl", {"result":result,"sweep_parameters":sweep_parameters_list})
    return {"result":result,"sweep_parameters":sweep_parameters_list}

#TODO: Write a function that also returns the partial sweeps instead of full sweeps.

def get_analysis_result(PATH_RESULTS, sweep_parameter_keys, fixed_parameters, function, save_analysis=False, **func_kwargs):
    """
    Looks into the analysis folder if the required analysis has been done already.
    If so, it returns the result. If not it returns None.
    """
    
    PATH = existing_analysis_path(PATH_RESULTS, sweep_parameter_keys, fixed_parameters, function, **func_kwargs)
    if PATH:
        return partial_analysis_result(PATH, PATH_RESULTS, sweep_parameter_keys, fixed_parameters)
    else:
        # if there is no available analysis already preformed, start a new one.
        return batch_analysis(PATH_RESULTS, sweep_parameter_keys, fixed_parameters, [{"func":function,"kwargs":func_kwargs}], save_analysis=save_analysis)

def existing_analysis_path(PATH_RESULTS, sweep_parameter_keys, fixed_parameters, function, **func_kwargs):
    """
    Checks if the analysis including the reqested analysis has already been performed.
    If yes it returns the path to where it is stored. If not it returns None.
    """
    
    result = None
    
    cmn.mkdir_p(PATH_RESULTS+"/analysis/")
    PATH = PATH_RESULTS+"/analysis/"+function.__name__
    analysis_folders = cmn.list_files(PATH_RESULTS+"/analysis/", folders=True)
    
    if function.__name__ in analysis_folders:
    
        potential = cmn.list_files(PATH, folders=True)
        kwargs = cmn.get_default_function_parameters(function)
        for kwarg_key in func_kwargs.keys():
            kwargs[kwarg_key]=func_kwargs[kwarg_key]
                      
        for pot in potential:
            md = cmn.load_metadata(PATH+"/"+pot)
                        
            if not function.__name__ == md["function_name"]:
                continue
            
            if not kwargs == md["func_kwargs"]:
                continue
            
            all_sweep_parameters1 = set(sweep_parameter_keys + fixed_parameters.keys())
            all_sweep_parameters2 = set(md["sweep_parameter_keys"] + md["fixed_parameters"].keys())
            if not all_sweep_parameters1 == all_sweep_parameters2:
                continue
                                    
            fp_check = [item in fixed_parameters.items() for item in md["fixed_parameters"].items()]
            if not np.all(fp_check):
                continue
            
            # UGLY and EXPENSIVE BUT WORKS. Perhpaps find a better way to check, without loading the result.
            tmp_res = cmn.load_pickle(PATH+"/"+pot+"/analysis.pkl")["sweep_parameters"]
            tmp_res_dict = {}
            for tmp in tmp_res:
                tmp_res_dict = dict(tmp_res_dict.items() + tmp.items())
            stay_in = True
            for item in fixed_parameters.items():
                if not (item in md["fixed_parameters"].items()):
                    if not item[1] in tmp_res_dict[item[0]]:
                        stay_in = False
                        break
            if not stay_in: 
                continue
            result = PATH+"/"+pot+"/"
            break
    
    return result

def partial_analysis_result(PATH, PATH_RESULTS, sweep_parameter_keys, fixed_parameters):
    full_result = cmn.load_pickle(PATH+"/analysis.pkl")
    params = {"res":full_result["result"], "sp_all":full_result["sweep_parameters"], "partial_result":[]}
    
    md = cmn.load_metadata(PATH)
    if md["fixed_parameters"] == fixed_parameters:
        return full_result
    
    def for_recursion(params):
        full_results = params["res"]
        sizes = []
        indices = []
        for sp in params["sp_all"]:
            key = sp.keys()[0]
            sizes.append(len(sp[key]))
            indices.append(np.where(sp[key]==params["chosen_sweep_params"][key])[0][0])
        try:
            full_results = np.array(full_results)
            full_results = full_results.reshape(sizes)
            params["partial_result"].append(full_results[tuple(indices)])
        except:
            print "*** Partial result retrevial is not available for this analysis. Will need to redo the analysis."
            # TODO: enable the partial retrevial fully.
            raise
    
    simulation_sweep_parameters = cmn.load_pickle(PATH_RESULTS+"/sweep_params/sweep_params.pkl")
    #print simulation_sweep_parameters
    sweep_parameters_list = []
    sweep_parameters = {}
    sweep_parameter_keys.sort()
    for key in sweep_parameter_keys:
        sweep_parameters_list.append({key:simulation_sweep_parameters[key]})
        sweep_parameters[key]=simulation_sweep_parameters[key]
    
    cmn.parameter_sweep(for_recursion, params, sweep_parameters, fixed_parameters)
    return {"result":params["partial_result"],"sweep_parameters":sweep_parameters_list}
    
    