'''
Created on Oct 9, 2012

@author: mpelko
'''

from neurovivo.cell.HayPassiveCell import HayPassiveCell

import matplotlib.pylab as pyl
import neurovivo.plotting as plt
from neurovivo.cell.synapse import exp2syn_GABA_synapse as syn_inh
from neurovivo.cell.synapse import exp2syn_AMPA_synapse as syn_exc

PLOT = False

def test_HayPassiveCell():
    cell = HayPassiveCell()
    
def test_simulate_simple_input():
    cell = HayPassiveCell()
    
    sim_params = {
                  "tstop": 400,
                  "dt": 0.1
                  }
    
    timings = [[200]]
    locations = [(cell.secs["soma"][0].name(),0.5)]
    sim_params["synaptic_input"] = [{"timings": timings, "locations": locations}]
    trace = cell.simulate(sim_params)[0]
    plt.plot_trace(trace)
    
    timings = [[200,300]]
    locations = [(cell.secs["soma"][0].name(),0.5)]
    sim_params["synaptic_input"] = [{"timings": timings, "locations": locations}]
    trace = cell.simulate(sim_params)[0]
    plt.plot_trace(trace)
    
    timings = [[200],[300, 325]]
    locations = [cell.secs["soma"][0].name(), cell.secs["basal"][0].name()]
    sim_params["synaptic_input"] = [{"timings": timings, "locations": locations}]
    trace = cell.simulate(sim_params)[0]

    if PLOT:
        plt.plot_trace(trace)
        pyl.show()

def test_simulate_current_injection():
    cell = HayPassiveCell()
    
    sim_params = {
                  "tstop": 400,
                  "dt": 0.1
                  }
    curr_inject_params1 = {"tstart":200, "duration":100, "magnitude":0.1, "location":cell.secs["soma"][0].name()}
    curr_inject_params2 = {"tstart":150, "duration":200, "magnitude":0.05, "location":cell.secs["soma"][0].name()}
    sim_params["current_injection"] = [curr_inject_params1,curr_inject_params2]     
    trace = cell.simulate(sim_params)[0]
    if PLOT:
        plt.plot_trace(trace)
        pyl.show()
    
def test_synaptic_input():
    cell = HayPassiveCell()
    
    sim_params = {
                  "tstop": 400,
                  "dt": 0.1
                  }
    
    synaptic_input_params={"timings":[[100,150,200]], 
                           "locations":[(cell.secs["soma"][0].name(),0.5)], 
                           "weights": 5*0.0003, 
                           "Syn_type":syn_inh.Exp2syn_GABA
                           }
    synaptic_input_params2={"timings":[[100,120,200]], 
                            "locations":[(cell.secs["soma"][0].name(),0.5)], 
                            "weights": 0.0003, 
                            "Syn_type":syn_exc.Exp2syn_AMPA}
    sim_params["synaptic_input"] = [synaptic_input_params,synaptic_input_params2]     
    trace = cell.simulate(sim_params)[0]
    if PLOT:
        plt.plot_trace(trace)
        pyl.show()

test_HayPassiveCell()        
test_simulate_simple_input()
test_simulate_current_injection()
test_synaptic_input()

