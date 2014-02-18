'''
Created on Mar 29, 2012

@author: mpelko
'''

from neurovivo.cell.HayCell import HayCell
import matplotlib.pylab as pyl
import neurovivo.plotting.plotting as plt
from neurovivo.cell.synapse import exp2syn_GABA_synapse as syn_inh
from neurovivo.cell.synapse import exp2syn_AMPA_synapse as syn_exc

PLOT = False

def test_simulate():
    
    cell = HayCell()
    trace = cell.simulate()[0]
    if PLOT:
        plt.plot_trace(trace)
        pyl.show()
    
def test_simulate_simple_input():
    cell = HayCell()
    
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

    
def test_synaptic_input():
    cell = HayCell()
    
    sim_params = {
                  "tstop": 400,
                  "dt": 0.1
                  }
    
    synaptic_input_params={"timings":[[100,150,200]], "locations":[(cell.secs["soma"][0].name(),0.5)], "weights": 5*0.0003, "Syn_type":syn_inh.Exp2syn_GABA}
    synaptic_input_params2={"timings":[[100,120,200]], "locations":[(cell.secs["soma"][0].name(),0.5)], "weights": 0.0003, "Syn_type":syn_exc.Exp2syn_AMPA}
    sim_params["synaptic_input"] = [synaptic_input_params,synaptic_input_params2]     
    trace = cell.simulate(sim_params)[0]
    if PLOT:
        plt.plot_trace(trace)
        pyl.show()

def test_additional_synaptic_parameters():
    cell = HayCell()
    
    sim_params = {
                  "tstop": 400,
                  "dt": 0.1
                  }
    
    timings = [[200,320], [100,350]]
    locations = [(cell.secs["soma"][0].name(),0.5),(cell.secs["soma"][0].name(),0.5)]
    asp = [{"tau1":5, "tau2":10}, {"tau1":0.1}]
    sim_params["synaptic_input"] = [{"timings": timings, "locations": locations, "additional_synaptic_parameters":asp}]
    trace = cell.simulate(sim_params)[0]
    print trace
    
    if PLOT:
        plt.plot_trace(trace)
        pyl.show()
    
    
test_simulate()
test_simulate_simple_input()
test_synaptic_input()
test_additional_synaptic_parameters()