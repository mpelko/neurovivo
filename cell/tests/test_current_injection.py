'''
Created on Mar 29, 2012

@author: mpelko
'''

from neurovivo.cell.HayCell import HayCell
import matplotlib.pylab as pyl
import neurovivo.plotting as plt
from neurovivo.cell.synapse import exp2syn_GABA_synapse as syn_inh
from neurovivo.cell.synapse import exp2syn_AMPA_synapse as syn_exc

PLOT = True
    

def test_simulate_current_injection():
    cell = HayCell()
    
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
    
test_simulate_current_injection()
