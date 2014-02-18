import numpy as np
from neuron import h
from neuron import gui
from networks.artificial_cell import ArtificialCell
from networks.targetLocation import TargetLocation
from networks.inputConnection import InputConnection

import matplotlib.pyplot as pyl

def test_simple_connection():
    
    soma = h.Section()
    ac = ArtificialCell()
    tl = TargetLocation(soma, 0.5)
    proc = h.ExpSyn 

    input_connections = []
     
    input_connections.append(InputConnection(ac, tl, proc))

def test_simple_synaptic_input():
    soma = h.Section()
    soma.insert("hh")
    soma.insert("pas")
    tl = TargetLocation(soma, 0.5)
    proc = h.ExpSyn 

    ac = ArtificialCell()
    fire_times = np.linspace(0,50,10)
    firetimes = h.Vector(fire_times)
    ac.setFireTimes(firetimes)
    
    ic = InputConnection(ac, tl, proc)

    v = h.Vector()
    v.record(soma(0.5)._ref_v)
   
    h.tstop = 100
    h.run()

    test_result = np.array(v)
    correct = np.load("test_simple_synaptic_input_data.npy")

    print sum(abs(test_result - correct))

    #pyl.plot(test_result)
    #pyl.plot(correct)
    #pyl.show()

    assert np.all(test_result == correct)
