import numpy as np
from neuron import h
#from neuron import gui

class InputConnection(object):
    """
    A class mapping the artificial cells providing the input spikes to the
    target locations on the cell.
    Creates the neuron NetCon object on the way.
    Requires the "source" artificial cell, the location on the target_cell
    (target_lcoation) and the process (Neuron POINT_PROCESS), specified in the
    .mod file or Neuron native.
    """
    rec = h.Vector()
    
    def __init__(self, artificial_cell, target_location, process_function,\
    weight=1):
        
        self.source = artificial_cell
        self.target = target_location
        self.process =\
        process_function(target_location.section(target_location.sectionpos),
        name = "syn", sec = target_location.section)
        self.weight = weight

        self.nc = h.NetCon(self.source.vecStim, self.process)
        self.nc.weight[0] = weight
        #self.nc.event(10)
        #self.nc.event(21)
        #self.nc.event(62)
        #self.nc.event(93)

