import numpy as np
from neuron import h
#from neuron import gui

class TargetLocation(object):

    def __init__(self, section, sectionpos):
        self.section = section
        self.sectionpos = sectionpos
        assert 0 <= sectionpos <= 1

    def distanceFromSoma(self):
        #TODO - implement
        assert False, "Not implemented"
