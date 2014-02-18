'''
Created on Mar 29, 2012

@author: mpelko
'''

class NeuronSynapse(object):
    """
    The template class for any Neuron synapse.
    """
    def __init__(self):
        self.name = "Template"
        self.strength = 0.
        self.delay = 0.
        self.synapse_for_neuron = None
        self.location = None
        self.section_name = None
    
    def attach(self, location, section):
        assert False, "Not inplemented" 