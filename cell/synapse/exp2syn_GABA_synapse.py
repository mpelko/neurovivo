'''
Created on Mar 29, 2012

@author: mpelko
'''
from neurovivo.cell.synapse.neuron_synapse import NeuronSynapse
from neuron import h

class Exp2syn_GABA(NeuronSynapse): 

    def __init__(self, additional_parameters={}):
        super(self.__class__, self).__init__()
        
        self.name = "exp2syn_GABA"
        self.strength = 0.0003
        self.tau1 = 1 # from Hausser and Roth 1997
        self.tau2 = 10 # from Ali et al. 2001
        self.additional_parameters = additional_parameters
        
    def attach(self, location, section_name):
        
        if not self.synapse_for_neuron == None:
            if location == self.location and section_name == self.section_name:  
                return self.synapse_for_neuron
            else:
                assert False, "This synapse is already assigned to another section_name/location."
        self.synapse_for_neuron = h.Exp2Syn(location, sec=section_name)
        self.synapse_for_neuron.e = -80.
        if "tau1" in self.additional_parameters:
            self.synapse_for_neuron.tau1 = self.additional_parameters["tau1"]
        else:
            self.synapse_for_neuron.tau1 = self.tau1
        if "tau2" in self.additional_parameters:
            self.synapse_for_neuron.tau2 = self.additional_parameters["tau2"]
        else:
            self.synapse_for_neuron.tau2 = self.tau2
        self.location = location
        self.section_name = section_name
