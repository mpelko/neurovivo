'''
Created on Mar 29, 2012

@author: mpelko
'''
from neurovivo.cell.synapse.neuron_synapse import NeuronSynapse
from neuron import h

class alpha_synapse(NeuronSynapse): 
    
    def __init__(self, additional_parameters={}):
        super(self.__class__, self).__init__()
        
        self.name = "alpha_synapse"
        self.strength = 0.0003
        self.tau = 2.5 # in ms
        self.e = 0.    # in mV
        self.additional_parameters = additional_parameters
        
    def attach(self, location, section_name):
        
        if not self.synapse_for_neuron == None:
            if location == self.location and section_name == self.section_name:  
                return self.synapse_for_neuron
            else:
                assert False, "This synapse is already assigned to another section_name/location."
        self.synapse_for_neuron = h.Exp2Syn(location, sec=section_name)
        if "e" in self.additional_parameters:
            self.synapse_for_neuron.e = self.additional_parameters["e"]
        else:
            self.synapse_for_neuron.e = self.e
        if "tau" in self.additional_parameters:
            self.synapse_for_neuron.tau1 = self.additional_parameters["tau"]
            self.synapse_for_neuron.tau2 = self.additional_parameters["tau"]
        else:
            self.synapse_for_neuron.tau1 = self.tau
            self.synapse_for_neuron.tau2 = self.tau
        self.location = location
        self.section_name = section_name