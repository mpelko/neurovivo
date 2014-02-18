'''
Created on Aug 12, 2011

@author: mpelko
'''

import brian as br
import brian.library.synapses as brsy
import numpy as np

class LIFCell(object):
    '''
    A class representing a LIF neuron parameters. Used as a default parameter collection.
    '''
    def __init__(self, N_cells):
        self.description = \
        """
        A default LIF cell with conducing synapses.
        """
        self.C = 605. * br.pF # based on HayCell surface with spines and 1uF/cm^2 relative capacitance
        self.taum = 24 * br.ms # based on Rm ~= 40 MOhm, number given by Paolo 
        self.gL = self.C / self.taum
        self.EL = -80 * br.mV # leak reversal potential Hay actually uses -90
        self.Eres = -65 * br.mV # the reset potential
        self.VT = -55 * br.mV # threshold
        self.Tref = 5 * br.ms # refractory period
        
        # Synapse parameters
        self.Ee = 0 * br.mV
        self.Ei = -80 * br.mV
        self.taue1 = 0.2 * br.ms  # from Hausser and Roth 2007
        self.taue2 = 1.7 * br.ms  # from Hausser and Roth 2007
        self.taui1 = 3 * br.ms
        self.taui2 = 3.01 * br.ms # can't be the same as taui1 - division by zero
        
        self.eqs = 'dvm/dt=(self.gL*(self.EL-vm)+ge*(self.Ee-vm)+gi*(self.Ei-vm)) / self.C :volt'
        self.eqs += brsy.biexp_conductance('gi', self.Ei, tau1=self.taui1, tau2=self.taui2) # from library.synapses
        self.eqs += brsy.biexp_conductance('ge', self.Ee, tau1=self.taue1, tau2=self.taue2) # from library.synapses
        
        self.conductances = ["ge", "gi"]
        
        self.group = br.NeuronGroup(N_cells, model=self.eqs, threshold=self.VT, reset=self.Eres, refractory=self.Tref)