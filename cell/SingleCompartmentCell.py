'''
Created on 9 Sep 2011

@author: mpelko
'''
from neurovivo.cell.NeuronCell import NeuronCell
from neuron import h

class SingleCompartmentCell(NeuronCell):
    '''
    Single compartment cell with HH and pas mechanism.
    The other parameters are the neuron defaults. 
    '''
    def __init__(self):
        
        super(self.__class__, self).__init__()
        self.name = "SingleCompartment"

    def set_up(self):

        soma = h.Section()
        soma.insert("hh")
        soma.insert("pas")
        soma.L = 100
        soma.diam = 100
        soma(0.5).pas.e = -80
        
        self.secs["soma"].append(soma)
        
    def destroy(self):
        """
        Takes care that neuron objects and python objects reffering to the cell
        are destroyed.

        At this point all the cell sections should be destroyed.
        
        Make sure all the references to the sections are deleted.
        """
        for key in self.secs.keys():
            self.secs[key] = []
            self.sec_names[key] = []