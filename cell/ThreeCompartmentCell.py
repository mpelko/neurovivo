'''
Created on 9 Sep 2011

@author: mpelko
'''
from neurovivo.cell.NeuronCell import NeuronCell
from neuron import h

class ThreeCompartmentCell(NeuronCell):
    '''
    The cell with soma and two dendrites of same dimensions hanging out of the soma.
    We call them trunk and basal, but they have same properties.
    '''
    def __init__(self):
        
        super(self.__class__, self).__init__()
        self.name = "ThreeCompartment"

    def set_up(self):

        # Simplified morphology based on the morphology of HayCell 
        soma = h.Section(name="soma")
        soma.L = 23
        soma.diam = 13.5
        
        trunk = h.Section(name="trunk1")
        trunk.L = 200
        trunk.nseg = 10
        trunk.diam = 2
                
        basal = h.Section(name="basal")
        basal.L = 200
        basal.nseg = 10
        basal.diam = 2
                
        trunk.connect(soma)
        basal.connect(soma,0,0)
        
        all_secs = [soma, basal, trunk]
        
        for sec in all_secs:
            sec.insert("pas")
            sec.Ra = 200
            sec.cm = 1
            sec.e_pas = -70
            sec.v = -70

        # structuring the sections
        self.secs["basal"].append(basal)
        self.secs["soma"].append(soma)
        self.secs["trunk"].append(trunk)
                
        self.bifurcation_info = ("trunk", 1)
        
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