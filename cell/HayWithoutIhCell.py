'''
Created on Aug 13, 2011

@author: mpelko
'''
from neurovivo.cell.NeuronCell import NeuronCell
from neuron import h
import quantities as pq
from itertools import chain


class HayWithoutIhCell(NeuronCell):
    '''
    Hay model:
    http://senselab.med.yale.edu/ModelDb/ShowModel.asp?model=139653
    '''
    def __init__(self):
        
        super(self.__class__, self).__init__()
        self.name = "HayWithoutIh"
        self.tuft_info = None

    def set_up(self):

        h.load_file(1, 'NEURON_stuff/HayStuff/set_up.hoc')

        # structuring the sections
        for sec in h.L5PC.basal:
            self.secs["basal"].append(sec)
            sec.gIhbar_Ih = 0.
            
        for sec in h.L5PC.somatic:
            self.secs["soma"].append(sec)
            sec.gIhbar_Ih = 0.

        for sec in h.L5PC.axonal:
            self.secs["axon"].append(sec)
        
        # spliting tuft and trunk part of apical    
        hoc_tuft = h.SectionList()
        hoc_tuft.subtree(sec=h.L5PC.apic[36])
        hoc_trunk = h.SectionList()
        for sec in h.L5PC.apical:
            sec.gIhbar_Ih = 0.
            hoc_trunk.append(sec=sec)
        
        for sec in hoc_tuft:
            if sec.name() != h.L5PC.apic[36].name():
                self.secs["tuft"].append(sec)

        for sec in self.secs["tuft"]:
            hoc_trunk.remove(sec=sec)

        for sec in hoc_trunk:
            self.secs["trunk"].append(sec)
            
        hoc_tuft = None     # making sure the object gets destroyed.
        
        self.bifurcation_info = (self.sections("trunk")[36].name(), 1)

    def distance_to_soma(self, section, position=0.5):
        soma = self.soma
        h.distance(0,0.5,sec=soma)
        distance = h.distance(position, sec=section)
        return distance * pq.um
    
    def distance_to_main_bifurcation(self, section, position=0.5):
        dend = self.section(self.bifurcation_info[0])
        h.distance(0,self.bifurcation_info[1],sec=dend)
        distance = h.distance(position, sec=section)
        return distance * pq.um

    def destroy(self):
        """
        Takes care that neuron objects and python objects reffering to the cell
        are destroyed.

        At this point it actually destroys all the sections.
        """
        super(self.__class__, self).destroy()