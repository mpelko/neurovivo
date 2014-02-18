'''
Created on Aug 13, 2011

@author: mpelko

Simplified morphology based on the HayCell. 
Passive, only includes the pas mechanism, the rest is same as SimplifiedL5Cell.
'''

from neurovivo.cell.NeuronCell import NeuronCell
from neuron import h
import quantities as pq

class SimplifiedL5CellPassive(NeuronCell):
    '''
    A Simplified model of L5 cell. Passive.
    '''
    def __init__(self):
        
        super(self.__class__, self).__init__()
        self.name = "SimplifiedL5CellPassive"

    def set_up(self):

        # Simplified morphology based on the morphology of HayCell 
        soma = h.Section(name="soma")
        soma.L = 23
        soma.diam = 13.5
        
        trunk1 = h.Section(name="trunk1")
        trunk1.L = 20
        trunk1.nseg = 2
        trunk1.diam = 5.5
        
        trunk2 = h.Section(name="trunk2")
        trunk2.L = 605
        trunk2.nseg = 30
        trunk2.diam = 2.6
        
        trunk_sections = [trunk1, trunk2]

        basal = h.Section(name="basal")
        basal.L = 750
        basal.diam = 6.5
        
        tuft = h.Section(name="tuft")
        tuft.L = 1200
        tuft.diam = 2.8
        
        trunk1.connect(soma)
        trunk2.connect(trunk1)
        basal.connect(soma,0,0)
        tuft.connect(trunk2)
        
        all_secs = [soma, basal, tuft]
        all_secs.extend(trunk_sections)
        
        for sec in all_secs:
            sec.insert("pas")
            sec.Ra = 100
            sec.cm = 1
            sec.e_pas = -70

        basal.cm = 2
        basal.g_pas = 0.0000467 
        
        soma.g_pas = 0.0000338 

        apical_secs = trunk_sections + [tuft]
        for sec in apical_secs:
            sec.cm = 2.
            sec.g_pas = 0.0000589

        #for sec in trunk_sections:
        #    sec.cm=3

        # structuring the sections
        self.secs["basal"].append(basal)
        self.secs["soma"].append(soma)
        self.secs["trunk"] = trunk_sections = [trunk1, trunk2]
        self.secs["tuft"].append(tuft)
                
        self.bifurcation_info = ("trunk2", 1)

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