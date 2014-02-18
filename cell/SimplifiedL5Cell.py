'''
Created on Aug 13, 2011

@author: mpelko

Simplified morphology based on the morphology of HayCell, excludes all the channels except passive, 
calibrated to be
'''
from neurovivo.cell.NeuronCell import NeuronCell
from neuron import h
import numpy as np
import quantities as pq

class SimplifiedL5Cell(NeuronCell):
    '''
    A Simplified model of L5 cell. soma, basal, trunk and tuft compartments only.
    '''
    def __init__(self):
        
        super(self.__class__, self).__init__()
        self.name = "SimplifiedL5Cell"

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
            sec.e_pas = -90
            
        basal.insert('Ih')
        basal.gIhbar_Ih = 0.0002
        basal.cm = 2
        basal.g_pas = 0.0000467 

        soma.insert('Ca_LVAst') 
        soma.insert('Ca_HVA') 
        soma.insert('SKv3_1') 
        soma.insert('SK_E2') 
        soma.insert('K_Tst') 
        soma.insert('K_Pst') 
        soma.insert('Nap_Et2') 
        soma.insert('NaTa_t')
        soma.insert('CaDynamics_E2')
        soma.insert('Ih')
        soma.ek = -85
        soma.ena = 50
        soma.gIhbar_Ih = 0.0002
        soma.g_pas = 0.0000338 
        soma.decay_CaDynamics_E2 = 460.0 
        soma.gamma_CaDynamics_E2 = 0.000501 
        soma.gCa_LVAstbar_Ca_LVAst = 0.00343 
        soma.gCa_HVAbar_Ca_HVA = 0.000992 
        soma.gSKv3_1bar_SKv3_1 = 0.693 
        soma.gSK_E2bar_SK_E2 = 0.0441 
        soma.gK_Tstbar_K_Tst = 0.0812 
        soma.gK_Pstbar_K_Pst = 0.00223 
        soma.gNap_Et2bar_Nap_Et2 = 0.00172 
        soma.gNaTa_tbar_NaTa_t = 2.04 
        
        apical_secs = trunk_sections + [tuft]
        for sec in apical_secs:
            sec.insert('Ca_LVAst') 
            sec.insert('Ca_HVA') 
            sec.insert('SKv3_1') 
            sec.insert('SK_E2') 
            sec.insert('NaTa_t')
            sec.insert('CaDynamics_E2')
            sec.cm = 2.
            sec.insert('Ih')
            sec.insert('Im')  
            sec.ek = -85
            sec.ena = 50
            sec.decay_CaDynamics_E2 = 122 
            sec.gamma_CaDynamics_E2 = 0.000509 
            sec.gSK_E2bar_SK_E2 = 0.0012 
            sec.gSKv3_1bar_SKv3_1 = 0.000261 
            sec.gNaTa_tbar_NaTa_t = 0.0213 
            sec.gImbar_Im = 0.0000675 
            sec.g_pas = 0.0000589
            sec.gCa_LVAstbar_Ca_LVAst=0.000187
            sec.gCa_HVAbar_Ca_HVA=0.0000555
        
        def fill_ih_in_section(sec, value_at_zero_dist, value_at_max_dist, max_dist):
            min_dist = self.distance_to_soma(sec, position=0)
            max_dist = self.distance_to_soma(sec, position=1)
            distances = min_dist + np.arange(sec.L/(2.*sec.nseg),sec.L,sec.L/(sec.nseg)) * pq.um
            for seg, distance in zip(sec, distances):
                seg.gIhbar_Ih = value_at_zero_dist + 1.*value_at_max_dist*distance/max_dist
        
        
        self.soma = soma
        for sec in trunk_sections:
            fill_ih_in_section(sec, value_at_zero_dist=0.00025, value_at_max_dist=0.0022, max_dist=618*pq.um)
            sec.cm=3

        tuft.gIhbar_Ih = 0.0022

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