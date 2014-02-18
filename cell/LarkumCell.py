'''
Created on Aug 12, 2011

@author: mpelko
'''
from neurovivo.cell.NeuronCell import NeuronCell
import neuron
from neuron import h
import quantities as pq
from itertools import chain


class LarkumCell(NeuronCell):
    '''
    Larkum model:
    http://senselab.med.yale.edu/ModelDb/showmodel.asp?model=124043
    '''
    def __init__(self):
        
        super(self.__class__, self).__init__()
        self.name = "Larkum"

    def set_up(self):

        h.load_file(1, 'NEURON_stuff/LarkumStuff/as_for_python.hoc')

        # structuring the sections
        self.secs["basal"] = [h.dend[i] for i in xrange(70)]
        for sec in h.tuft:
            if sec.name() != h.apic[10].name():
                self.secs["tuft"].append(sec)

        hoc_trunk_full = h.SectionList()
        hoc_trunk_full.subtree(sec=h.apic[0])
        for sec in self.secs["tuft"]:
            hoc_trunk_full.remove(sec=sec)
        for sec in hoc_trunk_full:
            self.secs["trunk"].append(sec)

        hoc_trunk_full = None     # making sure the object gets destroyed.

        self.secs["axon"].append(h.axon)
        self.secs["axon"].append(h.axon1)
        self.secs["soma"].append(h.soma)

    def distance_to_soma(self, section):
        if section.name() == h.soma.name():
            return 0 * pq.um
        h.distance(sec=h.soma)
        distance = h.distance(0.5, sec=section)
        # make sure the distance is to 0.5 in soma
        if list(chain(self.sec_names["axon"], self.sec_names["basal"])).__contains__(section.name()):
            distance = distance + (h.soma.L/2)
        else: 
            distance = distance - (h.soma.L/2)
        return distance * pq.um
    
    def distance_to_dend(self, section):
        dend = h.apic[10]
        if section.name() == dend.name():
            return 0 * pq.um 
        h.distance(sec=dend)
        distance = h.distance(0.5, sec=section)
        # make sure the distance is to 0.5 in dend section
        if self.secs_names["tuft"].__contains__(section.name()):
            distance = distance - (dend.L/2)
        else:
            distance = distance + (dend.L/2)
        return distance * pq.um

    def destroy(self):
        """
        Takes care that neuron objects and python objects reffering to the cell
        are destroyed.

        At this point it actually destroys all the sections.
        """
        neuron.h('objref tuft')
        super(self.__class__, self).destroy()
