'''
Created on Aug 12, 2011

@author: mpelko
'''

import neuron
from neuron import h
from itertools import chain
import neurovivo.common as cmn
from neurovivo.cell.synapse import exp2syn_AMPA_synapse as syn
from neurovivo.networks.artificial_cell import ArtificialCell
from neurovivo.common import Trace

import numpy as np
import quantities as pq

def destroy_previous():
    """
    Only one cell can exist at a time. Destroy all previous cells (sections).
    """
    neuron.h('forall delete_section()')
        
    for sec in h.allsec():
        assert False, "destroy_previous did not destroy all the sections. Revise. {}".format(sec.name())

class NeuronCell(object):
    '''
    A class representing a neuron. Used as a tempalte for various cells.
    '''
    def __init__(self):
        
        # We only allow for one cell to exist at a time
        destroy_previous()
        
        # Containers for hoc or python Neuron sections.
        # All but soma can stay empty after the set_up.
        self.name = "Template"
        self.secs = {
                         "axon":[],
                         "soma":[],
                         "basal":[],
                         "trunk":[],
                         "tuft":[],
                         }
        
        self.set_up()

        # Containers for names - needed as the 
        # sections can not be compared for some reason.
        self.sec_names = {
                         "axon":[],
                         "soma":[],
                         "basal":[],
                         "trunk":[],
                         "tuft":[],
                         }
        for sec in self.secs["basal"]:
            self.sec_names["basal"].append(sec.name())
        for sec in self.secs["trunk"]:
            self.sec_names["trunk"].append(sec.name())
        for sec in self.secs["axon"]:
            self.sec_names["axon"].append(sec.name())
        for sec in self.secs["soma"]:
            self.sec_names["soma"].append(sec.name())
        for sec in self.secs["tuft"]:
            self.sec_names["tuft"].append(sec.name())
         
        assert len(self.secs["soma"]) > 0,\
             "Check your {0} cell. No soma after setup.".format(self.name)
        self.soma = self.secs["soma"][0]
             
    def random_section(self, region=None, np_rnd = None):
        """
        Randomly selects a section. The probablities are weighted according to
        section lengths.
        """
        assert region == None or region in self.secs.keys(), "Invalid region: {0}".format(region)

        if not np_rnd:
            np_rnd = np.random
        
        sections = self.sections(region)
        lengths = [section.L for section in sections]
        
        return cmn.random_choice(sections, weights=lengths, np_rnd = np_rnd)
         

    def set_up(self):
        """
        Sets up the cell by creating the sections, filling it with channels, ...
        """
        assert False, "Not implemented error - you should write a set_up method for your cell"
        
        
    def distance_to_soma(self, section):
        """
        returns the distance from the middle of the section to soma.
        """
        assert False, "Not implemented error - you should write the method for your cell"

    def destroy(self):
        """
        Takes care that neuron objects and python objects reffering to the cell
        are destroyed.

        At this point all the cell sections should be destroyed.
        
        Make sure all the references to the sections are deleted
        """
        for key in self.secs.keys():
            self.secs[key] = []
            self.sec_names[key] = []
        self.artificial_cells = []
        self.synapses = []
        self.net_cons = []
        print "deleting all sections"
        neuron.h('forall delete_section()')
    
    def display(self):
        assert False, "Not implemented error."
        
    def sections(self, region=None):
        """
        Returns the sections in the region of a cell. 
        If region is not specified it returns all the sections in the cell.
        """
        if region == None:
            return list(chain(*[self.secs[key] for key in self.secs.keys()]))
        else:
            return self.secs[region]
        
    def section(self, sec_name):
        """
        Returns the first section with the section name sec_name.
        Doesn't care if there are several such sections, returns 
        the first one - could be random.
        """
        for sec in self.sections():
            if sec.name() == sec_name:
                return sec
        assert False, "No section with name {}.".format(sec_name)
        
    def classify_location(self, section):
        if section.name() in self.sec_names["basal"]:
            return "basal"
        if section.name() in self.sec_names["tuft"]:
            return "tuft"
        if section.name() in self.sec_names["trunk"]:
            return "trunk"
        if section.name() in self.sec_names["soma"]:
            return "basal"
        if section.name() in self.sec_names["axon"]:
            return "axon"
    
    def section_area(self, sec):
        """
        Returns the surface area of the section sec.
        """
        points = np.arange(1./(2*sec.nseg),1.,1./sec.nseg)
        return np.sum([h.area(p,sec=sec) for p in points])    
    
    def surface_area(self):
        """
        Returns the surface area of the cell.
        """
        return np.sum([self.section_area(sec) for sec in self.sections()])
        
    def surface_area_with_spines(self):
        """
        Returns the surface area of the cell. If the membrane capaticances 
        are changed in the dendrites to account for the spines, we take this
        into account. 
        """
        total_area = 0
        cm_soma = self.secs["soma"][0].cm
        for sec in self.sections():
            cm = sec.cm
            spine_correction_factor = 1.*cm/cm_soma 
            total_area += self.section_area(sec)*spine_correction_factor
        return total_area
    
    def simulate(self, simulation_parameters={}):
        """
        The grand function for simulation. Returns the trace.
        simulation_parameters can contain:
        - synaptic_input
            * timings - required
            * locations - required either as section names only or list of tuples (section_name, location)
            * Syn_type - default syn.Exp2syn_AMPA
            * weights
        - current_injection (only a box-type pulse allowed)
            * location
            * tstart
            * duration
            * magnitude
        - dt (in ms) - default whatever is Neuron simulator default
        - recording (a tuple of (section_name, location) - only the voltage can be recorded)
        - tstop (in ms) - default 500 ms
        - init_volage (in mV)
        """
        
        results = []
        recordings = []
        
        sp = simulation_parameters
        keys = sp.keys()
        
        if "synaptic_input" in keys:
            syns = []
            acs = []
            ncs = []
            firetimes_hoc=[]
            
            for ig_ind, input_group in enumerate(sp["synaptic_input"]):
                # Unfortunatelly the following line doesn't work - some Neuron-python memory probles
                # it seems to be related to netCon objects which somehow get screwed up as soon as they change the namespace 
                #self.add_syn_input(**input_group)
                
                # so we have to spell it out in this function - internal need to speak to Neuron it seems
                timings = input_group["timings"]
                locations = input_group["locations"]
                
                #assert len(timings) > 0 and len(locations)>0

                # setting the default location on the section, if not provided
                if not isinstance(locations[0], (list, tuple)): 
                    locations = [(loc, 0.5) for loc in locations]
                
                if not "Syn_type" in input_group.keys():
                    Syn_type=syn.Exp2syn_AMPA
                else:
                    Syn_type=input_group["Syn_type"]

                if not "additional_synaptic_parameters" in input_group.keys():
                    additional_synaptic_parameters = [{} for _ in xrange(len(timings))]
                else:
                    additional_synaptic_parameters=input_group["additional_synaptic_parameters"]
                
                if not "weights" in input_group.keys():
                    weights=None
                else:
                    weights=input_group["weights"]
                    
                # making sure we have timings and locations of same size
                assert len(timings) == len(locations), "len_timings = {}, len_locations={} - uncool".format(len(timings), len(locations))
                
                #synapse = Syn_type()
                
                # making weights into a list of length timings
                if weights==None:
                    synapse = Syn_type()
                    weights = [synapse.strength for _ in xrange(len(timings))]
                try:
                    iter(weights)
                except:
                    weights = [weights for _ in xrange(len(timings))]
                
                # transforming the timings in Neuron-friendly vector format
                firetimes_hoc.append([])
                for ft in timings:
                    firetimes_hoc[ig_ind].append(h.Vector(ft))
         
                for i, inp in enumerate(xrange(len(firetimes_hoc[ig_ind]))):
                    synapse = Syn_type()
                    synapse.additional_parameters = additional_synaptic_parameters[i] 
                    section = self.section(locations[inp][0])
                    acs.append(ArtificialCell())
                    acs[-1].setFireTimes(firetimes_hoc[ig_ind][inp])
                    synapse.attach(locations[inp][1], section_name=section)
                    syns.append(synapse)
                    weight = weights[inp]
                    nc = h.NetCon(acs[-1].vecStim, synapse.synapse_for_neuron, 10, 0, weight)
                    ncs.append(nc)

        if "current_injection" in keys:
            ics = []
            for curr_inj in sp["current_injection"]:
                location = curr_inj["location"]
                if not isinstance(location, (list, tuple)): 
                    location = (location, 0.5)
                #location = (self.section(curr_inj["location"]), 0.5)
                ic = h.IClamp(self.section(location[0])(location[1]))
                ics.append(ic)
                ics[-1].delay = curr_inj["tstart"]
                ics[-1].dur = curr_inj["duration"]
                ics[-1].amp = curr_inj["magnitude"]


        if "zap" in keys:
            # zap model taken from (http://paynesnotebook.net/Research/NeuronSimulation/Notes/10/index.html)
            zaps = []
            for zap in sp["zap"]:
                location = zap["location"]
                if not isinstance(location, (list, tuple)):
                    location = (location, 0.5)
                #location = (self.section(curr_inj["location"]), 0.5)
                zapy = h.Izap(self.section(location[0])(location[1]))
                zaps.append(zapy)
                zaps[-1].Ioff = zap.get("Ioff",0)
                zaps[-1].Astart = zap.get("Astart",0)
                zaps[-1].Astop = zap.get("Astop",0)
                zaps[-1].ttstart = zap.get("ttstart",0)
                zaps[-1].ttstop = zap.get("ttstop",0)
                zaps[-1].Fstart = zap.get("Fstart",0)
                zaps[-1].Fstop = zap.get("Fstop",0)

        if "dt" in keys:
            h.dt = sp["dt"]
            
        if "recording" in keys:
            for i,rec in enumerate(sp["recording"]):
                recordings.append(h.Vector())
                section = self.section(rec[0])
                recordings[i].record(section(rec[1])._ref_v)
        else:
            recordings.append(h.Vector())
            section = self.secs["soma"][0]
            recordings[0].record(section(0.5)._ref_v)
                
        if "tstop" in keys:
            tstop = sp["tstop"]
        else:
            tstop = 500
        
        #import time
        #time.sleep(2)
                    
        #if "init_voltage" in keys:
        #    h.finitialize(sp["init_voltage"])
        #else:
        #    h.finitialize(-65)
        
        #print "Starting the run."
        neuron.init()
        neuron.run(tstop)
                
        # need to do some cleaning up for further runs and to prevent 
        # segmentation fault 
        nc = None
        ncs = []
        
        for i,rec in enumerate(recordings):
            voltage = np.array(rec)
            trace = Trace(np.arange(len(voltage))*h.dt*pq.ms, voltage*pq.mV)
            results.append(trace)
        
        return results

    def print_statistics(self):
        print "Cell name: {0}".format(self.name)
        print "Sections:"
        print "somatic: {0}".format(len(self.sections("soma")))
        print "axonal: {0}".format(len(self.sections("axon")))
        print "basal: {0}".format(len(self.sections("basal")))
        print "trunk: {0}".format(len(self.sections("trunk")))
        print "tuft: {0}".format(len(self.sections("tuft")))
        print "total surface: {0}".format(self.surface_area())