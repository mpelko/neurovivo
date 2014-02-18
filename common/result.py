from morphforge.core import LocMgr, WriteToFile, FilterExpectSingle
import pickle

class Result(object):
    
    """ traces is a list of Trace Objects"""
    def __init__(self, traces, name = "Undefined simulation"):
   
        self.traces = traces
        self.name = name
        self.tStart = None
        self.tStop = None
    
    def setSimulationTime(self, tStart, tStop):
        self.tStart = tStart
        self.tStop = tStop
    
    def getTrace(self, name):
        return FilterExpectSingle(self.traces, lambda s: s.name == name) 
    
    # Loading & Saving:
    def saveToFile(self, filename):
        resString = pickle.dumps(self)
        return WriteToFile(resString, filename=filename, filedirectory=LocMgr.getSimulationTmpDir())
         
    @classmethod
    def loadFromFile(cls, filename):
        return pickle.load(open(filename))
