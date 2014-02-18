
#import numpy
#from numpy import *
import numpy

class SpikeFinderThreshCross(object):
    
    
    def __init__(self, trace, crossingthresh=0, firingthres=None):
        self.trace = trace
        self.crossingthresh = crossingthresh
        
        #Get the crossing times:
        threshIndices = self.findThresCrossings()
    
        # Make a spike for each one:
        self.spikes = [ Spike(self.trace, threshInd, firingthres=firingthres) for threshInd in threshIndices]
        
    def NumSpikes(self):
        return len(self.spikes)   
        
        
    def findThresCrossings(self):
        #t = self.trace.time
        d = self.trace._data.rescale("mV").magnitude
        
        aboveZero = numpy.zeros(d.shape, dtype=int)
        aboveZero[d > self.crossingthresh] = 1
        
        crossings = aboveZero - numpy.roll(aboveZero, 1)
        risingEdgeInd = numpy.where(crossings == 1)[0]  
        fallingEdgeInd = numpy.where(crossings == -1)[0]
        
        assert len(risingEdgeInd) == len(fallingEdgeInd)
        assert d[0] < self.crossingthresh
        
        threshIndices = zip(risingEdgeInd, fallingEdgeInd) 


        return threshIndices
        
        


class Spike(object):
    
    def getPeakSize(self):
        return self.trace._data[self.peakIndex]
    
    def __init__(self, trace, timeIndices, firingthres=None):
        self.trace = trace
        self.thresIndices = timeIndices
        self.firingthres = firingthres 
        
        self.init_getPeak()
        self.init_getDuration()
        
    def init_getPeak(self):
        d = numpy.copy(self.trace._data)
        d[ 0:self.thresIndices[0] ] = 0
        d[ self.thresIndices[1]:-1] = 0
        self.peakIndex = numpy.argmax(d) 
        
    
    
    def init_getDuration(self):
        
        self.fiftyPCLine = (self.trace._data.rescale("mV").magnitude[ self.peakIndex] + self.firingthres) / 2.0  
        
        d = numpy.copy(self.trace._data.rescale("mV").magnitude)
        
        d[ 0:self.thresIndices[0] ] = 0
        d[ self.thresIndices[1]:-1] = 0
        
        above50PC = numpy.zeros(d.shape, dtype=int)
        above50PC[d > self.fiftyPCLine] = 1
        
        crossings = above50PC - numpy.roll(above50PC, 1)
        risingEdgeInd = numpy.where(crossings == 1)  
        fallingEdgeInd = numpy.where(crossings == -1)  
        
        assert len(risingEdgeInd) == len(fallingEdgeInd) == 1
        
        self.durInd = risingEdgeInd, fallingEdgeInd
        self.duration = self.trace._time[fallingEdgeInd] - self.trace._time[risingEdgeInd]
        self.duration = self.duration.rescale("ms").magnitude
        
    def addToAxes(self, ax):
        t = self.trace._time
        d = self.trace._data.rescale("mV").magnitude
        
        # 50% Line:
        ax.plot((t[self.durInd[0]], t[self.durInd[1]]), (self.fiftyPCLine, self.fiftyPCLine), 'k--')

        # Peak Line:
        ax.plot((t[self.peakIndex], t[self.peakIndex]), (self.fiftyPCLine, d[self.peakIndex]), 'k:')

        
        # Annotate Plot
        print (t[self.peakIndex], d[self.peakIndex], self.duration, self.firingthres)
        print
        annotStr = """Time:%2.2fms \nPeak: %2.2f mV \nDur: %2.2f ms\n(ThresVoltage: %2.2f)""" % (t[self.peakIndex], d[self.peakIndex], self.duration, self.firingthres)
        ax.annotate(annotStr, xy=(t[self.peakIndex] + 2, d[self.peakIndex] - 10), xytext=None, xycoords='data', textcoords='data', arrowprops=None)

