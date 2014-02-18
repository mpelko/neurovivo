
from trace import Trace

import numpy
#import scipy.signal
from scipy.interpolate.interpolate import interp1d



class TraceSplicer(object):
    
    @classmethod
    def RebaseTime(cls, trace, newTimebase, **kwargs):
        bounds_error = False
        
        interpolator = interp1d(trace._time, trace.rawdata, "linear", bounds_error=bounds_error)
        newData = interpolator(newTimebase)
        rebasedTrace = Trace(newTimebase, newData, timeUnit=trace.timeUnit, dataUnit=trace.dataUnit, **kwargs)
        return rebasedTrace
    
    
    
    @classmethod
    def WindowTrace(cls, trace, timeWindow, **kwargs):
        
        if not isinstance(trace, Trace): raise ValueError()
        if timeWindow[0] - trace._time[0] < 0:
            print  "timeWindow[0]", timeWindow[0].rescale("s")
            print  "trace.time[0]", trace._time[0].rescale("s")
            print
            raise ValueError("Windowing outside of trace (min) WindowMin/TraceMin: %f %f  " % (timeWindow[0], trace.time[0]))
        
        #if timeWindow[1] > trace._time[-1]:
        if timeWindow[1] - trace._time[-1] > 0:
            print  "timeWindow[1]", timeWindow[1].rescale("s")
            print  "trace.time[-1]", trace._time[-1].rescale("s")
            print
            
            raise ValueError("Windowing outside of trace (max)")
        
        timeIndices1 = numpy.nonzero(trace._time > timeWindow[0])
        timeTraceNew = trace._time[timeIndices1]
        traceNew = trace._data[timeIndices1]
        
        timeIndices2 = numpy.nonzero(timeTraceNew < timeWindow[1])
        timeTraceNew = timeTraceNew[timeIndices2]
        traceNew = traceNew[timeIndices2]
        return Trace(timeTraceNew, traceNew, **kwargs)
    
    
    @classmethod
    def ShiftTrace(cls, trace, offset, **kwargs):
        return Trace(trace._time + offset, trace._data, **kwargs) 
    
    
    @classmethod
    def WindowAndShiftTrace(cls, trace, timeWindow, **kwargs):
        
        windowedTrace = cls.WindowTrace(trace, timeWindow)
        return cls.ShiftTrace(windowedTrace, -1.0 * windowedTrace._time[0], **kwargs)
    
    
    
