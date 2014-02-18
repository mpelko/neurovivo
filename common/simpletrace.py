'''
A lean version of the trace Class.
Data is always assumed to be in mV.
Time is always assumed to be in ms.

'''
import numpy as np

class SimpleTrace(object):
  
    def __init__(self, dt, data, name = "Unnamed"):
        self._name = name
        self._dt = dt
        self._data = data
   
    def getDT(self):
        return self._dt
    
    def getN(self):
        return len(self._data)
     
    def __str__(self):
        return "SimpleTrace object: " + self._name
        
    def __getitem__(self, time):
        assert time >= 0, "Simple trace starts with t = 0."
        assert time < len(self._data)*self._dt, "SimpleTrace not defined at requested time."
       
        from scipy.interpolate import interp1d         
        interpolator = interp1d(self._dt*np.arange(0,len(self._data)), self._data)
        
        dMag = interpolator(time)
        
        return dMag
    
    def mean(self):
        return np.mean(self._data)
       
    def std(self):
        return np.std( self._data)
    
    def cut_spikes(self, threshold = -45):
        self._data = self._data[self._data<-45]
    
    def __add__(self, rhs):
        return SimpleTrace( self._dt,  self._data + rhs)

    def __sub__(self, rhs):
        return SimpleTrace( self._dt,  self._data - rhs)
    
    def __div__(self, rhs):
        return SimpleTrace( self._dt,  self._data / rhs)
    
    def __mul__(self, rhs):
        return SimpleTrace( self._dt,  self._data * rhs)
    
    