'''
Created on Nov 29, 2009

@author: michael
'''
import numpy
import quantities as pq

class Trace(object):
  
    def __init__(self, time, data, name="Unnamed", label=None, tags=None, timeUnit=None, dataUnit=None):
        #self.label = label if label else "UnknownSrc"
        self.label = label
        
        if not isinstance(time, pq.quantity.Quantity): raise ValueError()
        if not isinstance(data, (pq.quantity.Quantity,pq.Dimensionless) ): raise ValueError()
        
        if not time.shape == data.shape: raise ValueError()
        
        self.tags = [] if tags == None else tags
        self._time = time
        self._data = data
        self.name = name
   
    def getDT(self):
        return float(self._time[1] - self._time[0])
    
    def getN(self):
        return len(self.time)
     
    def __str__(self):
        name = self.name
        if name == None: 
            name = "Unnamed"
        return "TraceObject: "+ name + " Shape:" + str(self._time.shape)
        
    def __getitem__(self, time):
        assert isinstance(time, pq.quantity.Quantity)
       
        # Rebase the Time:
        time.units = self._time.units
        
        from scipy.interpolate import interp1d         
        interpolator = interp1d(self._time.magnitude, self._data.magnitude)
        
        dMag = interpolator( time.magnitude)
        
        return dMag * self._data.units
    
    def mean(self):
        return numpy.mean( self._data )
       
    def stdDev(self):
        u = self._data.units
        
        stddev =  numpy.std( self._data.magnitude )
        return stddev * u

    def fft(self):
        return numpy.fft.fft(self._data.magnitude)
    
    def __add__(self, rhs):
        return Trace( self._time,  self._data + rhs)
    def __sub__(self, rhs):
        return Trace( self._time,  self._data - rhs)
    
    def __div__(self, rhs):
        return Trace( self._time,  self._data / rhs)
    
    def __mul__(self, rhs):
        return Trace( self._time,  self._data * rhs)
    

            
        
"""        
#Trace: 
#Trace Time-Unit: s
#Trace Data-Unit: V
#N-DataPoints:
0    0
6e-05    -0.0796464
0.00012    -0.0793053
0.00018    -0.0789763
0.00024    -0.0786589
0.0003    -0.0783527
0.00036    -0.0780574
0.00042    -0.0777725
0.00048    -0.0774976
0.00054    -0.0772326
"""
