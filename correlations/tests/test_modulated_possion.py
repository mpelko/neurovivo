"""
Testing suite for creating moduated poisson train
"""

import numpy as np
import neurovivo.correlations.correlated_input as ci

DEBUG = True

def test_get_sin_mod_parameters():
    freq = 3
    total_time = 10000
    mean_rate = 5
    mod_str = 0.5
    func, stim_count = ci.get_sin_mod_parameters(freq, total_time, mean_rate, mod_str)
    print stim_count
    assert stim_count == 67
    
    if DEBUG:
        import matplotlib.pylab as pyl
        time = np.arange(total_time)
        y= [func(i) for i in time]
        pyl.plot(time, y)
        pyl.show()
        
def test_modulated_poisson():
    freq = 3
    total_time = 1000
    mean_rate = 100
    mod_str = 0.9
    func, spike_tries = ci.get_sin_mod_parameters(freq, total_time, mean_rate, mod_str)
    spike_times = ci.modulated_poisson(func, spike_tries, total_time)
    if DEBUG:
        import matplotlib.pylab as pyl
        pyl.plot(spike_times,0.5*np.ones(len(spike_times)),"*")
        time = np.arange(total_time)
        y= [func(i) for i in time]
        pyl.plot(time,y)
        pyl.show()
        print 1000.*len(spike_times)/total_time, mean_rate
    
if __name__ == "__main__":
    test_get_sin_mod_parameters()
    test_modulated_poisson()