import numpy as np
import neurovivo.analysis.trace_analysis as ta
import quantities as pq
from common.trace import Trace

DEBUG = True 

def test_spike_times_fast():
    voltage1 = [2,3,1,-1,-4,-6,-2,3,5,-3,-1,2]*pq.mV
    voltage_trace1 = Trace(np.array(range(len(voltage1)))*pq.ms, voltage1)
    voltage2 = [-2,3,1,-1,-4,-6,-2,3,5,-3,-1,2]*pq.mV
    voltage_trace2 = Trace(np.array(range(len(voltage2)))*pq.ms, voltage2)
    voltage3 = [-2,-2]*pq.mV
    voltage_trace3 = Trace(np.array(range(len(voltage3)))*pq.ms, voltage3)
    voltage4 = []*pq.mV
    voltage_trace4 = Trace(np.array(range(len(voltage4)))*pq.ms, voltage4)
    voltage5 = [-1,-1, 1]*pq.mV
    voltage_trace5 = Trace(np.array(range(len(voltage5)))*pq.ms, voltage5)
    voltage6 = [1,-1, 1]*pq.mV
    voltage_trace6 = Trace(np.array(range(len(voltage6)))*pq.ms, voltage6)
    voltage7 = [1,-1, -1]*pq.mV
    voltage_trace7 = Trace(np.array(range(len(voltage7)))*pq.ms, voltage7)
    voltage8 = [-1,-1, -1]*pq.mV
    voltage_trace8 = Trace(np.array(range(len(voltage8)))*pq.ms, voltage8)
    voltage9 = [1,-1]*pq.mV
    voltage_trace9 = Trace(np.array(range(len(voltage9)))*pq.ms, voltage9)
    
    true_result1 = np.array([7,11])*pq.ms
    true_result2 = np.array([1,7,11])*pq.ms
    true_result3 = np.array([])*pq.ms
    true_result4 = np.array([])*pq.ms
    true_result5 = np.array([2])*pq.ms
    true_result6 = np.array([2])*pq.ms
    true_result7 = np.array([1])*pq.ms
    true_result8 = np.array([])*pq.ms
    true_result9 = np.array([])*pq.ms

    if DEBUG:
        print ta.spike_times_fast(voltage_trace1)
        print ta.spike_times_fast(voltage_trace2)
        print ta.spike_times_fast(voltage_trace3)
        print ta.spike_times_fast(voltage_trace4)
        print ta.spike_times_fast(voltage_trace5)
        print ta.spike_times_fast(voltage_trace6)
        print ta.spike_times_fast(voltage_trace7)
        print ta.spike_times_fast(voltage_trace8)
        print ta.spike_times_fast(voltage_trace9)

    assert np.all(ta.spike_times_fast(voltage_trace1) == true_result1)
    assert np.all(ta.spike_times_fast(voltage_trace2) == true_result2)
    assert np.all(ta.spike_times_fast(voltage_trace3) == true_result3)
    assert np.all(ta.spike_times_fast(voltage_trace4) == true_result4)
    assert np.all(ta.spike_times_fast(voltage_trace5) == true_result5)
    assert np.all(ta.spike_times_fast(voltage_trace6) == true_result6)
    assert np.all(ta.spike_times_fast(voltage_trace7) == true_result7)
    assert np.all(ta.spike_times_fast(voltage_trace8) == true_result8)

    assert np.all(ta.spike_times_fast(voltage_trace1) ==\
                 ta.spike_times(voltage_trace1))
    assert np.all(ta.spike_times_fast(voltage_trace2) ==\
                 ta.spike_times(voltage_trace2))
    assert np.all(ta.spike_times_fast(voltage_trace3) ==\
                 ta.spike_times(voltage_trace3))
    assert np.all(ta.spike_times_fast(voltage_trace4) ==\
                 ta.spike_times(voltage_trace4))
    assert np.all(ta.spike_times_fast(voltage_trace5) ==\
                 ta.spike_times(voltage_trace5))
    assert np.all(ta.spike_times_fast(voltage_trace6) ==\
                 ta.spike_times(voltage_trace6))
    assert np.all(ta.spike_times_fast(voltage_trace7) ==\
                 ta.spike_times(voltage_trace7))
    assert np.all(ta.spike_times_fast(voltage_trace8) ==\
                 ta.spike_times(voltage_trace8))

def test_spike_times():
    voltage1 = [2,3,1,-1,-4,-6,-2,3,5,-3,-1,2]*pq.mV
    voltage_trace1 = Trace(np.array(range(len(voltage1)))*pq.ms, voltage1)
    true_result = np.array([7,11])*pq.ms
    
    assert np.all(ta.spike_times(voltage_trace1) == true_result)

def test_was_spike_elicited():
    voltage1 = [2,3,1,-1,-4,-6,-2,3,5,-3,-1,2]*pq.mV
    voltage_trace1 = Trace(np.array(range(len(voltage1)))*pq.ms, voltage1)
    voltage2 = [-2,3,1,-1,-4,-6,-2,3,5,-3,-1,2]*pq.mV
    voltage_trace2 = Trace(np.array(range(len(voltage2)))*pq.ms, voltage2)
    voltage3 = [-2,-2]*pq.mV
    voltage_trace3 = Trace(np.array(range(len(voltage3)))*pq.ms, voltage3)
    voltage4 = []*pq.mV
    voltage_trace4 = Trace(np.array(range(len(voltage4)))*pq.ms, voltage4)
    voltage5 = [-1,-1, 1]*pq.mV
    voltage_trace5 = Trace(np.array(range(len(voltage5)))*pq.ms, voltage5)
    voltage6 = [1,-1, 1]*pq.mV
    voltage_trace6 = Trace(np.array(range(len(voltage6)))*pq.ms, voltage6)
    voltage7 = [1,-1, -1]*pq.mV
    voltage_trace7 = Trace(np.array(range(len(voltage7)))*pq.ms, voltage7)
    voltage8 = [-1,-1, -1]*pq.mV
    voltage_trace8 = Trace(np.array(range(len(voltage8)))*pq.ms, voltage8)
    voltage9 = [1,-1]*pq.mV
    voltage_trace9 = Trace(np.array(range(len(voltage9)))*pq.ms, voltage9)
    
    true_result1 = True 
    true_result2 = True 
    true_result3 = False 
    true_result4 = False 
    true_result5 = True 
    true_result6 = True 
    true_result7 = False 
    true_result8 = False 
    true_result9 = False 

    if DEBUG:
        print ta.was_spike_elicited(voltage_trace1)
        print ta.was_spike_elicited(voltage_trace2)
        print ta.was_spike_elicited(voltage_trace3)
        print ta.was_spike_elicited(voltage_trace4)
        print ta.was_spike_elicited(voltage_trace5)
        print ta.was_spike_elicited(voltage_trace6)
        print ta.was_spike_elicited(voltage_trace7)
        print ta.was_spike_elicited(voltage_trace8)
        print ta.was_spike_elicited(voltage_trace9)

    assert np.all(ta.was_spike_elicited(voltage_trace1) == true_result1)
    assert np.all(ta.was_spike_elicited(voltage_trace2) == true_result2)
    assert np.all(ta.was_spike_elicited(voltage_trace3) == true_result3)
    assert np.all(ta.was_spike_elicited(voltage_trace4) == true_result4)
    assert np.all(ta.was_spike_elicited(voltage_trace5) == true_result5)
    assert np.all(ta.was_spike_elicited(voltage_trace6) == true_result6)
    assert np.all(ta.was_spike_elicited(voltage_trace7) == true_result7)
    assert np.all(ta.was_spike_elicited(voltage_trace8) == true_result8)
