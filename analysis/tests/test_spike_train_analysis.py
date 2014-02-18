'''
Created on Dec 29, 2011

@author: mpelko
'''

import neurovivo.analysis.spike_train_analysis as sta
import numpy as np
from neurovivo.common.spike_train import SpikeTrain
import matplotlib.pylab as pyl

st1 = SpikeTrain([250,600,1200,3300],start_time = 200, total_time=5200)
st2 = SpikeTrain([360,710,1300,3000],start_time = 200, total_time=5200)
st3 = SpikeTrain([460,510,1000,2202],start_time = 200, total_time=5200)
st4 = SpikeTrain([310,750,1120,1834],start_time = 200, total_time=5200)
st5 = SpikeTrain([333,888,1700,2800],start_time = 200, total_time=5200)

sts = [st1,st2,st3,st4,st5]

def test_correlation():
    print sta.correlation(st1, st2, dt=1, window=600)
    print sta.correlation(st1, st2, dt=1, window=10)
    print sta.correlation(st1, st2, dt=1, window=10)
    print sta.correlation(st1, st2, dt=1, window=2000)
    print sta.correlation(st1, st2, dt=4, window=600)
    print sta.correlation(st1, st2, dt=4, window=2000)
    windows = np.arange(1,1000,40)
    corrs = [sta.correlation(st1, st2, dt=1, window=w) for w in windows]
    pyl.plot(windows,corrs)

def test_average_pairwise_correlation_full():
    print sta.average_pairwise_correlation_full(sts,dt=1,window=20)
    windows = np.arange(1,1000,50)
    corrs = [sta.average_pairwise_correlation_full(sts, dt=1, window=w) for w in windows]
    pyl.plot(windows,corrs)

def test_average_3rd_order_correlation_full():
    print sta.average_3rd_order_correlation_full(sts, dt=1, window=20)
    windows = np.arange(1,200,5)
    corrs = [sta.average_3rd_order_correlation_full(sts, dt=1, window=w) for w in windows]
    pyl.plot(windows,corrs)
    
def test_average_nth_order_correlation():
    print sta.average_nth_order_correlation(sts, 3, dt=1, window=20, nr_combinations=1000)
    windows = np.arange(1,200,5)
    corrs = [sta.average_nth_order_correlation(sts, 3, dt=1, window=w, nr_combinations=1000) for w in windows]
    pyl.plot(windows,corrs)

def main():
    #test_correlation()
    #test_average_pairwise_correlation_full()
    test_average_3rd_order_correlation_full()
    test_average_nth_order_correlation()
    
if __name__ == "__main__":
    main()
    pyl.show()