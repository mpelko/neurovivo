"""
Creating correlated input trains
"""

import numpy as np
import random
from neurovivo.correlations.correlated_input import *
from neurovivo.correlations.correlation_analysis import *
import cProfile

DEBUG = False

def test_population_pairwise_correlation():
    trains = np.array([[1,0,1],
                       [0,0,1]])

    true_result = [0.5, 0]
    calculated_result = population_pairwise_pearson_correlation(trains)

    if DEBUG:
        print "test_population_pairwise_correlation:"
        print "Calc: " + str(calculated_result)
        print "True: " + str(true_result)

    assert np.all(calculated_result == true_result)

def test_population_pairwise_correlation_destexhe():
    trains = correlated_trains_destexhe(60,10,5,1000)[0]
    res = population_pairwise_pearson_correlation(trains)
    if DEBUG:
        print res

def main():
    test_population_pairwise_correlation()
    test_population_pairwise_correlation_destexhe()

if __name__ == "__main__":
    if DEBUG:
        cProfile.run("main()")
    else:
        main()
