'''
Place for various mathematical functions I am often using.

Created on Jan 13, 2012

@author: mpelko
'''

import numpy as np

def expexp_fun(x,tau1,tau2):
    """
    Returns the normalised result exp-exp function (normalised so that the max of 
    the function is 1). Used in modelling synaptic currents following a 
    presynaptic spike.
    
    Calculations in ws/philipp/calculations, explanations in books on modelling 
    synapses.
    """
    if tau1==tau2:
        tau1*=1.
        factor = np.e/tau1
        x = np.array(x)
        result=factor*x*(np.exp(-x/tau1))
    else:
        tau1, tau2 = 1.*tau1, 1.*tau2
        factor = 1/((tau2/tau1)**(tau2/(tau1-tau2))-(tau2/tau1)**(tau1/(tau1-tau2)))
        x = np.array(x)
        result=factor*(np.exp(-x/tau1)-np.exp(-x/tau2))
    return result