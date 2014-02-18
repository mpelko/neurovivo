try:
    import neuron
    neuron_loaded = True
except ImportError:
    neuron_loaded = False 

try:
    import brian
    brian_loaded = True
except ImportError:
    brian_loaded = False 

if neuron_loaded:
    from LarkumCell import *
    #from HayCell import *
    #from HayPassiveCell import *
    from SingleCompartmentCell import *
    from SimplifiedL5Cell import *

if brian_loaded:
    from LIFCell import *