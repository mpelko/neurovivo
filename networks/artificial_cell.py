import numpy as np
from neuron import h

class ArtificialCell(object):
    vecStim = None

    def __init__(self, fire_times=h.Vector()):
        assert type(fire_times).__name__ == "HocObject",\
            "fire_times must be a hoc Vector"
        self.fire_times = np.array(fire_times)
        self.vecStim = h.VecStim()
        if not fire_times.size() == 0:
            firetimes = h.Vector(self.fire_times)
            self.vecStim.play(firetimes)
        else:
            self.vecStim.play()

    def setFireTimes(self, fire_times):
        assert type(fire_times).__name__ == "HocObject",\
            "fire_times must be a hoc Vector"
        if not fire_times.size() == 0:
            self.vecStim.play(fire_times)
            self.fire_times = np.array(fire_times)
        else:
            self.vecStim.play()
