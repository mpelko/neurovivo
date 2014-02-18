"""
Conclusion:
Both O-U process generations seem to work properly. The var of the process in the long run is 
sigma**2/(2*theta)
"""


import numpy as np
import neurovivo.common as cmn
import PyProcess.pyprocess as pyp
import matplotlib.pylab as pyl

theta = 0.1
mu = 1
sigma = 2
IC = mu
dt = 0.1

times = np.arange(10000) * dt

OU_parameters= {
                "theta":theta,
                "mu":mu,
                "sigma":sigma, 
                }
time_space_constraints = {'startTime': 0, 'startPosition':IC}

OU1 = pyp.OU_process(OU_parameters, time_space_constraints)

pairs1 = []
pairs2 = []
for _ in xrange(100):
    OU1_sample = np.array(OU1.generate_sample_path(times)).T[1]
    OU2_sample = cmn.OU_process(theta, mu, sigma, times, IC)
    pairs1.append([np.mean(OU1_sample),np.std(OU1_sample)])
    pairs2.append([np.mean(OU2_sample),np.std(OU2_sample)])

pairs1=np.array(pairs1).T
pairs2=np.array(pairs2).T

pyl.plot(pairs1[0],pairs1[1],"bo")
pyl.plot(pairs2[0],pairs2[1],"go")
pyl.plot([np.mean(pairs1[0])],[np.mean(pairs1[1])],"ro")
pyl.plot([np.mean(pairs2[0])],[np.mean(pairs2[1])],"ro")

print np.mean(pairs1[0]),np.mean(pairs1[1])
print np.mean(pairs2[0]),np.mean(pairs2[1])
pyl.show()