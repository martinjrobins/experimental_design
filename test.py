from exp_design import map,vector,e_surface
import pylab as plt
from numpy import amax
import numpy as np


params = map()
params['k0'] = 1.0
params['alpha'] = 0.5
params['Cdl'] = 0.0037
params['Ru'] = 0.002
params['Estart'] = -11.7217
params['Ereverse'] = 39.0725
params['E0'] = 0.5*(params['Ereverse'] - params['Estart'])
params['omega'] = 14.4728

Itot = vector()
t = vector()
e_surface(params,Itot,t)

plt.plot(t,Itot)
plt.show()
