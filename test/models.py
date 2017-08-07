import PyGran.Simulator as Si
from PyGran.Materials import stearicAcid
from numpy import random
import matplotlib.pylab as plt

models = [Si.models.SpringDashpot, Si.models.HertzMindlin]
radius = 1e-4

# Compute the displacement vs time curve for e = 0.9
for model in models:
	model = model(material=stearicAcid)

	time, delta = model.displacement(radius)
	force = model.normalForce(delta[:,0], radius) + model.dissForce(delta[:,0], delta[:,1], radius)

	deltan = delta[force > 0,0]
	force = force[force > 0]
	plt.plot(deltan, force)

plt.legend(['SpringDashpot', 'HertzMindlin'])
plt.show()
