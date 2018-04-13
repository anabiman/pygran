import PyGran.Simulator as Si
from PyGran.Materials import stearicAcid
import matplotlib.pylab as plt

# Create a list of contact models to analyze
models = [Si.models.SpringDashpot, Si.models.HertzMindlin]

# Set particle radius to 100 microns
stearicAcid['radius'] = 1e-4

# Compute the force-displacement curve
for model in models:
	model = model(material=stearicAcid)

	time, delta, force = model.displacement()

	deltan = delta[force > 0,0]
	force = force[force > 0]
	plt.plot(deltan * 1e6, force * 1e3)

plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(['SpringDashpot', 'HertzMindlin'], loc='best')
plt.xlabel(r'$\delta$ $(\mu m)$')
plt.ylabel('Force (mN)')
plt.grid(linestyle=':')
plt.show()
