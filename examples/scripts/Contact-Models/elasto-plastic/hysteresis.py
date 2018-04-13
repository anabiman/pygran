import PyGran.Simulator as Sim
from PyGran.Materials import organic
from numpy import arange, fabs
import matplotlib.pylab as plt
import matplotlib as mpl

# Setup matplotlib params
mpl.rc('text', usetex = True)
plt.rcParams.update({'font.size':18})

cModel = Sim.models.ThorntonNing

COR = []

# Create a yield pressure array to study
pressure = arange(1, 6, 0.1) * organic['youngsModulus'] * 0.01

# Set particle radius to 100 microns
organic['radius'] = 1e-4

for yieldPress in pressure:

	organic['yieldPress'] = yieldPress
	model = cModel(material=organic)

	time, disp, force = model.displacement()
	deltav = disp[:,1]

	COR.append(fabs(deltav[-1] / deltav[0]))

fig = plt.figure(figsize=(16,18), dpi=80)
plt.plot(pressure / 1e6, COR, 'o-')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel('Yield Pressure (MPa)')
plt.ylabel('Coefficient of Restitution')
plt.grid(linestyle=':')
plt.show()

