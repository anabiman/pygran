import PyGran.Simulator as Sim
from numpy import arange, fabs

cModel = Sim.models.ThorntonNing

# Define material properties
powderX = {
	'youngsModulus': 1e8,
	'poissonsRatio': 0.25,
	'coefficientRestitution': 0.9,
	'characteristicVelocity': 0.1,
	'density': 997.164,
	'radius': 1e-4
}

# Initialize variables
COR = []
pressure = arange(1e6, 4e6, 1e5)

for yieldPress in pressure:

	powderX['yieldPress'] = yieldPress
	model = cModel(material=powderX)

	time, disp, force = model.displacement()
	deltav = disp[:,1]

	COR.append(fabs(deltav[-1] / deltav[0]))
