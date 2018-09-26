import PyGran.Simulator as Sim
from PyGran.Params import cohesionless
from numpy import arange, fabs, array, sqrt
import matplotlib.pylab as plt
import matplotlib as mpl

# Setup matplotlib params
mpl.rc('text', usetex = True)
plt.rcParams.update({'font.size':18})

cModel = Sim.models.ThorntonNing

COR, yieldVel = [], []

# Create a yield pressure array to study
pressure = arange(1, 6, 0.1) * cohesionless['youngsModulus'] * 0.01

# Set particle radius to 100 microns
cohesionless['radius'] = 1e-4
cohesionless['coefficientRestitution'] = 0.5

for yieldPress in pressure:

	yeff = cohesionless['youngsModulus'] * 0.5 / (1. - cohesionless['poissonsRatio']**2)

	cohesionless['yieldPress'] = yieldPress
	model = cModel(material=cohesionless)

	time, disp, force = model.displacement(dt=2.33322158839e-07)

	deltav = disp[:,1]

	COR.append(fabs(deltav[-1] / cohesionless['characteristicVelocity']))
	yieldVel.append( model.yieldVel / cohesionless['characteristicVelocity'] )

ratio = array(yieldVel)
ratio[ratio > 1] = 1

COR_anal = sqrt(6. * sqrt(3.0) / 5.0) * sqrt(1.0 - (1.0/6.0) * (ratio)**2) * \
    (ratio / (ratio + 2 * sqrt( 6.0/5.0 - (1.0/5.0) * ratio**2 ) ) )**(0.25)


fig = plt.figure(figsize=(16,18), dpi=80)

#print(COR_anal * cohesionless['characteristicVelocity'], fabs(deltav[-4::]), deltav[0])

#plt.plot( deltav / cohesionless['characteristicVelocity']  )
plt.plot(pressure, COR, '-o')
plt.plot(pressure, COR_anal)

plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel('Yield Pressure (MPa)')
plt.ylabel('Coefficient of Restitution')
plt.grid(linestyle=':')
plt.legend(['Numerical', 'Analytical'])
plt.show()

