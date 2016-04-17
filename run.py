
import matplotlib.pylab as plt
from Liggghts import DEM
import numpy as np

if __name__ == '__main__':

	# Create a dictionary of simulation parameters
	params = {
			  'units': 'si',
			  'dim': 3,
			  'style': 'granular', # spherical deformable particles
			  'boundary': ('f','f','f'), # reflective BCs
			  'model': ('gran', 'model', 'hooke'),
			  'restart': (10**4, 'restart/restart.*'),
			  'dt': 10**-6,
			  'nSS': 1,  # number of components / subsystems
			  'idSS': [1],
			  'vel': [(0.001,-0.001,001)],
			  'insertRate': [50000],
			  'insertFreq': [1000],
			  'radius': [('constant', 0.00224)],
			  'gravity': (0, 0, -1, 0), # apply gravitional force in the negative direction along the y-axis
			  'box': (-2.1, 2.1, -0.1, 3.1, -2.1, 2.1), # simulation box size
			  'Natoms': [10000],
			  'print': ('time', 'atoms', 'ke', 'fmax'), # print the time, atom number, and kinetic energy
			  'density': [2500.0], # kg /m^3
			  'freq': 10**4, # frequency of saving/printing output
			  'insertionRun': 10**6,
			  'productionRun': 10**5, 
			  'mesh': 'hopper.stl',
			  'scaleMesh': 0.1
			  }

	# Create an instance of the Reaction class
	sim = DEM(**params)

	# Define the domain, create atoms, and initialize masses, velocities, etc.
	sim.initialize()

	# Setup material properties
	sim.lmp.command('fix m1 all property/global youngsModulus peratomtype 5.e6 5.e6')
	sim.lmp.command('fix m2 all property/global poissonsRatio peratomtype 0.45 0.45')
	sim.lmp.command('fix m3 all property/global coefficientRestitution peratomtypepair 2 0.9 0.9 0.9 0.9')
	sim.lmp.command('fix m4 all property/global coefficientFriction peratomtypepair 2 0.05 0.05 0.05 0.05')
	sim.lmp.command('fix m5 all property/global characteristicVelocity scalar 2.0 2.0')

	# Import mesh
	sim.importMesh(var='hopper')

	# Setup mesh as a wall
	sim.setupWall(var='hopper', wtype='mesh')

	# Setup a piston wall
	sim.setupWall(var='stopper', wtype='primitive', plane = 'yplane', peq = 0.0)

	# Print output specified in 'print' every 'freq' steps
	sim.printSetup(freq=params['freq'])

	# Setup integration method
	sim.setupIntegrate(name='integrator')

	# Monitor temperature as a function of time
	sim.monitor(name='globKE', group='all', var='ke')

	# Insertion stage
	sim.integrate(steps=params['insertionRun'])

	# Plot temperature vs time, then save the figure as a pdf
	plt.rc('text', usetex=True)
	time = np.array(range(len(sim.vars))) * params['freq'] * params['dt']
	#plt.plot(time, sim.vars)
	#plt.xlabel(r"Time (s)")
	#plt.ylabel("KE (Joules)")
	#plt.savefig("KE.pdf")

	# Write 'all' coordinates to 'traj.xyz' file every 'freq' steps
	sim.dumpSetup(sel='all', freq=params['freq'], traj='traj.xyz')

	# Equilibration stage
	sim.integrate(steps=params['productionRun'])

