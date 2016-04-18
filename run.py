
from Liggghts import DEM
import os

if __name__ == '__main__':

	# Create a dictionary of simulation parameters
	params = {
			  'units': 'si',
			  'dim': 3,
			  'style': 'granular', # spherical deformable particles
			  'boundary': ('p','p','p'), # periodic BCs
			  'model': ('gran', 'model', 'hooke'),
			  'restart': (10**4, 'restart', 'restart.*', False),
			  'traj': ('all', 5 * 10**2, 'traj', 'traj.xyz'),
			  'dt': 0.8 * 10**-6,
			  'nSS': 1,  # number of components / subsystems
			  'idSS': [1],
			  'vel': [(0.0,-0.01,0)],
			  'insertRate': [10**6],
			  'insertFreq': [10**3], # frequency of inserting particles
			  'radius': [('constant', 0.00224)],
			  'gravity': (9.81, 0, -1, 0), # apply gravitional force in the negative direction along the y-axis
			  'box': (-0.21, 0.21, -0.01, 0.31, -0.21, 0.21), # simulation box size
			  'Natoms': [30000],
			  'print': ('time', 'atoms', 'ke', 'fmax'), # print the time, atom number, avg. kinetic energy, and max force
			  'density': [2500.0], # kg /m^3
			  'insertionRun': 5 * 10**4,
			  'productionRun': 0.2 * 10**5,
			  'freq': 10**4, # print output every freq steps
			  'mesh': 'hopper.stl',
			  'scaleMesh': 0.01
			  }

	# Create an instance of the DEM class
	sim = DEM(**params)

	# Define the domain, create atoms, and initialize masses, velocities, etc.
	sim.initialize()

	# Setup material properties
	sim.lmp.command('fix m1 all property/global youngsModulus peratomtype 5.e6 5.e6')
	sim.lmp.command('fix m2 all property/global poissonsRatio peratomtype 0.45 0.45')
	sim.lmp.command('fix m3 all property/global coefficientRestitution peratomtypepair 2 0.9 0.9 0.9 0.9')
	sim.lmp.command('fix m4 all property/global coefficientFriction peratomtypepair 2 0.05 0.05 0.05 0.05')
	sim.lmp.command('fix m5 all property/global characteristicVelocity scalar 2.0 2.0')

	# Import mesh hopper.stl file
	sim.importMesh(var='hopper')

	# Setup mesh as a wall
	sim.setupWall(var='hopper', wtype='mesh')

	# Setup a stopper wall along the xoz plane (y = 0.1)
	sim.setupWall(var='stopper', wtype='primitive', plane = 'yplane', peq = 0.1)

	# Print output specified in 'print' every 'freq' steps
	sim.printSetup(freq=params['freq'])

	# Setup integration method for insertion
	sim.setupIntegrate(name='intInsert')

	# Insertion stage
	sim.integrate(steps=params['insertionRun'])

	# Monitor KE as a function of time
	sim.monitor(name='globKE', group='all', var='ke')

	# Write 'all' coordinates to 'traj.xyz' file every 'freq' steps
	sim.dumpSetup()

	# Setup integration method for equilibration
	sim.setupIntegrate(name='intProd', dt=10**-5)

	# Equilibration stage
	sim.integrate(steps=params['productionRun'])

	# Plot KE vs time, then save the figure as a pdf
	sim.plot(name='globKE', xlabel='Time (s)', ylabel='KE (J)')