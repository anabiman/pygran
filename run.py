
from Liggghts import DEM

if __name__ == '__main__':

	# Create a dictionary of simulation parameters
	params = {
			  'units': 'si',
			  'dim': 3,
			  'style': 'granular', # spherical deformable particles
			  'boundary': ('s','s','s'), # shrink-wrapping BCs
			  'model': ('gran', 'model', 'hooke'),
			  'restart': (10**4, 'restart', 'restart.*', False),
			  'traj': ('all', 500, 'traj', 'traj.xyz'),
			  'nSim': 4,
			  'output': 'out',
			  'dt': 10**-5,
			  'nSS': 1,  # number of components / subsystems
			  'idSS': [1],
			  'vel': [(0.0,-0.01,0)],
			  'insertRate': [10**7],
			  'insertFreq': [10**3], # frequency of inserting particles
			  'radius': [('constant', 0.00224)],
			  'gravity': (9.81, 0, -1, 0), # apply gravitional force in the negative direction along the y-axis
			  'box': (-0.84, 0.84, -0.01, 1.2, -0.84, 0.84), # simulation box size
			  'Natoms': [1000],
			  'print': ('time', 'atoms', 'ke', 'fmax'), # print the time, atom number, avg. kinetic energy, and max force
			  'density': [2500.0], # kg /m^3
			  'insertionRun': 10**4,
			  'productionRun': 2 * 10**4,
			  'freq': 10**4, # print output every freq steps
			  'mesh': 'hopper.stl',
			  'scaleMesh': 0.04
			  }

	# Create an instance of the DEM class
	sim = DEM(**params)

	# Define the domain, create atoms, and initialize masses, velocities, etc.
	sim.initialize()

	# Setup material properties
	sim.createProperty('mYmod', *('youngsModulus', 'peratomtype', '5.e6', '5.e6'))
	sim.createProperty('mPratio', *('poissonsRatio', 'peratomtype', '0.45', '0.45'))
	sim.createProperty('mCfric', *('coefficientFriction', 'peratomtypepair', '2', '0.05', '0.05', '0.05', '0.05'))
	sim.createProperty('mCvel', *('characteristicVelocity', 'scalar', '2.0', '2.0'))

	# Modify material property ~ useful for sensitivity analysis / series of simulations
	coeffRest = [('coefficientRestitution', 'peratomtypepair', '2', '0.9', '0.9', '0.9', '0.9'),
		('coefficientRestitution', 'peratomtypepair', '2', '1.1', '1.1', '1.1', '1.1'),
		('coefficientRestitution', 'peratomtypepair', '2', '1.3', '1.3', '1.3', '1.3'),
		('coefficientRestitution', 'peratomtypepair', '2', '1.5', '1.5', '1.5', '1.5')]

	sim.createProperty(name='mCrest', values=coeffRest)

	# Import mesh hopper.stl file
	sim.importMesh(name='hopper')

	# Setup mesh as a wall
	sim.setupWall(name='hopper', wtype='mesh')

	# Setup a stopper wall along the xoz plane (y = 0.1)
	sim.setupWall(name='stopper', wtype='primitive', plane = 'yplane', peq = 0.1)

	# Print output specified in 'print' every 'freq' steps
	sim.printSetup(freq=params['freq'])

	# Setup integration method for insertion
	sim.setupIntegrate(name='intInsert')

	# Insertion stage
	sim.integrate(steps=params['insertionRun'])

	# Write 'all' coordinates to 'traj.xyz' file every 'freq' steps
	sim.dumpSetup()

	# Remove stopper
	sim.remove(name='stopper')

	# Setup integration method for production run
	sim.setupIntegrate(name='intProd', dt=10**-5)

	# Production run
	sim.integrate(steps=params['productionRun'])	
