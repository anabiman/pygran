# !/usr/bin/python
# -*- coding: utf8 -*-

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
			  'nSim': 1,
			  'output': 'out',
			  'SS': ({'id':1, 'natoms': 5000, 'density': 2500, 'insert':True, 'rate':10**6, 'freq':10**3, 'region': \
				('sphere' , 0, 0.25, 0, 0.05)}, ), # rate of particles inserted = rate x dt x freq
			  'nSS': 2, 
			  'vel': ((0,0.0,0),),
			  'radius': (('constant', 0.00224),),
			  'gravity': (9.81, 0, -1, 0), # apply gravitional force in the negative direction along the y-axis
			  'box':  (-0.42, 0.42, -0.01, 0.62, -0.42, 0.42), # simulation box size
			  'print': (10**4, 'time', 'atoms', 'fmax', 'ke'), # print the time, atom number, avg. kinetic energy, and max force
			  'insertionRun':  {'steps':  10**4, 'dt': 10**-5},
			  'productionRun': {'steps': 5 * 10**4, 'dt': 10**-5},

			  'surfMesh': {
					'hopper': {'file': 'hopper.stl', 'scale':0.02},
				      },
			  }

	# Create an instance of the DEM class
	sim = DEM(**params)

	# Define the domain, create atoms, and initialize masses, velocities, etc.
	sim.initialize()

	# Setup material properties
	sim.createProperty('mYmod', *('youngsModulus', 'peratomtype', '5.e6', '5.e6'))
	sim.createProperty('mPratio', *('poissonsRatio', 'peratomtype', '0.45', '0.45'))
	sim.createProperty('mCfric', *('coefficientFriction', 'peratomtypepair', '2', '1.0', '1.0', '1.0', '1.0'))
	sim.createProperty('mCvel', *('characteristicVelocity', 'scalar', '2.0', '2.0'))

	# For 4 simulations, we need four tuples for the coefficient of restitution
	coeffRest = ('coefficientRestitution', 'peratomtypepair', '2', '0.051', '0.051', '0.051', '0.051')
			#('coefficientRestitution', 'peratomtypepair', '2', '1.0', '1.0', '1.0', '1.0'))

	# Overloaded function 'createProperty' will partition coeffRest based on MPI's coloring split scheme
	sim.createProperty('mCrest', *coeffRest)

	# Import all meshes
	for mesh in params['surfMesh'].keys():
		sim.importMesh(name=mesh, **params['surfMesh'][mesh])

	# Setup hopper mesh as a wall
	sim.setupWall(name='hopper', wtype='mesh')

	# Setup a stopper wall along the xoz plane (y = 0.0)
	sim.setupWall(name='stopper', wtype='primitive', plane = 'yplane', peq = 0.0)

	# Print output specified in 'print' every 'freq' steps
	sim.printSetup()

	# Create an NVE (micro canonical) integrator
	sim.setupIntegrate(name='intMicro')

	# Insertion run
	sim.integrate(**params['insertionRun'])

	# Write 'all' coordinates to 'traj.xyz' file every 'freq' steps
	sim.dumpSetup()

	# Remove stopper
	# sim.remove(name='stopper')

	# Monitor rotational KE
    	sim.monitor(name='ke', group='all', var='globKE', file='ke.dat')
	sim.monitor(name='erotate/sphere', group='all', var='globRE', file='re.dat')

	# Production run
	sim.integrate(**params['productionRun'])	

	# Save rotational KE to dat file
	sim.plot(fname='ke.dat', xlabel='Time (s)', ylabel='KE (J/atom)', output='ke.pdf', xscale = params['productionRun']['dt'])
	sim.plot(fname='re.dat', xlabel='Time (s)', ylabel='RE (J/atom)', output='re.pdf', xscale = params['productionRun']['dt'])
