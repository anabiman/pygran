# !/usr/bin/python
# -*- coding: utf8 -*-

from Liggghts import DEM

if __name__ == '__main__':

	# Create a dictionary of simulation parameters
	params = {

			  # Define the system
			  'units': 'si',
			  'dim': 3,
			  'style': 'granular', # spherical deformable particles
			  'boundary': ('f','f','f'), # fixed BCs
			  'model': ('gran', 'model', 'hooke', 'tangential', 'history'),
			  'box':  (-0.032, 0.032, -0.032, 0.032, -0.01, 0.122), # simulation box size

			  # Define component(s)
			  'SS': ({'id':1, 'natoms': 3000, 'density': 2500, 'insert': True, 'rate':10**6, 'freq': 10**3}, ), # rate of particles inserted = rate x dt x freq
			  'nSS': 2, 
			  'vel': ((0,0.0,0), ),
			  'radius': (('constant', 2.24e-3),),

			  # Apply gravitional force in the negative direction along the y-axis
			  'gravity': (9.81, 0, 0, -1), 

			  # I/O parameters
			  'restart': (10**4, 'restart', 'restart.*', False),
			  'traj': {'sel':'all', 'freq':100, 'dir':'traj', 'file':'traj.custom', 'args': ('x', 'y', 'z', 'omegax', 'omegay', 'omegaz', 'radius')},
			  'nSim': 1,
			  'output': 'out-3000',
			  'print': (10**4, 'time', 'atoms', 'fmax', 'ke', 'cpu', 'cu'), # print the time, atom number, avg. kinetic energy, and max force

			  # Stage runs
			  'insertion':  {'steps':  2 * 10**4, 'dt': 10**-4},
			  'flow': {'steps': 10**5, 'dt': 10**-4},

			  # Meshes
			  'surfMesh': {
					'hopper': {'file': 'hopper.stl', 'scale':0.001},
				      },

			  # Nearest neighbor searching paramers
			  'nns_skin': 1e-3,
			  'nns_type': 'bin',

			  }

	# Create an instance of the DEM class
	sim = DEM(**params)

	# Define the domain, create atoms, and initialize masses, velocities, etc.
	sim.initialize()

	# Setup material properties
	sim.createProperty('mYmod', *('youngsModulus', 'peratomtype', '7.1e7', '7.1e7'))
	sim.createProperty('mPratio', *('poissonsRatio', 'peratomtype', '0.22', '0.22'))
	sim.createProperty('mCfric', *('coefficientFriction', 'peratomtypepair', '2', '0.1', '0.1', '0.1', '0.1'))
	sim.createProperty('mCvel', *('characteristicVelocity', 'scalar', '0.01', '0.01'))

	# For 'm' simulations, we need 'm' tuples for the coefficient of restitution
	coeffRest = (('coefficientRestitution', 'peratomtypepair', '2', '0.9', '0.9', '0.9', '0.9'))

	# Overloaded function 'createProperty' will partition coeffRest based on MPI's coloring split scheme
	sim.createProperty('mCrest', *coeffRest)

	# Import all meshes
	for mesh in params['surfMesh'].keys():
		sim.importMesh(name=mesh, **params['surfMesh'][mesh])

	# Setup hopper mesh as a wall
	sim.setupWall(name='hopper', wtype='mesh')

	# Setup a stopper wall along the xoz plane (y = 0.0)
	sim.setupWall(name='stopper', wtype='primitive', plane = 'zplane', peq = 0.0)

	# Print output specified in 'print' every 'freq' steps
	sim.printSetup()

	# Create an NVE (micro canonical) integrator
	sim.setupIntegrate(name='intMicro')

	# Write 'all' coordinates to 'traj.xyz' file every 'freq' steps
	sim.dumpSetup()

	# Insert particles for stage 1
	lowerSphere = sim.insertParticles('lowerSphere', *('sphere', 0, 0, 0.03, 0.02))
	sim.integrate(**params['insertion'])
	sim.remove(lowerSphere)

	# Insert particles for stage 2
	midSphere = sim.insertParticles('midSphere', *('sphere', 0, 0, 0.06, 0.02))
	sim.integrate(**params['insertion'])
	sim.remove(midSphere)

	# Insert particles for stage 3
	upperSphere = sim.insertParticles('upperSphere', *('sphere', 0, 0, 0.09, 0.02))
	sim.integrate(**params['insertion'])
	sim.remove(upperSphere)

	# Remove stopper
	sim.remove(name='stopper')

	sim.integrate(**params['flow'])

	# Monitor translational and rotational <KE>
   	sim.monitor(name='ke', group='all', var='globKE', file='ke.dat')
	sim.monitor(name='erotate/sphere', group='all', var='globRE', file='re.dat')

	# Save rotational KE to dat file
	sim.plot(fname='ke.dat', xlabel='Time (s)', ylabel='KE (J/atom)', output='ke.pdf', xscale = params['flow']['dt'])
	sim.plot(fname='re.dat', xlabel='Time (s)', ylabel='RE (J/atom)', output='re.pdf', xscale = params['flow']['dt'])
