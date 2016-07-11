# !/usr/bin/python
# -*- coding: utf8 -*-

'''
Created on May 6, 2016
@author: Andrew Abi-Mansour

Center for Materials Sci. & Eng.,
Merck Inc., West Point
'''

from DEM import DEMSi

if __name__ == '__main__':

	# Create a dictionary of simulation parameters
	params = {

			  # Use LIGGGHTS as the DEM computational engine
			  'modName': 'liggghts',

			  # Define the system
			  'units': 'si',
			  'dim': 3,
			  'style': 'granular', # spherical deformable particles
			  'boundary': ('f','f','f'), # fixed BCs
			  'model': ('gran', 'model', 'hooke', 'tangential', 'history', 'rolling_friction', 'cdt', 'tangential_damping', 'on', 'limitForce', 'on'), # the order matters here
			  'box':  (-0.016, 0.016, -0.016, 0.016, -0.01, 0.061), # simulation box size

			  # Define component(s)
			  'SS': ({'id':1, 'natoms': 2.5e5, 'density': 2500, 'insert': True, 'rate':10**6, 'freq': 10**3}, ), # rate of particles inserted = rate x dt x freq
			  'nSS': 2, 
			  'vel': ((0,0,0), ),
			  'radius': (('gaussian number', 4.0e-4, 4.0e-5),),

			  # Apply gravitional force in the negative direction along the z-axis
			  'gravity': (9.81, 0, 0, -1), 

			  # I/O parameters
			  'restart': (5000, 'restart', 'restart.*', False),
			  'traj': {'sel': 'all', 'freq': 1000, 'dir': 'traj', 'file': 'traj.custom', 'args': ('id', 'x', 'y', 'z', 'radius', 'omegax', 'omegay', 'omegaz', 'vx', 'vy', 'vz', 'fx', 'fy', 'fz')},
			  'dump_modify': ('append', 'yes'),
			  'nSim': 1,
			  'output': 'out-radius-400-gaussian-roll-0.9-0.9',
			  'print': (10**4, 'time', 'atoms', 'fmax', 'ke', 'cpu', 'cu', 'density'), # print the time, atom number, avg. kinetic energy, and max force

			  # Stage runs
			  'insertion':  {'steps':  5e4, 'dt': 1e-5},
			  'flow': {'steps': 3e5, 'dt': 2e-5},

			  # Meshes
			  'surfMesh': {
					'hopper': {'file': 'hopper-2cm-6cm.stl', 'scale':5e-4},
				      },

			  # Nearest neighbor searching paramers
			  'nns_skin': 1e-3,
			  'nns_type': 'bin',

			  # grab shared libraries from here
			  'path': '/usr/lib64/',
			  }

	# Create an instance of the DEM class
	sim = DEMSi.DEM(**params)

	# Define the domain, create atoms, and initialize masses, velocities, etc.
	sim.initialize()

	# Setup material properties
	sim.createProperty('mYmod', *('youngsModulus', 'peratomtype', '7.1e7', '7.1e7'))
	sim.createProperty('mPratio', *('poissonsRatio', 'peratomtype', '0.22', '0.22'))
	sim.createProperty('mCfric', *('coefficientFriction', 'peratomtypepair', '2', '0.1', '0.9', '0.9', '0.9'))
	sim.createProperty('mRfric', *('coefficientRollingFriction', 'peratomtypepair', '2', '0.9', '0.9', '0.9', '0.9'))
	sim.createProperty('mCvel', *('characteristicVelocity', 'scalar', '0.01', '0.01'))

	# For 'm' simulations, we need 'm' tuples for the coefficient of restitution
	coeffRest = (('coefficientRestitution', 'peratomtypepair', '2', '0.9', '0.9', '0.9', '0.9'))

	# Overloaded function 'createProperty' will partition coeffRest based on MPI's coloring split scheme
	sim.createProperty('mCrest', *coeffRest)

	# Import and setup all meshes as rigid wall
	for mesh in params['surfMesh'].keys():
		sim.importMesh(name=mesh, **params['surfMesh'][mesh])
		sim.setupWall(name=mesh + 'Wall', wtype='mesh', meshName=mesh)

	# Setup a stopper wall along the xoz plane (y = 0.0)
	sim.setupWall(name='stopper', wtype='primitive', plane = 'zplane', peq = 0.0)

	# Print output specified in 'print' every 'freq' steps
	sim.printSetup()

	# Create an NVE (micro canonical) integrator
	sim.setupIntegrate(name='intMicro')

	# Write 'all' coordinates to 'traj.xyz' file every 'freq' steps
	sim.dumpSetup()

	if not params['restart'][3]:
		# Insert particles for stage 1
		cylinder = sim.insertParticles('void', *('cylinder', 'z', 0, 0, 0.014, 0.005, 0.055))
		sim.integrate(**params['insertion'])
		sim.remove(cylinder)

	# Remove stopper
	sim.remove(name='stopper')

	sim.integrate(**params['flow'])
