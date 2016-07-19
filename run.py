# !/usr/bin/python
# -*- coding: utf8 -*-

from PyDEM import Simulator, Analyzer, Visualizer
from PyDEM.Materials import glass

if __name__ == '__main__':

	# Create a dictionary of physical parameters
	params = {

			'model': Analyzer.models.HertzMindlin,
            'engine': Simulator.engines.liggghts,

			# Define the system
			'boundary': ('f','f','f'), # fixed BCs
			'box':  (-0.016, 0.016, -0.016, 0.016, -0.01, 0.061), # simulation box size

			# Define component(s)
			'SS': ({'id':1, 'natoms': 3e4, 'density': 2500.0, 'insert': True, 'rate': 1e6, 'freq': 1e3}, ), # number of particles inserted = rate * dt * freq every freq steps
			'radius': (('gaussian number', 4e-4, 4e-5),),

			# Material properties
			'materials': glass,

			# Apply gravitional force in the negative direction along the z-axis
			'gravity': (9.81, 0, 0, -1),

			# Stage runs
			'insertion': 1e6,
			'flow': 0e6,

			# I/O parameters
			'traj': {'sel': 'all', 'freq': 1000, 'dir': 'traj', 'style': 'custom', 'file': 'traj.dump', \
			       'args': ('id', 'x', 'y', 'z', 'radius', 'omegax', 'omegay', 'omegaz', \
			                'vx', 'vy', 'vz', 'fx', 'fy', 'fz')},
			'output': 'out-radius-400-repos',
			'print': (10**4, 'time', 'atoms', 'fmax', 'ke', 'cpu', 'cu', 'density'),

			# Meshes
			'surfMesh': {
				'hopper': {'file': 'hopper-2cm-6cm.stl', 'scale': 5e-4},
			      },

			# grab shared libraries from here
			'path': '/usr/lib64/',
		  }

	# Instantiate a linear spring dashpot class
	params['model'] = params['model'](**params)
	cModel = params['model']

	for ss in cModel.params['SS']:
		print 'Attempting to insert {} particles every {} steps for component {}'.format(ss['rate'] * ss['freq'] * cModel.params['dt'], ss['freq'], ss['id'])

	# Create an instance of the DEM class
	sim = Simulator.DEM(**cModel.params)

	# Define the domain, create atoms, and initialize masses, velocities, etc.
	sim.initialize()

	# Setup material properties
	for item in cModel.params['materials'].keys():
		# Overloaded function 'createProperty' will partition coeffRest based on MPI's coloring split scheme
		sim.createProperty(item, *cModel.params['materials'][item])

	# Import and setup all meshes as rigid wall
	for mesh in cModel.params['surfMesh'].keys():
		sim.importMesh(name=mesh, **cModel.params['surfMesh'][mesh])
		sim.setupWall(name=mesh + 'Wall', wtype='mesh', meshName=mesh)

	# Setup a stopper wall along the xoz plane (y = 0.0)
	sim.setupWall(name='stopper', wtype='primitive', plane = 'zplane', peq = 0.0)

	# Print output specified in 'print' every 'freq' steps
	sim.printSetup()

	# Create an NVE (micro canonical) integrator
	sim.setupIntegrate(name='intMicro')

	# Write 'all' coordinates to 'traj.xyz' file every 'freq' steps
	sim.dumpSetup()

	# Insert particles if not restarting/resuming sim
	cylinder = sim.insertParticles('void', *('cylinder', 'z', 0, 0, 0.001, 0.02, 0.045))
	sim.integrate(cModel.params['insertion'], cModel.params['dt'])
	sim.remove(cylinder)

	# Remove stopper
	sim.remove(name='stopper')

	sim.integrate(cModel.params['flow'], cModel.params['dt'])