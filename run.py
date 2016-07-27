# !/usr/bin/python
# -*- coding: utf8 -*-

from PyDEM import Simulator, Analyzer, Visualizer
from PyDEM.Materials import glass

if __name__ == '__main__':

	# Create a dictionary of physical parameters
	params = {

			'model': Simulator.models.HertzMindlin,
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
			'stages': {'insertion': 1e6},

			# Meshes
			'mesh': {
				'hopper': {'file': 'hopper-2cm-6cm.stl', 'scale': 5e-4},
			      },
		  }

	# Instantiate a linear spring dashpot class
	params['model'] = params['model'](**params)
	cModel = params['model']

	# Create an instance of the DEM class
	sim = Simulator.DEM(**cModel.params)

	# Setup a stopper wall along the xoz plane (y = 0.0)
	sim.setupWall(name='stopper', wtype='primitive', plane = 'zplane', peq = 0.0)

	# Create an NVE (micro canonical) integrator
	sim.setupIntegrate(name='intMicro')

	# Insert particles if not restarting/resuming sim
	cylinder = sim.insertParticles('void', *('cylinder', 'z', 0, 0, 0.001, 0.02, 0.045))
	sim.integrate(cModel.params['stages']['insertion'], cModel.params['dt'])
	sim.remove(cylinder)

	# Remove stopper
	sim.remove(name='stopper')