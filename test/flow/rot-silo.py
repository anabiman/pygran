'''
Created on April 22, 2017
@author: Andrew Abi-Mansour
'''

# !/usr/bin/python
# -*- coding: utf8 -*-

from PyGran import Simulator
from PyGran.Materials import glass, stearicAcid

pDict = {

		# Setup model + DEM engine
		'model': Simulator.models.SpringDashpot,
		'engine': Simulator.engines.liggghts,

		# Define the system
		'boundary': ('f','f','f'),
		'box':  (-1e-3, 1e-3, -1e-3, 1e-3, 0, 4e-3),

		# Define component(s)
		'SS': ({'insert': 'by_pack', 'material': stearicAcid, 'natoms': 1000, 'freq': 'once', 'radius': ('gaussian number', 5e-5, 5e-6), 'vol_lim': 1e-16}, 
		      ),

		# Setup I/O params
		'traj': {'freq':1000, 'pfile': 'traj.dump', 'mfile': 'mesh*.vtk'},

		# Define computational parameters
		'nns_skin': 1e-3,
		'dt': 1e-6,

		# Apply a gravitional force in the negative direction along the z-axis
		'gravity': (9.81, 0, 0, -1),

		# Import hopper mesh
		 'mesh': {
			'hopper': {'file': 'silo.stl', 'mtype': 'mesh/surface', 'import': True, 'material': glass, 'args': ('scale 1e-3',)},
			'valve': {'file': 'valve.stl', 'mtype': 'mesh/surface', 'import': True, 'material': glass, 'args': ('move 0 0 1.0', 'scale 1e-3',)},
		      },

		# Stage runs
		'stages': {'insertion': 1e4, 'run': 5e4},
	  }

# Instantiate a class based on the selected model
pDict['model'] = pDict['model'](**pDict)

# Create an instance of the DEM class
sim = Simulator.DEM(**pDict['model'].params)

# Setup a stopper wall along the xoy plane
sim.setupWalls(name='stopper', wtype='primitive', plane = 'zplane', peq = 0.0)

sim.moveMesh('valve', 'rotate origin 0. 0. 0.', 'axis  0. 0. 1.', 'period 5e-2')

# Insert particles in a cubic region
insert = sim.insert('cubic', 1, *('block', -5e-4, 5e-4, -5e-4, 5e-4, 2e-3, 3e-3))
sim.run(pDict['stages']['insertion'], pDict['dt'])
sim.remove(insert)

# Run equilibration run
sim.run(pDict['stages']['run'], pDict['dt'])

# Remove stopper and monitor flow
sim.remove('stopper')
sim.run(pDict['stages']['run'], pDict['dt'])
