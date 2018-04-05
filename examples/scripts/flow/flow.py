'''
Created on April 22, 2017
@author: Andrew Abi-Mansour
'''

# !/usr/bin/python
# -*- coding: utf8 -*-

from PyGran import Simulator, Visualizer
from PyGran.Materials import organic, glass

params = {
	# Define the system
	'boundary': ('f','f','f'),
	'box':  (-1e-3, 1e-3, -1e-3, 1e-3, 0, 4e-3),

	# Define component(s)
	'SS': ({'material': organic, 'radius': ('constant', 5e-5)},),

	# Setup I/O params
	'traj': {'freq':1000, 'style': 'custom/vtk', 'pfile': 'traj*.vtk', 'mfile': 'mesh*.vtk'},

	# Define computational parameters
	'dt': 1e-6,

	# Apply a gravitional force in the negative direction along the z-axis
	'gravity': (9.81, 0, 0, -1),

	# Import hopper mesh
	 'mesh': {
		'hopper': {'file': 'mesh/silo.stl', 'mtype': 'mesh/surface', 'material': glass, 'args': ('scale 1e-3',)},
		'impeller': {'file': 'mesh/valve.stl', 'mtype': 'mesh/surface', 'material': glass, 'args': ('move 0 0 1.0', 'scale 1e-3',)},
	},

	# Stage runs
	'stages': {'insertion': 1e5, 'run': 1e5},
}

# Create an instance of the DEM class
sim = Simulator.DEM(**params)

# Setup a stopper wall along the xoy plane
stopper = sim.setupWall(species=1, wtype='primitive', plane = 'zplane', peq = 0.0)

# Insert particles in a cubic region
insert = sim.insert(species=1, region=('block', -5e-4, 5e-4, -5e-4, 5e-4, 2e-3, 3e-3), mech='volumefraction_region', value=1, freq=1e4)
sim.run(params['stages']['insertion'], params['dt'])
sim.remove(insert)

# Rotate the impeller
rotImp = sim.moveMesh('impeller', 'rotate origin 0. 0. 0.', 'axis  0. 0. 1.', 'period 5e-2')

# Blend the system
sim.run(params['stages']['run'], params['dt'])
sim.remove(rotImp)

# Remove stopper and monitor flow
sim.remove(stopper)
sim.run(params['stages']['run'], params['dt'])
