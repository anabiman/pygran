'''
Created on July 9, 2016
@author: Andrew Abi-Mansour
'''

# !/usr/bin/python
# -*- coding: utf8 -*-

from PyDEM import Simulator, Analyzer, Visualizer
from PyDEM.Materials import glass

if __name__ == '__main__':

	# Create a dictionary of physical parameters
	pDict = {

			'model': Simulator.models.Hysteresis,
            'engine': Simulator.engines.liggghts,

			# Define the system
			'boundary': ('f','f','f'), # fixed BCs
			'box':  (-0.016, 0.016, -0.016, 0.016, -0.01, 0.061), # simulation box size

			# 'model-args': ('gran', 'model', 'hooke', 'limitForce', 'on', 'ktToKnUser', 'on'),

			# Define component(s)
			'SS': ({'id':1, 'natoms': 10000, 'density': 2500.0, 'insert': True, 'rate': 1e6, \
				'freq': 100, 'radius': ('gaussian number', 4e-4, 4e-5)}, 
			      ), 
				# number of particles inserted = rate * dt * freq every freq steps

			'traj': {'sel': 'all', 'freq': 1000, 'dir': 'traj', 'style': 'custom', 'file': 'traj.dump', \
                       'args': ('id', 'x', 'y', 'z', 'radius', 'omegax', 'omegay', 'omegaz', \
                       'vx', 'vy', 'vz', 'fx', 'fy', 'fz')},

			'dt': 1e-05,

			# Material properties
			'materials': glass,

			# Apply gravitional force in the negative direction along the z-axis
			'gravity': (9.81, 0, 0, -1),

			# Stage runs
			'stages': {'insertion': 2e4, 'run': 0e6},

			# Meshes
			'mesh': {
				'hopper': {'file': 'hopper-2cm-6cm.stl', 'scale': 5e-4},
			      },
		  }

	# Instantiate a class based on the selected model
	pDict['model'] = pDict['model'](**pDict)

	# Create an instance of the DEM class
	sim = Simulator.DEM(**pDict['model'].params)

	# Setup a stopper wall along the xoz plane (y = 0.0)
	sim.setupWall(name='stopper', wtype='primitive', plane = 'zplane', peq = 0.0)

	sim.insert('void', 1, *('cylinder', 'z', 0, 0, 0.01, 0.01, 0.04001))

	sim.run(pDict['stages']['insertion'], pDict['dt'])
