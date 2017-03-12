'''
Created on July 9, 2016
@author: Andrew Abi-Mansour
'''

# !/usr/bin/python
# -*- coding: utf8 -*-

from PyGran import Simulator, Analyzer, Visualizer
from PyGran.Materials import glass

if __name__ == '__main__':

	# Create a dictionary of physical parameters
	pDict = {

			'model': Simulator.models.Hysteresis,
            		'engine': Simulator.engines.liggghts,

			# Define the system
			'boundary': ('f','f','f'), # fixed BCs
			'box':  (-0.016, 0.016, -0.016, 0.016, -0.01, 0.061), # simulation box size

			# Define component(s)
			'SS': ({'id':1, 'natoms': 1, 'density': 2500.0, 'insert': True, 'rate': 1e6, \
				'freq': 100, 'radius': ('gaussian number', 20e-4, 4e-5)}, 
				{'id':2, 'natoms': 1, 'density': 2500.0, 'insert': True, 'rate': 1e6, \
				'freq': 100, 'radius': ('gaussian number', 20e-4, 4e-5)}
			      ), 
				# number of particles inserted = rate * dt * freq every freq steps

			'traj': {'sel': 'all', 'freq': 1, 'dir': 'traj', 'style': 'custom', 'file': 'traj.dump', \
                       'args': ('id', 'x', 'y', 'z', 'radius', 'omegax', 'omegay', 'omegaz', \
                       'vx', 'vy', 'vz', 'fx', 'fy', 'fz')},

			'dt': 3e-05,

			# Material properties
			'materials': glass,

			# Apply gravitional force in the negative direction along the z-axis
			'gravity': (9.81, 0, 0, 1),

			# Stage runs
			'stages': {'insertion': 50, 'run': 0e6},

			# Meshes
			'mesh': {
				'hopper': {'file': 'hopper-2cm-6cm.stl', 'scale': 5e-4},
			      },
		  }

	# Instantiate a class based on a user-defined model
	pDict['model'] = pDict['model'](**pDict)

	# Create an instance of the DEM class
	sim = Simulator.DEM(**pDict['model'].params)

	# Setup a stopper wall along the xoz plane (y = 0.0)
	sim.setupWall(name='stopper', wtype='primitive', plane = 'zplane', peq = 0.0)

	sim.insert('void', 1, *('cylinder', 'z', 0, 0, 0.01, 0.01, 0.04001))
	sim.insert('void2', 2, *('cylinder', 'z', 0, 0, 0.01, 0.0, 0.00001))

	sim.run(pDict['stages']['insertion'], pDict['dt'])
