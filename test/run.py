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

			# Define component(s)
			'SS': ({'id':1, 'natoms': 1e4, 'density': 2500.0, 'insert': True, 'rate': 1e6, 'freq': 1e3}, ), # number of particles inserted = rate * dt * freq every freq steps
			'radius': (('gaussian number', 4e-4, 4e-5),),

			'traj': {'sel': 'all', 'freq': 10000, 'dir': 'traj', 'style': 'custom', 'file': 'traj.dump', \
                       'args': ('id', 'x', 'y', 'z', 'radius', 'omegax', 'omegay', 'omegaz', \
                       'vx', 'vy', 'vz', 'fx', 'fy', 'fz')},

			# Material properties
			'materials': glass,

			# Apply gravitional force in the negative direction along the z-axis
			# 'gravity': (9.81, 0, 0, -1),

			# Assign initial velocity along the z-direction. TODO: Must assign velocities AFTER particles are created.
			#'velocity': (('all', 'set', 0, -1e-3, 0),),

			# Stage runs
			'stages': {'insertion': 1e6, 'run': 1e7},

			# Meshes
			'mesh': {
				'hopper': {'file': 'hopper-2cm-6cm.stl', 'scale': 5e-4},
			      },
		  }

	# Instantiate a linear spring dashpot class
	pDict['model'] = pDict['model']('thorn', **pDict)

	# Create an instance of the DEM class
	sim = Simulator.DEM(**pDict['model'].params)

	# Setup a stopper wall along the xoz plane (y = 0.0)
	sim.setupWall(name='stopper', wtype='primitive', plane = 'zplane', peq = 0.0)

	sim.insert('void', *('cylinder', 'z', 0, 0, 0.001, 0.001, 0.005))
	sim.run()
	sim.remove(name='stopper')