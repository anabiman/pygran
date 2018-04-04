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
			'box':  (-0.0011, 0.0011, -0.0011, 0.0011, -0.0001, 0.004), # simulation box size
			'nns_skin': 0.0001,

			# Define component(s)
			'SS': ({'by_pack':'', 'id':1, 'natoms': 40000, 'density': 2500.0, 'vol_lim': 1e-16, 'insert': True, 'rate': 1e6, \
				'freq': 'once', 'radius': ('gaussian number', 10.e-6, 1.e-6)}, 
			      ), 

			'traj': {'sel': 'all', 'freq': 1, 'dir': 'traj', 'style': 'custom', 'file': 'traj.dump', \
                       'args': ('id', 'x', 'y', 'z', 'radius', 'omegax', 'omegay', 'omegaz', \
                       'vx', 'vy', 'vz', 'fx', 'fy', 'fz')},

			'dt': 1e-7,

			# Material properties
			'materials': glass,

			# Apply gravitional force in the negative direction along the z-axis
			# 'gravity': (9.81, 0, 0, -1),

			# Stage runs
			'stages': {'insertion': 1, 'run': 0},

			# Meshes
			'mesh': {
				'hopper': {'file': 'cylinder.stl', 'scale': 1e-3},
			      },
		  }

	# Instantiate a class based on the selected model
	pDict['model'] = pDict['model'](**pDict)

	# Create an instance of the DEM class
	sim = Simulator.DEM(**pDict['model'].params)

	# Setup a stopper wall along the xoy plane
	# sim.setupWall(name='stopper', wtype='primitive', plane = 'zplane', peq = 0.0)

	sim.insert('void', 1, *('cylinder', 'z', 0, 0, 0.0009, 0.0, 0.002))
	sim.run(pDict['stages']['insertion'], pDict['dt'])

	#sim.remove('stopper')
	sim.run(pDict['stages']['run'], pDict['dt'])

