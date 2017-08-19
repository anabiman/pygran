'''
Created on July 9, 2016
@author: Andrew Abi-Mansour
'''

# !/usr/bin/python
# -*- coding: utf8 -*-

from PyGran import Simulator, Analyzer, Visualizer
from PyGran.Materials import glass

# Create a dictionary of physical parameters
pDict = {

		'model': Simulator.models.SpringDashpot,
        	'engine': Simulator.engines.liggghts,

		# Define the system
		'boundary': ('p','p','f'), # fixed BCs
		'box':  (-0.001, 0.001, -0.001, 0.001, -0.0001, 0.004), # simulation box size
		'nns_skin': 1e-3,

		# Define component(s)
		'SS': ({'insert':'by_pack', 'id':1, 'natoms': 800, 'material':glass, 'vol_lim': 1e-16, \
			'freq': 'once', 'radius': ('gaussian number', 70e-6, 10e-6)}, 
		      ), 

		'traj': {'sel': 'all', 'freq': 100, 'dir': 'traj', 'style': 'custom', 'file': 'traj.dump', \
                   'args': ('id', 'x', 'y', 'z', 'radius', 'omegax', 'omegay', 'omegaz', \
                   'vx', 'vy', 'vz', 'fx', 'fy', 'fz')},

		'dt': 2e-6,

		# Apply gravitional force in the negative direction along the z-axis
		'gravity': (9.81, 0, 0, -1),

		# Stage runs
		'stages': {'insertion': 1.5e3, 'run': 1.5e3},

		# Meshes
		 'mesh': {
			'hopper': {'file': 'square.stl', 'mtype': 'mesh/surface/stress', 'import': True, 'material': glass, \
			'args': ('scale', '1e-3')},
		      },
	  }

if __name__ == '__main__':

	# Instantiate a class based on the selected model
	pDict['model'] = pDict['model'](**pDict)

	# Create an instance of the DEM class
	sim = Simulator.DEM(**pDict['model'].params)

	# Setup a stopper wall along the xoy plane
	sim.setupWalls(name='stopper', wtype='primitive', plane = 'zplane', peq = 0.0)

	high = 0.5e-3
	scale = 10.0

	for i in range(2):

		factorLow = scale * (i / (i+1.))**2.0 * high - i * high
		factorHigh = scale * ((i+1.) / (i+2.))**2.0 * high - i * high

		print 'factorLow = {}, factorHigh = {}'.format(factorLow, factorHigh) 
		insert = sim.insert('void{}'.format(i), 1, *('block', -1e-3, 1e-3, -1e-3, 1e-3, 5e-5 + factorLow, 5e-5 + factorHigh))
		sim.run(pDict['stages']['insertion'], pDict['dt'])
		sim.remove(insert)

	sim.run(pDict['stages']['run'], pDict['dt'])
