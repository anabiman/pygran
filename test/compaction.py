'''
Created on July 9, 2016
@author: Andrew Abi-Mansour
'''

# !/usr/bin/python
# -*- coding: utf8 -*-

from PyDEM import Simulator, Analyzer, Visualizer
from PyDEM.Materials import glass



realGlass = {
	'youngsModulus': ('youngsModulus', 'peratomtype', '6.3e6'),
	'poissonsRatio': ('poissonsRatio', 'peratomtype', '0.24'),
	'coefficientFriction': ('coefficientFriction', 'peratomtypepair', '0.5'),
	'coefficientRollingFriction': ('coefficientRollingFriction', 'peratomtypepair', '5e-1'),
	'cohesionEnergyDensity': ('cohesionEnergyDensity', 'peratomtypepair', '0.1'),
	'yieldPress': ('yieldPress', 'peratomtype', '6.2e4'),
	}

# Create a dictionary of physical parameters
pDict = {

		'model': Simulator.models.Hysteresis,
        'engine': Simulator.engines.liggghts,

		# Define the system
		'boundary': ('p','p','f'), # fixed BCs
		'box':  (-0.001, 0.001, -0.001, 0.001, -0.0001, 0.004), # simulation box size
		'nns_skin': 1e-3,

		# Define component(s)
		'SS': ({'by_pack':'', 'id':1, 'natoms_local': 800, 'density': 2500.0, 'vol_lim': 1e-16, 'insert': True, 'rate': 1e6, \
			'freq': 'once', 'radius': ('gaussian number', 70e-6, 10e-6)}, 
		      ), 

		'traj': {'sel': 'all', 'freq': 100, 'dir': 'traj', 'style': 'custom', 'file': 'traj.dump', \
                   'args': ('id', 'x', 'y', 'z', 'radius', 'omegax', 'omegay', 'omegaz', \
                   'vx', 'vy', 'vz', 'fx', 'fy', 'fz')},

		'dt': 2e-6,

		# Material properties
		'materials': realGlass,

		# Apply gravitional force in the negative direction along the z-axis
		'gravity': (9.81, 0, 0, -1),

		# Stage runs
		'stages': {'insertion': 1.5e3, 'run': 1.5e3},

		# Meshes
		 'mesh': {
			'hopper': {'file': 'square.stl', 'mtype': 'mesh/surface/stress/servo', 'import': True, \
			'args': {'scale': 1e-3, 'com': '0 0 0', 'ctrlPV': 'force',  'axis': '0 0 -1', 'target_val': 1.76, 'vel_max': 20.0}},
		      },
	  }

if __name__ == '__main__':

	# Instantiate a class based on the selected model
	pDict['model'] = pDict['model'](**pDict)

	# Create an instance of the DEM class
	sim = Simulator.DEM(**pDict['model'].params)

	# Setup a stopper wall along the xoy plane
	sim.setupWall(name='stopper', wtype='primitive', plane = 'zplane', peq = 0.0)

	high = 0.5e-3
	scale = 10.0

	sim.command('dump meshDump all mesh/vtk 100 traj/mesh*.vtk id')

	for i in range(2):

		factorLow = scale * (i / (i+1.))**2.0 * high - i * high
		factorHigh = scale * ((i+1.) / (i+2.))**2.0 * high - i * high

		print 'factorLow = {}, factorHigh = {}'.format(factorLow, factorHigh) 
		insert = sim.insert('void{}'.format(i), 1, *('block', -1e-3, 1e-3, -1e-3, 1e-3, 5e-5 + factorLow, 5e-5 + factorHigh))
		sim.run(pDict['stages']['insertion'], pDict['dt'])
		sim.remove(insert)

	sim.run(pDict['stages']['run'], pDict['dt'])
