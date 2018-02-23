'''
Created on July 9, 2016
@author: Andrew Abi-Mansour
'''

# !/usr/bin/python
# -*- coding: utf8 -*-

from PyGran import Simulator, Analyzer, Visualizer
from PyGran.Materials import stearicAcid

ms = ('nspheres 3', 'ntry 1000000', 'spheres', '0. 0. 0. 4e-4', '8e-4 0 0 4e-4', '16e-4 0 0 4e-4', 'type 1')

# Create a dictionary of physical parameters
pDict = {

		'model': Simulator.models.SpringDashpot,
        	'engine': Simulator.engines.liggghts,

		'style': 'sphere',

		# Define the system
		'boundary': ('p','p','p'), # fixed BCs
		'box':  (-0.001, 0.001, -0.001, 0.001, -0.001, 0.001), # simulation box size
		'nns_skin': 1e-4,

		# Define component(s)
		'SS': ({'material':stearicAcid, 'vol_lim': 1e-16, 'style': 'sphere', 'density': stearicAcid['density'], 'radius': ('constant', 4e-4)},
		       {'material':stearicAcid, 'vol_lim': 1e-16, 'style': 'multisphere', 'density': stearicAcid['density'], 'args': ms, 'natoms':3, 'insert': 'by_pack', 'all_in':'no', 'freq': 'once'},
		      ),

		'dt': 1e-7,

		# Setup I/O
		'traj': {'args': ('id', 'type', 'x', 'y', 'z', 'radius')}, #, 'quat1', 'quat2', 'quat3', 'quat4', 'shapex', 'shapey', 'shapez')},

		# Stage runs
		'stages': {'insertion': 1.5e3, 'run': 1e5},
	  }

if __name__ == '__main__':

	# Instantiate a class based on the selected model
	pDict['model'] = pDict['model'](**pDict)

	# Create an instance of the DEM class
	sim = Simulator.DEM(**pDict['model'].params)

	# Create two spherical particles
	sim.createParticles(1, 'single', *(0,0,-9e-4))
	sim.set(*('atom 1', 'type 1', 'vz 0.1', 'diameter 4e-4'))

	sim.createParticles(1, 'single', *(0,0,9e-4))
        sim.set(*('atom 2', 'type 1', 'vz -0.1', 'diameter 4e-4'))

	# Create MS spheres
	insert = sim.insert('void0', 2, *('block', 2e-4, 2.8e-4, -4e-5, 4e-5, -9e-4, 9e-4, 'units box', 'volume_limit 1e-30'))

	# Create two SQ particles
	# sim.createParticles(1, 'single', *(0,0,-1e-3))
	#sim.set(*('atom 1', 'type 1', 'shape 0.002 0.002 0.004', 'blockiness 20 10', 'vz 0.1', 'density {}'.format(stearicAcid['density']), 'quat/random 81728623'))

	#sim.createParticles(2, 'single', *(0,0,1e-3))
        #sim.set(*('atom 2', 'type 2', 'shape 0.002 0.002 0.004', 'blockiness 20 10', 'vz -0.1', 'density {}'.format(stearicAcid['density']), 'quat/random 81728623'))


	sim.createGroup(*('spheres', 'id', '1 2'))
	sim.createGroup(*('multi', 'id', '3'))

	sim.run(pDict['stages']['run'], pDict['dt'], ['nve/sphere', 'multisphere'], ['spheres', 'multi'])
