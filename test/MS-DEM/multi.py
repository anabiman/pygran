'''
Created on July 9, 2016
@author: Andrew Abi-Mansour
'''

# !/usr/bin/python
# -*- coding: utf8 -*-

from PyGran import Simulator, Analyzer, Visualizer
from PyGran.Materials import stearicAcid

stearicAcid['coefficientRestitution'] = 0.051
stearicAcid['density'] = 1.0

ms1 = ('nspheres 3', 'ntry 1000000', 'spheres', '0 0 0 1.5e-4', '2e-4 0 0 1.5e-4', '4e-4 0 0 1.5e-4', 'type 1')
ms2 = ('nspheres 6', 'ntry 1000000', 'spheres', '0 0 0 1.5e-4', '0.8e-4 0 0 1.5e-4', '1.6e-4 0 0 1.5e-4', '2.4e-4 0 0 1.5e-4', '3.2e-4 0 0 1.5e-4', '4e-4 0 0 1.5e-4', 'type 1')
ms3 = ('nspheres 12', 'ntry 1000000', 'spheres', '0 0 0 1.5e-4', '0.363e-4 0 0 1.5e-4', '0.727e-4 0 0 1.5e-4', '1.091e-4 0 0 1.5e-4', '1.454e-4 0 0 1.5e-4', '1.818e-4 0 0 1.5e-4',  \
'2.182e-4 0 0 1.5e-4', '2.545e-4 0 0 1.5e-4', '2.909e-4 0 0 1.5e-4', '3.273e-4 0 0 1.5e-4', '3.636e-4 0 0 1.5e-4', '4e-4 0 0 1.5e-4', 'type 1')

# Create a dictionary of physical parameters
pDict = {

		'model': Simulator.models.SpringDashpot,
        	'engine': Simulator.engines.liggghts,

		# Define the system
		'boundary': ('p','p','p'), # fixed BCs
		'box':  (-0.002, 0.004, -0.001, 0.001, -0.001, 0.001), # simulation box size
		'nns_skin': 1e-3,

		# Define component(s)
		'SS': (
		#	{'material':stearicAcid, 'vol_lim': 1e-16, 'style': 'sphere', 'density': stearicAcid['density'], 'radius':('constant', 1.5e-4)},
		       {'material':stearicAcid, 'vol_lim': 1e-16, 'style': 'multisphere', 'density': stearicAcid['density'], 'args': ms3},
		      ),

		'dt': 1e-8,

		# Setup I/O
		'traj': {'args': ('id', 'type', 'x', 'y', 'z', 'radius')}, #, 'quat1', 'quat2', 'quat3', 'quat4', 'shapex', 'shapey', 'shapez')},

		# Stage runs
		'stages': {'insertion': 1.5e3, 'run': 15e5},
	  }

if __name__ == '__main__':

	# Instantiate a class based on the selected model
	pDict['model'] = pDict['model'](**pDict)

	# Create an instance of the DEM class
	sim = Simulator.DEM(**pDict['model'].params)

	# Create two spherical particles
	# sim.createParticles(1, 'single', *(0,0,-7e-4))
	# sim.set(*('atom 1', 'type 1', 'vz 0.1', 'diameter 4e-4'))

	# sim.createParticles(1, 'single', *(0,0,7e-4))
        # sim.set(*('atom 2', 'type 1', 'vz -0.1', 'diameter 4e-4'))

	# Create spheres
	#insert = sim.insert(name='void0', species=1, value=1, vel=(0,0,0.1), all_in= 'no', region=('block', -1e-6, 1e-6, -1e-6, 1e-6, -6.5e-4, -6.4e-4, 'units box', 'volume_limit 1e-30'))
	#insert = sim.insert(name='void1', species=1, value=1, vel=(0,0,-0.1), all_in= 'no', region=('block', -1e-6, 1e-6, -1e-6, 1e-6, 6.4e-4, 6.5e-4, 'units box', 'volume_limit 1e-30'))

	# Create multispheres
	insert = sim.insert(name='void2', species=1, value=1, vel=(0,0,0.1), all_in= 'no', region=('block', 2e-3, 2.001e-3, -1e-7, 1e-7, -1e-3, -6e-4, 'units box', 'volume_limit 1e-30'))
	insert = sim.insert(name='void3', species=1, value=1, vel=(0,0,-0.1), all_in= 'no', region=('block', 2e-3, 2.001e-3, -1e-7, 1e-7,  6e-4, 1e-3, 'units box', 'volume_limit 1e-30'))

	#insert = sim.insert(name='void4', species=3, value=1, vel=(0,0,0.1), region=('block', 2e-3, 2.01e-3, -1e-5, 1e-5, -1e-3, -6e-4, 'units box', 'volume_limit 1e-30'))
        #insert = sim.insert(name='void5', species=3, value=1, vel=(0,0,-0.1), region=('block', 2e-3, 2.01e-3, -1e-5, 1e-5, 6e-4, 1e-3, 'units box', 'volume_limit 1e-30'))

	# Create two SQ particles
	# sim.createParticles(1, 'single', *(0,0,-1e-3))
	#sim.set(*('atom 1', 'type 1', 'shape 0.002 0.002 0.004', 'blockiness 20 10', 'vz 0.1', 'density {}'.format(stearicAcid['density']), 'quat/random 81728623'))

	#sim.createParticles(2, 'single', *(0,0,1e-3))
        #sim.set(*('atom 2', 'type 2', 'shape 0.002 0.002 0.004', 'blockiness 20 10', 'vz -0.1', 'density {}'.format(stearicAcid['density']), 'quat/random 81728623'))

	sim.run(pDict['stages']['run'], pDict['dt'])
