'''
Created on March 3, 2018
@author: Andrew Abi-Mansour
'''

# !/usr/bin/python
# -*- coding: utf8 -*-

from PyGran import simulation
from PyGran.modules.params import organic

# Create a dictionary of physical parameters
params = {

	# Define the system
	'boundary': ('f','f','f'), # fixed BCs
	'box':  (-1, 1, -1 , 1, -1, 1), # simulation box size

	# Define component(s)
	'species': ({'material': organic, 'style': 'multisphere', 'nspheres': 12, 'radius': 2e-2, 'length': 1e-1},
				),

	# Set skin distance to be 1/4 particle diameter 
	'nns_skin': 5e-3,

	# Timestep
	'dt': 2e-7,
 
	# Apply gravitional force in the negative direction along the z-axis
	'gravity': (9.81, 0, 0, -1),

	# Setup I/O
	'traj': {'pfile': 'particles.dump', 'mfile': 'tumbler*.vtk'},

	# Stage runs [optional]
         'stages': {'insertion': 1e6, 'rotation': 1e6},

	# Define mesh for rotating drum (trumbler)
	'mesh': {
              'tumbler': {'file': 'mesh/tumbler.stl', 'mtype': 'mesh/surface/stress', 'material': organic, 'args': {'scale': 1e-3}},
	}
  }

def run(**params):

	# Create an instance of the DEM class
	sim = simulation.DEM(**params)

	# Insert 800 particles in a cylindrical region
	insert = sim.insert(species=1, value=1, region=('cylinder', 'y', 0, 0, 0.7, -0.4, 0.4)) #, args={'orientation': 'random'})

	# Add viscous force
	air_resistance = sim.add_viscous(species=1, gamma=0.1)

	# Run insertion stage
	sim.run(params['stages']['insertion'], params['dt'])

	# Delete insertion fix
	sim.remove(insert)

	# Rotate mesh along the xoz plane
	sim.moveMesh(name='tumbler', rotate=('origin', 0, 0, 0), axis=(0, 1, 0), period=5e-1)

	# Run rotation stage
	sim.run(params['stages']['rotation'], params['dt'])

if __name__ == '__main__':
	run(**params)
