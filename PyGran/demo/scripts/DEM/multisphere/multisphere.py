'''
Created on March 3, 2018
@author: Andrew Abi-Mansour
'''

# !/usr/bin/python
# -*- coding: utf8 -*-

from PyGran import simulation
from PyGran.params import organic

def template_tablet(nspheres, radius, length):
	""" This function creates a multi-sphere tablet (cylinder) of
	radius "radius" and height "length" constituting nspheres spheres.
	"""

	delta = (2 * radius * nspheres - length) / (nspheres-1)
	ms = ('nspheres {}'.format(nspheres), 'ntry 1000000 spheres')

	for i in range(nspheres):
		ms = ms + ((2 * radius - delta ) * i, 0, 0, radius,)

	ms = ms + ('type 1',)

	return ms


# Create a dictionary of physical parameters
params = {

	# Define the system
	'boundary': ('f','f','f'), # fixed BCs
	'box':  (-1, 1, -1 , 1, -1, 1), # simulation box size

	# Define component(s)
	'species': ({'material': organic, 'style':'multisphere', 'args': template_tablet(12, 2e-2, 1e-1)},
				),

	# Set skin distance to be 1/4 particle diameter 
	'nns_skin': 5e-3,

	# Timestep
	'dt': 2e-7,
 
	# Apply gravitional force in the negative direction along the z-axis
	'gravity': (9.81, 0, 0, -1),

	# Setup I/O
	'traj': {'pfile': 'rotating.dump', 'mfile': 'particles*.vtk'},

	# Stage runs [optional]
         'stages': {'insertion': 1e6, 'rotation': 1e6},

	# Define mesh for rotating drum (trumbler)
	'mesh': {
              'tumbler': {'file': 'mesh/tumbler.stl', 'mtype': 'mesh/surface/stress', 'material': organic, 'args': {'scale': 1e-3}},
	}
  }

# Create an instance of the DEM class
sim = simulation.DEM(**params)

# Insert 800 particles in a cylindrical region
insert = sim.insert(species=1, value=800, region=('cylinder', 'y', 0, 0, 0.7, -0.4, 0.4), args={'orientation': 'random'})

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
