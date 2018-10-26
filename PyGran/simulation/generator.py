# !/usr/bin/python
# -*- coding: utf8 -*- 
# -----------------------------------------------------------------------
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 2 of the License, or
#   (at your option) any later version.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

# -------------------------------------------------------------------------

'''
Created on March 30, 2016
@author: Andrew Abi-Mansour
'''

from PyGran import analysis
import numpy
from numpy import random
import collections

def rand(natoms, radius, overlapping=True, factor=1.0):
	""" Creates a set of randomly placed non-overlapping particles 

	@natoms: number of particles to create
	@radius: a float or array (of size Natoms)
	@[overlapping]: boolean variable which ensures all particles are allowed to overlap or natoms
	@[factor]: float to scale noise (w.r.t particle radius) by when correcting for overlapping particles

	Returns a Particles object.
	"""

	scale = 2.0 * natoms**(1.0/3)

	x,y,z = random.rand(natoms) * radius * scale, random.rand(natoms) * radius * scale, random.rand(natoms) * radius * scale

	if type(radius) is float or type(radius) is int:
		r = numpy.ones(natoms) * radius
	elif len(r) == natoms:
		r = radius
	else:
		raise ValueError('Input radius can be either a scalar (float or int) or a list/array of length natoms.')

	data = collections.OrderedDict()
	data['x'] = x
	data['y'] = y
	data['z'] = z
	data['radius'] = r
	
	if not overlapping:

		Particles = analysis.Particles(**data)
		indices_tot = numpy.arange(len(Particles))

		while(True):
			
			Neigh = analysis.Neighbors(Particles)
		
			if len(Neigh.overlaps):
				indices = numpy.array(Neigh.overlaps[:,1:], 'int')
				indices.reshape(len(indices)*2)
				indices = numpy.unique(indices)
 
 				Particles_overlap = Particles[indices]

 				# Find all particle indices not overlapping
 				# Is there a faster way to do this?
 				indices_rem = []
 				for i in indices_tot:
 					if i not in indices:
 						indices_rem.append(i)

 				Particles = Particles[indices_rem]

				Particles_overlap.perturb(factor * Particles_overlap.radius)

				Particles += Particles_overlap
				
				print len(Neigh.overlaps)

			else:
				return Particles

	else:
		return analysis.Particles(**data)

def hcp(natoms, radius):
	""" Generates a Hexagonal Close Packed structure
		
	@natoms: number of particles to create
	@radius: a float or array (of size Natoms)

	Returns a Particles object.
	"""
	data = collections.OrderedDict()
	N = int(natoms**(1/3.0))
	count = 0

	data['x'], data['y'], data['z'] = numpy.zeros(natoms), numpy.zeros(natoms), numpy.zeros(natoms)

	for i in range(N):
		for j in range(N):
			for k in range(N):
				data['x'][count] = (2.0 * i + ((j + k) % 2)) * radius
				data['y'][count] = (numpy.sqrt(3.) * (j + 1.0/3.0 * (k % 2))) * radius
				data['z'][count] = (2.0 * numpy.sqrt(6.0) / 3.0 * k) * radius
				count += 1

	data['radius'] = numpy.ones(natoms) * radius

	return analysis.Particles(**data)