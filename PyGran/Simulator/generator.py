from PyGran import Analyzer
from numpy import random, ones, zeros, sqrt
from sys import argv
import collections


def random(natoms, radius, overlapping=True):
	""" Creates a set of randomly placed non-overlapping particles 

	@natoms: number of particles to create
	@radius: a float or array (of size Natoms)
	@[overlapping]: boolean variable which ensures all particles are allowed to overlap or natoms

	Returns a Particles object.
	"""

	scale = 2.0 * natoms**(1.0/3)

	x,y,z = random.rand(natoms) * radius * scale, random.rand(natoms) * radius * scale, random.rand(natoms) * radius * scale

	if type(radius) is float or type(radius) is int:
		r = ones(natoms) * radius
	elif len(r) == natoms:
		r = radius
	else:
		raise ValueError('Input radius can be either a scalar (float or int) or a list/array of length natoms.')

	data = collections.OrderedDict()
	data['x'] = x
	data['y'] = y
	data['z'] = z
	data['radius'] = r
	
	return Analyzer.Particles(**data)

def hcp(natoms, radius):
	""" Generates a Hexagonal Close Packed structure
		
	@natoms: number of particles to create
	@radius: a float or array (of size Natoms)

	Returns a Particles object.
	"""
	data = collections.OrderedDict()
	N = int(natoms**(1/3.0))
	count = 0

	data['x'], data['y'], data['z'] = zeros(natoms), zeros(natoms), zeros(natoms)

	for i in range(N):
		for j in range(N):
			for k in range(N):
				data['x'][count] = (2.0 * i + ((j + k) % 2)) * radius
				data['y'][count] = (sqrt(3.) * (j + 1.0/3.0 * (k % 2))) * radius
				data['z'][count] = (2.0 * sqrt(6.0) / 3.0 * k) * radius
				count += 1

	data['radius'] = ones(natoms) * radius

	return Analyzer.Particles(**data)