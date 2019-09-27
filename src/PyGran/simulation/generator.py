'''
A set of routines for generating particle configurations

Created on March 30, 2016

Author: Andrew Abi-Mansour

This is the::

  ██████╗ ██╗   ██╗ ██████╗ ██████╗  █████╗ ███╗   ██╗
  ██╔══██╗╚██╗ ██╔╝██╔════╝ ██╔══██╗██╔══██╗████╗  ██║
  ██████╔╝ ╚████╔╝ ██║  ███╗██████╔╝███████║██╔██╗ ██║
  ██╔═══╝   ╚██╔╝  ██║   ██║██╔══██╗██╔══██║██║╚██╗██║
  ██║        ██║   ╚██████╔╝██║  ██║██║  ██║██║ ╚████║
  ╚═╝        ╚═╝    ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═══╝

DEM simulation and analysis toolkit
http://www.pygran.org, support@pygran.org

Core developer and main author:
Andrew Abi-Mansour, andrew.abi.mansour@pygran.org

PyGran is open-source, distributed under the terms of the GNU Public
License, version 2 or later. It is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
received a copy of the GNU General Public License along with PyGran.
If not, see http://www.gnu.org/licenses . See also top-level README
and LICENSE files.
'''

from PyGran import analysis
import numpy
from numpy import random
import collections

def rand(natoms, radius, overlapping=True, factor=1.0):
	""" Creates a set of randomly placed non-overlapping particles 

	:param natoms: number of particles to create
	:type natoms: int
	:param radius: particle radius
	:type radius: float or array of length natoms
	:param overlapping:  ensures all particles are allowed to overlap or not
	:type overlapping: bool
	:param factor: number to scale noise with (w.r.t particle radius) when correcting for overlapping particles
	:type factor: float
	:returns: a new Particles object
	:rtype: Particles
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
				
			else:
				return Particles

	else:
		return analysis.Particles(**data)

def hcp(natoms, radius):
	""" Generates a Hexagonal Close Packed structure
		
	:param natoms: number of particles to create
	:type natoms: int
	:param radius: particle radius
	:type radius: float or array of length natoms
	:returns: a new Particles object
	:rtype: Particles
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