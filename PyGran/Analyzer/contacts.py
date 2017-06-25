'''
Created on March 10, 2017
@author: Andrew Abi-Mansour
'''

# !/usr/bin/python
# -*- coding: utf8 -*- 
# -------------------------------------------------------------------------
#
#   Python module for determining contacts/overlaps in an N-particle system
#
# --------------------------------------------------------------------------
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

# -------------------------------------------------------------------------

import numpy
from scipy import spatial
from PyGran.Simulator.models import SpringDashpot

class Neighbors(object):
	""" A dynamic class that contains all the particle-particle (and optionally particle-wall)
	neighbors from which contacts, overlaps, force chains, etc. can be determined.
	"""
	def __init__(self, Particles, material = None, cutoff = None):

		self._Particles = Particles
		self._coords = numpy.array([Particles.x, Particles.y, Particles.z]).T
		self._tree = spatial.cKDTree(self._coords)
		
		if not cutoff:
			cutoff = 2.0 * Particles.radius.max()

		self._pairs = self._tree.query_pairs(cutoff)
		self._distances = numpy.zeros(len(self._pairs))
		self._overlaps = numpy.zeros((len(self._pairs),3))

		count = 0

		for pair in self._pairs:
			i, j = pair
			self._distances[count] = numpy.sqrt((Particles.x[i] - Particles.x[j])**2.0 + \
									(Particles.y[i] - Particles.y[j])**2.0 + \
									(Particles.z[i] - Particles.z[j])**2.0)

			if self._distances[count] <= Particles.radius[i] + Particles.radius[j]:
				self._overlaps[count] = (Particles.radius[i] + Particles.radius[j] - self._distances[count], i, j)

			count += 1 

		self._overlaps = self._overlaps[self._overlaps[:,0] > 0,:]

		if material:
			self._material = {}
			self._material['SS'] = ({'material':material},)
			self._model = SpringDashpot(**self._material)

	@property
	def distances(self):
	    return self._distances

	@property
	def pairs(self):
	    return self._pairs

	@property
	def overlaps(self):
	    return self._overlaps

	def forceChain(self, axis = (0,2)):
		""" Computes the force chain based on an algorithm published in Phys. Rev. E. 72, 041307 (2005): Characterization of force chains in granular material
		"""
		stress = numpy.zeros((self._Particles.natoms, 3, 3))
		stress_prin = numpy.zeros((self._Particles.natoms, 4)) # 3 stress components + angle = 4 dims

		for contact in self._overlaps:

			i, j = int(contact[1]), int(contact[2])
			overlap = self._Particles.radius[i] + self._Particles.radius[j] - numpy.array([self._Particles.x[i] - self._Particles.x[j], \
				self._Particles.y[i] - self._Particles.y[j], self._Particles.z[i] - self._Particles.z[j]])
			
			vi = 4.0/3.0 * numpy.pi * self._Particles.radius[i]**3.0
			vj = 4.0/3.0 * numpy.pi * self._Particles.radius[j]**3.0

			stress_i = numpy.outer(self._model.normalForce(overlap, self._Particles.radius[i]), overlap)
			stress_j = numpy.outer(self._model.normalForce(overlap, self._Particles.radius[j]), overlap)

			stress[i] += stress_i
			stress[j] += stress_j

		# Faster to loop over all particles
		for i in range(self._Particles.natoms):

			# Compute principal stress
			stress_i = stress[i]
			stress_p = 0.5 * (stress_i[axis[0], axis[0]] + stress_i[axis[1], axis[1]]) - numpy.sqrt( (0.5 * (stress_i[axis[0], axis[0]] - \
				stress_i[axis[1], axis[1]]))**2.0 + stress_i[axis[0],axis[1]] )

			stress_prin[i,:-1] = stress_p
			stress_prin[i,-1] = 0.5 * numpy.arctan( 2.0 * stress_i[axis[0],axis[1]] / (stress_i[axis[0], axis[0]] - stress_i[axis[1], axis[1]]) )

		return stress_prin

	def findWithin(self, coords , r):
		""" Find all points within distance r of point(s) coords. 
		TOD: Support walls aligned arbitrarily in space """

		indices =  numpy.arange(self._Particles.natoms) #self._tree.query_ball_point(coords, r)
		#indices = [item for i in indices for item in i]
		indices = list(numpy.unique(indices))
		
		if len(indices):
			# calculate distance along the z-axis (must take other axes into account)
			lengths = numpy.fabs(self._coords[indices,-1] - coords[:,-1][0])

			parts = self._Particles[indices]
			return parts[numpy.where(lengths <= parts.radius)]
	
		return None