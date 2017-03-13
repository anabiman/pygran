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

class Neighbors(object):
	""" A dynamic class that contains all the particle-particle (and optionally particle-wall)
	neighbors from which contacts, overlaps, force chains, etc. can be determined.
	"""
	def __init__(self, Particles, cutoff = None):

		coords = numpy.array([Particles.x, Particles.y, Particles.z]).T
		tree = spatial.cKDTree(coords, leafsize=100)
		
		if not cutoff:
			cutoff = 2.0 * Particles.radius.max()

		self._pairs = tree.query_pairs(cutoff)
		self._distances = numpy.zeros(len(self._pairs))
		self._overlaps = numpy.zeros((len(self._pairs),3))

		count = 0

		for pair in self._pairs:
			i, j = pair
			self._distances[count] = numpy.sqrt((Particles.x[i] - Particles.x[j])**2.0 + \
									(Particles.y[i] - Particles.y[j])**2.0 + \
									(Particles.z[i] - Particles.z[j])**2.0)

			if self._distances[count] <= Particles.radius[i] + Particles.radius[j]:
				self._overlaps[count] = (self._distances[count], i, j)

			count += 1 

		self._overlaps = self._overlaps[self._overlaps[:,0] > 0,:]

	@property
	def distances(self):
	    return self._distances

	@property
	def pairs(self):
	    return self._pairs

	@property
	def overlaps(self):
	    return self._overlaps

	def forceChain(self):
		pass
	
	   
	
