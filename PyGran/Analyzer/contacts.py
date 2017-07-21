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
from scipy.linalg import norm
from PyGran.Simulator.models import SpringDashpot
import matplotlib.pylab as plt

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
		self._cutoff = cutoff

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

	def forceChain(self, axis = (0,2), alpha = 45):
		""" Computes the force chain based on an algorithm published in Phys. Rev. E. 72, 041307 (2005):
		'Characterization of force chains in granular material'.

		@ axis: a tuple of size 2 that specifies the two axis to use for computing the force chain, e.g. axis=(0,1) -> (x,y)
		@ alpha: the angle (in degrees) that controls the deviation of the force chain. A value of 0 means a perfectly linear chain.
		Thus, alpha is a measure of the 'curvature' of the force chain. See page 5 of the paper cited above.
		"""

		stress = numpy.zeros((self._Particles.natoms, 2, 2))
		stress_prin = numpy.zeros((self._Particles.natoms, 3)) # 2 stress components + angle = 4 dims

		coords = numpy.array([self._Particles.x, self._Particles.y, self._Particles.z])
		x,y = coords[axis[0]], coords[axis[1]]

		# Compute net stresses on all particles
		for contact in self._overlaps:

			i, j = int(contact[1]), int(contact[2])
			overlap = self._Particles.radius[i] + self._Particles.radius[j] - numpy.array([x[i] - x[j], y[i] - y[j]])

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

			stress_p = 0.5 * (stress_i[0, 0] + stress_i[1, 1]) - numpy.sqrt( (0.5 * (stress_i[0, 0] - \
				stress_i[1, 1]))**2.0 + stress_i[0, 1] )

			if stress_i[0,0] - stress_i[1,1]: # otherwise this is prolly a boundary particle

				# Compute angle
				stress_prin[i,2] = 0.5 * numpy.arctan( 2.0 * stress_i[0,1] / (stress_i[0,0] \
								  - stress_i[1,1]) ) * 360.0
				# compute principal stress
				stress_prin[i,0] = stress_p * numpy.cos(stress_prin[i,2])
				stress_prin[i,1] = stress_p * numpy.sin(stress_prin[i,2])

		stress_norm = norm(stress_prin[:,:-1], axis=1)

		# Section V. THE ALGORITHM (page 4)
		# Step 1: filter out particles with less than the mean principal stress
		indices  = numpy.where(stress_norm >= stress_norm.mean(axis=0))[0]
		#stress_prin = stress_prin[indices,:]

		# Step 2: filter out particles which are in contact with 1 or less 'highly stressed' particles
		# Construct an nns list
		coords = self._coords[indices,:][:,axis] 
		tree = spatial.cKDTree(coords)

		nns = tree.query_ball_point(coords, self._Particles[indices].radius.max() * 2.0)
		self._indices = indices
		self._nns = [[] for i in range(len(nns))]

		# Update the nns indices so that only overlapping neighbors remain
		count = 0
		saxis = numpy.array(axis)

		for ns in nns:
			pind = indices[count]
			for index in ns:
				if index != pind:

					dist = numpy.sqrt(((self._coords[index,saxis] - self._coords[pind,saxis])**2.0).sum())

					if dist <= (self._Particles.radius[index] + self._Particles.radius[pind]):
						self._nns[count].append(index)

			count += 1

		count = 0
		cosAlpha = numpy.cos(alpha)
		chain = [[] for i in range(len(self._nns))]

		# Step 3: compute the force chain per particle
		for nns in self._nns:
			if len(nns) > 1: # select only particles surrounded by 2 or more highly stressed particles
				# Begin loop over each neighboring particle
				for ni in nns:
					# check for angle compliance
					angle_i = stress_prin[indices[count],-1]
					angle_j =  stress_prin[ni,-1]

					dist = self._coords[indices[count], saxis] - self._coords[ni,saxis]

					term1 = norm(dist * stress_prin[indices[count],:-1])
					term2 = norm(dist) * norm(stress_prin[indices[count],:-1])

					term3 = norm(- dist * stress_prin[ni,:-1])
					term4 = norm(dist) * norm(stress_prin[ni,:-1])

					# Search in the forward direction
					if (term1 <= term2) and (cosAlpha * term2) < term1:
						if (term3 <= term4) and (cosAlpha * term4) < term3:
							chain[count].append(ni)

					# Search in the reverse direction
					if (term1 <= term2) and (-cosAlpha * term2) > term1: # TODO: double check the math on this one
						if (term3 <= term4) and (-cosAlpha * term4) > term3: # TODO: double check the math on this one
							chain[count].append(ni)

			count += 1

		parts = self._Particles[indices]
		plt.scatter(parts.x *  1e6, parts.z * 1e6, s=parts.radius * 1e6, facecolors='none')

		count = 0
		for ind in chain:
			for i in ind:
				plt.plot([stress_prin[count,0], stress_prin[count,1]] ,'r')

		plt.show()

		#import networkx as nx
		#import numpy as np
		#import matplotlib.pyplot as plt

		#G = nx.Graph()
		#G.add_edges_from([(1,2), (1,3)])

		#pos = {1: (40, 20), 2: (20, 30), 3: (40, 30)}

		#nx.draw_networkx_nodes(G, pos, cmap=plt.get_cmap('jet'))
		#nx.draw_networkx_edges(G, pos, edgelist=((1,2),(1,3)), width=(1.0,5.0), edge_color='black', arrows=True)

		#plt.show()


		return chain

	def findWithin(self, coords , r):
		""" Find all points within distance r of point(s) coords.
		TODO: Support walls aligned arbitrarily in space """

		indices =  numpy.arange(self._Particles.natoms) #self._tree.query_ball_point(coords, r)
		#indices = [item for i in indices for item in i]
		indices = list(numpy.unique(indices))

		if len(indices):
			# calculate distance along the z-axis (must take other axes into account)
			lengths = numpy.fabs(self._coords[indices,-1] - coords[:,-1][0])

			parts = self._Particles[indices]
			return parts[numpy.where(lengths <= parts.radius)]

		return None
