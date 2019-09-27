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
#   the Free Software Foundation, either version 2 of the License, or
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
from ..simulation.models import SpringDashpot

class Neighbors(object):
	""" A dynamic class that contains all the particle-particle (and optionally particle-wall)
	neighbors from which contacts, overlaps, force chains, etc. can be determined.

	@Particles = a derivative of SubSystem to be analyzed
	@[material]: Python dictionary for material params
	@[cutoff]: max radius by default
	@[binary]: False by default. Set to True when analyzing 2-component systems.

	"""
	def __init__(self, Particles, material = None, cutoff = None, binary=False):

		self._Particles = Particles
		self._coords = numpy.array([Particles.x, Particles.y, Particles.z]).T
		self._tree = spatial.cKDTree(self._coords)

		if not cutoff:
			cutoff = 2.0 * Particles.radius.max()

		self._binary = binary
		if binary:
			if not hasattr(self._Particles, 'type'):
				raise RuntimeError('Binary Particles object must contain type attribute!')

		self._neigh = self._tree.query_ball_point(self._coords, cutoff)

		self._pairs = self._tree.query_pairs(cutoff)
		self._distances = numpy.zeros(len(self._pairs))
		self._overlaps = numpy.zeros((len(self._pairs),3), dtype=object)  # Can we use mixed dtypes for improved mem allocation?
		self._cutoff = cutoff

		count = 0

		for pair in self._pairs:
			i, j = pair
			self._distances[count] = numpy.sqrt((Particles.x[i] - Particles.x[j])**2.0 + \
									(Particles.y[i] - Particles.y[j])**2.0 + \
									(Particles.z[i] - Particles.z[j])**2.0)

			if self._distances[count] <= Particles.radius[i] + Particles.radius[j]:
				self._overlaps[count] = numpy.array([Particles.radius[i] + Particles.radius[j] - self._distances[count], i, j], dtype=object)

			count += 1

		self._overlaps = numpy.array(self._overlaps, dtype=object) # Can we use mixed dtypes for improved mem allocation?
		self._overlaps = self._overlaps[self._overlaps[:,0] > 0,:]

		if material:
			self._material = {}
			self._material['SS'] = ({'material':material},)
			self._model = SpringDashpot(**self._material)

	@property
	def distances(self):
		" Returns all neighbor pair-wise distance "
		return self._distances

	@property
	def pairs(self):
		" Returns the pair indices for each overlap "
		return self._pairs

	def coon(self, type1=1, type2=2):
		""" Computes the coordination number per particle. For binary mixtures
		the default coordinations numbers are returned as two cross coordination 
		number arrays (tuple). For auto (self) coordination numbers, specify type1
		and type2 accordingly.

		[type1]: which type to compute 1st coon for in a binay mixture
		[type2]:  which type to compute 2nd coon for in a binay mixture

		TODO: support multi-body entities for single component systems
		TODO: support tertiary systems 

		Returns an array of coon for all particles of size natoms x ntypes
		"""
		
		if self._binary:
			
			partsA = self._Particles[self._Particles.type == type1]
			partsB = self._Particles[self._Particles.type == type2]

			Na, Nb = len(partsA), len(partsB)

			typeA, typeB = numpy.zeros(Na, dtype='int64'), numpy.zeros(Nb, dtype='int64')
			count = 0

			# Must take into account the case where both A and B are molecules!!!

			# Begin calculating the C.N. of type A
			for i, cn in enumerate(self._neigh):
				if self._Particles[i].type == type1:
					if hasattr(self._Particles, 'mol'):
						if partsB.mol.max() > 0:
							# B is a molecule but A is a sphere
							typeA[count] = len(numpy.unique((self._Particles[cn][self._Particles[cn].type == type2]).mol))
						else:
							# A and B are spheres / *not* molecules
							if type1 == type2:
								typeA[count] = sum(self._Particles[cn].type == type2) - 1
							else:
								typeA[count] = sum(self._Particles[cn].type == type2)
					else:
						# no molecules , just spheres!
						if type1 == type2:
							typeA[count] = sum(self._Particles[cn].type == type2) - 1
						else:
							typeA[count] = sum(self._Particles[cn].type == type2)

					count += 1

			count = 0
			for i, cn in enumerate(self._neigh):
				if self._Particles[i].type == type2:
					if hasattr(self._Particles, 'mol'):
						if partsA.mol.max() > 0:
							typeB[cpunt] = len(numpy.unique((self._Particles[cn][self._Particles[cn].type == type1]).mol))
						else:
							# A and B are spheres / *not* molecules
							if type1 == type2:
								typeB[count] = sum(self._Particles[cn].type == type1)  - 1
							else:
								typeB[count] = sum(self._Particles[cn].type == type1)
					else:
						if type1 == type2:
							typeB[count] = sum(self._Particles[cn].type == type1) - 1
						else:
							typeB[count] = sum(self._Particles[cn].type == type1)

					count += 1

			# See if a molecule is defined (i.e. multi-body entity)
			if hasattr(self._Particles, 'mol'):
				if partsA.mol.max() > 0:
					tmpA = numpy.zeros(int(partsA.mol.max()), dtype='int64')

					for i, partA in enumerate(partsA):
						tmpA[int(partA.mol) - 1] += typeA[i] # mols are numbered from 1 ... N in LAMMPS

					typeA = tmpA
				
				if partsB.mol.max() > 0:
					tmpB = numpy.zeros(int(partsB.mol.max()), dtype='int64')

					for i, partB in enumerate(partsB):
						tmpB[int(partB.mol) - 1] += typeB[i] # mols are numbered from 1 ... N in LAMMPS

					typeB = tmpB

			return numpy.array([typeA, typeB]).T
		else:

			if hasattr(self._Particles, 'mol'):
				if self._Particles.mol[0] > 0:
					coon = numpy.zeros(len(self._Particles.mol))
					for i, cn in enumerate(self._neigh):
						coon[int(self._Particles.mol[i])] += len(numpy.unique(self._Particles[cn].mol)) - 1

					return coon
				else:
					return numpy.array([len(cn) - 1 for cn in self._neigh])
			else:
				return numpy.array([len(cn) - 1 for cn in self._neigh])

	@property
	def overlaps(self):
		""" Returns all overlapping distances and their indices in the form of an N x 3 numpy array
		where N is the number of particles. 1st column contains the distances, 2nd and 3rd colum the corresponding
		indices """
		return self._overlaps

	def filter(self, percent):
		""" Returns a Particles object in which particles overlap no more than `percent' of their effective radius

		@[percent]: filter all particles overlapping by a certain percentage.

		Returns a new Particles class
		"""

		percent /= 100.0
		indices = []

		Particles = self._Particles
		indices = list(range(len(Particles)))

		if percent:
			for count, pair in enumerate(self._pairs):
				i, j = pair
				reff = Particles.radius[i] * Particles.radius[j] / (Particles.radius[i] + Particles.radius[j])

				if (Particles.radius[i] + Particles.radius[j] - self._distances[count] > percent *  reff):

					# this must be slow ..... gotta optimize this for large systems with high overlaps
					if i in indices:
						indices.remove(i)

					if j in indices:
						indices.remove(j)

		indices = numpy.unique(indices)
		indices = numpy.array(indices, dtype='int64')

		return Particles[indices]

	def forceChain(self, axes = (0,2), alpha = numpy.pi/4, plot_stress=True, plot_parts=False, peters=True, threshold=1):
		""" Computes the force chain based on an algorithm published in Phys. Rev. E. 72, 041307 (2005):
		'Characterization of force chains in granular material'.

		@ [axes]: a tuple of size 2 that specifies the two axes to use for computing the force chain, e.g. axes=(0,1) -> (x,y)
		@ [alpha]: the angle (in radians) that controls the deviation of the force chain. A value of 0 means a perfectly linear chain.
		Thus, alpha is a measure of the 'curvature' of the force chain. See page 5 of the paper cited above.
		"""

		stress = numpy.zeros((self._Particles.natoms, 2, 2))
		stress_prin = numpy.zeros((self._Particles.natoms, 3)) # 2 stress components + angle = 4 dims

		coords = numpy.array([self._Particles.x, self._Particles.y, self._Particles.z])
		x,y = coords[axes[0]], coords[axes[1]]

		# Compute net stresses on all particles
		for contact in self._overlaps:

			i, j = int(contact[1]), int(contact[2])
			overlap = self._Particles.radius[i] + self._Particles.radius[j] + numpy.array([x[i] - x[j], y[i] - y[j]])

			vi = 4.0/3.0 * numpy.pi * self._Particles.radius[i]**3.0
			vj = 4.0/3.0 * numpy.pi * self._Particles.radius[j]**3.0

			stress_i = numpy.outer(self._model.normalForce(overlap, self._Particles.radius[i]), overlap) / vi
			stress_j = numpy.outer(self._model.normalForce(overlap, self._Particles.radius[j]), overlap) / vj

			stress[i] += stress_i
			stress[j] += stress_j

		# Faster to loop over all particles
		for i in range(self._Particles.natoms):

			# Compute principal stress
			stress_i = stress[i]

			stress_p = 0.5 * (stress_i[0, 0] + stress_i[1, 1]) + numpy.sqrt( (0.5 * (stress_i[0, 0] - \
				stress_i[1, 1]))**2.0 + stress_i[0, 1] )

			if stress_i[0,0] - stress_i[1,1]: # otherwise this is prolly a boundary particle

				# Compute angle
				stress_prin[i,2] = 0.5 * numpy.arctan( 2.0 * stress_i[0,1] / (stress_i[0,0] \
								  - stress_i[1,1]) )

				# compute principal stress
				stress_prin[i,0] = stress_p * numpy.cos(stress_prin[i,2])
				stress_prin[i,1] = stress_p * numpy.sin(stress_prin[i,2])

		stress_norm = norm(stress_prin[:,:-1], axis=1)

		if not peters: # plot all particle stresses
			
			import matplotlib.pylab as plt # imported here to make sure PyGran works without matplotlib

			indices  = numpy.where(stress_norm > 0)[0]
			Particles = self._Particles[indices]

			coords = numpy.array([Particles.x, Particles.y, Particles.z])
			x,y = coords[axes[0]], coords[axes[1]]

			stress_prin = stress_prin[indices,:]
			stress_norm = stress_norm[indices]

			stress_norm /= stress_norm.max() * 0.3

			plt.axes()

			for i in range(coords.shape[1]):
				theta = stress_prin[i,-1]
				xi = [x[i] - Particles.radius[i] * numpy.cos(theta), x[i] + Particles.radius[i] * numpy.cos(theta)]
				yi = [y[i] - Particles.radius[i] * numpy.sin(theta), y[i] + Particles.radius[i] * numpy.sin(theta)]

				if plot_stress:
					plt.plot(xi, yi, linewidth=stress_norm[i], color='black')

				if plot_parts:
					circle = plt.Circle((x[i], y[i]), radius=Particles.radius[i], fill=False, color='blue', linewidth=0.5, linestyle ='dotted')
					plt.gca().add_patch(circle)

					circle = plt.Circle((x[i], y[i]), radius=Particles.radius[i], fill=False, color='blue', linewidth=0.5, linestyle ='dotted')
					plt.gca().add_patch(circle)

			plt.axis('scaled')
			plt.show()

			return stress_prin

		else: # use Peters' algorithm for computing force chain propagation
			# Section V. THE ALGORITHM (page 4)
			# Step 1: filter out particles with less than the mean principal stress
			indices  = numpy.where(stress_norm >= threshold * stress_norm.mean(axis=0))[0]

			stress_norm = stress_norm[indices]

			# Step 2: filter out particles which are in contact with 1 or less 'highly stressed' particles
			# Construct an nns list
			coords = self._coords[indices,:][:,axes] 
			stress_prin = stress_prin[indices,:]
			radii = self._Particles.radius[indices]

			tree = spatial.cKDTree(coords)

			nns = tree.query_ball_point(coords, self._Particles[indices].radius.max() * 2.0)
			self._nns = [[] for i in range(len(nns))]

			# Update the nns indices so that only overlapping neighbors remain
			count = 0

			for ns in nns:
				for index in ns:
					if index != count:

						dist = norm(coords[index] - coords[count])

						if dist <= (radii[index] + radii[count]):
							self._nns[count].append(index)

				self._nns[count].insert(0, count)
				count += 1

			cosAlpha = numpy.cos(alpha)
			chain = []

			# Step 3: Remove overlapping particles with 2 or less neighbors 
			remInd = []
			for nns in self._nns:
				if len(nns) <= 2: # select only particles surrounded by 2 or more highly stressed particles
					remInd.append(nns)

			for ind in remInd:
				self._nns.remove(ind)

			# Step 4: compute the force chain per particle and make plots if requested by the user
			plt.axes()

			for nns in self._nns:
				
				index = nns[0]
				# Begin loop over each neighboring particle
				for ni in nns[1:]:

					# check for angle compliance
					angle_i = stress_prin[index,-1]
					angle_j =  stress_prin[ni,-1]

					dist = coords[ni] - coords[index]

					term1 = numpy.fabs(numpy.dot(dist, stress_prin[index,:-1]))
					term2 = norm(dist) * stress_norm[index]

					term3 = numpy.fabs(numpy.dot(-dist, stress_prin[ni,:-1]))
					term4 = norm(-dist) * stress_norm[ni]

					# Search in the forward direction
					if (term1 <= term2) and (cosAlpha * term2) < term1:
						if (term3 <= term4) and (cosAlpha * term4) < term3:
							chain.append((index,ni))
					
							x = [coords[index,0], coords[ni,0]]
							y = [coords[index,1], coords[ni,1]]

							if plot_stress:
								plt.plot(x, y, linewidth=max(1.2 * stress_norm[i] / stress_norm.max(), 0.01), color='black')

							if plot_parts:
								circle = plt.Circle((coords[index,0], coords[index,1]), radius=radii[index], fill=False, color='blue', linewidth=0.5, linestyle ='dotted')
								plt.gca().add_patch(circle)

								circle = plt.Circle((coords[ni,0], coords[ni,1]), radius=radii[ni], fill=False, color='blue', linewidth=0.5, linestyle ='dotted')
								plt.gca().add_patch(circle)

					# Search in the reverse direction
					# I dont understand what Im doing here
					#if (term1 <= term2) and (-cosAlpha * term2) > term1: # TODO: double check the math on this one
						#if (term3 <= term4) and (-cosAlpha * term4) > term3: # TODO: double check the math on this one
							#chain.append((index,ni))
							#print (index,ni)
			
			print('Natoms highly stressed = {} out of {} particles'.format(len(chain), self._Particles.natoms))
				
			plt.axis('scaled')
			plt.show()

			return chain, stress_prin
