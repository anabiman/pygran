#  -*- coding: utf8 -*-

'''
  Created on March 10, 2017
  @author: Andrew Abi-Mansour

  This is the 
   __________         ________                     
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

 -------------------------------------------------------------------------
  Python module for determining contacts/overlaps in an N-particle system
 -------------------------------------------------------------------------

 '''

import numpy
from scipy import spatial
from scipy.linalg import norm
from .models import SpringDashpot
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
	    return self._distances

	@property
	def pairs(self):
	    return self._pairs

	@property
	def overlaps(self):
	    return self._overlaps

	def filter(self, percent=None):
		""" Returns a non-overlapping Particles class from the given configuration 
		
		@[percent]: filter all particles overlapping by a certain percentage
		@
		"""

		percent /= 100.0

		if percent:
			indices = numpy.array(self._overlaps[:,1:], 'int')
			ri, rj = self._Particles.radius[indices[:,0]], self._Particles.radius[indices[:,1]]
			radii = ri * rj / (ri + rj)

			indices = numpy.unique(numpy.array(indices[self._overlaps[:,0] <= percent * radii, 1:]))

			return self._Particles[indices]
		
		count = numpy.zeros(len(self._Particles))

		for ind in indices:
			i,j = ind
			count[i] += 1
			count[j] += 1
			
		indices = numpy.where(count <= 1)[0]

		return self._Particles[indices]

	def forceChain(self, axis = (0,2), alpha = numpy.pi/4, plot_stress=True, plot_parts=False, peters=True, threshold=1):
		""" Computes the force chain based on an algorithm published in Phys. Rev. E. 72, 041307 (2005):
		'Characterization of force chains in granular material'.

		@ axis: a tuple of size 2 that specifies the two axis to use for computing the force chain, e.g. axis=(0,1) -> (x,y)
		@ alpha: the angle (in radians) that controls the deviation of the force chain. A value of 0 means a perfectly linear chain.
		Thus, alpha is a measure of the 'curvature' of the force chain. See page 5 of the paper cited above.
		"""

		stress = numpy.zeros((self._Particles.natoms, 2, 2))
		stress_prin = numpy.zeros((self._Particles.natoms, 3)) # 2 stress components + angle = 4 dims

		coords = numpy.array([self._Particles.x, self._Particles.y, self._Particles.z])
		x,y = coords[axis[0]], coords[axis[1]]

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
			
			indices  = numpy.where(stress_norm > 0)[0]
			Particles = self._Particles[indices]

			coords = numpy.array([Particles.x, Particles.y, Particles.z])
			x,y = coords[axis[0]], coords[axis[1]]

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
			coords = self._coords[indices,:][:,axis] 
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
