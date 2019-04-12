#  -*- coding: utf8 -*-

'''
  Created on July 11, 2016
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
 '''

import numpy

class Temporal(object):
	def __init__(self, System):
		self.System = System

	def series(self, attr):
		out = []
		self.System.rewind()

		for ts in self.System:
			if attr in self.System.Particles.data:
				out.append(self.System.Particles.__getattribute__(attr))

		self.System.rewind()

		return numpy.array(out)
		
	def flow(self, density=None, dt=1):
		"""
		Computes flow rate
		@[density]: true mass density. When supplied, the computed rate is the mass per unit of time, otherwise it is the number of particles per unit of time
		@[dt]: timestep in units of time, default 1.

		TODO: Make it more efficient
		"""

		self.System.rewind()

		N0 = self.System.Particles.natoms
		t0 = self.System.Particles.timestep

		mass = []

		count = 0

		for ts in self.System:
			
			if self.System.Particles.natoms:

				time = (self.System.Particles.timestep - t0) * dt
				mass.append( numpy.sum(- density * 4.0 / 3.0 * numpy.pi * (self.System.Particles.natoms - N0) * (self.System.Particles.radius**3.0) / time) )
				
		self.System.rewind()

		return numpy.array(mass)


	def wallCollision(self, **boundary):
		"""
		Computes the frequency of particle-wall collisions

		@[xmin]: float defining the left-hand wall pependicular to the y-z plane
		@[xmax]: float defining the right-hand wall pependicular to the y-z plane
		@[ymin]: float defining the right-hand wall pependicular to the x-z plane
		@[ymax]: float defining the right-hand wall pependicular to the x-z plane
		@[zmin]: float defining the right-hand wall pependicular to the x-y plane
		@[zmax]: float defining the right-hand wall pependicular to the x-y plane
		"""

		if not boundary:
			raise ValueError("Bounding walls must be specified with xmin, xmax, ... zmax.")

		collisions = 0

		if 'xmin' in boundary:
			collisions += len(self.System.Particles.x - self.System.radius <= boundary['xmin'])

		if 'xmax' in boundary:
			collisions += len(self.System.Particles.x - self.System.radius >= boundary['xmax'])

		if 'ymin' in boundary:
			collisions += len(self.System.Particles.y - self.System.radius <= boundary['ymin'])

		if 'ymax' in boundary:
			collisions += len(self.System.Particles.y - self.System.radius >= boundary['ymax'])

		if 'zmin' in boundary:
			collisions += len(self.System.Particles.z - self.System.radius <= boundary['zmin'])

		if 'zmax' in boundary:
			collisions += len(self.System.Particles.z - self.System.radius >= boundary['zmax'])

		return collisions