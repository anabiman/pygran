'''
Created on July 11, 2016
@author: Andrew Abi-Mansour
'''

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