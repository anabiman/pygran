'''
Created on July 11, 2016
@author: Andrew Abi-Mansour
'''

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

class Temporal(object):
	def __init__(self, System):
		self.System = System

	def timeSeries(self, att):
		out = []
		for ts in self.System:
			out.append(self.System.Particles.__getattribute__(att))

		return numpy.array(out)
		
	def computeFlow(self, density, dt):
		"""
		Computes flow rate: N/t for a selection *sel*
		@ data: list of dictionaries containing simulation and particle data (box size, x,y,z, etc.)

		TODO: Get this working for a particle distribution
		"""

		natoms = self.timeSeries('natoms')
		time = self.timeSeries('timestep') * dt

		N0 = natoms[0]
		t0 = time[0]
		mass = [None for i in range(len(natoms))]

		count = 0

		for ts in self.System:
			
			mass[count] = density * 4.0 / 3.0 * np.pi * (natoms - N0) * (self.Particles.radius**3.0)
			count += 1
			
		return - numpy.array(mass) / (time - t0)