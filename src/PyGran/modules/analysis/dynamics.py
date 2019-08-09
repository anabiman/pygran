'''
This module provides classes for time-dependent analysis

Created on July 11, 2016

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

import numpy

class Temporal(object):
	""" A generic class that enables time-dependent analysis of DEM data """

	def __init__(self, System):
		self.System = System

	def series(self, attr):
		"""
		Extracts a time series for a supplied *attr* stored in self.System

		:param attr: attribute to extract time series for
		:type attr: string
		:return: time series
		:rtype: numpy array of dimensions nframes x dim_attr
		"""
		out = []
		self.System.rewind()

		for ts in self.System:
			if attr in self.System.Particles.data:
				out.append(self.System.Particles.__getattribute__(attr))

		self.System.rewind()

		return numpy.array(out)
		
	def flow(self, density=None, dt=1):
		"""
		Computes flow rate. When *density* is supplied, the computed rate is the mass per unit of time, otherwise it is the number of particles per unit of time.

		:param density: true mass density
		:type density: float 
		:param dt: timestep in units of time, default 1.
		:type dt: float
		:return: flow rate
		:rtype: 1D numpy array

		.. todo:: 
			Make this routine more efficient
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

		:param xmin: bound that specifies the left-hand wall pependicular to the y-z plane
		:type xmin: float
		:param xmax: bound that specifies the right-hand wall pependicular to the y-z plane
		:type xmax: float
		:param ymin: bound that specifies the right-hand wall pependicular to the x-z plane
		:type ymin: float
		:param ymax: bound that specifies the right-hand wall pependicular to the x-z plane
		:type ymax: float
		:param zmin: bound that specifies the right-hand wall pependicular to the x-y plane
		:type zmin: float
		:param zmax: bound that specifies the right-hand wall pependicular to the x-y plane
		:type zmax: float
		:return: wall collision frequency
		:rtype: float

		.. note:: This routines supports only 3D systems
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