'''
Created on July 10, 2016
@author: Andrew Abi-Mansour

Center for Materials Sci. & Eng.,
Merck Inc., West Point
'''

# !/usr/bin/python
# -*- coding: utf8 -*- 
# -------------------------------------------------------------------------
#
#   Python module for creating the basic DEM (Granular) object for analysis
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

import numpy as np
from xlrd import open_workbook
from numbers import Number

class Granular(object):
	"""The Granular class contains all the information describing a ganular system.
	A system always requires a trajectory file to read. A trajectory is a (time) 
	series corresponding to the coordinates of all particles in the system. It can 
	also contain other variables such as momenta, angular velocities, forces, radii,
	etc. """

	def __init__(self, fname):

		self._fname = fname
		self._fp = open(fname, 'r')

		self.nFrames = None
		self.frame = 0
		self.data = {}

		self.__next__()

	def __iter__(self):
		return self

	def extract(self, key):

		if key in self.data:
			return self.data[key]
		else:
			return None

	def goto(self, frame):
		""" Go to a specific frame in the trajectory """

		if frame == self.frame:
			return 0

		# rewind if necessary (better than reading file backwads?)
		if frame < self.frame:
			self.rewind()


		# find the right frame number
		while self.frame < frame:

			line = self._fp.readline()

			if not line:
				raise StopIteration

			if line.find('TIMESTEP') >= 0:
				self.frame += 1

		# assert self.frame == frame else something's wrong
		if self.frame == frame:

			timestep = int(self._fp.readline())
			self.data['TIMESTEP'] = timestep

			while True:

				line = self._fp.readline()

				if not line:
					raise StopIteration

				if line.find('NUMBER OF ATOMS') >= 0:
					natoms = int(self._fp.readline())
					self.data['NATOMS'] = natoms

				if line.find('BOX') >= 0:
					boxX = self._fp.readline().split()
					boxY = self._fp.readline().split()
					boxZ = self._fp.readline().split()

					boxX = [float(i) for i in boxX]
					boxY = [float(i) for i in boxY]
					boxZ = [float(i) for i in boxZ]

					self.data['BOX'] = (boxX, boxY, boxZ)
					break

			line = self._fp.readline()

			if not line:
					raise StopIteration

			keys = line.split()[2:] # remove ITEM: and ATOMS keywords

			for key in keys:
				self.data[key] = np.zeros(natoms)

			for i in range(natoms):
				var = self._fp.readline().split()

				for j, key in enumerate(keys):
					self.data[key][i] = float(var[j])
		else:
			print 'Cannot find frame {} in current trajectory'.format(frame)

	def rewind(self):
		"""Read trajectory from the beginning"""
		self._fp.close()
		self._fp = open(self._fname)
		self.frame = 0
		self.__next__()

	def __next__(self):
		"""Forward one step to next frame when using the next builtin function."""
		return self.next()

	def next(self):
		while True:

			line = self._fp.readline()

			if not line:
				raise StopIteration
			
			if line.find('TIMESTEP') >= 0:
				timestep = int(self._fp.readline())
				self.data['TIMESTEP'] = timestep

			if line.find('NUMBER OF ATOMS') >= 0:
				natoms = int(self._fp.readline())
				self.data['NATOMS'] = natoms

			if line.find('BOX') >= 0:
				boxX = self._fp.readline().split()
				boxY = self._fp.readline().split()
				boxZ = self._fp.readline().split()

				boxX = [float(i) for i in boxX]
				boxY = [float(i) for i in boxY]
				boxZ = [float(i) for i in boxZ]

				self.data['BOX'] = (boxX, boxY, boxZ)
				break

		line = self._fp.readline()
		if not line:
					raise StopIteration

		self.frame += 1

		keys = line.split()[2:] # remove ITEM: and ATOMS keywords

		for key in keys:
			self.data[key] = np.zeros(self.data['NATOMS'])

		for i in range(self.data['NATOMS']):
			var = self._fp.readline().split()

			for j, key in enumerate(keys):
				self.data[key][i] = float(var[j])

		return timestep

	@property
	def granular(self):
		return self

	def __del__(self):
		self._fp.close()



def select(data, *region):	
	"""
	Create a selection of particles based on a subsystem defined by region
	@ region: a tuple (xmin, xmax, ymin, ymax, zmin, zmax). If undefined, all particles are considered.
	"""

	try:
		if not len(region):
			return np.arange(data['NATOMS'], dtype='int')
		else:
			if len(region) != 6:
				print 'Length of region must be 6: (xmin, xmax, ymin, ymax, zmin, zmax)'
				raise

		xmin, xmax, ymin, ymax, zmin, zmax = region

		x, y, z = data['x'], data['y'], data['z']

		# This is so hackish!
		if len(x) > 0:

			indices = np.intersect1d(np.where(x > xmin)[0], np.where(x < xmax)[0])
			indices = np.intersect1d(np.where(y > ymin)[0], indices)
			indices = np.intersect1d(np.where(y < ymax)[0], indices)

			indices = np.intersect1d(np.where(z > zmin)[0], indices)
			indices = np.intersect1d(np.where(z < zmax)[0], indices)

		else:
			indices = []

		return indices

	except:
		raise