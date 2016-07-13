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
import types

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

		# Read frame 0 to initialize function getters
		self.__next__()
		self.constructGetters()

	def __iter__(self):
		return self

	def constructGetters(self):
		""" Constructs dynamic functions (getters) for all keys found in the trajectory """
		for key in self.keys:
			getter = 'def get_{}(self): return self.data["{}"] \n'.format(key, key) + \
					 'self.get_{} = types.MethodType(get_{}, self)'.format(key, key)
			exec(getter)

	def _chkGetters(self):
		""" Checks if the trajectory file supports reduction in key getters """

		# TODO: make sure difference styles are supported depending on the trajectory extension
		keys = self.keys

		if self._fname.split('.')[-1] == 'dump':
			if 'x' in keys and 'y' in keys and 'z' in keys:
				self.data['positions'] = np.array([self.data['x'], self.data['y'], self.data['z']]).T

			if 'vx' in keys and 'vy' in keys and 'vz' in keys:
				self.data['velocities'] = np.array([self.data['vx'], self.data['vy'], self.data['vz']]).T

			if 'omegax' in keys and 'omegay' in keys and 'omegaz' in keys:
				self.data['angVelocities'] = np.array([self.data['omegax'], self.data['omegay'], self.data['omegaz']]).T

			if 'fx' in keys and 'fy' in keys and 'fz' in keys:
				self.data['forces'] = np.array([self.data['fx'], self.data['fy'], self.data['fz']]).T

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
			self.data['timestep'] = timestep

			while True:

				line = self._fp.readline()

				if not line:
					raise StopIteration

				if line.find('NUMBER OF ATOMS') >= 0:
					natoms = int(self._fp.readline())
					self.data['natoms'] = natoms

				if line.find('BOX') >= 0:
					boxX = self._fp.readline().split()
					boxY = self._fp.readline().split()
					boxZ = self._fp.readline().split()

					boxX = [float(i) for i in boxX]
					boxY = [float(i) for i in boxY]
					boxZ = [float(i) for i in boxZ]

					self.data['box'] = (boxX, boxY, boxZ)
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

			# for convenience, concatenate the x,y,z components
			self._updateSystem()

		else:
			raise NameError('Cannot find frame {} in current trajectory'.format(frame))

	@property
	def keys(self):
		""" returns the key objects found in the trajctory files """
		return self.data.keys()

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
		""" This method updates the system attributes! """
		while True:
			line = self._fp.readline()

			if not line:
				raise StopIteration
			
			if line.find('TIMESTEP') >= 0:
				timestep = int(self._fp.readline())
				self.data['timestep'] = timestep

			if line.find('NUMBER OF ATOMS') >= 0:
				natoms = int(self._fp.readline())
				self.data['natoms'] = natoms

			if line.find('BOX') >= 0:
				boxX = self._fp.readline().split()
				boxY = self._fp.readline().split()
				boxZ = self._fp.readline().split()

				boxX = [float(i) for i in boxX]
				boxY = [float(i) for i in boxY]
				boxZ = [float(i) for i in boxZ]

				self.data['box'] = (boxX, boxY, boxZ)
				break

		line = self._fp.readline()
		if not line:
					raise StopIteration

		self.frame += 1

		keys = line.split()[2:] # remove ITEM: and ATOMS keywords

		for key in keys:
			self.data[key] = np.zeros(self.data['natoms'])

		for i in range(self.data['natoms']):
			var = self._fp.readline().split()

			for j, key in enumerate(keys):
				self.data[key][i] = float(var[j])

		self._updateSystem()

		return timestep

	def _updateSystem(self):
		""" Makes sure the system is aware of any update in its attribvtes caused by
		a frame change. """
		self._chkGetters()

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