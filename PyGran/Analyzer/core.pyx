'''
Created on July 10, 2016
@author: Andrew Abi-Mansour
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
cimport numpy as np
import types
from random import choice
from string import ascii_uppercase
import glob
import re, sys
import collections
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from PyGran.Tools import siToMicro, microToSi

class SubSystem(object):
	""" The SubSystem class stores all particle properties and the methods that operate on \
these properties """

	def __init__(self, sel = None, units = 'si', **args):


		self._fname = args['fname']
		self._singleFile = None
		self._ftype = None
		self._fp = None

		# If data is already coped from Particles, do nothing
		if not hasattr(self, 'data'):
			self.data = collections.OrderedDict() # am ordered dict that contains either arrays
			#(for storing pos, vels, forces, etc.), scalars (natoms, ) or tuples (box size). 
			# ONLY arrays can be slices based on user selection.

		self._index = -1 # used for for loops only
		self._units = units

		
	def _metaget(self, key):
		"""A meta function for returning dynamic class attributes treated as lists (for easy slicing) 
		and return as numpy arrays for numerical computations / manipulations"""
		return self.data[key]

	def _constructAttributes(self, sel=None):
		""" Constructs dynamic functions (getters) for all keys found in the trajectory """

		for key in self.data.keys():
			if type(self.data[key]) == np.ndarray:

				if sel is not None:
					self.data[key] = self.data[key][sel]
					
		# Checks if the trajectory file supports reduction in key getters
		# It's important to construct a (lambda) function for each attribute individually
		# so that each dynamically created (property) function would have its own unique
		# key. For instance, placing the for loop below in  _constructAttributes would make
		# all property functions return the same (last) key variable.


			# We cannot know the information for any property function until that property is created, 
			# so we define either define metaget function and particularize it only later with a 
			# lambda function, thus, permanently updating the SubSystem class, or update the attributes
			# of this particular instance of SubSystem, which is the approach adopted here.
			setattr(self, key, self._metaget(key))

	def __getitem__(self, sel):
		""" SubSystem can be sliced with this function """
		
		# Get the type of the class (not necessarily SubSystem for derived classes)
		cName = eval(type(self).__name__)

		return cName(sel, self._units, **self.data)

	def __and__(self, ):
		""" Boolean logical operator on particles """ 
				

	def copy(self):
		""" Returns a hard copy of the SubSystem """
		data = self.data.copy()

		for key in self.data.keys():
			if type(data[key]) == np.ndarray:
				data[key] = data[key].copy()
			else:
				data[key] = data[key]

		return eval(type(self).__name__)(None, self._units, **data)

	def conversion(self, factors):
		""" Convesion factors from S.I. to micro and vice versa """

		for key in self.data.keys():

			if key == 'x' or key == 'y' or key == 'z' or key == 'radius':
				self.data[key] *= factors['distance']
			elif key == 'vx' or key == 'vy' or key == 'vz':
				self.data[key] *= factors['distance'] / factors['time']
			elif key == 'omegax' or key == 'omegay' or key == 'omegaz':
				self.data[key] /= factors['time']
			elif key == 'fx' or key == 'fy' or key == 'fz':
				self.data[key] *= factors['mass'] * factors['distance'] / factors['time']**2.0
			else:
				pass

	def units(self, units):
		""" Sets the unit system """
		if self._units == 'si':
			if units == 'micro':
				self.conversion(siToMicro)
		elif self._units == 'micro':
			if units == 'si':
				self.conversion(microToSi)

		self._units = units

	def __iadd__(self, obj):
		""" Adds attributes of obj to current SubSystem """

		if type(obj) is type(self):
			for key in self.data.keys():
				if key in obj.data:
					if type(self.data[key]) is np.ndarray:
						self.data[key] = np.concatenate((self.data[key], obj.data[key]))

		self.__init__(None, self.units, **self.data)

		return self

	def __add__(self, obj):
		""" Adds two classes together, or operates scalars/vectors on particle radii/positions
		TODO: get this working with meshes """

		if type(obj) is type(self):

			if len(obj) == len(self):

				if len(self.data.keys())  < len(obj.data.keys()):
					data = self.data.copy()
					datac = obj.data
				else:
					data = obj.data.copy()
					datac = self.data

				for key in datac.keys():
					if key in data:
						if type(datac[key]) is np.array:
							data[key] = data[key].copy()
							
						data[key] += datac[key]
					else:
						if type(datac[key]) is np.array:
							data[key] = datac[key].copy()
						else:
							data[key] = datac[key]

				return eval(type(self).__name__)(None, self._units, **data)
			else:
				print('Two subsystems with different numbers of elements cannot be added.')
				return None

	def __sub__(self, obj):
		""" Subtracts scalars/vectors from particle radii/positions
		TODO: get this working with meshes """

		if type(obj) is tuple:
			obj, att = obj
			if type(obj) is type(self):
				if att is 'all':
					pass
				elif len(obj) == len(self): # gotta make sure both classes being subtracted have the same number of elements
					self.data[att] -= obj.data[att]
				else:
					print 'Two subsystems with different number of elements cannot be subtracted.'

		return self

	def __mul__(self, obj):
		""" Subtracts scalars/vectors from particle radii/positions
		TODO: get this working with meshes """

		if type(obj) is tuple:
			obj, att = obj
			if type(obj) is type(self):
				if att is 'all':
					pass
				elif len(obj) == len(self): # gotta make sure both classes being multiplied have the same number of elements
					self.data[att] *= obj.data[att]
				else:
					print 'Two subsystems with different number of elements cannot be multiplied.'

		return self

	def __div__(self, obj):
		""" Subtracts scalars/vectors from particle radii/positions
		TODO: get this working with meshes """

		if type(obj) is tuple:
			obj, att = obj
			if type(obj) is type(self):
				if att is 'all':
					pass
				elif len(obj) == len(self): # gotta make sure both classes being divided have the same number of elements
					self.data[att] /= obj.data[att]
				else:
					print 'Two subsystems with different number of elements cannot be divided.'

		return self

	def __or__(self, att):
		""" Returns a (temporary) new subsystem with only a single attribute. the user can of course make this not *temporary* but 
		that is not what should be used for. In principle, the modulus is a reductionist operator that serves as a temporary metastate
		for binary operations. """

		# Make a hard copy of the class to make sure we preserve its state
		obj = self.copy()
		data = obj.data

		for key in self.data.keys():
			if key != att and key != 'natoms':
				del data[key]

		return eval(type(self).__name__)(None, obj._units, **data)

	def __iter__(self):
		self._index = -1
		return self

	def __next__(self):
		"""Forward one step to next frame when using the next builtin function."""
		return self.next()

	def next(self):
		""" This method is invoked by __iter__ in a for loop """

		if self._index < len(self) - 1:
			self._index += 1
			return self[self._index]
		else:
			raise StopIteration

	def __len__(self):
		""" Returns the number of elements (e.g. particles or nodes) in a SubSystem """
		for key in self.data.keys():
			if type(self.data[key]) == np.ndarray:
				return len(self.data[key])

		return -1

	def __del__(self):
		pass

class Mesh(SubSystem):
	"""  The Mesh class stores a list of meshes and their associated attributes / methods.
	"""
	def __init__(self, fname, units='si', vtk_type=None):

		self._vtk = vtk_type
		self._units = units

		if vtk_type == 'poly':
			self._reader = vtk.vtkPolyDataReader()
		else:
			self._reader = vtk.vtkUnstructuredGridReader()

		super(Mesh, self).__init__(None, self._units, fname=fname)

		self._reader.SetFileName(fname)
		self._reader.Update() # Needed if we need to call GetScalarRange
		self._output = self._reader.GetOutput()

		# Assert mesh input filname is VTK
		if self._mfname:
			if self._mfname.split('.')[-1] != 'vtk':
				raise IOError('Input mesh must be of VTK type.')

			self._mfname = sorted(glob.glob(self._mfname), key=numericalSort)
			self._mesh = self._mfname[0]


		try:
			points = self._output.GetPoints().GetData()
		except:
			points = None

		try:
			cells = self._output.GetCells().GetData()
		except:
			try:
				cells = self._output.GetCellData()
			except:
				cells = None

		if points:
			self.data[points.GetName()] = vtk_to_numpy(points)

		if cells:
			try:
				self.data['Cells'] = vtk_to_numpy(cells)
			except:
				pass

		if points:
			if cells:

				np_pts = np.zeros((self.nCells(), self._output.GetCell(0).GetPoints().GetNumberOfPoints(), self.data[points.GetName()].shape[1]))

				# Area/volume computation for elements is deprecated in VTK (only support for tetrahedra), so we do the computation ourselves here for rectangular elements.
				if self._output.GetCell(0).GetNumberOfFaces() == 6: # we're dealing with a rectangular element
					self.data['CellVol'] = np.zeros(self.nCells())
				elif self._output.GetCell(0).GetNumberOfFaces() == 0: # 2D mesh
					self.data['CellArea'] = np.zeros(self.nCells())

				for i in range(self.nCells()):
					pts = self._output.GetCell(i).GetPoints()
					np_pts[i] = vtk_to_numpy(pts.GetData())

				# Find the axis of alignment for all cells using 1 (arbitrary) node
				if np.diff(np_pts[:,0,0]).any() == .0:
					indices = 1,2
				elif np.diff(np_pts[:,0,1]).any() == .0:
					indices = 0,2
				else:
					indices = 0,1

				for i in range(self.nCells()):
					if 'CellVol' in self.data:
						# Need to update this to compoute correct the cell volume (assumed here to be a box)
						self.data['CellVol'][i] = (np_pts[i][:,0].max() - np_pts[i][:,0].min()) * (np_pts[i][:,1].max() - np_pts[i][:,1].min()) * (np_pts[i][:,2].max() - np_pts[i][:,2].min())

					if 'CellArea' in self.data:
						# Shoelace  algorithm for computing areas of polygons
						self.data['CellArea'][i] = 0.5 * np.abs( np.dot(np_pts[i][:,indices[0]] ,np.roll(np_pts[i][:,indices[1]],1)) - np.dot(np_pts[i][:,indices[1]],np.roll( np_pts[i][:,indices[0]],1)) )

				self.data['CellsPos'] = np_pts

		index = 0
		while True:
			key = self._output.GetCellData().GetArrayName(index)
			if key:
				self.data[key] = vtk_to_numpy(self._output.GetCellData().GetArray(key))
				index += 1
			else:
				break

	def _updateSystem(self):
		self.__init__(self._mesh, self._units, self._vtk)

	def nCells(self):
		return self._output.GetNumberOfCells()

	def nPoints(self):
		return self._output.GetNumberOfPoints()

	def _goto(self, frame):

		if self._fname:
			if frame >= len(self._fname):
				print 'Input frame exceeds max number of frames'
			else:
				if frame == self.frame:
					pass
				else:
					if frame == -1:
						frame = len(self._fname) - 1

					self._mesh = self._fname[frame]
					self.frame = frame

	def rewind(self):
		if self._fname:
			self._mesh = self._fname[self.frame]

	def __del__(self):

		SubSystem.__del__(self)

class Particles(SubSystem):
	""" The Particle class stores all particle properties and the methods that operate on \
	these properties """

	def __init__(self, sel = None, units = 'si', **data):


		if 'Particles' in data:
			self = data['Particles'] # do a hard copy here?

			if 'natoms' not in self.data: 
				self.data['natoms'] = len(self)
		else:
			self._units = units
			fname = data['fname']

			# Python 3.X: just do super()
			super(Particles, self).__init__(sel, units, fname=fname)

			if self._fname:
				self._ftype = fname.split('.')[-1]

				if self._ftype == 'dump': # need a way to figure out this is a LIGGGHTS/LAMMPS file

					if self._fname.split('.')[:-1][0].endswith('*'):
						self._singleFile = False
						self._files = sorted(glob.glob(self._fname), key=numericalSort)
						self._fp = open(self._files[0], 'r')
					else:
						self._singleFile = True
						self._fp = open(fname, 'r')

					self._params = None

					# Do some checking here on the traj extension to make sure
					# it's supported
					self._format = fname.split('.')[-1]

					self._readFile()
					
				else:
					raise IOError('Input trajectory must be a valid LAMMPS/LIGGHTS (dump), ESyS-Particle (txt), Yade (?), or DEM-Blaze file (?)')


			self._constructAttributes(sel)

	def _updateSystem(self):

		# If we are updating (looping over traj) we wanna keep the id (ref) of the class constant 
		# (soft copy)
		self.__init__(None, self._units, **self.data)
		self._constructAttributes()

	def rog(self):
		""" Computes the radius of gyration (ROG) for an N-particle system:
		ROG = <\sqrt(\sum_i (r_i - rm)^2)> where rm is the mean position of all
		particles, and <...> is the ensemble average. Alternatively, one can
		compute ROG as \sum_i <r_i^T r_i> - rm ^T rm
		"""
		positions = np.array([self.x, self.y, self.z]).T
		rm = positions.mean(axis=0)
		N = len(positions)

		dr = positions - rm
		rog = .0

		for pos in dr:
			rog += np.dot(pos,pos)

		return np.sqrt(rog/N)

	def computeRadius(self, N = 100):
		""" Computes the maximum radius of an N-particle (spherical) system
		by sorting the radial components and returning the average of the sqrt
		of the radius of the first N max data points. 
		"""
		positions = np.array([self.x, self.y, self.z]).T
		x,y,z = self.x, self.y, self.z
		rm = positions.mean(axis=0)

		r = ((x - rm[0])**2.0 + (y - rm[1])**2.0 + (z - rm[2])**2.0)
		r.sort()

		return np.sqrt(r[-N:]).mean()

	def rdf(self, dr = None, center = True, rMax=None):
		""" Computes the three-dimensional radial distribution function for a set of
	    spherical particles contained in a cube with side length S.  This simple
	    function finds reference particles such that a sphere of radius rMax drawn
	    around the particle will fit entirely within the cube, eliminating the need
	    to compensate for edge effects.  If no such particles exist, an error is
	    returned.  Try a smaller rMax...or write some code to handle edge effects! ;)

	    Arguments:

	        S               length of each side of the cube in space
	        rMax            outer diameter of largest spherical shell
	        dr              increment for increasing radius of spherical shell
	    Implicit arguments:
	    	x               an array of x positions of centers of particles
	        y               an array of y positions of centers of particles
	        z               an array of z positions of centers of particles

	    Returns a tuple: (g, radii, interior_indices)
	        g(r)            a numpy array containing the correlation function g(r)
	        radii           a numpy array containing the radii of the
	                        spherical shells used to compute g(r)
	        reference_indices   indices of reference particles
		"""
	    
		x, y, z = self.x, self.y, self.z
		
		# center positions around 0
		if center:
			x -= x.mean()
			y -= y.mean()
			z -= z.mean()

		S = min(x.max(), y.max(), z.max())

		if rMax is None:
			rMax = S / 2.0

		if dr is None:
			dr = rMax / 100

		# Find particles which are close enough to the cube center that a sphere of radius
		# rMax will not cross any face of the cube

		print 'Constructing a cube of length {} and a circumscribed sphere of radius {}'.format(S * 2.0, rMax)
		print 'Resolution chosen is {}'.format(dr)

		bools1 = x > rMax - S
		bools2 = x < (S - rMax)
		bools3 = y > rMax - S
		bools4 = y < (S - rMax)
		bools5 = z > rMax - S
		bools6 = z < (S - rMax)

		interior_indices, = np.where(bools1 * bools2 * bools3 * bools4 * bools5 * bools6)    
		num_interior_particles = len(interior_indices)

		if num_interior_particles < 1:
			raise  RuntimeError ("No particles found for which a sphere of radius rMax will lie entirely within a cube of side length S.  Decrease rMax or increase the size of the cube.")

		edges = np.arange(0., rMax + 1.1 * dr, dr)
		num_increments = len(edges) - 1
		g = np.zeros([num_interior_particles, num_increments])
		radii = np.zeros(num_increments)
		numberDensity = num_interior_particles / S**3

		# Compute pairwise correlation for each interior particle
		for p in range(num_interior_particles):
			index = interior_indices[p]
			d = np.sqrt((x[index] - x)**2 + (y[index] - y)**2 + (z[index] - z)**2)
			d[index] = 2.0 * rMax
			(result, bins) = np.histogram(d, bins=edges, normed=False)
			g[p,:] = result

		# Average g(r) for all interior particles and compute radii
		g_average = np.zeros(num_increments)
		for i in range(num_increments):
			radii[i] = (edges[i] + edges[i+1]) / 2.
			rOuter = edges[i + 1]
			rInner = edges[i]
			g_average[i] = np.mean(g[:, i]) / (4.0 / 3.0 * np.pi * (rOuter**3 - rInner**3))

		print numberDensity
		return (g_average / numberDensity, radii, interior_indices)
		# Number of particles in shell/total number of particles/volume of shell/number density
		# shell volume = 4/3*pi(r_outer**3-r_inner**3)

	def angleRepose(self):
		"""
		Computes the angle of repos theta = arctan(h_max/L)
		in a sim box defined by [-Lx, Lx] x [-Ly, Ly] x [0, Lz]
		"""
		x, y, z = self.x, self.y, self.z
		dL = 0.25 * (x.max() - x.min()) + 0.25 * (y.max() - y.min())
		z_max = z.max() - z.min() - self.radius.mean()

		return np.arctan(z_max / dL) * 180.0 / np.pi
		
	def density(self, bdensity, shape = 'box'):
		"""
		Computes the bulk density for a selection of particles from their *true* density. 
		The volume is determined approximately by constructing a box/cylinder/cone 
		embedding the particles. Particles are assumed to be spherical in shape.
		"""

		if(self.natoms > 0):

			radius = self.radius
			volume = self.volume(shape)
			mass = np.sum(bdensity * 4.0 / 3.0 * np.pi * (radius**3.0))

			return mass / volume

		return 0

	def densityLocal(self, bdensity, dr, axis):
		"""" Computes a localized density at a series of discretized regions of thickness 'dr'
		along an axis specified by the user """
		
		if axis == 'x':
			r = self.x
		elif axis == 'y':
			r = self.y
		elif axis == 'z':
			r = self.z
		else:
			raise ValueError("Axis can be only x, y, or z.")

		thick = np.arange(r.min(), r.max(), dr)
		odensity = []

		for i in range(len(thick) - 1):
			parts = self[r <= thick[i+1]]

			if axis == 'x':
				denLoc = parts[parts.x >= thick[i]].density(bdensity)
			elif axis == 'y':
				denLoc = parts[parts.y >= thick[i]].density(bdensity)
			elif axis == 'z':
				denLoc = parts[parts.z >= thick[i]].density(bdensity)

			odensity.append( denLoc )

		return thick, odensity

	def volume(self, shape = 'box'):
		""" Computes the volume of a granular system based on a simple geometry """

		if(self.natoms > 0):

			x, y, z = self.x, self.y, self.z
			xmin, xmax = min(x), max(x)
			ymin, ymax = min(y), max(y)
			zmin, zmax = min(z), max(z)

			if shape == 'box':
				volume = (xmax - xmin) * (ymax - ymin) * (zmax - zmin)

			elif shape == 'cylinder-z':
				height = zmax - zmin
				radius = (ymax - ymin) * 0.25 + (xmax - xmin) * 0.25
				volume = np.pi * radius**2.0 * height

			elif shape == 'cylinder-y':
				height = ymax - ymin
				radius = (zmax - zmin) * 0.25 + (xmax - xmin) * 0.25
				volume = np.pi * radius**2.0 * height

			elif shape == 'cylinder-x':
				height = xmax - xmin
				radius = (ymax - ymin) * 0.25 + (zmax - zmin) * 0.25
				volume = np.pi * radius**2.0 * height

			return volume

		return 0

	def scale(self, value, attr=('x','y','z')):
		""" Scales all particles by a float or int 'value'

		@[attr]: attribute to scale (positions by default)
		"""
		for at in attr:
			if at in self.data:
				self.data[at] *= value
			
		self._constructAttributes()

	def translate(self, value, attr=('x','y','z')):
		""" Translates all particles by a float or int 'value'

		@[attr]: attribute to translate (positions by default)
		"""
		for at in attr:
			if at in self.data:
				self.data[at] += value
				
		self._constructAttributes()

	def perturb(self, percent, attr=('x','y','z')):
		""" Adds white noise to all particles by a certain 'percent'

		@[attr]: attribute to perturb (positions by default)
		"""
		for at in attr:
			if at in self.data:
				self.data[at] += percent * (1.0 - 2.0 * np.random.rand(len(self.data[at])))
				
		self._constructAttributes()

	def rewind(self):
		
		if self._fp:
			self._fp.close()

			if self._singleFile:
				self._fp = open(self._fname)
			else:
				self._fp = open(self._files[0])
		

	def _goto(self, frame):
		""" This function assumes we're reading a non-const N trajectory. 
		"""

		# find the right frame number
		if self._singleFile:
			# We must be reading a single particle trajectory file (i.e. not a mesh)
			while self.frame < frame or frame == -1:

				line = self._fp.readline()

				if not line and frame >= 0:
					raise StopIteration('End of file reached.')
				elif not line and frame == -1:
					break

				if line.find('TIMESTEP') >= 0:
					self.frame += 1

			# assert self.frame == frame else something's wrong
			# or the user wants to go to the last frame
			if self.frame == frame:

				timestep = int(self._fp.readline())
				self.data['timestep'] = timestep
				self._readDumpFile()

			else:
				if frame == -1:
					# this is very low and inefficient -- can't we move the file pointer backwards?
					tmp = self.frame
					self.rewind()
					self.goto(tmp)
				else:
					raise NameError('Cannot find frame {} in current trajectory'.format(frame))

		else: # no need to find the input frame, just select the right file if available
			
			if self._fp:
				if frame >= len(self._files):
					print 'Input frame exceeds max number of frames'
				else:
					if frame == self.frame:
						pass
					else:
						if frame == -1:
							frame = len(self._files) - 1

						self._fp.close()
						self._fp = open(self._files[frame], 'r')
						self.frame = frame
						self._readDumpFile()

		return self.frame

	def _readDumpFile(self):
		""" Reads a single dump file"""
		count = 0
		timestep = None

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

			count += 1

		line = self._fp.readline()

		if not line:
			raise StopIteration

		keys = line.split()[2:] # remove ITEM: and ATOMS keywords

		for key in keys:
			self.data[key] = np.zeros(natoms)

		for i in range(self.data['natoms']):
			var = self._fp.readline().split()
			for j, key in enumerate(keys):
				self.data[key][i] = float(var[j])

		count += self.data['natoms'] + 1

		return timestep

	def write(self, filename):
		""" Write a single output file """
		ftype = filename.split('.')[-1]
		if ftype == 'dump':
			self._writeDumpFile(filename)
		else:
			raise NotImplementedError

	def _writeDumpFile(self, filename):
		""" Writes a single dump file"""
		with  open(filename ,'a') as fp:

			if 'timestep' not in self.data:
				self.data['timestep'] = -1

			if 'box' not in self.data:
				maxR = self.data['radius'].max()
				self.data['box'] = ([self.data['x'].min() - maxR, self.data['x'].max() + maxR], \
									[self.data['y'].min() - maxR, self.data['y'].max() + maxR], \
									[self.data['z'].min() - maxR, self.data['z'].max() + maxR], \
									)

			fp.write('ITEM: TIMESTEP\n{}\n'.format(self.data['timestep']))
			fp.write('ITEM: NUMBER OF ATOMS\n{}\n'.format(self.data['natoms']))
			fp.write('ITEM: BOX BOUNDS\n')
			for box in self.data['box']:
				fp.write('{} {}\n'.format(box[0], box[1]))

			var = 'ITEM: ATOMS '
			for key in self.data.keys():
				if key != 'timestep' and key != 'natoms' and key != 'box':
					var = var + '{} '.format(key)

			fp.write(var)
			fp.write('\n')

			for i in range(self.data['natoms']):
				var = ()
				for key in self.data.keys():
					if key != 'timestep' and key != 'natoms' and key != 'box':
						if key == 'id':
							var += (int(self.data[key][i]),)
						else:	
							var += (self.data[key][i],)

				nVars = len(var)
				fp.write(('{} ' * nVars).format(*var))
				fp.write('\n')

	def _readFile(self):
		""" Read a particle trajectory file """
		if self._singleFile:
			if self._ftype == 'dump':
				timestep = self._readDumpFile()
		else:
			if self._ftype == 'dump':

				if self.frame > len(self._files) - 1:
					raise StopIteration
					
				self._fp.close()
				self._fname = self._files[self.frame]
				self._fp = open(self._fname, 'r')
				timestep = self._readDumpFile()


class SystemFactory(object):
	"""The system class contains all the information describing a DEM system.
	A system always requires at least one trajectory file (for Particles or Mesh)
	to read. A trajectory is a (time) series corresponding to the coordinates of 
	all particles/meshes in the system. System handles the time frame and controls 
	i/o operations. It contains one or more subclasses ('Paticles'/'Mesh') which store 
	all the particle/mesh attributes read from the trajectory file (variables such as 
	positions, momenta, angular velocities, forces, stresses, radii, etc.).

	@fname: filename (or list of filenames) for the particle trajectory
	@mfname: filename (or list of filenames) for the mesh trajectory
	@dname (optional): the name of the python script file used to run a PyGran simu,

	"""

	def __init__(self, units='si', **args):

		self.frame = 0

		for ss in args:
			if(self._str_to_class(ss)):
				obj = self._str_to_class(ss)(fname=args[ss])	
				setattr(self, ss, obj)

	def __iter__(self):
		self.frame = -1
		return self

	def units(self, units):
		""" Change unit system to:  si (S.I.), or micro (microns) """
		if units == 'si' or 'micro':
			for ss in self.__dict__:
				if hasattr(self.__dict__[ss], 'units'):
					self.__dict__[ss].units(units)

			self._units = units
		else:
			print 'Only S.I. and micro units currently supported'

	def goto(self, frame):
		""" Go to a specific frame in the trajectory. If frame is -1
		then this function will read the last frame """

		# nothing to do
		if frame == self.frame:
			return 0

		# rewind if necessary (better than reading file backwads?)
		if frame < self.frame and frame >= 0:
			self.rewind()

		for ss in self.__dict__:
			if hasattr(self.__dict__[ss], '_goto'):
				self.__dict__[ss]._goto(frame)

		self._updateSystem()

	@property
	def keys(self):
		""" returns the key objects found in the trajctory files """
		return self.data.keys()

	def rewind(self):
		"""Read trajectory from the beginning"""

		self.frame = 0

		for ss in self.__dict__:
			if hasattr(self.__dict__[ss], 'rewind'):
				self.__dict__[ss].rewind()

	def __next__(self):
		"""Forward one step to next frame when using the next builtin function."""
		return self.next()

	def next(self):
		""" This method updates the system attributes! """
		self.frame += 1
		timestep = self.frame

		self.frame = 0

		if self._fname:
			try:
				self._readFile()
			except:
				self.frame -=1
				print('End of trajectory file reached.')
				raise StopIteration

		if self._mfname:
			try:
				self._mesh = self._mfname[self.frame]
			except:
				self.frame -=1
				print('End of trajectory file reached.')
				raise StopIteration

		self._updateSystem()

		return timestep

	def _updateSystem(self):
		""" Makes sure the system is aware of any update in its attributes caused by
		a frame change. """

		# A generic way of invoking _updateSystem
		for ss in self.__dict__:
			if hasattr(self.__dict__[ss], '_updateSystem'):
				self.__dict__[ss]._updateSystem()

	@property
	def system(self):
		return self

	def _str_to_class(self, str):
		try:
			return getattr(sys.modules[__name__], str)
		except:
			return None

def numericalSort(value):
	""" A sorting function by numerical numbers for glob.glob """

	numbers = re.compile(r'(\d+)')
	parts = numbers.split(value)
	parts[1::2] = map(int, parts[1::2])
	return parts

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
