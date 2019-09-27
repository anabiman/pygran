'''
A module that provides the fundamental classes used in the analysis module

Created on July 10, 2016

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

import numpy as np

import types
from random import choice
from string import ascii_uppercase
import glob
import re, sys, os
import collections

try:
	import vtk
	from vtk.util.numpy_support import vtk_to_numpy
except:
	pass
	
from ..tools import convert
from scipy.stats import binned_statistic
import importlib
import numbers

class SubSystem(object):
	""" The SubSystem is an abstract class the implementation of which stores all DEM object properties and the methods that operate on \
these properties. This class is iterable but NOT an iterator. 

	:param units: unit system
	:type units: string 

	:note: see `link <https://github.com/Andrew-AbiMansour/PyGran/blob/master/src/PyGran/modules/tools/__init__.py>`_ for available unit systems
	"""

	def __init__(self, **args):

		self._units = 'si'
		self._fname = None
		self._args = args

		if type(self).__name__ in args:

			# copy constructor
			# delete data if it exists
			for key in list(self.__dict__.keys()):
				delattr(self, key)

			if hasattr(args[type(self).__name__], 'data'):
				self.data = args[type(self).__name__].data
			else:
				raise IOError('data attribute not found in {} object'.format(type(self).__name__))

			if hasattr(args[type(self).__name__], '_units'):
				self._units = args[type(self).__name__]._units

			if hasattr(args[type(self).__name__], '__module__'):
				self.__module__ = args[type(self).__name__].__module__

			self._constructAttributes()

		else: # otherwise, we're instantiating a new object
			if 'units' in args:
				self._units = args['units']

			if 'data' in args:
				self.data = args['data']

				if 'units' in args:
					self._units = args['units']
				else:
					raise RuntimeError('units must be supplied when creating SubSystem from a data dictionary.')

			if 'fname' in args:
				self._fname = args['fname']

		if '__module__' in args:
			self.__module__ = args['__module__']

		# If data is already copied from Particles, do nothing
		if not hasattr(self, 'data'):
			self.data = collections.OrderedDict() # am ordered dict that contains either arrays
			#(for storing pos, vels, forces, etc.), scalars (natoms, ) or tuples (box size).
			# ONLY arrays can be slices based on user selection.

	@property
	def keys(self):
		""" Returns all stored attribute keynames """
		if hasattr(self, 'data'):
			return self.data.keys()
		else:
			return None

	def _metaget(self, key):
		"""A meta function for returning dynamic class attributes treated as lists (for easy slicing)
		and return as numpy arrays for numerical computations / manipulations"""

		if isinstance(self.data[key], np.ndarray):
			self.data[key].flags.writeable = False

		return self.data[key]

	def _resetSubSystem(self):
		""" Mainly used by System for rewinding traj back to frame 0 """

		self._constructAttributes()

		return 0

	def _constructAttributes(self, sel=None, mesh=False):
		""" Constructs dynamic functions (getters) for all keys found in the trajectory """

		for key in self.data.keys():

			if hasattr(self, key):
				delattr(self, key)

			if isinstance(self.data[key], np.ndarray):

				if sel is not None:
					if isinstance(self.data[key], np.ndarray):
						self.data[key].flags.writeable = True
					
					self.data[key] = self.data[key][sel]

					if isinstance(self.data[key], np.ndarray):
						self.data[key].flags.writeable = False

		# Checks if the trajectory file supports reduction in key getters
		# It's important to construct a (lambda) function for each attribute individually
		# so that each dynamically created (property) function would have its own unique
		# key. For instance, placing the for loop below in  _constructAttributes would make
		# all property functions return the same (last) key variable.


			# We cannot know the information for any property function until that property is created,
			# so we define either define metaget function and particularize it only later with a
			# lambda function, thus, permanently updating the SubSystem class, or update the attributes
			# of this particular instance of SubSystem, which is the approach adopted here.

			if mesh:
				if '.' in key:
					obj, newkey = key.split('.')
					if not hasattr(self, obj):
						setattr(self, obj, lambda: None)

					setattr(getattr(self, obj), newkey, self._metaget(key))

			setattr(self, key, self._metaget(key))
			
			#Factory.addprop(self, key, lambda x: self.data[key])
			#Why do we need this?

		self.units(self._units)

	def __setattr__(self, name, value):

		if hasattr(self, name):
			# Make sure variable is not internal
			if not name.startswith('_') and not callable(value):
				raise Exception("{} property is read-only".format(name))

		self.__dict__[name] = value

	def __getitem__(self, sel):
		""" SubSystem can be sliced with this function """

		if hasattr(self, '__module__'):
			# Could not find cName, search for it in cwd (if it is user-defined)
			module = importlib.import_module(self.__module__)
		else:
			# Get the type of the class (not necessarily SubSystem for derived classes)
			module = importlib.import_module(__name__)
			
		cName = getattr(module, type(self).__name__)

		return cName(sel=sel, units=self._units, data=self.data.copy())

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

		return eval(type(self).__name__)(units=self._units, data=data)

	def conversion(self, factors):
		""" Convesion factors from S.I., micro, cgs, or nano, and vice versa 

		TODO: support all possible variables (is this possible???) """

		for key in self.data.keys():

			if key == 'x' or key == 'y' or key == 'z' or key == 'radius':

				if(type(self.data[key]) == np.ndarray):
					self.data[key].flags.writeable = True
				
				self.data[key] *= factors['distance'][0]

				if(type(self.data[key]) == np.ndarray):
					self.data[key].flags.writeable = False

			elif key == 'vx' or key == 'vy' or key == 'vz':

				if(type(self.data[key]) == np.ndarray):
					self.data[key].flags.writeable = True

				self.data[key] *= factors['distance'][0] / factors['time'][0]

				if(type(self.data[key]) == np.ndarray):
					self.data[key].flags.writeable = False

			elif key == 'omegax' or key == 'omegay' or key == 'omegaz':

				if(type(self.data[key]) == np.ndarray):
					self.data[key].flags.writeable = True

				self.data[key] /= factors['time'][0]
				
				if(type(self.data[key]) == np.ndarray):
					self.data[key].flags.writeable = False

			elif key == 'fx' or key == 'fy' or key == 'fz':

				if(type(self.data[key]) == np.ndarray):
					self.data[key].flags.writeable = True
				
				self.data[key] *= factors['mass'][0] * factors['distance'][0] / factors['time'][0]**2.0
				
				if(type(self.data[key]) == np.ndarray):
					self.data[key].flags.writeable = False
			else:
				pass

	def units(self, units=None):
		""" Change unit system to units ore return units if None
		units:  si, micro, cgs, or nano """

		if not units:
			return self._units

		self.conversion(convert(self._units, units))
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
		TODO: get this working with meshes."""

		data = collections.OrderedDict()

		if type(obj) is type(self):

			for key in self.keys:

				if key in obj.data:
					if isinstance(self.data[key], np.ndarray):
						data[key] = np.concatenate((self.data[key], obj.data[key]))
					elif isinstance(self.data[key], numbers.Number):
						data[key] = self.data[key] + obj.data[key]
					else:
						# what to do with tuples / lists such as box size ???
						pass

			module = importlib.import_module(__name__)
			cName = getattr(module, type(self).__name__)

			return cName(sel=None, units=self._units, data=data)
		else:
			raise RuntimeError("Cannot add two objects of different types: {} and {}".format(type(obj), type(self)))

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
					print('Two subsystems with different number of elements cannot be subtracted.')

		return self

	def __mul__(self, obj):
		""" Multiplies scalars/vectors from particle radii/positions
		TODO: get this working with meshes """

		if type(obj) is not type(self):
			raise RuntimeError('Two subsystems with different types cannot be multiplied.')
		
		if len(obj) != len(self): # gotta make sure both classes being multiplied have the same number of elements
			raise RuntimeError('Two subsystems with different number of elements cannot be multiplied.')

		data = collections.OrderedDict()

		for key in self.keys:

			if key in obj.data:
				if isinstance(self.data[key], np.ndarray):
					data[key] = np.sqrt(np.outer(self.data[key], obj.data[key])).flatten()
				else:
					# what to do with tuples / lists such as box size ???
					pass

		module = importlib.import_module(__name__)
		cName = getattr(module, type(self).__name__)

		return cName(sel=None, units=self._units, data=data)

	def __div__(self, obj):
		""" Divides scalars/vectors from particle radii/positions
		TODO: get this working with meshes """

		if type(obj) is tuple:
			obj, att = obj
			if type(obj) is type(self):
				if att is 'all':
					pass
				elif len(obj) == len(self): # gotta make sure both classes being divided have the same number of elements
					self.data[att] /= obj.data[att]
				else:
					print('Two subsystems with different number of elements cannot be divided.')

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

		for i in range(len(self)):
			yield self[i]

	def scale(self, value, attr=('x','y','z')):
		""" Scales all ss elements by a float or int 'value'

		@[attr]: attribute to scale (positions by default)
		"""
		for at in attr:
			if at in self.data:

				if type(self.data[at]) is np.ndarray:
					self.data[at].flags.writeable = True

				self.data[at] *= value

				if type(self.data[at]) is np.ndarray:
					self.data[at].flags.writeable = False

		self._constructAttributes()

	def translate(self, value, attr=('x','y','z')):
		""" Translates all ss elements by a tuple of integers or floats

		@value: tuple of float or int by which to translate the system
		@[attr]: tuple of attributes to translate system by (positions by default)
		"""

		if len(value) != len(attr):
			raise ValueError('The length of values must be equal to that of attributes.')

		for i, at in enumerate(attr):
			if at in self.data:

				if isinstance(self.data[at], np.ndarray):
					self.data[at].flags.writeable = True

				self.data[at] += value[i]

				if isinstance(self.data[at], np.ndarray):
					self.data[at].flags.writeable = False

		self._constructAttributes()

	def noise(self, sigma, attr=('x','y','z')):
		""" Adds white noise of standard deviation `sigma' to all elements with attribute 'attr'.

		@sigma: standard deviation of the Gaussian (white) noise
		@[attr]: attribute to perturb (positions by default)
		"""

		if sigma < 0:
			raise ValueError('Standard deviation must be positive.')

		for at in attr:
			if at in self.data:

				if type(self.data[at]) is np.ndarray:
					self.data[at].flags.writeable = True

					self.data[at] += sigma * (1.0 - 2.0 * np.random.rand(len(self.data[at])))

					self.data[at].flags.writeable = True


		self._constructAttributes()

	def __len__(self):
		""" Returns the number of elements (e.g. particles or nodes) in a SubSystem """
		for key in self.data.keys():
			if type(self.data[key]) == np.ndarray:
				return len(self.data[key])

		# Assume we're dealing with single particle object
		return 1

	def __del__(self):
		pass

class Mesh(SubSystem):
	"""  The Mesh class stores a list of meshes and their associated attributes / methods.
	This class is iterable but NOT an iterator.

	@fname: mesh filename
	"""
	def __init__(self, fname, **args):

		if 'avgCellData' not in args:
			args['avgCellData'] = False

		if 'vtk_type' in args:
			self._vtk = args['vtk_type']
		else:
			self._vtk = None

		if self._vtk == 'poly':
			# Try polydata otherwise unstructured reader ... need to make this better
			# The try-except control makes self._vtk kinda useless, right?
			try:
				self._reader = vtk.vtkPolyDataReader()
			except:
				self._reader = vtk.vtkUnstructuredGridReader()
		else:
			try:
				self._reader = vtk.vtkUnstructuredGridReader()
			except:
				self._reader = vtk.vtkUnstructuredGridReader()			

		super(Mesh, self).__init__(fname=fname, **args)

		# Assert mesh input filname is VTK
		if not hasattr(self, '_mesh'):
			if self._fname:
				if self._fname.split('.')[-1] != 'vtk':
					raise IOError('Input mesh must be of VTK type.')

				self._fname = sorted(glob.glob(self._fname), key=numericalSort)
				self._mesh = self._fname[0]

		self._reader.SetFileName(self._mesh)
		self._reader.ReadAllVectorsOn()
		self._reader.ReadAllScalarsOn()

		self._reader.Update() # Needed if we need to call GetScalarRange
		self._output = self._reader.GetOutput()

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

				newkey = 'cells.' + key
				self.data[newkey] = vtk_to_numpy(self._output.GetCellData().GetArray(key))

				if args['avgCellData']:
					if 'CellArea' in self.data:
						self.data[newkey] = (self.data['CellArea'] * self.data[newkey].T).T.sum(axis=0) / self.data['CellArea'].sum()

					elif 'CellVol' in self.data:
						self.data[newkey] = (self.data['CellVol'] * self.data[newkey].T).T.sum(axis=0) / self.data['CellVol'].sum()
						
				index += 1
			else:
				break

		index = 0
		while True:
			key = self._output.GetPointData().GetArrayName(index)
			if key:
				newkey = 'points.' + key
				self.data[newkey] = vtk_to_numpy(self._output.GetPointData().GetArray(key))
				index += 1
			else:
				break

		self._constructAttributes(mesh=True)
		
	def _updateSystem(self):
		""" Class function for updating the state of a Mesh """
		# Must make sure fname is passed in case we're looping over a trajectory
		self._args['mesh'] = self._mesh
		self._args['units'] = self._units
		self._args['vtk_type'] = self._vtk
		self._args['fname'] = self._fname

		self.__init__(**self._args) # it's imp to pass args so that any subsystem-specific args are
		# passed along to the next frame (e.g. vtk_type, etc.)
		self._constructAttributes()

	def nCells(self):
		return self._output.GetNumberOfCells()

	def nPoints(self):
		return self._output.GetNumberOfPoints()

	def _goto(self, iframe, frame):

		if self._fname:
			if frame >= len(self._fname):
				print('Input frame exceeds max number of frames')
			else:
				if frame == iframe:
					pass
				elif frame == -1:
					frame = len(self._fname) - 1

				self._mesh = self._fname[frame]

		return frame

	def _readFile(self, frame):
		""" Reads a mesh file """

		if 'skip' in self._args:
			skip = self._args['skip']
		else:
			skip = 0

		self._mesh = self._fname[frame + frame * skip]

		return frame + 1

	def _resetSubSystem(self):

		self.__init__(**self._args) 
		super(Mesh, self)._resetSubSystem()

		return 0

	def __del__(self):

		SubSystem.__del__(self)

class Particles(SubSystem):
	""" The Particle class stores all particle properties and the methods that operate on
	these properties. This class is iterable but NOT an iterator. """

	def __init__(self, **args):

		if 'sel' in args:
			sel = args['sel']
		else:
			sel = None

		super(Particles, self).__init__(**args)

		if not hasattr(self, '_fp'): # Make sure a file is already not open
			if hasattr(self, '_fname'):
				if self._fname:
					self._ftype = self._fname.split('.')[-1]

					if self._ftype == 'dump': # need a way to figure out this is a LIGGGHTS/LAMMPS file

						if self._fname.split('.')[:-1][0].endswith('*'):
							self._files = sorted(glob.glob(self._fname), key=numericalSort)
							self._fp = open(self._files[0], 'r')
						else:
							self._fp = open(self._fname, 'r')

						self._params = None

						# Do some checking here on the traj extension to make sure
						# it's supported
						self._format = self._fname.split('.')[-1]

						# Read 1st frame
						self._readFile(0)

					elif self._ftype == 'vtk':
						if self._fname.split('.')[:-1][0].endswith('*'):
							self._files = sorted(glob.glob(self._fname), key=numericalSort)

							for filen in self._files:
								if 'boundingBox' in filen:
									self._files.remove(filen)

							self._fp = open(self._files[0], 'r')
						else:
							self._fp = open(self._fname, 'r')
					else:
						raise IOError('Input trajectory must be a valid LAMMPS/LIGGGHTS (dump) or vtk file.')

		self._constructAttributes(sel)
		self.data['natoms'] = len(self)

		# Make sure natoms is updated ~ DUH
		self._constructAttributes()

		# see if multi-spheres are present
		# TODO: support polydisperse multisphere
		if 'mol' in args:
			data = collections.OrderedDict()
			nmols = int(self.mol.max())
			nspheres = self.natoms / nmols

			for key in self.data.keys():
				if key != 'mol': # elminate recursive call (we dont need mol anyway)
					if type(self.data[key]) == np.ndarray:
						data[key] = np.zeros(nmols)

						for atom, prop in enumerate(self.data[key]):
							data[key][int(self.mol[atom])-1] += prop

						data[key] /= nspheres
						# we average all atomic properties and assign them to the molecules ~ makes no sense to atomic IDs?
					else:
						data[key] = self.data[key] # must be natoms

			data['natoms'] = nmols

			# Compute length of multispheres
			dim = 3 # ~ HACKISH!
			data['length'] = np.zeros((nmols,dim,2))

			for atom in range(self.natoms):
				data['length'][int(self.mol[atom])-1,0,0] = min(self.x[atom], data['length'][int(self.mol[atom])-1,0,0])
				data['length'][int(self.mol[atom])-1,1,0] = min(self.y[atom], data['length'][int(self.mol[atom])-1,1,0])
				data['length'][int(self.mol[atom])-1,2,0] = min(self.z[atom], data['length'][int(self.mol[atom])-1,2,0])

				data['length'][int(self.mol[atom])-1,0,1] = max(self.x[atom], data['length'][int(self.mol[atom])-1,0,1])
				data['length'][int(self.mol[atom])-1,1,1] = max(self.y[atom], data['length'][int(self.mol[atom])-1,1,1])
				data['length'][int(self.mol[atom])-1,2,1] = max(self.z[atom], data['length'][int(self.mol[atom])-1,2,1])

			tmp = np.zeros(nmols)

			for dim in range(dim):
				tmp += (data['length'][:,dim,1] - data['length'][:,dim,0])**2

			data['length'] = np.sqrt(tmp)

			# Make sure molecules can be updated if reading a trajectory, so we delete it
			if hasattr(self, 'molecules'):
				del self.molecules

			self.molecules = Particles(units=self._units, data=data)

	def _updateSystem(self):
		""" Class function for updating the state of Particles """
		# Must make sure fname is passed in case we're looping over a trajectory

		# If we are updating (looping over traj) we wanna keep the id (ref) of the class constant
		# (soft copy)
		self.__init__(sel=None, units=self._units, fname=self._fname, **self.data)
		self._constructAttributes()

	def _resetSubSystem(self):

		self.__init__(**self._args) 
		super(Particles, self)._resetSubSystem()

		return 0

	def computeROG(self):
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

	def computeCOM(self):
		""" Returns center of mass """
		vol = 4.0/3.0 * np.pi * self.radius**3

		return np.array([np.dot(self.x, vol), np.dot(self.y, vol), np.dot(self.z, vol)]) / vol.sum()

	def computeGCOM(self):
		""" Returns the geometric center of mass """
		vol = 4.0/3.0 * np.pi * self.radius**3
		r = len(vol)

		return np.array([(self.x * vol)**r, (self.y, vol)**r, (self.z, vol)**r]).sum(axis=1) / vol.sum()


	def computeRadius(self, N=100):
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

	def computeRDF(self, dr = None, center = True, rMax=None):
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

		if not (self.natoms > 0):
			raise RuntimeError('No Particles found.')

		x, y, z = self.x.copy(), self.y.copy(), self.z.copy()

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

		print('Constructing a cube of length {} and a circumscribed sphere of radius {}'.format(S * 2.0, rMax))
		print('Resolution chosen is {}'.format(dr))

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

		return (g_average / numberDensity, radii, interior_indices)
		# Number of particles in shell/total number of particles/volume of shell/number density
		# shell volume = 4/3*pi(r_outer**3-r_inner**3)

	def computeAngleRepose(self):
		"""
		Computes the angle of repos theta = arctan(h_max/L)
		in a sim box defined by [-Lx, Lx] x [-Ly, Ly] x [0, Lz]
		"""
		x, y, z = self.x, self.y, self.z
		dL = 0.25 * (x.max() - x.min()) + 0.25 * (y.max() - y.min())
		z_max = z.max() - z.min() - self.radius.mean()

		return np.arctan(z_max / dL) * 180.0 / np.pi
		
	def computeMass(self, tdensity):
		""" Computes the mass of all particles 

		@tdensity: true density of the powder
		returns the summation of the mass of all particles """

		if self.natoms > 0:

			radius = self.radius
			mass = tdensity * 4.0 / 3.0 * np.pi * (radius**3.0)

			return mass
		else:
			return None

	def computeIntensitySegregation(self, resol=None):
		""" Computes the intensity of segregation for binary mixture
		as defined by Danckwerts:

		I = sigma_a**2 / (mean_a (1 - mean_a))

		@[resol]: bin size for grid construction (default 3 * diameter)
		returns scalar intensity 
		"""

		if not resol:
			resol = self.radius.min() * 3

		if len(np.unique(self.type)) != 2:
			raise ValueError("Intensity of segergation can be computed only for a binary system.")
			# should I support tertiary systems or more ???

		nx = int((self.x.max() - self.x.min()) / resol) + 1
		ny = int((self.y.max() - self.y.min()) / resol) + 1
		nz = int((self.z.max() - self.z.min()) / resol) + 1

		indices = np.zeros((nx,ny,nz), dtype='float64')

		for sn, ctype in enumerate(np.unique(self.type)):

			parts = self[self.type==ctype]

			parts.translate(value=(-parts.x.min(), -parts.y.min(), -parts.z.min()), attr=('x','y','z'))

			x = np.array(parts.x / resol, dtype='int64')
			y = np.array(parts.y / resol, dtype='int64')
			z = np.array(parts.z / resol, dtype='int64')

			for i in range(parts.natoms):
				indices[x[i],y[i],z[i]] += parts.radius[i]**3

			if sn == 0:
				indices_a = indices.copy()

		indices_a[indices > 0] /= indices[indices > 0]
		aMean = indices_a[indices > 0].mean()
		aStd = indices_a[indices > 0].std()

		return aStd**2 / (aMean * (1.0 - aMean)), indices_a, indices

	def computeScaleSegregation(self, nTrials=1000, resol=None, Npts=50, maxDist=None):
		""" Computes the correlation coefficient as defined by Danckwerts:
		R(r) = a * b / std(a)**2

		This is done via a Monte Carlo simulation.

		@[resol]: bin size for grid construction (default min radius)
		@[nTrials]: number of Monte Carlo trials (sample size)
		@[Npts]: number of bins for histogram construction
		@[maxDist]: maximum distance (in units of grid size) to sample

		Returns the coefficient of correlation R(r) and separation distance (r)
		"""

		if not resol:
			resol = self.radius.min()
	
		_, a, total = self.intensitySegregation(resol)

		if not maxDist:
			maxDim = max(a.shape)
			maxDist = int(np.sqrt(3 * maxDim**2)) + 1
		
		volMean = a[total > 0].mean()
		volVar = a[total > 0].std()**2

		corr = np.zeros(nTrials)
		dist = np.zeros(nTrials)
		count = 0

		# Begin Monte Carlo simulation
		while count < nTrials:

			i1, i2 = np.random.randint(0, a.shape[0], size=2)
			j1, j2 = np.random.randint(0, a.shape[1], size=2)
			k1, k2 = np.random.randint(0, a.shape[2], size=2)

			# Make sure we are sampling non-void spatial points
			if total[i1, j1, k1] > 0 and total[i2, j2, k2] > 0:

				distance = np.sqrt((i2 - i1)**2 + (j2 - j1)**2 + (k2 - k1)**2)

				if distance <= maxDist:

					corr[count] = ((a[i1,j1,k1] - volMean) * (a[i2,j2,k2] - volMean)) / volVar
					dist[count] = distance
					count += 1

		corrfunc, distance, _ = binned_statistic(dist, corr, 'mean', Npts)

		return corrfunc[np.invert(np.isnan(corrfunc))], distance[0:-1][np.invert(np.isnan(corrfunc))] * resol

	def computeDensity(self, tdensity, shape = 'box', bounds=None):
		"""
		Computes the bulk density for a selection of particles from their *true* density.
		The volume is determined approximately by constructing a box/cylinder/cone
		embedding the particles. Particles are assumed to be spherical in shape.

		@tdensity: true powder density
		@[shape]: box, cylinder-x, cylinder-y, or cylinder-z

		"""
		
		return self.computeMass(tdensity).sum() / self.computeVolume(shape)

	def computeDensityLocal(self, bdensity, dr, axis):
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

	def computeVolume(self, shape = 'box'):
		""" Computes the volume of a granular system based on a simple geometry 
		@[shape]: box, cylinder-x, cylinder-y, or cylinder-z"""

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

	def _goto(self, iframe, frame):
		""" This function assumes we're reading a non-const N trajectory.
		"""

		# find the right frame number
		if self._singleFile:
			# We must be reading a single particle trajectory file (i.e. not a mesh)
			while frame < iframe or iframe == -1:

				line = self._fp.readline()

				if not line and iframe >= 0:
					raise StopIteration('End of file reached.')
				elif not line and iframe == -1:
					break

				if line.find('TIMESTEP') >= 0:
					frame += 1

			# assert self.frame == frame else something's wrong
			# or the user wants to go to the last frame
			if iframe == frame:

				ts = int(self._fp.readline())
				self.data['timestep'] = ts
				self._constructAttributes()

				self._readFile(frame)

			elif iframe == -1:
				tmp = frame
				frame = 0

				del self._fp
				self._readFile(0)
				return self._goto(tmp, 1)

			else:
				raise NameError('Cannot find frame {} in current trajectory'.format(frame))

		else: # no need to find the input frame, just select the right file if available

			if self._fp:
				if iframe >= len(self._files):
					print('Input frame exceeds max number of frames')
				else:
					if frame == iframe:
						pass
					else:
						if iframe == -1:
							iframe = len(self._files) - 1

						self._fp.close()
						self._fp = open(self._files[iframe], 'r')
						frame = iframe
						self._readDumpFile()

		return frame

	def _readDumpFile(self):
		""" Reads a single dump file"""
		count = 0
		ts = None

		while True:
			line = self._fp.readline()

			if not line:
				raise StopIteration

			if line.find('TIMESTEP') >= 0:
				ts = int(self._fp.readline())
				self.data['timestep'] = ts

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

		return ts

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

	def _readFile(self, frame):
		""" Read a particle trajectory file 

		TODO: Support skip for single dump file
		"""

		if 'skip' in self._args:
			skip = self._args['skip']
		else:
			skip = 0

		# We are opening the traj file for the 1st time
		if not hasattr(self, '_fp'):
			self._fp = open(self._fname, 'r')

		if not hasattr(self, '_singleFile'):

			if self._fname.split('.')[:-1][0].endswith('*'):
				self._singleFile = False
			else:
				self._singleFile = True

		if self._singleFile: # must support skip
			if self._ftype == 'dump': 
				ts = self._readDumpFile()
				self.data['timestep'] = ts
				self._constructAttributes()

			else:
				raise IOError('{} format is not a supported trajectory file.'.format(self._ftype))

		else:
			if self._ftype == 'dump':

				if frame > len(self._files) - 1:
					raise StopIteration

				self._fp.close()
				self._fname = self._files[frame + skip * frame]
				self._fp = open(self._fname, 'r')
				ts = self._readDumpFile()
				self._constructAttributes()

			elif self._ftype == 'vtk':

				if frame > len(self._files) - 1:
					raise StopIteration

				self._fp.close()
				self._fname = self._files[frame + skip * frame]

				self._reader = vtk.vtkPolyDataReader()

				self._reader.SetFileName(self._fname)
				self._reader.Update() # Needed if we need to call GetScalarRange
				self._output = self._reader.GetOutput()

				pos = vtk_to_numpy(self._output.GetPoints().GetData())
				self.data['x'] = pos[:,0]

				if pos.shape[1]  >= 1:
					self.data['y'] = pos[:,1]

				if pos.shape[1] >= 2:
					self.data['z'] = pos[:,2]

				index = 0
				while True:
					key = self._output.GetCellData().GetArrayName(index)
					if key:
						self.data[key] = vtk_to_numpy(self._output.GetCellData().GetArray(key))
						index += 1
					else:
						break

				index = 0
				while True:
					key = self._output.GetPointData().GetArrayName(index)
					if key:
						self.data[key] = vtk_to_numpy(self._output.GetPointData().GetArray(key))
						index += 1
					else:
						break

				# This doesnot work for 2D / 1D systems
				for key in self.data.keys():
					if isinstance(self.data[key], np.ndarray):
						if len(self.data[key].shape) > 1:
							if self.data[key].shape[1]  == 3:
								self.data[key + 'x'] = self.data[key][:,0]
								self.data[key + 'y'] = self.data[key][:,1]
								self.data[key + 'z'] = self.data[key][:,2]

								del self.data[key]

				self._constructAttributes()

		return frame + 1

class Factory(object):
	"""A factory for system class. It creates subclasses of SubSystems. Its only two methods
	are static, thus no need to instantiate this class.

	@SubSystem: filename (or list of filenames) for the subsystem trajectory
	@[units](si): unit system can be either 'si' or 'micro'
	"""

	def factory(**args):
		""" Returns a list of Subsystems 

		System(Obj1=[(value11, args11), (value12, args12), ... (value1N, args1N)], ...)

		args11, args12, etc. are dictionaries of keywords args for objects of type obj1.
		"""
		obj = []

		if 'module' in args:
			module = args['module']
		else:
			module = None

		for ss in args:

			sclass = Factory._str_to_class(ss, module=module)

			if not sclass:
				# make sure selected class definition exists, otherwise, this had better be a list
				raise ValueError('System takes only keywords of classes that are defined. Class type {} not found.'.format(args[ss]))

			if isinstance(args[ss], list) or isinstance(args[ss], tuple): # we need to create a list of objects

				objs = []

				for item in args[ss]:

					if isinstance(item, tuple): # we must pass additional args

						try:
							obj_fname, obj_args = item
						except:
							raise ValueError('SubSystem arguments must be passed as a dictionary, i.e. System(SubSystem=[(path1, arg1), ...]) ')
						else:
							objs.append(sclass(fname=obj_fname, **obj_args))

					elif isinstance(item, str):
						objs.append(sclass(fname=item))
					elif isinstance(item, sclass):
						objs.append(item)
					else:
						raise ValueError('Incorrect keyarg supplied to System: {} when creating {} SubSystem'.format(item, ss))

				if objs:
					obj.append([ss, objs])

			elif isinstance(args[ss], str):

				fname = args[ss]
				obj.append([ss, sclass(fname=fname)])

			elif isinstance(args[ss], tuple):

				obj_fname, obj_args = args[ss]
				obj.append([ss, sclass(fname=obj_fname, **obj_args)])

			elif isinstance(args[ss], sclass):
				obj.append([ss, args[ss]])
			else:
				raise ValueError('Incorrect keyarg supplied to System: {} when creating {} SubSystem'.format(args[ss], ss))
			
		return obj

	def _str_to_class(string, module=None):

		# See if the object exists within the current module
		try:
			if module:
				module = importlib.import_module(module)
				return getattr(module, string)
			else:
				return getattr(sys.modules[__name__], string)
		except:
			pass # string has no class definition in the specified module

	def addprop(inst, name, method):

		cls = type(inst)

		if not hasattr(cls, '__perinstance'):
			cls = type(cls.__name__, (cls,), {})
			cls.__perinstance = True
			inst.__class__ = cls

		setattr(cls, name, property(method))

	factory = staticmethod(factory)
	addprop = staticmethod(addprop)

	_str_to_class = staticmethod(_str_to_class)


class System(object):
	"""A System contains all the information describing a DEM system.
	A meaningful system always requires at least one trajectory file to read. A trajectory is a (time)
	series corresponding to the coordinates of all particles/meshes/objects in the system.
	System handles the time frame and controls i/o operations. It contains one or more SubSystem
	derivatives (e.g. 'Particles', 'Mesh') which store all the particle and mesh attributes
	read from the trajectory file (variables such as positions, momenta, angular velocities,
	forces, stresses, radii, etc.). Objects that can be created by this class must be of the 
	'Particles' or 'Mesh' type. Multiple objects can be created by this class if a list of 
	filenames are passed to its constructors.

	@[Particles]: filename (or list of filenames) for the particle trajectory
	@[Mesh]: filename (or list of filenames) for the mesh trajectory
	@[units](si): unit system can be either 'si' or 'micro'
	@[module]: for user-defined dervied SubSystems, 'module' must be a string that specifies the name of 
	the module in which the user-defined class is defined.

	How time stepping works: when looping over System (traj file), the frame is controlled only by System
	through methods defined in a SubSystem sublass (read/write functions).

	This class is an iterator.
	"""

	def __init__(self, **args):

		self.frame = 0
		self.args = args
		objs = Factory.factory(**args)

		for ss, obj in objs:
			setattr(self, ss, obj)

		if 'units' in args:
			self.units(args['units'])
		else:
			self.units('si')

	def __iter__(self):
		self.rewind()
		return self

	def units(self, units=None):
		""" Change unit system to units ore return units if None
		units:  si, micro, cgs, or nano """

		if not units:
			return self._units

		if (units == 'si') or (units == 'micro'):
			for ss in self.__dict__:
				if hasattr(self.__dict__[ss], 'units'):
					self.__dict__[ss].units(units)

			self._units = units
		else:
			print('Only S.I. and micro units currently supported')

	def goto(self, frame):
		""" Go to a specific frame in the trajectory. If frame is -1
		then this function will read the last frame """

		# nothing to do
		if frame == self.frame:
			return 0
		else:
			newFrame = None

		# rewind if necessary (better than reading file backwards?)
		if frame < self.frame and frame > 0:
			self.rewind()
			self.goto(frame)
		elif frame == 0: # rewind
			for ss in self.__dict__:
				if isinstance(self.__dict__[ss], list):
					for item in self.__dict__[ss]:
						if hasattr(item, '_resetSubSystem'):
							self.frame  = item._resetSubSystem()

				elif hasattr(self.__dict__[ss], '_resetSubSystem'):
					self.frame  = self.__dict__[ss]._resetSubSystem()

		else:
			for ss in self.__dict__:
				if isinstance(self.__dict__[ss], list):
					for item in self.__dict__[ss]:
						if hasattr(item, '_goto'):
							newFrame = item._goto(self.frame, frame)

				elif hasattr(self.__dict__[ss], '_goto'):
					newFrame = self.__dict__[ss]._goto(self.frame, frame)

			if newFrame:
				self.frame = newFrame

			# Rewind already updates the system, so we call _updateSystem only if
			# the frame is moving forward
			self._updateSystem()

	def skip(self):
		""" Skips all empty frames i.e. moves the trajectory to the 1st frame containing
		non-zero elements """
		forward = True

		while forward:
			for ss in self.__dict__: 
				if isinstance(self.__dict__[ss], list):
					for item in self.__dict__[ss]:
						if hasattr(item, '_constructAttributes'):
							if len(item):
								forward = False
				elif hasattr(self.__dict__[ss], '_constructAttributes'):
					if len(self.__dict__[ss]):
						forward = False

			if not forward:
				break
			else:
				self.frame = self.next()
	@property
	def keys(self):
		""" returns the key objects found in the trajctory files """
		return self.data.keys()

	def rewind(self):
		"""Read trajectory from the beginning"""

		self.goto(0)

	def __next__(self):
		"""Forward one step to next frame when using the next builtin function."""
		return self.next()

	def next(self):
		""" This method updates the system attributes! """
		update_frame = False

		try:
			for ss in self.__dict__:

				if 'module' in self.args:
					module = self.args['module']
				else:
					module = None

				if Factory._str_to_class(ss, module):
					if isinstance(self.__dict__[ss], list):

						for obj in self.__dict__[ss]:
							if hasattr(obj, '_readFile'):
								obj._readFile(self.frame + 1)
								update_frame = True

					elif hasattr(self.__dict__[ss], '_readFile'):
						self.__dict__[ss]._readFile(self.frame + 1)
						update_frame = True
		except:
			self.rewind()
			raise StopIteration
		else:
			# we update the fame only after _readFile is called. If the latter
			# fails, the frame is not updated (due to raise), which is exactly
			# the behavior we want. Ding!
			if update_frame:
				self.frame += 1 
				self._updateSystem()

		return self.frame

	def _updateSystem(self):
		""" Makes sure the system is aware of any update in its attributes caused by
		a frame change. """

		# A generic way of invoking _updateSystem

		if 'module' in self.args:
			module = self.args['module']
		else:
			module = None

		for ss in self.__dict__:
			if Factory._str_to_class(ss, module): # Make sure the obj exists
				if isinstance(self.__dict__[ss], list):
					for obj in self.__dict__[ss]:
						# Make sure this is a list of defined SubSystems
						if hasattr(obj, '_updateSystem'):
							obj._updateSystem()

				elif hasattr(self.__dict__[ss], '_updateSystem'):
					self.__dict__[ss]._updateSystem()

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
				print('Length of region must be 6: (xmin, xmax, ymin, ymax, zmin, zmax)')
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
