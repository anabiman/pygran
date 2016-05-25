#   Author: Andrew Abi-Mansour
#   Python wrapper for LIGGGHTS library via ctypes
#
# ----------------------------------------------------------------------
#
#   Modified from  LAMMPS source code
#   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
#   http://lammps.sandia.gov, Sandia National Laboratories
#   Steve Plimpton, sjplimp@sandia.gov
#
#   Copyright (2003) Sandia Corporation.  Under the terms of Contract
#   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
#   certain rights in this software.  This software is distributed under 
#   the GNU General Public License.
#
#   See the README file in the top-level LAMMPS directory.
#
# -------------------------------------------------------------------------

import sys,traceback,types
from ctypes import *
from os.path import dirname, abspath, join
from inspect import getsourcefile

import numpy as np
import itertools
from numpy.linalg import norm
from scipy import spatial
from mpi4py import MPI
import os
import matplotlib.pylab as plt
import glob
import sys
from importlib import import_module

def find(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)

class liggghts:
  # detect if Python is using version of mpi4py that can pass a communicator
  
  has_mpi4py_v2 = False
  try:
    from mpi4py import MPI
    from mpi4py import __version__ as mpi4py_version
    if mpi4py_version.split('.')[0] == '2':
      has_mpi4py_v2 = True
  except:
    pass

  # create instance of LIGGGHTS
 
  def __init__(self,name="libliggghts.so",cmdargs=None,ptr=None,comm=None):

    # determine module location
    print "Looking for {} ...".format(name)
    foundliggghts = find(name, "/")
 
    if foundliggghts:
      print "Using " + foundliggghts
    else:
      print "Make sure a " + name + " is installed on your system"
      sys.exit()

    self.lib = CDLL(foundliggghts, RTLD_GLOBAL)

    # if no ptr provided, create an instance of LIGGGHTS
    #   don't know how to pass an MPI communicator from PyPar
    #   but we can pass an MPI communicator from mpi4py v2.0.0 and later
    #   no_mpi call lets LIGGGHTS use MPI_COMM_WORLD
    #   cargs = array of C strings from args
    # if ptr, then are embedding Python in LIGGGHTS input script
    #   ptr is the desired instance of LIGGGHTS
    #   just convert it to ctypes ptr and store in self.lmp
    
    if not ptr:
      # with mpi4py v2, can pass MPI communicator to LIGGGHTS
      # need to adjust for type of MPI communicator object
      # allow for int (like MPICH) or void* (like OpenMPI)
      
      if liggghts.has_mpi4py_v2 and comm != None:
        if liggghts.MPI._sizeof(liggghts.MPI.Comm) == sizeof(c_int):
          MPI_Comm = c_int
        else:
          MPI_Comm = c_void_p

        narg = 0
        cargs = 0
        if cmdargs:
          cmdargs.insert(0,"liggghts.py")
          narg = len(cmdargs)
          cargs = (c_char_p*narg)(*cmdargs)
          self.lib.lammps_open.argtypes = [c_int, c_char_p*narg, \
                                           MPI_Comm, c_void_p()]
        else:
          self.lib.lammps_open.argtypes = [c_int, c_int, \
                                           MPI_Comm, c_void_p()]

        self.lib.lammps_open.restype = None
        self.opened = 1
        self.lmp = c_void_p()
        comm_ptr = liggghts.MPI._addressof(comm)
        comm_val = MPI_Comm.from_address(comm_ptr)
        self.lib.lammps_open(narg,cargs,comm_val,byref(self.lmp))

      else:
        self.opened = 1
        if cmdargs:
          cmdargs.insert(0,"liggghts.py")
          narg = len(cmdargs)
          cargs = (c_char_p*narg)(*cmdargs)
          self.lmp = c_void_p()
          self.lib.lammps_open_no_mpi(narg,cargs,byref(self.lmp))
        else:
          self.lmp = c_void_p()
          self.lib.lammps_open_no_mpi(0,None,byref(self.lmp))
          # could use just this if LIGGGHTS lib interface supported it
          # self.lmp = self.lib.lammps_open_no_mpi(0,None)
          
    else:
      self.opened = 0
      # magic to convert ptr to ctypes ptr
      pythonapi.PyCObject_AsVoidPtr.restype = c_void_p
      pythonapi.PyCObject_AsVoidPtr.argtypes = [py_object]
      self.lmp = c_void_p(pythonapi.PyCObject_AsVoidPtr(ptr))

  def __del__(self):
    if self.lmp and self.opened: self.lib.lammps_close(self.lmp)

  def close(self):
    if self.opened: self.lib.lammps_close(self.lmp)
    self.lmp = None

  def version(self):
    return self.lib.lammps_version(self.lmp)

  def file(self,file):
    self.lib.lammps_file(self.lmp,file)

  def command(self,cmd):
    self.lib.lammps_command(self.lmp,cmd)

  def extract_global(self,name,type):
    if type == 0:
      self.lib.lammps_extract_global.restype = POINTER(c_int)
    elif type == 1:
      self.lib.lammps_extract_global.restype = POINTER(c_double)
    else: return None
    ptr = self.lib.lammps_extract_global(self.lmp,name)
    return ptr[0]

  def extract_atom(self,name,type):
    if type == 0:
      self.lib.lammps_extract_atom.restype = POINTER(c_int)
    elif type == 1:
      self.lib.lammps_extract_atom.restype = POINTER(POINTER(c_int))
    elif type == 2:
      self.lib.lammps_extract_atom.restype = POINTER(c_double)
    elif type == 3:
      self.lib.lammps_extract_atom.restype = POINTER(POINTER(c_double))
    else: return None
    ptr = self.lib.lammps_extract_atom(self.lmp,name)
    return ptr

  def extract_compute(self,id,style,type):
    if type == 0:
      if style > 0: return None
      self.lib.lammps_extract_compute.restype = POINTER(c_double)
      ptr = self.lib.lammps_extract_compute(self.lmp,id,style,type)
      return ptr[0]
    if type == 1:
      self.lib.lammps_extract_compute.restype = POINTER(c_double)
      ptr = self.lib.lammps_extract_compute(self.lmp,id,style,type)
      return ptr
    if type == 2:
      self.lib.lammps_extract_compute.restype = POINTER(POINTER(c_double))
      ptr = self.lib.lammps_extract_compute(self.lmp,id,style,type)
      return ptr
    return None

  # in case of global datum, free memory for 1 double via lammps_free()
  # double was allocated by library interface function
  
  def extract_fix(self,id,style,type,i=0,j=0):
    if style == 0:
      self.lib.lammps_extract_fix.restype = POINTER(c_double)
      ptr = self.lib.lammps_extract_fix(self.lmp,id,style,type,i,j)
      result = ptr[0]
      self.lib.lammps_free(ptr)
      return result
    elif (style == 1) or (style == 2):
      if type == 1:
        self.lib.lammps_extract_fix.restype = POINTER(c_double)
      elif type == 2:
        self.lib.lammps_extract_fix.restype = POINTER(POINTER(c_double))
      else:
        return None
      ptr = self.lib.lammps_extract_fix(self.lmp,id,style,type,i,j)
      return ptr
    else:
      return None

  # free memory for 1 double or 1 vector of doubles via lammps_free()
  # for vector, must copy nlocal returned values to local c_double vector
  # memory was allocated by library interface function
  
  def extract_variable(self,name,group,type):
    if type == 0:
      self.lib.lammps_extract_variable.restype = POINTER(c_double)
      ptr = self.lib.lammps_extract_variable(self.lmp,name,group)
      result = ptr[0]
      self.lib.lammps_free(ptr)
      return result
    if type == 1:
      self.lib.lammps_extract_global.restype = POINTER(c_int)
      nlocalptr = self.lib.lammps_extract_global(self.lmp,"nlocal")
      nlocal = nlocalptr[0]
      result = (c_double*nlocal)()
      self.lib.lammps_extract_variable.restype = POINTER(c_double)
      ptr = self.lib.lammps_extract_variable(self.lmp,name,group)
      for i in xrange(nlocal): result[i] = ptr[i]
      self.lib.lammps_free(ptr)
      return result
    return None

  # set variable value
  # value is converted to string
  # returns 0 for success, -1 if failed
  
  def set_variable(self,name,value):
    return self.lib.lammps_set_variable(self.lmp,name,str(value))

  # return total number of atoms in system
  
  def get_natoms(self):
    return self.lib.lammps_get_natoms(self.lmp)

  # return vector of atom properties gathered across procs, ordered by atom ID

  def gather_atoms(self,name,type,count):
    natoms = self.lib.lammps_get_natoms(self.lmp)
    if type == 0:
      data = ((count*natoms)*c_int)()
      self.lib.lammps_gather_atoms(self.lmp,name,type,count,data)
    elif type == 1:
      data = ((count*natoms)*c_double)()
      self.lib.lammps_gather_atoms(self.lmp,name,type,count,data)
    else: return None
    return data

  # scatter vector of atom properties across procs, ordered by atom ID
  # assume vector is of correct type and length, as created by gather_atoms()

  def scatter_atoms(self,name,type,count,data):
    self.lib.lammps_scatter_atoms(self.lmp,name,type,count,data)

class DEMPy:
  """A class that implements a python interface for DEM computations"""

  def __init__(self, sid, split, units, dim, style, **pargs):
    """ Initialize some settings and specifications 
    @ units: unit system (si, cgs, etc.)
    @ dim: dimensions of the problem (2 or 3)
    # style: granular, atom, or ...
    """
    self.rank = split.Get_rank()

    self.pargs = pargs
    self.monitorList = []
    self.vars = {}
    self.path = os.getcwd()
    self.nSS = len(self.pargs['SS'])

    self.output = self.pargs['output'] if self.pargs['nSim'] == 1 else (self.pargs['output'] + '{}'.format(sid))

    if not self.rank:

      if not os.path.exists(self.output):
        os.makedirs(self.output)

      global logging

      os.chdir(self.output)
      logging = import_module(name='logging')

      logging.basicConfig(filename='dem.log', format='%(asctime)s:%(levelname)s: %(message)s', level=logging.DEBUG)

      logging.info('Creating i/o directories')

      if not os.path.exists(self.pargs['traj']['dir']):
        os.makedirs(self.pargs['traj']['dir'])

      if not os.path.exists(self.pargs['restart'][1]):
        os.makedirs(self.pargs['restart'][1])

      logging.info('Instantiating LIGGGHTS object')

    self.lmp = liggghts(comm=split, cmdargs=['-log', 'liggghts.log'])

    if not self.rank:
      logging.info('Setting up problem dimensions and boundaries')

    self.lmp.command('units {}'.format(units))
    self.lmp.command('dimension {}'.format(dim))
    self.lmp.command('atom_style {}'.format(style))
    self.lmp.command('atom_modify map array') # array is faster than hash in looking up atomic IDs, but the former takes more memory
    self.lmp.command('boundary {} {} {}'.format(*pargs['boundary']))
    self.lmp.command('newton off') # turn off newton's 3rd law ~ should lead to better scalability
    self.lmp.command('communicate single vel yes') # have no idea what this does, but it's imp for ghost atoms
    self.lmp.command('processors * * *') # let LIGGGHTS handle DD

  def createDomain(self):
    """ Define the domain of the simulation
    @ nsys: number of subsystems
    @ pos: 6 x 1 tuple that defines the boundaries of the box 
    """
    if not self.rank:
      logging.info('Creating domain')

    self.lmp.command('region domain block {} {} {} {} {} {} units box'.format(*self.pargs['box']))
    self.lmp.command('create_box {} domain'.format(self.pargs['nSS']))

  def insertParticles(self):
    """ Create atoms in a pre-defined region
    @ N: max total number of particles to be inserted
    @ density: initial density of the particles
    @ vel: 3 x 1 tuple of initial velocities of all particles
    @ args: dictionary of params
    """
    if not self.rank:
      logging.info('Inserting particles')

    for i, ss in enumerate(self.pargs['SS']):

      if ss['insert'] == True:
        radius = self.pargs['radius'][i]

        self.lmp.command('fix pts all particletemplate/sphere 1 atom_type {id} density constant {density} radius'.format(**ss) + (' {}' * len(radius)).format(*radius))
        self.lmp.command('fix pdd all particledistribution/discrete 63243 1 pts 1.0')

        self.lmp.command('region factory ' + ('{} ' * len(ss['region'])).format(*ss['region']) + 'units box')
        self.lmp.command('fix ins all insert/rate/region seed 123481 distributiontemplate pdd nparticles {natoms} particlerate {rate} insert_every {freq} overlapcheck yes vel constant'.format(**ss) + ' {} {} {} region factory ntry_mc 1000'.format(*self.pargs['vel'][i] ))

  def importMesh(self, name, file, scale = None):
    """
    TODO: fix type for mesh
    """
    fname = self.path + '/' + file

    if not self.rank:
      logging.info('Importing mesh from {}'.format(fname))

    if scale == None:
	self.lmp.command('fix {} all mesh/surface file {} type 2'.format(name, fname))
    else:
    	self.lmp.command('fix {} all mesh/surface file {} type 2 scale {}'.format(name, fname, scale))

  def setupWall(self, name, wtype, plane = None, peq = None):
    """
    Creates a wall
    @ name: name of the variable defining a wall or a mesh
    @ wtype: type of the wall (primitive or mesh)
    @ plane: x, y, or z plane for primitive walls
    @ peq: plane equation for primitive walls
    """

    if wtype == 'mesh':
      self.lmp.command('fix myMesh all wall/gran model hooke {} n_meshes 1 meshes {}'.format(wtype, name))
    elif wtype == 'primitive':
      self.lmp.command('fix {} all wall/gran model hooke {} type 1 {} {}'.format(name, wtype, plane, peq))
    else:
      raise ValueError('Wall type can be either primitive or mesh')
 
  def remove(self, name):
    """
    Deletes a specified variable
    """
    self.lmp.command('unfix {}'.format(name))

  def createGroup(self, group = None):
    """ Create groups of atoms 
    """
    if not self.rank:
      logging.info('Creating atom group {}'.format(group))

    if group == None:
      for idSS in self.pargs['idSS']:
        self.lmp.command('group group{} type {}'.format(idSS, idSS))

  def setupNeighbor(self):
    """
    """
    if not self.rank:
      logging.info('Setting up nearest neighbor searching parameters')

    self.lmp.command('neighbor 0.001 bin')
    self.lmp.command('neigh_modify delay 0')

  def createProperty(self, name, *args):
    """
    Material and interaction properties required
    """
    if not self.rank:
      logging.info('Creating property {} with args'.format(name) + (' {}' * len(args)).format(*args))

    self.lmp.command('fix {} all property/global'.format(name) + (' {}' * len(args)).format(*args))

  def setupPhysics(self):
    """
    Specify the interation forces
    """
    if not self.rank:
      logging.info('Setting up interaction parameters')

    args = self.pargs['model']

    self.lmp.command('pair_style ' + (' {}' * len(args)).format(*args))
    self.lmp.command('pair_coeff * *')

  def setupGravity(self):
    """
    Specify in which direction the gravitational force acts
    """
    self.lmp.command('fix myGravity all gravity {} vector {} {} {}'.format(*self.pargs['gravity']))

  def initialize(self):
    """
    """

    self.lmp.command('restart {} {}/{}'.format(*self.pargs['restart']))

    if self.pargs['restart'][-1] == False:

      self.createDomain()
      #self.createGroup()
      self.setupPhysics()
      self.setupNeighbor()
      self.insertParticles()
      self.setupGravity()

    else:
      self.resume()
      self.setupPhysics()
      self.setupNeighbor()
      self.setupGravity()

  def setupIntegrate(self, name):
    """
    Specify how Newton's eqs are integrated in time. 
    @ name: name of the fixed simulation ensemble applied to all atoms
    @ dt: timestep
    @ ensemble: ensemble type (nvt, nve, or npt)
    @ args: tuple args for npt or nvt simulations
    """
    if not self.rank:
      logging.info('Setting up integration scheme parameters')

    self.lmp.command('fix {} all nve/sphere'.format(name))

  def integrate(self, steps, dt = None):
    """
    Run simulation in time
    """
    if not self.rank:
      logging.info('Integrating the system for {} steps'.format(steps))

    for tup in self.monitorList:
      self.lmp.command('compute {} {} {}'.format(*tup))

    if dt is not None:
      self.lmp.command('timestep {}'.format(dt))

    self.lmp.command('run {}'.format(steps))

  def printSetup(self):
    """
    Specify which variables to write to file, and their format
    """
    if not self.rank:
      logging.info('Setting up printing options')

    freq, args = self.pargs['print'][0], self.pargs['print'][1:]

    self.lmp.command('thermo_style custom' + (' {}' * len(args)).format(*args))
    self.lmp.command('thermo {}'.format(freq))
    self.lmp.command('thermo_modify norm no lost ignore')

  def dumpSetup(self):
    """
    """
    if not self.rank:
      logging.info('Setting up trajectory i/o')

    self.lmp.command('dump dump {sel} custom {freq} {dir}/{file}'.format(**self.pargs['traj']) + (' {} ' * len(self.pargs['traj']['args'])).format(*self.pargs['traj']['args']))

  def extractCoords(self, coords):
    """
    Extracts atomic positions from a certian frame and adds it to coords
    """
    # Extract coordinates from liggghts
    self.lmp.command('variable x atom x')
    x = Rxn.lmp.extract_variable("x", "group1", 1)

    self.lmp.command('variable y atom y')
    y = Rxn.lmp.extract_variable("y", "group1", 1)

    self.lmp.command('variable z atom z')
    z = Rxn.lmp.extract_variable("z", "group1", 1)

    for i in range(Rxn.lmp.get_natoms()):
      coords[i,:] += x[i], y[i], z[i]

    self.lmp.command('variable x delete')
    self.lmp.command('variable y delete')
    self.lmp.command('variable z delete')

    return coords

  def monitor(self, name, group, var, file):
    """
    """
    self.lmp.command('compute {} all {}'.format(var, name))
    self.lmp.command('fix my{} {} ave/time 1 1 1 c_{} file {}'.format(var, group, var, file))

  def plot(self, fname, xlabel, ylabel, output=None, xscale=None):
    """
    """
    if not self.rank:
      try:
     	#plt.rc('text', usetex=True)
     	data = np.loadtxt(fname, comments='#')

	time = data[:,0]

	if xscale is not None:
		time *= xscale

      	plt.plot(time, data[:,1])
     	plt.xlabel(r"{}".format(xlabel))
      	plt.ylabel(ylabel)
	plt.grid()

      	if output:
          plt.savefig(output)
      except:
	print "Unexpected error:", sys.exc_info()[0]
	raise	

  def saveas(self, name, fname):
    """
    """
    if not self.rank:

      try:
      	np.savetxt(fname, np.array(self.vars[name]))
      except:
	print "Unexpected error:", sys.exc_info()[0]
        raise

  def resume(self):
    """
    """
    rdir = '{}/*'.format(self.pargs['restart'][1])
    rfile = max(glob.iglob(rdir), key=os.path.getctime)

    self.lmp.command('read_restart {}'.format(rfile))

  def __del__(self):
    """ Destructor
    """
    self.lmp.close()

class DEM:
  """A class that handles communication for the DEM object"""

  def __init__(self, **pargs):
    """ Initialize COMM and partition proccesors based on user input """

    self.comm = MPI.COMM_WORLD
    self.rank = self.comm.Get_rank()
    self.tProcs = self.comm.Get_size()
    self.nSim = pargs['nSim']

    if self.nSim > self.tProcs:
      print "Number of simulations ({}) cannot exceed number of available processors ({})".format(self.nSim, self.tProcs)
      sys.exit(0)

    self.nPart = self.tProcs // self.nSim

    for i in range(self.nSim):
      if self.rank < self.nPart * (i + 1):
        self.color = i
        self.split = self.comm.Split(self.color, key=0)

        self.dem = DEMPy(i, self.split, **pargs) # logging module imported here      
        break

    if not self.rank:
      logging.info("Initializing MPI for a total of {} procs".format(self.tProcs))

      if self.nSim > 1:
        logging.info('Running {} simulations: multi-mode on'.format(self.nSim))

  def initialize(self):

    for i in range(self.nSim):
      if self.rank < self.nPart * (i + 1):
        self.dem.initialize()
        break
   
  def createProperty(self, name, *args):
    """
    Material and interaction properties required
    """
    for i in range(self.nSim):
      if self.rank < self.nPart * (i + 1):
        if type(args[0]) is tuple:
          self.dem.createProperty(name, *args[i])
        else:
          self.dem.createProperty(name, *args)
        break

  def importMesh(self, name, file, scale = None):
    """
    Imports a mesh file (STL or VTK)
    """
    for i in range(self.nSim):
      if self.rank < self.nPart * (i + 1):
        self.dem.importMesh(name, file, scale)
        break

  def setupWall(self, name, wtype, plane = None, peq = None):
    """
    Creates a wall
    @ name: name of the variable defining a wall or a mesh
    @ wtype: type of the wall (primitive or mesh)
    @ plane: x, y, or z plane for primitive walls
    @ peq: plane equation for primitive walls
    """
    for i in range(self.nSim):
      if self.rank < self.nPart * (i + 1):
        self.dem.setupWall(name, wtype, plane, peq)
        break

  def printSetup(self):
    """
    Specify which variables to write to file, and their format
    """
    for i in range(self.nSim):
      if self.rank < self.nPart * (i + 1):
        self.dem.printSetup()
        break

  def dumpSetup(self):
    """
    """
    for i in range(self.nSim):
      if self.rank < self.nPart * (i + 1):
        self.dem.dumpSetup()
        break

  def setupIntegrate(self, name):
    """
    Specify how Newton's eqs are integrated in time. 
    @ name: name of the fixed simulation ensemble applied to all atoms
    @ dt: timestep
    @ ensemble: ensemble type (nvt, nve, or npt)
    @ args: tuple args for npt or nvt simulations
    """
    for i in range(self.nSim):
      if self.rank < self.nPart * (i + 1):
        self.dem.setupIntegrate(name)
        break

  def integrate(self, steps, dt):
    """
    Run simulation in time
    """
    for i in range(self.nSim):
      if self.rank < self.nPart * (i + 1):
        self.dem.integrate(steps, dt)
        break

  def remove(self, name):
    """
    Delete variable/object
    """
    for i in range(self.nSim):
      if self.rank < self.nPart * (i + 1):
        self.dem.remove(name)
        break

  def monitor(self, name, group, var, file):
    """
    """
    for i in range(self.nSim):
      if self.rank < self.nPart * (i + 1):
        self.dem.monitor(name, group, var, file)
        break

  def plot(self, fname, xlabel, ylabel, output=None, xscale=None):
    """
    """
    for i in range(self.nSim):
      if self.rank < self.nPart * (i + 1):
        self.dem.plot(fname, xlabel, ylabel, output, xscale)
        break

  def saveas(self, name, fname):
    """
    """
    for i in range(self.nSim):
      if self.rank < self.nPart * (i + 1):
        self.dem.saveas(name, fname)
        break

  def __del__(self):

    MPI.Finalize()
	
