# ----------------------------------------------------------------------
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
# -------------------------------------------------------------------------

# Python wrapper on LAMMPS library via ctypes

import sys,traceback,types
from ctypes import *

import numpy as np
import itertools
import logging
from numpy.linalg import norm
from scipy import spatial
from mpi4py import MPI

logging.basicConfig(filename='dem.log', format='%(asctime)s:%(levelname)s: %(message)s', level=logging.DEBUG)

class liggghts:
  def __init__(self,name="",cmdargs=None):

    # load liblammps.so by default
    # if name = "g++", load liblammps_g++.so
    
    try:
      if not name: self.lib = CDLL("Liggghts/libliggghts.so",RTLD_GLOBAL)
      else: self.lib = CDLL("libliggghts_%s.so" % name,RTLD_GLOBAL)
    except:
      type,value,tb = sys.exc_info()
      traceback.print_exception(type,value,tb)
      raise OSError, "Could not load LIGGGHTS dynamic library"

    # create an instance of LAMMPS
    # don't know how to pass an MPI communicator from PyPar
    # no_mpi call lets LAMMPS use MPI_COMM_WORLD
    # cargs = array of C strings from args
    
    if cmdargs:
      cmdargs.insert(0,"liggghts.py")
      narg = len(cmdargs)
      cargs = (c_char_p*narg)(*cmdargs)
      self.lmp = c_void_p()
      self.lib.lammps_open_no_mpi(narg,cargs,byref(self.lmp))
    else:
      self.lmp = c_void_p()
      self.lib.lammps_open_no_mpi(0,None,byref(self.lmp))
      # could use just this if LAMMPS lib interface supported it
      # self.lmp = self.lib.lammps_open_no_mpi(0,None)

  def __del__(self):
    if self.lmp: self.lib.lammps_close(self.lmp)

  def close(self):
    self.lib.lammps_close(self.lmp)
    self.lmp = None

  def file(self,file):
    self.lib.lammps_file(self.lmp,file)

  def command(self,cmd):
    self.lib.lammps_command(self.lmp, cmd)

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
    if type == 0:
      if style > 0: return None
      self.lib.lammps_extract_fix.restype = POINTER(c_double)
      ptr = self.lib.lammps_extract_fix(self.lmp,id,style,type,i,j)
      result = ptr[0]
      self.lib.lammps_free(ptr)
      return result
    if type == 1:
      self.lib.lammps_extract_fix.restype = POINTER(c_double)
      ptr = self.lib.lammps_extract_fix(self.lmp,id,style,type,i,j)
      return ptr
    if type == 2:
      self.lib.lammps_extract_fix.restype = POINTER(POINTER(c_double))
      ptr = self.lib.lammps_extract_fix(self.lmp,id,style,type,i,j)
      return ptr
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

class DEM:
  """A class that implement  apython interface for DEM simulations"""

  def __init__(self, units, dim, style, **pargs):
    """ Initialize some settings and specifications 
    @ units: unit system (si, cgs, etc.)
    @ dim: dimensions of the problems (2 or 3)
    """
    self.comm = MPI.COMM_WORLD
    logging.info("Initialing MPI for a total of %d procs" % (self.comm.Get_size()))

    logging.info('Initializing LAMMPS object')

    self.lmp = liggghts()
    self.pargs = pargs
    self.monitorList = []
    self.vars = []

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
    logging.info('Creating domain')

    self.lmp.command('region domain block {} {} {} {} {} {} units box'.format(*self.pargs['box']))
    self.lmp.command('create_box {} domain'.format(self.pargs['nSS'] + 1))

  def insertParticles(self):
    """ Create atoms in a pre-defined region
    @ N: max total number of particles to be inserted
    @ density: initial density of the particles
    @ vel: 3 x 1 tuple of initial velocities of all particles
    @ args: dictionary of params
    """
    logging.info('Inserting particles')

    for ss in range(self.pargs['nSS']):
    
      radius = self.pargs['radius'][ss]
      density = self.pargs['density'][ss]

      self.lmp.command('fix pts all particletemplate/sphere 1 atom_type 1 density constant {} radius'.format(density) + ' %s' * len(radius) % (radius))
      self.lmp.command('fix pdd all particledistribution/discrete 63243 1 pts 1.0')
      self.lmp.command('region factory sphere 0 2.0 0 0.5 units box')
      
      self.lmp.command('fix ins all insert/rate/region seed 123481 distributiontemplate pdd nparticles {} particlerate {} insert_every {} overlapcheck yes vel constant {} {} {} region factory ntry_mc 1000'.format(self.pargs['Natoms'][ss], self.pargs['insertRate'][ss], self.pargs['insertFreq'][ss], *self.pargs['vel'][ss]))
      #self.lmp.command('fix myInsRate all insert/rate/region seed 123481 distributiontemplate pdd \
       #nparticles {} particlerate {} insert_every {} \
       #overlapcheck yes vel constant {} region factory ntry_mc 10000'.format(self.pargs['Natoms'][ss], self.pargs['insertRate'][ss], self.pargs['insertFreq'][ss], \
       #*self.pargs['vel'][ss]))

  def importMesh(self, var):
    """
    """
    logging.info('importing mesh')
    fname = self.pargs['mesh']

    logging.info('importing mesh from {}'.format(fname))
    self.lmp.command('fix {} all mesh/surface file {} type 2 scale {}'.format(var, fname, self.pargs['scaleMesh']))

  def setupWall(self, var, wtype, mesh = None, plane = None, peq = None):
    """
    Use the imported mesh as granular wall
    """

    if wtype == 'mesh':
      self.lmp.command('fix myMesh all wall/gran model hooke {} n_meshes 1 meshes {}'.format(wtype, var))
    elif wtype == 'primitive':
      self.lmp.command('fix {} all wall/gran model hooke {} type 1 {} {}'.format(var, wtype, plane, peq))
    else:
      raise ValueError('Wall type can be either primitive or mesh')

  def createGroup(self, group = None):
    """ Create groups of atoms 
    """
    logging.info('Creating atom group {}'.format(group))

    if group is None:
      for idSS in self.pargs['idSS']:
        self.lmp.command('group group{} type {}'.format(idSS, idSS))

  def setupNeighbor(self):
    """
    """
    logging.info('Setting up nearest neighbor searching parameters')
    self.lmp.command('neighbor 0.03 bin')
    self.lmp.command('neigh_modify delay 0')

  def createProperty(self, var, prop, type, valueProp, valueType = None):
    """
    Material and interaction properties required
    """
    logging.info('Creating proprety')
    self.lmp.command('fix {} all property/global {} {} {}'.format(var, prop, type, valueProp))

  def setupPhysics(self):
    """
    Specify the interation forces
    """
    logging.info('Setting up interaction parameters')

    self.lmp.command('pair_style ' + ' %s '* len(self.pargs['model']) % self.pargs['model'])
    self.lmp.command('pair_coeff * *')

  def setupGravity(self):
    """
    Specify in which direction the gravitational force acts
    """
    self.lmp.command('fix myGravity all gravity {} vector {} {} {}'.format(*self.pargs['gravity']))

  def initialize(self):
    """
    """
    self.lmp.command('restart {} {}'.format(*self.pargs['restart']))
    self.createDomain()
    #self.createGroup()
    self.setupNeighbor()
    self.setupPhysics()
    self.insertParticles()
    self.setupGravity()

  def setupIntegrate(self, name, dt = None):
    """
    Specify how Newton's eqs are integrated in time. 
    @ name: name of the fixed simulation ensemble applied to all atoms
    @ dt: timestep
    @ ensemble: ensemble type (nvt, nve, or npt)
    @ args: tuple args for npt or nvt simulations
    """

    logging.info('Setting up integration scheme parameters')

    self.lmp.command('fix {} all nve/sphere'.format(name))

    if dt is None:
      self.lmp.command('timestep {}'.format(self.pargs['dt']))

  def integrate(self, steps):
    """
    Run simulation in time
    """
    logging.info('Integrating the system for {} steps'.format(steps))

    self.lmp.command('run {}'.format(steps))

    for tup in self.monitorList:
      self.lmp.command('compute {} {} {}'.format(*tup))
      self.vars.append(self.lmp.extract_compute(tup[0], 0, 0))
      self.lmp.command('uncompute {}'.format(tup[0]))

  def printSetup(self, freq):
    """
    Specify which variables to write to file, and their format
    """
    logging.info('Setting up printing options')

    self.lmp.command('thermo_style custom' + ' %s '*len(self.pargs['print']) % self.pargs['print'])
    self.lmp.command('thermo {}'.format(freq))
    self.lmp.command('thermo_modify norm no lost ignore')

  def dumpSetup(self, sel, freq, traj):
    """
    """
    logging.info('Setting up trajectory i/o')
    traj, trajFormat = traj.split('.')

    self.lmp.command('dump dump {} xyz {} {}.{}'.format(sel, freq, traj, trajFormat))

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

  def monitor(self, name, group, var):
    """
    """
    self.monitorList.append((name, group, var))

  def __del__(self):
    """ Destructor
    """
    self.lmp.close()
    MPI.Finalize()
