# !/usr/bin/python
# -*- coding: utf8 -*- 
#
# ----------------------------------------------------------------------
#
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
'''
Created on March 30, 2016
@author: Andrew Abi-Mansour
'''

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
import glob
import sys
from importlib import import_module

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
 
  def __init__(self, library=None, style = 'granular', dim = 3, units = 'si', path=None, cmdargs=[], ptr=None, comm=None):

    comm = MPI.COMM_WORLD

    if library:
      if not comm.Get_rank():
        print "Using " + library + " as a shared library for DEM computations"
    else:
      if not comm.Get_rank():
        print "Make sure " + library + " is properly installed on your system"
      else:
        print 'Catastrophic FAILURE: library {} detected by one processor but not found by another'.format(library)
      sys.exit()

    self.lib = CDLL(library, RTLD_GLOBAL)

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
      self.lib.lammps_gather_atoms(self.lmp, name, type, count, data)
    elif type == 1:
      data = ((count*natoms)*c_double)()
      self.lib.lammps_gather_atoms(self.lmp, name, type, count, data)
    else: return None
    return data

  # scatter vector of atom properties across procs, ordered by atom ID
  # assume vector is of correct type and length, as created by gather_atoms()

  def scatter_atoms(self,name,type,count,data):
    self.lib.lammps_scatter_atoms(self.lmp,name,type,count,data)

class DEMPy:
  """A class that implements a python interface for DEM computations"""
  # TODO: This class should be generic (not specific to liggghts), must
  # handle all I/O, garbage collection, etc. and then moved to DEM.py

  def __init__(self, sid, split, library, units, dim, style, **pargs):
    """ Initialize some settings and specifications 
    @ units: unit system (si, cgs, etc.)
    @ dim: dimensions of the problem (2 or 3)
    # style: granular, atom, or ...
    """
      
    if 'print' not in pargs:
      pargs['print'] = (10**4, 'time', 'atoms')

    self.rank = split.Get_rank()
    self.split = split
    self.pargs = pargs
    self.monitorList = []
    self.vars = {}
    self.path = os.getcwd()
    self.nSS = len(self.pargs['SS'])
    self.output = self.pargs['output']

    if not self.rank:
      global logging

    os.chdir(self.output)

    if not self.rank:
      logging = import_module(name='logging')

      logging.basicConfig(filename='dem.log', format='%(asctime)s:%(levelname)s: %(message)s', level=logging.DEBUG)

      logging.info("Working in {}".format(self.path))
      logging.info('Creating i/o directories')

      if not os.path.exists(self.pargs['traj']['dir']):
        os.makedirs(self.pargs['traj']['dir'])

      if not os.path.exists(self.pargs['restart'][1]):
        os.makedirs(self.pargs['restart'][1])

      logging.info('Instantiating LIGGGHTS object')

    self.lmp = liggghts(comm=split, library=library.strip(), cmdargs=['-log', 'liggghts.log'])

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

    if not self.rank:
      from sys import argv
      scriptFile = argv[0]
      logging.info('Backing up {} file'.format(scriptFile))
      os.system('cp {}/{} {}'.format(self.path, scriptFile, scriptFile.split('.')[0] + '-bk.py'))

  def createDomain(self):
    """ Define the domain of the simulation
    @ nsys: number of subsystems
    @ pos: 6 x 1 tuple that defines the boundaries of the box 
    """
    if not self.rank:
      logging.info('Creating domain')

    if 'box' in self.pargs:
      self.lmp.command('region domain block {} {} {} {} {} {} units box'.format(*self.pargs['box']))
    elif 'cylinder' in self.pargs:
      self.lmp.command('region domain cylinder {} {} {} {} {} {} units box'.format(*self.pargs['cylinder']))

    self.lmp.command('create_box {} domain'.format(self.pargs['nSS']))

  def setupParticles(self):
    """ Setup particle for insertion if requested by the user
    """

    for ss in self.pargs['SS']:

      # Make sure we are setting up particles, not walls (so we check for id existence)
      if 'id' in ss:
        if not self.rank:
          logging.info('Setting up particles for group{id}'.format(**ss))

        if 'insert' in ss:
          radius = ss['radius']

          randName = np.random.randint(10**5,10**8)
          pddName = 'pdd' + '{}'.format(np.random.randint(10**5,10**8))

          if 'vol_lim' not in ss:
            ss['vol_lim'] = 1e-12

          self.lmp.command('group group{id} type {id}'.format(**ss))
          self.lmp.command('fix {} '.format(randName) + 'group{id} particletemplate/sphere 15485867 volume_limit {vol_lim} atom_type {id} density constant {density} radius'.format(**ss) + (' {}' * len(radius)).format(*radius))
          self.lmp.command('fix {} '.format(pddName) + 'group{id} particledistribution/discrete 67867967 1'.format(**ss) + ' {} 1.0'.format(randName))

          #Do NOT unfix randName! Will cause a memory corruption error
          self.pddName.append(pddName)

  def insert(self, name, species, *region):
    """
    This function inserts particles, and assigns particle velocities if requested by the user. 
    """
    if not self.pddName:
      print 'Probability distribution not set for particle insertion. Exiting ...'
      sys.exit()

    def insert_loc(self, ss, i, name, *region):
      if 'insert' in ss:
        if not self.rank:
          logging.info('Inserting particles for species {}'.format(i))
    
        if 'natoms_local' in ss:
          natoms = ss['natoms_local']
        else:
          natoms = ss['natoms'] - self.lmp.get_natoms()

        if natoms < 0:
          if not self.rank: 
            print 'Too many particles requested for insertion. Increase the total number of particles in your system.'
          raise

        if natoms > 0:
          randName = 'insert' + '{}'.format(np.random.randint(0,10**6))
          self.lmp.command('region {} '.format(name) + ('{} ' * len(region)).format(*region) + 'units box')

          if ss['insert'] == 'by_rate':
            self.lmp.command('fix {} group{} insert/rate/region seed 123481 distributiontemplate {} nparticles {}'.format(randName, ss['id'], self.pddName[i], natoms) + \
              ' particlerate {rate} insert_every {freq} overlapcheck yes vel constant'.format(**ss) \
              + ' {} {} {}'.format(*self.pargs['vel'][i])  + ' region {} ntry_mc 1000'.format(name) )
          elif ss['insert'] == 'by_pack':
            self.lmp.command('fix {} group{} insert/pack seed 123481 distributiontemplate {}'.format(randName, ss['id'], self.pddName[i]) + \
              ' insert_every {freq} overlapcheck yes vel constant'.format(**ss) \
              + ' {} {} {}'.format(*self.pargs['vel'][i])  + ' particles_in_region {} region {} ntry_mc 1000'.format(natoms, name) )
          else:
            print 'WARNING: Insertion mechanism not specified by user. Assuming insertion by rate ...'
            self.lmp.command('fix {} group{} insert/rate/region seed 123481 distributiontemplate {} nparticles {}'.format(randName, ss['id'], self.pddName[i], natoms) + \
              ' particlerate {rate} insert_every {freq} overlapcheck yes vel constant'.format(**ss) \
              + ' {} {} {}'.format(*self.pargs['vel'][i])  + ' region {} ntry_mc 1000'.format(name) )
        else:
          if not self.rank:
            print 'WARNING: no more particles to insert. Ignoring user request for more insertion ...'
          raise

        return randName
      else:
        return None

    if species != 'all':
      i, ss = species - 1, self.pargs['SS'][species - 1]
      randName = insert_loc(self, ss, i, name, *region)
    else:
      randName = []
      for i, ss in enumerate(self.pargs['SS']):
        randName.append(insert_loc(self, ss, i, name, *region))

    # Check if the user has supplied any initial velocities
    if 'velocity' in self.pargs:
      for comp in self.pargs['velocity']:
        self.velocity(*comp)

    return randName

  def run(self, nsteps, dt=None):
    """ runs a simulation for number of steps specified by the user """
      
    self.integrator = self.setupIntegrate(name=np.random.randint(0,10**6))

    if not dt:
      if 'dt' in self.pargs:
        dt = self.pargs['dt']
      else:
        if not self.rank:
          print 'Could not find dt in user-supplied dictionary. Aborting ...'
        sys.exit()

    # check timestep
    self.lmp.command('fix ts_check all check/timestep/gran 1000 0.25 0.25')


    self.integrate(nsteps, dt)

  def moveMesh(self, name, *args):
    self.lmp.command('fix moveMesh all move/mesh mesh {} '.format(name) + ('{} ' * len(args)).format(*args))

  def importMesh(self, name, file, mtype, *args):
    """
    TODO: fix type for mesh
    """
    fname = self.path + '/' + file

    if not self.rank:
      logging.info('Importing mesh from {}'.format(fname))
      
    self.lmp.command('fix {} all {} file {} type 2 '.format(name, mtype, fname) + ('{} ' * len(args)).format(*args))
    
  def setupWalls(self, name, wtype, meshName = None, plane = None, peq = None):
    """
    Creates a wall
    @ name: name of the variable defining a wall or a mesh
    @ wtype: type of the wall (primitive or mesh)
    @ plane: x, y, or z plane for primitive walls
    @ peq: plane equation for primitive walls

    This function can be called only ONCE for setting up all walls (restriction from LIGGGHTS)
    """

    gran = 'gran' # VERY HACKISH
    model = []
    modelExtra = []

    for item in self.pargs['model-args']:
      if item != 'gran' and item != 'tangential_damping' and item != 'on' and item != 'limitForce' and item != 'ktToKnUser' \
        and item != 'off' and item != 'radiusGrowth':
        model.append(item)
      elif item != 'gran':
        modelExtra.append(item)

    model = tuple(model)

    if wtype == 'mesh':
      nMeshes = len(self.pargs['mesh'])
      meshName = tuple(self.pargs['mesh'].keys())
      self.lmp.command('fix walls all wall/{} '.format(gran) + ('{} ' * len(model)).format(*model) + ' {} n_meshes {} meshes'.format(wtype, nMeshes) \
        + (' {} ' * len(meshName)).format(*meshName)  + ('{} ' * len(modelExtra)).format(*modelExtra))
    elif wtype == 'primitive':
      self.lmp.command('fix {} all wall/{} '.format(name, gran) + ('{} ' * len(model)).format(*model) +  '{} type 1 {} {}'.format(wtype, plane, peq))
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

  def setupNeighbor(self, **params):
    """
    """
    if not self.rank:
      logging.info('Setting up nearest neighbor searching parameters')

    self.lmp.command('neighbor {nns_skin} {nns_type}'.format(**params))
    self.lmp.command('neigh_modify delay 0 every 100')

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

    args = self.pargs['model-args']

    self.lmp.command('pair_style ' + (' {}' * len(args)).format(*args))
    self.lmp.command('pair_coeff * *')

  def velocity(self, *args):
    """
    Assigns velocity to selected particles.
    """
    self.lmp.command('velocity' + (' {}' * len(args)).format(*args))

  def setupGravity(self):
    """
    Specify in which direction the gravitational force acts
    """
    if 'gravity' in self.pargs:
      self.lmp.command('fix myGravity all gravity {} vector {} {} {}'.format(*self.pargs['gravity']))

  def initialize(self, **params):
    """
    """

    self.lmp.command('restart {} {}/{}'.format(*self.pargs['restart'][:-1]))
    self.pddName = []
    self.integrator = None

    if self.pargs['restart'][3] == False and self.pargs['read_data'] == False:

      self.createDomain()
      #self.createGroup()
      self.setupPhysics()
      self.setupNeighbor(**self.pargs)
      self.setupParticles()
      self.setupGravity()

    elif self.pargs['read_data']:
      self.createDomain()
      self.readData()
      self.setupPhysics()
      self.setupNeighbor(**self.pargs)
      self.setupParticles()
      self.setupGravity()

    else:
      self.resume()
      self.setupPhysics()
      self.setupNeighbor(**self.pargs)
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

    if self.integrator:
      self.remove(self.integrator)

    self.lmp.command('fix {} all nve/sphere'.format(name))

    return name

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

    self.lmp.command('dump dump {sel} {style} {freq} {dir}/{pfile}'.format(**self.pargs['traj']) + (' {} ' * len(self.pargs['traj']['args'])).format(*self.pargs['traj']['args']))

    self.lmp.command('dump_modify dump ' +  (' {} ' * len(self.pargs['dump_modify'])).format(*self.pargs['dump_modify']))

    if 'mfile' in self.pargs['traj']:
      self.lmp.command('dump meshDump all mesh/vtk {freq} {dir}/{mfile} id stress stresscomponents vel '.format(**self.pargs['traj']))

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

  def command(self, cmd):
    """
    Passes a specific command to LIGGGHTS 
    """
    self.lmp.command(cmd)

  def resume(self):
    """
    """
    rdir = '{}/*'.format(self.pargs['restart'][1])

    if self.pargs['restart'][-1]:
      rfile = self.pargs['restart'][1] + '/' + self.pargs['restart'][-1]
    else:
      rfile = max(glob.iglob(rdir), key=os.path.getctime)

    self.lmp.command('read_restart {}'.format(rfile))

  def readData(self):
    """
    """
    args = self.pargs['read_data']

    self.lmp.command('read_dump '  + (' {}' * len(args)).format(*args))

  def __del__(self):
    """ Destructor
    """
    self.lmp.close()
