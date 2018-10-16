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
# -----------------------------------------------------------------------
#
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

'''
Created on March 30, 2016
@author: Andrew Abi-Mansour
'''

import sys,traceback,types
import ctypes
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
from PyGran.Tools import find

class liggghts:
  # detect if Python is using version of mpi4py that can pass a communicator

  try:
    from mpi4py import MPI
  except:
    pass

  # create instance of LIGGGHTS
  def __init__(self, library=None, style = 'spherical', dim = 3, units = 'si', path=None, cmdargs=[], ptr=None, comm=None):

    if not comm:
      comm = MPI.COMM_WORLD

    if library:
      if not comm.Get_rank():
        print("Using " + library + " as a shared library for DEM computations")
    else:
      if not comm.Get_rank():
        if library:
          print("Make sure " + library + " is properly installed on your system")
        else:
          print("No liggghts library supplied. Exiting ...")
      else:
        print('Catastrophic FAILURE: library {} detected by one processor but not found by another'.format(library))
      sys.exit()

    try:
      self.lib = ctypes.CDLL(library, ctypes.RTLD_GLOBAL)
    except:
      etype,value,tb = sys.exc_info()
      traceback.print_exception(etype,value,tb)
      raise RuntimeError("Could not load LIGGGHTS dynamic library")

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

        narg = 0
        cargs = 0
        if cmdargs:
          cmdargs.insert(0,"liggghts.py")
          narg = len(cmdargs)
          cargs = [st.encode('utf-8') for st in cmdargs]
          cargs = (ctypes.c_char_p*narg)(*cargs)
          self.lmp = ctypes.c_void_p()
          self.lib.lammps_open_no_mpi(narg,cargs,ctypes.byref(self.lmp))
        else:
          self.lmp = ctypes.c_void_p()
          self.lib.lammps_open_no_mpi(0,None,ctypes.byref(self.lmp))

        self.opened = True
    else:
      self.opened = False
      # magic to convert ptr to ctypes ptr
      pythonapi.PyCObject_AsVoidPtr.restype = ctypes.c_void_p
      pythonapi.PyCObject_AsVoidPtr.argtypes = [py_object]
      self.lmp = ctypes.c_void_p(pythonapi.PyCObject_AsVoidPtr(ptr))

  def __del__(self):
    if hasattr(self, 'lmp') and self.opened: 
      self.close()

  def close(self):
    if self.opened: 
        self.lib.lammps_close(self.lmp)
        self.lmp = None

  def file(self,file):
    self.lib.lammps_file(self.lmp,file)

  def command(self,cmd):
    """ For python 3, I had to encode the string as an 8 character utf """
    self.lib.lammps_command(self.lmp,cmd.encode('utf-8'))

  def extract_global(self,name,type):
    if type == 0:
      self.lib.lammps_extract_global.restype = ctypes.POINTER(ctypes.c_int)
    elif type == 1:
      self.lib.lammps_extract_global.restype = ctypes.POINTER(ctypes.c_double)
    else: return None
    ptr = self.lib.lammps_extract_global(self.lmp,name)
    return ptr[0]

  def extract_atom(self,name,type):
    if type == 0:
      self.lib.lammps_extract_atom.restype = ctypes.POINTER(ctypes.c_int)
    elif type == 1:
      self.lib.lammps_extract_atom.restype = ctypes.POINTER(ctypes.POINTER(ctypes.c_int))
    elif type == 2:
      self.lib.lammps_extract_atom.restype = ctypes.POINTER(ctypes.c_double)
    elif type == 3:
      self.lib.lammps_extract_atom.restype = ctypes.POINTER(ctypes.POINTER(ctypes.c_double))
    else: return None
    ptr = self.lib.lammps_extract_atom(self.lmp,name)
    return ptr

  def extract_compute(self,id,style,type):
    if type == 0:
      if style > 0: return None
      self.lib.lammps_extract_compute.restype = ctypes.POINTER(ctypes.c_double)
      ptr = self.lib.lammps_extract_compute(self.lmp,id,style,type)
      return ptr[0]
    if type == 1:
      self.lib.lammps_extract_compute.restype = ctypes.POINTER(ctypes.c_double)
      ptr = self.lib.lammps_extract_compute(self.lmp,id,style,type)
      return ptr
    if type == 2:
      self.lib.lammps_extract_compute.restype = ctypes.POINTER(ctypes.POINTER(ctypes.c_double))
      ptr = self.lib.lammps_extract_compute(self.lmp,id,style,type)
      return ptr
    return None

  # in case of global datum, free memory for 1 double via lammps_free()
  # double was allocated by library interface function

  def extract_fix(self,id,style,type,i=0,j=0):
    if style == 0:
      self.lib.lammps_extract_fix.restype = ctypes.POINTER(ctypes.c_double)
      ptr = self.lib.lammps_extract_fix(self.lmp,id,style,type,i,j)
      result = ptr[0]
      self.lib.lammps_free(ptr)
      return result
    elif (style == 1) or (style == 2):
      if type == 1:
        self.lib.lammps_extract_fix.restype = ctypes.POINTER(ctypes.c_double)
      elif type == 2:
        self.lib.lammps_extract_fix.restype = ctypes.POINTER(ctypes.POINTER(ctypes.c_double))
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
      self.lib.lammps_extract_variable.restype = ctypes.POINTER(ctypes.c_double)
      ptr = self.lib.lammps_extract_variable(self.lmp,name,group)
      result = ptr[0]
      self.lib.lammps_free(ptr)
      return result
    if type == 1:
      self.lib.lammps_extract_global.restype = ctypes.POINTER(ctypes.c_int)
      nlocalptr = self.lib.lammps_extract_global(self.lmp,"nlocal")
      nlocal = nlocalptr[0]
      result = (ctypes.c_double*nlocal)()
      self.lib.lammps_extract_variable.restype = ctypes.POINTER(ctypes.c_double)
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
      data = ((count*natoms)*ctypes.c_int)()
      self.lib.lammps_gather_atoms(self.lmp, name, type, count, data)
    elif type == 1:
      data = ((count*natoms)*ctypes.c_double)()
      self.lib.lammps_gather_atoms(self.lmp, name, type, count, data)
    else: return None
    return data

  # scatter vector of atom properties across procs, ordered by atom ID
  # assume vector is of correct type and length, as created by gather_atoms()

  def scatter_atoms(self,name,type,count,data):
    self.lib.lammps_scatter_atoms(self.lmp, name, type, count, data)

class DEMPy:
  """A class that implements a python interface for DEM computations"""
  # TODO: This class should be generic (not specific to liggghts), must
  # handle all I/O, garbage collection, etc. and then moved to DEM.py

  def __init__(self, sid, split, library, style, **pargs):
    """ Initialize some settings and specifications
    @ units: unit system (si, cgs, etc.)
    @ dim: dimensions of the problem (2 or 3)
    # style: granular, atom, or ...
    """

    if 'print' not in pargs:
      pargs['print'] = (10**4, 'time', 'dt', 'atoms')

    self.rank = split.Get_rank()
    self.split = split
    self.pargs = pargs
    self.monitorList = []
    self.vars = {}
    self.path = os.getcwd()
    self.nSS = len(self.pargs['species'])
    self.output = self.pargs['output']
    self._dir, _ = __file__.split(__name__.split('PyGran.Simulator.')[-1] +'.py')

    if '__version__' in pargs:
      self.__version__ = self.pargs['__version__']

    if not self.rank:
      global logging

      logging = import_module(name='logging')

      logging.basicConfig(filename='dem.log', format='%(asctime)s:%(levelname)s: %(message)s', level=logging.DEBUG)

      logging.info("Working in {}".format(self.path))
      logging.info('Creating i/o directories')

      if not os.path.exists(self.pargs['traj']['dir']):
        os.makedirs(self.pargs['traj']['dir'])

      if self.pargs['restart']:
        if not os.path.exists(self.pargs['restart'][1]):
          os.makedirs(self.pargs['restart'][1])

      logging.info('Instantiating LIGGGHTS object')

    self.lmp = liggghts(comm=self.split, library=library.strip(), cmdargs=['-log', 'liggghts.log'])

    if not self.rank:
      logging.info('Setting up problem dimensions and boundaries')

    self.lmp.command('units {}'.format(self.pargs['units']))

    if hasattr(self, '__version__'): 
      if self.__version__ >= 3.6:
        self.lmp.command('hard_particles yes')
    else:
      # Get version from version_liggghts.txt. TODO: find a faster way to do this.
      with open(find('version_liggghts.txt', '/'), 'r+') as fp:
        major, minor, _ = fp.readline().rstrip().split('.')
        self.__version__ = float(major + '.' + minor)

        # Write version to config file if it exists
        if not self.rank:
          if os.path.isfile(self._dir + '../.config'):
            with open(self._dir + '../.config', 'a+') as fp:
              fp.write('\nversion={}'.format(self.__version__))

      if self.__version__ >= 3.6:
        self.lmp.command('hard_particles yes')

    self.lmp.command('dimension {}'.format(self.pargs['dim']))
    self.lmp.command('atom_style {}'.format(style))
    self.lmp.command('atom_modify map array') # array is faster than hash in looking up atomic IDs, but the former takes more memory
    self.lmp.command('boundary ' + ('{} ' * len(pargs['boundary'])).format(*pargs['boundary']))
    self.lmp.command('newton off') # turn off newton's 3rd law ~ should lead to better scalability
    self.lmp.command('communicate single vel yes') # have no idea what this does, but it's imp for ghost atoms
    self.lmp.command('processors * * *') # let LIGGGHTS handle DD

  def get_natoms(self):
    return self.lmp.get_natoms()

  def scatter_atoms(self,name,type,count,data):
    return self.lmp.scatter_atoms(name,type,count,data)

  def gather_atoms(self,name,type,count):
    return self.lmp.gather_atoms(name,type,count)

  def extract_global(self,name,type):
    return self.lmp.extract_global(name, type)

  def extract_compute(self,id,style,type):
    return self.lmp.extract_compute(id,style,type)

  def createParticles(self, type, style, *args):
    self.lmp.createParticles(type, style, *args)

  def createDomain(self):
    """ Define the domain of the simulation
    @ nsys: number of subsystems
    @ pos: 6 x 1 tuple that defines the boundaries of the box
    """
    if not self.rank:
      logging.info('Creating domain')

    if 'box' in self.pargs:
      self.lmp.command('region domain block ' + ('{} ' * len(self.pargs['box'])).format(*self.pargs['box']) + ' units box volume_limit 1e-20')
    elif 'cylinder' in self.pargs:
      self.lmp.command('region domain cylinder ' + ('{} ' * len(self.pargs['cylinder'])).format(*self.pargs['cylinder']) + ' units box volume_limit 1e-20') 

    self.lmp.command('create_box {} domain'.format(self.pargs['nSS']))

  def setupParticles(self):
    """ Setup particle for insertion if requested by the user
    """

    for ss in self.pargs['species']:

      # Make sure we are setting up particles, not walls (so we check for id existence)
      if 'id' in ss and 'wall' not in ss:
        if not self.rank:
          logging.info('Setting up particles for group{id}'.format(**ss))

        randName = np.random.randint(10**5,10**8)
        pddName = 'pdd' + '{}'.format(np.random.randint(10**5,10**8))

        if 'vol_lim' not in ss:
          ss['vol_lim'] = 1e-20

        id = ss['id'] - 1
        self.lmp.command('group group{} type {}'.format(id, ss['id']))

        if 'args'in ss:
          args = ss['args']
        else:
          args = ()

        if 'radius' in ss:
          radius = ss['radius']
          self.lmp.command('fix {} '.format(randName) + 'group{}'.format(id) + ' particletemplate/{style} 15485867 volume_limit {vol_lim} atom_type {id} density constant {density} radius'.format(**ss) + (' {}' * len(radius)).format(*radius) \
          + (' {}' * len(args)).format(*args))
        else:
          self.lmp.command('fix {} '.format(randName) + 'group{}'.format(id) + ' particletemplate/{style} 15485867 volume_limit {vol_lim} atom_type {id} density constant {density}'.format(**ss) + (' {}' * len(args)).format(*args))
        
        self.lmp.command('fix {} '.format(pddName) + 'group{}'.format(id) + ' particledistribution/discrete 67867967 1'.format(**ss) + ' {} 1.0'.format(randName))

        if ss['style'] is 'multisphere':
          itype = ss['style']
        else:
          itype = 'nve/{style}'.format(**ss)

        #Do NOT unfix randName! Will cause a memory corruption error
        self.pddName.append(pddName)

  def insert(self, species, value, **args):
    """
    This function inserts particles, and assigns particle velocities if requested by the user. If species is 'all',
    all components specified in SS are inserted. Otherwise, species must be the id of the component to be inserted.

    region: tuple of the form ('shape', (xmin, xmax, ymin, ymax, zmin, zmax)) or ('shape', xmin, xmax, ymin, ymax, zmin, zmax)

    TODO: support insertion of all or multiple species at the same time for multiple regions. 
    """
    if not self.pddName:
      print('Probability distribution not set for particle insertion. Exiting ...')
      sys.exit()

    if 'region' in args:
      region = args['region']
    else:
      # Default region is sim box
      if 'cylinder' in self.pargs:
        region = ('cylinder', self.pargs['cylinder'])
      else:
        region = ('block', self.pargs['box'])
        
      region = tuple([region[0]] + [i for i in region[1:][0]])
      args['region'] = region

    # I think this is for creating tuples of lists, corresponding to many regions
    # This is prolly for inserting many species at the same time in different regions
    if isinstance(region[1], tuple):
      targs = list(region[1])
      targs.insert(0, region[0])

      tmp = region[1]
      for i in range(len(tmp)):
        targs[i+1] = tmp[i]

      region = tuple(targs)

    def insert_loc(self, id, value, vel, vel_type, region, mech, **ss):
      """ For multi-component system, this function can cause REAL *trouble*. For now, make sure components
      are inserted consecutively or all at once.

      TODO: let the user override volume_limit
      """

      if not self.rank:
        logging.info('Inserting particles for species {}'.format(id+1))

      seed = 32452843
      name = np.random.randint(0,1e8)

      randName = 'insert' + '{}'.format(np.random.randint(0,10**6))
      self.lmp.command('region {} '.format(name) + ('{} ' * len(region)).format(*region) + 'units box volume_limit 1e-20')

      if 'args' not in ss:
        ss['args'] = ()

      if 'freq' not in ss:
        ss['freq'] = 'once'

      if 'all_in' not in ss:
        ss['all_in'] = 'yes'

      if 'insert' not in ss:
        ss['insert'] = 'by_pack'

      if ss['insert'] == 'by_rate':
        if not mech:
          mech = 'nparticles'

        if mech is 'nparticles':
          value += self.lmp.get_natoms()

        self.lmp.command('fix {} group{} insert/rate/region seed 123481 distributiontemplate {} {} {}'.format(randName, id, self.pddName[id], mech, value) + \
          ' particlerate {rate} insert_every {freq} overlapcheck yes all_in {all_in}'.format(**ss) + ' vel {}'.format(vel_type) \
          + (' {}' * len(vel)).format(*vel)  + (' {}' * len(ss['args'])).format(*ss['args']) + ' region {} ntry_mc 10000'.format(name) )
      elif ss['insert'] == 'by_pack':
        if not mech:
          mech = 'particles_in_region'

        if mech is 'particles_in_region':
          value += self.lmp.get_natoms()

        self.lmp.command('fix {} group{} insert/pack seed {} distributiontemplate {}'.format(randName, id, seed, self.pddName[id]) + \
          ' insert_every {freq} overlapcheck yes all_in {all_in}'.format(**ss) + ' vel {}'.format(vel_type) \
          + (' {}' * len(vel)).format(*vel)  + (' {}' * len(ss['args'])).format(*ss['args']) + ' {} {} region {} ntry_mc 10000'.format(mech, value, name) )
      else:
        print('WARNING: Insertion mechanism {insert} not found. Assuming insertion by rate ...'.format(**ss))

        if mech is 'nparticles':
          value += self.lmp.get_natoms()

        value += self.lmp.get_natoms()

        self.lmp.command('fix {} group{} insert/rate/region seed 123481 distributiontemplate {} {} {}'.format(randName, id, self.pddName[id], mech, value) + \
          ' {rate_type} {rate} insert_every {freq} overlapcheck yes all_in {all_in}'.format(**ss) + ' vel {}'.format(vel_type) \
          + (' {}' * len(vel)).format(*vel)  + (' {}' * len(ss['args'])).format(*ss['args']) + ' region {} ntry_mc 10000'.format(name))

      return randName

    if species != 'all':

      species = int(species)

      ss = self.pargs['species'][species - 1]

      if 'vel' not in args:
        args['vel'] = (0,0,0)

      if 'vel_type' not in args:
        args['vel_type'] = 'constant'

      if 'mech' not in args:
        args['mech'] = None

      randName = insert_loc(self, species - 1, value, **args)
    else:
      raise RuntimeError('Insertion species {} not supported in PyGran.'.format(species))

    return randName

  def run(self, nsteps, dt=None, itype=None):
    """ Runs a simulation for number of steps specified by the user
     @itype = sphere (rotational motion on) or rigid_sphere (rotational motion off)
     @dt = timestep"""
    
    name = self.setupIntegrate(itype=itype)

    if not dt:
      if 'dt' in self.pargs:
        dt = self.pargs['dt']
      else:
        if not self.rank:
          print('Could not find dt in user-supplied dictionary. Aborting ...')
        sys.exit()

    self.integrate(nsteps, dt)

    return name

  def moveMesh(self, name, *args):
    
    randName = 'moveMesh' + str(np.random.randint(10**5,10**8))

    self.lmp.command('fix {} all move/mesh mesh {} '.format(randName, name) + ('{} ' * len(args)).format(*args))

    return randName

  def importMeshes(self, name=None):
    """ Imports all meshes and sets them up as walls. Can import only one mesh specified by the 'name' keyword.
    @file: mesh filename
    @mtype: mesh type
    @args: additional args
    """
    wall = False

    if 'mesh' in self.pargs:
      for mesh in self.pargs['mesh'].keys():

        if 'file' in self.pargs['mesh'][mesh]:
            if name:
              if mesh == name:
                self.pargs['mesh'][mesh]['import'] = True
                self.importMesh(mesh, self.pargs['mesh'][mesh]['file'], self.pargs['mesh'][mesh]['mtype'], self.pargs['mesh'][mesh]['id'], *self.pargs['mesh'][mesh]['args'])  
                wall = True

            elif 'import' in self.pargs['mesh'][mesh]:
              if self.pargs['mesh'][mesh]['import']:
                self.importMesh(mesh, self.pargs['mesh'][mesh]['file'], self.pargs['mesh'][mesh]['mtype'], self.pargs['mesh'][mesh]['id'], *self.pargs['mesh'][mesh]['args'])  
                wall = True
              
      if wall:
        self.setupWall(wtype='mesh')
    

  def importMesh(self, name, file, mtype, material, *args):
    """
    Imports a specific surface mesh requested by the user
    """
    fname = self.path + '/../' + file
    
    if not self.rank:
      logging.info('Importing mesh from {}'.format(fname))

    self.lmp.command('fix {} all {} file {} type {} '.format(name, mtype, fname, material) + ('{} ' * len(args)).format(*args))

  def setupWall(self, wtype, species = None, plane = None, peq = None):
    """
    Creates a wall
    @ name: name of the variable defining a wall or a mesh
    @ wtype: type of the wall (primitive or mesh)
    @ plane: x, y, or z plane for primitive walls
    @ peq: plane equation for primitive walls

    This function can be called only ONCE for setting up all mesh walls (restriction from LIGGGHTS)

    TODO: support additional keywords (shear, etc.) for primitive walls
    """

    gran = 'gran' # VERY HACKISH
    model = []
    modelExtra = []

    name = np.random.randint(0,1e8)

    for item in self.pargs['model-args']:
      if item.startswith('model'):
        model.append(item)
      elif item.startswith('cohesion') or item.startswith('tangential ') or item.startswith('rolling_friction'):
        modelExtra.append(item)

    # Replace any user-specified model args for all mesh walls
    for i, key in enumerate(modelExtra):
      if wtype == 'mesh':
        kname = key.split()[0]
        if kname in self.pargs['mesh']:
          if isinstance(self.pargs['mesh'][kname], str): # make sure this is an actual mesh keyword, not a mesh defined with a keyname same as a mesh arg!
            modelExtra[i] = kname + ' ' + self.pargs['mesh'][kname]

    model = tuple(model)
    modelExtra = tuple(modelExtra)

    # Can we take model args into account for walls???

    if wtype == 'mesh':
      meshName = tuple([mname for mname in self.pargs['mesh'].keys() if 'file' in self.pargs['mesh'][mname]])
      nMeshes = len(meshName)

      self.lmp.command('fix walls all wall/{} '.format(gran) + ('{} ' * len(model)).format(*model) + ('{} ' * len(modelExtra)).format(*modelExtra) + \
      '{} n_meshes {} meshes'.format(wtype, nMeshes) + (' {} ' * nMeshes).format(*meshName))
    elif wtype == 'primitive':
      self.lmp.command('fix {} all wall/{} '.format(name, gran) + ('{} ' * len(model)).format(*model) +  ('{} ' * len(modelExtra)).format(*modelExtra) + \
        '{} type {} {} {}'.format(wtype, species, plane, peq))
    else:
      raise ValueError('Wall type can be either primitive or mesh')

    return name

  def remove(self, name):
    """
    Deletes a specified fix. If the fix is for a mesh, we must unfix it and re-import all meshes again and setup them
    up as walls. Very tedious!
    """
    # Remove any DUMP-IDS 1st in case the user wants to move a mesh
    if 'mesh' in self.pargs:
      if name in self.pargs['mesh']:
        # must delete all meshes / dumps in order to re-import remaining meshes
        for dump in self.pargs['traj']['dump_mname']:
          self.lmp.command('undump {}'.format(dump))

        self.lmp.command('unfix walls')

        for i, mesh in enumerate(self.pargs['mesh'].keys()):
          self.lmp.command('unfix {}'.format(mesh))

        if 'mfile' in self.pargs['traj']:
          if isinstance(self.pargs['traj']['mfile'], list):
            raise RuntimeError('mfile cannot be a list. Something is not setup correctly.')
          elif self.pargs['traj']['mfile']: # the user has requested all mesh(es) be written as one file
            pass
          else: # self.pargs['traj']['mfile'] had better be None
            assert(self.pargs['traj']['mfile'] is None)

        del self.pargs['mesh'][name]

        # Re-import any remaining meshes
        self.importMeshes()

        # Create new dump setups, leaving particle dumps intact
        self.dumpSetup(only_mesh=True)

        return 0
    
    # Otherwise, we are just unfixing a non-mesh fix
    self.lmp.command('unfix {}'.format(name))

  def createGroup(self, *group):
    """ Create groups of atoms. If group is empty, groups{i} are created for every i species.
    """
    if not self.rank:
      logging.info('Creating atom group {}'.format(group))

    if not len(group):
      for idSS in self.pargs['idSS']:
        self.lmp.command('group group{} type {}'.format(idSS, idSS))
    else:
      self.lmp.command('group ' + ('{} ' * len(group)).format(*group))

  def createParticles(self, type, style, *args):
    """
    Creates particles of type 'type' (1,2, ...) using style 'style' (box or region or single or random)
    @[args]: 'basis' or 'remap' or 'units' or 'all_in' 
    """
    if not self.rank:
      logging.info('Creating particles {} with args'.format(type) + (' {}' * len(args)).format(*args))

    self.lmp.command('create_atoms {} {}'.format(type, style) +  (' {}' * len(args)).format(*args))

  def set(self, *args):
    """ Set group/atom attributes """
    self.lmp.command('set ' + (' {}' * len(args)).format(*args))

  def setupNeighbor(self, **params):
    """
    Sets up NNS list parameters
    """
    if not self.rank:
      logging.info('Setting up nearest neighbor searching parameters')

    if 'nns_freq' not in params:
      params['nns_freq'] = 10

    if 'nns_skin' not in params:
      radius = 0

      for ss in params['species']:
        if 'radius' in ss:
          radius = max(radius, ss['radius'][1])

      params['nns_skin'] = radius * 4

    self.lmp.command('neighbor {nns_skin} {nns_type}'.format(**params))
    self.lmp.command('neigh_modify delay 0 every {nns_freq} check yes'.format(**params))

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

    if self.pargs['restart']:
      self.lmp.command('restart {} {}/{}'.format(*self.pargs['restart'][:-1]))
    else:
      # create dummy restart tuple to pass below
      self.pargs['restart'] = (None, None, None, False)

    self.pddName = []
    self.integrator = []

    if not self.pargs['restart'][3] and not self.pargs['read_data']:

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
      self.setupParticles()
      self.setupGravity()

    self.setupIntegrate()
    self.importMeshes()
    
    # Write output to trajectory by default unless the user specifies otherwise
    if 'dump' in self.pargs:
      if self.pargs['dump'] == True:
        self.dumpSetup()
    else:
      self.dumpSetup()

  def setupIntegrate(self, itype=None, group=None):
    """
    Specify how Newton's eqs are integrated in time. MUST BE EXECUTED ONLY ONCE.
    TODO: extend this to SQ particles
    """
    if not self.rank:
      logging.info('Setting up integration scheme parameters')

    spheres = []
    multi = []

    if not self.integrator:

      # check timestep ~ do this only ONCE
      self.lmp.command('fix ts_check all check/timestep/gran 1000 0.5 0.5')

      # Find which components (types) are spheres, multi-spheres, QS, etc.
      for i, ss in enumerate(self.pargs['species']):
        if 'id' in ss and 'wall' not in ss: # dont count mesh wall(s)
          if ss['style'] == 'sphere':
            spheres.append('{}'.format(i+1))
          elif ss['style'] == 'multisphere':
            multi.append('{}'.format(i+1))

      if len(spheres):
        #self.createGroup(*('spheres type', (' {}' * len(spheres)).format(*spheres)))

        for sphere in spheres:
          name = 'sphere_' + str(np.random.randint(0,10**6))
          if not itype:
            self.lmp.command('fix {} group{} nve/sphere'.format(name, int(sphere[0]) -1))
          else:
            self.lmp.command('fix {} group{} {}'.format(name, int(sphere[0]) -1, itype))
          self.integrator.append(name)

      # LIGGGHTS does not permit more than one multisphere group to exist / integrated
      # So we will reject any MS groups beyond the 1st
      if len(multi) > 1:
        raise RuntimeError("LIGGGHTS (3.x) does not currently support more than one multisphere group.")
      elif len(multi): # must be of length 1

        # When LIGGGHTS supports multiple multisphere groups, I should uncomment this
        #self.createGroup(*('multi type', (' {}' * len(multi)).format(*multi)))

        ms = True
        for integ in self.integrator:
          if integ.startswith('multisphere'):
            ms = False

        if ms:
          name = 'multisphere_' + str(np.random.randint(0,10**6))
          self.lmp.command('fix {} group{} multisphere'.format(name, int(multi[0])-1))
          self.integrator.append(name)

    return self.integrator

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

  def dumpSetup(self, only_mesh=False, name=None):
    """
    This creates dumps for particles and meshes in the system. In LIGGGHTS, all meshes must be declared once, so if a mesh is removed during
    the simulation, this function has to be called again, usually with only_mesh=True to keep the particle dump intact.
    """
    if not self.rank:
      logging.info('Setting up trajectory i/o')

    # Make sure the user did not request no particles be saved to a traj file, or we're not just re-initializing the meshes
    if not only_mesh and self.pargs['traj']['pfile']:

      if hasattr(self, 'dump'):
        if self.dump:
          self.lmp.command('undump dump')

      if not name:
        name = 'dump'

      self.lmp.command('dump {} '.format(name) + ' {sel} {style} {freq} {dir}/{pfile}'.format(**self.pargs['traj']) + (' {} ' * len(self.pargs['traj']['args'])).format(*self.pargs['traj']['args']))
      self.lmp.command('dump_modify {} '.format(name) +  (' {} ' * len(self.pargs['dump_modify'])).format(*self.pargs['dump_modify']))

    self.pargs['traj']['dump_mname'] = []

    # Make sure meshes are defined so we can dump them if requested (or not)
    if 'mesh' in self.pargs:
      if hasattr(self, 'dump'):
        if self.dump:
          for dname in self.pargs['traj']['dump_mname']:
            self.lmp.command('undump ' + dname)

      if 'mfile' not in self.pargs['traj']:
        for mesh in self.pargs['mesh'].keys():

          if 'file' in self.pargs['mesh'][mesh]:
          # Make sure only mehs keywords supplied with files are counter, otherwise, they're args to the mesh wall!

            if self.pargs['mesh'][mesh]['import']:
              args = self.pargs['traj'].copy()
              args['mfile'] = mesh + '-*.vtk'
              args['mName'] = mesh
              name = 'dump' + str(np.random.randint(0,1e8))
              self.pargs['traj']['dump_mname'].append(name)

              self.lmp.command('dump ' + name + ' all mesh/vtk {freq} {dir}/{mfile} id stress stresscomponents vel {mName}'.format(**args))
      elif not isinstance(self.pargs['traj']['mfile'], list):
        if self.pargs['traj']['mfile']:
          name = ''

          # see if we have many meshes to dump if the name of one mfile supplied by the user
          for mesh in self.pargs['mesh'].keys():
            if 'file' in self.pargs['mesh'][mesh]:
              if self.pargs['mesh'][mesh]['import']:
                name += mesh + ' '

          if len(name):
            name = name[:-1] # remove last space to avoid fake (empty) mesh IDs

            args = self.pargs['traj'].copy()
            args['mfile'] = self.pargs['traj']['mfile']
            args['mName'] = name
            
            dname = 'dump' + str(np.random.randint(0,1e8))
            self.pargs['traj']['dump_mname'] = [dname]

            self.lmp.command('dump ' + dname + ' all mesh/vtk {freq} {dir}/{mfile} id stress stresscomponents vel '.format(**args) + name)

    self.dump = True

    return name

  def extractCoords(self):
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

    coords = np.zeros((self.lmp.get_natoms(),3))

    for i in range(self.lmp.get_natoms()):
      coords[i,:] = x[i], y[i], z[i]

    self.lmp.command('variable x delete')
    self.lmp.command('variable y delete')
    self.lmp.command('variable z delete')

    return coords

  def monitor(self, name, group, var, file, species='all'):
    """
    """
    self.lmp.command('compute {} {} {}'.format(var, species, name))
    self.lmp.command('fix my{} {} ave/time 1 1 1 c_{} file {}'.format(var, group, var, file))

  def add_viscous(self, **args):
    """ Adds a viscous damping force: F = - gamma * v for each particle
    @species = 1,2, ... or all
    @gamma: real number (viscosity coefficient)
    @[scale]: tuple (species, ratio) to scale gamma with
    """
    if 'scale' not in args:
      args['scale'] = (args['species'], 1)

    if 'species' in args:
      if isinstance(args['species'], int):
        args['species'] = 'group' + str(args['species']-1)
    else:
      raise RuntimeError('Species must be specified (1,2,..., or "all") for which the viscous force applies.')

    name = np.random.randint(0, 1e8)

    self.lmp.command('fix {}'.format(name) + ' {species} viscous {gamma} scale '.format(**args) + (' {} ' * len(args['scale'])).format(*args['scale']) )

    return name

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
        raise("Unexpected error:", sys.exc_info()[0])

  def saveas(self, name, fname):
    """
    """
    if not self.rank:

      try:
        np.savetxt(fname, np.array(self.vars[name]))
      except:
         raise("Unexpected error:", sys.exc_info()[0])

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
    pass
