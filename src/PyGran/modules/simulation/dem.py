# !/usr/bin/python
# -*- coding: utf8 -*- 
# -----------------------------------------------------------------------
#
#   Python interface for running DEM simulations
#
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
Created on April 25, 2016
@author: Andrew Abi-Mansour
'''

from mpi4py import MPI
from importlib import import_module
from datetime import datetime
import os, sys
from PyGran.modules.tools import find, _setConfig
from . import models
import shutil
import PyGran

class DEM:
  """A *generic* class that handles communication for a DEM object independent of the engine used"""

  def __init__(self, **pargs):
    """ Initializes COMM and partition proccesors based on user input """

    # Instantiate contact model and store it in pargs
    if 'model' not in pargs:
      pargs['model'] = models.SpringDashpot

    pargs = pargs['model'](**pargs).params

    self.comm = MPI.COMM_WORLD
    self.rank = self.comm.Get_rank()
    self.tProcs = self.comm.Get_size()
    self.nSim = pargs['nSim']
    self.model = str(pargs['model']).split("'")[1].split('.')[-1]
    self.pargs = pargs
    self.library = None
    self._dir,_ = os.path.abspath(__file__).split(os.path.basename(__file__))
    
    # Check if .config files eixsts else create it
    # Only one process needs to do this
    if not self.rank:
      self.library, src, version = _setConfig(wdir=self._dir, engine=self.pargs['engine'].split('engine_')[1])

      if version:
        self.pargs['__version__'] = version

      if src:
        self.pargs['liggghts_src'] = src

      for slave in range(1,self.tProcs):
          self.comm.send(self.library, dest=slave, tag=0)
          self.comm.send('__version__' in self.pargs, dest=slave, tag=1)

          if '__version__' in self.pargs:
            self.comm.send(self.pargs['__version__'], dest=slave, tag=2)
    else:
      self.library = self.comm.recv(source=0, tag=0)
      if self.comm.recv(source=0, tag=1):
        self.pargs['__version__'] = self.comm.recv(source=0, tag=2)

    if 'output' not in self.pargs:
      # The idea is to create a unique output name that depends on the current time. Since the processes are not in sunc, it's safer
      # to create the output name on the master processor and then send it to the slaves.
      if not self.rank:
        time = datetime.now()
        self.pargs['output'] = 'out-{}-{}:{}:{}-{}.{}.{}'.format(self.model, time.hour, time.minute, time.second, time.day, time.month, time.year)
        
        for slave in range(1,self.tProcs):
          self.comm.send(self.pargs['output'], dest=slave, tag=3)

      else:
        self.pargs['output'] = self.comm.recv(source=0, tag=3)

    if self.nSim > self.tProcs:
      print("Number of simulations ({}) cannot exceed number of available processors ({})".format(self.nSim, self.tProcs))
      sys.exit(0)

    self.pProcs = self.tProcs // self.nSim

    for i in range(self.nSim):
      if self.rank < self.pProcs * (i + 1):
        self.color = i
        break
      else:
        # In case of odd number of procs, place the one left on the last communicator 
        self.color = self.nSim

    self.split = self.comm.Split(color=self.color, key=self.rank)

    # update rank locally for each comm
    self.rank = self.split.Get_rank()

    module = import_module('PyGran.modules.simulation.' + self.pargs['engine'])
    output = self.pargs['output'] if self.nSim == 1 else (self.pargs['output'] + '-multi-mode-' + str(self.color))

    if not self.split.Get_rank():
      if os.path.exists(output):
        print('WARNING: output dir {} already exists. Proceeding ...'.format(output))
      else:
        os.mkdir(output)

    # Make sure output in self.pargs is updated before instantiating dem class
    self.pargs['output'] = output

    self.split.barrier() # Synchronize all procs
        
    os.chdir(self.pargs['output'])

    self.dem = module.DEMPy(self.split, self.library, **self.pargs) # logging module imported here  

    if not self.rank:
      global logging

      logging = import_module(name='logging')
      logging.basicConfig(filename='dem.log', format='%(asctime)s:%(levelname)s: %(message)s', level=logging.DEBUG)

      logging.info('Initializing simulation with PyGran version {}'.format(PyGran.__version__))
      logging.info("Initializing MPI for a total of {} cores".format(self.tProcs))

      if self.nSim > 1:
        logging.info('Running {} simulations: multi-mode on'.format(self.nSim))

      if self.pProcs > 0:
        logging.info('Using {} cores per simulation'.format(self.pProcs))

      from sys import argv
      scriptFile = argv[0]

      if not scriptFile.endswith('__main__.py'): # user is importing their script as a module, dont back up:
        if scriptFile.endswith('.py'):
            logging.info('Backing up {} file'.format(scriptFile))
            shutil.copyfile('{}/../{}'.format(os.getcwd(), scriptFile), '{}'.format(scriptFile.split('.')[0] + '-bk.py'))
      else:
        logging.info('Input script run as a module. Not backing up file')

    # All I/O done ~ phew! Now initialize DEM
    # Import and setup all meshes as rigid walls
    self.initialize()

    # Setup material properties
    if 'materials' in self.pargs:
      for item in self.pargs['materials'].keys():
        # Overloaded function 'createProperty' will partition material propreties based on MPI's coloring split scheme
        # Do we even need this?
        if isinstance(self.pargs['materials'][item], tuple): # Make sure we're not reading user-defined scalars (e.g. density)
          self.createProperty(item, *self.pargs['materials'][item])

    self.printSetup()

    # Create links to the particle/mesh files (easily accessible to the user)
    if 'pfile' in self.pargs['traj']:
      if self.pargs['traj']['pfile']:
        self.pfile = self.pargs['output'] + '/traj/' + self.pargs['traj']['pfile']
      else:
        self.pfile = None

    if 'mfile' in self.pargs['traj']:
      if self.pargs['traj']['mfile']:
        self.mfile = self.pargs['output'] + '/traj/' + self.pargs['traj']['mfile']
      else:
        self.mfile = None

  def scatter_atoms(self,name,type,count,data):
    for i in range(self.nSim):
      if self.rank < self.pProcs * (i + 1):
        return self.dem.scatter_atoms(name,type,count,data)

  def createParticles(self, type, style, *args):
    for i in range(self.nSim):
      if self.rank < self.pProcs * (i + 1):
        self.dem.createParticles(type, style, *args)
        break

  def createGroup(self, *group):
    for i in range(self.nSim):
      if self.rank < self.pProcs * (i + 1):
        self.dem.createGroup(*group)
        break

  def set(self, *args):
    for i in range(self.nSim):
      if self.rank < self.pProcs * (i + 1):
        self.dem.set(*args)
        break

  def gather_atoms(self,name,type,count):
    for i in range(self.nSim):
      if self.rank < self.pProcs * (i + 1):
        return self.dem.gather_atoms(name,type,count)

  def get_natoms(self):
    for i in range(self.nSim):
      if self.rank < self.pProcs * (i + 1):
        return self.dem.get_natoms()

  def extract_global(self,name,type):
    for i in range(self.nSim):
      if self.rank < self.pProcs * (i + 1):
        return self.dem.extract_global(name, type)
        
  def extract_compute(self,id,style,type):
    for i in range(self.nSim):
      if self.rank < self.pProcs * (i + 1):
        return self.dem.extract_compute(id,style,type)

  def extract_fix(self,id,style,type,i=0,j=0):
    for i in range(self.nSim):
      if self.rank < self.pProcs * (i + 1):
        return self.dem.extract_fix(id,style,type,i,j)

  def initialize(self):
    for i in range(self.nSim):
      if self.rank < self.pProcs * (i + 1):
        self.dem.initialize()
        break
  
  def velocity(self, *args):
      for i in range(self.nSim):
        if self.rank < self.pProcs * (i + 1):
          self.dem.velocity(*args)
          break

  def add_viscous(self, **args):
    """ Adds a viscous damping force: F = - gamma * v for each particle
    @species = 1,2, ... or all
    @gamma: real number (viscosity coefficient)
    @[scale]: tuple (species, ratio) to scale gamma with
    """
    for i in range(self.nSim):
        if self.rank < self.pProcs * (i + 1):
          return self.dem.add_viscous(**args)

  def insert(self, species, value, **args):
    for i in range(self.nSim):
      if self.rank < self.pProcs * (i + 1):
        return self.dem.insert(species, value, **args)

  def run(self, nsteps, dt=None, itype=None):
    for i in range(self.nSim):
      if self.rank < self.pProcs * (i + 1):
        return self.dem.run(nsteps, dt, itype)
        
  def setupParticles(self):

    for i in range(self.nSim):
      if self.rank < self.pProcs * (i + 1):
        self.dem.setupParticles()
        break

  def createProperty(self, name, *args):
    """
    Material and interaction properties required
    """

    for i in range(self.nSim):
      if self.rank < self.pProcs * (i + 1):
        if type(args[0]) is tuple:
          self.dem.createProperty(name, *args[i])
        else:
          self.dem.createProperty(name, *args)
        break

  def importMeshes(self, name=None):
    """
    An internal function that is called during DEM initialization for importing meshes.
    This function imports all meshes and sets them up as walls. Can import only one mesh specified by the 'name' keyword.
    """
    for i in range(self.nSim):
      if self.rank < self.pProcs * (i + 1):
        self.dem.importMeshes(name)
        break

  def importMesh(self, name, file, mtype, **args):
    """
    Imports a mesh file (STL or VTK)
    """
    for i in range(self.nSim):
      if self.rank < self.pProcs * (i + 1):
        self.dem.importMesh(name, file, mtype, **args)
        break

  def setupWall(self, wtype, species = None, plane = None, peq = None):
    """
    Creates a wall
    @ name: name of the variable defining a wall or a mesh
    @ wtype: type of the wall (primitive or mesh)
    @ plane: x, y, or z plane for primitive walls
    @ peq: plane equation for primitive walls
    """
    for i in range(self.nSim):
      if self.rank < self.pProcs * (i + 1):
        return self.dem.setupWall(wtype, species, plane, peq)

  def printSetup(self):
    """
    Specify which variables to write to file, and their format
    """
    for i in range(self.nSim):
      if self.rank < self.pProcs * (i + 1):
        self.dem.printSetup()
        break

  def writeSetup(self, only_mesh=False, name=None):
    """
    """
    for i in range(self.nSim):
      if self.rank < self.pProcs * (i + 1):
        dumpID = self.dem.writeSetup(only_mesh, name)

        # Create or update links to the particle/mesh files (easily accessible to the user)
        if 'pfile' in self.lmp.pargs['traj']:
          self.pfile = self.lmp.pargs['output'] + '/traj/' + self.lmp.pargs['traj']['pfile']

        if 'mfile' in self.lmp.pargs['traj']:
          self.mfile = self.lmp.pargs['output'] + '/traj/' + self.lmp.pargs['traj']['mfile']

        return dumpID

  def setupIntegrate(self, name, itype='nve/sphere', group='all'):
    """
    Specify how Newton's eqs are integrated in time. 
    @ name: name of the fixed simulation ensemble applied to all atoms
    @ dt: timestep
    @ ensemble: ensemble type (nvt, nve, or npt)
    @ args: tuple args for npt or nvt simulations
    """
    for i in range(self.nSim):
      if self.rank < self.pProcs * (i + 1):
        self.dem.setupIntegrate(name, itype, group)
        break

  def integrate(self, steps, dt, itype=None):
    """
    Run simulation in time
    """
    for i in range(self.nSim):
      if self.rank < self.pProcs * (i + 1):
        self.dem.integrate(steps, dt, itype)
        break

  def remove(self, name):
    """
    Delete variable/object
    """
    for i in range(self.nSim):
      if self.rank < self.pProcs * (i + 1):
        self.dem.remove(name)
        break

  def monitor(self, **args):
    """
    """
    for i in range(self.nSim):
      if self.rank < self.pProcs * (i + 1):
        return self.dem.monitor(**args)

  def plot(self, fname, xlabel, ylabel, output=None, xscale=None):
    """
    """
    for i in range(self.nSim):
      if self.rank < self.pProcs * (i + 1):
        self.dem.plot(fname, xlabel, ylabel, output, xscale)
        break

  def moveMesh(self, name, **args):
    """
    """
    for i in range(self.nSim):
      if self.rank < self.pProcs * (i + 1):
        return self.dem.moveMesh(name, **args)

  def saveas(self, name, fname):
    """
    """
    for i in range(self.nSim):
      if self.rank < self.pProcs * (i + 1):
        self.dem.saveas(name, fname)
        break

  def command(self, cmd):
    """
    """
    for i in range(self.nSim):
      if self.rank < self.pProcs * (i + 1):
        self.dem.command(cmd)
        break

  def close(self):
    """
    Deletes allocated memory and changes directory back to cwd. 
    """
    # Dont call this since the user might be running multiple simulations in one script
    #MPI.Finalize()
    for i in range(self.nSim):
      if self.rank < self.pProcs * (i + 1):
        self.dem.close()
        break

    os.chdir('..')    

  def __enter__(self):
    return self

  def __exit__(self, exc_type, exc_value, traceback):
    os.chdir('..')
