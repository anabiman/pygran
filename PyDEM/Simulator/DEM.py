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
#   the Free Software Foundation, either version 3 of the License, or
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
import os

def find(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)

class DEM:
  """A *generic* class that handles communication for a DEM object independent of the engine used"""

  def __init__(self, **pargs):
    """ Initializes COMM and partition proccesors based on user input """

    self.comm = MPI.COMM_WORLD
    self.rank = self.comm.Get_rank()
    self.tProcs = self.comm.Get_size()
    self.nSim = pargs['nSim']
    self.model = str(pargs['model']).split('.')[-1]
    self.pargs = pargs
    self.library = None
    self._dir, _ = __file__.split('DEM.py')
    self.pargs['output'] = None
    
    # Check if .config files eixsts else create it
    # Only one process needs to do this
    if not self.rank:
      if os.path.isfile(self._dir + '../.config'):
        with open(self._dir + '../.config', 'r+') as fp:
          self.library = fp.read().split('=')[-1]

          # Make sure the library exists; else, find it somewhere else 
          if not os.path.isfile(self.library):
            self.library = find('lib' + self.pargs['engine'].split('engine_')[1] + '.so', '/')
            fp.seek(0,0)
            fp.write('library=' + self.library)
            print 'WARNING: Could not find user-specified library. Will use {} instead ...'.format(self.library)
      else:
        with open(self._dir + '../.config', 'w') as fp:
          self.library = find('lib' + self.pargs['engine'].split('engine_')[1] + '.so', '/')
          print 'WARNING: No config file found. Creating one for {}'.format(self.library)
          fp.write('library=' + self.library)

      for slave in range(1,self.tProcs):
          self.comm.send(self.library, dest=slave)
    else:
      self.library = self.comm.recv(source=0)

    if 'out' not in self.pargs:
      # The idea is to create a unique output name that depends on the current time. Since the processes are not in sunc, it's safer
      # to create the output name on the master processor and then send it to the slaves.
      if not self.rank:
        time = datetime.now()
        self.pargs['output'] = 'out-{}-{}:{}:{}-{}.{}.{}'.format(self.model, time.hour, time.minute, time.second, time.day, time.month, time.year)
        
        for slave in range(1,self.tProcs):
          self.comm.send(self.pargs['output'], dest=slave)

      else:
        self.pargs['output'] = self.comm.recv(source=0)

    if self.nSim > self.tProcs:
      print "Number of simulations ({}) cannot exceed number of available processors ({})".format(self.nSim, self.tProcs)
      sys.exit(0)

    self.nPart = self.tProcs // self.nSim

    for i in range(self.nSim):
      if self.rank < self.nPart * (i + 1):
        self.color = i
        self.split = self.comm.Split(self.color, key=0)

        module = import_module('PyDEM.Simulator.' + self.pargs['engine'])
        self.output = self.pargs['output'] if self.nSim == 1 else (self.pargs['output'] + '{}'.format(i))

        if not self.split.Get_rank():
          if os.path.exists(self.output):
            print 'WARNING: output dir {} already exists. Proceeding ...'.format(self.output)
          else:
            os.mkdir(self.output)

        self.split.barrier() # Synchronize all procs
        self.dem = module.DEMPy(i, self.split, self.library, **self.pargs) # logging module imported here      
        break

    if not self.rank:
      global logging

      logging = import_module(name='logging')
      logging.basicConfig(filename='dem.log', format='%(asctime)s:%(levelname)s: %(message)s', level=logging.DEBUG)

      logging.info("Initializing MPI for a total of {} procs".format(self.tProcs))

      if self.nSim > 1:
        logging.info('Running {} simulations: multi-mode on'.format(self.nSim))

    # All I/O done ~ phew! Now initialize DEM
    self.initialize()

    # Setup material properties
    if 'materials' in self.pargs:
      for item in self.pargs['materials'].keys():
        # Overloaded function 'createProperty' will partition coeffRest based on MPI's coloring split scheme
        self.createProperty(item, *self.pargs['materials'][item])

    # Import and setup all meshes as rigid walls
    if 'mesh' in self.pargs:
      for mesh in self.pargs['mesh'].keys():
        self.importMesh(name=mesh, **self.pargs['mesh'][mesh])
        self.setupWall(name=mesh + 'Wall', wtype='mesh', meshName=mesh)

    self.printSetup()

    # Write output to trajectory by default unless the user specifies otherwise
    if 'dump' in self.pargs:
      if self.pargs['dump'] == True:
        self.dumpSetup()
    else:
      self.dumpSetup()

  def initialize(self):

    for i in range(self.nSim):
      if self.rank < self.nPart * (i + 1):
        self.dem.initialize()
        break
  
  def velocity(self, *args):
      for i in range(self.nSim):
        if self.rank < self.nPart * (i + 1):
          self.dem.velocity(*args)
          break

  def insert(self, name, *args):
    for i in range(self.nSim):
      if self.rank < self.nPart * (i + 1):
        self.dem.insert(name, *args)
        break

  def run(self, nsteps=None, dt=None):
    for i in range(self.nSim):
      if self.rank < self.nPart * (i + 1):
        self.dem.run(nsteps, dt)
        break
        
  def setupParticles(self):

    for i in range(self.nSim):
      if self.rank < self.nPart * (i + 1):
        self.dem.setupParticles()
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

  def setupWall(self, name, wtype, meshName = None, plane = None, peq = None):
    """
    Creates a wall
    @ name: name of the variable defining a wall or a mesh
    @ wtype: type of the wall (primitive or mesh)
    @ plane: x, y, or z plane for primitive walls
    @ peq: plane equation for primitive walls
    """
    for i in range(self.nSim):
      if self.rank < self.nPart * (i + 1):
        self.dem.setupWall(name, wtype, meshName, plane, peq)
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
    """
    """
    MPI.Finalize()