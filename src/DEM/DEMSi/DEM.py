from mpi4py import MPI
from importlib import import_module

class DEM:
  """A class that handles communication for the DEM object"""

  def __init__(self, **pargs):
    """ Initializes COMM and partition proccesors based on user input """

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


        module = import_module(pargs['modName'])

        self.dem = module.DEMPy(i, self.split, **pargs) # logging module imported here      
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
  
  def insertParticles(self, name, *args):

    for i in range(self.nSim):
      if self.rank < self.nPart * (i + 1):
        return self.dem.insertParticles(name, *args)

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

    MPI.Finalize()