from liggghts import liggghts

class DEM:

	def __init__(self, units = si, dim = 3, **args):
		""" Initialize some settings and specifications 
		@ lmp: lammps object
		@ units: unit system (si, cgs, etc.)
		@ dim: dimensions of the problems (2 or 3)
		"""

		self.lmp = liggghts(**args)
		self.lmp.command('units {}'.format(units))
		self.lmp.command('dimension {}'.format(dim))
		self.lmp.command('atom_style granular')
		self.lmp.command('atom_modify map array') # array is faster than hash in looking up atomic IDs, but the former takes more memory
		self.lmp.command('boundary f f f')
		self.lmp.command('newton off')
		self.lmp.command('processors * * *') # let LAMMPS handle DD

	def createDomain(self, nss, *pos):
		""" Define the domain of the simulation
		@ lmp: lammps object
		@ nsys: number of subsystems
		@ pos: 6 x 1 tuple that defines the boundaries of the box 
		"""

		self.lmp.command('region domain block {} units box'.format(pos)
		self.lmp.command('create_box {} domain'.format(nss))

	def insertParticle(self, N, density = 1000, *vel, **args):
		""" Insert particles in a pre-defined region
		@ lmp: lammps object
		@ N: max total number of particles to be inserted
		@ density: initial density of the particles
		@ vel: 3 x 1 tuple of initial velocities of all particles
		@ args: dictionary of params
		"""
	
		radius = args['radius']

		if radius is 'uniform number':
			radiusMin = args['radiusMin']
			radiusMax = args['radius']

		self.lmp.command('fix {} all particletemplate/sphere 1 atom_type 1 density constant {} radius {} {} {}'.format(var1, density, radius, radiusMin, radiusMax))

		self.lmp.command('fix {} all particledistribution/discrete 63243 1 pts 1.0'.format(var2))
		self.lmp.command('region factory sphere 0 1.5 0 0.5 units box'.format(var3))
		self.lmp.command('fix {} all insert/rate/region seed 123481 distributiontemplate pdd & \
			nparticles {} particlerate {} insert_every {} & \
			overlapcheck yes vel constant {} region factory ntry_mc 10000'.format(var4, N, args['particleRate'], args['freq'], vel))

	def setupNeighbor(self):
		"""
		"""
		self.lmp.command('neighbor 0.03 bin')
		self.lmp.command('neigh_modify delay 0')

	def createProperty(self, var, property, type, valueProp, valueType):
		"""
		Material and interaction properties required
		"""
		self.lmp.command('fix {} all property/global {} {} {} {}'.format(var, property, type, valueProp, valueType))

	def importMesh(self, var, fname):
		"""
		"""
		self.lmp.command('fix {} all mesh/surface file {} type 2 scale {}'.format(var, fname))

	def setupWall(self, var, wtype, mesh = None, plane = None, peq = None):
		"""
		Use the imported mesh as granular wall
		"""

		if wtype == 'mesh':
			lmp.command('fix {} all wall/gran model hertz tangential history {} n_meshes 1 meshes {}'.format(var, wtype, mesh))
		elif wtype == 'primitive':
			lmp.command('fix {} all wall/gran model hertz tangential history {} type 2 {} {}'.format(var, wtype, plane, peq))
		else:
			raise ValueError('Wall type can be either primitive or mesh')

	def setupPhysics(self):
		"""
		Specify the interation forces
		"""
		self.lmp.command('pair_style gran model hertz tangential history')
		self.lmp.command('pair_coeff * *')

	def setupGravity(self, var, *gravity):
		"""
		Specify in which direction the gravitational force acts
		"""
		self.lmp.command('fix {} all gravity 9.81 vector {}'.format(var, gravity))

	def integrate(self, var, dt):
		"""
		Integrate Newton's eqs in time. 
		TODO: find out if we can do NVT simulations
		"""
		self.lmp.command('fix {} all nve/sphere'.format(var))
		self.lmp.command('timestep {}'.format(dt))

	def printSetup(self, freq, *args=('step', 'atoms', 'ke')):
		"""
		Specify which variables to write to file, and their format
		"""
		self.lmp.command('thermo_style custom {}'.format(args))
		self.lmp.command('thermo {}'.format(freq))
		self.lmp.command('thermo_modify norm no lost ignore')
	