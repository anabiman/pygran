import numpy as np

"""
TODO: compute mass flow rate, density, nns distribution (radial distribution function)
- Compute BEVERLO coefficients
- 

"""
def computeAngleRepos(data, *args):
	"""
	Computes the angle of repos theta = arctan(h_max/L)
	in a sim box defined by [-Lx, Lx] x [-Ly, Ly] x [0, Lz]
	"""
	Lx, Ly = args
	x, y, z = data['x'], data['y'], data['z']
	r = data['radius'].mean()

	h_max = z.max()

	# Select all particles close to the walls (within r distance)
	z = z[x**2.0 + y**2.0 >= (0.5 * (Lx + Ly) - r)**2.0]

	print len(z)

	if len(z):
		zm = z.max()

		dzMin = zm * 0.9
		dzMax = zm 

		z = z[z >= dzMin]
		z = z[z <= dzMax]
		h_max -= z.mean()

		print np.arctan(h_max / Lx) * 180.0 / np.pi
		return np.arctan(h_max / Lx)
	else:
		return 0

def computeFlow(data, density, t0 = 0, N0 = 0, sel = None, dt = 1e-4):
	"""
	Computes flow rate: N/t for a selection *sel*
	@ data: list of dictionaries containing simulation and particle data (box size, x,y,z, etc.)

	TODO: Get this working for a particle distribution
	"""

	if N0 == None or t0 == None:
		return 0
	else:
		mass = density * 4.0 / 3.0 * np.pi * (len(sel) - N0) * np.mean(data['radius'][sel])**3.0
		return - mass / ((data['TIMESTEP'] - t0) * dt)

def computeDensity(data, density, shape = 'box', sel = None):
	"""
	Computes the bulk density for a selection of particles from their true *density*. 
	The volume is determined approximately by constructing a box/cylinder/cone 
	embedding the particles. Particles are assumed to be spherical in shape.
	"""
	x, y, z = data['x'][sel], data['y'][sel], data['z'][sel]

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

	mass = np.sum(density * 4.0 / 3.0 * np.pi * (data['radius'][sel])**3.0)

	return mass / volume

def computeHeight(data, axis):
	"""
	Computes the mean max height of an N-particle system along the x, y, or z axis.
	"""
	height = data[axis].max()
	hmin = height * 0.99
	hmax = height * 1.01

	if axis == 'x':
		region = (hmin, hmax, -np.inf, np.inf, -np.inf, np.inf)
	elif axis == 'y':
		region = (-np.inf, np.inf, hmin, hmax, -np.inf, np.inf)
	elif axis == 'z':
		region = (-np.inf, np.inf, -np.inf, np.inf, hmin, hmax)
	else:
		print "axis must be x, y, or z"
		raise

	return data[axis][select(data, *region)].mean()

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

def readCustomTraj(fname, height = False, flow = False, density = None, shape = None, angle = None, region = (), dt = 1e-4):
	"""
	transforms a LAMMPS/LIGGGHTS custom dump file(s) to a python trajectory
	"""
	dicList = [] # Keep this empty for now to avoid large memory buffers
	t0, N0 = None, None

	if flow:
		flowRate = []

	if shape:
		bDensity = []

	if angle:
		rAngle = []

	if height:
		hdist = []
		axis = shape[-1]

	with open(fname,'r') as fp:
		while True:

			dic = {}

			while True:
				line = fp.readline()

				if not line: 

					if flow:
						flowRate = np.array(flowRate)
						np.savetxt('flow.dat', flowRate)

					if shape:
						bDensity = np.array(bDensity)
						np.savetxt('density.dat', bDensity)

					if height:
						hdist = np.array(hdist)
						np.savetxt('height.dat', hdist)

					if angle:
						rAngle = np.array(rAngle)
						np.savetxt('angle.dat', rAngle)

					return dicList

				if line.find('TIMESTEP') >= 0:
					timestep = int(fp.readline())
					dic['TIMESTEP'] = timestep

				if line.find('NUMBER OF ATOMS') >= 0:
					natoms = int(fp.readline())
					dic['NATOMS'] = natoms

				if line.find('BOX') >= 0:
					boxX = fp.readline().split()
					boxY = fp.readline().split()
					boxZ = fp.readline().split()

					boxX = [float(i) for i in boxX]
					boxY = [float(i) for i in boxY]
					boxZ = [float(i) for i in boxZ]

					dic['BOX'] = (boxX, boxY, boxZ)

				if line.find('ITEM: ATOMS') >= 0:
					break

			keys = line.split()[2:] # remove ITEM: and ATOMS keywords

			for key in keys:
				dic[key] = np.zeros(natoms)

			for i in range(natoms):
				var = fp.readline().split()

				for j, key in enumerate(keys):
					dic[key][i] = float(var[j]) 

			sel = select(dic, *region)

			if len(region):
				print "Found {} particles based on user-supplied selection".format(len(sel))

			if len(sel): # Make sure the frame contains NATOMS > 0 else dont bother compute anything

				if flow:
					flowRate.append(computeFlow(dic, density, t0, N0, sel, dt))

				if shape:
					bDensity.append(computeDensity(dic, density, shape, sel))

				if height:
					hdist.append(computeHeight(dic, axis))

				if angle:
					rAngle.append(computeAngleRepos(dic, *angle))

			t0 = timestep
			N0 = len(sel)

	if flow:
		flowRate = np.array(flowRate)
		np.savetxt('flow.dat', flowRate)

	if shape:
		bDensity = np.array(bDensity)
		np.savetxt('density.dat', bDensity)

	if height:
		hdist = np.array(hdist)
		np.savetxt('height.dat', hdist)

	if angle:
		rAngle = np.array(rAngle)
		np.savetxt('angle.dat', rAngle)