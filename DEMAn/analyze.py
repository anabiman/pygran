import numpy as np
from xlrd import open_workbook
from numbers import Number
from numpy import zeros, sqrt, where, pi, mean, arange, histogram

"""
TODO: compute mass flow rate, density, nns distribution (radial distribution function)
- Compute BEVERLO coefficients
- 

"""

def pairCorrelationFunction(x, y, z, S, rMax, dr):
    """Compute the three-dimensional pair correlation function for a set of
    spherical particles contained in a cube with side length S.  This simple
    function finds reference particles such that a sphere of radius rMax drawn
    around the particle will fit entirely within the cube, eliminating the need
    to compensate for edge effects.  If no such particles exist, an error is
    returned.  Try a smaller rMax...or write some code to handle edge effects! ;)

    Arguments:
        x               an array of x positions of centers of particles
        y               an array of y positions of centers of particles
        z               an array of z positions of centers of particles
        S               length of each side of the cube in space
        rMax            outer diameter of largest spherical shell
        dr              increment for increasing radius of spherical shell

    Returns a tuple: (g, radii, interior_indices)
        g(r)            a numpy array containing the correlation function g(r)
        radii           a numpy array containing the radii of the
                        spherical shells used to compute g(r)
        reference_indices   indices of reference particles
    """
    
    # center positions around 0
    x -= x.min()
    y -= y.min()
    z -= z.min()

    # Find particles which are close enough to the cube center that a sphere of radius
    # rMax will not cross any face of the cube
    bools1 = x > rMax
    bools2 = x < (S - rMax)
    bools3 = y > rMax
    bools4 = y < (S - rMax)
    bools5 = z > rMax
    bools6 = z < (S - rMax)

    interior_indices, = where(bools1 * bools2 * bools3 * bools4 * bools5 * bools6)
    num_interior_particles = len(interior_indices)

    if num_interior_particles < 1:
        raise  RuntimeError ("No particles found for which a sphere of radius rMax\
        	will lie entirely within a cube of side length S.  Decrease rMax\
                or increase the size of the cube.")

    edges = arange(0., rMax + 1.1 * dr, dr)
    num_increments = len(edges) - 1
    g = zeros([num_interior_particles, num_increments])
    radii = zeros(num_increments)
    numberDensity = len(x) / S**3

    # Compute pairwise correlation for each interior particle
    for p in range(num_interior_particles):
        index = interior_indices[p]
        d = sqrt((x[index] - x)**2 + (y[index] - y)**2 + (z[index] - z)**2)
        d[index] = 2 * rMax

        (result, bins) = histogram(d, bins=edges, normed=False)
        g[p,:] = result / numberDensity

    # Average g(r) for all interior particles and compute radii
    g_average = zeros(num_increments)
    for i in range(num_increments):
        radii[i] = (edges[i] + edges[i+1]) / 2.
        rOuter = edges[i + 1]
        rInner = edges[i]
        g_average[i] = mean(g[:, i]) / (4.0 / 3.0 * pi * (rOuter**3 - rInner**3))

    return (g_average, radii, interior_indices)
    # Number of particles in shell/total number of particles/volume of shell/number density
    # shell volume = 4/3*pi(r_outer**3-r_inner**3)


def readExcel(fname):

	"""
	reads an excel sheet(s) and appends each data column to 
	a dictionary
	"""
	book = open_workbook(fname, on_demand=True)
	data = {}

	for name in book.sheet_names():
		sheet = book.sheet_by_name(name)

		for coln in range(sheet.ncols):
			for celln, cell in enumerate(sheet.col(coln)):
				if celln == 0:
					dname = cell.value
					data[dname] = []
				else:
					if isinstance(cell.value, Number):
						data[dname].append(cell.value)
					else:
						print cell.value, ' ignored'

	return data


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