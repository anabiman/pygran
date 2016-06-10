import numpy as np

"""
TODO: compute mass flow rate, density, nns distribution (radial distribution function)
- Compute BEVERLO coefficients
- 

"""

def computeFlow(data, t0 = 0, N0 = 0, sel = None):
	"""
	Computes flow rate: N/t for a selection *sel*
	@ data: list of dictionaries containing simulation and particle data (box size, x,y,z, etc.)
	"""
	
	if data['TIMESTEP'] - t0 <= 0:
		return 0
	else:
		return - (np.float(len(sel)) - N0) / (data['TIMESTEP'] - t0)

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

	elif shape == 'cylinderZ':
		height = zmax - zmin
		radius = (ymax - ymin) * 0.25 + (xmax - xmin) * 0.25
		volume = np.pi * radius**2.0 * height

	elif shape == 'cylinderY':
		height = ymax - ymin
		radius = (zmax - zmin) * 0.25 + (xmax - xmin) * 0.25
		volume = np.pi * radius**2.0 * height

	elif shape == 'cylinderX':
		height = xmax - xmin
		radius = (ymax - ymin) * 0.25 + (zmax - zmin) * 0.25
		volume = np.pi * radius**2.0 * height

	mass = np.sum(density * 4.0 / 3.0 * np.pi * (data['radius'][sel])**3.0)

	return mass / volume

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
		indices = np.where(x > xmin)
		indices = indices + np.where(x < xmax)
		indices = np.unique(indices)

		indices = indices + np.where(y > ymin)
		indices = indices + np.where(y < ymax)
		indices = np.unique(indices)

		indices = indices + np.where(z > zmin)
		indices = indices + np.where(z < zmax)
		indices = np.unique(indices)

		return indices

	except:
		raise

def readCustomTraj(fname, flow = False, density = None, shape = 'box', region = ()):
	"""
	transforms a LAMMPS/LIGGGHTS custom dump file(s) to a python trajectory
	"""
	dicList = []
	foundT0 = False
	foundN0 = False

	if flow:
		flowRate = []

	if density:
		bDensity = []

	with open(fname,'r') as fp:
		while True:

			dic = {}

			while True:
				line = fp.readline()

				if not line: 

					if flow:
						flowRate = np.array(flowRate)
						np.savetxt('flow.dat', flowRate)

					if density:
						bDensity = np.array(bDensity)
						np.savetxt('density.dat', bDensity)

					return dicList

				if line.find('TIMESTEP') >= 0:
					timestep = int(fp.readline())
					dic['TIMESTEP'] = timestep

					if not foundT0:
						t0 = timestep
						foundT0 = True

				if line.find('NUMBER OF ATOMS') >= 0:
					natoms = int(fp.readline())
					dic['NATOMS'] = natoms

					if not foundN0:
						N0 = natoms
						foundN0 = True

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

			sel = select(dic)
			
			if flow:
				flowRate.append(computeFlow(dic, t0, N0, sel))

			if density:
				bDensity.append(computeDensity(dic, density, shape, sel))

			dicList.append(dic)

	if flow:
		flowRate = np.array(flowRate)
		np.savetxt('flow.dat', flowRate)

	if density:
		bDensity = np.array(bDensity)
		np.savetxt('density.dat', bDensity)

	return dicList