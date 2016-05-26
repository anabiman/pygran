import numpy as np

"""
TODO: compute mass flow rate, density, nns distribution (radial distribution function)
"""

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

def readCustomTraj(fname):
	"""
	transforms a LAMMPS/LIGGGHTS custom dump file(s) to a python trajectory
	"""
	dicList = []

	with open(fname,'r') as fp:
		while True:

			dic = {}

			while True:
				line = fp.readline()

				if not line: 
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

			dicList.append(dic)

	return dicList