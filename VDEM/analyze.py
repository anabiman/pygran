import numpy as np

"""
TODO: compute mass flow rate, density, nns distribution (radial distribution function)
"""

def computeFlow(data, sel):
	"""
	Computes flow rate: N/t for a selection *sel*
	@ data: list of dictionaries containing simulation and particle data (box size, x,y,z, etc.)
	"""
	
	if data['TIMESTEP'] <= 0:
		return 0
	else:
		return np.float(len(sel)) / data['TIMESTEP']

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

def readCustomTraj(fname, flow = False, region = ()):
	"""
	transforms a LAMMPS/LIGGGHTS custom dump file(s) to a python trajectory
	"""
	dicList = []

	if flow:
		flowRate = []

	with open(fname,'r') as fp:
		while True:

			dic = {}

			while True:
				line = fp.readline()

				if not line: 

					if flow:
						flowRate = np.array(flowRate)
						np.savetxt('flow.dat', flowRate)

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

			if flow:
				sel = select(dic)
				flowRate.append(computeFlow(dic, sel))

			dicList.append(dic)

	if flow:
		flowRate = np.array(flowRate)
		np.savetxt('flow.dat', flowRate)

	return dicList