import numpy as np

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