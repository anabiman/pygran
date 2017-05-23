from numbers import Number
from PIL import Image
from numpy import random, array, linspace, sqrt, fabs

def readExcel(fname):
	"""
	reads an excel sheet(s) and appends each data column to 
	a dictionary
	"""
	from xlrd import open_workbook

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

	for key in data.keys():
		data[key] = array(data[key])

	return data

def genImg(Particles, zmin, zmax, dz, output=None, imgShow = None):
	"""
	Generates a 2D image from a slice (limited by 'zmin/zmax' and of thickness 'dz') 
	of a 3D config in the Particles class. The scale is the number of microns per pixel.
	"""
	scale = 1.0/dz

	maxRad = Particles.radius.max()

	length, width = max(array(Particles.y * scale,'int')), max(array(Particles.x *scale, 'int'))

	Particles = Particles[Particles.z >= zmin - maxRad]
	Particles = Particles[Particles.z <= zmax + maxRad]

	Particles = Particles[fabs(Particles.z - zmax) <= Particles.radius]
	Particles = Particles[fabs(Particles.z - zmin) <= Particles.radius]
	
	N = Particles.natoms

	img = Image.new('RGB', (length, width), "black") # create a new black image
	pixels = img.load() # create the pixel map


	if N > 0:
		# map particle positions to pixels
		x = Particles.x
		x = array(x * scale, 'int')

		y = Particles.y
		y = array(y * scale, 'int')

		z = Particles.z

		# uncommnet the line below for alpha channel estimation
		# trans = array(z * 255 / z.max(), 'int')

		zmean = (zmax + zmin) * 0.5

		r = sqrt(Particles.radius**2.0 - (z - zmean)**2.0)
		r = array(r * scale, 'int')
		
		for n in range(N):

			i, j = x[n], y[n]
			radius = r[n]

			if (i + radius < length and i - radius > 0) and (j + radius < width and j - radius > 0):
				for ic in range(-radius, radius+1):
					for jc in range(-radius, radius+1):
						if (ic)**2 + (jc)**2 <= radius**2:
				        		pixels[i+ic,j+jc] = (255, 255, 255) #  add trans[n] for transparency (e.g. png files) and then set the colour accordingly
			else:
				pass

	if imgShow:
		img.show()

	if output:
		img.save(output)