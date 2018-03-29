import scipy.spatial
from PyGran import Analyzer
from numpy import arange, array
import os

wdir = '/home/abimanso/Desktop/Papers/Manuscript-1/Void-analysis/Stearic-acid/fractal/DEM-sim/traj/'

Granl = Analyzer.SystemFactory([wdir + 'traj.dump']) #([wdir + 'Sim4-nc/traj/traj-nc-40Hz.dump'])
system = 0

for Gran in Granl.System:

	Gran.goto(-1)
	system += 1
	os.system('mkdir system{}'.format(system))

	Particles = Gran.Particles

	print Particles.natoms

	#Particles = Gran.Particles[Gran.Particles.z >= 0.5e-3]
	Particles = Particles[Particles.z <= 1.5e-3]

	Particles = Particles[Particles.x >= -0.5e-3]
	Particles = Particles[Particles.x <= 0.5e-3]


	Particles = Particles[Particles.y >= -0.5e-3]
	Particles = Particles[Particles.y <= 0.5e-3]

	print Particles.natoms

	resol = 1.24e-6 # microns/pixel
	count = 0

	Particles.x -= Particles.x.min()
	Particles.y -= Particles.y.min()
	Particles.z -= Particles.z.min()

	for z in arange(.0, Particles.z.max() + resol, resol):
		print count
		zmin, zmax = z, z + resol
		count += 1
		Analyzer.xrct.genImg(Particles, zmin, zmax, resol, 'system{}/poured{}.bmp'.format(system, count))
