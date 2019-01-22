import scipy.spatial
from PyGran import Analyzer
from numpy import arange, array
import os

Gran = Analyzer.System(Particles='traj.dump')

# Go to last frame
Gran.goto(-1)

# Create a new class containing particles between  z=0 and z=1e-3
Particles = Gran.Particles
Particles = Particles[Particles.z <= 1e-3 & Particles.z >= 0]

# Set image resolution and size
resol = 1.24e-6 # microns/pixel
size = (512, 512)

for i, z in enumerate(arange(0, Particles.z.max() + resol, resol)):
	zmin, zmax = z, z + resol
	output =  'output/poured{}.bmp'.format(i)
	Analyzer.imaging.slice(Particles, zmin, zmax, 'z', size, resol, output)
