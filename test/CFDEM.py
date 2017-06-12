from PyGran import Analyzer
from matplotlib.mlab import find
from matplotlib import pylab as plt 
from scipy.linalg import norm
from numpy import array

nProc = 4

# Create empty lists on all processors
pressDrop = [[] for i in range(nProc)]
velInlet = [[] for i in range(nProc)]
voidFraction = [[] for i in range(nProc)]

# CFDEM working dir
wdir = '/home/abimanso/Disk-1/CFDEM/cfd-dem/tutorials/FluidBed/'

# Create a particle trajectory file and go to last frame
pTraj = Analyzer.System(fname=wdir+'/DEM/post/CFDEM*.dump')
pTraj.goto(-1)

for i in range(nProc):

	# Create a mesh trajectory file
	fluid = wdir + 'CFD/processor{}/VTK/processor{}_*.vtk'.format(i,i)
	mTraj = Analyzer.System(mfname=fluid)

	# Compute center of each mesh cell
	cellCenter = mTraj.Mesh.CellsPos.mean(axis=1)

	inlet = find(cellCenter[:,2] == cellCenter[:,2].min())
	outlet = find(cellCenter[:,2] == cellCenter[:,2].max())

	# Loop over mesh trajectory and compute the pressure drop
	for ts in mTraj:

		# Find the boundaries of the bed
		bed = find(cellCenter[:,2] <= pTraj.Particles.z.max())

		# Compute the variables of interest
		pressDrop[i].append(mTraj.Mesh.p[inlet].mean() - mTraj.Mesh.p[outlet].mean())
		velInlet[i].append(norm(mTraj.Mesh.U[inlet,:].mean(axis=0)))
		voidFraction[i].append(mTraj.Mesh.voidfraction[bed].mean(axis=0))

# Average computed variables across all processors
pressDrop = array(pressDrop).mean(axis=0)
velInlet = array(velInlet).mean(axis=0)
voidFraction = array(voidFraction).mean(axis=0)

# Plot pressure drop
plt.plot(pressDrop)
plt.xlabel('Frame')
plt.ylabel('Pressure drop (Pa)')
plt.grid()
plt.show()

# Plot inlet velocity
plt.plot(velInlet)
plt.xlabel('Frame')
plt.ylabel('Inlet velocity (m/s)')
plt.grid()
plt.show()

# Plot void fraction
plt.plot(velInlet, voidFraction)
plt.xlabel('Inlet velocity (m/s)')
plt.ylabel('Mean void fraction')
plt.grid()
plt.show()
