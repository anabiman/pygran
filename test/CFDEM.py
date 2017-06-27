from PyGran import Analyzer
from matplotlib.mlab import find
from matplotlib import pylab as plt 
from scipy.linalg import norm
from numpy import array

nProc = 4

# Create empty lists on all processors
pressInlet = [[] for i in range(nProc)]
pressOutlet = [[] for i in range(nProc)]
velInlet = [[] for i in range(nProc)]
voidFraction = [[] for i in range(nProc)]

# CFDEM working dir
wdir = '/home/abimanso/Disk-1/CFDEM/cfd-dem/tutorials/FluidBed/'

# Create a particle trajectory file and go to last frame
pTraj = Analyzer.System(fname=wdir+'/DEM/post/CFDEM*.dump')
pTraj.goto(-1)

for i in range(nProc):

	# Create a mesh trajectory file for the bed
	domain = wdir + 'CFD/processor{}/VTK/processor{}_*.vtk'.format(i,i)
	bTraj = Analyzer.System(mfname=domain)

	# Create a mesh trajectory file for the inlet and outlet
	inlet = wdir + 'CFD/processor{}/VTK/inlet/inlet_*.vtk'.format(i)
	iTraj = Analyzer.System(mfname=inlet, vtk_type='poly')

	outlet = wdir + 'CFD/processor{}/VTK/outlet/outlet_*.vtk'.format(i)
	oTraj =	Analyzer.System(mfname=outlet, vtk_type='poly')

	# Compute center of each mesh cell
	cellCenter = bTraj.Mesh.CellsPos.mean(axis=1)

	# Loop over bed trajectory and compute mean void fraction
	for ts in bTraj:

		bed = find(cellCenter[:,2] <= pTraj.Particles.z.max())

		# compute the weighted-average mean void fraction
		voidFrac = (bTraj.Mesh.voidfraction[bed] * bTraj.Mesh.CellVol[bed]).sum() / bTraj.Mesh.CellVol[bed].sum()
		voidFraction[i].append(voidFrac)

	# Loop over inlet trajectory and compute the inlet pressure + vel
	for ts in iTraj:
		# Compute the weighted-average pressure drop
		iPress = (iTraj.Mesh.p * iTraj.Mesh.CellArea).sum() / iTraj.Mesh.CellArea.sum()
		pressInlet[i].append(iPress)

		# compure the weighted-average inlet velocity
		iVel = (iTraj.Mesh.U.T * iTraj.Mesh.CellArea).sum(axis=1) / iTraj.Mesh.CellArea.sum()
		velInlet[i].append(norm(iVel))

	# Loop over outlet trajectory and compute the outlet vel
	for ts in oTraj:
		oPress = (oTraj.Mesh.p * oTraj.Mesh.CellArea).sum() / oTraj.Mesh.CellArea.sum()
		pressOutlet[i].append(oPress)

# Average computed variables across all processors
pressInlet = array(pressInlet).mean(axis=0)
pressOutlet = array(pressOutlet).mean(axis=0)
velInlet = array(velInlet).mean(axis=0)
voidFraction = array(voidFraction).mean(axis=0)

# Plot pressure drop
plt.plot(pressInlet - pressOutlet)
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
