from PyGran.Analyzer import System
from scipy.linalg import norm

# Create a mesh trajectory file for the inlet & outlet files
inlet, outlet = 'CFD/inlet_*.vtk', 'CFD/outlet_*.vtk'
Traj = System(Particles='DEM/*.dump', Mesh=[inlet, outlet], vtk_type='poly')

# Compute mesh surface areas
iMesh, oMesh = Traj.Mesh
iArea, oArea = iMesh.CellArea.sum(), oMesh.CellArea.sum()

# Loop over inlet trajectory and compute the inlet pressure & vel
for i, timestep in enumerate(Traj):

	# Compute the weghted-average pressure inlet + outlet
	iPress[i].append((iMesh.p * iMesh.CellArea).sum() / iArea)
	oPress[i].append((oMesh.p * oMesh.CellArea).sum() / oArea)

	# Compute the weighted-average inlet velocity
	iVel[i].append(norm((iMesh.U.T * iMesh.CellArea).sum(axis=1) / iArea))

	# Compute mean particle position along the z-azis
	zMean[i].append(Traj.Particles.z.mean())
