from PyGran import Analyzer
from matplotlib.mlab import find
from scipy.linalg import norm

Fluid = Analyzer.System(mfname='/home/levnon/Desktop/VTK/processor0_*.vtk')

# Compute center of each mesh cell
cellCenter = Fluid.Mesh.CellsPos.mean(axis=1)

inlet = find(cellCenter[:,2] == cellCenter[:,2].min())
outlet = find(cellCenter[:,2] == cellCenter[:,2].max())

# Loop over mesh trajectory and compute the pressure drop
pressDrop = []
velInlet = []

for ts in Fluid:
	pressDrop.append(Fluid.Mesh.p[inlet].mean() - Fluid.Mesh.p[outlet].mean())
	velInlet.append(norm(Fluid.Mesh.U[inlet,:].mean(axis=0)))
