from PyGran import Analyzer
from matplotlib.mlab import find

Fluid = Analyzer.System(mfname='CFD/processor0/VTK/processor0_*.vtk')

# Go to last frame to find inlet / outlet
Fluid.goto(-1)
inlet = find(Fluid.Mesh.p == Fluid.Mesh.p.min()).min()
outlet = find(Fluid.Mesh.p == Fluid.Mesh.p.max()).max()

# Loop over mesh trajectory and compute the pressure drop
pressure = []

for ts in Fluid:
	pressure.append(Fluid.p[outlet] - Fluid.p[inlet])