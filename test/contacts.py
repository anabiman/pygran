from PyGran import Analyzer
from PyGran.Materials import glass

System = Analyzer.System('traj.dump')

# Go to last frame
System.goto(-1)

Neigh = Analyzer.Neighbors(System.Particles, material=glass)

