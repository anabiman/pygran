from PyGran import Analyzer
from PyGran.Materials import glass

System = Analyzer.System('traj.dump')
System.goto(-1)

Neigh = Analyzer.Neighbors(System.Particles)
overlaps, i, j = Neigh.overlaps


