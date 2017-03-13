from PyGran import Analyzer
G = Analyzer.Granular('traj.dump')
G.goto(2)
#print G.data['x'].shape, G.Particles.x.shape

N = Analyzer.Neighbors(G.Particles)
