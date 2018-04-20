from PyGran import Analyzer
import matplotlib.pylab as plt

# Create a granular object from a LIGGGHTS dump file
Sys = Analyzer.System(Particles='/home/levnon/Desktop/flow/output/traj/traj*.dump', units='micro')

# Go to last frame
Sys.goto(-1)

# Compute the radial distribution function
g, r, _ = Sys.Particles.rdf()

# Plot rdf vs radial distance
plt.plot(r, g)
plt.show()

# Construct a class for nearest neighbor searching
Neigh = Analyzer.equilibrium.Neighbors(Sys.Particles)

# Extract coordination number per particle
coon = Neigh.coon
