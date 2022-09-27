import matplotlib.pylab as plt
from pygran import analysis

# Create a granular object from a LIGGGHTS dump file
Sys = analysis.System(Particles="traj*.dump", units="micro")

# Go to last frame
Sys.goto(-1)

# Compute the radial distribution function
g, r, _ = Sys.Particles.rdf()

# Plot rdf vs radial distance
plt.plot(r, g)
plt.show()

# Construct a class for nearest neighbor searching
Neigh = analysis.equilibrium.Neighbors(Sys.Particles)

# Extract coordination number per particle
coon = Neigh.coon
