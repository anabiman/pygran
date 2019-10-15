from PyGran import analysis
import sys

# Create a granular object from a LIGGGHTS dump file
Sys = analysis.System(Particles=sys.argv[1], units='micro')

# Go to last frame
Sys.goto(-1)

# Compute the radial distribution function
g, r, _ = Sys.Particles.rdf()

# Construct a class for nearest neighbor searching
Neigh = analysis.equilibrium.Neighbors(Sys.Particles)

# Extract coordination number per particle
coon = Neigh.coon
