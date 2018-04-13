from PyGran import Analyzer

# Create a granular object from a LIGGGHTS dump file
Sys = Analyzer.System(Particles='traj.dump', units='micro')

# Compute the radial distribution function
g, r, _ = Sys.Particles.rdf()

# Construct a class for nearest neighbor searching
Neigh = Analyzer.equilibrium.Neighbors(Sys.Particles)

# Extract coordination number per particle
coon = Neigh.coon

