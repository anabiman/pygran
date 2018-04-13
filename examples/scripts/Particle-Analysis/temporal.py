from PyGran import Analyzer

# Material propreties and simulation parameters for glass
tDensity, timestep = 2500.0, 1e-6

# Create a PyGran System from a series of ESyS (trajectory) files
Gran = Analyzer.System(Particles='/home/levnon/Desktop/compaction/out-SpringDashpot/traj/traj.dump')

# Each Granular object contains a Particle object that can be sliced
Parts = Gran.Particles

# Compute the mass flow rate using the 1st 10 frames
Temp = Analyzer.dynamics.Temporal(Gran)
rate = Temp.computeFlow(density=tDensity, dt=timestep)

# Compute the bulk density as a time series by looping over the Gran trajectory
density = []

for ts in Gran:
	density.append(Parts.density(tDensity))
