
from PyGran import simulation as sim
from PyGran.params import stearicAcid, glass
import numpy as np

# Create a dictionary of physical parameters
params = {

	# Define the system
	'boundary': ('p','p','f'),
	'box':  (-0.00025, 0.00025, -0.00025, 0.00025, 0, 0.002),

	# Define component(s)
	'species': ({'material': stearicAcid, 'radius': ('constant', 25.00e-6)},
		{'material': glass, 'radius': ('constant', 25.00e-6)}
		),

	# Setup I/O
	'traj': {'freq': 100000, 'pfile': 'traj.dump', 'mfile':'mesh-*.vtk'},
	'output': 'Binary',

	# Write a restart file (restart.binary.*) every 5000 steps to 'restart' dir 
	'restart': (5000, 'restart', 'restart.binary', False, None),

	# Define (optional) timestep to be used later
	'dt': 2e-9,

	# Apply gravitional force in the negative direction along the z-axis
	'gravity': (9.81, 0, 0, -1),

	# Stage runs
	'stages': {'insertion': 1e8}
}

# Create an instance of the DEM class
sim = sim.DEM(**params)

# Insert all particles throughout the sim box
insert = sim.insert(species=1, value=100)

# Run the simulation
sim.run(params['stages']['insertion'], params['dt'])

# Remove insertion
sim.remove(insert)

sim.run(params['stages']['insertion'], params['dt'])
