
from PyGran import simulation as sim
from PyGran.modules.params import stearicAcid, glass
import numpy as np

def run(**params):

	# Create an instance of the DEM class
	sim = sim.DEM(**params)

	# Insert all particles throughout the sim box
	for species in [1,2]:
		insert = sim.insert(species, value=100)
		sim.run(1, params['dt'])

		# Remove insertion
		sim.remove(insert)

	# Run the simulation
	sim.run(params['stages']['insertion'], params['dt'])


	# Setup params for vibrating mesh
	freq = 40 * 2 * np.pi
	nTaps, period = 100, 1.0 / freq
	nSteps = period / (pDict['dt'])

	for i in range(nTaps):
		moveMesh = sim.moveMesh(name='wall', viblin=('axis', 0, 0, 1), order=1, amplitude=12.5e-6, phase=0, period=period)
		sim.run(nSteps, params['dt'])

		sim.remove(moveMesh)
		# Relax the system if requested by the user

	if params['relax']:
		sim.run(nSteps, params['dt'])

if __name__ == '__main__':

	# Create a dictionary of physical parameters
	params = {

		# Define the system
		'boundary': ('p','p','f'),
		'box':  (-0.00025, 0.00025, -0.00025, 0.00025, 0, 0.0005),

		# Define component(s)
		'species': ({'material': stearicAcid, 'radius': ('constant', 25.00e-6)},
			{'material': glass, 'radius': ('constant', 25.00e-6)}
			),

		# Setup I/O
		'traj': {'freq': 100000, 'pfile': 'traj.dump'},
		'output': 'Binary',

		# Write a restart file (restart.binary.*) every 5000 steps to 'restart' dir 
		'restart': (5000, 'restart', 'restart.binary', False, None),

		# Define (optional) timestep to be used later
		'dt': 2e-9,

		# Apply gravitional force in the negative direction along the z-axis
		'gravity': (9.81, 0, 0, -1),

		# Stage runs
		'stages': {'insertion': 1e8, 'run': 1e7, 'relax': 0e5},

		# Meshes
		#'mesh': {
		#	'wall': {'file': 'square.stl', 'mtype': 'mesh/surface/stress', 'material': stearicAcid, \
		#		 'args': {'scale': 2.5e-4, 'move': (0, 0, -8e-4)}
		#	},
		#},
	}

	run(**params)

	

