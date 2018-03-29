from PyGran import Simulator, Analyzer, Visualizer
from PyGran.Materials import glass

glass['youngsModulus'] = 1e7

# Create a dictionary of physical parameters
pDict = {

		# Define the system
		'boundary': ('p','p','p'), # fixed BCs
		'box':  (-0.001, 0.001, -0.001, 0.001, 0, 0.004), # simulation box size

		# Define component(s)
		'SS': ({'material': glass, 'radius': ('constant', 2e-4)}, 
		      ),

		# Timestep
		'dt': 1e-6,

		# Apply gravitional force in the negative direction along the z-axis
		'gravity': (9.81, 0, 0, -1),

		'output': 'MS-DEM',

		# Number of simulation steps (non-PyGran variable)
		'nsteps': 2.5e4,

		# Import surface mesh
		 'mesh': {
			'wallZ': {'file': 'mesh/square.stl', 'mtype': 'mesh/surface/stress', 'material': glass, 'args': ('scale 1e-3',)}
		      },
	  }

if __name__ == '__main__':

	# Create an instance of the DEM class
    sim = Simulator.DEM(**pDict)

	# Setup a primitive wall along the xoy plane at z=0
	sim.setupWall(species=1, wtype='primitive', plane = 'zplane', peq = 0.0)

	# Insert the particles
	insert = sim.insert(species=1, value=100, region=('block', -1e-3,1e-3, -1e-3, 1e-3, 0, 3e-3))
	sim.run(1e4, pDict['dt'])
	sim.remove(insert)

    # Move wall at constant speed
    moveZ = sim.moveMesh('wallZ', *('linear', '0 0 -0.01'))
    sim.run(pDict['nsteps'], pDict['dt'])
    sim.remove(moveZ)

    # Relax the system
    moveZ = sim.moveMesh('wallZ', *('linear', '0 0 0.01'))
    sim.run(pDict['nsteps'], pDict['dt'])
