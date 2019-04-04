from PyGran import simulation
from PyGran.modules.params import organic

xdim = 1e-3

organic['cohesionEnergyDensity'] = 1e5
organic['youngsModulus'] = 6e7

wall = organic.copy()
wall['cohesionEnergyDensity'] = 0.0

PID = wall.copy()
wall['coefficientFriction'] = 0

# Create a dictionary of physical parameters
params = {

	# Define the system
	'boundary': ('p','f','f'), # fixed BCs
	'box':  (-xdim, xdim, -xdim, xdim, 0, xdim*12),

	'model': simulation.models.HertzMindlin,

	# Define component(s)
	'species': ({'material': organic, 'radius': ('constant', 2e-4)},
		),

	# Timestep
	'dt': 1e-6,

	'nns_skin': 1e-4,

	# Apply gravitional force in the negative direction along the z-axis
	'gravity': (9.81, 0, 0, -1),

	# Number of simulation steps (non-PyGran variable)
	'run': 1e5,

	'traj': {'mfile': 'mesh*.vtk', 'freq':1000, 'pfile': 'traj.dump'},

	'output': 'output',

	# Import surface mesh
	'mesh': {
		'PID': {'file': 'mesh/wall-2D.stl', 'mtype': 'mesh/surface/stress/servo', 'material': PID, 'args': {'scale':9.8e-4, 'move': (0, 0, 6e-3), 'com': (0, 0, 0), 'axis': (0, 0, -1), 'vel_max': 10.0, 'ctrlPV': 'force', 'target_val': 5.0}},
		'wallX1': {'file': 'mesh/wall-2D-side.stl', 'mtype': 'mesh/surface/stress', 'material': wall, 'args': {'scale':1e-3, 'rotate': ('axis', 0, 1, 0, 'angle', 90), 'move': (0, 0, 2e-3)}}, 
		'wallX2': {'file': 'mesh/wall-2D-side.stl', 'mtype': 'mesh/surface/stress', 'material': wall, 'args': {'scale':1e-3, 'rotate': ('axis', 0, 1, 0, 'angle', -90), 'move': (0, 0, 2e-3) }},
		'wallY1': {'file': 'mesh/wall-2D-side.stl', 'mtype': 'mesh/surface/stress', 'material': wall, 'args': {'scale':1e-3, 'rotate': ('axis', 1, 0, 0, 'angle', 90, 'rotate', 'axis', 0, 1, 0, 'angle', 90), 'move': (0, 0, 2e-3)}},
		'wallY2': {'file': 'mesh/wall-2D-side.stl', 'mtype': 'mesh/surface/stress', 'material': wall, 'args': {'scale':1e-3, 'rotate': ('axis', 1, 0, 0, 'angle', -90, 'rotate', 'axis', 0, 1, 0, 'angle', 90), 'move': (0, 0, 2e-3)}}
		}
}

def run(**params):
	# Create an instance of the DEM class
	sim = simulation.DEM(**params)

	# Setup a primitive wall along the xoy plane at z=0
	sim.setupWall(species=1, wtype='primitive', plane = 'zplane', peq = 0.0)

	# Insert particles by volume fraction
	sim.insert(species=1, value=1.0, mech='volumefraction_region', region=('block', '-{} {} -{} {} 0 {}'.format(xdim, xdim, xdim, xdim, xdim*8)))

	# Run the simulation
	sim.run(params['run'], params['dt'])

if __name__ == '__main__':
	run(**params)
