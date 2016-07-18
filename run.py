# !/usr/bin/python
# -*- coding: utf8 -*-

from PyDEM import Simulator as Simulator, Analyzer, Visualizer


if __name__ == '__main__':

	# Create a dictionary of physical parameters
	physParams = {

			  # Define the system
			  'boundary': ('f','f','f'), # fixed BCs
			  'model': ('gran', 'model', 'hooke', 'tangential', 'history', 'rolling_friction', \
                        'cdt', 'tangential_damping', 'on', 'limitForce', 'on'), # the order matters here
			  'box':  (-0.016, 0.016, -0.016, 0.016, -0.01, 0.061), # simulation box size

			  # Define component(s)
			  'SS': ({'id':1, 'natoms': 1e5, 'density': 2500.0, 'insert': True, 'rate': 1e6, 'freq': 1e3}, ), # number of particles inserted = rate * dt * freq every freq steps
			  'nSS': 2,
			  'vel': ((0,0,0), ),
			  'radius': (('gaussian number', 4e-4, 4e-5),),

			  # Material properties
			  'materials': {
			  				'yMod': ('youngsModulus', 'peratomtype', '7.1e7', '7.1e7'),
			  				'pRatio': ('poissonsRatio', 'peratomtype', '0.22', '0.22'),
			  				'cFric': ('coefficientFriction', 'peratomtypepair', '2', '0.5', '0.5', '0.5', '0.5'),
			  				'cRollFric': ('coefficientRollingFriction', 'peratomtypepair', '2', '5e-4', '5e-4', '5e-4', '5e-4'),
			  				'cVel': ('characteristicVelocity', 'scalar', '0.1', '0.1'),
			  				'cRest': (('coefficientRestitution', 'peratomtypepair', '2', '1.0', '1.0', '1.0', '1.0'))
			  				},

			  # Apply gravitional force in the negative direction along the z-axis
			  'gravity': (9.81, 0, 0, -1), 
			}
			
	# Create a dictionary of dynamical parameters for a linear spring model
	coreParams = {

				'model': DEMAn.models.LinearSpring(**physParams),
                'engine': DEMSi.engines.liggghts,

				# Stage runs
				'insertion': 1e6,
				'flow': 0e6
				}
    
    # Estimate allowed sim timestep
	dynamicParams['dt'] = (dynamicParams['alpha'] * dynamicParams['lsModel'].contactTime()).min()

	for ss in physParams['SS']:
		print 'Attempting to insert {} particles every {} steps for component {}'.format(ss['rate'] * ss['freq'] * dynamicParams['dt'], ss['freq'], ss['id'])

	compParams = {

				# Use LIGGGHTS as the DEM computational engine
				'modName': 'liggghts',

				# I/O parameters
				'restart': (5000, 'restart', 'restart.binary', False),
				'traj': {'sel': 'all', 'freq': 1000, 'dir': 'traj', 'style': 'custom', 'file': 'traj.dump', \
				       'args': ('id', 'x', 'y', 'z', 'radius', 'omegax', 'omegay', 'omegaz', \
				                'vx', 'vy', 'vz', 'fx', 'fy', 'fz')},
				'dump_modify': ('append', 'yes'),
				'nSim': 1,
				'output': 'out-radius-400-repos',
				'print': (10**4, 'time', 'atoms', 'fmax', 'ke', 'cpu', 'cu', 'density'),

				# Meshes
				'surfMesh': {
					'hopper': {'file': 'hopper-2cm-6cm.stl', 'scale': 5e-4},
				      },

				# Nearest neighbor searching params
				'nns_skin': 1e-3,
				'nns_type': 'bin',

				# grab shared libraries from here
				'path': '/usr/lib64/',
			  }

	# Create an instance of the DEM class
	sim = DEMSi.DEM(**dict(physParams, **compParams))

	# Define the domain, create atoms, and initialize masses, velocities, etc.
	sim.initialize()

	# Setup material properties
	for item in physParams['materials'].keys():
		# Overloaded function 'createProperty' will partition coeffRest based on MPI's coloring split scheme
		sim.createProperty(item, *physParams['materials'][item])

	# Import and setup all meshes as rigid wall
	for mesh in compParams['surfMesh'].keys():
		sim.importMesh(name=mesh, **compParams['surfMesh'][mesh])
		sim.setupWall(name=mesh + 'Wall', wtype='mesh', meshName=mesh)

	# Setup a stopper wall along the xoz plane (y = 0.0)
	sim.setupWall(name='stopper', wtype='primitive', plane = 'zplane', peq = 0.0)

	# Print output specified in 'print' every 'freq' steps
	sim.printSetup()

	# Create an NVE (micro canonical) integrator
	sim.setupIntegrate(name='intMicro')

	# Write 'all' coordinates to 'traj.xyz' file every 'freq' steps
	sim.dumpSetup()

	if not compParams['restart'][3]:
		# Insert particles if not restarting/resuming sim
		cylinder = sim.insertParticles('void', *('cylinder', 'z', 0, 0, 0.001, 0.02, 0.045))
		sim.integrate(dynamicParams['insertion'], dynamicParams['dt'])
		sim.remove(cylinder)

	# Remove stopper
	sim.remove(name='stopper')

	sim.integrate(dynamicParams['flow'], dynamicParams['dt'])
