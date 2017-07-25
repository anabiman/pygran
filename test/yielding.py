from PyGran import Analyzer
from PyGran.Materials import stearicAcid
from numpy import pi, sqrt

stearicAcid['yieldPress'] = 2e6

def checkYield(reff, **material):
	""" Solves the cubic equation x^3 - b*x - c = 0 for yielding contact_radius = sqrt(x)
	based on Thornton's elasto-plastic cohesive model:

	b = py * pi * reff / (2 * YoungEff)
	c = reff * sqrt(gamma * pi / (2 * YoungEff))

	@reff: effective radius
	@py: yielding pressure
	@YoungEff: Young's effective modulus
	@gamma: cohesion energy density
	"""
	# Extract material params from supplied database
	py = material['yieldPress']
	poiss = material['poissonsRatio']
	gamma = material['cohesionEnergyDensity']
	Young = material['youngsModulus']
	YoungEff = Young * 0.5 / (1.0  - poiss)

	# Compute the 'b' and 'c' coefficients
	b = py * pi * reff / (2 * YoungEff)
	c = reff * sqrt(gamma * pi / (2 * YoungEff))

	# Solve the algebraic equation symbolically
	frac = 0.333333333333333333333333333333333333333333333333333
	common = (9*c + sqrt(3*(27*c**2 - 4*b**3)))**frac

	print 3*(27*c**2 - 4*b**3), 9*c

	#print 27*c**2 - 4*b**3
	x = b * (2.0/3.0)**frac /  common + common / (18**frac)

	# Compute contact yield radius
	ay = x**2

	# Return yield overlap
	return ay**2 / reff - sqrt(2 * pi * gamma * ay / YoungEff)


System = Analyzer.System('../../compaction/out-SpringDashpot/traj/traj.dump')
System.goto(-1)

Neigh = Analyzer.Neighbors(System.Particles)
overlaps, indices = Neigh.overlaps[:,0], Neigh.overlaps[:,1:]

# Extract radii of all particles
radii = System.Particles.radius

for contact, index in enumerate(indices):

	# Get the two particle (in contact) indices
	i,j = index

	# Compute reff
	reff = 1.0 / (1.0 / radii[i] + 1.0 / radii[j])

	# Compute yielding overlap
	delta_y = checkYield(reff, **stearicAcid)

	print overlaps[contact], delta_y
