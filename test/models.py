import PyDEM.Simulator as Si
from numpy import sqrt, array, arange, pi

radius = 1e-5
mass = 2500.0 * 4.0/3.0 * pi * radius**3.0

cohDensity =  arange(0.1, 1.1, 0.1) * 10.0
JKRRad = []

for den in cohDensity:
	Hertz = Si.HertzMindlin(**{'poissonsRatio': 0.4, 'youngsModulus': 7.1e7, 'mass':mass, 'cohesionEnergyDensity': den, 'radius': radius, 'characteristicVelocity':0.1})
	_, disp = Hertz.displacement()

	disp = disp.mean()
	print disp

	yEff = Hertz.materials['youngsModulus'] * 0.5 /  (1.0 - Hertz.materials['poissonsRatio'])

	HertzRad = sqrt(disp * radius)
	JKRRad.append([Hertz.contactRadius(disp)[0], HertzRad + sqrt(pi * den / (2.0 * yEff * HertzRad)) * radius / 2.0])

JKRRad = array(JKRRad)
