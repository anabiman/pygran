'''
  This is the 
   __________         ________                     
  ██████╗ ██╗   ██╗ ██████╗ ██████╗  █████╗ ███╗   ██╗
  ██╔══██╗╚██╗ ██╔╝██╔════╝ ██╔══██╗██╔══██╗████╗  ██║
  ██████╔╝ ╚████╔╝ ██║  ███╗██████╔╝███████║██╔██╗ ██║
  ██╔═══╝   ╚██╔╝  ██║   ██║██╔══██╗██╔══██║██║╚██╗██║
  ██║        ██║   ╚██████╔╝██║  ██║██║  ██║██║ ╚████║
  ╚═╝        ╚═╝    ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═══╝
                                                      
  DEM simulation and analysis toolkit
  http://www.pygran.org, support@pygran.org

  Core developer and main author:
  Andrew Abi-Mansour, andrew.abi.mansour@pygran.org

  PyGran is open-source, distributed under the terms of the GNU Public
  License, version 2 or later. It is distributed in the hope that it will
  be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
  received a copy of the GNU General Public License along with PyGran.
  If not, see http://www.gnu.org/licenses . See also top-level README
  and LICENSE files.
 '''

from PyGran import Analyzer
from PyGran.Materials import stearicAcid
from numpy import pi, sqrt, array
from scipy import optimize
import matplotlib.pylab as plt

stearicAcid['cohesionEnergyDensity'] = 0.033
stearicAcid['yieldPress'] = 2.2e6

def checkYieldNum(reff, **material):
	""" Solves numerically the cubic equation x^3 - b*x - c = 0 for yielding contact_radius = sqrt(x)
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

	def eq(x, *args):

		py, gamma, YoungEff = args

		b = py * pi * reff / (2 * YoungEff)
		c = reff * sqrt(gamma * pi / (2 * YoungEff))

		return x**3 - b * x - c

	x0 = sqrt(py * pi * reff / (2. * YoungEff))
	x = optimize.fsolve(func=eq, x0=x0, args=(py, gamma, YoungEff), xtol=1e-16)
	ay = x * x

	return ay*ay/reff - sqrt(2. * pi * gamma * ay / YoungEff)

def checkYield(reff, **material):
	""" Solves symbolically the cubic equation x^3 - b*x - c = 0 for yielding contact_radius = sqrt(x)
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
	b = py * pi * reff / (2. * YoungEff)
	c = reff * sqrt(gamma * pi / (2. * YoungEff))

	# Solve the algebraic equation symbolically
	frac = 0.333333333333333333333333333333333333333333333333333	
	common = (9*c + sqrt(3*(27*c**2 - 4*b**3)))**frac

	x = b * (2.0/3.0)**frac /  common + common / (18**frac)

	# Compute contact yield radius
	ay = x*x

	# Return yielding contact radius
	return ay*ay/reff - sqrt(2. * pi * gamma * ay / YoungEff)


System = Analyzer.System('traj-tapping.dump')
data = []

for ts in System:
	Neigh = Analyzer.Neighbors(System.Particles)
	overlaps, indices = Neigh.overlaps[:,0], Neigh.overlaps[:,1:]

	# Extract radii of all particles
	radii = System.Particles.radius

	ny = 0

	for contact, index in enumerate(indices):

		# Get the two particle (in contact) indices
		i,j = index

		# Compute reff
		reff =  (radii[i] * radii[j]) / (radii[i] + radii[j])

		# Get overlap
		delta = overlaps[contact]

		# Compute yield contact radius symbolically and numerically
		# ay = checkYield(reff, **stearicAcid)
		deltay = checkYieldNum(reff, **stearicAcid)

		if delta >= deltay: 
			ny += 1
		
	data.append([ts, ny])	

data = array(data)
#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#plt.plot(data, '.'); plt.legend(['Hertzian', 'Numerical', 'Symbolic'])
#plt.grid(lineStyle=':')
#plt.ylabel('Yield contact radius (m)')
#plt.xlabel('Contact #')
#plt.show()
