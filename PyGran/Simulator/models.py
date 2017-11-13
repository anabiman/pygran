'''
Created on July 1, 2016
@author: Andrew Abi-Mansour
'''

# !/usr/bin/python
# -*- coding: utf8 -*-
# -------------------------------------------------------------------------
#
#   Python module for analyzing contact models for DEM simulations
#
# --------------------------------------------------------------------------
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

# -------------------------------------------------------------------------

# TODO: Support 2-particle analysis by replacing mass radius, etc. with reduces mass, radius, etc.
# i.e. 1/m_ij = 1/m_i + 1/m_j

import numpy as np
from scipy.integrate import ode
from scipy.optimize import fsolve
from PyGran import Materials
import math

class Model(object):
	def __init__(self, **params):

		self.params = params
		self.params['nSS'] = 0

		if 'SS' in self.params:
			self.params['nSS'] += len(self.params['SS'])

		idc = 1

		if 'SS' in self.params:
			for ss in self.params['SS']:

				if 'id' not in ss:
					ss['id'] = idc
					idc += 1

		# Treat any mesh as an additional component
		if 'mesh' in self.params:
			self.params['nSS'] += len(self.params['mesh'])
			for mesh in self.params['mesh']:
				self.params['SS'] += ({'material':self.params['mesh'][mesh]['material']},)

				if 'id' not in self.params['mesh'][mesh]:
					self.params['mesh'][mesh]['id'] = idc
					idc += 1

				if 'args' not in self.params['mesh'][mesh]:
					self.params['mesh'][mesh]['args'] = ()

		if 'style' not in self.params:
			self.params['style'] = 'granular'

		if 'units' not in self.params:
			self.params['units'] = 'si'

		if 'dim' not in self.params:
			self.params['dim'] = 3

		if 'vel' not in self.params:
			self.params['vel'] = ()

			if 'SS' in self.params:
				# the user did not specify the initial velocities, so we assume they're zero
				for comp in range(len(self.params['SS'])):
					self.params['vel'] += ((0,0,0),)

		if 'nns_skin' not in self.params:
			self.params['nns_skin'] = 1e-3

		if 'nns_type' not in self.params:
			self.params['nns_type'] = 'bin'

		if 'nns_freq' not in self.params:
                        self.params['nns_freq'] = 100

		if 'restart' not in self.params:
			self.params['restart'] = (5000, 'restart', 'restart.binary', False, None)

		if 'dump_modify' not in self.params:
			self.params['dump_modify'] = ('append', 'yes')

		if 'nSim' not in self.params:
			self.params['nSim'] =  1

		if 'read_data' not in self.params:
			self.params['read_data'] = False

		# Default traj I/O args
		traj = {'sel': 'all', 'freq': 1000, 'dir': 'traj', 'style': 'custom', 'pfile': 'traj.dump', \
                   'args': ('id', 'x', 'y', 'z', 'radius', \
                   'vx', 'vy', 'vz', 'fx', 'fy', 'fz')}

		if 'traj' in self.params:
			for key in self.params['traj']:
				traj[key] = self.params['traj'][key]

		self.params['traj'] = traj

        # Compute mean material properties
		self.materials = {}

		# For analysis
		if 'material' in params:
			self.materials = params['material']

		# Expand material properties based on number of components
		if 'SS' in self.params:
			for ss in self.params['SS']:
				ss['material'] = Materials.LIGGGHTS(**ss['material'])

			# Use 1st component to find all material params ~ hackish!!! 
			ss = self.params['SS'][0]

			if 'material' in ss:

				for item in ss['material']:
					if type(ss['material'][item]) is not float:
						# register each material proprety then populate per number of components
						if ss['material'][item][1] == 'peratomtype':
							self.materials[item] = ss['material'][item][:2]
						elif ss['material'][item][1] == 'peratomtypepair':
							self.materials[item] = ss['material'][item][:2] + ('{}'.format(self.params['nSS']),)
						elif ss['material'][item][1] == 'scalar':
							self.materials[item] = ss['material'][item][:2]
				
			for ss in self.params['SS']:
				for item in ss['material']:
					if type(ss['material'][item]) is float:
						
						# This is for running DEM sim
						ss[item] = ss['material'][item]


			for item in self.materials:
				if type(ss['material'][item]) is not float:

					for ss in self.params['SS']:

						if ss['material'][item][1] == 'peratomtype':
							self.materials[item] =  self.materials[item] + (('{}').format(ss['material'][item][2]),)

						elif ss['material'][item][1] == 'peratomtypepair':
							# assume the geometric mean suffices for estimating binary properties
							for nss in range(self.params['nSS']):

								prop = np.sqrt(float(ss['material'][item][2]) * float(self.params['SS'][nss]['material'][item][2]))
								self.materials[item] =  self.materials[item] + (('{}').format(prop),)

						# we should set this based on species type
						elif ss['material'][item][1] == 'scalar':
							self.materials[item] = self.materials[item] + (('{}').format(ss['material'][item][2]),)

						else:
							print('Error: Material database flawed.')
							sys.exit()

			self.params['materials'] = self.materials

		if 'dt' not in self.params:
			# Estimate the allowed sim timestep
			try:
				self.params['dt'] = (0.25 * self.contactTime()).min()
			except:
				self.params['dt'] = 1e-6

				if 'model' in self.params:
					print('Model {} does not yet support estimation of contact period. Using a default value of {}'.format(self.params['model'], self.params['dt']))

		self.setupProps()

		if hasattr(self,'density'):
			self.mass = 4.0 * self.density * 4.0/3.0 * np.pi * self.radius**3.0

	def setupProps(self):
		""" Creates class attributes for all material properties """
		for prop in self.materials:
			setattr(self, prop, self.materials[prop])

	def contactTime(self):
		""" Computes the characteristic collision time assuming for a spring dashpot model """

		rest = self.coefficientRestitution
		poiss = self.poissonsRatio
		yMod = self.youngsModulus
		radius = self.radius
		mass = self.mass

		# Create SpringDashpot instance to etimate contact time
		SD = SpringDashpot(material=self.materials)

		kn = SD.springStiff(radius)

		return np.sqrt(mass * (np.pi**2.0 + np.log(rest)) / kn)

	def displacement(self):
		""" Generator that computes (iteratively) the contact overlap as a function of time """

		y0 = np.array([0, self.characteristicVelocity])
		t0 = .0

		inte = ode(self.numericalForce)
		inte.set_f_params(*())
		inte.set_integrator('dopri5')
		inte.set_initial_value(y0, t0)

		Tc = self.contactTime()
		dt = Tc / 100.

		self.end = False
		time, delta, force = [], [], []

		def generator():
			while inte.successful() and (inte.t <= Tc):
				inte.integrate(inte.t + dt)
				yield inte.t + dt, inte.y, self.normalForce(inte.y[0]) + self.dissForce(inte.y[0], inte.y[1])

		for t, soln, f in generator():
			if self.end:
				break
			else:
				time.append(t)
				delta.append(soln)
				force.append(f)

				# for hysteretic models
				if hasattr(self, 'maxDisp'):
					if soln[0] >= self.maxDisp:
						self.maxDisp = soln[0]
					else:
						self.unloading = True

				# for hysteretic models
				if hasattr(self,'maxForce'):
					self.maxForce = max(f, self.maxForce)
						
		return np.array(time), np.array(delta), np.array(force)

	def contactRadius(self, delta):
		""" Returns the contact radius based on Hertzian or JKR models"""

		if type(delta) == np.ndarray:
			if (delta < 0).any():
				self.end = True
		elif delta < 0:
			self.end = True

		radius = self.radius
		contRadius = np.sqrt(delta * radius)
		
		if False:
			Gamma = self.cohesionEnergyDensity

			poiss = self.poissonsRatio
			yMod = self.youngsModulus
			yMod /= 2.0 * (1.0  - poiss )

			def jkr_disp(a, *args):
				delta, Gamma, yMod, radius = args
				return delta - a**2.0/radius + np.sqrt(2.0 * np.pi * Gamma * a / yMod)

			def jkr_jacob(a, *args):
				_, Gamma, yMod, radius = args
				return - 2.0 * a /radius + np.sqrt(np.pi * Gamma / (a * 2.0 * yMod))

			output = fsolve(jkr_disp, x0 = np.sqrt(delta * radius), args = (delta, Gamma, yMod, radius), full_output=True, fprime = jkr_jacob)
			contRadius = output[0]
			info = output[1]

			print(info)

		return contRadius

	def dissCoef(self):
		raise NotImplementedError('Not yet implemented')

	def springStiff(self):
		raise NotImplementedError('Not yet implemented')

	def normalForce(self):
		raise NotImplementedError('Not yet implemented')

	def dissForce(self):
		raise NotImplementedError('Not yet implemented')

	def numericalForce(self, time, delta):
		""" Returns the force used for numerical solvers """

		radius = self.radius

		Fn = self.normalForce(float(delta[0]))
		Fd = self.dissForce(float(delta[0]), float(delta[1]))
		
		Force = Fn + Fd

		mass = self.mass

		return np.array([delta[1], - Force / mass])

	def tangForce(self):
		raise NotImplementedError('Not yet implemented')

class SpringDashpot(Model):
	"""
	A class that implements the linear spring model for granular materials
	"""

	def __init__(self, **params):

		super(SpringDashpot, self).__init__(**params)

		if 'model-args' not in self.params:
			self.params['model-args'] = ('gran', 'model', 'hooke', 'tangential', 'history', 'rolling_friction', \
						'cdt', 'tangential_damping', 'on', 'limitForce', 'on', 'ktToKnUser', 'on') # the order matters here
		else:
			self.params['model-args'] = self.params['model-args']

	def springStiff(self, delta = None):
		""" Computes the spring constant kn for F = - kn * \delta
		"""
		poiss = self.poissonsRatio
		yMod = self.youngsModulus
		radius = self.radius
		mass = self.mass

		yMod /= 2.0 * (1.0  - poiss )

		v0 = self.characteristicVelocity

		return 16.0/15.0 * np.sqrt(radius) * yMod * (15.0 * mass \
			* v0 **2.0 / (16.0 * np.sqrt(radius) * yMod))**(1.0/5.0)

	def dissCoef(self):

		rest = self.coefficientRestitution
		poiss = self.poissonsRatio
		yMod = self.youngsModulus
		radius = self.radius

		mass = self.mass
		yMod /= 2.0 * (1.0  - poiss )
		v0 = self.characteristicVelocity

		kn = self.springStiff()
		loge = np.log(rest)

		return loge * np.sqrt(4.0 * mass * kn / (np.pi**2.0 + loge**2.0))

	def dissForce(self, delta, deltav):
		""" Returns the dissipative force """
		radius = self.radius

		return - self.dissCoef() * deltav

	def displacementAnalytical(self, dt = None):
		""" Computes the displacement based on an analytical solution """

		rest = self.coefficientRestitution
		poiss = self.poissonsRatio
		yMod = self.youngsModulus
		radius = self.radius
		mass = self.mass

		if dt is None:
			dt = self.contactTime()

		v0 = self.characteristicVelocity
		kn = self.springStiff()
		cn = self.dissCoef()

		const = np.sqrt(4.0 * mass * kn - cn**2.0) / mass

		return time, np.exp(- 0.5 * cn * time / mass) * 2.0 * v0 / const * np.sin(const * time / 2.0)

	def normalForce(self, delta):
		""" Returns the normal force based on Hooke's law: Fn = kn * delta """

		poiss = self.poissonsRatio
		yMod = self.youngsModulus
		radius = self.radius
		mass = self.mass

		kn = self.springStiff(radius)

		return kn * delta

class HertzMindlin(Model):
	"""
	A class that implements the linear spring model for granular materials
	"""

	def __init__(self, **params):
		super(HertzMindlin, self).__init__(**params)

		if 'model-args' not in self.params:
			self.params['model-args'] = ('gran', 'model', 'hertz', 'tangential', 'history', 'cohesion', 'sjkr', \
			'tangential_damping', 'on', 'limitForce', 'on') # the order matters here
		else:
			self.params['model-args'] = self.params['model-args']

	def springStiff(self, delta):
		""" Computes the spring constant kn for
			F = - kn * delta
		"""
		poiss = self.poissonsRatio
		yMod = self.youngsModulus
		radius = self.radius
		yEff = yMod * 0.5 / (1.0  - poiss )

		return 4.0 / 3.0 * yEff * np.sqrt(self.radius * delta)

	def normalForce(self, delta):
		""" Returns the Hertzian normal force"""

		return self.springStiff(delta) * delta

	def dissCoef(self, delta):
		""" Returns the dissipative force coefficient """
		rest = self.coefficientRestitution
		yMod = self.youngsModulus
		poiss = self.poissonsRatio
		yEff = yMod * 0.5 / (1.0  - poiss )

		radius = self.radius
		mass = self.mass

		contRadius = self.contactRadius(delta)

		return 2.0 * np.sqrt(5.0/6.0) * np.log(rest) / np.sqrt(np.log(rest)**2 + np.pi**2) * \
					np.sqrt(mass * 2 * yEff * contRadius)

	def dissForce(self, delta, deltav):
		""" Returns the dissipative force """
		return - self.dissCoef(delta) * deltav

class ThorntonNing(Model):
	"""
	A basic class that implements the Thornton elasto-plastic model
	"""

	def __init__(self, **params):

		super(ThorntonNing, self).__init__(**params)

		if 'model-args' not in self.params:
			self.params['model-args'] = ('gran', 'model', 'hysteresis_coh/thorn', \
					'tangential', 'history')
		else:
			self.params['model-args'] = self.params['model-args']

		# We check for the radius 1st since it can change in this model
		if hasattr(self, 'radius'):
			self.deltay = self.computeYieldDisp()
			self.maxDisp = 0
			self.maxForce = 0
			self.yieldForce = self.springStiff(self.deltay) * self.deltay
			self.unloading = False
			self.noCheck = False

	def computeYieldDisp(self):
		""" Computes the contact radius at the yield point """

		poiss = self.poissonsRatio
		yMod = self.youngsModulus
		yEff = self.youngsModulus / (2.0 * (1. - poiss))
		py = self.yieldPress

		return (0.5 * py * np.pi / yEff)**2 * self.radius

	def springStiff(self, delta):
		""" Computes the spring constant kn for
			F = - kn * delta
		"""
		poiss = self.poissonsRatio
		yMod = self.youngsModulus
		radius = self.radius	
		mass = self.mass
		yEff = yMod * 0.5 / (1.0  - poiss )

		return 4.0 / 3.0 * yEff * np.sqrt(radius * delta)

	def normalForce(self, delta):
		""" Returns the Hertzian-like normal force"""

		poiss = self.poissonsRatio
		yMod = self.youngsModulus
		yEff = yMod * 0.5 / (1.0  - poiss)

		if self.unloading:
			if not self.noCheck:
				self.noCheck = True
				self.radius = self.springStiff(self.maxDisp) * self.maxDisp / self.maxForce * self.radius
				self.deltap = self.maxDisp - (self.maxForce * 3 / (4. * yEff * np.sqrt(self.radius)))**(2.0/3.0)

			return self.springStiff(delta - self.deltap) * (delta - self.deltap)

		if delta < self.deltay:
			return self.springStiff(delta) * delta
		else:
			return self.springStiff(self.deltay) * self.deltay
	
	def dissForce(self, delta, deltav = None):
		""" Computes the piece-wise defined normal force based on Thornton's model """

		def compute(delta):
			py = self.yieldPress

			if self.unloading:
				return 0

			if delta > self.deltay:
				return  np.pi * py * self.radius * (delta - self.deltay)
			else:
				return 0

		if type(delta) == np.ndarray:
			out = []
			for disp in delta:
				out.append(compute(disp))

			return np.array(out)
		else:
			return compute(delta)
				
	def yieldVel(self):
		""" Returns the minimum velocity required for a colliding particle to undergo plastic deformation """

		poiss = self.poissonsRatio
		yEff = self.youngsModulus / (2.0 * (1. - poiss))
		density = self.density
		py = self.yieldPress

		return 1.56 * np.sqrt(py**5 / (yEff**4 * density))
