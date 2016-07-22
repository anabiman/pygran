'''
Created on July 1, 2016
@author: Andrew Abi-Mansour

Center for Materials Sci. & Eng.,
Merck Inc., West Point
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

class Model:
	def __init__(self, **params):

		self.params = params
		self.params['nSS'] = 0

		if 'SS' in self.params:
			self.params['nSS'] += len(self.params['SS'])

		if 'mesh' in self.params:
			self.params['nSS'] += len(self.params['mesh'])

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

		if 'restart' not in self.params:
			self.params['restart'] = (5000, 'restart', 'restart.binary', False)

		if 'dump_modify' not in self.params:
			self.params['dump_modify'] = ('append', 'yes')

		if 'nSim' not in self.params:
			self.params['nSim'] =  1

		self.materials = {}

		if 'materials' in self.params:
			for item in self.params['materials']:
				if self.params['materials'][item][1] == 'peratomtype':
					self.materials[self.params['materials'][item][0]] = np.array([np.float(it) for \
						it in self.params['materials'][item][3:]]).mean()
				else:	
					self.materials[params['materials'][item][0]] = np.array([np.float(it) for \
						it in params['materials'][item][2:]]).mean()

		if 'SS' in self.params:
			self.SS = params['SS']
			self.radius = np.array([radius[1] for radius in params['radius']])
			self.mass = np.array([4.0/3.0 * np.pi * self.radius[i]**3 * self.SS[i]['density'] for \
									i in range(len(self.SS))])

			# Estimate the allowed sim timestep
			self.params['dt'] = (0.25 * self.contactTime()).min()
		else:
			print 'Warning: no components found in your supplied dictionary!'

	def contactTime(self):
		raise NotImplementedError('Not yet implemented')

	def overlap(self):
		raise NotImplementedError('Not yet implemented')

	def dissCoef(self):
		raise NotImplementedError('Not yet implemented')

	def springStiff(self):
		raise NotImplementedError('Not yet implemented')

	def normalForce(self):
		raise NotImplementedError('Not yet implemented')

	def tangForce(self):
		raise NotImplementedError('Not yet implemented')

class SpringDashpot(Model):
	"""
	A class that implements the linear spring model for granular materials
	"""

	def __init__(self, **params):

		params['materials']['cVel'] = ('characteristicVelocity', 'scalar', '0.2', '0.2')

		Model.__init__(self, **params)

		if 'model-args' not in self.params:
			self.params['model-args'] = ('gran', 'model', 'hooke', 'tangential', 'history', 'rolling_friction', \
						'cdt', 'tangential_damping', 'on', 'limitForce', 'on', 'ktToKnUser', 'on') # the order matters here
		else:
			self.params['model-args'] = self.params['model-args']

	def springStiff(self, radius = None, yMod = None, poiss = None, mass = None, v0 = None):
		""" Computes the spring constant kn for F = - kn * \delta
		"""
		if poiss is None:
			poiss = self.materials['poissonsRatio']

		if yMod is None:
			yMod = self.materials['youngsModulus']
			yMod /= 2.0 * (1.0  - poiss )

		if v0 is None:
			v0 = self.materials['characteristicVelocity']

		if mass is None:
			mass = self.mass

		if radius is None:
			radius = self.radius

		return 16.0/15.0 * np.sqrt(radius) * yMod * (15.0 * mass \
			* v0 **2.0 / (16.0 * np.sqrt(radius) * yMod))**(1.0/5.0)

	def dissCoef(self, radius = None, yMod = None, poiss = None, mass = None, v0 = None, rest = None):

		if rest is None:
			rest = self.materials['coefficientRestitution']

		if mass is None:
			mass = self.mass

		kn = self.springStiff(radius, yMod, poiss, mass, v0)
		loge = np.log(rest)

		return loge * np.sqrt(4.0 * mass * kn / (np.pi**2.0 + loge**2.0))

	def contactTime(self, mass = None, rest = None, kn = None):

		if kn is None:
			 kn = self.springStiff()

		if mass is None:
			mass = self.mass

		if rest is None:
			rest = self.materials['coefficientRestitution']

		return np.sqrt(mass * (np.pi**2.0 + np.log(rest)) / kn) 

	def overlap(self, delta = None, radius = None, yMod = None, poiss = None, mass = None, v0 = None, rest = None):
		""" Computes the overlap distance """
		kn = self.springStiff(radius, yMod, poiss, mass, v0)
		cn = self.dissCoef(radius, yMod, poiss, mass, v0, rest)
		dt = self.contactTime(mass, rest, kn)

		time = np.arange(0, dt*100, dt)

		if v0 is None:
			v0 = self.materials['characteristicVelocity']

		if mass is None:
			mass = self.mass

		if radius is None:
			radius = self.radius

		const = np.sqrt(4.0 * mass * kn - cn**2.0) / mass
		return np.exp(- 0.5 * cn * time / mass) * 2.0 * v0 / const * np.sin(const * time / 2.0)

	def normalForce(self, delta = None, radius = None, yMod = None, poiss = None, mass = None, v0 = None):
		if delta is None:
			delta = self.overlap()

		return - self.springStiff(radius, yMod, poiss, mass, v0) * delta

class HertzMindlin(Model):
	"""
	A class that implements the linear spring model for granular materials
	"""

	def __init__(self, **params):
		Model.__init__(self, **params)

		if 'model-args' not in self.params:
			self.params['model-args'] = ('gran', 'model', 'hertz', 'tangential', 'history', 'rolling_friction', \
						'cdt', 'tangential_damping', 'on', 'limitForce', 'on') # the order matters here
		else:
			self.params['model-args'] = self.params['model-args']

	def springStiff(self, radius = None, yMod = None, mass = None, v0 = None):
		""" Computes the spring constant kn for 
			F = - kn * \delta
		"""
		if yMod is None:
			yMod = self.materials['youngsModulus']

		if mass is None:
			mass = self.mass

		if radius is None:
			radius = self.radius

		pass

	def contactTime(self, mass = None, cR = None, kn = None):
		return 2.5e-5 * np.ones(2)

	def normalForce(self, ):
		pass