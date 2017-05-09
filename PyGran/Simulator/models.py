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

class Model:
	def __init__(self, **params):

		self.params = params
		self.params['nSS'] = 0

		if 'SS' in self.params:
			self.params['nSS'] += len(self.params['SS'])

		# Treat any mesh as an additional component
		if 'mesh' in self.params:
			self.params['nSS'] += len(self.params['mesh'])
			for mesh in self.params['mesh']:
				self.params['SS']  += ({'material': self.params['mesh'][mesh]['material']},)

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

		# Expand material properties based on number of components
		if 'SS' in self.params:
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
					else:
						ss[item] = ss['material'][item]

			for item in self.materials:
				if type(ss['material'][item]) is not float:
					for ss in self.params['SS']:
						if ss['material'][item][1] == 'peratomtype': 
							self.materials[item] =  self.materials[item] + (('{}').format(ss['material'][item][2]),)

						elif ss['material'][item][1] == 'peratomtypepair':
							# assume the arithmetic mean suffices for estimating binary properties
							for nss in range(self.params['nSS']):
								prop = 0.5 * (float(ss['material'][item][2]) + float(self.params['SS'][nss]['material'][item][2]))
								self.materials[item] =  self.materials[item] + (('{}').format(prop),)

						# we should set this based on species type
						elif ss['material'][item][1] == 'scalar':
							self.materials[item] = self.materials[item] + (('{}').format(ss['material'][item][2]),)

						else:
							print('Error: Material database flawed.')
							sys.exit()

			self.params['materials'] = self.materials

		else:
			print('Warning: No components found in your supplied dictionary! Proceeding ...')

		if 'dt' not in self.params:
			# Estimate the allowed sim timestep
			try:
				self.params['dt'] = (0.25 * self.contactTime()).min()
			except:
				self.params['dt'] = 1e-6

				if 'model' in self.params:
					print('Model {} does not yet support estimation of contact period. Using a default value of {}'.format(self.params['model'], self.params['dt']))

	def contactTime(self):
		raise NotImplementedError('Not yet implemented')

	def displacement(self, v0 = None, dt = None):

		if v0 is None:
			v0 = self.materials['characteristicVelocity']

		if dt is None:
			dt = 1e-6

		y0 = np.array([0, v0])
		t0 = .0

		inte = ode(self.normalForce)
		inte.set_f_params(*(v0,))
		inte.set_integrator('dopri5')
		inte.set_initial_value(y0, t0)

		time, soln = [], []
		Tc = self.contactTime()

		while inte.successful() and inte.t <= Tc:
			inte.integrate(inte.t + dt)
			time.append(inte.t + dt)
			soln.append(inte.y)

		return np.array(time), np.array(soln)

	def contactRadius(self, delta):
		""" Returns the contact radius based on purely Hertzian or the JKR models"""

		radius = self.materials['radius']
		ontRadius = np.sqrt(delta * radius)

		if 'cohesionEnergyDensity' in self.materials:
			Gamma = self.materials['cohesionEnergyDensity']

			poiss = self.materials['poissonsRatio']
			yMod = self.materials['youngsModulus']
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

	def tangForce(self):
		raise NotImplementedError('Not yet implemented')

class SpringDashpot(Model):
	"""
	A class that implements the linear spring model for granular materials
	"""

	def __init__(self, **params):

		Model.__init__(self, **params)

		if 'model-args' not in self.params:
			self.params['model-args'] = ('gran', 'model', 'hooke', 'tangential', 'history', 'rolling_friction', \
						'cdt', 'tangential_damping', 'on', 'limitForce', 'on', 'ktToKnUser', 'on') # the order matters here
		else:
			self.params['model-args'] = self.params['model-args']

	def springStiff(self, v0 = None):
		""" Computes the spring constant kn for F = - kn * \delta
		"""
		poiss = self.materials['poissonsRatio']
		yMod = self.materials['youngsModulus']
		radius = self.materials['radius']
		mass = self.materials['mass']

		yMod /= 2.0 * (1.0  - poiss )

		if v0 is None:
			v0 = self.materials['characteristicVelocity']

		return 16.0/15.0 * np.sqrt(radius) * yMod * (15.0 * mass \
			* v0 **2.0 / (16.0 * np.sqrt(radius) * yMod))**(1.0/5.0)

	def dissCoef(self, v0 = None):

		rest = self.materials['coefficientRestitution']
		poiss = self.materials['poissonsRatio']
		yMod = self.materials['youngsModulus']
		radius = self.materials['radius']
		mass = self.materials['mass']

		yMod /= 2.0 * (1.0  - poiss )

		if v0 is None:
			v0 = self.materials['characteristicVelocity']

		kn = self.springStiff(v0)
		loge = np.log(rest)

		return loge * np.sqrt(4.0 * mass * kn / (np.pi**2.0 + loge**2.0))

	def contactTime(self):
		""" Computes the characteristic collision time """

		rest = self.materials['coefficientRestitution']
		poiss = self.materials['poissonsRatio']
		yMod = self.materials['youngsModulus']
		radius = self.materials['radius']
		mass = self.materials['mass']

		kn = self.springStiff()

		return np.sqrt(mass * (np.pi**2.0 + np.log(rest)) / kn) 

	def displacementExact(self, v0 = None, dt = None):
		""" Computes the displacement based on an analytical solution """

		rest = self.materials['coefficientRestitution']
		poiss = self.materials['poissonsRatio']
		yMod = self.materials['youngsModulus']
		radius = self.materials['radius']
		mass = self.materials['mass']

		kn = self.springStiff(v0)
		cn = self.dissCoef(v0)
		
		if dt is None:
			dt = self.contactTime()

		if v0 is None:
			v0 = self.materials['characteristicVelocity']

		const = np.sqrt(4.0 * mass * kn - cn**2.0) / mass

		return time, np.exp(- 0.5 * cn * time / mass) * 2.0 * v0 / const * np.sin(const * time / 2.0)

	def normalForce(self, time, delta, v0 = None):
		""" Returns the normal force based on Hooke's law: Fn = kn * delta """

		poiss = self.materials['poissonsRatio']
		yMod = self.materials['youngsModulus']
		radius = self.materials['radius']
		mass = self.materials['mass']

		kn = self.springStiff(v0)

		if len(delta) > 1:
			return np.array([delta[1], - kn * delta[0] / mass])
		else:
			return - kn * delta

class HertzMindlin(Model):
	"""
	A class that implements the linear spring model for granular materials
	"""

	def __init__(self, **params):
		Model.__init__(self, **params)

		if 'model-args' not in self.params:
			self.params['model-args'] = ('gran', 'model', 'hertz', 'tangential', 'history', 'cohesion', 'sjkr', \
			'rolling_friction', 'cdt', 'tangential_damping', 'on', 'limitForce', 'on') # the order matters here
		else:
			self.params['model-args'] = self.params['model-args']

	def springStiff(self, delta):
		""" Computes the spring constant kn for 
			F = - kn * \delta
		"""
		poiss = self.materials['poissonsRatio']
		yMod = self.materials['youngsModulus']
		radius = self.materials['radius']
		mass = self.materials['mass']
		
		yEff = yMod * 0.5 / (1.0  - poiss )

		return 4.0 / 3.0 * yEff * np.sqrt(radius * delta) 

	def contactTime(self):
		""" Estimate contact time based on a spring model """

		poiss = self.materials['poissonsRatio']
		yMod = self.materials['youngsModulus'] / (2.0 * (1 - poiss))
		radius = self.materials['radius']
		mass = self.materials['mass']
		v0 = self.materials['characteristicVelocity']

		kn = 16.0/15.0 * np.sqrt(radius) * yMod * (15.0 * mass \
			* v0 **2.0 / (16.0 * np.sqrt(radius) * yMod))**(1.0/5.0)

		return np.sqrt(mass * (np.pi**2.0 / kn)) 

	def normalForce(self, time, delta, v0 = None):
		""" Computes the Hertzian normal force"""

		poiss = self.materials['poissonsRatio']
		yMod = self.materials['youngsModulus']
		radius = self.materials['radius']
		mass = self.materials['mass']
		
		yEff = yMod * 0.5 / (1.0  - poiss )

		if delta[0]:
			contRadius = self.contactRadius(delta[0])[0]
		else:
			contRadius = .0

		if len(delta) > 1:
			Fn = np.array([delta[1], - 4.0/3.0 * yEff / mass * contRadius**3.0 / radius])
		else:
			Fn = - 4.0/3.0 * yEff * contRadius**3.0 / radius

		if 'cohesionEnergyDensity' in self.materials:
			Gamma = self.materials['cohesionEnergyDensity']
			Fn -= np.sqrt(8.0 * np.pi * yEff * Gamma * contRadius**3.0)

		return Fn

class Hysteresis(Model):
	"""
	A basic class that implements the Thornton elasto-plastic model
	"""

	def __init__(self, **params):
		Model.__init__(self, **params)

		if 'name' not in params:
			name = 'thorn'
		else:
			name = params['name']

		if 'model-args' not in self.params:
			self.params['model-args'] = ('gran', 'model', 'hysteresis_coh/{}'.format(name), \
					'tangential', 'history', 'rolling_friction', 'epsd')
		else:
			self.params['model-args'] = self.params['model-args']

		if 'mass' in params: # very hackish way to detect if we're doing analysis or running simulation
			py = self.materials['yieldPress']
			poiss = self.materials['poissonsRatio']

			yMod = self.materials['youngsModulus'] / (2.0 * (1. - poiss))
			radius = self.materials['radius']

			self.deltay = (np.pi * py / (2.0 * yMod))**2.0 * radius
			self.deltam = .0

	def springSitff(self, delta):

		poiss = self.materials['poissonsRatio']
		yMod = self.materials['youngsModulus']
		radius = self.materials['radius']
		mass = self.materials['mass']
		
		yEff = yMod * 0.5 / (1.0  - poiss )

		return 4.0 / 3.0 * yEff * np.sqrt(radius * delta)

	def normalForce(self, time, delta, v0 = None):
		""" Computes the piece-wise defined normal force based on Thornton's model """

		deltan, deltav = delta

		poiss = self.materials['poissonsRatio']
		yEff = self.materials['youngsModulus'] / (2.0 * (1. - poiss))
		radius = self.materials['radius']
		mass = self.materials['mass']
		py = self.materials['yieldPress']

		if deltam < deltay:
			if len(delta) > 1:
				return np.array([deltav, -self.springSitff(deltan) * deltan / mass]) 
			else:
				return -self.springSitff(deltan) * deltan
		else:
			Fy = - self.springSitff(deltay) * deltay 

			if deltav > 0:
				if deltan >= deltam:
					self.Fmax = Fy + np.pi * py * radius * (deltan - deltay)
					return self.Fmax

				else:
					pass
