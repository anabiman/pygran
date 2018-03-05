'''
Created on July 18, 2016
@author: Andrew Abi-Mansour
'''

# !/usr/bin/python
# -*- coding: utf8 -*-
# -------------------------------------------------------------------------
#
#   Python module for storing material properties
#
# --------------------------------------------------------------------------
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 2 of the License, or
#   (at your option) any later version.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.


# -------------------------------------------------------------------------

glass = {
		'youngsModulus': 63e9,
		'poissonsRatio': 0.24,
		'coefficientFriction': 0.5,
		'coefficientRollingFriction': .0,
		'cohesionEnergyDensity': 0.05,
		'coefficientRestitution': 0.9,
		'coefficientRollingViscousDamping': 0.1,
		'yieldPress': 62e9,
		'characteristicVelocity': 0.1,
		'density': 2500.0
		}

stearicAcid = {
		'youngsModulus': 4.15e7,
		'poissonsRatio': 0.25,
		'coefficientFriction': 0.5,
		'coefficientRollingFriction': 0.0,
		'cohesionEnergyDensity': 0.033,
		'coefficientRestitution': 0.9,
		'coefficientRollingViscousDamping': 0.1,
		'yieldPress': 2.2e6,
		'characteristicVelocity': 0.1,
		'density': 997.164
		}


organic = {
		'youngsModulus': 1e7,
		'poissonsRatio': 0.25,
		'coefficientFriction': 0.5,
		'coefficientRollingFriction': 0.0,
		'cohesionEnergyDensity': 0.05,
		'coefficientRestitution': 0.9,
		'coefficientRollingViscousDamping': 0.1,
		'yieldPress': 2.2e6,
		'characteristicVelocity': 0.1,
		'density': 1000.0
		}

def LIGGGHTS(**material):
	""" Transform a PyGran material database into a LIGGGHTS material dictionary """

	for key in material:
		if key is 'youngsModulus' or key is 'poissonsRatio' or key is 'yieldPress':
			material[key] = (key, 'peratomtype', str(material[key]))
		elif key is 'coefficientFriction' or key is 'coefficientRollingFriction' or key is 'cohesionEnergyDensity' \
			or key is 'coefficientRestitution' or key is 'coefficientRollingViscousDamping':
			material[key] = (key, 'peratomtypepair', str(material[key]))
		elif key is 'characteristicVelocity':
			material[key] = (key, 'scalar', str(material[key]))

	return material