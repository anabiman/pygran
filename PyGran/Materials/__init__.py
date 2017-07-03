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
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.


# -------------------------------------------------------------------------

glass = {
		'youngsModulus': ('youngsModulus', 'peratomtype', '63e9'),
		'poissonsRatio': ('poissonsRatio', 'peratomtype', '0.24'),
		'coefficientFriction': ('coefficientFriction', 'peratomtypepair', '0.5'),
		'coefficientRollingFriction': ('coefficientRollingFriction', 'peratomtypepair', '0.1'),
		'cohesionEnergyDensity': ('cohesionEnergyDensity', 'peratomtypepair', '0.05'),
		'coefficientRestitution': ('coefficientRestitution', 'peratomtypepair', '0.9'),
		'coefficientRollingViscousDamping': ('coefficientRollingViscousDamping', 'peratomtypepair', '0.1'),
		'yieldPress': ('yieldPress', 'peratomtype', '62e9'),
		'characteristicVelocity': ('characteristicVelocity', 'scalar', '0.1'),
		'density': 2500.0
		}

stearicAcid = {
		'youngsModulus': ('youngsModulus', 'peratomtype', '4.15e7'),
		'poissonsRatio': ('poissonsRatio', 'peratomtype', '0.25'),
		'coefficientFriction': ('coefficientFriction', 'peratomtypepair', '0.5'),
		'coefficientRollingFriction': ('coefficientRollingFriction', 'peratomtypepair', '0.1'),
		'cohesionEnergyDensity': ('cohesionEnergyDensity', 'peratomtypepair', '0.033'),
		'coefficientRestitution': ('coefficientRestitution', 'peratomtypepair', '0.9'),
		'coefficientRollingViscousDamping': ('coefficientRollingViscousDamping', 'peratomtypepair', '0.1'),
		'yieldPress': ('yieldPress', 'peratomtype', '2.2e6'),
		'characteristicVelocity': ('characteristicVelocity', 'scalar', '0.1'),
		'density': 997.164
		}
