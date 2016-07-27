'''
Created on July 18, 2016
@author: Andrew Abi-Mansour

Center for Materials Sci. & Eng.,
Merck Inc., West Point
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
		'yMod': ('youngsModulus', 'peratomtype', '7.1e7', '7.1e7'),
		'pRatio': ('poissonsRatio', 'peratomtype', '0.22', '0.22'),
		'cFric': ('coefficientFriction', 'peratomtypepair', '2', '0.5', '0.5', '0.5', '0.5'),
		'cRollFric': ('coefficientRollingFriction', 'peratomtypepair', '2', '5e-4', '5e-4', '5e-4', '5e-4'),
		'cEnerDen': ('cohesionEnergyDensity', 'peratomtypepair', '2', '1.0', '1.0', '1.0', '1.0'),
		'cRest': (('coefficientRestitution', 'peratomtypepair', '2', '0.9', '0.9', '0.9', '0.9'))
		}