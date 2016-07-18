# !/usr/bin/python
# -*- coding: utf8 -*- 
# -----------------------------------------------------------------------
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

'''
Created on March 30, 2016
@author: Andrew Abi-Mansour

Center for Materials Sci. & Eng.,
Merck Inc., West Point
'''

from DEM import *
from glob import glob as _glob

class _findEngines:
	def __init__(self):

		ext = __file__.split('__init__.py')[0]
		pyFiles = _glob(ext + '*.py')

		for file in pyFiles:
			if file != ext + 'DEM.py' and file != ext + '__init__.py':
				
				engine, _ = file.split(ext)[1].split('.py')
				setattr(self, engine, engine)

engines = _findEngines()