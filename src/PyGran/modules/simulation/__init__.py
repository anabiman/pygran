'''
Created on March 30, 2016
@author: Andrew Abi-Mansour
'''

# !/usr/bin/python
# -*- coding: utf8 -*- 
# -----------------------------------------------------------------------
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

from .dem import *
from glob import glob as _glob
from .model_liggghts import *

class _findEngines:
	def __init__(self):
		# Any engine module *must* follow the naming convention: engine_foo.py
		# If the engine is found, it will be linked via setattr to be imported
		# in DEM.py as PyGran.simulation.engine_foo. The engine is set by the user
		# as DEM.simulation.engines.foo

		_dir, _ = __file__.split('__init__.py')
		pyFiles = _glob(_dir + '*.py')

		for file in pyFiles:
			_, fileName = file.split(_dir)

			if fileName.startswith('engine_'):	
				engine, _ = fileName.split('.py')
				setattr(self, engine.split('engine_')[1], engine)

engines = _findEngines()
