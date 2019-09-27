'''
  Created on March 30, 2016
  @author: Andrew Abi-Mansour

  This is the 

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
