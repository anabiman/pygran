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

from PyGran import Simulator, Analyzer, Tools
import os

try:
	from PyGran import Visualizer
except:
	print("WARNING: Visualization not supported. Make sure wxwidgets and vtk libraries are properly installed on this system.")

def run(program):
	""" Launches an executable program available in the PATH env variable """
	paths = os.environ['PATH']

	for path in paths.split(':'):
		found = Tools.find(program, path)

		if found:
			print('Launching {}'.format(found))
			os.system(found + ' &')
			return 0

	print('Could not find {} in {}'.format(program, paths))
	return 1
