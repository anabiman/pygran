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
# -----------------------------------------------------------------------

'''
Created on April 25, 2016
@author: Andrew Abi-Mansour
'''

import os

# Conversion unit systems
# si.: meters, seconds, kilograms
# cgs: centimeters, seconds, grams
# micro: microns, seconds, micrograms
# nano: nanometers, seconds, nanograms 

conversion = {'si': {'distance': [1, 'm'], 'time': [1, 's'], 'mass': [1, 'kg']}, \
			'cgs': {'distance': [1e-2, 'cm'],'time': [1, 's'], 'mass': [1e-3, 'g']}, \
			'micro': {'distance': [1e-6, '$\mu m$'], 'time': [1, 's'], 'mass': [1e-9, '$\mu g$']}, \
			'nano': {'distance': [1e-9, 'nm'], 'time': [1, 's'], 'mass': [1e-12, 'ng']}
			}

def convert(unitso, unitsf):

	if unitso in conversion:
		if unitsf in conversion:
			conv = conversion[unitso]

			for key in conv:
				conv[key][0] /= conversion[unitsf][key][0]

			return conv
		else:
			raise ValueError('Input unit system not supported: {}'.format(unitsf))
	else:
		raise ValueError('Input unit system not supported: {}'.format(unitso))

def find(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)

    return None

def run(program):
	""" Unix only: launches an executable program available in the PATH environment variable.
	Returns 0 if successful and 1 otherwise. """
	paths = os.environ['PATH']

	for path in paths.split(':'):
		found = Tools.find(program, path)

		if found:
			print('Launching {}'.format(found))
			os.system(found + ' &')
			return 0

	print('Could not find {} in {}'.format(program, paths))
	return 1