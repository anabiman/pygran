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

def dictToTuple(**args):
	""" Converts a dictionary (args) to a tuple 

	e.g. if arg = {'key_1': value1, 'key_2': value2}
	this function returns ('key_1', value1, 'key_2', value2)

	if value1 is a tuple, it is broken into a string
	e.g. if arg = {'key_1': (1,2,3)}
	this function returns ('key_1', '1 2 3')

	"""

	keys = args.keys()
	vals = args.values()

	tup = ()
	for pair in zip(keys, vals):
		if isinstance(pair[1], tuple):
			pair = (pair[0], ('{} ' * len(pair[1])).format(*pair[1]).strip())
		tup = tup + pair

	return tup

def convert(unitso, unitsf):
	""" Generic function that converts length/time/mass from one unit system to another

	@unitso: string specifying unit system to convert from
	@unitsf: string specifying unit system to convert to

	returns the new unit system as a dictionary
	"""

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

def find(fname, path):
	""" Finds a filename (fname) along the path `path' 

	@fname: string specifying filename
	@path: string specifying the search path 

	returns the absolute path of the fname (if found) as a string
	"""
	for root, dirs, files in os.walk(path):
		if fname in files:
			return os.path.join(root, fname)

	return None

def run(program):
	""" Unix only: launches an executable program available in the PATH environment variable.
	
	@program: string specifying the executable to search for

	returns 0 if successful and 1 otherwise. """
	paths = os.environ['PATH']

	for path in paths.split(':'):
		found = Tools.find(program, path)

		if found:
			print('Launching {}'.format(found))
			os.system(found + ' &')
			return 0

	print('Could not find {} in {}'.format(program, paths))
	return 1

def configure(path, version=None, src=None):
	""" configures PyGran to use a specific LIGGGHTS library """
	_setLIGGGHTS(path, version, src)

def _setLIGGGHTS(path, version=None, src=None):
	""" Write libliggghts path to .config file 

	@[version]: a set of numbers &/or characters indicating the version of the library, e.g. 1.5a
	"""

	wdir, _ = os.path.abspath(__file__).split(os.path.basename(__file__))

	with open(wdir + '../.config', 'w') as fp:

		fp.seek(0,0)
		fp.write('library=' + path)

		if src:
			fp.write('\nsrc=' + src)

		if version:
			fp.write('\nversion=' + version)

def _findEngines(engine):
	""" Searches for and lists all available libraries for a specific engine """

	engines = [os.path.join(root, engine) for root, dirs, files in os.walk('/') if engine in files]

	if engines:
		print('Engine(s) found:')
		for engine in engines:
			print(engine)
	else:
		print('No engines found.')

	return engines

def _setConfig(wdir, engine):
	""" Reads/writes libliggghts to .config file """

	library, src, __version__ = None, None, None
	file = wdir + '../.config'

	if os.path.isfile(file):

		if os.stat(file).st_size: # file not empty

			with open(wdir + '../.config', 'r+') as fp: # r+ is for reading and writing
				for line in fp.readlines():
					if 'library=' in line:
						library = line.split('=')[-1].rstrip()
					elif 'src' in line:
						src = line.split('=')[-1].rstrip()
					elif 'version' in line:
						__version__ = float(line.split('=')[-1].rstrip())

			# Make sure the library exists; else, find it somewhere else
			if not os.path.isfile(library):
				library = find('lib' + engine+ '.so', '/')
				_setLIGGGHTS(library)
				print('WARNING: Could not find user-specified library. Will use {} instead ...'.format(library))

			return library, src, __version__

	with open(wdir + '../.config', 'w') as fp:
		library = find('lib' + engine + '.so', '/')

		if library:
			print('No config file found. Creating one for {} in {}'.format(library, os.path.abspath(file)))
			_setLIGGGHTS(library)
		else:
			print('No installation of {} was found. Make sure your selected DEM engine is properly installed first.'.format(engine))
			print('PyGran looked for ' + 'lib' + engine + '.so' + '. If the file exists, make sure it can be executed by the user.')
			sys.exit()

	return library, src, __version__
