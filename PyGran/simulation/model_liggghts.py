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
#   the Free Software Foundation, either version 2 of the License, or
#   (at your option) any later version.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

# ----------------------------------------------

from PyGran.tools import find
try:
	from mpi4py import MPI
except:
	pass # no MPI, no problem
import os, glob

def _find_number_models(src_dir, mtype='normal'):
	""" Finds the total number of contact models available in the liggghts src dir 

	@src_dir: directory to search the contact model header files in
	@[mtype]: 'normal' (default) or 'tangential' contact models to search for

	"""

	nModels = []

	for file in glob.glob(src_dir + '{}_*'.format(mtype)):

		with open(file, 'r') as fp:
			for line in fp:
				if '{}_MODEL'.format(mtype.upper()) in line and (not line.startswith('#')):
					nModels.append(line.rstrip().split(',')[-1].split(')')[0])
					break
	return nModels

def register(**args):
	""" Generates a c++ header file for a contact model and compiles it during runtime. """

	# Make sure everything is done on a single core
	try:
		rank = MPI.COMM_WORLD.Get_rank()
	except:
		rank = 0

	if not rank:

		args['name_lower'] = args['name']
		args['name'] = args['name'].upper()

		if 'mtype' not in args:
			args['mtype'] = 'normal'

		_dir, _ = __file__.split(__name__.split('PyGran.simulation.')[-1] +'.py')

		if 'src_dir' not in args:
			if os.path.isfile(_dir + '../.config'):
				with open(_dir + '../.config', 'r+') as fp:
					fp.readline(); fp.readline();
					args['src_dir'] = fp.readline().split('=')[-1].rstrip()
			else:
				raise ValueError('Could not find a LIGGGHTS source code directory. Specify this with src_dir=path/to/LIGGGHTS/src.')

		# find contact model number
		nModels = _find_number_models(mtype=args['mtype'], src_dir=args['src_dir'])

		for number in range(100): # more than 100 contact models? WTF! TODO: Make this better automated
			if str(number) not in nModels:
				args['number'] = number
				print(number, nModels)
				break

		with open(_dir + 'model_template.h', 'r') as fp:
			lines = fp.readlines()

		with open('{mtype}_model_{name_lower}.h'.format(**args), 'w') as fp:

			for line in lines:
				if '{name}' in line or '{name_lower}' in line or '{stiffness}' in line or '{viscosity}' in line or '{number}' in line:
					line = line.format(**args)

				fp.write(line)