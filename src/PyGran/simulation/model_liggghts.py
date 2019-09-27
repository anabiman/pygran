'''
A module for creating contact models on the fly for LIGGGHTS

Created on July 1, 2016

Author: Andrew Abi-Mansour

This is the::

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

from ..tools import find
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

def _parse(exp):

	keywords = {'Yeff': 'Yeff[itype][jtype]', 'Geff': 'Geff[itype][jtype]', 'restLogChosen': 'coeffRestLogChosen', \
	'restLog': 'coeffRestLog', 'meff': 'meff', 'charVel': 'charVel', 'reff': 'reff', 'deltan': 'sidata.deltan', \
	'vn': 'sidata.vn', 'PI': 'M_PI'}

	for keyword in keywords:
		exp = exp.replace(keyword, "{{{}}}".format(keyword))

	return exp.format(**keywords)

def register(**args):
	""" Generates a c++ header file for a contact model and compiles it during runtime.

	@stiffness: analytical form (str) of the stiffness = force / deltan
	@viscosity: analytical form (str) of the viscosity term (force = viscosity * vn).

	Material parameters that can be used:
	Yeff: effective Young's modulus
	Geff: effective Shear modulus
	meff: effective mass
	reff: effective radius
	charVel: characteristic impact velocity
	restLog: log of the coefficient of restitution
	deltan: normal displacement
	vn: velocity of the normal displacement
	kn: stiffness
	PI: constant (3.14 ...)

	e.g. register(name='my_model',
		stiffness = '6./15.*sqrt(reff)*(Yeff)*pow(15.*meff*charVel*charVel/(16.*sqrt(reff)*Yeff),0.2)',
		viscosity = 'sqrt(4.*meff*kn*restLogChosen*restLogChosen/(restLogChosen*restLogChosen+PI*PI))')

		produces a template header file for the spring-dashpot model called 'my_model'.

	"""

	# Make sure everything is done on the master processor
	try:
		rank = MPI.COMM_WORLD.Get_rank()
	except:
		rank = 0

	if not rank:

		if 'stiffness' not in args:
			raise('An analytical equation for stiffness must be specified.')
		else:
			args['stiffness'] = _parse(args['stiffness'])

		if 'viscosity' not in args:
			raise('An analytical equation for viscosity must be specified.')
		else:
			args['viscosity'] = _parse(args['viscosity'])

		args['name_lower'] = args['name']
		args['name'] = args['name'].upper()

		if 'ktToKn' not in args:
			args['ktToKn'] = 2.0/7.0

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
