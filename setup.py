'''
Setup file for PyGran

Created on July 09, 2016

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

import os, sys, shutil
import subprocess
from setuptools import setup, find_packages
import glob, shutil
from distutils.command.install import install
from distutils.command.clean import clean

# Extract metadata from simulation._version
with open(os.path.join('src', 'PyGran', '_version.py'), 'r') as fp:
        for line in fp.readlines():
                if '__version__' in line:
                        __version__ = line.split('=')[-1].strip().strip("''")
                elif '__email__' in line:
                        __email__ = line.split('=')[-1].strip().strip("''")
                elif '__author__' in line:
                        __author__ = line.split('=')[-1].strip().strip("''")

try:
	from Cython.Build import cythonize
	import numpy
	optimal_list = cythonize("src/PyGran/analysis/PyGranAnalysis/core.pyx", compiler_directives={'language_level' : sys.version_info[0]})
	include_dirs = [numpy.get_include()]
except:
	print('Could not cythonize. Make sure Cython is properly installed.')
	optimal_list = []
	include_dirs = []

class Track(install):
	""" An install class that enables the tracking of installation/compilation progress """

	def execute(self, cmd, cwd='.'):
		popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True, cwd=cwd, shell=True)

		for stdout_line in iter(popen.stdout.readline, ""):
			yield stdout_line

		popen.stdout.close()
		return_code = popen.wait()

		if return_code:
			raise subprocess.CalledProcessError(return_code, cmd)

	def print_progress(self, iteration, prefix='', suffix='', decimals=1, total = 100):
		"""
		Call in a loop to create terminal progress bar
		@params:
			iteration   - Required  : current iteration (Int)
			total       - Required  : total iterations (Int)
			prefix      - Optional  : prefix string (Str)
			suffix      - Optional  : suffix string (Str)
			decimals    - Optional  : positive number of decimals in percent complete (Int)
			bar_length  - Optional  : character length of bar (Int)
		"""
		str_format = "{0:." + str(decimals) + "f}"
		percents = str_format.format(100 * (iteration / float(total)))
		sys.stdout.write('\r %s%s %s' % (percents, '%', suffix))
		sys.stdout.flush()

	def run(self):
		self.do_pre_install_stuff()

	def do_pre_install_stuff(self):
		raise NotImplementedError

class LIGGGHTS(Track):
		""" A class that enables the compilation of LIGGGHTS-PUBLIC from github """

		def do_pre_install_stuff(self):

			if os.path.exists('LIGGGHTS-PUBLIC'):
				print('Deleting ' + 'LIGGGHTS-PUBLIC')
				shutil.rmtree('LIGGGHTS-PUBLIC')

			self.spawn(cmd=['git', 'clone', 'https://github.com/CFDEMproject/LIGGGHTS-PUBLIC.git'])

			files = glob.glob('LIGGGHTS-PUBLIC/src/*.cpp')

			count = 0
			os.chdir('LIGGGHTS-PUBLIC/src')
			self.spawn(cmd=['make', 'clean-all'])

			print('Compiling LIGGGHTS as a shared library\n')

			for path in self.execute(cmd='make auto'):
				count +=1
				self.print_progress(count, prefix = 'Progress:', suffix = 'Complete', total = len(files) * 2.05)

			self.spawn(cmd=['make', '-f', 'Makefile.shlib', 'auto'])
			sys.stdout.write('\nInstallation of LIGGGHTS-PUBLIC complete\n')
			os.chdir(os.path.join('..','..'))

class Clean(clean):

	def run(self):
		for ddir in ['build', 'dist', 'PyGran.egg-info']: 
			if os.path.isdir(ddir):
				print('Deleting ' + os.path.abspath(ddir))
				shutil.rmtree(ddir)

		super().run()

setup(
	name = "PyGran",
	version = __version__,
	author = __author__,
	author_email = __email__,
	description = ("A DEM toolkit for rapid quantitative analysis of granular/powder systems"),
	license = "GNU v2",
	keywords = "Discrete Element Method, Granular Materials",
	url = "https://github.com/Andrew-AbiMansour/PyGran",
	packages=find_packages('src'),
	package_dir={'PyGran': os.path.join('src','PyGran')},
	include_package_data=True,
	install_requires=['numpy', 'scipy'],
	extras_require={'extra': ['vtk', 'mpi4py', 'Pillow', 'pytest', 'pytest-cov', 'codecov']},
	long_description='A DEM toolbox for rapid quantitative analysis of granular/powder systems. See http://www.pygran.org.',
	classifiers=[
			"Development Status :: 4 - Beta",
			"Topic :: Utilities",
			"License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
			"Programming Language :: Python :: 2.7",
			"Programming Language :: Python :: 3",
			"Programming Language :: Python :: 3.4",
			"Programming Language :: Python :: 3.5",
			"Programming Language :: Python :: 3.6",
			"Programming Language :: Python :: 3.7",
			"Programming Language :: C",
			"Operating System :: POSIX :: Linux"
	],

	cmdclass={'build_liggghts': LIGGGHTS, 'clean': Clean},
	zip_safe=False,
	ext_modules=optimal_list,
	include_dirs=include_dirs
)
