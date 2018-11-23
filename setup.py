#  -*- coding: utf8 -*-
'''
  Created on July 9, 2016
  @author: Andrew Abi-Mansour

  This is the 
   __________         ________                     
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

 -------------------------------------------------------------------------

 -------------------------------------------------------------------------

 '''

import os, sys
import subprocess
from setuptools import setup, find_packages
import glob, shutil, re

__version__ = 1.2

try:
	from Cython.Build import cythonize
	import numpy
	optimal_list = cythonize("PyGran/analysis/core.pyx")
	include_dirs = [numpy.get_include()]
except:
	optimal_list = []
	include_dirs = []

from distutils.command.install import install


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

def link(dir):
	txt="""[Desktop Entry]
Version=1.0
Type=Application
Name=PyGran
Comment=graphical User Interface for PyGran
Exec=python -m PyGran
Icon={}/Gran.png
Terminal=false
MimeType=image/x-foo
NotShowIn=KDE""".format(dir)

	with open('PyGran.desktop', 'w') as fp:
		fp.write(txt)

	os.system('chmod +x PyGran.desktop')
	os.system('mv PyGran.desktop ~/.local/share/applications')

class Test(install):
    """ A class that performs basics tests on PyGran's simulator """

    def run(self):
        pass


class LIGGGHTS(install):
    """ A class that enables the compilation of LIGGGHTS-PUBLIC from github """

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

    def do_pre_install_stuff(self):
        
        if os.path.exists('LIGGGHTS-PUBLIC'):
            shutil.rmtree('LIGGGHTS-PUBLIC')

        self.execute(cmd='git clone https://github.com/CFDEMproject/LIGGGHTS-PUBLIC.git')
        
        files = glob.glob('LIGGGHTS-PUBLIC/src/*.cpp')

        count = 0
        self.execute(cmd='make clean-all', cwd='LIGGGHTS-PUBLIC/src')

        print('Compiling LIGGGHTS as a shared library\n')

        for path in self.execute(cmd='make auto', cwd='LIGGGHTS-PUBLIC/src'):
            count +=1
            self.print_progress(count, prefix = 'Progress:', suffix = 'Complete', total = len(files) * 2.05)

        self.execute(cmd='make -f Makefile.shlib auto', cwd='LIGGGHTS-PUBLIC/src')
        sys.stdout.write('\nInstallation of LIGGGHTS-PUBLIC complete\n')

    def run(self):
        self.do_pre_install_stuff()

example_files = []
for pdir in glob.glob('PyGran/demo/scripts/*'):
    for psdir in glob.glob(pdir+'/*'):
        example_files.append(psdir + '/*.py') #re.sub('/scripts', '', pyfile))

setup(
    name = "PyGran",
    version = __version__,
    author = "Andrew Abi-Mansour",
    author_email = "support@pygran.org",
    description = ("A DEM toolbox for rapid quantitative analysis of granular/powder systems"),
    license = "GNU v2",
    keywords = "Discrete Element Method, Granular Materials",
    url = "https://github.com/Andrew-AbiMansour/PyGran",
    packages=find_packages(),
    package_dir={'PyGran': 'PyGran'},
    package_data={'simulation': ['model_template.h'], 'gui': ['Icons/*.png'], 'demo': ['*']},
    include_package_data=True,
    install_requires=['numpy', 'scipy', 'vtk', 'pytool', 'cython', 'mpi4py', 'pillow'],
    long_description='A DEM toolbox for rapid quantitative analysis of granular/powder systems. See https://andrew-abimansour.github.io/PyGran.',
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Utilities",
        "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
	"Programming Language :: Python :: 2.7",
	"Programming Language :: Python :: 3.6"
    ],

    cmdclass={'build_liggghts': LIGGGHTS, 'test': Test},
    zip_safe=False,
    ext_modules=optimal_list,
    include_dirs=include_dirs
)
