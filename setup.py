'''
Created on July 9, 2016
@author: Andrew Abi-Mansour
'''

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

import os, sys, numpy, site
from setuptools import setup, find_packages
from Cython.Build import cythonize
from distutils.command.install import install

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

def link(dir):
	txt="""[Desktop Entry]
Version=1.0
Type=Application
Name=PyDEM
Comment=graphical User Interface for PyDEM	     
Exec=python -m PyDEM
Icon={}/PyDEM.png
Terminal=false
MimeType=image/x-foo
NotShowIn=KDE""".format(dir)

	with open('PyDEM.desktop', 'w') as fp:
		fp.write(txt)

	os.system('chmod +x PyDEM.desktop')
	os.system('mv PyDEM.desktop ~/.local/share/applications')

class LIGGGHTS(install):
    def do_pre_install_stuff(self):
        os.chdir('engines/liggghts/src')
        os.system('make makeshlib')
        os.sytem('make -f Makefile.shlib openmpi')

    def run(self):
        self.do_pre_install_stuff()
        install.run(self)
        self.do_post_install_stuff()

setup(
    name = "PyDEM",
    version = "0.0.1",
    author = "Andrew Abi-Mansour",
    author_email = "andrew.gaam@gmail.com",
    description = ("A set of tools for running, analyzing, and visualizing DEM simulations"),
    license = "GNU v3",
    keywords = "Discrete Element Method, Granular Materials",
    url = "https://github.com/Andrew-AbiMansour/PyDEM",
    packages=find_packages(),
    package_dir={'PyDEM': 'PyDEM'},
    package_data={'': ['.config', 'GUI/Icons/*.png']},
    install_requires=['numpy', 'pyvtk', 'pytools>=2011.2', 'pygmsh'],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 1 - Alpha",
        "Topic :: Utilities",
        "License :: GNU License",
    ],
    cmdclass={'build_liggghts': LIGGGHTS},
    zip_safe=False,
    ext_modules=cythonize("PyDEM/Analyzer/*.pyx"),
	include_dirs=[numpy.get_include()]
)

if sys.argv[1] == 'install':

	sys.path.remove(os.getcwd()) # remove current path to make sure PyDEM is imported from elsewhere
	sys.path.append('/home/abimanso/.local/lib/python2.7/site-packages')
	print sys.path
	print 'Verifying installation ....................................................'
	print '...........................................................................'
	print '...........................................................................'

	try:
		import PyDEM
		link(PyDEM.__path__[0] + '/GUI/Icons')
		print 'PyDEM successfully installed'
	except:
		print 'PyDEM installation failed ...'
		raise
