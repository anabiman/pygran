"""
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
"""

import os, sys, shutil
import subprocess
from setuptools import setup, find_packages
import shutil, pathlib
from distutils.command.install import install
from distutils.command.clean import clean
import versioneer
import fileinput


short_description = __doc__.split("\n")

try:
    with open("README.md", "r") as handle:
        long_description = handle.read()
except Exception:
    long_description = "\n".join(short_description[2:])

try:
    from Cython.Build import cythonize
    import numpy

    optimal_list = cythonize(
        "pygran/analysis/pygran_analysis/core.pyx",
        compiler_directives={"language_level": sys.version_info[0]},
    )
    include_dirs = [numpy.get_include()]
except Exception:
    print("Could not cythonize. Make sure Cython is properly installed.")
    optimal_list = []
    include_dirs = []


class Track(install):
    """An install class that enables the tracking of installation/compilation progress"""

    def execute(self, cmd, cwd="."):
        popen = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, universal_newlines=True, cwd=cwd, shell=True
        )

        for stdout_line in iter(popen.stdout.readline, ""):
            yield stdout_line

        popen.stdout.close()
        return_code = popen.wait()

        if return_code:
            raise subprocess.CalledProcessError(return_code, cmd)

    def print_progress(self, iteration, total, prefix="", suffix="", decimals=1):
        """
        Call in a loop to create terminal progress bar
        @params:
                iteration   - Required  : current iteration (Int)
                total       - Required  : total iterations (Int)
                prefix      - Optional  : prefix string (Str)
                suffix      - Optional  : suffix string (Str)
                decimals    - Optional  : positive number of decimals in percent complete (Int)
        """
        str_format = "{0:." + str(decimals) + "f}"
        percents = str_format.format(100 * iteration / total)
        sys.stdout.write(f"\r{prefix} {percents}% {suffix}")
        sys.stdout.flush()

    def run(self):
        self.do_pre_install_stuff()

    def do_pre_install_stuff(self):
        raise NotImplementedError


class LIGGGHTS(Track):
    """A class that enables the compilation of LIGGGHTS-PUBLIC from github"""

    def fix_mkfile(self, filename: str, search_txt: str, insert_txt: str):
        with fileinput.input(files=(filename,), inplace=True, backup=".bak") as fileobj:
            for line in fileobj:
                if search_txt in line:
                    line = line.strip() + " " + insert_txt + "\n"
                print(line, end="")

    def do_pre_install_stuff(self):

        if os.path.exists("LIGGGHTS-PUBLIC"):
            print("Deleting " + "LIGGGHTS-PUBLIC")
            shutil.rmtree("LIGGGHTS-PUBLIC")

        self.spawn(
            cmd=["git", "clone", "https://github.com/CFDEMproject/LIGGGHTS-PUBLIC.git"]
        )

        path = pathlib.Path("LIGGGHTS-PUBLIC", "src")
        files = list(path.glob("*.cpp"))

        count = 0
        os.chdir(path)
        self.spawn(cmd=["make", "clean-all"])

        print("Configuring makefile ...")
        mkfile = pathlib.Path("MAKE", "Makefile.auto")
        self.fix_mkfile(mkfile, "-funroll-loops", "-w")

        print("Compiling LIGGGHTS as a executable ...\n")
        for path in self.execute(cmd="make auto"):
            count += 1
            self.print_progress(
                count,
                total=2*len(files)*1.0215, # makefile prints 2 lines per file
                prefix="Progress:",
                suffix="complete",
            )

        print("\nGenerating shared object file ...\n")
        self.spawn(cmd=["make", "-f", "Makefile.shlib", "auto"])
        sys.stdout.write("\nInstallation of LIGGGHTS-PUBLIC complete\n")
        os.chdir(os.path.join("..", ".."))


cmdclass = versioneer.get_cmdclass()
cmdclass["build_liggghts"] = LIGGGHTS

setup(
    name="pygran",
    author="Andrew Abi-Mansour",
    author_email="support@pygran.org",
    description=(
        "A DEM toolkit for rapid quantitative analysis and simulation of granular/powder systems"
    ),
    license="GNU v2",
    keywords="Discrete Element Method, Granular Materials",
    url="https://github.com/Andrew-AbiMansour/PyGran",
    include_package_data=True,
    install_requires=["numpy", "scipy"],
    packages=find_packages(),
    extras_require={
        "extra": ["vtk", "mpi4py", "Pillow", "pytest", "pytest-cov", "codecov"]
    },
    long_description="A DEM toolbox for rapid quantitative analysis and simulation of granular/powder systems. See http://www.pygran.org.",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Topic :: Utilities",
        "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: C",
        "Operating System :: POSIX :: Linux",
    ],
    version=versioneer.get_version(),
    cmdclass=cmdclass,
    zip_safe=False,
    ext_modules=optimal_list,
    include_dirs=include_dirs,
)
