import os
from setuptools import setup, find_packages

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "PyDEM",
    version = "0.0.1",
    author = "Andrew Abi-Mansour",
    author_email = "andrew.abi.mansour@merck.com",
    description = ("A set of tools for running and analyzing DEM simulations"),
    license = "GNU",
    keywords = "Discrete Element Method, Granular Materials",
    url = "https://github.com/Andrew-AbiMansour/PyDEM",
    packages=find_packages(),
    package_dir={'PyDEM': 'PyDEM'},
    install_requires=['numpy', 'pyvtk', 'pytools>=2011.2'],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 1 - Alpha",
        "Topic :: Utilities",
        "License :: GNU License",
    ],
    zip_safe=True,
)