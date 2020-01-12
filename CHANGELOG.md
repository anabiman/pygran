# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and as of v1.2 this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0]
### Added

- New dynamic contact model module for the LIGGGHTS-PUBLIC engine. Currently this module supports only contact models derived from the spring-dashpot model.
- Support for dynamic keyword parsing in the LIGGGHTS-PUBLIC engine.
- New tutorial for installing PyGran/LIGGGHTS on Ubuntu 18.04 LTS.

### Changed
- Renamed Simulator, Analyzer, and Material modules to simulation, analysis, and params, respectively.
- Fixed bug in the contact duraction computation of the numerical solver in the simulation module.

## [1.2.0]
### Added
- New `configure' function for specifying which LIGGGHTS version to use
- New `increment' option for System frame propagation. Certain frames can now be skipped when reading multiple traj files.
- Cell and point attributes are exposed separately via: Mesh.points.attr and Mesh.cells.attr
- New version navigation
- Online documentation based on Sphinx
- Pytest scripts for validation

### Removed
- Manual in pdf format

### Changed
- SubSystem objects are now constructed only once (factory) and thus retain the same memory address from frame to frame
- Mesh SubSystem now takes the optional arg `avgCellData' which transforms cell arrays to their weighted-avg values
- System can now be constructed with additional arguments passed for each SubSystem type
- System rewinds itself after the last frame is reached and it starts reading frame 1 when it is looped over.
- Fixed bug with System frame propagation when both Particles and Mesh were stored
- Fixed bug in SSMP mode when the slave processors did not produce any output 
- Fixed bug with user-created scripts importing themselves (leading to double execution)
- Renamed method `add_viscous' to `addViscous' in engine_liggghts and dem modules
- DEM engine config files are stored in ~/.config/engine_name.ini
- System method 'units' always applies to the currently loaded frame
- System method 'units' always returns current unit system (str) 

## [1.2.1]
### Added
- A random prime number generator for LIGGGHTS (in engine_ligggghts.py)
- Support for polydisperse systems (normal psd, lognormal psd, and discrete distribution)
- PyGran.simulation version printed to pygran.log file

### Changed
- Fixed bug in Particles._goto() method

## [1.2.2]
### Added
- New interface for user-defined functions used in multisphere particles
- version, author, and email available from __init__ module
- New polydisperse example for DEM simulation

### Changed
- Fixed bugs in intensitySegergation routines in analysis.core.py*
- Improved docstrings in analysis.core.py* and analysis.equilibrium.py
- Renamed Particles.molecules to Particles.Molecules
- Fixed typo in compaction.html file (missing quotes for wallZ)

### Removed
- Method computeGCOM() in analysis.core.py*
- Attribute 'length' for Particles.molecule
