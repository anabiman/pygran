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
