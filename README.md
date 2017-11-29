# Welcome to the PyGran webpage!

PyGran is a toolkit primarily designed for analyzing DEM simulation data. In addition to performing basic and custom post-processing, PyGran enables running DEM simulation with [LIGGGHTS](https://www.cfdem.com/liggghtsr-open-source-discrete-element-method-particle-simulation-code) in Python.

## Quick Installation
Installing PyGran is quite straight forward on a Unix/Unix-like machine. Just fire up a terminal and run from the PyGran source directory:
```bash
python setup.py install --user
```
For more options and information on setting up PyGran, see chapter I in the manual.
## Basic Usage
### Contact Analysis
Using PyGran is also quite straight forward. Computing particle overlaps for instance can be done in few lines of code:

```python
from PyGran import Analyzer
from PyGran.Materials import glass

Gran = Analyzer.System(Particles='granular.dump')
NNS = Analyzer.Neighbors(Gran.Particles, material=glass)
overlaps = NNS.overlaps
```
<p style="text-align:center;"><img src="images/overlap-hist.png" width="600"></p>
### Simulating Granular Flow
PyGran also provides an interface for running DEM simulation with [LIGGGHTS](https://www.cfdem.com/liggghtsr-open-source-discrete-element-method-particle-simulation-code). A sample script that simulates flow in a hopper is shown below.

```python
from PyGran import Simulator
from PyGran.Materials import stearicAcid, steel

pDict = {

	'model': Simulator.models.SpringDashpot,
	'boundary': ('f','f','f'),
	'box':  (-1e-3, 1e-3, -1e-3, 1e-3, 0, 4e-3),

	'SS': ({'insert': 'by_pack', 'material': stearicAcid,'natoms': 1000, \
		'freq': 'once', 'radius': ('constant', 5e-5),}, 
		),
		
	'dt': 1e-6,
	'gravity': (9.81, 0, 0, -1),

	'mesh': {
		'hopper': {'file': 'silo.stl', 'mtype': 'mesh/surface', 'material': steel},
		},
	'stages': {'insertion': 1e4},
}

pDict['model'] = pDict['model'](**pDict)
sim = Simulator.DEM(**pDict['model'].params)
insert = sim.insert('cubic', 1, *('block', pDict['box'])
sim.run(pDict['stages']['insertion'], pDict['dt'])
```
<p style="text-align:center;"><img src="images/hopper.png" width="600"></p>
For more examples on using PyGran for running DEM simulation, see chapter II in the manual.
