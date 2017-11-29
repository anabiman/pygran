# Welcome to the PyGran webpage!

PyGran is a toolkit primarily designed for analyzing DEM simulation data. In addition to performing basic and advanced post-processing, PyGran enables running DEM simulation with [LIGGGHTS](https://www.cfdem.com/liggghtsr-open-source-discrete-element-method-particle-simulation-code) in Python.

Installing PyGran is quite straight forward on a Unix/Unix-like machine. Just fire up a terminal and run from the PyGran source directory:

```bash
python setup.py install --user
```

Using PyGran is also quite straight forward:

```python
from PyGran import Analyzer
from PyGran.Materials import glass

Gran = Analyzer.System(Particles='granular.dump')
NNS = Analyzer.Neighbors(Gran.Particles, material=glass)
overlaps = NNS.overlaps
```

<center><img src="images/overlap-hist.png" width="600" alt="Smiley face"></center>
