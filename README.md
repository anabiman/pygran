# Welcome to the PyGran webpage!

PyGran is a toolkit for analyzing DEM simulation data. In addition to performing routine and custom post-processing, PyGran enables running DEM simulation with LIGGGHTS in Python.

Installing PyGran is quite straight forward on a Uni/Unix-like machine. Just fire up a terminal and run from the PyGran directory:

```python
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

![alt text](images/overlap-hist.pdf "Logo Title Text 1")
