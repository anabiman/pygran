import PyGran.Simulator as Si
from PyGran.Materials import stearicAcid
from numpy import random

SD = Si.models.SpringDashpot(material=stearicAcid)
radius = 1e-4

time, delta = SD.displacement(radius)
