# Import the analysis module from PyGran
from PyGran import analysis
import sys
import pytest


def test_foo(trajf):
    # Create a System object from a LIGGGHTS dump file. They keyword 'Particles' is mandatory, since it instructs
    # System to create an object of type 'Particles' which can be assessed from the instantiated System object.

    Sys = analysis.System(Particles=trajf)

    # Go to last frame
    Sys.goto(-1)
    assert Sys.frame == 30

    # Create a reference to Sys.Particles
    Particles = Sys.Particles

    # Any changes made to Particles is reflected in Sys.Particles. To avoid that, a hard copy of Particles should be
    # created instead:
    Particles = Sys.Particles.copy()

    len_parts = len(Particles)
    assert len_parts > 0

    nparts = 0
    # Looping over Particles yields a Particles class of length 1
    for part in Particles:
        nparts += len(part)

    assert nparts == len_parts

    # Slice Particles into a new class containing the 1st 10 particles
    Slice = Particles[:10]
    assert len(Slice) == 10

    # Slice Particles into a new class containing particles 1, 2, and 10
    Slice = Particles[[1, 2, 10]]
    assert len(Slice) == 3

    # Slice Particles into a new class containing the 1st 10 particles and the last 10 particles
    Slice = Particles[:10] + Particles[-10:]

    # More sophisticated slicing can be done with 1 or boolean expressions. For example:

    # Create a Particles class containing particles smaller than 25% of the mean particle diameter
    SmallParts = Particles[Particles.radius <= Particles.radius.mean() * 0.25]
    assert len(SmallParts) == 0

    # Create a Particles class containing particles in the positive xy plane
    SmallParts = Particles[(Particles.y >= 0) & (Particles.x >= 0)]
    assert len(SmallParts) > 0
