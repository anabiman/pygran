from PyGran import analysis
import sys


def test_foo(trajf):
    # Create a granular object from a LIGGGHTS dump file
    Sys = analysis.System(Particles=trajf, units="si")

    # Go to last frame
    Sys.goto(-1)

    # Switch to micro unit system
    Sys.units("micro")

    # Compute the radial distribution function
    g, r, _ = Sys.Particles.computeRDF()

    assert r.max() > 4.0

    # Construct a class for nearest neighbor searching
    Neigh = analysis.equilibrium.Neighbors(Sys.Particles)

    # Extract coordination number per particle
    coon = Neigh.coon
