from pygran import analysis


def test_foo(trajf):
    # Create a granular object from a LIGGGHTS dump file
    sys = analysis.System(Particles=trajf, units="si")

    # Go to last frame
    sys.goto(-1)

    # Switch to micro unit system
    sys.units("micro")

    # Compute the radial distribution function
    g, r, _ = sys.Particles.computeRDF(optimize=False)

    assert r.max() > 4.0

    # Construct a class for nearest neighbor searching
    neigh = analysis.equilibrium.Neighbors(sys.Particles)

    # Extract coordination number per particle
    coon = neigh.coon
