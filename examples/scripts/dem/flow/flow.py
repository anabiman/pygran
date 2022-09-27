"""
Created on April 22, 2017
@author: Andrew Abi-Mansour
"""

from pygran import simulation
from pygran.params import glass, organic

params = {
    # Define the system
    "boundary": ("f", "f", "f"),
    "box": (-1e-3, 1e-3, -1e-3, 1e-3, 0, 4e-3),
    # Define component(s)
    "species": ({"material": organic, "radius": ("constant", 5e-5)},),
    # Setup I/O params
    "traj": {"freq": 1000, "pfile": "particles*.dump", "mfile": "mesh*.vtk"},
    # Output dir name
    "output": "DEM_flow",
    # Define computational parameters
    "dt": 1e-6,
    # Apply a gravitional force in the negative direction along the z-axis
    "gravity": (9.81, 0, 0, -1),
    # Import hopper mesh
    "mesh": {
        "hopper": {
            "file": "mesh/silo.stl",
            "mtype": "mesh/surface",
            "material": glass,
            "args": {"scale": 1e-3},
        },
        "impeller": {
            "file": "mesh/valve.stl",
            "mtype": "mesh/surface",
            "material": glass,
            "args": {"move": (0, 0, 1), "scale": 1e-3},
        },
    },
    # Stage runs
    "stages": {"insertion": 1e5, "run": 1e5},
}


def run(**params):

    # Create an instance of the DEM class
    sim = simulation.DEM(**params)

    # Setup a stopper wall along the xoy plane
    stopper = sim.setupWall(species=1, wtype="primitive", plane="zplane", peq=0.0)

    # Monitor the time average ensemble kinetic energy
    ke = sim.monitor(
        var="ke",
        species="all",
        name="ensemble_ke",
        file="ke.dat",
        nevery=100,
        nfreq=1000,
        nrepeat=10,
    )

    # Insert particles in a cubic region
    insert = sim.insert(
        species=1,
        region=("block", -5e-4, 5e-4, -5e-4, 5e-4, 2e-3, 3e-3),
        mech="volumefraction_region",
        value=1,
        freq=1e4,
    )
    sim.run(params["stages"]["insertion"], params["dt"])
    sim.remove(insert)

    # Rotate the impeller
    rotImp = sim.moveMesh(
        name="impeller", rotate=("origin", 0, 0, 0), axis=(0, 0, 1), period=5e-2
    )

    # Blend the system
    sim.run(params["stages"]["run"], params["dt"])
    sim.remove(rotImp)

    # Remove stopper and monitor flow
    sim.remove(stopper)
    sim.run(params["stages"]["run"], params["dt"])


if __name__ == "__main__":
    run(**params)
