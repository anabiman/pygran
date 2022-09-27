import matplotlib as mpl
import matplotlib.pylab as plt
from numpy import arange, array, fabs, sqrt
from pygran import simulation as sim
from pygran.params import cohesionless


def run():

    # Setup matplotlib params
    mpl.rc("text", usetex=True)
    plt.rcParams.update({"font.size": 18})

    cModel = sim.models.ThorntonNing

    COR, yieldVel = [], []

    # Create a yield pressure array to study
    pressure = arange(1, 6, 0.1) * cohesionless["youngsModulus"] * 0.01

    # Set particle radius to 100 microns
    cohesionless["radius"] = 1e-4

    for yieldPress in pressure:

        cohesionless["yieldPress"] = yieldPress
        model = cModel(material=cohesionless)

        time, disp, force = model.displacement()

        deltav = disp[:, 1]

        COR.append(fabs(deltav[-1] / cohesionless["characteristicVelocity"]))
        yieldVel.append(model.yieldVel / cohesionless["characteristicVelocity"])

    ratio = array(yieldVel)
    ratio[ratio > 1] = 1

    COR_anal = (
        sqrt(6.0 * sqrt(3.0) / 5.0)
        * sqrt(1.0 - (1.0 / 6.0) * (ratio) ** 2)
        * (ratio / (ratio + 2 * sqrt(6.0 / 5.0 - (1.0 / 5.0) * ratio**2))) ** (0.25)
    )

    fig = plt.figure()

    plt.plot(pressure, COR, "-o")
    plt.plot(pressure, COR_anal)

    plt.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    plt.xlabel("Yield Pressure (MPa)")
    plt.ylabel("Coefficient of Restitution")
    plt.grid(linestyle=":")
    plt.legend(["Numerical", "Analytical"])
    plt.show()


if __name__ == "__main__":

    run()
