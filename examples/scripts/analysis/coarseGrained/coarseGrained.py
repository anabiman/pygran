import os

from pygran import analysis


class CoarseParticles(analysis.Particles):
    def __init__(self, **args):
        super().__init__(**args)

        if "scale" in args and "percent" in args:
            self.scale(args["scale"], ("radius",))
            CG = analysis.equilibrium.Neighbors(self).filter(percent=args["percent"])

            self.__init__(CoarseParticles=CG)


if __name__ == "__main__":
    Traj = analysis.System(
        CoarseParticles="traj.dump", scale=3, percent=10.0, module="coarseGrained"
    )
    Traj.CoarseParticles.write("CG.dump")
