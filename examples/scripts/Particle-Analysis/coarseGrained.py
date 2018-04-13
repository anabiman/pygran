from PyGran import Analyzer

class CoarseParticles(Analyzer.Particles):
	def __init__(self, **args):
		super(CoarseParticles, self).__init__(**args)

		if 'scale' in args and 'percent' in args:
			self.scale(args['scale'], ('radius',))
			CG = Analyzer.equilibrium.Neighbors(self).filter(percent=args['percent'])

			self.__init__(CoarseParticles=CG)

if __name__ == '__main__':
	Traj = Analyzer.System(CoarseParticles='traj.dump', units='micro', scale=3, percent=25.0)
	Traj.CoarseParticles.write('CG.dump')
