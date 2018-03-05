from PyGran.Analyzer import System, Particles, Neighbors

class CoarseGrained(Particles):
	def __init__(self, **args):
		super(CoarseGrained, self).__init__(**args)

		if 'scale' in args and 'percent' in args:
			self.scale(args['scale'], ('radius',))
			NNS = Neighbors(self)

			self.__init__(CoarseGrained=NNS.filter(percent=args['percent']))

if __name__ == '__main__':
	Traj = System(CoarseGrained='traj.dump', units='micro', scale=2.0, percent=1.0)
