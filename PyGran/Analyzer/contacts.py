import scipy, numpy

class Contacts(object):
	""" A dynamic class that contains all the particle-particle (and optionally particle-wall)
	contacts
	"""
	def __init__(self, Particles):

		coords = numpy.array([Particles.x, Particles.y, Particles.z]).T
		tree = scipy.spatial.cKDTree(coords, leafsize=100)
		maxDist = 2.0 * Particles.radius.max()

		self.pairs = tree.query_pairs(maxDist)
		self.overlaps = numpy.zeros(len(self.pairs))
		count = 0

		for pair in self.pairs:
			self.overlaps[count] = numpy.sqrt((Particles.x[pair[0]] - Particles.x[pair[1]])**2.0 + \
									(Particles.y[pair[0]] - Particles.y[pair[1]])**2.0 + \
									(Particles.z[pair[0]] - Particles.z[pair[1]])**2.0)

			count += 1 
