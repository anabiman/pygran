import numpy as np

class C(np.ndarray):

    def __init__(self, obj):
	if type(obj) != np.ndarray:
		raise TypeError('Error in input. Only ndarray can be used to instantiate this class')
	else:
		print('pass')
