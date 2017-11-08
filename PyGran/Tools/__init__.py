import os

# Conversion unit systems
# S.I.: m, s, Kg
# CGS: cm, s, g
# Micro: micron, micro s, micro g

conversion = {'si': {'distance': 1, 'time': 1, 'mass': 1}, \
			'cgs': {'distance': 1e-2,'time': 1, 'mass': 1e-3}, \
			'micro': {'distance': 1e-6, 'time': 1e-6, 'mass': 1e-9}, \
			'nano': {'distance': 1e-9, 'time': 1e-9, 'mass': 1e-12}
			}

def convert(unitso, unitsf):

	if unitso in conversion:
		if unitsf in conversion:
			conv = conversion[unitso]

			for key in conv:
				conv[key] /= conversion[unitsf][key]

			return conv
		else:
			raise ValueError('Input unit system not supported: {}'.format(unitsf))
	else:
		raise ValueError('Input unit system not supported: {}'.format(unitso))

def find(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)

    return None