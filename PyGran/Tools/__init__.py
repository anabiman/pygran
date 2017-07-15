import os

# Conversion systems
# S.I. to microns
siToMicro = {'distance': 1e6, \
			'time': 1e6, \
			'mass': 1e6, 
}

microToSi = {'distance': 1e-6, \
			'time': 1e-6, \
			'mass': 1e-6, 
}

def find(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)

    return None