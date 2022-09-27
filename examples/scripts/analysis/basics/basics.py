# Basics
# ======
# In this tutorial, a general overview of the Particles object in pygran is presented. Class instantiation, manipulation (slicing, indexing, etc.), and basic propreties are covered.

import sys

# Import the analysis module from pygran
from pygran import analysis

# Read particle trajectory input filename
pfname = sys.argv[1]

# Create a System object from a LIGGGHTS dump file. They keyword 'Particles' is mandatory, since it instructs
# System to create an object of type 'Particles' which can be assessed from the instantiated System object.
sys = analysis.System(Particles=pfname)

# go to last frame
sys.goto(-1)

# Particles Class
# ===============
# The code below shows how a Particles object and its dynamic attributes can be accessed.

# Create a reference to Sys.Particles
particles = sys.Particles

# Any changes made to Particles is reflected in Sys.Particles. to avoid that, a hard copy of Particles should be
# created instead:
particles = sys.Particles.copy()

# The length of Particles is the number of particles contained in this class
print("Number of particles stored = ", len(particles))

nparts = 0
# Looping over Particles yields a Particles class of length 1
for part in particles:
    nparts += len(part)

print("Number of particles calculated = ", nparts)


# Print all attributes stored in Particles
print(particles.keys)

# Particle slicing
# =================
# Particles can be sliced using a numpy-like syntax as shown below. A sliced particles orject is always created as a hard copy of the target object.


# Slice particles into a new class containing the 1st 10 particles
slice = particles[:10]

print("number of particles in 1st slice is", len(slice))

# Slice particles into a new class containing particles 1, 2, and 10
slice = particles[[1, 2, 10]]

print("number of particles in 2nd slice is", len(slice))

# Slice particles into a new class containing the 1st 10 particles and the last 10 particles
slice = particles[:10] + particles[-10:]

print("number of particles in 3rd slice is", len(slice))

# More sophisticated slicing can be done with 1 or boolean expressions. For example:

# Create a Particles object containing particles smaller than 25% of the mean particle diameter
small_parts = particles[particles.radius <= particles.radius.mean() * 0.25]

# Create a Particles object containing particles below the mean height (z-direction) and in the positive portion
# of the domain along the x-axis
small_parts = particles[(particles.z <= particles.z.mean()) & (particles.x >= 0)]

# Modifying and Writing Particles
# ===============================

# 'Particles' provides many methods that change its state. For example, all particles are translated along
# the z-direction by -z.min() with the following command:
particles.translate(("z",), (-particles.z.min(),))

# Particles can be written as a new dump file, which by default is appended to an existing file.
particles.write("translated_system.dump")

print("Test completed successfully")
