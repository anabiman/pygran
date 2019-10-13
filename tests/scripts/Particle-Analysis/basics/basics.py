# Basics
# ======
# In this tutorial, a general overview of the Particles object in PyGran is presented. Class instantiation, manipulation (slicing, indexing, etc.), and basic propreties are covered.

# Import the analysis module from PyGran
from PyGran import analysis
import sys

# Read particle trajectory input filename
pfname = sys.argv[1]

# Create a System object from a LIGGGHTS dump file. They keyword 'Particles' is mandatory, since it instructs 
# System to create an object of type 'Particles' which can be assessed from the instantiated System object.
Sys = analysis.System(Particles=pfname)

# go to last frame
Sys.goto(-1)

# Particles Class
# ===============
# The code below shows how a Particles class and its dynamic attributes can be accessed.

# Create a reference to Sys.Particles
Particles = Sys.Particles

# Any changes made to Particles is reflected in Sys.Particles. to avoid that, a hard copy of Particles should be 
# created instead:
Particles = Sys.Particles.copy()

# The length of Particles is the number of particles contained in this class
print("Number of particles stored = ", len(Particles))

nparts = 0
# Looping over Particles yields a Particles class of length 1
for part in Particles:
    nparts += len(part)

print("Number of particles calculated = ", nparts)


# Print all attributes stored in Particles
print(Particles.keys)

# Particle Slicing
# =================
# Particles can be sliced using a numpy-like syntax as shown below. A sliced Particles is always created as a hard copy of the target object.


# Slice Particles into a new class containing the 1st 10 particles
Slice = Particles[:10]

print("number of particles in 1st slice is", len(Slice))

# Slice Particles into a new class containing particles 1, 2, and 10
Slice = Particles[[1,2,10]]

print("number of particles in 2nd slice is", len(Slice))

# Slice Particles into a new class containing the 1st 10 particles and the last 10 particles
Slice = Particles[:10] + Particles[-10:]

print("number of particles in 3rd slice is", len(Slice))

# More sophisticated slicing can be done with 1 or boolean expressions. For example:

# Create a Particles class containing particles smaller than 25% of the mean particle diameter
SmallParts = Particles[Particles.radius <= Particles.radius.mean() * 0.25]

# Create a Particles class containing particles below the mean height (z-direction) and in the positive portion
# of the domain along the x-axis
SmallParts = Particles[(Particles.z <= Particles.z.mean()) & (Particles.x >= 0)]

# Modifying and Writing Particles
# ===============================

# 'Particles' provides many methods that change its state. For example, all particles are translated along 
# the z-direction by -z.min() with the following command:
Particles.translate(('z',), (-Particles.z.min(),))

# Particles can be written as a new dump file, which by default is appended to an existing file.
Particles.write('translated_system.dump')

print("Test completed successfully")
