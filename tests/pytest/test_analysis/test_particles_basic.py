# Import the analysis module from PyGran
from PyGran import analysis
import sys

traj_name = sys.argv[1]

# Create a System object from a LIGGGHTS dump file. They keyword 'Particles' is mandatory, since it instructs 
# System to create an object of type 'Particles' which can be assessed from the instantiated System object.
Sys = analysis.System(Particles=traj_name)

# Create a reference to Sys.Particles
Particles = Sys.Particles

# Any changes made to Particles is reflected in Sys.Particles. to avoid that, a hard copy of Particles should be 
# created instead:
Particles = Sys.Particles.copy()

nparts = 0
# Looping over Particles yields a Particles class of length 1
for part in Particles:
    nparts += len(part)

# Slice Particles into a new class containing the 1st 10 particles
Slice = Particles[:10]

# Slice Particles into a new class containing particles 1, 2, and 10
Slice = Particles[[1,2,10]]

# Slice Particles into a new class containing the 1st 10 particles and the last 10 particles
Slice = Particles[:10] + Particles[-10:]

# More sophisticated slicing can be done with 1 or boolean expressions. For example:

# Create a Particles class containing particles smaller than 25% of the mean particle diameter
SmallParts = Particles[Particles.radius <= Particles.radius.mean() * 0.25]

# Create a Particles class containing particles below the mean height (z-direction) and in the positive portion
# of the domain along the x-axis
SmallParts = Particles[(Particles.z <= Particles.z.mean()) & (Particles.x >= 0)]
