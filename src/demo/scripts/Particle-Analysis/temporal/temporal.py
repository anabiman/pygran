from PyGran import analysis
import matplotlib.pylab as plt
from numpy import array

# Create a PyGran System from a dump file (default si units)
System = analysis.System(Particles='traj*.dump')

# Timestep used in simulation (s)
dt = 1e-6

# Skip empty frames
System.skip()

# Extract bed height + timesteps
data = array([[ts * dt,System.Particles.z.max()] for ts in System])

# Plot bed height (mm) vs time (ms)
plt.plot(data[:,0] * 1e3, data[:,1] * 1e3, '-o')

plt.xlabel('Time (ms)')
plt.ylabel('Height (mm)')
plt.grid(linestyle=':')
plt.show()
