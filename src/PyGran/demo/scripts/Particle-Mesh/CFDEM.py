from PyGran import analysis
from scipy.linalg import norm
from numpy import array
import matplotlib.pylab as plt

# Create a mesh trajectory file for the inlet & outlet files
inlet, outlet, inter = 'CFD/inlet_*.vtk', 'CFD/outlet_*.vtk', 'CFD/inter_*.vtk'
inlet_args = {'vtk_type': 'poly'}
outlet_args = {'vtk_type': 'poly'}
inter_args = {'avgCellData': False}

Meshes = [(inlet, inlet_args), (outlet, outlet_args), (inter, inter_args)]

Bed = analysis.System(Particles='DEM/particles_*.dump', Mesh=Meshes)

# Compute bed height
bed_height = Bed.Particles.z.max()

print(Bed.__dict__)

# Loop over inlet trajectory and compute the pressure drop
for ts in Bed:
	print(ts)

press_drop = [iMesh.p - oMesh.p for ts in Bed]

# Go back to frame 0 ~ do we need this?
Bed.rewind()

# Loop over inlet trajectory and compute the inlet velocity
ivel = [norm(Bed.Mesh[0].U) for ts in Bed]

Bed.rewind()
voidFraction = array([Bed.Mesh[2].voidfraction])

# Compute pressure drop based on Ergun equation
density, viscous, diameter = 10.0, 1.5e-4, 1e-3

term1 = 150.0 * viscous / diameter**2 * ivel * (1 - voidFraction)**2 / voidFraction**3
term2 = 1.75 * density / diameter * ivel * (1 - voidFraction) / voidFraction**3 * ivel

ergun_dp = bed_height * 2 * (term1 + term2)

# Plot pressure drop vs inlet velocity

plt.plot(ivel, ergun_dp * 1e-3, '--')
plt.plot(ivel, press_drop, '*')

plt.xlabel('Inlet Velocity (m/s)')
plt.ylabel('Pressure drop (Pa)')
plt.grid()
plt.show()
