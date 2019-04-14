from PyGran import analysis
from scipy.linalg import norm
from numpy import array
import matplotlib.pylab as plt

# Create a mesh trajectory file for the inlet & outlet files
inlet, outlet, inter = 'CFD/inlet_*.vtk', 'CFD/outlet_*.vtk', 'CFD/inter_*.vtk'
inlet_args = {'vtk_type': 'poly', 'avgCellData': True} #, 'skip': 4}
outlet_args = {'vtk_type': 'poly', 'avgCellData': True} #, 'skip': 4}
inner_args = {} #{'skip': 4}

Meshes = (inlet, inlet_args), (outlet, outlet_args), (inter, inner_args)

#Bed = analysis.System(Particles='DEM/particles_*.dump', 
Bed = analysis.System(Mesh=Meshes)

Inlet, Outlet, Inner = Bed.Mesh

# Compute bed height
Parts = analysis.System(Particles='DEM/particles_*.dump')
bed_height = Parts.Particles.z.max() # Outlet.Points[:,-1].max() - Inlet.Points[:,-1].min()

# Loop over inlet trajectory and compute the pressure drop
press_drop = [Inlet.cells.p - Outlet.cells.p for ts in Bed]

# Loop over inlet trajectory and compute the inlet velocity
ivel = [norm(Inlet.cells.U) for ts in Bed]

voidFraction = [Inner.points.voidfraction[Inner.points.voidfraction < 1].mean() for ts in Bed]

# Compute pressure drop based on Ergun equation
density, viscous, diameter = 10.0, 1.5e-4, Parts.Particles.radius.mean() * 2

# Convert lists to arrays
press_drop = array(press_drop)
ivel = array(ivel)
voidFraction = array(voidFraction)

term1 = 150.0 * viscous / diameter**2 * ivel * (1 - voidFraction)**2 / voidFraction**3
term2 = 1.75 * density / diameter * ivel * (1 - voidFraction) / voidFraction**3 * ivel

ergun_dp = 1.8 * bed_height * (term1 + term2)

# Plot pressure drop vs inlet velocity
plt.plot(ivel * 1e3, press_drop, '*')
plt.plot(ivel * 1e3, ergun_dp, '--')

plt.xlabel('Inlet Velocity (mm/s)')
plt.ylabel('Pressure drop (Pa)')
plt.legend(['Sim', 'Anal'])
plt.grid()
plt.show()
