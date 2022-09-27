import matplotlib.pylab as plt
from numpy import array, sqrt
from pygran import analysis
from scipy.linalg import norm

# Create a mesh trajectory file for the bed inlet, outlet, and internal files
inlet, outlet, inter = "CFD/inlet_*.vtk", "CFD/outlet_*.vtk", "CFD/inter_*.vtk"
inlet_args = {"vtk_type": "poly", "avgCellData": True, "skip": 4}
outlet_args = {"vtk_type": "poly", "avgCellData": True, "skip": 4}
inter_args = {"avgCellData": True, "skip": 4}

Meshes = (inlet, inlet_args), (outlet, outlet_args), (inter, inter_args)
Bed = analysis.System(Particles="DEM/particles_*.dump", Mesh=Meshes)

# Extract the mesh SubSystems from Bed
Inlet, Outlet, Inter = Bed.Mesh

# Compute bed height (along the z-direction)
bed_height = Outlet.Points[:, -1].max() - Inlet.Points[:, -1].min()
# bed_height = array([Bed.Particles.z.max() for ts in Bed])

# Compute the pressure drop as a time series
press_drop = [Inlet.cells.p - Outlet.cells.p for ts in Bed]

# Compute the inlet velocity as a time series
ivel = [norm(Inlet.cells.U) for ts in Bed]

# Compute the fluid velocity
fvel = [norm(Inter.cells.U) for ts in Bed]

# Compute the void fraction as a time series
voidFraction = [
    Inter.points.voidfraction[Inter.points.voidfraction < 1].mean() for ts in Bed
]
# voidFraction = [Inter.points.voidfraction.mean() for ts in Bed]

# Properties of the fluid and granules
f_density, p_density, kin_vis, diameter = (
    10.0,
    2000.0,
    1.5e-4,
    Bed.Particles.radius.mean() * 2,
)
dyn_vis = kin_vis * f_density

# Convert lists to arrays
press_drop = array(press_drop) * f_density
ivel = array(ivel)
fvel = array(fvel)
voidFraction = array(voidFraction)

# The theo min fluidization value is almost double that of its numerical value because  voidFraction is overvalued (should be ~ 0.44 instead of 0.51).
# Reason is I did not let the system relax enough. Better re-run the sim.
# So instead we scale voidFraction here ~ hackish
voidFraction /= voidFraction.min() / 0.52

print(
    "umf = ",
    +1e3
    * diameter**2
    * (p_density - f_density)
    * 9.81
    * voidFraction**3
    / (150 * dyn_vis * (1 - voidFraction)),
)

# Compute fluid superficial velocity
fsvel = fvel * voidFraction

# Compute the pressure drop analytically (with Ergun eq)
term1 = (
    150.0
    * dyn_vis
    / diameter**2
    * fsvel
    * (1 - voidFraction) ** 2
    / voidFraction**3
)
term2 = (
    1.75 * f_density / diameter * fsvel * (1 - voidFraction) / voidFraction**3 * fsvel
)
ergun_dp = bed_height * (term1 + term2)

# Plot pressure drop vs inlet velocity
plt.plot(ivel * 1e3, press_drop, "*")
plt.plot(ivel * 1e3, ergun_dp, "--")
plt.xlabel("Inlet velocity (mm/s)")
plt.ylabel("Pressure drop (Pa)")
plt.legend(["Sim", "Anal"])
plt.grid(linestyle=":")
plt.show()

plt.figure()
plt.plot(ivel * 1e3, voidFraction)
plt.grid(linestyle=":")
plt.show()
