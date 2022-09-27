"""
  Created on Sep 25, 2021
  Author: Andrew Abi-Mansour

  This is the::

  ██████╗ ██╗   ██╗ ██████╗ ██████╗  █████╗ ███╗   ██╗
  ██╔══██╗╚██╗ ██╔╝██╔════╝ ██╔══██╗██╔══██╗████╗  ██║
  ██████╔╝ ╚████╔╝ ██║  ███╗██████╔╝███████║██╔██╗ ██║
  ██╔═══╝   ╚██╔╝  ██║   ██║██╔══██╗██╔══██║██║╚██╗██║
  ██║        ██║   ╚██████╔╝██║  ██║██║  ██║██║ ╚████║
  ╚═╝        ╚═╝    ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═══╝

  DEM simulation and analysis toolkit
  http://www.pygran.org, support@pygran.org

  Core developer and main author:
  Andrew Abi-Mansour, andrew.abi.mansour@pygran.org

  PyGran is open-source, distributed under the terms of the GNU Public
  License, version 2 or later. It is distributed in the hope that it will
  be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
  received a copy of the GNU General Public License along with PyGran.
  If not, see http://www.gnu.org/licenses . See also top-level README
  and LICENSE files.
"""

glass = {
    "youngsModulus": 63e9,
    "poissonsRatio": 0.24,
    "coefficientFriction": 0.5,
    "coefficientRollingFriction": 0.0,
    "cohesionEnergyDensity": 0.05,
    "coefficientRestitution": 0.9,
    "coefficientRollingViscousDamping": 0.1,
    "yieldPress": 62e9,
    "characteristicVelocity": 0.1,
    "density": 2500.0,
}

organic = {
    "youngsModulus": 1e7,
    "poissonsRatio": 0.25,
    "coefficientFriction": 0.5,
    "coefficientRollingFriction": 0.0,
    "cohesionEnergyDensity": 0.0,
    "coefficientRestitution": 0.9,
    "coefficientRollingViscousDamping": 0.1,
    "yieldPress": 2.2e6,
    "characteristicVelocity": 0.1,
    "density": 1000.0,
}

stearicAcid = {
    "youngsModulus": 4.15e7,
    "poissonsRatio": 0.25,
    "coefficientFriction": 0.5,
    "coefficientRollingFriction": 0.0,
    "cohesionEnergyDensity": 0.033,
    "coefficientRestitution": 0.9,
    "coefficientRollingViscousDamping": 0.1,
    "yieldPress": 2.2e6,
    "characteristicVelocity": 0.1,
    "density": 997.164,
}

cohesionless = {
    "youngsModulus": 1e7,
    "poissonsRatio": 0.25,
    "coefficientFriction": 0.5,
    "coefficientRollingFriction": 0.0,
    "coefficientRestitution": 0.9,
    "coefficientRollingViscousDamping": 0.1,
    "yieldPress": 2.2e6,
    "characteristicVelocity": 0.1,
    "density": 1000.0,
}

cohesive = {
    "youngsModulus": 1e7,
    "poissonsRatio": 0.25,
    "coefficientFriction": 0.5,
    "coefficientRollingFriction": 0.0,
    "cohesionEnergyDensity": 2e5,
    "coefficientRestitution": 0.9,
    "coefficientRollingViscousDamping": 0.1,
    "yieldPress": 2.2e6,
    "characteristicVelocity": 0.1,
    "density": 1000.0,
}

steel = {
    "youngsModulus": 2e11,
    "poissonsRatio": 0.3,
    "coefficientFriction": 0.5,
    "coefficientRollingFriction": 0.0,
    "cohesionEnergyDensity": 0.0,
    "coefficientRestitution": 0.9,
    "coefficientRollingViscousDamping": 0.1,
    "yieldPress": 2e10,
    "characteristicVelocity": 0.1,
    "density": 8050,
}

__all__ = ["glass", "stearicAcid", "cohesionless", "cohesive", "steel", "organic"]
