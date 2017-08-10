'''
Created on July 11, 2016
@author: Andrew Abi-Mansour
'''

#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

# -------------------------------------------------------------------------


class Temporal(object):
	def __init__(self):
		pass
		
def computeFlow(data, density, t0 = 0, N0 = 0, sel = None, dt = 1e-4):
	"""
	Computes flow rate: N/t for a selection *sel*
	@ data: list of dictionaries containing simulation and particle data (box size, x,y,z, etc.)

	TODO: Get this working for a particle distribution
	"""

	if N0 == None or t0 == None:
		return 0
	else:
		mass = density * 4.0 / 3.0 * np.pi * (len(sel) - N0) * np.mean(data['radius'][sel])**3.0
		return - mass / ((data['TIMESTEP'] - t0) * dt)