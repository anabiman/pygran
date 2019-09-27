'''
Created on February 28, 2018
@author: Andrew Abi-Mansour
'''

# !/usr/bin/python
# -*- coding: utf8 -*- 
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published byO
#   the Free Software Foundation, either version 2 of the License, or
#   (at your option) any later version.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

# -------------------------------------------------------------------------

import subprocess
import os
from mpi4py import MPI

def visualize(**args):

	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()

	if rank == 0:
		cmd = ['ovito']
		traj = args['traj']
		output = os.path.abspath('.')

		if 'pfile' in traj:
			cmd.append('{}/traj/{}'.format(output, traj['pfile']))

		# Can we even call meshes as args to ovito from the terminal? maybe we dont need this
		if 'mesh' in traj:
			if isinstance(traj['mfile'], list):
				for mfile in traj['mfile']:
					cmd.append(' {}/traj/{}'.format(output, mfile))
			else:
				cmd.append(' {}/traj/{}'.format(output, traj['mfile']))

		return subprocess.Popen(cmd)