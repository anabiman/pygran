'''
Created on March 06, 2018
@author: Andrew Abi-Mansour
'''

# !/usr/bin/python
# -*- coding: utf8 -*-
# -------------------------------------------------------------------------
#
#   Python module for creating the basic DEM (Granular) object for analysis
#
# --------------------------------------------------------------------------
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 2 of the License, or
#   (at your option) any later version.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

# -------------------------------------------------------------------------
import matplotlib.pylab as plt
import matplotlib
import matplotlib.ticker as ticker
import types

import numpy as np
from PyGran.Tools import conversion

font = {'size': 13}

matplotlib.rc('font', **font)
matplotlib.rc('text', usetex=True)

def _fmt(x, pos):
	""" An internal function that i used to create sci notation for colorbars """
	a, b = '{:.2e}'.format(x).split('e')
	b = int(b)
	return r'${} \times 10^{{{}}}$'.format(a, b)

def _initialize(Particles, fig, value, subplot, axes, **args):
	"""
	An internal function that initializes fig and extracts target value.
	"""

	z = None

	if 'figsize' not in args:
		figsize = (8, 6)
	else:
		figsize = args['figsize']

	if 'dpi' not in args:
		dpi = 80

	if not fig:
		fig = plt.figure(figsize=figsize, dpi=dpi)

	if value:
		if isinstance(value, tuple):
			value, z = value
		elif isinstance(value, str):
			if len(value) > 1:
				z = Particles.data['{}'.format(value)]
			elif len(value) == 1:
				z = Particles.data['{}{}'.format(value,axes[0])], Particles.data['{}{}'.format(value, axes[1])]
		else:
			raise IOError('value must be a string or a tuple')

	ax = fig.add_subplot(subplot)

	x,y = Particles.data['{}'.format(axes[0])], Particles.data['{}'.format(axes[1])]

	return fig, ax, value, x, y, z

def quiver(Particles, value=None, axes='xy', title=None, color='k', units='xy', scale=None, cmap='seismic', subplot=111, fig=None, **args):
	"""
	Plots a 2D field of arrows for a set of particles

	@Particles: PyGran.analysis.SubSystem object

	@[value]: a string attribute such as 'v' (velocity), 'f' (force), 't' (torque), or ... any attribute contained in Particles.
	or a tuple ('attr', array) with array being a list or numpy array of length natoms 
	@[axes]: a 2-letter string specifying which axes to plot ('xy', 'yz', or 'xz') 
	@[title]: string specifying the title of the plot
	@[color]: string specifying the color of the arrows (default black)
	@[radius]: radii of the particles  (numpy array)
	"""

	fig, ax, value, x, y, z = _initialize(Particles, fig, value, subplot, axes, **args)

	if value:
		vx, vy = z
		v = np.sqrt(vx**2+vy**2)
		Q = ax.quiver(x, y, vx, vy, v, units=units, cmap=cmap)
	else:
		Q = None

	if scale:
		ax.scatter(x, y, color=color, s=Particles.radius * scale)
	else:
		ax.scatter(x, y, color=color, s=Particles.radius)

	rmax = Particles.radius.max()

	format(Particles, axes, ax, fig, title, value, x.min() - rmax, x.max() + rmax, y.min() - rmax, y.max() + rmax, Q, **args)
	fig.show()

	return fig

def pcolor(Particles, value, axes='xy', title=None, cmap='autumn', subplot=111, fig=None, **args):
	"""
	Plots a 2D pcolor for a set of particles

	@Particles: PyGran.analysis.SubSystem object
	@value: a string attribute such as 'vx' (x-comp velocity), 'fy' (y-comp force), 'tz' (z-comp torque), or ... any attribute contained in Particles.
	or a tuple ('attr', array) with array being a list or numpy array of length natoms

	@[axes]: a 2-letter string specifying which axes to plot ('xy', 'yz', or 'xz') 
	@[title]: string specifying the title of the plot
	@[color]: string specifying the color of the arrows (default black)
	@[radius]: radii of the particles  (numpy array)
	"""


	fig, ax, value, x, y, z = _initialize(Particles, fig, value, subplot, axes, **args)

	ymin, ymax = y.min(), y.max()
	xmin, xmax = x.min(), x.max()

	dr = Particles.radius.min()
	X,Y = np.mgrid[slice(xmin, xmax + dr, dr), slice(ymin, ymax + dr, dr)]

	Z = X * 0
	xi = np.array((x - xmin) / dr, dtype='int64')
	yi = np.array((y - ymin) / dr, dtype='int64')

	for i in range(Particles.natoms):
		Z[xi[i],yi[i]] = z[i]

	Q = ax.pcolor(X, Y, Z, cmap=cmap, vmin=z.min(), vmax=z.max())

	format(Particles, axes, ax, fig, title, value, xmin, xmax, ymin, ymax, Q, **args)
	fig.show()

	return fig

def format(Particles, axes, ax, fig, title, value, xmin, xmax, ymin, ymax, cbar, **args):

	if title:
		ax.set_title(title)

	ax.set_xlabel('${}$ ({})'.format(axes[0], conversion[Particles.units()]['distance'][1]), fontsize=16)
	ax.set_ylabel('${}$ ({})'.format(axes[1], conversion[Particles.units()]['distance'][1]), fontsize=16)

	ax.set_xlim(xmin, xmax)
	ax.set_ylim(ymin, ymax)

	if cbar:
		cbar = fig.colorbar(cbar, extend='max')

		if 'cbar_title' in args:
			cbar_title = args['cbar_title']
		else:
			if len(value) > 1:
				cbar_title = '${}_{}$'.format(value[0], value[1])
			else:
				cbar_title = '${}$'.format(value)

		cbar.ax.set_ylabel(cbar_title, fontsize=16)

		cbar.formatter.set_powerlimits((0, 0))
		cbar.update_ticks()

	ax.grid(linestyle=':')

	ax.set_xticks(np.linspace(xmin, xmax, 8))
	ax.set_yticks(np.linspace(ymin, ymax, 8))

	ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
	ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

def timePlot(System, attr, framei=0, framef=None, metric=None, title=None, xlabel=None, ylabel=None, xscale=None, yscale=None, **args):

	time, var = _timeExtract(System, attr, framei, framef, metric)

	if xscale:
		time = time * xscale

	if yscale:
		var = var * yscale

	if 'figsize' not in args:
		figsize = (8, 6)
	else:
		figsize = args['figsize']

	if 'dpi' not in args:
		dpi = 80
	else:
		dpi = args['dpi']

	if 'lineStyle' in args:
		lineStyle = args['lineStyle']
	else:
		lineStyle = '-'

	if 'marker' in args:
		marker = args['marker']
	else:
		marker = 'o'

	fig = plt.figure(figsize=figsize, dpi=dpi)

	ax = fig.gca()

	if not xlabel:
		xlabel = 'Time'

	if not ylabel:
		ylabel = attr

	ax.set_xlabel('Time ({})'.format(conversion[System.Particles.units()]['time'][1]), fontsize=16)
	ax.set_ylabel(ylabel, fontsize=16)

	plt.plot(time, var, lineStyle=lineStyle, marker=marker)

	ax.grid(linestyle=':')

	#ax.set_xticks(np.linspace(time.min(), time.max(), 8))
	#ax.set_yticks(np.linspace(var.min(), var.max(), 8))

	ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
	ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

	fig.show()

	return fig

def _timeExtract(System, attr, framei=0, framef=None, metric=None):

	Particles = System.Particles

	if not hasattr(Particles, attr):
		raise IOError('{} not found in Particles class.'.format(attr))

	System.goto(framei)
	System.skip()
	var = []
	time = []

	if not metric:
		if not isinstance(getattr(Particles, attr), types.MethodType):
			if isinstance(Particles.data[attr], np.ndarray):
				metric = 'mean'
			# otherwise, the user requested a function call

	for ts in System:

		if framef:
			if ts >= framef:
				return np.array(var)

		if metric:
			method_to_call = getattr(np, metric)
			var.append(method_to_call(System.Particles.data[attr]))
		else:
			if isinstance(getattr(System.Particles, attr), types.MethodType):
				var.append(getattr(System.Particles, attr)())
			else:
				var.append(System.Particles.data[attr])

		time.append(ts)

	return np.array(time), np.array(var)