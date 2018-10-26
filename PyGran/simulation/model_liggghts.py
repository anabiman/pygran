'''
Created on July 1, 2016
@author: Andrew Abi-Mansour
'''

# !/usr/bin/python
# -*- coding: utf8 -*-
# -------------------------------------------------------------------------
#
#   Python module for analyzing contact models for DEM simulations
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

# ----------------------------------------------

model_specs = {'name':'MY_MODEL',\
			   'type':'NORMAL', \
			   'params': {'Yeff': 'double **', 'Geff': 'double **' }, \
			   'registerOnOff': {'tangential_damping': 'bool'}, \

			   }

model_specs['sNAME'] = model_specs['name'].lower()
model_specs['mNumber'] = 4

initial_params = {}

# Read core params 
for param in model_specs['params']:
	dtype = model_specs['params'][param].split()
	if len(dtype) == 1:
		initial_params[param] = (model_specs['params'][param], 0)
	elif dtype[-1] == '*':
		initial_params[param] = (model_specs['params'][param], 'NULL')
	elif dtype[-1] == '**':
		initial_params[param] = (model_specs['params'][param], 'NULL')
	elif dtype[-1] == '***':
		initial_params[param] = (model_specs['params'][param], 'NULL')

# Read registerOnOff settings
for param in model_specs['registerOnOff']:
	dtype = model_specs['registerOnOff'][param].split()
	if len(dtype) == 1:
		initial_params[param] = (model_specs['registerOnOff'][param], 0)
	else:
		raise TypeError('Input for registerOnOff should be a boolean variable.')

# Write C++ code for initial_params

key_params = initial_params.keys()
key_code = ''
protected_vars = ''

# Concetenate key params with curly brackets to be used to fill in the code below
for i in range(len(key_params)):
	key_code = key_code + '          {}({}),\n'.format(key_params[i], initial_params[key_params[i]][-1])
	protected_vars = protected_vars + '    {} {};,\n'.format(initial_params[key_params[i]][0], key_params[i])

register_code = ''

# Write C++ code for registerOnOff
for param in model_specs['registerOnOff']:
	register_code = register_code + '	settings.registerOnOff("{}", {}, false);\n'.format(param, param)

code = """
#ifdef NORMAL_MODEL
{type}_MODEL({name},{sNAME},{mNumber})
#else
#ifndef {type}_MODEL_{name}_H_
#define {type}_MODEL_{name}_H_
#include "global_properties.h"
#include <math.h>

namespace LIGGGHTS {{

namespace ContactModels
{{
  template<>
  class NormalModel<{name}> : protected Pointers
  {{
  public:
    static const int MASK = CM_REGISTER_SETTINGS | CM_CONNECT_TO_PROPERTIES | CM_SURFACES_INTERSECT;

   	NormalModel(LAMMPS * lmp, IContactHistorySetup* hsetup, class ContactModelBase *c) : Pointers(lmp), \n""".format(**model_specs) + \
key_code + \
   	"""          cmb(c)
    {

    void registerSettings(Settings & settings)
    {
    """ + \
register_code + \
	"""
    }

  protected:
    class ContactModelBase *cmb; \n""" + \
protected_vars + \
"""
  }};

}}

}}
#endif
#endif
""".format(**model_specs)



