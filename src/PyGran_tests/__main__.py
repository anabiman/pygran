#  -*- coding: utf8 -*-
'''
  Created on November 22, 2018
  @author: Andrew Abi-Mansour

  This is the 
   __________         ________                     
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

 -------------------------------------------------------------------------
  Main module for running the various demo scripts in PyGran
 -------------------------------------------------------------------------

 '''

import sys, os
from importlib import import_module
import glob
import logging

try:
  from mpi4py import MPI
except:
  class MPI(object):
    class COMM_WORLD(object):
      def Get_rank(): return 0

  logging.warning("Failed to import MPI from mpi4py.")

if __name__ == '__main__':

  sdir, _ = __file__.split(__name__.split('PyGran_tests.')[-1] +'.py')

  # In case you forgot how to use this ... run -h
  if sys.argv[1] == '-h' or sys.argv[1] == '--help':

    if not MPI.COMM_WORLD.Get_rank():
      possible_dirs = glob.glob(os.path.join(sdir, 'scripts', '*'))
      actual_dirs = []

      actual_dirs = [pdir.split('/scripts/')[-1] for pdir in possible_dirs if os.path.isdir(pdir)]

      print('Available options:\n==================')
      for dtype in actual_dirs:
        py_dirs = glob.glob(os.path.join(sdir, 'scripts', dtype, '*'))
        py_dirs = [py_dir.split(dtype + '/')[-1] for py_dir in py_dirs if os.path.isdir(py_dir)]
        print( '{}: '.format(dtype), ('{}, ' * len(py_dirs)).format(*py_dirs)[:-2] ) # get rid of the last comma

    sys.exit()

  if sys.argv[1] == '--all' or sys.argv[1] == '-a':
    print('All tests not yet supported.')
    raise NotImplementedError

  if len(sys.argv) != 3:
    if not MPI.COMM_WORLD.Get_rank():
      raise ValueError('PyGran.demo takes only 2 input arguments: type demo_name. See help (-h) for options. ')
  else:
    try:
      demo = import_module('PyGran_tests.scripts.' + sys.argv[1] + '.' + sys.argv[2] + '.' + sys.argv[2])

      # see if demo contains a mesh file; adjust its location if so
      if hasattr(demo, 'params'):
        if 'mesh' in demo.params:
          for key in demo.params['mesh']:
            demo.params['mesh'][key]['file'] = os.path.join(sdir, 'scripts', sys.argv[1], sys.argv[2], demo.params['mesh'][key]['file'])

        # now run the demo~!
        demo.run(**demo.params)

      else: # we're just running a simple script (no DEM sim)
        demo.run()
    except:
      raise

  sys.exit()
