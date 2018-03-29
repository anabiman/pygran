from mpi4py import MPI
from PyGran.Simulator.engine_liggghts import liggghts

def foo():

	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()
	tProcs = comm.Get_size()

      	if rank < 1: #pProcs * (i + 1):

		color = 0
	else:
		color = 1

        split = MPI.COMM_WORLD.Split(color=color, key=rank)
        split.barrier()

	lmp = liggghts(comm=split, library='/home/levnon/Downloads/LIGGGHTS-PUBLIC/src/libliggghts.so', cmdargs=['-log', 'liggghts{}.log'.format(color)])
	lmp.lib.lammps_command(lmp.lmp, 'shell pwd')

foo()
