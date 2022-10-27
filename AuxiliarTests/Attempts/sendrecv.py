# send_array
from mpi4py import MPI
import numpy 

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

if rank == 0:
  data = numpy.arange(5,dtype=numpy.float64)
  comm.Send([data,3,MPI.DOUBLE],1,11)
elif rank == 1:
  data = 10.*numpy.arange(5,dtype=numpy.float64)
  print('on task',rank,'before Recv:   data = ',data)
  comm.Recv(data,source=0,tag=11)
  print('on task',rank,'after Recv:    data = ',data)