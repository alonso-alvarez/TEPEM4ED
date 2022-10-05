import numpy
from mpi4py import MPI
import time

comm = MPI.COMM_WORLD
my_rank = comm.Get_rank()
p = comm.Get_size()

if my_rank != 0:
	message = "Hello from "+str(my_rank)
	comm.send( message, dest=0 )
else:
	for procid in range( 1,p ):
		message = comm.recv(source=procid)
		print("process 0 receives message from process",procid,":", message)

