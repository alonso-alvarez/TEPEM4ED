import numpy as np
from mpi4py import MPI
import scipy.sparse as sp
import math as math
from petsc4py import PETSc


world_comm = MPI.COMM_WORLD
my_rank = world_comm.Get_rank()
if my_rank == 0:
	n=10

	A = sp.csr_matrix((n,n),dtype=float)
	print(A.shape)
	A[1:5,:]=1+5*math.pi

	p1=A.indptr
	p2=A.indices
	p3=A.data

	print( A.shape )
	print( p1, type(p1) )
	print( p2, type(p2) )
	print( p3, type(p3) )
	petsc_mat = PETSc.Mat().createAIJ(size=A.shape,csr=(p1,p2,p3))
	print( petsc_mat )