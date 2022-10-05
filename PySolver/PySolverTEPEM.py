import numpy as np
from mpi4py import MPI
import petsc4py
petsc4py.init()
from petsc4py import PETSc
import time
import os

def FindString( io, string, loc ):
	io.seek( 0 )
	while True:
		str0 =  io.readline().split()
		if len( str0 ) > loc and str0[loc] == string:
			return 0
		if len( str0 ) > loc and str0[loc] == '*END':
			return 1
	return 1

if __name__ == "__main__":
	world_comm = MPI.COMM_WORLD
	world_size = world_comm.Get_size()
	my_rank = world_comm.Get_rank()

	xData = np.zeros((5,3))

	if my_rank == 0:
		print(' ')
		print('**********************************************************************')
		print('************************ Initializing solver *************************')
		meshFile = open('Mesh.txt', 'r')

		str0 = '*NODAL'; iError = FindString( meshFile, str0, 0); idofs = int( meshFile.readline().split()[0] )
		str0 = '*DIMEN'; iError = FindString( meshFile, str0, 0); nDim  = int( meshFile.readline().split()[0] )
		str0 = '*COORDINATES'; iError = FindString( meshFile, str0, 0); 
		noCoor = int( meshFile.readline().split()[0] )
		xData = np.zeros((noCoor,nDim))
		for icont in range( noCoor ):
			xData[ icont ] = np.float_( meshFile.readline().split() )

		str0 = 'GROUPS'; iError = FindString( meshFile, str0, 1); 
		meshFile.readline().split()[0]
		noElement = int( meshFile.readline().split()[1] )

		str0 = 'TYPE'; iError = FindString( meshFile, str0, 1); noTypeElement = 5 * [0]
		# noTypeElement -> [ volume, flow, pressure, resistance, no-slip ]
		for icont in range( noElement ):
			idElement = int(meshFile.readline().split()[0])
			noTypeElement[ idElement - 1 ] = noTypeElement[ idElement - 1 ] + 1 

		print('\n*Solving in parallel with %d processors'%world_size)

		print('\n*Mesh: ')
		print('\t Coordinates:   %d'%noCoor)
		print('\t Dimension:     %d'%nDim)
		print('\t DoFs per node: %d'%idofs)

		print('\n*Elements: ')
		print('\t Fluid volume:  %d'%noTypeElement[0])
		print('\t BC flow:       %d'%noTypeElement[1])
		print('\t BC pressure:   %d'%noTypeElement[2])
		print('\t BC resistance: %d'%noTypeElement[3])
		print('\t BC no-slip:    %d'%noTypeElement[4])


		#idElement = PETSc.Vec().createSeq( noElement )
		#for icont in range( noElement ):
		#	idElement.setValues( icont,)
		# idElement.getValues( index ) 

		#
		#print( meshFile.readline() )
		meshFile.close()

	world_comm.Barrier()

	PID = os.getpid()
	print( my_rank, PID )

	if my_rank == 0:
		localData = np.array([0.])
		world_comm.Scatter( xData[0:5][:], localData, root = 0)
		print( my_rank, localData )


	#localData = np.array([0.])
	#world_comm.Scatter( xData[0:5][:], localData, root = 0)
	#localData = world_comm.bcast( xData[my_rank], root = 0 )
	#localData = world_comm.sendrecv( xData[my_rank], dest = my_rank, source = 0)
	#print( my_rank, localData )

	world_comm.Barrier()	
	if my_rank == 0:
		print( xData[0:5])
		del( xData )
	#print("World Size: " + str(world_size) + "   " + "Rank: " + str(my_rank))