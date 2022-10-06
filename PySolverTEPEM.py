import numpy as np
from mpi4py import MPI
import petsc4py
petsc4py.init()
from petsc4py import PETSc
from scipy.sparse import coo_matrix, csr_matrix
from timeit import default_timer as timer
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

	system_size = None; noLocalEl = None; incidence = None; xData = None
	if my_rank == 0:
		print(' ')
		print('**********************************************************************')
		print('************************ Initializing solver *************************')
		meshFile = open('Mesh.txt', 'r')

		str0 = '*NODAL'; iError = FindString( meshFile, str0, 0); idofs = int( meshFile.readline().split()[0] )
		str0 = '*DIMEN'; iError = FindString( meshFile, str0, 0); nDim  = int( meshFile.readline().split()[0] )
		str0 = '*COORDINATES'; iError = FindString( meshFile, str0, 0); 

		idofs = 4 ### Modificando para resolver NS sem postprocessamento

		noCoor = int( meshFile.readline().split()[0] )
		xData = np.zeros((noCoor,nDim))
		for icont in range( noCoor ):
			xData[ icont ] = np.float_( meshFile.readline().split() )

		str0 = 'GROUPS'; iError = FindString( meshFile, str0, 1); 
		meshFile.readline().split()[0]
		noElement = int( meshFile.readline().split()[1] )

		str0 = 'TYPE'; iError = FindString( meshFile, str0, 1); noTypeElement = 5 * [0]
		for icont in range( noElement ):
			idElement = int(meshFile.readline().split()[0])
			noTypeElement[ idElement - 1 ] = noTypeElement[ idElement - 1 ] + 1 

		str0 = '*INCIDENCE'; iError = FindString( meshFile, str0, 0); 
		incidence = np.zeros((sum(noTypeElement),27), dtype = int)
		for icont in range( sum(noTypeElement) ):
			inc_local = np.float_( meshFile.readline().split() )
			incidence[ icont ][:len(inc_local)] = inc_local

		print('\n*Solving in parallel with %d processors'%world_size)
		print('\n*Mesh: ')
		print('\t Coordinates:   %d'%noCoor)
		print('\t Dimension:     %d'%nDim)
		print('\t DoFs per node: %d'%idofs)

		print('\n*Elements: ')
		print('\t Fluid volume:  %d'%noTypeElement[0])  # sempre divisivel por 48 !!!
		print('\t BC flow:       %d'%noTypeElement[1])
		print('\t BC pressure:   %d'%noTypeElement[2])
		print('\t BC resistance: %d'%noTypeElement[3])
		print('\t BC no-slip:    %d'%noTypeElement[4])

		noLocalEl = noTypeElement[0]//world_comm.size
		system_size = noCoor*idofs
		print('\n*Elements per process: %d'%noLocalEl)

		meshFile.close()
	world_comm.Barrier()

	noLocalEl = world_comm.bcast( noLocalEl, root = 0 )
	system_size = world_comm.bcast( system_size, root = 0 )

	row = np.array([]) #np.zeros(noCoor*(idofs+1))
	col = np.array([]) #np.zeros(noCoor*(idofs+1))
	for icont in range( 1, 2 ): #noLocalEl+1 ):
		if my_rank == 0:
			for ilocal in range( world_size ):
				local_index = ilocal * noLocalEl + icont 
				world_comm.send( incidence[ local_index - 1 ], dest = ilocal, tag = 0 )
		local_inc = world_comm.recv( source = 0, tag = 0 )  ## [1, 2, 3, ..., p ] \in N^27
		print( local_inc, local_inc[0], local_inc[1] )
		np.append( row, np.array([local_inc[0]]) )
		np.append( col, [local_inc[1]] )

	np.append( row,1)
	np.append( row,2)
	np.append( row,3)
	np.append( row,4)
	print( my_rank, row, col )


	#tic = timer()
	#AE = csr_matrix( ( system_size,system_size), dtype = float )
	#toc = timer()
	#print( my_rank, AE[2,0], AE.get_shape(), toc - tic)

	#row = np.array([0, 0, 1, 2, 2, 2])
	#col = np.array([0, 2, 2, 0, 1, 2])
	#data = np.array([1, 2, 3, 0, 5, 6])
	#local = csr_matrix((data, (row, col)), shape=(3, 3))
	#local[2,0] = 11
	#print( my_rank, local[2,0], local[2,2])

	#row = np.array( system_size )
	#col = np.array( system_size )

	local_a = my_rank * noLocalEl +  1
	local_b = my_rank * noLocalEl + noLocalEl
	#print( my_rank, local_a, local_b )





	BE = np.zeros( system_size, dtype = float )


	#print( my_rank, system_size )
	#AE = PETSc.Mat().createAIJ([system_size,system_size]); AE.setUp(); AE.assemble()
	#BE = PETSc.Vec().createSeq( system_size )
	#AE.setValue(0,0,0.); AE.assemble()
	#world_comm.Barrier()
	#print( '***', my_rank, AE.getValues(0,0))

	for icont in range( 1,noLocalEl+1 ):
		if my_rank == 0:
			for ilocal in range( world_size ):
				local_index = ilocal * noLocalEl + icont
				world_comm.send( incidence[ local_index - 1 ], dest = ilocal, tag = 0 )
				world_comm.send( xData[ incidence[ local_index - 1 ] - 1 ], dest = ilocal, tag = 1 )
		local_inc = world_comm.recv( source = 0, tag = 0 )
		local_coor= world_comm.recv( source = 0, tag = 1 )

		## Realizar calculo local 
		## FluidVolume( local_inc, local_coor, localA(27*idoft,27*idoft), localB(27*idoft)) )

	#AE.setValue( my_rank, my_rank, 2*my_rank + 1); AE.assemble()
	
	#	#for icont in range( 10 ): #system_size ):
	#	icont = 0
	#	AE.setValue( icont,icont,icont )
	#	AE.assemble()
	#	


	world_comm.Barrier()	
	#AE = world_comm.reduce( AE, root = 0 )
	#BE = world_comm.reduce( BE, root = 0 )
	if my_rank == 0:
		local_a = noLocalEl * world_comm.size + 1; local_b = noTypeElement[0]
		print( '**', my_rank, local_a, local_b )
		#for icont in range( local_a, local_b + 1 ):
		#	print( '**', my_rank, icont, local_a, local_b )

		#usol = PETSC.Vec().createSeq( system_size )
		#ksp = PETSC.KSP().create(); ksp.setOperators( AE )
		#ksp.setFromOptions()
		#ksp.solve( BE, usol )
		#print( usol.getArray() )

	#PID = os.getpid(); print( my_rank, PID )

	world_comm.Barrier()
	del( xData ); del( incidence ); del( noLocalEl )
	#del( AE ); del( BE ); del( usol )
