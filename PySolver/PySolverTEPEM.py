import numpy as np
from mpi4py import MPI
import petsc4py
petsc4py.init()
from petsc4py import PETSc
from scipy.sparse import coo_matrix, csr_matrix
from timeit import default_timer as timer
import os
import sys

sys.path.append('./PySolver/')
from FluidFlowElements import *

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
	idofs = None; nDim = None; local_count = None; row_global = None; col_global = None
	if my_rank == 0:
		print(' ')
		print('**********************************************************************')
		print('************************ Initializing solver *************************')
		tic = timer()
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

		meshFile.close()

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

		noGlobalEl = sum( noTypeElement )
		noLocalEl = noTypeElement[0]//world_comm.size
		system_size = noCoor*idofs
		print('\n*Elements per process: %d'%noLocalEl)
		print('*Time reading input files: %4.2f sec'%(timer() - tic))

		print('\n**********************************************************************')
		print('**************** PreProcessing - allocating space ********************')

		local_count, *_ = SymbolicElement_CSR( idofs, elem_type = 0, inc = [*range(27)] );
		row_global = np.zeros( noTypeElement[0] * local_count, dtype = int )
		col_global = np.zeros( noTypeElement[0] * local_count, dtype = int )

	idofs = world_comm.bcast( idofs, root = 0 )
	nDim = world_comm.bcast( nDim, root = 0 )
	incidence = world_comm.bcast( incidence, root = 0 )

	noLocalEl = world_comm.bcast( noLocalEl, root = 0 )
	system_size = world_comm.bcast( system_size, root = 0 )
	local_count = world_comm.bcast( local_count, root = 0 )
	
	row_local = np.zeros( noLocalEl * local_count, dtype = int )
	col_local = np.zeros( noLocalEl * local_count, dtype = int )

	range_a = my_rank * noLocalEl 
	range_b = my_rank * noLocalEl + noLocalEl

	tic = timer()
	local_count = 0
	for iBC in range( 1 ):
		for ielem in range( range_a,range_b ):
			globalInc = incidence[ielem] - 1
			
			_, row, col = SymbolicElement_CSR( idofs, elem_type = iBC, inc = globalInc );
			for irow in range( len(row) ):
				row_local[ local_count ] = row[irow]; col_local[ local_count ] = col[irow]
				local_count = local_count + 1
	del( row,col )

	if my_rank == 0:
		local_count = 0
		for irank in range(world_size):
			if irank > 0:
				world_comm.Recv( row_local, source = irank, tag = 0 )
				world_comm.Recv( col_local, source = irank, tag = 1 )

			print( 'Receiving from %d: %d %d %d %d'%( irank,row_local[0],row_local[-1],col_local[0],col_local[-1]) ); 
			for icount in range( len(row_local) ):
				row_global[ local_count ] = row_local[ icount ]
				col_global[ local_count ] = col_local[ icount ]
				local_count = local_count + 1

		AE_global = csr_matrix((np.zeros(len(row_global), dtype = float), (row_global,col_global) ))
		BE_global = np.zeros( system_size, dtype = float )
		print('\n*Time assembling symbolic structure: %4.2f sec'%(timer() - tic))
		print('Symbolic matrix size: %4.2f bytes'%(sys.getsizeof(AE_global) ) )
		print('Symbolic vector size: %4.2f Mb\n'%(sys.getsizeof(BE_global)/(1024*1024) ) )
	else:
		#print( 'Sending from %d to 0: %d %d %d %d'%(my_rank,row_local[0],row_local[-1],col_local[0],col_local[-1]) ); 
		world_comm.Send( [row_local,MPI.INT], dest = 0, tag = 0 );
		world_comm.Send( [col_local,MPI.INT], dest = 0, tag = 1 );

	world_comm.Barrier()
	del( row_local, col_local, row_global, col_global )
	del( xData, incidence )

	#for icont in range( 1, noLocalEl+1 ):
		#if my_rank == 0:
		#	for ilocal in range( world_size ):
		#		local_index = ilocal * noLocalEl + icont 
		#		world_comm.send( incidence[ local_index - 1 ], dest = ilocal, tag = 0 )
		#local_inc = world_comm.recv( source = 0, tag = 0 )  ## [1, 2, 3, ..., p ] \in N^27

	# 	## Build symbolic structure of AE - count elements (i,j) = 1
	# 	tic = timer()
	# 	local_count, *_ = SymbolicElement_CSR( idofs, elem_type = 0, inc = [*range(27)] );

	# 	print( '****', local_count)

	# 	row = noTypeElement[0] * local_count * [0]
	# 	col = noTypeElement[0] * local_count * [0]

	# 	print( 'Symbolic matrix size: %4.2f Mb'%(sys.getsizeof(row)/(1024*1024) ) )
	# 	local_count = 0
	# 	for iBC in range( 1 ):
	# 		ieleminit = sum( noTypeElement[0:iBC])
	# 		_, AElocal = SymbolicElement_CSR( idofs, elem_type = iBC, inc = [*range(27)] );

	# 		rowrange, colrange = np.where( AElocal == 1 )
	# 		for ielem in range( 1000 ): #ieleminit, ieleminit + noTypeElement[iBC] ):
	# 			globalInc = incidence[ielem] - 1
				
	# 			for irow, icol in zip( rowrange,colrange ):
	# 				ipRow = irow//idofs; m = idofs*( irow/idofs - irow//idofs)
	# 				ipCol = icol//idofs; n = idofs*( icol/idofs - icol//idofs)
	# 				row[ local_count ] = globalInc[ ipRow ]*idofs + m
	# 				col[ local_count ] = globalInc[ ipCol ]*idofs + n 
	# 				local_count = local_count + 1

	# 		#print( AElocal.shape)
	# 		#for ielem in range( 1000 ): #ieleminit, ieleminit + noTypeElement[iBC] ):
	# 		#	globalInc = incidence[ielem] - 1
	# 		#	
	# 		#	for irow in range( len(row_local) ):
	# 		#		row[ local_count ] = row_local[irow]; col[ local_count ] = col_local[irow]
	# 		#		local_count = local_count + 1

	# 	print( 'Symbolic matrix size: %4.2f Mb'%(sys.getsizeof(row)/(1024*1024) ) )

	# 	print('\n*Time assembling symbolic structure: %4.2f sec'%(timer() - tic))
	# 	exit()


	# 		# 	for icol in range( icolmax ):
	# 		# 		globalindex = incidence[ielem][icol] - 1
	# 		# 		for irow in range( idofs ):
	# 		# 			BE[ globalindex*idofs + irow ] = BElocal[ idofs*icol + irow ]

	# 	del( row_local ); del( col_local )

	# 	BE = system_size * [0]; BElocal = 27 * idofs * [ 1 ]
	# 	for inode in range( 1,28 ):
	# 		if( inode == 1 or inode == 3 or inode == 7 or inode == 9 ):
	# 			continue
	# 		if( inode == 19 or inode == 21 or inode == 25 or inode == 27 ):
	# 			continue
	# 		BElocal[ idofs*( inode - 1 ) + idofs - 1 ] = 0 

	# 	tic = timer()
	# 	icolmax = 27
	# 	for iBC in range( 4 ):
	# 		ieleminit = sum( noTypeElement[0:iBC])

	# 		if iBC > 0:
	# 			icolmax = 10; BElocal = 10 * idofs * [ 1 ]
	# 			for itemp in range( 9 ):
	# 				BElocal[ itemp*idofs + 0 ] = 1
	# 				BElocal[ itemp*idofs + 1 ] = 1
	# 				BElocal[ itemp*idofs + 2 ] = 1
	# 			if iBC == 1 or iBC == 3:
	# 				BElocal[ -1 ] = 1

	# 		for ielem in range( ieleminit, ieleminit + noTypeElement[iBC] ):
	# 			for icol in range( icolmax ):
	# 				globalindex = incidence[ielem][icol] - 1
	# 				for irow in range( idofs ):
	# 					BE[ globalindex*idofs + irow ] = BElocal[ idofs*icol + irow ]

		

	# 	del( BElocal ); del( AElocal )
	# 	print( system_size, sum(BE) )
	# 	AEglobal = csr_matrix((len(BE)*[0], (BE, BE)), shape=(3, 3))
	# 	BEglobal = len(BE) * [ 0.0 ]
	# 	del( BE )

	# world_comm.Barrier()

	# noLocalEl = world_comm.bcast( noLocalEl, root = 0 )
	# system_size = world_comm.bcast( system_size, root = 0 )

	#for icont in range( 1, noLocalEl+1 ):
		#if my_rank == 0:
		#	for ilocal in range( world_size ):
		#		local_index = ilocal * noLocalEl + icont 
		#		world_comm.send( incidence[ local_index - 1 ], dest = ilocal, tag = 0 )
		#local_inc = world_comm.recv( source = 0, tag = 0 )  ## [1, 2, 3, ..., p ] \in N^27


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

	#local_a = my_rank * noLocalEl +  1
	#local_b = my_rank * noLocalEl + noLocalEl
	#print( my_rank, local_a, local_b )


	#BE = np.zeros( system_size, dtype = float )


	#print( my_rank, system_size )
	#AE = PETSc.Mat().createAIJ([system_size,system_size]); AE.setUp(); AE.assemble()
	#BE = PETSc.Vec().createSeq( system_size )
	#AE.setValue(0,0,0.); AE.assemble()
	#world_comm.Barrier()
	#print( '***', my_rank, AE.getValues(0,0))

	# for icont in range( 1,noLocalEl+1 ):
	# 	if my_rank == 0:
	# 		for ilocal in range( world_size ):
	# 			local_index = ilocal * noLocalEl + icont
	# 			world_comm.send( incidence[ local_index - 1 ], dest = ilocal, tag = 0 )
	# 			world_comm.send( xData[ incidence[ local_index - 1 ] - 1 ], dest = ilocal, tag = 1 )
	# 	local_inc = world_comm.recv( source = 0, tag = 0 )
	# 	local_coor= world_comm.recv( source = 0, tag = 1 )

		## Realizar calculo local 
		## FluidVolume( local_inc, local_coor, localA(27*idoft,27*idoft), localB(27*idoft)) )

	#AE.setValue( my_rank, my_rank, 2*my_rank + 1); AE.assemble()
	
	#	#for icont in range( 10 ): #system_size ):
	#	icont = 0
	#	AE.setValue( icont,icont,icont )
	#	AE.assemble()
	#	


	#world_comm.Barrier()	
	#AE = world_comm.reduce( AE, root = 0 )
	#BE = world_comm.reduce( BE, root = 0 )
	#if my_rank == 0:
	#	local_a = noLocalEl * world_comm.size + 1; local_b = noTypeElement[0]
	#	print( '**', my_rank, local_a, local_b )
		#for icont in range( local_a, local_b + 1 ):
		#	print( '**', my_rank, icont, local_a, local_b )

		#usol = PETSC.Vec().createSeq( system_size )
		#ksp = PETSC.KSP().create(); ksp.setOperators( AE )
		#ksp.setFromOptions()
		#ksp.solve( BE, usol )
		#print( usol.getArray() )

	#PID = os.getpid(); print( my_rank, PID )

	# world_comm.Barrier()
	# del( xData ); del( incidence ); del( noLocalEl )
	# if my_rank == 0:
	# 	del( AEglobal )
	# 	del( BEglobal )
	# #del( AE ); del( BE ); del( usol )
