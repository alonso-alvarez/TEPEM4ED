import numpy as np
from mpi4py import MPI
import petsc4py
import sys

petsc4py.init(sys.argv)

from petsc4py import PETSc
from scipy.sparse import coo_matrix, csr_matrix
from timeit import default_timer as timer
from scipy.sparse.linalg import gmres
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

	nlocal = None
	nlocal = PETSc.Options().getInt( "n" )
	print( nlocal )

	local_count = 6
	row_local = np.zeros( local_count, dtype = 'int32' )
	col_local = np.zeros( local_count, dtype = 'int32' )
	val_local = np.zeros( local_count, dtype = float )

	local_count = 0
	if my_rank == 0:
		incidence = [ 1, 4, 6 ]
		row_local[ :4 ] = [ 0, 3, 3, 5 ]
		col_local[ :4 ] = [ 0, 3, 5, 3 ]
		val_local[ :4 ] = [ 1, 1, 1, 1 ]
		local_count = 4
	if my_rank == 1:
		incidence = [ 1, 6, 3 ]
		row_local[ :1 ] = [ 5 ]
		col_local[ :1 ] = [ 5 ]
		val_local[ :1 ] = [ 1 ]
		local_count = 1
	if my_rank == 2:
		incidence = [ 4, 2, 6 ]
		row_local[ :6 ] = [ 3, 3, 3, 1, 1, 5 ]
		col_local[ :6 ] = [ 3, 1, 5, 3, 1, 3 ]
		val_local[ :6 ] = [ 1., 1, 1, 1, 1, 1 ]
		local_count = 6
	if my_rank == 3:
		incidence = [ 5, 4, 1 ]
		row_local[ :3 ] = [ 4, 4, 0 ]
		col_local[ :3 ] = [ 4, 0, 4 ]
		val_local[ :3 ] = [ 4., 1., 1. ]
		local_count = 3

	world_comm.Barrier()
	global_count = world_comm.reduce( local_count, root = 0)

	pA = PETSc.Mat().createAIJ( [6,6] ); pA.setUp(); pA.assemble()	
	[ range_a, range_b ] = pA.getOwnershipRange()
	Range_a = world_comm.gather( range_a, root = 0 )
	Range_b = world_comm.gather( range_b, root = 0 )

	pB = PETSc.Vec().create( comm = world_comm ); pB.setType( PETSc.Vec.Type.MPI )
	pB.setSizes( [range_b - range_a, 6] )
	pB.zeroEntries()

	AE_global = None; AE_local = None
	if my_rank == 0:
		print( '---------', Range_a )
		print( '---------', Range_b)

		row_global = np.zeros( global_count, dtype = int )
		col_global = np.zeros( global_count, dtype = int )
		val_global = np.zeros( global_count, dtype = float )

		global_count = 0
		for irank in range( world_size ):
			if irank > 0:
				world_comm.Recv( row_local, source = irank, tag = 0 )
				world_comm.Recv( col_local, source = irank, tag = 1 )
				world_comm.Recv( val_local, source = irank, tag = 2 )
				local_count = world_comm.recv( source = irank, tag = 3 )

			range_a = global_count
			range_b = range_a + local_count
			row_global[ range_a:range_b ] = row_local[:local_count]
			col_global[ range_a:range_b ] = col_local[:local_count]
			val_global[ range_a:range_b ] = val_local[:local_count]
			global_count = range_b

		AE_global = csr_matrix((val_global, (row_global, col_global)), shape=(6, 6), dtype = float )
		BE_global = np.zeros( 6, dtype = float ); BE_global[0] = 1.

		#print( AE_global )
		for irank in range( world_size ):
			range_a = Range_a[ irank ]; range_b = Range_b[ irank ]
			AE_local = AE_global[range_a:range_b][:]
			
			p1 = AE_local.indptr;  world_comm.Send( [ p1,MPI.INT ], dest = irank, tag = 10 )
			p2 = AE_local.indices; world_comm.Send( [ p2,MPI.INT ], dest = irank, tag = 11 )
			p3 = AE_local.data;    world_comm.Send( [ p3,MPI.FLOAT ], dest = irank, tag = 12 )

			Range_a[ irank ] = len( p1 )
			Range_b[ irank ] = len( p2 )

		pB.setValues( range(6),BE_global, addv = PETSc.InsertMode.INSERT_VALUES ) # addv = PETSC.InsertMode.ADD_VALUES
		del( AE_global )
		del( BE_global )
	else:
		world_comm.Send( [row_local,MPI.INT], dest = 0, tag = 0 );
		world_comm.Send( [col_local,MPI.INT], dest = 0, tag = 1 );
		world_comm.Send( [val_local,MPI.FLOAT], dest = 0, tag = 2 );
		world_comm.send( local_count, dest = 0, tag = 3 )

	world_comm.Barrier()

	pB.assemble()

	[ range_a, range_b ] = pA.getOwnershipRange()
	print( '>>>>>>>>', my_rank, range_a, range_b)
	for irange in range( range_a, range_b ):
		print( my_rank, irange )
		#pB.setValue( irange, my_rank + 1, addv = PETSc.InsertMode.ADD_VALUES )
		pB.setValue( irange, my_rank + 1 ) #, addv = PETSc.InsertMode.INSERT_VALUES )

	range_a = world_comm.scatter( Range_a, root = 0 )
	range_b = world_comm.scatter( Range_b, root = 0 )

	p1 = np.zeros( range_a, dtype = 'int32' ); world_comm.Recv( p1, source = 0, tag = 10 )
	p2 = np.zeros( range_b, dtype = 'int32' ); world_comm.Recv( p2, source = 0, tag = 11 )
	p3 = np.zeros( range_b, dtype = float );   world_comm.Recv( p3, source = 0, tag = 12 )

	pA = PETSc.Mat().createAIJ( size = [6,6], csr=(p1,p2,p3) ); pA.setUp(); pA.assemble()

	#pB.setValueLocal(0,-(my_rank+1) )
	#pB.assemble()
	pB.assemble()
	#pB.view()
	[ range_a, range_b ] = pA.getOwnershipRange()
	#print( my_rank, pA.getValues(range(range_a,range_b),range(6)) )
	#print( my_rank, pB.getValues(range(range_a,range_b)))

	x = PETSc.Vec().createMPI( 6 ); x.assemble()


	ksp = PETSc.KSP().create()
	ksp.setType('preonly')
	pc = ksp.getPC()
	pc.setType('none')
	ksp.setOperators( pA )
	ksp.setFromOptions()

	
	#print( '<<<', my_rank, pB.getArray())
	#ksp.solve( pB,x )
	pA.view()
	pB.view()
	#print( '>>>', my_rank, x.getArray())
	pA.mult( pB, x )
	x.view()

	pA.destroy()
	pB.destroy()
	x.destroy()

	exit()

	[ range_a, range_b ] = pA.getOwnershipRange()
	print( '--->', my_rank, range_a, range_b)
	AE_global = world_comm.bcast( AE_global[range_a:range_b][:], root = 0 )


	print( my_rank, AE_global )
	exit()



	AE_local = csr_matrix( (val_local, (row_local,col_local) ), shape = (8,8))
	p1 = AE_local.indptr
	p2 = AE_local.indices 
	p3 = AE_local.data

	print( AE_local )
	print( AE_local.indptr )
	print( AE_local.indices )
	print( AE_local.data )

	world_comm.Barrier()
	pA = PETSc.Mat().createAIJ( [8,20] )
	pA.setUp()
	pA.assemble()
	#print( my_rank, type( row_local), type(row_local[0]))
	#pA.setValue(row_local[0],col_local[0],val_local[0])
	pA.setValue(0,0,my_rank)
	print( '--->', my_rank, pA.getOwnershipRange())
	pA.assemble()

	world_comm.Barrier()
	#print( my_rank, '***',pA.getValues(0,0) )
	#print( my_rank, '***',pA.getValues(3,0) )

	# print(my_rank, 'chegou aqui')
	# pA = PETSc.Mat()
	# pA.create( PETSc.COMM_WORLD )
	# pA.setSizes([6,6])
	# pA.setType('aij')
	# #pA.setPreallocationNNZ(11)
	# #print( '--->', my_rank, pA.getOwnershipRange())
	# pA.setValues( row_local, col_local, val_local )
	# #pA.createAIJ( size=(6,6), csr = (p1,p2,p3) )
	# pA.assemblyBegin()
	# pA.assemblyEnd()
	# print( pA.getSizes() )


	#print( my_rank, local_count )
	global_count = world_comm.reduce( local_count, root = 0)

	if my_rank == 0:
		print( '>>>', my_rank, global_count )
		row_global = np.zeros( global_count, dtype = int )
		col_global = np.zeros( global_count, dtype = int )
		val_global = np.zeros( global_count, dtype = float )

		global_count = 0
		for irank in range( world_size ):
			if irank > 0:
				world_comm.Recv( row_local, source = irank, tag = 0 )
				world_comm.Recv( col_local, source = irank, tag = 1 )
				world_comm.Recv( val_local, source = irank, tag = 2 )
				local_count = world_comm.recv( source = irank, tag = 3 )

			print('--', irank)
			print('	', row_local)
			print('	', col_local)
			print('	', val_local)
			#print( 'Receiving from %d: %d %d %d %d'%( irank,row_local[0],row_local[local_count],col_local[0],col_local[local_count]) ); 
			range_a = global_count
			range_b = range_a + local_count
			row_global[ range_a:range_b ] = row_local[:local_count]
			col_global[ range_a:range_b ] = col_local[:local_count]
			val_global[ range_a:range_b ] = val_local[:local_count]
			global_count = range_b

		print('Final bcast: ')
		print('	', row_global )
		print('	', col_global )
		print('	', val_global )

		AE_global = csr_matrix((val_global, (row_global, col_global)), shape=(6, 6), dtype = float )
		BE_global = np.zeros( 6, dtype = float ); BE_global[0] = 1.

		print( AE_global )

		x, exitCode = gmres( AE_global, BE_global )
		print( x, exitCode )
		# p1 = AE_global.indptr
		# p2 = AE_global.indices 
		# p3 = AE_global.data
		# print( p1, type(p1) )
		# print( p2, type(p2) )
		# print( p3, type(p3) )		

		# print( 'Solving with petsc ')
		# print( 'one')
		# pA = PETSc.Mat().createAIJ( size = AE_global.shape, csr=(p1,p2,p3) )
		# #pB = PETSc.Vec().createSeq( 6 ); pB.setValue(0,1)
		# #pX = PETSc.Vec().createSeq( 6 )
		# print( 'two')

		#ksp = PETSc.KSP().create()
		#ksp.setOperators( pA )
		#print( 'three')

		#ksp.setFromOptions()
		#print("Solving with: ", ksp.getType() )
		#print( 'four')
		#ksp.solve( pB, pX )
		#print( 'five')

		#print( "Converged in ", ksp.getIterationNumber(), 'iterations')
		#print( pX )

    # ksp = PETSc.KSP()
    # ksp.create()
    # ksp.setOperators(pA)
    # ksp.setType('preonly')
    # pc = ksp.getPC()
    # pc.setType('lu')
    # pc.setFactorSolverPackage('mumps')
    # ksp.solve(rhv,x)


		#BE_global = np.zeros( system_size, dtype = float )
		#print('\n*Time assembling symbolic structure: %4.2f sec'%(timer() - tic))
		#print('Symbolic matrix size: %4.2f bytes'%(sys.getsizeof(AE_global) ) )
		#print('Symbolic vector size: %4.2f Mb\n'%(sys.getsizeof(BE_global)/(1024*1024) ) )
	else:
		#print( 'Sending from %d to 0: %d %d %d %d'%(my_rank,row_local[0],row_local[local_count],col_local[0],col_local[local_count]) ); 
		world_comm.Send( [row_local,MPI.INT], dest = 0, tag = 0 );
		world_comm.Send( [col_local,MPI.INT], dest = 0, tag = 1 );
		world_comm.Send( [val_local,MPI.FLOAT], dest = 0, tag = 2 );
		world_comm.send( local_count, dest = 0, tag = 3 )
		#world_comm.Send( [local_count,MPI_INT], dest = 0, tag = 3 );

	# world_comm.Barrier()
	# del( row_local, col_local, row_global, col_global )
	# del( xData, incidence )





	#if my_rank == 0:
#		del( AE_global, BE_global )
