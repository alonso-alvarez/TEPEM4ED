import numpy as np

def SymbolicFluidElement( idofs, globalInc, row_local, col_local, local_count ):
	indexP = [ 0, 2, 6, 8, 18, 20, 24, 26 ]

	for irow in range( 27 ):
		ipRow = globalInc[ irow ]*idofs 

		for icol in range( 27 ):
			ipCol = globalInc[ icol ]*idofs

			row_local[ local_count ] = ipRow + 0; col_local[ local_count ] = ipCol + 0; local_count +=1
			row_local[ local_count ] = ipRow + 0; col_local[ local_count ] = ipCol + 1; local_count +=1
			row_local[ local_count ] = ipRow + 0; col_local[ local_count ] = ipCol + 2; local_count +=1

			row_local[ local_count ] = ipRow + 1; col_local[ local_count ] = ipCol + 0; local_count +=1
			row_local[ local_count ] = ipRow + 1; col_local[ local_count ] = ipCol + 1; local_count +=1
			row_local[ local_count ] = ipRow + 1; col_local[ local_count ] = ipCol + 2; local_count +=1

			row_local[ local_count ] = ipRow + 2; col_local[ local_count ] = ipCol + 0; local_count +=1
			row_local[ local_count ] = ipRow + 2; col_local[ local_count ] = ipCol + 1; local_count +=1
			row_local[ local_count ] = ipRow + 2; col_local[ local_count ] = ipCol + 2; local_count +=1

		if irow < 8:
			ipRow = globalInc[ indexP[ irow ] ]*idofs 
			for icol in range( 8 ):
				ipCol = globalInc[ indexP[ icol ] ]*idofs

				row_local[ local_count ] = ipRow + 0; col_local[ local_count ] = ipCol + 3; local_count += 1
				row_local[ local_count ] = ipRow + 1; col_local[ local_count ] = ipCol + 3; local_count += 1
				row_local[ local_count ] = ipRow + 2; col_local[ local_count ] = ipCol + 3; local_count += 1

				row_local[ local_count ] = ipRow + 3; col_local[ local_count ] = ipCol + 0; local_count += 1
				row_local[ local_count ] = ipRow + 3; col_local[ local_count ] = ipCol + 1; local_count += 1
				row_local[ local_count ] = ipRow + 3; col_local[ local_count ] = ipCol + 2; local_count += 1

				row_local[ local_count ] = ipRow + 3; col_local[ local_count ] = ipCol + 3; local_count += 1	

	return local_count 

def SymbolicElement_CSR( idofs, elem_type, inc ):
	if elem_type == 0:
		return Symbolic_FluidElement( idofs, inc )

	return 0, 0 

def Symbolic_FluidElement_v0( idofs, inc ):
	AE = np.ones( (27*idofs,27*idofs), dtype = int )

	for irow in range( 27 ):
		ipRow = irow*idofs + 3 
		if irow == 0 or irow == 2 or irow == 6 or irow == 8:
			continue
		if irow == 18 or irow == 20 or irow == 24 or irow ==26:
			continue
		AE[ ipRow,: ] = 0; AE[ :,ipRow ] = 0

	return sum( sum( AE )), AE 

def Symbolic_FluidElement( idofs, incidence ):
	tempsize = 27*idofs*27*idofs
	
	local_count = 0; 
	row = tempsize * [ 0 ]; col = tempsize * [ 0 ]	
	for irow in range( 27 ):
		for icol in range( 27 ):
			for m in range( 3 ):
				for n in range( 3 ):
					row[ local_count ] = incidence[ irow ]*idofs + m
					col[ local_count ] = incidence[ icol ]*idofs + n 
					local_count = local_count + 1

	for irow in [ 0, 2, 6, 8, 18, 20, 24, 26 ]:
		for icol in range( 27 ):
			for m in range( 3 ):
				row[ local_count ] = incidence[ irow ]*idofs + 3 
				col[ local_count ] = incidence[ icol ]*idofs + m 
				local_count = local_count + 1

				row[ local_count ] = incidence[ icol ]*idofs + m 
				col[ local_count ] = incidence[ irow ]*idofs + 3 
				local_count = local_count + 1

		for icol in [ 0, 2, 6, 8, 18, 20, 24, 26 ]:
			row[ local_count ] = incidence[ irow ]*idofs + 3
			col[ local_count ] = incidence[ icol ]*idofs + 3
			local_count = local_count + 1

	return local_count, row[:local_count], col[:local_count]
