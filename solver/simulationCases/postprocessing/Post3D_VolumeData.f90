!#########################################################
!###   
!#########################################################
 IMPLICIT NONE
 INTEGER, PARAMETER :: Ndim = 3
 
 INTEGER :: iFile, oFile, dFile
 INTEGER :: nci, ncj, nck, nc0, pp, nod, ntotal
 INTEGER :: noPoints, noCells, itime
 INTEGER :: noCoor, idoft, noPipe, par_n, par_z
 INTEGER :: iFindStringInFile, iError, rewindfile 
 INTEGER :: iInc( 300 )
 INTEGER :: iType, pType( 20,3 )
 
 DOUBLE PRECISION :: SerPoints( 54 ), x0( NDim ), det, u0( 7 ), sol1( 300*7 )
 
 INTEGER, ALLOCATABLE :: iPar_Real( : ), nType( : ), nGroup( : )
 DOUBLE PRECISION, ALLOCATABLE :: dataIn( : ), dataOut(:), xlocal( :,: )
 
 CHARACTER*40 :: Str, Path, istr
 
 iFile = 101; oFile = 102; dFile = 103

  CALL GETARG( 1,str ); READ( str,'(I2)') par_n
  CALL GETARG( 2,str ); READ( str,'(I2)') par_z

  IF( par_n .EQ. 0 ) par_n = 3
  IF( par_z .EQ. 0 ) par_z = 3
  !par_n =  4; par_z =   2
 
 !!!!!! Parameters related to the interpolation
 OPEN( iFile, FILE = 'Basparam.txt' )
  Str = '#Data'; iError = iFindStringInFile ( Str, iFile, rewindfile )
  READ( iFile,* ) pp    !! Transversal order
  READ( iFile,* ) noPipe  !! Number of pipe elements
  READ( iFile,* ) noCoor  !! Number of coordinates
  READ( iFile,* ) idoft !! Degrees of freedom

  REWIND( iFile )
  Str = '*ElementLibraryControl'; iError = iFindStringInFile ( Str, iFile, rewindfile )
  READ( iFile,* ) nc0

  Str = '*SubStep'; iError = iFindStringInFile ( Str, iFile, rewindfile )
  READ( iFile,* ) Str
  iType = 0
  DO nck = 1, nc0
    READ( iFile,* ) ncj
    IF( ncj .EQ. 800 ) iType = iType + 1
    IF( nck .EQ. 1 ) nod = ncj
  ENDDO
  READ( iFile,* ) nci
  DO nck = 1, iType
  	pType(nck,1:3) = [nci,nci,nci]
  	WRITE(*,*) pType(nck,1:3)
  ENDDO
  !DO nck = 1, iType
  !  IF( nod .EQ. 806 .OR. nod .EQ. 807) THEN
  !    pType(nck,1:3) = [2,2,2]
  !  ELSE
  !    !READ( iFile,* ) pType( nck,1:3 )
  !    pType(nck,1:3) = [2,2,2]
  !  ENDIF
  !ENDDO
 CLOSE( iFile )

 !!!!!! Reading Param00 - Geometry
 ALLOCATE( iPar_Real( noPipe ) )
 OPEN( iFile, FILE = 'Param00.txt' )
  Str = '*Real'; iError = iFindStringInFile ( Str, iFile, rewindfile )
  READ( iFile,* ) iPar_Real(1:noPipe)
 CLOSE( iFile )

 !!!!!!! Reading Mesh - element type
 ALLOCATE( nType( noPipe ) )
 ALLOCATE( nGroup( noPipe ))
 OPEN( iFile, FILE = 'Mesh.txt')
  Str = '*ELEMENT TYPE'; iError = iFindStringInFile ( Str, iFile, rewindfile )
  READ( iFile,* ) nType( 1:noPipe )

  Str = '*ELEMENT GROUPS'; iError = iFindStringInFile ( Str, iFile, rewindfile )
  READ( iFile,* ) nc0
  READ( iFile,* ) nc0
  READ( iFile,* ) nGroup( 1:noPipe )
 CLOSE( iFile ) 

 
 !!!!!! Local data / per element
 noPoints = par_z * par_n**2; nc0 = 0
 noCells  = noPipe*( par_z - 1 )*( par_n - 1 )**2
 ALLOCATE( xlocal( noPoints,Ndim ) )
 DO nck = 1, par_z
 DO ncj = 1, par_n
 DO nci = 1, par_n
  nc0 = nc0 + 1
  xlocal( nc0,1 ) = ( -par_n - 1 + 2*nci )/( par_n - 1. )
  xlocal( nc0,2 ) = ( -par_n - 1 + 2*ncj )/( par_n - 1. )
  xlocal( nc0,3 ) = ( -par_z - 1 + 2*nck )/( par_z - 1. )
 ENDDO
 ENDDO
 ENDDO
 
 ALLOCATE( dataIn( noCoor*idoft) )
 ALLOCATE( dataOut( noPipe*noPoints*idoft) )

  OPEN( dFile, FILE ='DataOut.txt')
    READ( dFile,* ) iStr
    itime = 1; iError = iFindStringInFile ( iStr, dFile, rewindfile )
    DO WHILE( rewindfile .EQ. 0 )   
      itime = itime + 1
      iError = iFindStringInFile ( iStr, dFile, rewindfile )
    ENDDO
    ntotal = itime
  CLOSE( dFile )

 itime = 0
 OPEN( dFile, FILE = 'DataOut.txt' )
  Str = trim(iStr); iError = iFindStringInFile ( Str, dFile, rewindfile )
  Str = trim(iStr); iError = iFindStringInFile ( Str, dFile, rewindfile )
  DO WHILE( rewindfile .EQ. 0 ) 
    IF( iStr(1:3) .EQ. '*Ti') READ( dFile,* ) nci
    READ( dFile,* ) dataIn
    
    WRITE( str,'(i3.3)' ) itime
    itime = itime + 1
    
    Path = 'solution_VolumeData.vtk'
    WRITE(*,*) 'Writing ', itime, '/', ntotal 
  
    !!!!!! 3D Volume mesh
    OPEN( oFile, FILE = trim(Path) )
    WRITE( oFile, '(a,/,a)' ) '# vtk DataFile Version 4.2', 'TEPEM approximation'
    WRITE( oFile, '(a,/,a)' ) 'ASCII', 'DATASET UNSTRUCTURED_GRID'
    WRITE( oFile, '(a,i8,a)' ) 'POINTS ', noPoints*noPipe, ' float' 
    OPEN( iFile, FILE = 'Param00.txt' )
      Str = '*Parameter'; iError = iFindStringInFile ( Str, iFile, rewindfile )
      READ( iFile,* ) nck
      
      Str = '*Real'; iError = iFindStringInFile ( Str, iFile, rewindfile )
      READ( iFile,* ) ( nci, nc0 = 1, nck )
      DO nc0 = 1, noPipe
        READ( iFile,* ) SerPoints( 1:iPar_Real(nc0)-1 ), det
        
        DO nci = 1, noPoints
          x0 = xlocal(nci,:)
          CALL EvaluateMap( x0, SerPoints, x0, det, (iPar_Real(nc0)-1)/3 )
          WRITE( oFile, '(3(f16.8))' ) x0
        ENDDO
      ENDDO
    CLOSE( iFile )

    !!Geometry
    WRITE( oFile, '(/,a,i10,i10)' ) 'CELLS ', NoCells, 9*NoCells
    DO nc0 = 1, noPipe
    DO nck = 1, par_z - 1 
    DO ncj = 1, par_n - 1
    DO nci = 1, par_n - 1
      nod = ( nc0 - 1 ) * par_z * par_n ** 2
      WRITE( oFile, '(i2)', advance = 'no') 8
      WRITE( oFile,'(i10)', advance = 'no') nod + (nck-1)*par_n**2 + nci + (ncj - 1)*par_n - 1
      WRITE( oFile,'(i10)', advance = 'no') nod + (nck-1)*par_n**2 + nci + 1 + (ncj - 1)*par_n - 1
      WRITE( oFile,'(i10)', advance = 'no') nod + (nck-1)*par_n**2 + par_n + nci + 1 + (ncj - 1)*par_n - 1
      WRITE( oFile,'(i10)', advance = 'no') nod + (nck-1)*par_n**2 + par_n + nci + (ncj - 1)*par_n - 1
      WRITE( oFile,'(i10)', advance = 'no') nod + nck*par_n**2 + nci + (ncj - 1)*par_n - 1
      WRITE( oFile,'(i10)', advance = 'no') nod + nck*par_n**2 + nci + 1 + (ncj - 1)*par_n - 1
      WRITE( oFile,'(i10)', advance = 'no') nod + nck*par_n**2 + par_n + nci + 1 + (ncj - 1)*par_n - 1
      WRITE( oFile,'(i10)') nod + nck*par_n**2 + par_n + nci + (ncj - 1)*par_n - 1
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    WRITE( oFile, '(/,a,i10)' ) 'CELL_TYPES ', NoCells
    WRITE( oFile, '(i4)' ) ( 12, nci = 1, noCells )
    
    OPEN( iFile, FILE = 'Mesh.txt' )
      Str = '*INCIDENCE'; iError = iFindStringInFile ( Str, iFile, rewindfile )
      !nod = 3*(pp+1)**2
      nci = 0
      DO nc0 = 1, noPipe
        nod = nGroup( nc0 )
        READ( iFile,* ) iInc( 1:nod )
        DO nck = 1, nod
          ncj = iInc( nck )
          Sol1((nck-1)*idoft+1:nck*idoft) = dataIn((ncj-1)*idoft+1:ncj*idoft)
        ENDDO

        iInc(1:3) = pType( nType( nc0 ), 1:3)
        DO nck = 1, noPoints
          x0 = xlocal(nck,:)
          CALL EvaluateLocalVectorField( x0, Sol1, nod, iInc(1:3), idoft, u0, 0 )
          nci = nci + 1
          dataOut((nci-1)*idoft+1:nci*idoft) = u0
        ENDDO
      ENDDO
    CLOSE( iFile )
    
    WRITE( oFile, '(/,a,i10)' ) 'POINT_DATA ', noPoints*noPipe
    WRITE( oFile, '(a)' ) 'SCALARS Pressure float 1'
    WRITE( oFile, '(a)' ) 'LOOKUP_TABLE default'
    WRITE( oFile, '(e16.8)') ( dataOut( ( nci - 1 )*idoft + 4 ), nci = 1, noPoints*noPipe )
    
    WRITE( oFile, '(/,a)' ) 'FIELD FieldData 1'
    WRITE( oFile, '(a,i8,a)') 'Velocity 3 ', noPoints*noPipe, ' float'
    WRITE( oFile, '(e16.8)') ( dataOut( ( nci - 1 )*idoft + 1:( nci - 1 )*idoft + 3 ), nci = 1, noPoints*noPipe )
    
    WRITE( oFile, '(/,a)' ) 'FIELD FieldData 1'
    WRITE( oFile, '(a,i8,a)') 'WSS 3 ', noPoints*noPipe, ' float'
    WRITE( oFile, '(e16.8)') ( dataOut( ( nci - 1 )*idoft + 5:( nci - 1 )*idoft + 7 ), nci = 1, noPoints*noPipe )
    CLOSE( oFile )
    
    !Str = '*Time'; iError = iFindStringInFile ( Str, dFile, rewindfile )
    Str = trim(iStr); iError = iFindStringInFile ( Str, dFile, rewindfile )
   ENDDO
 CLOSE( dFile )
 
 DEALLOCATE( iPar_Real, nType, nGroup )
 DEALLOCATE( dataIn, dataOut, xlocal ) 
END

 
!     ------------------------------------------------------------------
Integer Function iFindStringInFile (Str, ioUnit, rewindfile )
!     ------------------------------------------------------------------
!     Busca un String en un archivo (Str), si no lo encuentra rebovina 
!     el archivo y lo busca nuevamente.
!
!     Str: String to find, ioUnit: Unit assigned to Input File;
!     iError: return value, 0: String Found, 1: string not found
!
      Character(*) Str, DummyString*120
      INTEGER :: rewindfile
!
      iError=0
      rewindfile = 0
      Leng = Len_Trim(Str)
!
      Do While (iError .eq. 0)
         Read (ioUnit, "(A120)", iostat = iError) DummyString
         If (iError .ne. 0) Exit
         If (DummyString(1:1) .ne. Str(1:1)) Cycle
         If (DummyString(1:Leng) .eq. Str(1:Leng)) Exit
      End Do
!
      If (iError .ne. 0) Then   ! Still not found
         iError = 0
         rewind (ioUnit)
         rewindfile = 1
         Do While (iError .eq. 0)
            Read (ioUnit, "(A120)", iostat=iError) DummyString
!
            If (iError .ne. 0) Exit
            If (DummyString(1:1) .ne. Str(1:1)) Cycle
            If (DummyString(1:Leng) .eq. Str(1:Leng)) Exit
         End Do
      End If
!
      iFindStringInFile = iError
      Return
      End Function    
