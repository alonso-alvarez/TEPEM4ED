!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EvaluateMap( coor, SerPoints, map, det, nodG )
  IMPLICIT NONE
  INTEGER, PARAMETER :: Ndim = 3

  INTEGER :: n, i, j, k, id, NN, nodG
  DOUBLE PRECISION :: map( Ndim )

  DOUBLE PRECISION :: SerPoints( NodG * Ndim )
  DOUBLE PRECISION :: x, y, z, det
  DOUBLE PRECISION :: coor( Ndim ), grad( Ndim, Ndim )
  DOUBLE PRECISION :: Geo( NodG ), dGeo( Ndim, NodG )

  CALL Pipe_GeoBasis3D( nodG/3, coor, Geo, dGeo )

  map  = 0.0d0; grad = 0.d0
  DO n = 1, NodG
  DO i = 1, Ndim
    map(i) = map(i) + SerPoints( (n-1)*Ndim + i )*Geo( n )
    DO j = 1, Ndim
           grad( i, j ) = grad( i, j ) + SerPoints( (n-1)*Ndim + i )*dGeo( j, n )
    ENDDO
  ENDDO
  ENDDO
  det = grad(1,1)*grad(2,2)*grad(3,3) + grad(2,1)*grad(3,2)*grad(1,3) + grad(1,2)*grad(3,1)*grad(2,3)
  det = det - grad(2,2)*grad(3,1)*grad(1,3) - grad(3,2)*grad(1,1)*grad(2,3) - grad(1,2)*grad(2,1)*grad(3,3) 
END SUBROUTINE 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Pipe_GeoBasis3D( nodG, x0, Geo, dGeo )
  IMPLICIT NONE
  INTEGER, PARAMETER :: Ndim = 3
  INTEGER, INTENT( IN ) :: nodG
 
  INTEGER :: nci, ncj, nck, i
  DOUBLE PRECISION :: x, y, z
  DOUBLE PRECISION :: x0( ndim ),  Geo( 3*nodG ), dGeo( ndim,3*nodG )
  DOUBLE PRECISION :: Fi( nodG ), dFi( ndim,nodG )

  x = x0( 1 ); y = x0( 2 ); z = x0( 3 )
  CALL Pipe_GeoBasis2D( nodG, x0, Fi, dFi)

  Geo  = 0.0d0; dGeo = 0.0d0
  DO i = 1, nodG
       Geo( 0*nodG + i ) = Fi( i ) * 0.5 * z * ( z - 1 )
       Geo( 1*nodG + i ) = Fi( i ) * ( 1. - z**2 )
       Geo( 2*nodG + i ) = Fi( i ) * 0.5 * z * ( z + 1 )
       
       dGeo( 1, 0*nodG + i ) = dFi( 1, i ) * 0.5 * z * ( z - 1 )
       dGeo( 2, 0*nodG + i ) = dFi( 2, i ) * 0.5 * z * ( z - 1 )
       dGeo( 3, 0*nodG + i ) = Fi( i ) * 0.5 * ( 2*z - 1 )
       
       dGeo( 1, 1*nodG + i ) = dFi( 1, i ) * ( 1. - z**2 )
       dGeo( 2, 1*nodG + i ) = dFi( 2, i ) * ( 1. - z**2 )
       dGeo( 3, 1*nodG + i ) = Fi( i ) * ( -2*z )
       
       dGeo( 1, 2*nodG + i ) = dFi( 1, i ) * 0.5 * z * ( z + 1 )
       dGeo( 2, 2*nodG + i ) = dFi( 2, i ) * 0.5 * z * ( z + 1 )
       dGeo( 3, 2*nodG + i ) = Fi( i ) * 0.5 * ( 2*z + 1 )
  ENDDO 
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Pipe_GeoBasis2D( nodG, x0, Fi, dFi )
  IMPLICIT NONE
  INTEGER, PARAMETER :: Ndim = 3
  INTEGER, INTENT( IN ) :: nodG

  DOUBLE PRECISION :: x, y, z, x0( ndim )
  DOUBLE PRECISION :: Fi( nodG ), dFi( ndim,nodG )

  x = x0( 1 ); y = x0( 2 ); z = x0( 3 )
  IF( nodG .EQ. 4 ) THEN
    Fi( 1 ) = 0.25*( 1 - x )*( 1 - y ); dFi(1,1) = -0.25*( 1 - y ); dFi(2,1) = -0.25*( 1 - x )
    Fi( 2 ) = 0.25*( 1 + x )*( 1 - y ); dFi(1,2) =  0.25*( 1 - y ); dFi(2,2) = -0.25*( 1 + x )
    Fi( 3 ) = 0.25*( 1 + x )*( 1 + y ); dFi(1,3) =  0.25*( 1 + y ); dFi(2,3) =  0.25*( 1 + x )
    Fi( 4 ) = 0.25*( 1 - x )*( 1 + y ); dFi(1,4) = -0.25*( 1 + y ); dFi(2,4) =  0.25*( 1 - x )
  ENDIF
  IF( nodG .EQ. 5 ) THEN
    Fi( 1 ) = 0.25*( 1 - x )*( 1 - y ); 
    Fi( 2 ) = 0.25*( 1 + x )*( y**2 - y )
    Fi( 3 ) = 0.25*( 1 + x )*( y**2 + y )
    Fi( 4 ) = 0.25*( 1 - x )*( 1 + y )
    Fi( 5 ) = 0.50*( 1 + x )*( 1 - y**2 )

    dFi(1,1) = -0.25*( 1 - y );   dFi(2,1) = -0.25*( 1 - x )
    dFi(1,2) = 0.25*( y**2 - y ); dFi(2,2) = 0.25*( 1 + x )*( 2*y - 1 )
    dFi(1,3) = 0.25*( y**2 + y ); dFi(2,3) = 0.25*( 1 + x )*( 2*y + 1 )
    dFi(1,4) = -0.25*( 1 + y );   dFi(2,4) = 0.25*( 1 - x )
    dFi(1,5) = 0.50*( 1 - y**2 ); dFi(2,5) = 0.50*( 1 + x )*( -2*y )
  ENDIF
  IF( nodG .EQ. 6 ) THEN
    Fi( 1 ) = 0.25*(x - 1)*(y - 1)
    Fi( 2 ) = -0.03125*(x + 1)*(3*y + 1)*(3*y - 1)*(y - 1)
    Fi( 3 ) =  0.03125*(x + 1)*(3*y + 1)*(3*y - 1)*(y + 1)
    Fi( 4 ) = -0.25*(x - 1)*(y + 1)
    Fi( 5 ) =  0.28125*(x + 1)*(3*y - 1)*(y + 1)*(y - 1)
    Fi( 6 ) = -0.28125*(x + 1)*(3*y + 1)*(y + 1)*(y - 1)  

    dFi(1,1) = 1/4.*y - 1/4.; dFi(2,1) = 1/4.*x - 1/4.
    dFi(1,2) = -1/32.*(3*y + 1)*(3*y - 1)*(y - 1);
    dFi(2,2) = -27/32.*x*y**2 + 9/16.*x*y - 27/32.*y**2 + 1/32.*x + 9/16.*y + 1/32.
    dFi(1,3) = 1/32.*(3*y + 1)*(3*y - 1)*(y + 1)
    dFi(2,3) = 27/32.*x*y**2 + 9/16.*x*y + 27/32.*y**2 - 1/32.*x + 9/16.*y - 1/32.
    dFi(1,4) = -1/4.*y - 1/4.; dFi(2,4) = -1/4.*x + 1/4.
    dFi(1,5) = 9/32.*(3*y - 1)*(y + 1)*(y - 1)
    dFi(2,5) = 81/32.*x*y**2 - 9/16.*x*y + 81/32.*y**2 - 27/32.*x - 9/16.*y - 27/32.
    dFi(1,6) = -9/32.*(3*y + 1)*(y + 1)*(y - 1)
    dFi(2,6) = -81/32.*x*y**2 - 9/16.*x*y - 81/32.*y**2 + 27/32.*x - 9/16.*y + 27/32.
  ENDIF
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EvaluatePlanarMap( coor, SerPoints, map, nodG )
  IMPLICIT NONE
  INTEGER, PARAMETER :: Ndim = 3

  INTEGER :: n, i, j, k, id, NN, nodG
  DOUBLE PRECISION :: map( Ndim )

  DOUBLE PRECISION :: SerPoints( NodG * Ndim )
  DOUBLE PRECISION :: x, y, z, det
  DOUBLE PRECISION :: coor( Ndim ), Geo( NodG ), dGeo( ndim,nodG )

  x = coor( 1 ); y = coor( 2 ); z = coor( 3 ) 
  CALL Pipe_GeoBasis2D( nodG, coor, Geo, dGeo )

  map  = 0.0d0
  DO n = 1, NodG
  DO i = 1, Ndim
    map(i) = map(i) + SerPoints( (n-1)*Ndim + i )*Geo( n )
  ENDDO
  ENDDO
END SUBROUTINE 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
SUBROUTINE CreateVoidSurface()
  IMPLICIT NONE
  INTEGER, PARAMETER :: Ndim = 3
  INTEGER :: iFile, oFile
  INTEGER :: iFindStringInFile, iError, rewindfile 

  INTEGER :: nci, ncj, nck, nc0, i, j, k, noCoor, i0
  INTEGER :: iPlane, iVolume, nPlane, nAxial, idet, cont
  DOUBLE PRECISION :: geonodes( 18*Ndim )

  CHARACTER*120 :: Str

  nPlane = 7; nAxial = 4
  iFile = 201; oFile = 202; iPlane = 0; iVolume = 0
  OPEN( iFile, FILE = 'Data_FluidElements.out' )
  OPEN( oFile, FILE = 'TEPEM_VoidSurface.vtk' )
  WRITE( oFile, '(a)' ) '# vtk DataFile Version 3.0'
  WRITE( oFile, '(a)' ) 'Volume Mesh'
  WRITE( oFile, '(a)' ) 'ASCII'
  WRITE( oFile, '(a)' ) 'DATASET UNSTRUCTURED_GRID'

  REWIND( iFile )
  Str = '*Begin'; iError = iFindStringInFile ( Str, iFile, rewindfile )   
  Str = '*PipeElement'; iError = iFindStringInFile ( Str, iFile, rewindfile ) 
  DO WHILE( rewindfile .EQ. 0 ) 
  	iVolume = iVolume + 1
  	
  	READ( iFile,* ) nci
  	IF( nci .EQ. 18 ) iPlane = iPlane + 1
    iError = iFindStringInFile ( Str, iFile, rewindfile )
  ENDDO

  noCoor = iPlane * Naxial * Nplane
  WRITE( oFile, '(a,i8,a)' ) 'POINTS ', noCoor, ' float'
  Str = '*Begin'; iError = iFindStringInFile ( Str, iFile, rewindfile )   
  DO nci = 1, iVolume    
    Str = '*PipeElement'; iError = iFindStringInFile ( Str, iFile, rewindfile )   

    READ( iFile, * ) i0
    READ( iFile, * ) GeoNodes( 1:i0*Ndim )
    IF( i0 .NE. 18 ) CYCLE
       
    CALL WriteLateralElement( oFile, GeoNodes, Naxial, Nplane )
  ENDDO
  CLOSE( iFile )
  WRITE( oFile, * ) ''
  nci =  iPlane*( Naxial - 1 )*( Nplane - 1 )
  WRITE( oFile, '(a, i12, i12)' ) 'CELLS ', nci, 5*nci
  DO cont = 1, iPlane
  DO k = 1, Naxial - 1
  DO i = 1, Nplane - 1
    nci = ( cont - 1 ) * Naxial * Nplane
    write( oFile,* ) 4, nci + (k-1)*Nplane + i - 1,  &
      nci + (k-1)*Nplane + i + 0,  nci + (k+0)*Nplane + i + 0,  &
      nci + (k+0)*Nplane + i - 1
  ENDDO
  ENDDO
  ENDDO

  WRITE( oFile, * ) ''
  nci = ( Nplane - 1 ) * ( Naxial - 1 ) * iPlane
  WRITE( oFile, '(a, i8)' ) 'CELL_TYPES ', nci
  DO i = 1, nci
    WRITE( oFile, *) 9
  ENDDO
  CLOSE( oFile ) 
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
SUBROUTINE WriteLateralElement( io, GeoNodes, NZ, NN )
  IMPLICIT NONE
  INTEGER, PARAMETER :: Ndim = 3
  INTEGER, INTENT( IN ) :: io

  INTEGER :: i, j, k, nci, NN, NZ, nodG
  DOUBLE PRECISION :: GeoNodes( 18*Ndim ), refcoor( ndim ), x0( ndim ), det

  DO k = 1, NZ
  DO j = 1, NN
    RefCoor( 1 ) =  1.
    RefCoor( 2 ) = -1. + 2.*( j - 1 )/( 1.*NN - 1 )
    RefCoor( 3 ) = -1. + 2.*( k - 1 )/( 1.*NZ - 1 )

    CALL EvaluateMap( RefCoor, GeoNodes, x0, det, 18 )
    WRITE( io, * ) x0
  ENDDO
  ENDDO 
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
SUBROUTINE WritePipeElement( io, nodG, GeoNodes, NZ, NN, idet, NodalDet, minD )
  IMPLICIT NONE
  INTEGER, PARAMETER :: Ndim = 3
  INTEGER, INTENT( IN ) :: io

  INTEGER :: i, j, k, nci
  INTEGER :: NN, NZ, nodG
  INTEGER :: idet, idet0

  DOUBLE PRECISION :: NodalDet, maxdet, mindet, minD
  DOUBLE PRECISION :: RefCoor( Ndim ), det, x0( Ndim ), det1, x1( Ndim )
  DOUBLE PRECISION :: GeoNodes( NodG*Ndim ), xi( 2 ), alpha, maxxi

  DIMENSION NodalDet( * ), minD(*)

  idet0 = idet
  maxdet = 0.d0
  mindet = 1.e+8
  DO k = 1, NZ
  DO j = 1, NN
  DO i = 1, NN
    RefCoor( 1 ) = ( -1.*NN + 2.*i - 1. )/( 1.*NN - 1 )
    RefCoor( 2 ) = ( -1.*NN + 2.*j - 1. )/( 1.*NN - 1 )
    RefCoor( 3 ) = ( -1.*NZ + 2.*k - 1. )/( 1.*NZ - 1 )

    CALL EvaluateMap( RefCoor, GeoNodes, x0, det, nodG )
    WRITE( io, * ) x0

    idet = idet + 1
    NodalDet( idet ) = det
    IF( det .GT. maxdet ) maxdet = det
    IF( det .LT. mindet ) mindet = det
  ENDDO
  ENDDO
  ENDDO 

  IF( maxdet .EQ. 0.0d0 ) maxdet = 1.e+9
  DO i = 1, NN**2 * NZ
    idet0 = idet0 + 1
    nodaldet(idet0) = nodaldet(idet0)/maxdet
    minD(idet0) = mindet
  ENDDO
END SUBROUTINE
