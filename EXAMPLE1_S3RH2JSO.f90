!********************************************************************************************   
! example 1: coordinate transformation from S3RH to JSO 
!********************************************************************************************                   
!   VARIABLE        INPUT/OUTPUT         DESCRIPTION 
!   --------        ------------         ----------------------------------------------------
!    time             INPUT              UTC, Universal Coordinate Time
!    IVECTOR          INPUT              Initial vector
!    tdb                                 TDB, Barycentric Dynamical Time 
!    X_S3RH                              Unit vector of S3RH X-axis expressed in JCRS
!    Y_S3RH                              Unit vector of S3RH Y-axis expressed in JCRS
!    Z_S3RH                              Unit vector of S3RH Z-axis expressed in JCRS
!    X_JSO                               Unit vector of JSO X-axis expressed in JCRS
!    Y_JSO                               Unit vector of JSO Y-axis expressed in JCRS
!    Z_JSO                               Unit vector of JSO Z-axis expressed in JCRS
!    TM                                  Transformation matrix from S3RH to JSO
!    FVECTOR          OUTPUT             Final vector
!
!In this example:
!    time=2021-027 00:00:31.478
!
!    IVECTOR=[-4001865.916,-5946601.666,-3637712.347]
!
!When this code is executed, the output (include some intermediate results) is:
!
!           Transformation from S3RH to JSO coordinate system
!  
!           Time(Julian date,TDB)
!             2459241.501165077
!  
!           Initial vector
!                     -4001865.916          -5946601.666          -3637712.347
!  
!           Transformation matrix
!             9.7956274636932655E-002    0.99515956648484361        -7.8743566690941402E-003
!            -0.99370632346307930        9.8239103991244317E-002    5.3822125157663492E-002
!             5.4335172482872729E-002    2.5525831417951284E-003    0.99851949072141832     
!  
!           Final vector
!                    -6281180.768           3196701.238          -3864947.950
!*******************************************************************************************
PROGRAM EXAMPLE1

USE MOD_PRECISION

USE MOD_TOP2013

USE MOD_TIME_TRANSFORMATION

USE MOD_CS_TRANSFORMATION

IMPLICIT NONE

!variables
REAL(KIND=dp_ud):: tdb,X_S3RH(3,1),Y_S3RH(3,1),Z_S3RH(3,1),&
                   X_JSO(3,1),Y_JSO(3,1),Z_JSO(3,1), &
                   IVECTOR(3,1),FVECTOR(3,1),TM(3,3)
                   
OPEN(8,FILE='example1.txt',STATUS='new')
 
! convert UTC to JDB
! time format: year, month, day, hour, minute, second
CALL WR_UTC2TDB(2021,1,27,0,0,31.478d0,tdb)

! convert UTC to JDB
! time format: year, day of year, hour, minute, second
!CALL WR_UTC2TDB2(2021,27,0,0,31.478d0,tdb)

!Initial position vector
IVECTOR=RESHAPE((/-4001865.916D0,-5946601.666D0,-3637712.347D0/),(/3,1/))

! unit vectors of S3RH coordinate axes expressed in JCRS
CALL S3RH_JCRS(tdb,X_S3RH,Y_S3RH,Z_S3RH)

! unit vectors of JSO coordinate axes expressed in JCRS
CALL JSO_JCRS(tdb,X_JSO,Y_JSO,Z_JSO)

!Transformation matrix from S3RH to JSO
CALL TRANSFMATRIX(X_S3RH,Y_S3RH,Z_S3RH,X_JSO,Y_JSO,Z_JSO,TM)

!Final vector
FVECTOR=MATMUL(TM,IVECTOR)

!Output
WRITE(8,*) 'Transformation from S3RH to JSO coordinate system'
WRITE(8,*) ' '
WRITE(8,*) 'Time(Julian date, TDB):'
WRITE(8,"(F20.9)") tdb
WRITE(8,*) ' '
WRITE(8,*) 'Initial vector:'
WRITE(8,"(3(2X,F20.3))") IVECTOR
WRITE(8,*) ' '
WRITE(8,*) 'Transformation matrix:'
WRITE(8,*) TM(1,1:3)
WRITE(8,*) TM(2,1:3)
WRITE(8,*) TM(3,1:3)
WRITE(8,*) ' '
WRITE(8,*) 'Final vector:'
WRITE(8,"(3(2X,F20.3))") FVECTOR


END PROGRAM EXAMPLE1