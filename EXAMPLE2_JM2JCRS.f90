!********************************************************************************************   
! example 2: coordinate transformation from JM to JCRS 
!********************************************************************************************                   
!   VARIABLE        INPUT/OUTPUT         DESCRIPTION 
!   --------        ------------         ----------------------------------------------------
!    time             INPUT              UTC, Universal Coordinate Time
!    IVECTOR          INPUT              Initial vector
!    tdb                                 TDB, Barycentric Dynamical Time 
!    X_JM                                Unit vector of JM X-axis expressed in JCRS
!    Y_JM                                Unit vector of JM Y-axis expressed in JCRS
!    Z_JM                                Unit vector of JM Z-axis expressed in JCRS
!    TM                                  Transformation matrix from JM to JCRS
!    FVECTOR          OUTPUT             Final vector
!
!
!In this example:
!    time=2021-027 00:00:31.478
!
!    IVECTOR=[-4001865.916,-5946601.666,-3637712.347]
!
!When this code is executed, the output (include some intermediate results) is:
!
!           Transformation from JM to JCRS coordinate system
!  
!           Time(Julian date,TDB):
!             2459241.501165077
!  
!           Initial vector:
!                    -4001865.916          -5946601.666          -3637712.347
!  
!           Transformation matrix:
!           -0.86490514103363658       0.47182443557849807      -0.17123317144939215     
!           -0.33712242623282507      -0.79880654947698582      -0.49825351603733331     
!           -0.37187036282080704      -0.37321548534821469       0.84995448981214150     
!
!           Final vector:
!                     1278379.452           7911806.072            615649.215
!
!
!*******************************************************************************************
PROGRAM EXAMPLE2

USE MOD_PRECISION

USE MOD_TOP2013

USE MOD_TIME_TRANSFORMATION

USE MOD_CS_TRANSFORMATION

IMPLICIT NONE

!variables
REAL(KIND=dp_ud):: tdb,X_JM(3,1),Y_JM(3,1),Z_JM(3,1),&
                   IVECTOR(3,1),FVECTOR(3,1),TM(3,3),M(3,3)
                   
OPEN(9,FILE='example2.txt',STATUS='new')

! convert UTC to JDB
! time format: year, month, day, hour, minute, second
CALL WR_UTC2TDB(2021,1,27,0,0,31.478d0,tdb)

! convert UTC to JDB
! time format: year, day of year, hour, minute, second
!CALL WR_UTC2TDB2(2021,27,0,0,31.478d0,tdb)

!Initial position vector
IVECTOR=RESHAPE((/-4001865.916D0,-5946601.666D0,-3637712.347D0/),(/3,1/))

! unit vectors of JM coordinate axes expressed in JCRS
! three Jovian magnetic field models:VIP4,JRM09,JRM33
CALL JM_JCRS(tdb,'JRM33',X_JM,Y_JM,Z_JM)

!Transformation matrix from JM to JCRS
TM=RESHAPE((/X_JM(1,1),X_JM(2,1),X_JM(3,1),Y_JM(1,1),Y_JM(2,1),Y_JM(3,1),&
             Z_JM(1,1),Z_JM(2,1),Z_JM(3,1)/),(/3,3/))
             
!Final vector
FVECTOR=MATMUL(TM,IVECTOR)

!Output
WRITE(9,*) 'Transformation from JM to JCRS coordinate system'
WRITE(9,*) ' '
WRITE(9,*) 'Time(Julian date, TDB):'
WRITE(9,"(F20.9)") tdb
WRITE(9,*) ' '
WRITE(9,*) 'Initial vector:'
WRITE(9,"(3(2X,F20.3))") IVECTOR
WRITE(9,*) ' '
WRITE(9,*) 'Transformation matrix:'
WRITE(9,*) TM(1,1:3)
WRITE(9,*) TM(2,1:3)
WRITE(9,*) TM(3,1:3)
WRITE(9,*) ' '
WRITE(9,*) 'Final vector:'
WRITE(9,"(3(2X,F20.3))") FVECTOR


END PROGRAM EXAMPLE2