! MOD_CS_TRANSFORMATION.f90
! MODULE OF COORDINATE SYSTEMS' TRANSFORMATION IS CREATED TO TRANSFORM BETWEEN
! ANY TWO TARGET COORDINATE SYSTEMS.

! THIS MODULE CALCULATES THE VALUES OF ALL JOVIAN COORDINATE AXES IN JCRS
! AND THE FINAL TRANSFORMATION MATRIX BETWEEN ANY TWO COORDINATE SYSTEMS

!
!CREATED BY YUJIE WANG ON 2023/04/02
!REVISED BY YABING WANG ON 2023/04/20
!REVISED BY YUJIE WANG ON 2023/06/25

MODULE MOD_CS_TRANSFORMATION
    
USE MOD_PRECISION

USE MOD_TIME_TRANSFORMATION

USE MOD_TOP2013

IMPLICIT NONE

PRIVATE

!********************************************************************************

PUBLIC  ROTATION_X

PUBLIC  ROTATION_Y

PUBLIC  ROTATION_Z

PUBLIC  CROSS_PRODUCT

PUBLIC  NORM

PUBLIC  ROTATIONAL_ELEMENTS

PUBLIC  TRANSFMATRIX

PUBLIC  S3RH_JCRS

PUBLIC  JEI_JCRS

PUBLIC  JSE_JCRS

PUBLIC  JSO_JCRS

PUBLIC  JH_JCRS

PUBLIC  JSM_JCRS

PUBLIC  JM_JCRS

PUBLIC  JSMAG_JCRS

!********************************************************************************
CONTAINS

!==================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
!==================================================================================
SUBROUTINE ROTATION_X (ANGLE,MATRIX)

! Calculate the rotation matrix about the X-axis
!------------------------------------------------------------------------------
! INPUT:
!       angle     unit: degree
! OUTPUT:
!       matrix    a 3*3 transformation matrix
!------------------------------------------------------------------------------
IMPLICIT NONE

REAL(KIND = dp_ud) ::  ANGLE, MATRIX(3,3)

    MATRIX = RESHAPE((/1.D0, 0.D0, 0.D0, 0.D0, DCOS(ANGLE*DEG2RAD), &
                    -DSIN(ANGLE*DEG2RAD), 0.D0, DSIN(ANGLE*DEG2RAD), &
                    DCOS(ANGLE*DEG2RAD)/), (/3,3/))

    RETURN
END SUBROUTINE ROTATION_X
!==================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
!==================================================================================
SUBROUTINE ROTATION_Y (ANGLE,MATRIX)

! Calculate the rotation matrix about the Y-axis
!------------------------------------------------------------------------------
IMPLICIT NONE

REAL(KIND = dp_ud) ::  ANGLE, MATRIX(3,3)

    MATRIX = RESHAPE((/DCOS(ANGLE*pi/180), 0.D0, DSIN(ANGLE*DEG2RAD), &
                      0.D0, 1.D0, 0.D0, -DSIN(ANGLE*DEG2RAD), 0.D0, &
                     DCOS(ANGLE*DEG2RAD)/), (/3,3/))

    RETURN
END SUBROUTINE ROTATION_Y

!==================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
!==================================================================================
SUBROUTINE ROTATION_Z (ANGLE,MATRIX)

! Calculate the rotation matrix about the Z-axis
!------------------------------------------------------------------------------
IMPLICIT NONE

REAL(KIND = dp_ud) ::  ANGLE,MATRIX(3,3)

    MATRIX = RESHAPE((/DCOS(ANGLE*DEG2RAD), -DSIN(ANGLE*DEG2RAD), &
                      0.D0, DSIN(ANGLE*DEG2RAD), DCOS(ANGLE*DEG2RAD), &
                      0.D0, 0.D0, 0.D0, 1.D0/), (/3,3/))

    RETURN
END SUBROUTINE ROTATION_Z

!==================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
!==================================================================================
SUBROUTINE  CROSS_PRODUCT(M1,M2,M)

! Vector cross product
!------------------------------------------------------------------------------
IMPLICIT NONE

REAL(KIND = dp_ud) ::  M1(3,1), M2(3,1), M(3,1)

    M(1,1) = M1(2,1) * M2(3,1) - M1(3,1) * M2(2,1)

    M(2,1) = M1(3,1) * M2(1,1) - M1(1,1) * M2(3,1)

    M(3,1) = M1(1,1) * M2(2,1) - M1(2,1) * M2(1,1)

    RETURN
END SUBROUTINE  CROSS_PRODUCT

!==================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
!==================================================================================
SUBROUTINE NORM(VECTOR,N)

! Magnitude of a vector
!------------------------------------------------------------------------------
IMPLICIT NONE

REAL(KIND = dp_ud) ::  VECTOR(3,1), N

    N = DSQRT( VECTOR(1,1)**2 + VECTOR(2,1)**2 + VECTOR(3,1)**2 )

    RETURN
END SUBROUTINE NORM

!==================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
!==================================================================================
SUBROUTINE  ROTATIONAL_ELEMENTS (tdb, RA, DE, W)
!------------------------------------------------------------------------------
! Rotational elements are from Archinal (2018)
!
! INPUT:
!       tdb   time, julian date in TDB
!
! OUTPUT:
!       RA      right ascension of Jupiter at time tdb
!       DE      declination of Jupiter at time tdb
!       W       angle related to rotation of Jupiter at time tdb
!------------------------------------------------------------------------------
IMPLICIT NONE

REAL(KIND = dp_ud) ::  T, A, B, C, D, E, RA, DE, W, tdb

    T = ( tdb - 2451545.D0 )/36525.D0

    A = ( 99.360714D0 + 4850.4046D0*T )*DEG2RAD

    B = ( 175.895369D0 + 1191.9605D0*T )*DEG2RAD

    C = ( 300.323162D0 + 262.5475D0*T )*DEG2RAD

    D = ( 114.012305D0 + 6070.2476D0*T )*DEG2RAD

    E = ( 49.511251D0 + 64.3D0*T )*DEG2RAD

    RA = 268.056595D0 - 0.006499D0*T + 0.000117D0*DSIN(A) + &
         0.000938D0*DSIN(B) + 0.001432D0*DSIN(C) + &
         0.00003D0*DSIN(D) + 0.00215D0*DSIN(E)

    DE = 64.495303D0 + 0.002413D0*T + 0.00005D0*DCOS(A) + &
         0.000404D0*DCOS(B) + 0.000617D0*DCOS(C) - &
         0.000013D0*DCOS(D) + 0.000926D0*DCOS(E)

    W = 284.95D0 + 870.536D0*T*36525.D0

    RETURN
END SUBROUTINE ROTATIONAL_ELEMENTS

!==================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
!==================================================================================
SUBROUTINE  TRANSFMATRIX (X_INI,Y_INI,Z_INI,X_FIN,Y_FIN,Z_FIN,TM)
!------------------------------------------------------------------------------
! Calculate the transformation matrix between two coordinate systems
!
! INPUT:
!       X_INI   X-axis REPRESENTATION OF THE INITIAL COORDINATE SYSTEM in JCRS
!       Y_INI   Y-axis ...
!       Z_INI   Z-axis ...
!
!       X_FIN   X-axis REPRESENTATION OF THE FINAL COORDINATE SYSTEM in JCRS
!       Y_FIN   Y-axis ...
!       Z_FIN   Z-axis ...
!
! OUTPUT:
!
!       TM      transformation matrix of two coordinate systems
!------------------------------------------------------------------------------
IMPLICIT NONE

REAL(KIND = dp_ud) ::  X_INI(3),Y_INI(3),Z_INI(3), X_FIN(3),Y_FIN(3),Z_FIN(3),TM(3,3)

    TM = RESHAPE((/DOT_PRODUCT(X_INI,X_FIN),DOT_PRODUCT(X_INI,Y_FIN),DOT_PRODUCT(X_INI,Z_FIN),&
            DOT_PRODUCT(Y_INI,X_FIN),DOT_PRODUCT(Y_INI,Y_FIN),DOT_PRODUCT(Y_INI,Z_FIN),&
            DOT_PRODUCT(Z_INI,X_FIN),DOT_PRODUCT(Z_INI,Y_FIN),DOT_PRODUCT(Z_INI,Z_FIN)/),(/3,3/))

    RETURN
END SUBROUTINE  TRANSFMATRIX

!==================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
!==================================================================================
SUBROUTINE  S3RH_JCRS(tdb, X_S3RH, Y_S3RH, Z_S3RH)
!------------------------------------------------------------------------------
! Calculate the expressions of System III (right hand) coordinate axes in JCRS
!
! INPUT:
!       tdb        time, julian date in TDB
! OUTPUT:
!       X_S3RH    expression of S3RH x-aixs in JCRS 3*1 matrix
!       Y_S3RH    expression of S3RH y-aixs in JCRS ...
!       Z_S3RH    expression of S3RH z-aixs in JCRS ...
!------------------------------------------------------------------------------
IMPLICIT NONE

REAL(KIND = dp_ud) ::  RA, DE, W, tdb, N

REAL(KIND = dp_ud) ::  X_S3RH(3,1),Y_S3RH(3,1),Z_S3RH(3,1),M1(3,3),M2(3,3),M3(3,3),Y_S3RH0(3,1)

!calculate rotational elements of Jupiter at tdb

    CALL ROTATIONAL_ELEMENTS(tdb,RA,DE,W)

!transform S3RH axes to JCRS: three angles

    CALL ROTATION_Z(-(90.D0+RA),M1)
    CALL ROTATION_X(-(90.D0-DE),M2)
    CALL ROTATION_Z(-W,M3)

    X_S3RH = MATMUL(MATMUL(M1,MATMUL(M2,M3)),RESHAPE((/1,0,0/),(/3,1/)))

    Z_S3RH = MATMUL(MATMUL(M1,MATMUL(M2,M3)),RESHAPE((/0,0,1/),(/3,1/)))

    CALL CROSS_PRODUCT(Z_S3RH,X_S3RH,Y_S3RH0)

    CALL NORM(Y_S3RH0,N)

    Y_S3RH = Y_S3RH0/N

    RETURN
END SUBROUTINE  S3RH_JCRS

!==================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
!==================================================================================
SUBROUTINE  JEI_JCRS(X_JEI,Y_JEI,Z_JEI)
!------------------------------------------------------------------------------
! Calculate the expressions of JEI coordinate axes in JCRS
!
! OUTPUT:
!       X_JEI    expression of JEI x-aixs in JCRS FOR J2000
!       Y_JEI    expression of JEI y-aixs in JCRS ...
!       Z_JEI    expression of JEI z-aixs in JCRS ...
!------------------------------------------------------------------------------
IMPLICIT NONE

REAL(KIND = dp_ud) ::  RA, DE, W, tdb, N

REAL(KIND = dp_ud) ::  X_JEI(3,1),Y_JEI(3,1),Z_JEI(3,1),M1(3,3),M2(3,3),M3(3,3),Y_JEI0(3,1)

    tdb = 2451545.D0 ! January/1/2000, 12h TDB

!calculate rotational elements of Jupiter at tdb

    CALL ROTATIONAL_ELEMENTS(tdb,RA,DE,W)

!transform JEI axes to JCRS: two angles (no W)

    CALL ROTATION_Z(-(90.D0+RA),M1)
    CALL ROTATION_X(-(90.D0-DE),M2)
    CALL ROTATION_Z(-W,M3)

    X_JEI=MATMUL(MATMUL(M1,MATMUL(M2,M3)),RESHAPE((/1,0,0/),(/3,1/)))

    Z_JEI=MATMUL(MATMUL(M1,MATMUL(M2,M3)),RESHAPE((/0,0,1/),(/3,1/)))

    CALL CROSS_PRODUCT(Z_JEI,X_JEI,Y_JEI0)

    CALL NORM(Y_JEI0,N)

    Y_JEI=Y_JEI0/N

    RETURN
END SUBROUTINE  JEI_JCRS

!==================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
!==================================================================================
SUBROUTINE  JSE_JCRS(tdb,X_JSE,Y_JSE,Z_JSE)
!------------------------------------------------------------------------------
! Calculate the expressions of JSE coordinate axes in JCRS
!
! INPUT:
!       tdb        time
! OUTPUT:
!       X_JSE    expression of JSE x-aixs in JCRS
!       Y_JSE    expression of JSE y-aixs in JCRS
!       Z_JSE    expression of JSE z-aixs in JCRS
!------------------------------------------------------------------------------
IMPLICIT NONE

REAL(KIND = dp_ud) ::  tdb,PV(6,1),R(3,1),RR,N1,N2,X_JSE0(3,1),Y_JSE0(3,1)

REAL(KIND = dp_ud) ::  X_JSE(3,1),Y_JSE(3,1),Z_JSE(3,1),JS(3,1),X_S3RH(3,1),Y_S3RH(3,1),Z_S3RH(3,1)

!z-axis of JSE is same with z-axis of S3RH

    CALL S3RH_JCRS(tdb,X_S3RH,Y_S3RH,Z_S3RH)

    Z_JSE=Z_S3RH

!position and velocity of JOVIAN BARYCENTER in heliocentric coordinate with ICRF (TOP2013)

    CALL JUPITERPV(tdb,PV)

    R = RESHAPE((/PV(1,1),PV(2,1),PV(3,1)/),(/3,1/))

    CALL NORM(R,RR)

!y-axis of JSE is determined by Jupiter-Sun vector and Jupiter spin axis

    JS = -R/RR

    CALL CROSS_PRODUCT(Z_JSE,JS,Y_JSE0)

    CALL NORM(Y_JSE0,N1)

    Y_JSE = Y_JSE0/N1

    CALL CROSS_PRODUCT(Y_JSE,Z_JSE,X_JSE0)

    CALL NORM(X_JSE0,N2)

    X_JSE = X_JSE0/N2

    RETURN
END SUBROUTINE  JSE_JCRS

!==================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
!==================================================================================
SUBROUTINE  JSO_JCRS(tdb, X_JSO, Y_JSO, Z_JSO)
!------------------------------------------------------------------------------
! Calculate the expressions of JSO coordinate axes in JCRS

! OUTPUT:
!       X_JSO    expression of JSO x-aixs in JCRS
!       Y_JSO    expression of JSO y-aixs in JCRS
!       Z_JSO    expression of JSO z-aixs in JCRS
!------------------------------------------------------------------------------
IMPLICIT NONE

REAL(KIND = dp_ud)::  tdb, PV(6,1), R(3,1), V(3,1), RR, VV, V_UNIT(3,1), &
                      X_JSO(3,1), Y_JSO(3,1), Z_JSO(3,1), Y_JSO0(3,1),&
                      Z_JSO0(3,1), N1, N2

!position and velocity of JOVIAN BARYCENTER in heliocentric coordinate with ICRF (TOP2013)

    CALL JUPITERPV(tdb,PV)

    R = RESHAPE((/PV(1,1),PV(2,1),PV(3,1)/),(/3,1/))

    V = RESHAPE((/PV(4,1),PV(5,1),PV(6,1)/),(/3,1/))

!x-axis of JSO is along with Jupiter-Sun vector
    CALL NORM(R,RR)

    X_JSO = -R/RR

!z-axis of JSO is perpendicular to the orbital plane determind by JOVAIN center's position and velocity vectors

    CALL NORM(V,VV)

    V_UNIT = V/VV

    CALL CROSS_PRODUCT(-X_JSO,V_UNIT,Z_JSO0)

    CALL NORM(Z_JSO0,N1)

    Z_JSO = Z_JSO0/N1

    CALL CROSS_PRODUCT(Z_JSO,X_JSO,Y_JSO0)

    CALL NORM(Y_JSO0,N2)

    Y_JSO = Y_JSO0/N2

    RETURN
END SUBROUTINE  JSO_JCRS

!==================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
!==================================================================================
SUBROUTINE JH_JCRS(tdb,X_JH,Y_JH,Z_JH)
!------------------------------------------------------------------------------
! Calculate the expressions of JH coordinate axes in JCRS

! OUTPUT:
!       X_JH    expression of JH x-aixs in JCRS
!       Y_JH    expression of JH y-aixs in JCRS
!       Z_JH    expression of JH z-aixs in JCRS
!------------------------------------------------------------------------------

IMPLICIT NONE

REAL(KIND = dp_UD):: tdb, R(3,1), RR, RA, DE, W, N1, N2, PV(6,1)

REAL(KIND = dp_ud)::  Z_SUN(3,1), X_JH(3,1), Y_JH(3,1), Z_JH(3,1), Y_JH0(3,1), Z_JH0(3,1)

REAL(KIND = dp_ud)::  M1(3,3), M2(3,3), M3(3,3)

! x_axis of JH is along with Jupiter-Sun vector

    CALL JUPITERPV(tdb,PV)

    R = RESHAPE((/PV(1,1),PV(2,1),PV(3,1)/),(/3,1/))

    CALL NORM(R,RR)

    X_JH=-R/RR

!z_axis of JH is related to the Sun spin aixs

    RA = 286.13D0

    DE = 63.87D0

    W = 84.176D0+14.1844D0*(tdb-2451545.D0)

!transform the Sun spin axis to JCRS

    CALL ROTATION_Z(-(90.D0+RA),M1)
    CALL ROTATION_X(-(90.D0-DE),M2)
    CALL ROTATION_Z(-W,M3)

    Z_SUN = MATMUL(MATMUL(M1,MATMUL(M2,M3)),RESHAPE((/0,0,1/),(/3,1/)))

!y_axis of JH is determind by its x-axis and the Sun spin axis

    CALL CROSS_PRODUCT(Z_SUN,X_JH,Y_JH0)

    CALL NORM(Y_JH0,N1)

    Y_JH = Y_JH0/N1

    CALL CROSS_PRODUCT(X_JH,Y_JH,Z_JH0)

    CALL NORM(Z_JH0,N2)

    Z_JH = Z_JH0/N2

    RETURN
END SUBROUTINE JH_JCRS

!==================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
!==================================================================================
SUBROUTINE  JSM_JCRS(tdb, model, X_JSM, Y_JSM, Z_JSM)
!------------------------------------------------------------------------------
! Calculate the expressions of JSM coordinate axes in JCRS
! Here we use JRM33 (Connerney 2021) to provide Jupiter's magnetic dipole moment

! INPUT:
!       tdb      time, julian date in TDB
! OUTPUT:
!       X_JSO    expression of JSO x-aixs in JCRS
!       Y_JSO    expression of JSO y-aixs in JCRS
!       Z_JSO    expression of JSO z-aixs in JCRS
!------------------------------------------------------------------------------

IMPLICIT NONE

REAL(KIND = dp_ud)::  PV(6,1), R(3,1), RR, RA, DE, W, tdb, N1, N2,COLAT_DIPOLE,ELONG_DIPOLE

REAL(KIND = dp_ud)::  X_JSM(3,1), Y_JSM(3,1), Z_JSM(3,1), Z_M(3,1), Y_JSM0(3,1)

REAL(KIND = dp_ud)::  M1(3,3), M2(3,3), M3(3,3), DM(3,1), Z_JSM0(3,1)

CHARACTER(LEN=*):: model

!x-axis of JSM is along with Jupiter-Sun vector

    CALL JUPITERPV(tdb,PV)

    R = RESHAPE((/PV(1,1),PV(2,1),PV(3,1)/),(/3,1/))

    CALL NORM(R,RR)

    X_JSM = -R/RR

!z-axis of JSM is related to Jupiter's magnetic dipole moment

!rectangular coordinates of dipole moment in S3RH
    IF ( model=='VIP4' ) THEN
        
        COLAT_DIPOLE=COLAT_DIPOLE_VIP4
        ELONG_DIPOLE=ELONG_DIPOLE_VIP4
        
    ELSE IF ( model=='JRM09' ) THEN
        
        COLAT_DIPOLE=COLAT_DIPOLE_JRM09
        ELONG_DIPOLE=ELONG_DIPOLE_JRM09
        
    ELSE IF ( model=='JRM33' ) THEN
        
        COLAT_DIPOLE=COLAT_DIPOLE_JRM33
        ELONG_DIPOLE=ELONG_DIPOLE_JRM33
        
    END IF

    Z_M(1,1)=DSIN(COLAT_DIPOLE*DEG2RAD)*DCOS(ELONG_DIPOLE*DEG2RAD)
    Z_M(2,1)=DSIN(COLAT_DIPOLE*DEG2RAD)*DSIN(ELONG_DIPOLE*DEG2RAD)
    Z_M(3,1)=DCOS(COLAT_DIPOLE*DEG2RAD)

!transform dipole moment from S3RH to JCRS

    CALL ROTATIONAL_ELEMENTS(tdb,RA,DE,W)

    CALL ROTATION_Z(-(90+RA),M1)
    CALL ROTATION_X(-(90-DE),M2)
    CALL ROTATION_Z(-W,M3)

    DM = MATMUL(MATMUL(M1,MATMUL(M2,M3)),Z_M)

    CALL CROSS_PRODUCT(DM,X_JSM,Y_JSM0)

    CALL NORM(Y_JSM0,N1)

    Y_JSM = Y_JSM0/N1

    CALL CROSS_PRODUCT(X_JSM,Y_JSM,Z_JSM0)

    CALL NORM(Z_JSM0,N2)

    Z_JSM=Z_JSM0/N2

    RETURN
END SUBROUTINE  JSM_JCRS

!==================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
!==================================================================================
SUBROUTINE  JM_JCRS(tdb, model, X_JM, Y_JM, Z_JM)
!------------------------------------------------------------------------------
! Calculate the expressions of JM coordinate axes in JCRS
!
! INPUT:
!       tdb      time, julian date in TDB
! OUTPUT:
!       X_JSO    expression of JSO x-aixs in JCRS
!       Y_JSO    expression of JSO y-aixs in JCRS
!       Z_JSO    expression of JSO z-aixs in JCRS
!------------------------------------------------------------------------------
IMPLICIT NONE

REAL(KIND = dp_ud)::  RA, DE, W, tdb, N1, N2,COLAT_DIPOLE,ELONG_DIPOLE

REAL(KIND = dp_ud)::  X_JM(3,1), Y_JM(3,1), Z_JM(3,1), X_S3RH(3,1), Y_S3RH(3,1), Z_S3RH(3,1)

REAL(KIND = dp_ud)::  M1(3,3), M2(3,3), M3(3,3), Y_M(3,1), Y_JM0(3,1), X_JM0(3,1), Z_M(3,1)

CHARACTER(LEN=*) :: model

!rectangular coordinates of dipole moment in S3RH

    IF ( model=='VIP4' ) THEN
        
        COLAT_DIPOLE=COLAT_DIPOLE_VIP4
        ELONG_DIPOLE=ELONG_DIPOLE_VIP4
        
    ELSE IF ( model=='JRM09' ) THEN
        
        COLAT_DIPOLE=COLAT_DIPOLE_JRM09
        ELONG_DIPOLE=ELONG_DIPOLE_JRM09
        
    ELSE IF ( model=='JRM33' ) THEN
        
        COLAT_DIPOLE=COLAT_DIPOLE_JRM33
        ELONG_DIPOLE=ELONG_DIPOLE_JRM33
        
    END IF

    Z_M(1,1) = DSIN(COLAT_DIPOLE*DEG2RAD)*DCOS(ELONG_DIPOLE*DEG2RAD)
    Z_M(2,1) = DSIN(COLAT_DIPOLE*DEG2RAD)*DSIN(ELONG_DIPOLE*DEG2RAD)
    Z_M(3,1) = DCOS(COLAT_DIPOLE*DEG2RAD)

!z-axis of JM

    CALL ROTATIONAL_ELEMENTS(tdb,RA,DE,W)

    CALL ROTATION_Z(-(90+RA),M1)
    CALL ROTATION_X(-(90-DE),M2)
    CALL ROTATION_Z(-W,M3)

    Z_JM=MATMUL(MATMUL(M1,MATMUL(M2,M3)),Z_M)
   
!y-axis of JM is determined by z axes of S3RH and dipole moment.  Rectangular coordinate in JCRS

    CALL S3RH_JCRS(tdb,X_S3RH,Y_S3RH,Z_S3RH)

    CALL CROSS_PRODUCT(Z_S3RH,Z_JM,Y_JM0)

    CALL NORM(Y_JM0,N1)

    Y_JM=Y_JM0/N1

!x_axis of JM

    CALL CROSS_PRODUCT(Y_JM,Z_JM,X_JM0)

    CALL NORM(X_JM0,N2)

    X_JM=X_JM0/N2

    RETURN
END SUBROUTINE  JM_JCRS

!==================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
!==================================================================================
SUBROUTINE  JSMAG_JCRS(tdb,model,X_JSMAG,Y_JSMAG,Z_JSMAG)

!------------------------------------------------------------------------------
! Calculate the expressions of JSMAG coordinate axes in JCRS

! OUTPUT:
!       X_JSMAG    expression of JSMAG x-aixs in JCRS
!       Y_JSMAG    expression of JSMAG y-aixs in JCRS
!       Z_JSMAG    expression of JSMAG z-aixs in JCRS
!------------------------------------------------------------------------------
IMPLICIT NONE

REAL(KIND = dp_ud)::  X_JM(3,1),Y_JM(3,1),Z_JM(3,1),tdb,N1,N2,PV(6,1),JS(3,1)

REAL(KIND = dp_ud)::  RR,R(3,1),X_JSMAG(3,1),Y_JSMAG(3,1),Z_JSMAG(3,1),Y_JSMAG0(3,1),X_JSMAG0(3,1)

CHARACTER(LEN=*):: model 

    CALL JUPITERPV(tdb,PV)

    R = RESHAPE((/PV(1,1),PV(2,1),PV(3,1)/),(/3,1/))

    CALL NORM(R,RR)

    JS = -R/RR

!z-axis of JSMAG is same with z-axis of JM

    CALL JM_JCRS(tdb,model,X_JM,Y_JM,Z_JM)

    Z_JSMAG = Z_JM

    CALL CROSS_PRODUCT(Z_JSMAG,JS,Y_JSMAG0)

    CALL NORM(Y_JSMAG0,N1)

    Y_JSMAG = Y_JSMAG0/N1

    CALL CROSS_PRODUCT(Y_JSMAG,Z_JSMAG,X_JSMAG0)

    CALL NORM(X_JSMAG0,N2)

    X_JSMAG=X_JSMAG0/N2

    RETURN
END SUBROUTINE  JSMAG_JCRS

!==================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
!==================================================================================

END MODULE MOD_CS_TRANSFORMATION

