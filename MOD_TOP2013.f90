! MOD_TOP2013.f90

! THIS MODULE IS REVISED BY YUJIE WANG ON 23/4/2023 TO 
! CALCULATE THE POSITION AND VELOCITY OF JOVAIN BARYCENTER.
! ORIGINAL F77 CODE IS CONVERTED TO F90
! USERS ARE REFERRED TO SEE ORIGINAL F77 FOR DETAILED INSTRUCTIONS.
! 
! NOTES TO USER:
! THE MODULE MOD_TOP2013 IS TAKEN FROM TOP2013 SOLUTIONS 
!       (https://ftp.imcce.fr/pub/ephem/planets/top2013/), and refers to Simon et al. (2013) 
!       'New analytical planetary theories VSOP2013 and TOP2013'.  
! -------------------------------------------------------------------------
MODULE MOD_TOP2013
  
USE MOD_PRECISION

IMPLICIT NONE

PRIVATE
!********************************************************************************

PUBLIC JUPITERPV

PUBLIC TOP2013

PUBLIC ELLXYZ

!********************************************************************************
CONTAINS
!==================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
!==================================================================================
SUBROUTINE JUPITERPV(t,r2)
! -----------------------------------------------------------------------------
! RETURN THE POSITION AND VELOCITY OF JOVIAN BARYCENTER
! t here is tdb
IMPLICIT NONE

INTEGER(KIND = i4_ud):: i,j,ip,nul,n,ierr
REAL(KIND = dp_ud):: t0,t,tj
REAL(KIND = dp_ud):: eps,phi,ceps,seps,cphi,sphi
REAL(KIND = dp_ud):: pi,dpi,dgrad,sdrad
CHARACTER(LEN = 40):: fich

CHARACTER(LEN = 8),DIMENSION(5):: body
REAL(KIND = dp_ud),DIMENSION(6):: el,r1,r2
REAL(KIND = dp_ud),DIMENSION(3,3):: rot(3,3)
REAL(KIND = dp_ud),PARAMETER:: t2000=2451545.d0

! --- TOP2013 FILE -----------------------------------------------------
    nul=10
    fich='TOP2013.dat'

! --- Rotation Matrix : Ecliptic -> Equator ----------------------------

    dpi   = 2.D0*pi
    sdrad = DEG2RAD/3600.D0

    eps  = (23.D0+26.D0/60.D0+21.41136D0/3600.d0)*DEG2RAD
    phi  = -0.05188D0*sdrad
    ceps = DCOS(eps)
    seps = DSIN(eps)
    cphi = DCOS(phi)
    sphi = DSIN(phi)

    rot(1,1) =  cphi
    rot(1,2) = -sphi*ceps
    rot(1,3) =  sphi*seps
    rot(2,1) =  sphi
    rot(2,2) =  cphi*ceps
    rot(2,3) = -cphi*seps
    rot(3,1) =  0.D0
    rot(3,2) =  seps
    rot(3,3) =  ceps
!
! --- Computation ------------------------------------------------------
!
    ip = 5 ! JUPITER
    tj = t-t2000
    CALL TOP2013 (tj,ip,fich,nul,el,ierr)
        IF (ierr /= 0) STOP 1
    CALL ELLXYZ (ip,el,r1,ierr)  !r1 is elliptic position
        IF (ierr /= 0) STOP 2
    DO i = 1,3
        r2(i) = 0.D0
        r2(i+3) = 0.D0
        DO j=1,3
            r2(i) = r2(i)+rot(i,j)*r1(j)
            r2(i+3) = r2(i+3)+rot(i,j)*r1(j+3)
        END DO
    END DO
    RETURN
END SUBROUTINE JUPITERPV
!==================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
!==================================================================================
SUBROUTINE TOP2013 (tdj,ipla,nFILE,IFILE,el,ierr)
! ----------------------------------------------------------------------
!
!     Reference : GF-JLS-2013
!
! --- Object -----------------------------------------------------------
!
!     Substitution of time in TOP2013 planetary solutions.
!     Variables :  Elliptic coordinates a, l, k, h, q, p.
!     Frame :      Dynamical equinox and ecliptic J2000.
!     Time scale : TDB, Temps Dynamique Barycentrique.
!     
! --- Input ------------------------------------------------------------
!
!     tdj      Julian DATE in dynamical time TDB from J2000 (REAL*8).
!
!     ipla     Planet index (INTEGER).
!              5: Jupiter
!              6: Saturn
!              7: Uranus
!              8: Neptune
!              9: Pluto
!
!     nFILE    Name of the FILE corresponding to the planet (CHARACTER*40)
!
!     IFILE    Logical unit index of the FILE (INTEGER).
!
! --- Output ----------------------------------------------------------- 
!
!     el(6)    Table of elliptic coordinates (REAL*8).
!              el(1): semi-major axis (au)
!              el(2): mean longitude (rd)
!              el(3): k = e*COS(pi) (rd)
!              el(4): h = e*SIN(pi) (rd)
!              el(5): q = SIN(i/2)*COS(omega) (rd)
!              el(6): p = SIN(i/2)*SIN(omega) (rd)
!                    e:     eccentricity
!                    pi:    perihelion longitude
!                    i:     inclination
!                    omega: ascENDing node longitude
!
!     ierr :   Error index (INTEGER).
!              ierr=0: no error.
!              ierr=1: planet index error.
!              ierr=2: FILE error (OPEN).
!              ierr=3: FILE error (READ).
!
! --- Declarations -----------------------------------------------------
!
IMPLICIT NONE

REAL(KIND = dp_ud):: tdj
INTEGER(KIND = i4_ud):: ipla
CHARACTER(LEN = 40):: nFILE
INTEGER(KIND = i4_ud):: IFILE
REAL(KIND = dp_ud),DIMENSION(6):: el
INTEGER(KIND = i4_ud):: ierr

INTEGER(KIND = i4_ud):: ip,it,iv,nt
INTEGER(KIND = i4_ud):: i,j,k,n,io,IFILE0
REAL(KIND = dp_ud):: dmu,DATE,arg,xl
CHARACTER(LEN = 26):: ca,sa
CHARACTER(LEN = 40):: nFILE0
!
! --- FILE PARAMETERs  -------------------------------------------------
!
LOGICAL:: fEXIST
INTEGER(KIND = i4_ud),DIMENSION(5,6,0:12):: limit
INTEGER(KIND = i4_ud),DIMENSION(5):: nterm
REAL(KIND = dp_ud),DIMENSION(336428):: c,s
INTEGER(KIND = i4_ud),DIMENSION(336428):: m
REAL(KIND = dp_ud),DIMENSION(0:12):: time
!
! --- Initialization ---------------------------------------------------
!
REAL(KIND = dp_ud),PARAMETER:: dpi = 2.D0*pi

REAL(KIND = dp_ud),DIMENSION(5)::  freq

    freq = (/0.5296909622785881d+03,&
             0.2132990811942489d+03,&
             0.7478166163181234d+02,&
             0.3813297236217556d+02,&
             0.2533566020437000d+02/)

    IFILE0 = 0
    nFILE0 = '*'

! --- Check planet index ----------------------------------------------- 

    ierr = 1 ! A WRONG PLANET INDEX
    IF (ipla < 5 .OR. ipla > 9) RETURN
!
! --- Time ------------------------------------------------------------- 
!
    DATE = tdj/365250.D0
    time(0) = 1.D0
    DO i = 1,12
        time(i) = time(i-1)*DATE
    END DO
!
! --- FILE READing ----------------------------------------------------- 
!
    IF (IFILE /= IFILE0 .OR. nFILE /= nFILE0) THEN
    !
        ierr = 2
        INQUIRE (FILE = nFILE, EXIST = fEXIST)
        IF ( .NOT. fEXIST) RETURN
        OPEN (IFILE, FILE=nFILE, STATUS='old', IOSTAT=io)
        IF (io /= 0) RETURN
            ierr = 3
            k = 0
        DO ip = 1,5
            nterm(ip) = 0
            DO iv = 1,6
                DO it = 0,12
                    limit(ip,iv,it) = 0
                END DO
            END DO
        END DO
    
        DO
            READ (IFILE,1001,IOSTAT = io) ip,iv,it,nt
            IF (io == -1) EXIT
            IF (io /= 0) RETURN
            IF (nt == 0) CYCLE
            limit(ip-4,iv,it) = nt
            DO n = 1,nt
                k = k+1
                nterm(ip-4) = nterm(ip-4)+1
                READ (IFILE,1002,IOSTAT=io) m(k),ca,sa
                IF (io /= 0) RETURN
                ca(23:23)='D'
                sa(23:23)='D'
                READ (ca,'(d26.16)') c(k)
                READ (sa,'(d26.16)') s(k)
            ENDDO
        ENDDO
    
        DO i = 5,2,-1
            nterm(i) = 0
            DO j = i-1,1,-1
                nterm(i) = nterm(i)+nterm(j)
            ENDDO
        ENDDO
        nterm(1) = 0
        CLOSE (IFILE)
        dmu = (freq(1)-freq(2))/880.D0
        IFILE0 = IFILE
        nFILE0 = nFILE
    ENDIF
!
! --- Substitution of time ---------------------------------------------                    
!
    k = nterm(ipla-4)
    DO iv = 1,6
        el(iv) = 0.D0
        DO it = 0,12
            IF (limit(ipla-4,iv,it) == 0) CYCLE
            DO nt = 1,limit(ipla-4,iv,it)
                k = k+1
                IF (iv == 2 .AND. it == 1 .AND. m(k) == 0) CYCLE
                arg = m(k)*dmu*time(1)
                el(iv) = el(iv)+time(it)*(c(k)*DCOS(arg)+s(k)*DSIN(arg))
            ENDDO
        ENDDO
    ENDDO
!
    xl = el(2)+freq(ipla-4)*time(1)
    xl = MOD(xl,dpi)
    IF (xl < 0.D0) xl = xl+dpi
    el(2) = xl
    ierr = 0
!
! --- FORMATs ----------------------------------------------------------
!
1001  FORMAT (21x,i2,12x,i2,7x,i2,2x,i6)
1002  FORMAT (1x,i8,2a26)
    
    RETURN
END SUBROUTINE TOP2013
!==================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
!==================================================================================
SUBROUTINE ELLXYZ (ibody,v,w,ierr)
!-----------------------------------------------------------------------
!
!     Reference : GF-JLS-1212
!
! --- Object -----------------------------------------------------------
!
!     Computation of planetary rectangular coordinates from elliptic 
!     variables.
!     
! --- Input ------------------------------------------------------------
!
!     ibody    Body index (INTEGER).
!              1: Mercury
!              2: Venus
!              3: Earth-Moon barycenter
!              4: Mars
!              5: Jupiter
!              6: Saturn
!              7: Uranus
!              8: Neptune
!              9: Pluto
!
!     v(6)     Elliptic variables reel (REAL*8).
!              v(1): semi-major axis (au)
!              v(2): mean longitude (rd)
!              v(3): k = e*COS(pi) (rd)
!              v(4): h = e*SIN(pi) (rd)
!              v(5): q = SIN(i/2)*COS(omega) (rd)
!              v(6): p = SIN(i/2)*SIN(omega) (rd)
!                    e:     eccentricity
!                    pi:    perihelion longitude
!                    i:     inclination
!                    omega: ascENDing node longitude
!
! --- Output ----------------------------------------------------------- 
!
!     w(6)     Rectangular coordinates (REAL*8)).
!              w(i),i=1,3 Positions  X, Y, Z (au)
!              w(i),i=4,5 Velocities X',Y',Z'(au/d)
!
!     ierr :   Error index (INTEGER).
!              ierr=0: no error.
!              ierr=1: ibody error.
!
! --- Declarations -----------------------------------------------------

IMPLICIT NONE

INTEGER(KIND = i4_ud):: ibody,ierr
REAL(KIND = dp_ud),DIMENSION(6):: v,w

COMPLEX(KIND = dp_ud):: DCMPLX
COMPLEX(KIND = dp_ud):: z,z1,z2,z3,zi,zeta,zto,zteta

REAL(KIND = dp_ud):: rgm,xa,xl,xk,xh,xq,xp,xfi,xki
REAL(KIND = dp_ud):: u,ex,ex2,ex3,gl,gm,e,dl,rsa
REAL(KIND = dp_ud):: xcw,xsw,xm,xr,xms,xmc,xn

REAL(KIND = dp_ud),PARAMETER:: dpi = 2.D0*pi
!
! --- Masses system (INPOP10A) -----------------------------------------
!
REAL(KIND = dp_ud),DIMENSION(9):: gmp

REAL(KIND = dp_ud),PARAMETER:: gmsol = 2.9591220836841438269D-04

    gmp = (/4.9125474514508118699D-11,&
            7.2434524861627027000D-10,&
            8.9970116036316091182D-10,&
            9.5495351057792580598D-11,&
            2.8253458420837780000D-07,&
            8.4597151856806587398D-08,&
            1.2920249167819693900D-08,&
            1.5243589007842762800D-08,&
            2.1886997654259696800D-12/)
!
! --- Initialization ---------------------------------------------------
!
    ierr = 11
    IF (ibody < 1 .OR. ibody > 9) RETURN
    ierr = 0
    !
    rgm = DSQRT(gmp(ibody)+gmsol)
    xa = v(1)
    xl = v(2)
    xk = v(3)
    xh = v(4)
    xq = v(5)
    xp = v(6)
    !
    ! --- Computation ------------------------------------------------------
    !
    xfi = DSQRT(1.D0-xk*xk-xh*xh)
    xki = DSQRT(1.D0-xq*xq-xp*xp)
    u = 1.D0/(1.D0+xfi)
    z = DCMPLX(xk,xh)
    ex = CDABS(z)
    ex2 = ex*ex
    ex3 = ex2*ex
    z1 = DCONJG(z)
!
    gl = DMOD(xl,dpi)
    gm = gl-DATAN2(xh,xk)
    e = gl+(ex-0.125D0*ex3)*DSIN(gm)&
        +0.5D0*ex2*DSIN(2.D0*gm)&
        +0.375D0*ex3*DSIN(3.D0*gm)
!
    DO
        z2 = DCMPLX(0.D0,e)
        zteta = CDEXP(z2)
        z3 = z1*zteta
        dl = gl-e+DIMAG(z3)
        rsa = 1.D0-DREAL(z3)
        e = e+dl/rsa
        IF (DABS(dl) < 1.d-15) EXIT
    END DO
!
    z1 = u*z*DIMAG(z3)
    z2 = DCMPLX(DIMAG(z1),-DREAL(z1))
    zto = (-z+zteta+z2)/rsa
    xcw = DREAL(zto)
    xsw = DIMAG(zto)
    xm = xp*xcw-xq*xsw
    xr = xa*rsa
    !
    w(1) = xr*(xcw-2.D0*xp*xm)
    w(2) = xr*(xsw+2.D0*xq*xm)
    w(3) = -2.D0*xr*xki*xm
    !
    xms = xa*(xh+xsw)/xfi
    xmc = xa*(xk+xcw)/xfi
    xn = rgm/xa**(1.5D0)
    !
    w(4) = xn*((2.D0*xp*xp-1.D0)*xms+2.D0*xp*xq*xmc)
    w(5) = xn*((1.D0-2.d0*xq*xq)*xmc-2.D0*xp*xq*xms)
    w(6) = 2.D0*xn*xki*(xp*xms+xq*xmc)
!
    RETURN
END SUBROUTINE ELLXYZ
!==================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
!==================================================================================
END MODULE MOD_TOP2013
