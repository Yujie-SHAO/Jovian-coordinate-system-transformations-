! MOD_PRECISION.f90
! THIS MODULE DEFINE THE PRECISION OF PARAMETERS AND SOME COMMON CONSTANTS. 
!
! CREATED BY YABING WANG ON 2022/06/07.
! COPYRIGHT 2022. ALL RIGHTS RESERVED.
! REVISED BY YABING WANG ON 2023/04/14 AND YUJIE WANG ON 06/25/2023.

MODULE MOD_PRECISION

	IMPLICIT NONE
                                                               ! user defined
        INTEGER, PARAMETER :: sp_ud = selected_real_kind(6,37) ! single precision real numbers
    	INTEGER, PARAMETER :: dp_ud = selected_real_kind(15,307) ! double precision
    	INTEGER, PARAMETER :: qp_ud = selected_real_kind(33,4931) ! quadruple precision

        INTEGER, PARAMETER :: i1_ud = selected_int_kind(2) ! kind = 1
        INTEGER, PARAMETER :: i2_ud = selected_int_kind(4) ! kind = 2
        INTEGER, PARAMETER :: i4_ud = selected_int_kind(9) ! kind = 4
        
        REAL(KIND=dp_ud),PARAMETER :: pi = DACOS(-1.D0), &
                                      COLAT_DIPOLE_JRM33 = 10.250332D0, & ! The colatitude and longitude of the dipole axis
                                      ELONG_DIPOLE_JRM33 = 163.620757D0, &      ! in S3RH coordinate system from JRM33 magnetic field model.
                                      COLAT_DIPOLE_JRM09 = 10.307870D0, & ! The colatitude and longitude of the dipole axis
                                      ELONG_DIPOLE_JRM09 = 163.388276D0, &      ! in S3RH coordinate system from JRM09 magnetic field model.
                                      COLAT_DIPOLE_VIP4 = 9.515259D0, & ! The colatitude and longitude of the dipole axis
                                      ELONG_DIPOLE_VIP4 = 159.225124D0, &      ! in S3RH coordinate system from VIP4 magnetic field model.
                                      DEG2RAD = pi/180.D0, RAD2DEG = 180.D0/pi  ! DEGREE TO RADIAN, RADIAN TO DEGREE


END MODULE MOD_PRECISION
