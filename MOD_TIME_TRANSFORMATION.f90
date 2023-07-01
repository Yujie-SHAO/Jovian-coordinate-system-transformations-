! MOD_TIME_TRANSFORMATION.f90

! THIS MODULE IS CREATED BY YABING WANG ON 07/07/2022 TO
! INCLUDE ALL RELEVANT SUBROUTINES IN SOFA TO TRANSFORM
! AMONG DIFFERENT TIME SYSTEMS.
! ORIGINAL F77 CODE ARE CONVERTED TO F90.
! USERS ARE REFERRED TO SEE ORIGINAL F77 FOR DETAILED INSTRUCTIONS.

! NOTES TO USER:
! USES ROUTINES AND COMPUTATIONS DREIVED BY US FROM SOFTWARE PROVIDED BY
! SOFA UNDER LICENSE TO US, AND DOSE NOT ITSELF CONSTITUTE SOFTWARE PROVIDED
! BY AND/OR ENDORSED BY SOFA.

! THE DATA PRECISION IS REVISED BY YABING WANG on 2023/04/14 AND YUJIE WANG ON 2023/06/25.
! ===============================================================================
MODULE  MOD_TIME_TRANSFORMATION

USE MOD_PRECISION

IMPLICIT NONE

PRIVATE

!********************************************************************************
PUBLIC  WR_CAL2JD

PUBLIC  WR_JD2CAL

PUBLIC  WR_DAT

PUBLIC  WR_DTF2D

PUBLIC  WR_UTCTAI

PUBLIC  WR_TAITT

PUBLIC  WR_TTTDB

PUBLIC  WR_UTC2TDB

PUBLIC  WR_UTC2TDB2
!********************************************************************************
!INTERFACE

CONTAINS

!==================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
!==================================================================================
SUBROUTINE WR_CAL2JD ( IY, IM, ID, DJM0, DJM, J )

! GREGORIAN CALENDAR TO JULIAN DATE.
!++++++++++++++++++++++INPUT++++++++++++++++++++++++++++
! IY, IM, AND ID: YEAR, MONTH, DAY IN GREGORIAN CALENDAR
!++++++++++++++++++++++OUTPUT++++++++++++++++++++++++++++
! DJM0:  MJD zero-point: always 2400000.5D0
! DJM:   Modified Julian Date for 0 hrs
! J: 0 = OK, -1 = bad year (JD not computed),
! -2 = bad month(JD not computed), -3 = bad day(JD computed)

IMPLICIT NONE

INTEGER(KIND=i4_ud):: IY, IM, ID
REAL(KIND=dp_ud)::  DJM0, DJM

INTEGER(KIND=i4_ud):: J, NDAYS, MY, IYPMY

!  Earliest year allowed (4800BC)
!  But this implementation rejects dates before -4799 January 1
INTEGER(KIND=i4_ud),PARAMETER:: IYMIN = -4799

! Month lengths in days
INTEGER(KIND=i4_ud):: MTAB(12)=(/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  Preset status.
    J = 0

    !  Validate year.
    IF ( IY < IYMIN ) THEN
       J = -1
    ELSE IF ( IM >= 1 .AND. IM <= 12 ) THEN ! Validate month.

        ! Days in current month.
        NDAYS = MTAB(IM)

        ! Allow for leap year.
        IF ( IM == 2 ) THEN
            IF ( MOD(IY,4) == 0 ) NDAYS = 29
            IF ( MOD(IY,100) == 0 .AND. MOD(IY,400) /= 0 ) NDAYS = 28
        END IF

        ! Validate day.
        IF ( ID < 1 .OR. ID > NDAYS ) J = -3

            ! Result.
            MY = ( IM - 14 ) / 12
            IYPMY = IY + MY
            DJM0 = 2400000.5D0
            DJM = DBLE( ( 1461 * ( IYPMY + 4800 ) ) / 4 &
                    + (  367 * ( IM-2 - 12*MY ) ) / 12 &
                    - (    3 * ( ( IYPMY + 4900 ) / 100 ) ) / 4 &
                    + ID - 2432076)

    ! Bad month
    ELSE
        J = -2
    END IF

    RETURN
END SUBROUTINE WR_CAL2JD

!==================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
!==================================================================================
SUBROUTINE WR_JD2CAL ( DJ1, DJ2, IY, IM, ID, FD, J )
! JULIAN DATE TO GREGORIAN YEAR, MONTH, DAY, AND FRACTION OF A DAY
!++++++++++++++++++++++INPUT++++++++++++++++++++++++++++
! DJ1, DJ2: Julian Date
!++++++++++++++++++++++OUTPUT++++++++++++++++++++++++++++
!     IY:     year
!     IM:     month
!     ID:     day
!     FD:     fraction of day
!     J:     status:
!                   0 = OK
!                   -1 = unacceptable date
!  Notes:
!
!  1) The earliest valid date is -68569.5 (-4900 March 1).  The
!     largest value accepted is 10^9.
!
!  2) The Julian Date is apportioned in any convenient way between
!     the arguments DJ1 and DJ2.  For example, JD=2450123.7 could
!     be expressed in any of these ways, among others:
!
!             DJ1            DJ2
!
!         2450123.7D0        0D0        (JD method)
!          2451545D0      -1421.3D0     (J2000 method)
!         2400000.5D0     50123.2D0     (MJD method)
!         2450123.5D0       0.2D0       (date & time method)
!
!     Separating integer and fraction uses the "compensated summation"
!     algorithm of Kahan-Neumaier to preserve as much precision as
!     possible irrespective of the JD1+JD2 apportionment.
IMPLICIT NONE

REAL(KIND=dp_ud):: DJ1, DJ2
INTEGER(KIND=i4_ud):: IY, IM, ID
REAL(KIND=dp_ud):: FD
INTEGER(KIND=i4_ud):: J

!  Smallest number such that 1D0+EPS /= 1D0
REAL(KIND=dp_ud),PARAMETER:: EPS = 2.2204460492503131D-16

!  Minimum and maximum allowed JD
REAL(KIND=dp_ud),PARAMETER:: DJMIN = -68569.5D0, DJMAX = 1.D9

INTEGER(KIND=i4_ud):: JD, I, L, N, K
REAL(KIND=dp_ud):: DJ, F1, F2, D, S, CS, V(2), X, T, C, F

!  Verify date is acceptable.
    DJ = DJ1 + DJ2

    IF ( DJ < DJMIN .OR. DJ > DJMAX ) THEN
       J = -1
    ELSE
       J = 0

    !     Separate day and fraction (where -0.5 <= fraction < 0.5).
       D = ANINT(DJ1)
       F1 = DJ1 - D
       JD = INT(D)
       D = ANINT(DJ2)
       F2 = DJ2 - D
       JD = JD + INT(D)

    !     Compute F1+F2+0.5 using compensated summation (Klein 2006).
       S = 0.5D0
       CS = 0.D0
       V(1) = F1
       V(2) = F2
        DO I=1,2
            X = V(I)
            T = S + X
            IF ( ABS(S) >= ABS(X) ) THEN
                 C = (S-T) + X
            ELSE
                 C = (X-T) + S
            END IF
                CS = CS + C
                S = T
            IF ( S >= 1.D0 ) THEN
                 JD = JD + 1
                 S = S - 1.D0
            END IF
        END DO
       F = S + CS
       CS = F - S

        !     Deal with negative F.
        IF ( F < 0.D0 ) THEN
            !        Compensated summation: assume that |S| <= 1.
            F = S + 1.D0
            CS = CS + ((1.D0-F) + S)
            S = F
            F = S + CS
            CS = F - S
            JD = JD - 1
        END IF

    !     Deal with F that is 1D0 or more (when rounded to double).
       IF ( (F-1.D0) >= -EPS/4.D0 ) THEN

    !        Compensated summation: assume that |S| <= 1. */
          T = S - 1.D0
          CS = CS + ((S-T) - 1.D0)
          S = T
          F = S + CS
          IF ( -EPS/2.D0 < F ) THEN
             JD = JD + 1
             F = MAX ( F, 0.D0 )
          END IF
       END IF

        !     Express day in Gregorian calendar.
        L = JD + 68569
        N = ( 4*L ) / 146097
        L = L - ( 146097*N + 3 ) / 4
        I = ( 4000 * (L+1) ) / 1461001
        L = L - ( 1461*I ) / 4 + 31
        K = ( 80*L ) / 2447
        ID = L - ( 2447*K ) / 80
        L = K / 11
        IM = K + 2 - 12*L
        IY = 100 * ( N-49 ) + I + L
        FD = F
    END IF

    RETURN
END SUBROUTINE WR_JD2CAL

!==================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
!==================================================================================
SUBROUTINE WR_DAT ( IY, IM, ID, FD, DELTAT, J )

! For a given UTC date, calculate Delta(AT) = TAI-UTC
!     :------------------------------------------:
!     :                                          :
!     :                 IMPORTANT                :
!     :                                          :
!     :  A new version of this routine must be   :
!     :  produced whenever a new leap second is  :
!     :  announced.  There are five items to     :
!     :  change on each such occasion:           :
!     :                                          :
!     :  1) The parameter NDAT must be           :
!     :     increased by 1.                      :
!     :                                          :
!     :  2) The set of DATA statements that      :
!     :     initialize the arrays IDAT and       :
!     :     DATS must be extended by one line.   :
!     :                                          :
!     :  3) The parameter IYV must be set to     :
!     :     the current year.                    :
!     :                                          :
!     :  4) The "Latest leap second" comment     :
!     :     below must be set to the new leap    :
!     :     second date.                         :
!     :                                          :
!     :  5) The "This revision" comment, later,  :
!     :     must be set to the current date.     :
!     :                                          :
!     :  Change (3) must also be carried out     :
!     :  whenever the routine is re-issued,      :
!     :  even if no leap seconds have been       :
!     :  added.                                  :
!     :                                          :
!     :  Latest leap second:  2016 December 31   :
!     :                                          :
!     :__________________________________________:
!++++++++++++++++++++++INPUT++++++++++++++++++++++++++++
! IY:    UTC:  year
! IM:          month
! ID:          day
! FD:          fraction of day
!++++++++++++++++++++++OUTPUT++++++++++++++++++++++++++++
! DELTAT:   TAI minus UTC, seconds
!     J :    status
!            1 = dubious year
!            0 = OK
!            -1 = bad year
!            -2 = bad month
!            -3 = bad day
!            -4 = bad fraction
!            -5 = internal error

!  Notes:
!
!  1) UTC began at 1960 January 1.0 (JD 2436934.5) and it is improper
!     to call the routine with an earlier date.  If this is attempted,
!     zero is returned together with a warning status.
!
!     Because leap seconds cannot, in principle, be predicted in
!     advance, a reliable check for dates beyond the valid range is
!     impossible.  To guard against gross errors, a year five or more
!     after the release year of the present routine (see parameter IYV)
!     is considered dubious.  In this case a warning status is returned
!     but the result is computed in the normal way.
!
!     For both too-early and too-late years, the warning status is J=+1.
!     This is distinct from the error status J=-1, which signifies a
!     year so early that JD could not be computed.
!
!  2) If the specified date is for a day which ends with a leap second,
!     the TAI-UTC value returned is for the period leading up to the
!     leap second.  If the date is for a day which begins as a leap
!     second ends, the TAI-UTC returned is for the period following the
!     leap second.
!
!  3) The day number must be in the normal calendar range, for example
!     1 through 30 for April.  The "almanac" convention of allowing
!     such dates as January 0 and December 32 is not supported in this
!     routine, in order to avoid confusion near leap seconds.
!
!  4) The fraction of day is used only for dates before the introduction
!     of leap seconds, the first of which occurred at the end of 1971.
!     It is tested for validity (0 to 1 is the valid range) even if not
!     used;  if invalid, zero is used and status J=-4 is returned.  For
!     many applications, setting FD to zero is acceptable;  the
!     resulting error is always less than 3 ms (and occurs only
!     pre-1972).
!
!  5) The status value returned in the case where there are multiple
!     errors refers to the first error detected.  For example, if the
!     month and day are 13 and 32 respectively, J=-2 (bad month) will be
!     returned.  The "internal error" status refers to a case that is
!     impossible but causes some compilers to issue a warning.
!
!  6) In cases where a valid result is not available, zero is returned.
!  References:
!
!  1) For dates from 1961 January 1 onwards, the expressions from the
!     file ftp://maia.usno.navy.mil/ser7/tai-utc.dat are used.
!
!  2) The 5ms timestep at 1961 January 1 is taken from 2.58.1 (p87) of
!     the 1992 Explanatory Supplement.

IMPLICIT NONE

INTEGER(KIND=i4_ud):: IY, IM, ID
REAL(KIND=dp_ud):: FD, DELTAT
INTEGER(KIND=i4_ud):: J

!  Release year for this version of WR_DAT
INTEGER(KIND=i4_ud),PARAMETER:: IYV = 2021

!  Number of Delta(AT) changes (increase by 1 for each new leap second)
INTEGER(KIND=i4_ud),PARAMETER :: NDAT = 42

!  Number of Delta(AT) expressions before leap seconds were introduced
INTEGER(KIND=i4_ud),PARAMETER :: NERA1 = 14

!  Dates (year, month) on which new Delta(AT) came into force
INTEGER(KIND=i4_ud):: IDAT(2,NDAT)

!  New Delta(AT) which came into force on the given dates
REAL(KIND=dp_ud):: DATS(NDAT)

!  Reference dates (MJD) and drift rates (s/day), pre leap seconds
REAL(KIND=dp_ud):: DRIFT(2,NERA1)

!  Miscellaneous local variables
LOGICAL:: MORE
INTEGER(KIND=i4_ud):: JS, M, N, IS
REAL(KIND=dp_ud):: DA, DJM0, DJM

!  Dates, Delta(AT)s, reference dates, and drift rates
DATA ((IDAT(M,N),M=1,2),DATS(N),(DRIFT(M,N),M=1,2),N=1,14)&
  / 1960,  1,  1.4178180D0, 37300D0, 0.001296D0, &
    1961,  1,  1.4228180D0, 37300D0, 0.001296D0, &
    1961,  8,  1.3728180D0, 37300D0, 0.001296D0, &
    1962,  1,  1.8458580D0, 37665D0, 0.0011232D0, &
    1963, 11,  1.9458580D0, 37665D0, 0.0011232D0, &
    1964,  1,  3.2401300D0, 38761D0, 0.001296D0, &
    1964,  4,  3.3401300D0, 38761D0, 0.001296D0, &
    1964,  9,  3.4401300D0, 38761D0, 0.001296D0, &
    1965,  1,  3.5401300D0, 38761D0, 0.001296D0, &
    1965,  3,  3.6401300D0, 38761D0, 0.001296D0, &
    1965,  7,  3.7401300D0, 38761D0, 0.001296D0, &
    1965,  9,  3.8401300D0, 38761D0, 0.001296D0, &
    1966,  1,  4.3131700D0, 39126D0, 0.002592D0, &
    1968,  2,  4.2131700D0, 39126D0, 0.002592D0 /

!  Dates and Delta(AT)s
DATA ((IDAT(M,N),M=1,2),DATS(N),N=15,30) &
  / 1972,  1, 10D0, &
    1972,  7, 11D0, &
    1973,  1, 12D0, &
    1974,  1, 13D0, &
    1975,  1, 14D0, &
    1976,  1, 15D0, &
    1977,  1, 16D0, &
    1978,  1, 17D0, &
    1979,  1, 18D0, &
    1980,  1, 19D0, &
    1981,  7, 20D0, &
    1982,  7, 21D0, &
    1983,  7, 22D0, &
    1985,  7, 23D0, &
    1988,  1, 24D0, &
    1990,  1, 25D0 /

DATA ((IDAT(M,N),M=1,2),DATS(N),N=31,NDAT) &
  / 1991,  1, 26D0, &
    1992,  7, 27D0, &
    1993,  7, 28D0, &
    1994,  7, 29D0, &
    1996,  1, 30D0, &
    1997,  7, 31D0, &
    1999,  1, 32D0, &
    2006,  1, 33D0, &
    2009,  1, 34D0, &
    2012,  7, 35D0, &
    2015,  7, 36D0, &
    2017,  1, 37D0 /

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  Initialize the result to zero and the status to OK.
    DA = 0.D0
    JS = 0

!  If invalid fraction of a day, set error status and give up.
    IF ( FD < 0.D0 .OR. FD > 1.D0 ) THEN
       JS = -4
       GO TO 9000
    END IF

!  Convert the date into an MJD.
    CALL WR_CAL2JD ( IY, IM, ID, DJM0, DJM, JS )

!  If invalid year, month, or day, give up.
    IF ( JS < 0 ) GO TO 9000

!  If pre-UTC year, set warning status and give up.
    IF ( IY < IDAT(1,1) ) THEN
       JS = 1
       GO TO 9000
    END IF

!  If suspiciously late year, set warning status but proceed.
    IF ( IY > IYV+5 ) JS = 1

!  Combine year and month.
    M = 12*IY+IM

!  Find the most recent table entry.
    IS = 0
    MORE = .TRUE.
    DO N=NDAT,1,-1
       IF ( MORE ) THEN
          IS = N
          MORE = M < ( 12*IDAT(1,N) + IDAT(2,N) )
       END IF
    END DO

!  Prevent underflow warnings.
    IF ( IS < 1 ) THEN
       JS = -5
       GO TO 9000
    END IF

!  Get the Delta(AT).
    DA = DATS(IS)

!  If pre-1972, adjust for drift.
    IF ( IS <= NERA1 ) DA = DA + &
                    ( DJM + FD - DRIFT(1,IS) ) * DRIFT(2,IS)

!  Return the Delta(AT) value and the status.
9000 CONTINUE

    DELTAT = DA

    J = JS

    RETURN
END SUBROUTINE WR_DAT

!==================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
!==================================================================================

SUBROUTINE WR_DTF2D(SCALE, IY, IM, ID, IHR, IMN, SEC, D1, D2, J)

!  Encode date and time fields into 2-part Julian Date (or in the case
!  of UTC a quasi-JD form that includes special provision for leap
!  seconds).

!++++++++++++++++++++++INPUT++++++++++++++++++++++++++++
! SCALE      time scale ID
! IY,IM,ID   year, month, day in Gregorian calendar
! IHR,IMN    hour, minute
! SEC        seconds
!++++++++++++++++++++++OUTPUT++++++++++++++++++++++++++++
! D1,D2      2-part Julian Date
! J          status: +3 = both of next two
!                    +2 = time is after end of day
!                    +1 = dubious year
!                     0 = OK
!                    -1 = bad year
!                    -2 = bad month
!                    -3 = bad day
!                    -4 = bad hour
!                    -5 = bad minute
!                    -6 = bad second (<0)

!  Notes:
!
!  1) SCALE identifies the time scale.  Only the value 'UTC' (in upper
!     case) is significant, and enables handling of leap seconds (see
!     Note 4).
!
!  2) For calendar conventions and limitations, see iau_CAL2JD.
!
!  3) The sum of the results, D1+D2, is Julian Date, where normally D1
!     is the Julian Day Number and D2 is the fraction of a day.  In the
!     case of UTC, where the use of JD is problematical, special
!     conventions apply:  see the next note.
!
!  4) JD cannot unambiguously represent UTC during a leap second unless
!     special measures are taken.  The SOFA internal convention is that
!     the quasi-JD day represents UTC days whether the length is 86399,
!     86400 or 86401 SI seconds.  In the 1960-1972 era there were
!     smaller jumps (in either direction) each time the linear UTC(TAI)
!     expression was changed, and these "mini-leaps" are also included
!     in the SOFA convention.
!
!  5) The warning status "time is after end of day" usually means that
!     the SEC argument is greater than 60D0.  However, in a day ending
!     in a leap second the limit changes to 61D0 (or 59D0 in the case of
!     a negative leap second).
!
!  6) The warning status "dubious year" flags UTCs that predate the
!     introduction of the time scale or that are too far in the future
!     to be trusted.  See iau_DAT for further details.
!
!  7) Only in the case of continuous and regular time scales (TAI, TT,
!     TCG, TCB and TDB) is the result D1+D2 a Julian Date, strictly
!     speaking.  In the other cases (UT1 and UTC) the result must be
!     used with circumspection;  in particular the difference between
!     two such results cannot be interpreted as a precise time
!     interval.
IMPLICIT NONE

CHARACTER(LEN = *):: SCALE
INTEGER(KIND=i4_ud):: IY, IM, ID, IHR, IMN, J
REAL(KIND=dp_ud):: SEC, D1, D2

!  Days to seconds
REAL(KIND=dp_ud),PARAMETER:: D2S = 86400.D0

INTEGER(KIND=i4_ud):: JS, IY2, IM2, ID2
REAL(KIND=dp_ud):: DJ, W, DAY, SECLIM, DAT0, DAT12, DAT24, DLEAP, TIME

! Today's Julian Day Number.

    CALL WR_CAL2JD ( IY, IM, ID, DJ, W, JS )
        IF ( JS /= 0 ) GO TO 9
            DJ = DJ + W

! Day length and final minute length in seconds (provisional).
            DAY = D2S
            SECLIM = 60.D0

!  Deal with the UTC leap second case.
        IF ( SCALE == 'UTC' ) THEN

!     TAI-UTC at 0h today.
         CALL WR_DAT ( IY, IM, ID, 0.D0, DAT0, JS )
         IF ( JS < 0 ) GO TO 9

!     TAI-UTC at 12h today (to detect drift).
         CALL WR_DAT ( IY, IM, ID, 0.5D0, DAT12, JS )
         IF ( JS < 0 ) GO TO 9

!    TAI-UTC at 0h tomorrow (to detect jumps).
         CALL WR_JD2CAL ( DJ, 1.5D0, IY2, IM2, ID2, W, JS )
         IF ( JS /= 0 ) GO TO 9
         CALL WR_DAT ( IY2, IM2, ID2, 0.D0, DAT24, JS )
         IF ( JS < 0 ) GO TO 9

!     Any sudden change in TAI-UTC between today and tomorrow.
         DLEAP = DAT24 - ( 2.D0 * DAT12 - DAT0 )

!     If leap second day, correct the day and final minute lengths.
         DAY = DAY + DLEAP
         IF ( IHR == 23 .AND. IMN == 59 ) SECLIM = SECLIM + DLEAP

!    End of UTC-specific actions.
        END IF

!  Validate the time.
      IF ( IHR >= 0 .AND. IHR <= 23 ) THEN
         IF ( IMN >= 0 .AND. IMN <= 59 ) THEN
            IF ( SEC >= 0D0 ) THEN
               IF ( SEC >= SECLIM ) THEN
                  JS = JS + 2
               END IF
            ELSE
               JS = -6
            END IF
         ELSE
            JS = -5
         END IF
      ELSE
         JS = -4
      END IF
      IF ( JS < 0 ) GO TO 9

!  The time in days.
      TIME = (60.D0*DBLE(60*IHR+IMN)+SEC) / DAY

!  Return the date and time.
      D1 = DJ
      D2 = TIME

! Return the status.
 9    CONTINUE
      J = JS

!  Finished.
    RETURN
END SUBROUTINE WR_DTF2D

!==================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
!==================================================================================

SUBROUTINE WR_UTCTAI ( UTC1, UTC2, TAI1, TAI2, J )
!  Time scale transformation:  Coordinated Universal Time, UTC, to
!  International Atomic Time, TAI.

!++++++++++++++++++++++INPUT++++++++++++++++++++++++++++
!  UTC1,UTC2    UTC as a 2-part quasi Julian Date
!++++++++++++++++++++++OUTPUT++++++++++++++++++++++++++++
!  TAI1,TAI2    TAI as a 2-part Julian Date (Note 5)
!  J            status: +1 = dubious year (Note 3)
!                        0 = OK
!                       -1 = unacceptable date

!  Notes:
!
!  1) UTC1+UTC2 is quasi Julian Date (see Note 2), apportioned in any
!     convenient way between the two arguments, for example where UTC1
!     is the Julian Day Number and UTC2 is the fraction of a day.
!
!  2) JD cannot unambiguously represent UTC during a leap second unless
!     special measures are taken.  The convention in the present routine
!     is that the JD day represents UTC days whether the length is
!     86399, 86400 or 86401 SI seconds.  In the 1960-1972 era there were
!     smaller jumps (in either direction) each time the linear UTC(TAI)
!     expression was changed, and these "mini-leaps" are also included
!     in the SOFA convention.
!
!  3) The warning status "dubious year" flags UTCs that predate the
!     introduction of the time scale or that are too far in the future
!     to be trusted.  See iau_DAT for further details.
!
!  4) The routine iau_DTF2D converts from calendar date and time of day
!     into 2-part Julian Date, and in the case of UTC implements the
!     leap-second-ambiguity convention described above.
!
!  5) The returned TAI1,TAI2 are such that their sum is the TAI Julian
!     Date.

IMPLICIT NONE

REAL(KIND = dp_ud) :: UTC1, UTC2, TAI1, TAI2
INTEGER(KIND = i4_ud) :: J

!  Days to seconds
REAL(KIND = dp_ud),PARAMETER :: D2S = 86400.D0

LOGICAL:: BIG1
INTEGER(KIND = i4_ud) :: IY, IM, ID, JS, IYT, IMT, IDT
REAL(KIND = dp_ud) :: U1, U2, FD, DAT0, DAT12, W, DAT24, DLOD, DLEAP, Z1, Z2, A2

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  Put the two parts of the UTC into big-first order.
    BIG1 = ABS(UTC1) >= ABS(UTC2)
    IF ( BIG1 ) THEN
       U1 = UTC1
       U2 = UTC2
    ELSE
       U1 = UTC2
       U2 = UTC1
    END IF

!  Get TAI-UTC at 0h today.
    CALL WR_JD2CAL ( U1, U2, IY, IM, ID, FD, JS )
    IF ( JS /= 0 ) GO TO 9
    CALL WR_DAT ( IY, IM, ID, 0.D0, DAT0, JS )
    IF ( JS < 0 ) GO TO 9

!  Get TAI-UTC at 12h today (to detect drift).
    CALL WR_DAT ( IY, IM, ID, 0.5D0, DAT12, JS )
    IF ( JS < 0 ) GO TO 9

!  Get TAI-UTC at 0h tomorrow (to detect jumps).
    CALL WR_JD2CAL ( U1+1.5D0, U2-FD, IYT, IMT, IDT, W, JS )
    IF ( JS /= 0 ) GO TO 9
    CALL WR_DAT ( IYT, IMT, IDT, 0.D0, DAT24, JS )
    IF ( JS < 0 ) GO TO 9

!  Separate TAI-UTC change into per-day (DLOD) and any jump (DLEAP).
    DLOD = 2.D0 * ( DAT12 - DAT0 )
    DLEAP = DAT24 - ( DAT0 + DLOD )

!  Remove any scaling applied to spread leap into preceding day.
    FD = FD * (D2S+DLEAP)/D2S

!  Scale from (pre-1972) UTC seconds to SI seconds.
    FD = FD * (D2S+DLOD)/D2S

!  Today's calendar date to 2-part JD.
    CALL WR_CAL2JD ( IY, IM, ID, Z1, Z2, JS )
    IF ( JS/=0 ) GO TO 9

!  Assemble the TAI result, preserving the UTC split and order.
    A2 = Z1 - U1
    A2 = ( A2 + Z2 ) + ( FD + DAT0/D2S )
    IF ( BIG1 ) THEN
       TAI1 = U1
       TAI2 = A2
    ELSE
       TAI1 = A2
       TAI2 = U1
    END IF

!  Status.
9    CONTINUE

    J = JS

    RETURN
END SUBROUTINE WR_UTCTAI
!==================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
!==================================================================================
SUBROUTINE WR_TAITT ( TAI1, TAI2, TT1, TT2, J )
!  Time scale transformation:  International Atomic Time, TAI, to
!  Terrestrial Time, TT.
!++++++++++++++++++++++INPUT++++++++++++++++++++++++++++
!  TAI1,TAI2    TAI as a 2-part Julian Date
!++++++++++++++++++++++OUTPUT++++++++++++++++++++++++++++
!  TT1,TT2      TT as a 2-part Julian Date
!  J            status:  0 = OK
!  Note:
!
!     TAI1+TAI2 is Julian Date, apportioned in any convenient way
!     between the two arguments, for example where TAI1 is the Julian
!     Day Number and TAI2 is the fraction of a day.  The returned
!     TT1,TT2 follow suit.

IMPLICIT NONE

REAL(KIND=dp_ud) :: TAI1, TAI2, TT1, TT2
INTEGER(KIND=i4_ud) :: J

!  TT minus TAI (days).
REAL(KIND=dp_ud),PARAMETER:: DTAT = 32.184D0/86400.D0

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  Result, safeguarding precision.
    IF ( ABS(TAI1) > ABS(TAI2) ) THEN
       TT1 = TAI1
       TT2 = TAI2 + DTAT
    ELSE
       TT1 = TAI1 + DTAT
       TT2 = TAI2
    END IF

!  Status (always OK).
    J = 0

    RETURN
END SUBROUTINE WR_TAITT

!==================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
!==================================================================================

SUBROUTINE WR_TTTDB ( TT1, TT2, TDB )
!  Time scale transformation:  Terrestrial Time, TT, to Barycentric
!  Dynamical Time, TDB.
!++++++++++++++++++++++INPUT++++++++++++++++++++++++++++
!  TT1,TT2      TT as a 2-part Julian Date
!++++++++++++++++++++++OUTPUT++++++++++++++++++++++++++++
!  TDB    TDB as a Julian Date

IMPLICIT NONE

    REAL(KIND=dp_ud) :: TT1, TT2, TDB, g

    g=6.24D0 + 0.017202D0 * ( TT1 + TT2 - 2451545.D0 )

    TDB = TT1 + TT2 + 0.001657D0 * DSIN ( g ) /24.D0/3600.D0 ! UNIT IN DAY

    RETURN
END SUBROUTINE WR_TTTDB

!==================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
!==================================================================================
SUBROUTINE WR_UTC2TDB(yr, mon, day, hr, min, ss, tdb)

!  Time scale transformation:  Universal Coordinate Time, UTC, to Barycentric
!  Dynamical Time, TDB.
!++++++++++++++++++++++INPUT++++++++++++++++++++++++++++
!  UTC     format: year, month, day, hour, minute, second 
!++++++++++++++++++++++OUTPUT++++++++++++++++++++++++++++
!  TDB     Julian date 
    
IMPLICIT NONE

REAL(KIND=dp_ud) ::  ss,tdb,UTC1,UTC2,TAI1,TAI2,TT1,TT2
INTEGER(KIND=i4_ud) :: yr,mon,day,hr,min,J    
    
CALL WR_DTF2D('UTC',yr,mon,day,hr,min,ss,UTC1,UTC2,J)

CALL WR_UTCTAI(UTC1,UTC2,TAI1,TAI2,J)

CALL WR_TAITT(TAI1,TAI2,TT1,TT2,J)

CALL WR_TTTDB(TT1,TT2,tdb)

RETURN

END SUBROUTINE WR_UTC2TDB
!==================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
!==================================================================================
SUBROUTINE WR_UTC2TDB2(yr, doy, hr, min, ss, tdb)

!  Time scale transformation:  Universal Coordinate Time, UTC, to Barycentric
!  Dynamical Time, TDB.
!++++++++++++++++++++++INPUT++++++++++++++++++++++++++++
!  UTC     format: year, day of year, hour, minute, second 
!++++++++++++++++++++++OUTPUT++++++++++++++++++++++++++++
!  TDB     Julian date 

IMPLICIT NONE

REAL(KIND=dp_ud) ::  ss,tdb,UTC1,UTC2,TAI1,TAI2,TT1,TT2
INTEGER(KIND=i4_ud) :: yr,doy,mon,day,hr,min,J,nday ,n   
    
INTEGER(KIND=i4_ud):: MTAB1(12)=(/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
INTEGER(KIND=i4_ud):: MTAB2(12)=(/ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

nday=0

! leap year
IF (((MOD(yr,4)==0) .AND. (MOD(yr,100)/=0)) .OR. (MOD(yr,400)==0)) THEN
    
    IF (doy>366) THEN
        
        WRITE(*,*) "bad day of year"
        
        STOP
        
    ELSE
    
        DO n=1,12
        
           nday=nday+MTAB2(n)
       
           IF(doy<=nday) EXIT
       
        END DO
    
        mon=n
    
        day=doy-(nday-MTAB2(n))
    
        IF (doy==nday) THEN
        
            mon=n
        
            day=MTAB2(n)
       
        END IF
    
    END IF

! not leap year    
ELSE  

    IF (doy>365) THEN
        
        WRITE(*,*) "bad day of year"
        
        STOP
        
    ELSE

        DO n=1,12
        
           nday=nday+MTAB1(n)
       
           IF(doy<=nday) EXIT
       
        END DO
    
        mon=n
    
        day=doy-(nday-MTAB1(n))
       
        IF (doy==nday) THEN
        
            mon=n
        
            day=MTAB1(n)
       
        END IF
        
    END IF
    
END IF

!convert UTC to TDB
CALL WR_DTF2D('UTC',yr,mon,day,hr,min,ss,UTC1,UTC2,J)

CALL WR_UTCTAI(UTC1,UTC2,TAI1,TAI2,J)

CALL WR_TAITT(TAI1,TAI2,TT1,TT2,J)

CALL WR_TTTDB(TT1,TT2,tdb)    
    
RETURN

END SUBROUTINE WR_UTC2TDB2    
!==================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
!==================================================================================



END MODULE MOD_TIME_TRANSFORMATION

