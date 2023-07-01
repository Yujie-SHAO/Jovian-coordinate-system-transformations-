# Jovian coordinate system transformations 

 This file is used to achieve the transformation between all Jovian coordinate system.  The **Universal Coordinated Time (UTC)** and **target vectors** are the input parameters.  DOI:10.5281/zenodo.8104077

 ***Before using these code, please pay more attention to announcements!***

 ## Jovian coordinate systems

 System III(1965) (S3RH)

 Jupiter Equatorial Inertial (JEI_J2000)

 Jupiter Solar Equatorial (JSE)

 Jupiter Solar Orbital (JSO)

 Jupiter Heliospheric (JH)

 Jupiter Solar Magnetic (JSM)

 Jupiter Magnetic (JM)

 Jupiter Solar MAGnetosphere (JSMAG) 

 ## Information of each file:

 **MOD_PRECISION.f90:** defines the precision of parameters and some common constants.

 **MOD_TIME_TRANSFORMATION.f90:** includes all revelant subroutines in SOFA to transform UTC to TDB.

 **MOD_CS_TRANSFORMATION.f90:** calculates all Jovian coordinate systems in JCRS and the final transformation matrix between any two coordinate systems.  In particular, we offer three Jovian magnetic field models (VIP4,JRM09, and JRM33) for magnetic field-related coordinate systems.

  **MOD_TOP2013.f90:** provides the relative motion between Jupiter and the Sun.

  **EXAPMLE1_S3RH2JSO.f90:** an example of transformation from S3RH to JSO.

  **EXAPMLE2_JM2JCRS.f90:** the other example of transformation from JM to JCRS.  
  
  **SOFA_announcements.txt:** includes the terms and conditions when using IAU SOFA and our codes.

  **TOP2013_announcements.txt:** includes the terms and conditions when using TOP2013.

  **example 1.txt:** output of *EXAPMLE1_S3RH2JSO.f90*; makes it easy for users to check.

   **example 2.txt:** output of *EXAPMLE2_JM2JCRS.f90*; makes it easy for users to check.

