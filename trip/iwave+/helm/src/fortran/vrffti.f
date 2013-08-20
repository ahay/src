Careful! Anything free comes with no guarantee.
      SUBROUTINE VRFFTI (N,WSAVE)
C***BEGIN PROLOGUE  VRFFTI
C***DATE WRITTEN   860701   (YYMMDD)
C***REVISION DATE  900509   (YYMMDD)
C***CATEGORY NO.  J1A1
C***KEYWORDS  FAST FOURIER TRANSFORM, REAL PERIODIC TRANSFORM,
C             MULTIPLE SEQUENCES
C***AUTHOR  SWEET, R.A. (NIST) AND LINDGREN, L.L. (NIST)
C***PURPOSE  Initialization for VRFFTF and VRFFTB.
C***DESCRIPTION
C
C  Subroutine VRFFTI initializes the array WSAVE which is used in
C  both VRFFTF and VRFFTB.  The prime factorization of N together with
C  a tabulation of certain trigonometric functions are computed and
C  stored in the array WSAVE.
C
C  Input Parameter
C
C  N       the length of the sequence to be transformed.  There is no
C          restriction on N.
C
C  Output Parameter
C
C  WSAVE   a work array which must be dimensioned at least N+15.
C          The same work array can be used for both VRFFTF and VRFFTB
C          as long as N remains unchanged.  Different WSAVE arrays
C          are required for different values of N.  The contents of
C          WSAVE must not be changed between calls of VRFFTF or VRFFTB.
C
C
C              * * * * * * * * * * * * * * * * * * * * *
C              *                                       *
C              *         PROGRAM SPECIFICATIONS        *
C              *                                       *
C              * * * * * * * * * * * * * * * * * * * * *
C
C
C     DIMENSION OF    R(MDIMR,N), RT(MDIMR,N), WSAVE(N+15)
C     ARGUMENTS
C
C     LATEST          AUGUST 1, 1985
C     REVISION
C
C     SUBPROGRAMS     VRFFTI, VRFTI1, VRFFTF, VRFTF1, VRADF2, VRADF3,
C     REQUIRED        VRADF4, VRADF5, VRADFG, VRFFTB, VRFTB1, VRADB2,
C                     VRADB3, VRADB4, VRADB5, VRADBG, PIMACH
C
C     SPECIAL         NONE
C     CONDITIONS
C
C     COMMON          NONE
C     BLOCKS
C
C     I/O             NONE
C
C     PRECISION       SINGLE
C
C     SPECIALIST      ROLAND SWEET
C
C     LANGUAGE        FORTRAN
C
C     HISTORY         WRITTEN BY LINDA LINDGREN AND ROLAND SWEET AT THE
C                     NATIONAL BUREAU OF STANDARDS (BOULDER).
C
C     ALGORITHM       A REAL VARIANT OF THE STOCKHAM AUTOSORT VERSION
C                     OF THE COOLEY-TUKEY FAST FOURIER TRANSFORM.
C
C     PORTABILITY     AMERICAN NATIONAL STANDARDS INSTITUTE FORTRAN 77.
C                     THE ONLY MACHINE DEPENDENT CONSTANT IS LOCATED IN
C                     THE FUNCTION PIMACH.
C
C     REQUIRED        COS,SIN
C     RESIDENT
C     ROUTINES
C
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C               Computations, (G. Rodrigue, ed.), Academic Press, 1982,
C               pp. 51-83.
C***ROUTINES CALLED  VRFTI1
C***END PROLOGUE  VRFFTI
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION       WSAVE(N+15)
C***FIRST EXECUTABLE STATEMENT  VRFFTI
      IF (N .EQ. 1) RETURN
      CALL VRFTI1 (N,WSAVE(1),WSAVE(N+1))
      RETURN
      END
