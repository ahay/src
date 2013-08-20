Careful! Anything free comes with no guarantee.
      SUBROUTINE VRFFTF (M,N,R,RT,MDIMR,WSAVE)
C***BEGIN PROLOGUE  VRFFTF
C***DATE WRITTEN   850801   (YYMMDD)
C***REVISION DATE  900509   (YYMMDD)
C***CATEGORY NO.  J1A1
C***KEYWORDS  FAST FOURIER TRANSFORM, REAL PERIODIC TRANSFORM, 
C             FOURIER ANALYSIS, FORWARD TRANSFORM, MULTIPLE SEQUENCES
C***AUTHOR  SWEET, R.A. (NIST) AND LINDGREN, L.L. (NIST)
C***PURPOSE  Forward real periodic transform, M sequences.
C***DESCRIPTION
C
C  Subroutine VRFFTF computes the Fourier coefficients (forward 
C  transform) of a number of real periodic sequences.  Specifically,
C  for each sequence the subroutine claculates the independent
C  Fourier coefficients described below at output parameter R.
C
C  The array WSAVE which is used by subroutine VRFFTF must be
C  initialized by calling subroutine VRFFTI(N,WSAVE).
C
C
C  Input Parameters
C
C  M       the number of sequences to be transformed.
C
C  N       the length of the sequences to be transformed.  The method
C          is most efficient when N is a product of small primes,
C          however n may be any positive integer.
C
C  R       areal two-dimensional array of size MDIMX x N containing the
C          the sequences to be transformed.  The sequences are stored
C          in the ROWS of R.  Thus, the I-th sequence to be transformed,
C          X(I,J), J=0,1,...,N-1, is stored as
C
C               R(I,J) = X(I,J-1) , J=1, 2, . . . , N.
C
C  RT      a real two-dimensional work array of size MDIMX x N.
C
C  MDIMR   the row (or first) dimension of the arrays R and RT exactly 
C          as they appear in the calling program.  This parameter is 
C          used to specify the variable dimension of these arrays.
C
C  WSAVE   a real one-dimensional work array which must be dimensioned
C          at least N+15.  The WSAVE array must be initialized by 
C          calling subroutine VRFFTI.  A different WSAVE array must be
C          used for each different value of N.  This initialization does
C          not have to be repeated so long as N remains unchanged.  The
C          same WSAVE array may be used by VRFFTF and VRFFTB.
C
C  Output Parameters
C
C  R       contains the Fourier coefficients F(K) for each of the M 
C          input sequences.  Specifically, row I of R, R(I,J), 
C          J=1,2,..,N, contains the independent Fourier coefficients
C          F(I,K), for the I-th input sequence stored as
C
C             R(I,1) = REAL( F(I,0) ),
C                    = SQRT(1/N)*SUM(J=0,N-1)[ X(I,J) ],
C
C             R(I,2*K) = REAL( F(I,K) )
C                      = SQRT(1/N)*SUM(J=0,N-1)[X(I,J)*COS(2J*K*PI/N)]
C
C             R(I,2*K+1) = IMAG( F(I,K) )
C                        =-SQRT(1/N)*SUM(J=0,N-1)[X(I,J)*SIN(2J*K*PI/N)]
C
C                   for K = 1, 2, . . . , M-1,
C
C              and, when N is even,
C
C              R(I,N) = REAL( F(I,N/2) ).
C                     = SQRT(1/N)*SUM(J=0,N-1)[ (-1)**J*X(I,J) ].
C
C  WSAVE   contains results which must not be destroyed between calls
C          to VRFFTF or VRFFTB.
C
C  -----------------------------------------------------------------
C
C  NOTE  -  A call of VRFFTF followed immediately by a call of
C           of VRFFTB will return the original sequences R.  Thus,
C           VRFFTB is the correctly normalized inverse of VRFFTF.
C
C  -----------------------------------------------------------------
C
C  VRFFTF is a straightforward extension of the subprogram RFFTF to
C  handle M simultaneous sequences.  RFFTF was originally developed
C  by P. N. Swarztrauber of NCAR.
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
C***ROUTINES CALLED  VRFTF1
C***END PROLOGUE  VRFFTF
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION       R(MDIMR,N)  ,RT(MDIMR,N)    ,WSAVE(N+15)
C***FIRST EXECUTABLE STATEMENT  VRFFTF
      IF (N .EQ. 1) RETURN
      CALL VRFTF1 (M,N,R,RT,MDIMR,WSAVE(1),WSAVE(N+1))
      RETURN
      END
