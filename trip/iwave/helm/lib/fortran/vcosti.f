Careful! Anything free comes with no guarantee.
      SUBROUTINE VCOSTI(N,WSAVE)
C***BEGIN PROLOGUE  VCOSTI
C***DATE WRITTEN   860701   (YYMMDD)
C***REVISION DATE  900509   (YYMMDD)
C***CATEGORY NO.  J1A3
C***KEYWORDS  FAST FOURIER TRANSFORM, COSINE TRANSFORM, MULTIPLE
C             SEQUENCES
C***AUTHOR  BOISVERT, R. F. (NIST)
C***PURPOSE  Initialize for VCOST.
C***DESCRIPTION
C
C  Subroutine VCOSTI initializes the array WSAVE which is used in
C  subroutine VCOST.  The prime factorization of N together with
C  a tabulation of the trigonometric functions are computed and
C  stored in WSAVE.
C
C  Input Parameter
C
C  N       the length of the sequence to be transformed.  The method
C          is most efficient when N-1 is a product of small primes.
C
C  Output Parameter
C
C  WSAVE   a work array which must be dimensioned at least 3*N+15.
C          Different WSAVE arrays are required for different values
C          of N.  The contents of WSAVE must not be changed between
C          calls of VCOST.
C
C  -----------------------------------------------------------------
C
C  VCOSTI is a straightforward extension of the subprogram COSTI to
C  handle M simultaneous sequences.  COSTI was originally developed
C  by P. N. Swarztrauber of NCAR.
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C               Computations, (G. Rodrigue, ed.), Academic Press, 1982,
C               pp. 51-83.
C***ROUTINES CALLED  VRFFTI
C***END PROLOGUE  VCOSTI
      DIMENSION       WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  VCOSTI
      PI = PIMACH(1.0)
      IF (N .LE. 3) RETURN
      NM1 = N-1
      NP1 = N+1
      NS2 = N/2
      DT = PI/REAL(NM1)
      FK = 0.
      DO 101 K=2,NS2
         FK = FK+1.
         WSAVE(K) = 2.*SIN(FK*DT)
  101 CONTINUE
      FK = 0.
      DO 102 K=2,NS2
         KC = NP1-K
         FK = FK+1.
         WSAVE(KC) = 2.*COS(FK*DT)
  102 CONTINUE
      CALL VRFFTI (NM1,WSAVE(N+1))
      RETURN
      END
