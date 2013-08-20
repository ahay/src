Careful! Anything free comes with no guarantee.
      SUBROUTINE VSINTI(N,WSAVE)
C***BEGIN PROLOGUE  VSINTI
C***DATE WRITTEN   860701   (YYMMDD)
C***REVISION DATE  900509   (YYMMDD)
C***CATEGORY NO.  J1A3
c***KEYWORDS  FAST FOURIER TRANSFORM, SINE TRANSFORM, MULTIPLE
C             SEQUENCES
C***AUTHOR  BOISVERT, R. F. (NIST)
C***PURPOSE  Initialize for VSINT.
C***DESCRIPTION
C
C  Subroutine VSINTI initializes the array WSAVE which is used in
C  subroutine SINT.  The prime factorization of N together with
C  a tabulation of the trigonometric functions are computed and
C  stored in WSAVE.
C
C  Input Parameter
C
C  N       the length of the sequence to be transformed.  The method
C          is most efficient when N+1 is a product of small primes.
C
C  Output Parameter
C
C  WSAVE   a work array with at least INT(2.5*N+15) locations.
C          Different WSAVE arrays are required for different values
C          of N.  The contents of WSAVE must not be changed between
C          calls of VSINT.
C
C  -----------------------------------------------------------------
C
C  VSINTI is a straightforward extension of the subprogram SINTI to
C  handle M simultaneous sequences.  SINTI was originally developed
C  P. N. Swarztrauber of NCAR.
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C               Computations, (G. Rodrigue, ed.), Academic Press, 1982,
C               pp. 51-83.
C***ROUTINES CALLED  VRFFTI
C***END PROLOGUE  VSINTI
      DIMENSION       WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  SINTI
      PI = PIMACH(1.0)
      IF (N .LE. 1) RETURN
      NP1 = N+1
      NS2 = N/2
      DT = PI/REAL(NP1)
      KS = 1
      KF = KS+NS2-1
      FK = 0.
      DO 101 K=KS,KF
         FK = FK+1.
         WSAVE(K) = 2.*SIN(FK*DT)
  101 CONTINUE
      CALL VRFFTI (NP1,WSAVE(KF+1))
      RETURN
      END
