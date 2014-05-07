Careful! Anything free comes with no guarantee.
      SUBROUTINE VCOST(M,N,X,XT,MDIMX,WSAVE)
C***BEGIN PROLOGUE  VCOST
C***DATE WRITTEN   860701   (YYMMDD)
C***REVISION DATE  900509   (YYMMDD)
C***CATEGORY NO.  J1A3
C***KEYWORDS  FAST FOURIER TRANSFORM, COSINE TRANSFORM, MULTIPLE
C             SEQUENCES
C***AUTHOR  BOISVERT, R. F. (NIST)
C***PURPOSE  Cosine transform of one or more real, even sequences.
C***DESCRIPTION
C
C  Subroutine VCOST computes the discrete Fourier cosine transform
C  of M even sequences X(J,I), J=1,...,M.  The transform is defined
C  below at output parameter X.
C
C  The array WSAVE which is used by subroutine VCOST must be
C  initialized by calling subroutine VCOSTI(N,WSAVE).
C
C  Input Parameters
C
C  M       the number of sequences to be transformed.
C
C  N       the length of the sequence to be transformed.  N must be
C          greater than 1.  The method is most efficient when N-1 is
C          is a product of small primes.
C
C  X       an array of size at least X(MDIMX,N) which contains the
C          the sequences to be transformed.  The sequences are stored
C          in the ROWS of X.  Thus, the Jth sequence is stored in
C          X(J,I), I=1,..,N.
C
C  XT      a work array of size at least XT(MDIMX,N-1).
C
C  MDIMX   the first dimension of the array X exactly as it appears in
C          the calling program.
C
C  WSAVE   a work array which must be dimensioned at least 3*N+15
C          in the program that calls VCOST.  The WSAVE array must be
C          initialized by calling subroutine VCOSTI(N,WSAVE), and a
C          different WSAVE array must be used for each different
C          value of N.  This initialization does not have to be
C          repeated so long as N remains unchanged.  Thus subsequent
C          transforms can be obtained faster than the first.
C
C  Output Parameters
C
C  X       For I=1,...,N and J=1,...,M
C
C             X(J,I) = ( X(J,1)+(-1)**(I-1)*X(J,N)
C
C               + the sum from K=2 to K=N-1
C
C                 2*X(J,K)*COS((K-1)*(I-1)*PI/(N-1)) )/SQRT(2*(N-1))
C
C  WSAVE   contains initialization calculations which must not be
C          destroyed between calls of VCOST.
C
C  -----------------------------------------------------------------
C
C  NOTE  -  A call of VCOST followed immediately by another call
C           of VCOST will return the original sequences X.  Thus,
C           VCOST is the correctly normalized inverse of itself.
C
C  -----------------------------------------------------------------
C
C  VCOST is a straightforward extension of the subprogram COST to
C  handle M simultaneous sequences.  The scaling of the sequences
C  computed by VCOST is different than that of COST.  COST was
C  originally developed by P. N. Swarztrauber of NCAR.
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C               Computations, (G. Rodrigue, ed.), Academic Press, 1982,
C               pp. 51-83.
C***ROUTINES CALLED  VRFFTF
C***END PROLOGUE  VCOST
      DIMENSION       X(MDIMX,*), XT(MDIMX,*), WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  VCOST
      IF (M .LE. 0)  GO TO 900
      IF (N .LE. 1)  GO TO 900
      IF (N .GT. 3)  GO TO 400
      IF (N .EQ. 3)  GO TO 300
C
C  CASE  N = 2
C
      SCALE = SQRT(0.50E0)
      DO 210 J=1,M
         X1H = SCALE*(X(J,1)+X(J,2))
         X(J,2) = SCALE*(X(J,1)-X(J,2))
         X(J,1) = X1H
  210 CONTINUE
      GO TO 900
C
C  CASE  N = 3
C
  300 CONTINUE
      SCALE = 0.50E0
      DO 310 J=1,M
         X1P3 = X(J,1)+X(J,3)
         TX2 = X(J,2)+X(J,2)
         X(J,2) = SCALE*(X(J,1)-X(J,3))
         X(J,1) = SCALE*(X1P3+TX2)
         X(J,3) = SCALE*(X1P3-TX2)
  310 CONTINUE
      GO TO 900
C
C  CASE  N .GT. 3
C
C     ... PREPROCESSING
C
  400 CONTINUE
      NM1 = N-1
      NP1 = N+1
      NS2 = N/2
      DO 410 J=1,M
         XT(J,1) = X(J,1)-X(J,N)
         X(J,1) = X(J,1)+X(J,N)
  410 CONTINUE
      DO 420 K=2,NS2
         KC = NP1-K
         DO 420 J=1,M
            T1 = X(J,K)+X(J,KC)
            T2 = X(J,K)-X(J,KC)
            XT(J,1) = XT(J,1)+WSAVE(KC)*T2
            T2 = WSAVE(K)*T2
            X(J,K) = T1-T2
            X(J,KC) = T1+T2
  420 CONTINUE
      MODN = MOD(N,2)
      IF (MODN .NE. 0) THEN
         DO 430 J=1,M
            X(J,NS2+1) = X(J,NS2+1)+X(J,NS2+1)
  430    CONTINUE
      ENDIF
      DO 435 J=1,M
         X(J,N) = XT(J,1)
  435 CONTINUE
C
C     ... REAL PERIODIC TRANSFORM
C
      CALL VRFFTF (M,NM1,X,XT,MDIMX,WSAVE(NP1))
C
C     ... POSTPROCESSING
C
      FACTOR = 1.0/SQRT(REAL(NM1))
      DO 440 J=1,M
         XT(J,1) = X(J,2)
         X(J,2) = FACTOR*X(J,N)
  440 CONTINUE
      DO 450 I=4,N,2
         DO 450 J=1,M
            XI = X(J,I)
            X(J,I) = X(J,I-2)-X(J,I-1)
            X(J,I-1) = XT(J,1)
            XT(J,1) = XI
  450 CONTINUE
      IF (MODN .NE. 0) THEN
         DO 460 J=1,M
            X(J,N) = XT(J,1)
  460    CONTINUE
      ENDIF
C
C     ... NORMALIZATION
C
      SCALE = SQRT(0.5)
      DO 490 I=1,N
         DO 490 J=1,M
            X(J,I) = SCALE*X(J,I)
  490 CONTINUE
C
C  EXIT
C
  900 CONTINUE
      RETURN
      END
