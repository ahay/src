Careful! Anything free comes with no guarantee.
      SUBROUTINE VSINT(M,N,X,XT,MDIMX,WSAVE)
C***BEGIN PROLOGUE  VSINT
C***DATE WRITTEN   860701   (YYMMDD)
C***REVISION DATE  900509   (YYMMDD)
C***CATEGORY NO.  J1A3
c***KEYWORDS  FAST FOURIER TRANSFORM, SINE TRANSFORM, MULTIPLE
C             SEQUENCES
C***AUTHOR  BOISVERT, R. F., (NIST)
C***PURPOSE  Sine transform of one or more real, odd sequences.
C***DESCRIPTION
C
C  Subroutine VSINT computes the discrete Fourier sine transform
C  of M odd sequences X(J,I), J=1,...,M.  The transform is defined
C  below at output parameter X.
C
C  The array WSAVE which is used by subroutine VSINT must be
C  initialized by calling subroutine VSINTI(N,WSAVE).
C
C  Input Parameters
C
C  M       the number of sequences to be transformed.
C
C  N       the length of the sequence to be transformed.  The method
C          is most efficient when N+1 is the product of small primes.
C
C  X       an array of size at least X(MDIMX,N+1) which contains the
C          the sequences to be transformed.  The sequences are stored
C          in the ROWS of X.  Thus, the Jth sequence is stored in
C          X(J,I), I=1,..,N.  The extra column of X is used as work
C          storage.
C
C  XT      a work array of size at least XT(MDIMX,N+1).
C
C  MDIMX   the first dimension of the array X exactly as it appears in
C          the calling program.
C
C  WSAVE   a work array with dimension at least INT(2.5*N+15)
C          in the program that calls VSINT.  The WSAVE array must be
C          initialized by calling subroutine VSINTI(N,WSAVE), and a
C          different WSAVE array must be used for each different
C          value of N.  This initialization does not have to be
C          repeated so long as N remains unchanged.
C
C  Output Parameters
C
C  X       for I=1,...,N and J=1,...,M
C
C               X(J,I)= the sum from k=1 to k=N
C
C                    2*X(J,K)*SIN(K*I*PI/(N+1))/SQRT(2*(N+1))
C
C  WSAVE   contains initialization calculations which must not be
C          destroyed between calls of VSINT.
C
C  -----------------------------------------------------------------
C
C  NOTE  -  A call of VSINT followed immediately by another call
C           of VSINT will return the original sequences X.  Thus,
C           VSINT is the correctly normalized inverse of itself.
C
C  -----------------------------------------------------------------
C
C  VSINT is a straightforward extension of the subprogram SINT to
C  handle M simultaneous sequences.  The scaling of the sequences
C  computed by VSINT is different than that of SINT.  SINT was
C  originally developed by P. N. Swarztrauber of NCAR.
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C               Computations, (G. Rodrigue, ed.), Academic Press, 1982,
C               pp. 51-83.
C***ROUTINES CALLED  VRFFTF
C***END PROLOGUE  VSINT
      DIMENSION       X(MDIMX,*), XT(MDIMX,*), WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  SINT
      IF (M .LE. 0)  GO TO 900
      IF (N .LE. 1)  GO TO 900
      IF (N .GT. 2)  GO TO 300
C
C  CASE   N = 2
C
      SQRTH = SQRT(0.50E0)
      DO 201 J=1,M
         XH = SQRTH*(X(J,1)+X(J,2))
         X(J,2) = SQRTH*(X(J,1)-X(J,2))
         X(J,1) = XH
  201 CONTINUE
      GO TO 900
C
C  CASE   N .GT. 2
C
C     ... PREPROCESSING
C
  300 CONTINUE
      NP1 = N+1
      NS2 = N/2
      DO 301 J=1,M
         XT(J,1) = 0.0
  301 CONTINUE
      DO 310 K=1,NS2
         KC = NP1-K
         DO 310 J=1,M
            T1 = X(J,K)-X(J,KC)
            T2 = WSAVE(K)*(X(J,K)+X(J,KC))
            XT(J,K+1) = T1+T2
            XT(J,KC+1) = T2-T1
  310 CONTINUE
      MODN = MOD(N,2)
      IF (MODN .NE. 0) THEN
         DO 320 J=1,M
            XT(J,NS2+2) = 4.0*X(J,NS2+1)
  320    CONTINUE
      ENDIF
C
C     ... REAL PERIODIC TRANSFORM
C
      NF = NS2+1
      CALL VRFFTF(M,NP1,XT,X,MDIMX,WSAVE(NF))
C
C     ... POSTPROCESSING
C
      DO 330 J=1,M
         X(J,1) = 0.5*XT(J,1)
  330 CONTINUE
      DO 350 I=3,N,2
         DO 340 J=1,M
            X(J,I-1) = -XT(J,I)
  340    CONTINUE
         DO 345 J=1,M
            X(J,I) = X(J,I-2)+XT(J,I-1)
  345    CONTINUE
  350 CONTINUE
      IF (MODN .EQ. 0) THEN
         DO 360 J=1,M
            X(J,N) = -XT(J,N+1)
  360    CONTINUE
      ENDIF
C
C     ... NORMALIZATION
C
      SCALE = SQRT(0.5)
      DO 370 I=1,N
         DO 370 J=1,M
            X(J,I) = SCALE*X(J,I)
  370 CONTINUE
C
C  EXIT
C
  900 CONTINUE
      RETURN
      END
