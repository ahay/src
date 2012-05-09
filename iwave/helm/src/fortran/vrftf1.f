      SUBROUTINE VRFTF1 (M,N,C,CH,MDIMC,WA,FAC)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION       CH(MDIMC,N) ,C(MDIMC,N)  ,WA(N)   ,FAC(15)
      NF = FAC(2)
      NA = 1
      L2 = N
      IW = N
      DO 111 K1=1,NF
         KH = NF-K1
         IP = FAC(KH+3)
         L1 = L2/IP
         IDO = N/L2
         IDL1 = IDO*L1
         IW = IW-(IP-1)*IDO
         NA = 1-NA
         IF (IP .NE. 4) GO TO 102
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IF (NA .NE. 0) GO TO 101
         CALL VRADF4 (M,IDO,L1,C,CH,MDIMC,WA(IW),WA(IX2),WA(IX3))
         GO TO 110
  101    CALL VRADF4 (M,IDO,L1,CH,C,MDIMC,WA(IW),WA(IX2),WA(IX3))
         GO TO 110
  102    IF (IP .NE. 2) GO TO 104
         IF (NA .NE. 0) GO TO 103
         CALL VRADF2 (M,IDO,L1,C,CH,MDIMC,WA(IW))
         GO TO 110
  103    CALL VRADF2 (M,IDO,L1,CH,C,MDIMC,WA(IW))
         GO TO 110
  104    IF (IP .NE. 3) GO TO 106
         IX2 = IW+IDO
         IF (NA .NE. 0) GO TO 105
         CALL VRADF3 (M,IDO,L1,C,CH,MDIMC,WA(IW),WA(IX2))
         GO TO 110
  105    CALL VRADF3 (M,IDO,L1,CH,C,MDIMC,WA(IW),WA(IX2))
         GO TO 110
  106    IF (IP .NE. 5) GO TO 108
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IX4 = IX3+IDO
         IF (NA .NE. 0) GO TO 107
      CALL VRADF5(M,IDO,L1,C,CH,MDIMC,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 110
  107 CALL VRADF5 (M,IDO,L1,CH,C,MDIMC,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 110
  108    IF (IDO .EQ. 1) NA = 1-NA
         IF (NA .NE. 0) GO TO 109
         CALL VRADFG (M,IDO,IP,L1,IDL1,C,C,C,CH,CH,MDIMC,WA(IW))
         NA = 1
         GO TO 110
  109    CALL VRADFG (M,IDO,IP,L1,IDL1,CH,CH,CH,C,C,MDIMC,WA(IW))
         NA = 0
  110    L2 = L1
  111 CONTINUE
      SCALE=SQRT(1./N)
      IF (NA .EQ. 1) GO TO 113
      DO 112 J=1,N
      DO 112 I=1,M
         C(I,J) = SCALE*CH(I,J)
  112 CONTINUE
      RETURN
  113 DO 114 J=1,N
      DO 114 I=1,M
         C(I,J)=SCALE*C(I,J)
  114 CONTINUE
      RETURN
      END
