      SUBROUTINE VRADF5 (MP,IDO,L1,CC,CH,MDIMC,WA1,WA2,WA3,WA4)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION  CC(MDIMC,IDO,L1,5)    ,CH(MDIMC,IDO,5,L1)     ,
     1           WA1(IDO)     ,WA2(IDO)     ,WA3(IDO)     ,WA4(IDO)
      ARG=2.*PIMACH(1.0)/5.
      TR11=COS(ARG)
      TI11=SIN(ARG)
      TR12=COS(2.*ARG)
      TI12=SIN(2.*ARG)
      DO 101 K=1,L1
         DO 1001 M=1,MP
         CH(M,1,1,K) = CC(M,1,K,1)+(CC(M,1,K,5)+CC(M,1,K,2))+
     1    (CC(M,1,K,4)+CC(M,1,K,3))
         CH(M,IDO,2,K) = CC(M,1,K,1)+TR11*(CC(M,1,K,5)+CC(M,1,K,2))+
     1    TR12*(CC(M,1,K,4)+CC(M,1,K,3))
         CH(M,1,3,K) = TI11*(CC(M,1,K,5)-CC(M,1,K,2))+TI12*
     1    (CC(M,1,K,4)-CC(M,1,K,3))
         CH(M,IDO,4,K) = CC(M,1,K,1)+TR12*(CC(M,1,K,5)+CC(M,1,K,2))+
     1    TR11*(CC(M,1,K,4)+CC(M,1,K,3))
         CH(M,1,5,K) = TI12*(CC(M,1,K,5)-CC(M,1,K,2))-TI11*
     1    (CC(M,1,K,4)-CC(M,1,K,3))
 1001    CONTINUE
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
            DO 1002 M=1,MP
            CH(M,I-1,1,K) = CC(M,I-1,K,1)+((WA1(I-2)*CC(M,I-1,K,2)+
     1       WA1(I-1)*CC(M,I,K,2))+(WA4(I-2)*CC(M,I-1,K,5)+WA4(I-1)*
     1       CC(M,I,K,5)))+((WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))+(WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*CC(M,I,K,4)))
            CH(M,I,1,K) = CC(M,I,K,1)+((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))+(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*
     1       CC(M,I-1,K,5)))+((WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4)))
            CH(M,I-1,3,K) = CC(M,I-1,K,1)+TR11*
     1      ( WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2)
     1       +WA4(I-2)*CC(M,I-1,K,5)+WA4(I-1)*CC(M,I,K,5))+TR12*
     1      ( WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3)
     1       +WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*CC(M,I,K,4))+TI11*
     1      ( WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2)
     1       -(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*CC(M,I-1,K,5)))+TI12*
     1      ( WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*CC(M,I-1,K,3)
     1       -(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*CC(M,I-1,K,4)))
            CH(M,IC-1,2,K) = CC(M,I-1,K,1)+TR11*
     1      ( WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2)
     1       +WA4(I-2)*CC(M,I-1,K,5)+WA4(I-1)*CC(M,I,K,5))+TR12*
     1     ( WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3)
     1      +WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*CC(M,I,K,4))-(TI11*
     1      ( WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2)
     1       -(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*CC(M,I-1,K,5)))+TI12*
     1      ( WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*CC(M,I-1,K,3)
     1       -(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*CC(M,I-1,K,4))))
            CH(M,I,3,K) = (CC(M,I,K,1)+TR11*((WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))+(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*
     1       CC(M,I-1,K,5)))+TR12*((WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))+(TI11*((WA4(I-2)*CC(M,I-1,K,5)+
     1       WA4(I-1)*CC(M,I,K,5))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2)))+TI12*((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4))-(WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))))
            CH(M,IC,2,K) = (TI11*((WA4(I-2)*CC(M,I-1,K,5)+WA4(I-1)*
     1       CC(M,I,K,5))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2)))+TI12*((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4))-(WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))))-(CC(M,I,K,1)+TR11*((WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))+(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*
     1       CC(M,I-1,K,5)))+TR12*((WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))
            CH(M,I-1,5,K) = (CC(M,I-1,K,1)+TR12*((WA1(I-2)*
     1       CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2))+(WA4(I-2)*
     1       CC(M,I-1,K,5)+WA4(I-1)*CC(M,I,K,5)))+TR11*((WA2(I-2)*
     1       CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3))+(WA3(I-2)*
     1       CC(M,I-1,K,4)+WA3(I-1)*CC(M,I,K,4))))+(TI12*((WA1(I-2)*
     1       CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2))-(WA4(I-2)*CC(M,I,K,5)-
     1       WA4(I-1)*CC(M,I-1,K,5)))-TI11*((WA2(I-2)*CC(M,I,K,3)-
     1       WA2(I-1)*CC(M,I-1,K,3))-(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))
            CH(M,IC-1,4,K) = (CC(M,I-1,K,1)+TR12*((WA1(I-2)*
     1       CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2))+(WA4(I-2)*
     1       CC(M,I-1,K,5)+WA4(I-1)*CC(M,I,K,5)))+TR11*((WA2(I-2)*
     1       CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3))+(WA3(I-2)*
     1       CC(M,I-1,K,4)+WA3(I-1)*CC(M,I,K,4))))-(TI12*((WA1(I-2)*
     1       CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2))-(WA4(I-2)*CC(M,I,K,5)-
     1       WA4(I-1)*CC(M,I-1,K,5)))-TI11*((WA2(I-2)*CC(M,I,K,3)-
     1       WA2(I-1)*CC(M,I-1,K,3))-(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))
            CH(M,I,5,K) = (CC(M,I,K,1)+TR12*((WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))+(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*
     1       CC(M,I-1,K,5)))+TR11*((WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))+(TI12*((WA4(I-2)*CC(M,I-1,K,5)+
     1       WA4(I-1)*CC(M,I,K,5))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2)))-TI11*((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4))-(WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))))
            CH(M,IC,4,K) = (TI12*((WA4(I-2)*CC(M,I-1,K,5)+WA4(I-1)*
     1       CC(M,I,K,5))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2)))-TI11*((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4))-(WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))))-(CC(M,I,K,1)+TR12*((WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))+(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*
     1       CC(M,I-1,K,5)))+TR11*((WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))
 1002       CONTINUE
  102    CONTINUE
  103 CONTINUE
      RETURN
      END
