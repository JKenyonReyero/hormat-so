      FUNCTION VCC(JX1,JX2,JX3,MX1,MX2)                                 00000010
      IMPLICIT REAL*8(A-H,O-Z)                                          00000020
      DIMENSION FACT(299)                                               00000030
      COMMON/FACTRL/FACLOG(300)                                         00000040
      EQUIVALENCE(FACT(1),FACLOG(2))                                    00000050
      SQRT(XXX)=DSQRT(XXX)                                              00000060
      EXP(XXX)=DEXP(XXX)                                                00000070
      FLOAT(III)=DFLOAT(III)                                            00000080
      VCC=0.d0                                                          00000090
      J1=JX1                                                            00000100
      J2=JX2                                                            00000110
      J3=JX3                                                            00000120
      M1=MX1                                                            00000130
      M2=MX2                                                            00000140
      IF(J1.LT.J2) GO TO 20                                             00000150
      IF(J3.LT.J2) GO TO 30                                             00000160
      ICNTR=0                                                           00000170
      GO TO 40                                                          00000180
   20 IF(J3.LT.J1) GO TO 30                                             00000190
      ICNTR=-1                                                          00000200
      IT=J1                                                             00000210
      J1=J2                                                             00000220
      J2=IT                                                             00000230
      IT=M1                                                             00000240
      M1=M2                                                             00000250
      M2=IT                                                             00000260
      GO TO 40                                                          00000270
   30 ICNTR=1                                                           00000280
      IT=J2                                                             00000290
      J2=J3                                                             00000300
      J3=IT                                                             00000310
      M2=-M1-M2                                                         00000320
   40 CONTINUE                                                          00000330
      JZ1=(J1+J2-J3)/2                                                  00000340
      IF(JZ1.LT.0) GO TO 150                                            00000350
      JZ2=(J1+J3-J2)/2                                                  00000360
      IF(JZ2.LT.0) GO TO 150                                            00000370
      JZ3=(J2+J3-J1)/2                                                  00000380
      IF(JZ3.LT.0) GO TO 150                                            00000390
      IF(J1-IABS(M1).LT.0) GO TO 150                                    00000400
      IF(J2-IABS(M2).LT.0) GO TO 150                                    00000410
      IF(J3-IABS(M1+M2).LT.0) GO TO 150                                 00000420
      JT1=(J1-J3+M2)/2                                                  00000430
      JT2=(J2-J3-M1)/2                                                  00000440
      NUMIN=MAX0 (JT1,JT2,0)                                            00000450
      JT3=(J1-M1)/2                                                     00000460
      JT4=(J2+M2)/2                                                     00000470
      NUMAX=MIN0 (JT3,JT4,JZ1)                                          00000480
      JT5=(J2-M2)/2                                                     00000490
      IF(NUMAX.LT.NUMIN) GO TO 150                                      00000500
      J4=J1/2                                                           00000510
      J5=J3/2                                                           00000520
      PHAS=PHASEF(NUMIN)                                                00000530
      DO 100 NU=NUMIN,NUMAX                                             00000540
      VCT=FACT(J4)+FACT(J5)-FACT(JT3-NU)-FACT(NU-JT2)-FACT(JT4-NU)      00000550
     1    -FACT(NU-JT1)-FACT(JZ1-NU)-FACT(NU)                           00000560
      VCC=VCC+PHAS*EXP(VCT)                                             00000570
      PHAS=-PHAS                                                        00000580
  100 CONTINUE                                                          00000590
      FCTOT=FACT((J1+M1)/2)+FACT(JT3)+FACT(JZ2)+FACT((J3+M1+M2)/2)      00000600
     1   +FACT((J3-M1-M2)/2)+FACT(JZ1)+FACT(JZ3)+FACT(JT4)+FACT(JT5)    00000610
     2    -FACT(J4)-FACT(J4)-FACT((J1+J2+J3)/2+1)-FACT(J5)-FACT(J5)     00000620
      FCTOR=FLOAT(J3+1)*EXP(FCTOT)                                      00000630
      VCC=SQRT(FCTOR)*VCC                                               00000640
      IF(ICNTR)120,150,110                                              00000650
  110 VCC=VCC*SQRT(FLOAT(J2+1)/FLOAT(J3+1))*PHASEF(JT3)                 00000660
      GO TO 150                                                         00000670
  120 VCC=VCC*PHASEF(JZ1)                                               00000680
  150 RETURN                                                            00000690
      END                                                               00000700
      FUNCTION RACAH(J1,J2,J3,J4,J5,J6)                                 00000710
      IMPLICIT REAL*8(A-H,O-Z)                                          00000720
      DIMENSION FACT(299)                                               00000730
      COMMON/FACTRL/FACLOG(300)                                         00000740
      EQUIVALENCE(FACT(1),FACLOG(2))                                    00000750
      SQRT(XXX)=DSQRT(XXX)                                              00000760
      EXP(XXX)=DEXP(XXX)                                                00000770
      RACAH=0.0                                                         00000780
      Z1=DELR(J1,J2,J5)                                                 00000790
      IF(Z1.EQ.0.0) GO TO 90                                            00000800
      Z1=DELR(J3,J4,J5)*Z1                                              00000810
      IF(Z1.EQ.0.0) GO TO 90                                            00000820
      Z2=DELR(J1,J3,J6)                                                 00000830
      IF(Z2.EQ.0.0) GO TO 90                                            00000840
      Z2=DELR(J2,J4,J6)*Z2                                              00000850
      IF(Z2.EQ.0.0) GO TO 90                                            00000860
      Z1=SQRT(Z1/Z2)*Z2                                                 00000870
      JT1=(J1+J2+J5)/2                                                  00000880
      JT2=(J3+J4+J5)/2                                                  00000890
      JT3=(J1+J3+J6)/2                                                  00000900
      JT4=(J2+J4+J6)/2                                                  00000910
      JZ1=(J1+J2+J3+J4)/2                                               00000920
      JZ2=(J1+J4+J5+J6)/2                                               00000930
      JZ3=(J2+J3+J5+J6)/2                                               00000940
      NUMIN=MAX0 (JT1,JT2,JT3,JT4)                                      00000950
      NUMAX=MIN0 (JZ1,JZ2,JZ3)                                          00000960
      IF(NUMAX.LT.NUMIN) GO TO 90                                       00000970
      PHASE=PHASEF(NUMIN+JZ1)*Z1                                        00000980
      DO 80 NU=NUMIN,NUMAX                                              00000990
      JY1=NU-JT1                                                        00001000
      JY2=NU-JT2                                                        00001010
      JY3=NU-JT3                                                        00001020
      JY4=JZ1-NU                                                        00001030
      JY5=JZ2-NU                                                        00001040
      JY6=JZ3-NU                                                        00001050
      FCTOT=FACT(JY1)+FACT(JY2)+FACT(JY3)+FACT(JY4)+FACT(JY5)+FACT(JY6) 00001060
     1   +FACT(NU-JT4)-FACT(NU+1)                                       00001070
      FCTOR=EXP(FCTOT)                                                  00001080
      RACAH=RACAH+PHASE/FCTOR                                           00001090
      PHASE=-PHASE                                                      00001100
   80 CONTINUE                                                          00001110
   90 RETURN                                                            00001120
      END                                                               00001130
      FUNCTION DELR(J1,J2,J3)                                           00001140
      IMPLICIT REAL*8(A-H,O-Z)                                          00001150
      DIMENSION FACT(299)                                               00001160
      COMMON/FACTRL/FACLOG(300)                                         00001170
      EQUIVALENCE(FACT(1),FACLOG(2))                                    00001180
      EXP(XXX)=DEXP(XXX)                                                00001190
      JZ1=(J1+J2-J3)/2                                                  00001200
      IF(JZ1.LT.0) GO TO 130                                            00001210
      JZ2=(J1-J2+J3)/2                                                  00001220
      IF(JZ2.LT.0) GO TO 130                                            00001230
      JZ3=(J2+J3-J1)/2                                                  00001240
      IF(JZ3.LT.0) GO TO 130                                            00001250
      JZ4=(J1+J2+J3)/2+1                                                00001260
      DELT=FACT(JZ1)+FACT(JZ2)+FACT(JZ3)-FACT(JZ4)                      00001270
      DELR=EXP(DELT)                                                    00001280
      RETURN                                                            00001290
  130 DELR=0.0                                                          00001300
      RETURN                                                            00001310
      END                                                               00001320
      FUNCTION WINEJ(J1,J2,J3,J4,J5,J6,J7,J8,J9)                        00001330
      IMPLICIT REAL*8(A-H,O-Z)                                          00001340
      FLOAT(III)=DFLOAT(III)                                            00001350
      WINEJ=0.0                                                         00001360
      MUMIN=MAX0(IABS(J1-J9),IABS(J2-J6),IABS(J4-J8))                   00001370
      MUMAX=MIN0(J1+J9,J2+J6,J4+J8)                                     00001380
      IF(MUMAX.LT.MUMIN) GO TO 40                                       00001390
       DO 20 MU=MUMIN,MUMAX,2                                           00001400
      PROD=RACAH(J1,J4,J9,J8,J7,MU)*RACAH(J2,J5,MU,J4,J8,J6)*           00001410
     1     RACAH(J9,MU,J3,J2,J1,J6)*FLOAT(MU+1)                         00001420
      WINEJ=WINEJ+PROD                                                  00001430
   20 CONTINUE                                                          00001440
      WINEJ=WINEJ*PHASEF((J1+J3+J5+J8)/2+J2+J4+J9)                      00001450
   40 RETURN                                                            00001460
      END                                                               00001470
      FUNCTION PHASEF(N)                                                00001480
      IMPLICIT REAL*8(A-H,O-Z)                                          00001490
      FLOAT(III)=DFLOAT(III)                                            00001500
      PHASEF= FLOAT(1-2*IABS(N-2*(N/2)))                                00001510
      RETURN                                                            00001520
      END                                                               00001530
      FUNCTION YXFCT(M,N)                                               00001540
      IMPLICIT REAL*8(A-H,O-Z)                                          00001550
C     COMPUTES NFACT/MFACT                                              00001560
      DIMENSION FACT(299)                                               00001570
      COMMON/FACTRL/FACLOG(300)                                         00001580
      EQUIVALENCE(FACT(1),FACLOG(2))                                    00001590
      EXP(XXX)=DEXP(XXX)                                                00001600
      YXFCT=1.0                                                         00001610
      NUMAX=M-N                                                         00001620
      IF(NUMAX.EQ.0) GO TO 100                                          00001630
      FCTOR=FACT(N)-FACT(M)                                             00001640
      YXFCT=EXP(FCTOR)                                                  00001650
  100 CONTINUE                                                          00001660
      RETURN                                                            
      END                                                              
