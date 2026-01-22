      FUNCTION TRB8(IN,IL,INN,ILL,LAM,IN1,IL1,IN2,IL2)
      implicit real*8(a-h,o-z)
C PROGRAM FOR CALCULATING OSCILLATOR BRACKETS
C THE MAIN ROUTINE OF A PROGRAM, ANY SUBROUTINE OF WHICH
C IS TO CALL TRB, MUST HAVE THE
C            CALL FACL(D)
C STATEMENT AT ITS VERY BEGINNING
      COMMON /PRAHA/ FL(400), GM(400), DG(50)
      COMMON /BRNO/A,TRM, DLN,D12,D21, Q1,Q2
      COMMON /REZ/N,L,NN,LL,N1,L1,N2,L2,IROO
C N,L (NN,LL) - RELATIVE (CENTRE-OF-MASS) QUANTUM NUMBERS
C N1,L1,N2,L2 PERTAIN TO INDIVIDUAL PARTICLES
C N,NN,N1,N2 = 0,1,2,...
C LAM IS TOTAL ANGULAR MOMENTUM,   IROO - ENERGY INDEX
C DLN, D12, D21 CONNECTED WITH MASS RATIO M1/M2
C SUM TO BE ACCUMULATED IN TRM,   PHASE FACTOR IN A
      INTEGER P1,P2,PD1,PD2,PH1,PH2,P1S,P2S
C STORE THE ARGUMENTS
      N = IN
      L = IL
      NN = INN
      LL = ILL
      N1 = IN1
      L1 = IL1
      N2 = IN2
      L2 = IL2
C TEST CONSERVATION LAWS, REARRANGE ARGUMENTS
      CALL SYMTR(LAM)
C ASSIGN TRB = 0. IF CONSERVATION LAWS VIOLATED
      IF (IROO) 800,100,100
  100 IF (L.NE.0) GO TO 150
C APPLY SIMPLIFIED FORMULAE
      CALL ZEROL
      GO TO 700
C
C                  GENERAL CASE  -  L.NE.0
C
  150 LS = L + 1
      DO 200 LA1S = 1,LS
      LA1 = LA1S - 1
      LA2 = L - LA1
      PD1 = IABS(L1-LA1) + 1
      IF ( (-1)**(PD1+LA1+L1).EQ.1)  PD1 = PD1 + 1
      PH1 = MIN0( L1+LA1, 2*N1+L1-LA1 ) + 1
      IF (PH1.LT.PD1) GO TO 200
      DO 210 P1S = PD1,PH1,2
      P1 = P1S - 1
      PD2 = MAX0( IABS(L2-LA2), IABS(LL-P1) ) + 1
      IF ( (-1)**(PD2+P1+LL).EQ.1) PD2 = PD2 + 1
      PH2 = MIN0( L2+LA2, 2*N2+L2-LA2, LL+P1 ) + 1
      IF (PH2.LT.PD2) GO TO 210
      MH1 = N1 - (LA1+P1-L1)/2 + 1
      DO 220 P2S = PD2,PH2,2
      P2 = P2S - 1
      MH2 = N2 - (LA2+P2-L2)/2 + 1
      A1 = 0.
      IF (MH1.LT.1) GO TO 220
      DO 230 M1S = 1,MH1
      M1 = M1S - 1
      A2 = 0.
      MD2 = MAX0( 0, NN - (P1+P2-LL)/2 - M1 ) + 1
      IF (MH2.LT.MD2) GO TO  230
      DO 240 M2S = MD2,MH2
      M2 = M2S - 1
      I1 = M1 + M2 + (P1 +P2 - LL)/2
      I3 = N2 - M2 - (LA2 + P2 - L2)/2
      I4 = I1 - NN
      I5 = (P1 + P2 + LL)/2 + M1 + M2 + 1
      I6 = M2 + P2 + 1
      T = FL(I1+1) - FL(M2S) - FL(I3+1) - FL(I4+1) + GM(I5) - GM(I6)
     1   - dFLOAT(M2) * Q2
  240 A2 = A2 + (-1)**M2 * dEXP(T)
      I2 = N1 - M1 - (LA1 + P1 - L1)/2
      I4 = M1 + P1 + 1
      T = GM(1) +dFLOAT(M1)*(Q1-Q2) - FL(M1S) - FL(I2+1) - GM(I4)
  230 A1 = A1 + (-1)**M1 * A2 * dEXP(T)
C SUMMATIONS OVER M1,M2 TERMINATED
C FUNCTION EIGHT CALCULATES THE S FACTOR WHICH IS INDEPENDENT
C                                                 OF M1 AND M2
      TRM = TRM + A1 * EIGHT(L1,P1,LA1,L2,P2,LA2,LAM,LL)
  220 CONTINUE
C SUMMATION OVER P2 TERMINATED
  210 CONTINUE
C SUMMATION OVER P1 TERMINATED
  200 CONTINUE
C SUMMATION OVER LA1 (AND LA2) TERMINATED
C
C CALCULATE THE NORMALIZATION FACTORS B0 AND B1
  700 CALL NORM(LAM,ALGB01)
      TRM = A * TRM * (-1.d0)**(LAM + N + N1 + N2) * dEXP(ALGB01)
  800 TRB8= TRM
      RETURN
      END
      SUBROUTINE SYMTR(LAM)
      implicit real*8(a-h,o-z)
      COMMON /BRNO/A,TRM, DLN,D12,D21, Q1,Q2
      COMMON /REZ/N,L,NN,LL,N1,L1,N2,L2,IROO
      TRM = 0.d0
      A = 1.d0
      Q1 = DLN
      Q2 = D12
      IR = 1
C TEST ENERGY AND PARITY CONSERVATION, TRIANGULAR INEQUALITIES,
C                       AND NON-NEGATIVITY OF ARGUMENTS
      IROO = 2 * (N + NN) + L + LL
      IRO = 2 * (N1 + N2) + L1 + L2
      IF (IROO.NE.IRO.OR.2*((L+LL+L1+L2)/2).NE.L+ LL+L1+L2) GO TO 800
      IF ( N.LT.0.OR. L.LT.0.OR.NN.LT.0.OR.LL.LT.0) GO TO 800
      IF (N1.LT.0.OR.L1.LT.0.OR.N2.LT.0.OR.L2.LT.0) GO TO 800
      IF (LAM.LT.IABS(L-LL).OR.LAM.GT.L+LL) GO TO 800
      IF (LAM.LT.IABS(L1-L2).OR.LAM.GT.L1+L2) GO TO 800
C   REARRANGE ARGUMENTS SO AS LNEW = MIN(L,LL,L1,L2)
      IF (LL.GE.L) GO TO 50
      I1 = N
      M1 = L
      N = NN
      L = LL
      NN = I1
      LL = M1
      A = A * (-1)**(L1-LAM)
      IR = - IR
   50 IF (L2.GE.L1) GO TO 55
      I1 = N1
      M1 = L1
      N1 = N2
      L1 = L2
      N2 = I1
      L2 = M1
      A = A * (-1)**(LL-LAM)
      IR = - IR
   55 IF (L1.GE.L) GO TO 60
      I1 = N1
      M1 = L1
      I2 = N2
      M2 = L2
      N1 = N
      L1 = L
      N2 = NN
      L2 = LL
      N = I1
      L = M1
      NN = I2
      LL = M2
      A = A * (-1)**(L2+LL)
   60 IF (IR) 70,80,80
   70 Q1 = - DLN
      Q2 = D21
   80 RETURN
  800 IROO = -10
      RETURN
      END
      SUBROUTINE ZEROL
      implicit real*8(a-h,o-z)
      COMMON /PRAHA/ FL(400), GM(400), DG(50)
      COMMON /BRNO/A,TRM, DLN,D12,D21, Q1,Q2
      COMMON /REZ/N,L,NN,LL,N1,L1,N2,L2,IROO
      IF (N.NE.0) GO TO 550
C L=0, N=0  THEN ALSO M1=N1, M2=N2
      I1 = NN + LL + 1
      I2 = N1 + L1 + 1
      I3 = N2 + L2 + 1
      I4 = (L1 + L2 + LL)/2
      I5 = (L1 + L2 - LL)/2
      I6 = (L1 + LL - L2)/2
      I7 = (L2 + LL - L1)/2
      M1 = L2 + LL - L1
      M2 = L1 + LL - L2
      M3 = L1 + L2 - LL
      M4 = L1 + L2 + LL + 1
      T = GM( 1) + GM(I1) + FL(NN+1) - FL(N1+1) - FL(N2+1)
     1  - GM(I2) - GM(I3) + FL(I4+1) - FL(I5+1) - FL(I6+1) - FL(I7+1)
     2  + 0.5d0*(FL(M1+1) + FL(M2+1) + FL(M3+1) - FL(M4+1) + DG(L1+1)
     3         + DG(L2+1) -dFLOAT(IROO)*Q2 +dFLOAT(2*N1+L1)*Q1  )
      TRM = (-1)**(NN + (L1 + L2 + LL)/2) * dEXP(T)
      RETURN
C L=0, N.NE.0 THEN P1=L1, P2=L2, LA1=0, LA2=0
  550 A1 = 0.
      MH1 = N1 + 1
      DO 560 M1S = 1,MH1
      M1 = M1S - 1
      A2 = 0.
      MD2 = MAX0( 0, NN-M1-(L1+L2-LL)/2 ) + 1
      MH2 = N2 + 1
      IF (MH2.LT.MD2) GO TO 560
      DO 570 M2S = MD2,MH2
      M2 = M2S - 1
      I1 = M1 + M2 + (L1 + L2 - LL)/2
      I2 = N2 - M2
      I3 = M1 + M2 - NN + (L1 + L2 - LL)/2
      I4 = (L1 + L2 + LL)/2 + M1 + M2 + 1
      I5 = M2 + L2 + 1
      T = FL(I1+1) - FL(M2S) - FL(I2+1) - FL(I3+1) + GM(I4) - GM(I5)
     1  -dFLOAT(M2) * Q2
  570 A2 = A2 + (-1)**M2 * dEXP(T)
      I1 = N1 - M1
      I2 = M1 + L1 + 1
      T = - FL(M1+1) - FL(I1+1) + GM(1) - GM(I2) -dFLOAT(M1)*(Q2-Q1)
  560 A1 = A1 + (-1)**M1 * A2 * dEXP(T)
      I1 = (L1 + L2 + LL)/2
      I2 = (L1 + L2 - LL)/2
      I3 = (L1 + LL - L2)/2
      I4 = (L2 + LL - L1)/2
      I5 = L1 + L2 - LL
      I6 = L1 + LL - L2
      I7 = L2 + LL - L1
      I8 = L1 + L2 + LL + 1
      T =      FL(I1+1) - FL(I2+1) - FL(I3+1) - FL(I4+1)
     1  + 0.5d0*(FL(I5+1) + FL(I6+1) + FL(I7+1) -dFLOAT(L1+L2)*Q2
     2       - FL(I8+1) + DG(L1+1) + DG(L2+1) +dFLOAT(L1)*Q1    )
      TRM = (-1)**(L2-L1) * A1 * dEXP(T)
      RETURN
      END
      FUNCTION EIGHT(J1,J2,J12,J3,J4,J34,J13,J24)
      implicit real*8(a-h,o-z)
C SUBROUTINE IS BASED ON A FORMULA DUE TO
C JUCYS AND BANDZAITIS FOR SPECIAL 9J SYMBOL
      COMMON /PRAHA/ FL(400), GM(400), DG(50)
      COMMON /BRNO/A,TRM, DLN,D12,D21, Q1,Q2
      R = 0.
      IXH = J2 + J4 - J24 + 1
      IF (IXH.LT.1) GO TO 30
      DO 10 IXS = 1,IXH
      IX = IXS - 1
      IYD = MAX0( 0, IX+J3+J12-J2-J13 ) + 1
      IYH = MIN0( IX+J3+J24-J2-J34, J1+J3-J13 ) + 1
      R1 = 0.
      IF (IYH.LT.IYD) GO TO 10
      DO 15 IYS = IYD,IYH
      IY = IYS - 1
      I1 = 2 * J3 - IY
      I2 = J1 - J3 + J13 + IY
      I3 = J1 + J3 - J13 - IY
      I4 = J2 + J13 - J3 - J12 + IY - IX
      I5 = J3 + J24 - J2 - J34 + IX - IY
      T = FL(I1+1) + FL(I2+1) - FL(I3+1) - FL(I4+1) - FL(I5+1) - FL(IYS)
   15 R1 = R1 + (-1)**IY * dEXP(T)
      I1 = 2 * J2 - IX
      I2 = J24 + J4 - J2 + IX
      I4 = J2 + J4 - J24 - IX
      T = FL(I1+1) + FL(I2+1) - FL(IXS) - FL(I4+1)
      R = R + (-1)**IX * R1 * dEXP(T)
   10 CONTINUE
      I1 = J1 + J2 - J12
      I2 = J1 + J2 + J12 + 1
      I3 = J2 + J4 - J24
      I4 = J2 + J4 + J24 + 1
      I5 = J3 + J4 - J34
      I6 = J3 + J4 + J34 + 1
      I7 = (J1 + J2 + J12)/2
      I8 = (J1 + J12 - J2)/2
      I9 = (J1 + J2 - J12)/2
      I0 = (J2 + J12 - J1)/2
      T1 = FL(I1+1) - FL(I2+1) + FL(I3+1) - FL(I4+1) + FL(I5+1)
     1   - FL(I6+1) + FL(I7+1) - FL(I8+1) - FL(I9+1) - FL(I0+1)
      I1 = (J3 + J4 + J34)/2
      I2 = (J3 + J4 - J34)/2
      I3 = (J3 + J34 - J4)/2
      I4 = (J4 + J34 - J3)/2
      I5 = (J2 + J4 + J24)/2
      I6 = (J2 + J4 - J24)/2
      I7 = (J2 + J24 - J4)/2
      I8 = (J4 + J24 - J2)/2
      T2 = FL(I1+1) - FL(I2+1) - FL(I3+1) - FL(I4+1) + FL(I5+1)
     1   - FL(I6+1) - FL(I7+1) - FL(I8+1) + DG(J2+1) + DG(J4+1)
      T3 = 0.5d0*(dFLOAT(J2+J34)*Q1 -dFLOAT(J2+J4)*Q2 )
      T = (-1)**(J34 + (J1+J2+J12)/2 + (J3+J4+J34)/2)
      EIGHT = R * T * dEXP(T1+T2+T3)
      RETURN
   30 EIGHT = 0.
      RETURN
      END
      SUBROUTINE NORM(LAM,T)
      implicit real*8(a-h,o-z)
      COMMON /PRAHA/ FL(400), GM(400), DG(50)
      COMMON /BRNO/A,TRM, DLN,D12,D21, Q1,Q2
      COMMON /REZ/N,L,NN,LL,N1,L1,N2,L2,IROO
      I1 = N1 + L1 + 1
      I2 = N2 + L2 + 1
      I3 =  N +  L + 1
      I4 = NN + LL + 1
      T = 0.5d0*(FL(N1+1)+FL(N2+1) + FL(N+1) - FL(NN+1) + GM(I1)
     1                 + GM(I2) - GM(I3) - GM(I4) ) + DG(L+1)
      IF (L.EQ.0) RETURN
      I1 = L1 + L2 - LAM
      I2 = L  + LL + LAM + 1
      I3 = L  + LL - LAM
      I4 = L  - LL + LAM
      I5 = LL -  L + LAM
      I6 = L1 + L2 + LAM + 1
      I7 = L1 - L2 + LAM
      I8 = L2 - L1 + LAM
      T =T+0.5d0*(FL(I1+1)+FL(I2+1) + FL(I3+1) + FL(I4+1) + FL(I5+1)
     1        - FL(I6+1) - FL(I7+1) - FL(I8+1) + DG(L1+1) + DG(L2+1)
     2        + DG(LL+1) - DG( L+1) -dFLOAT(L)*Q2 )
      RETURN
      END
      SUBROUTINE FACL(D)
      implicit real*8(a-h,o-z)
c     COMMON /PRAHA/ FL(200), GM(200), DG(25)
      COMMON /PRAHA/ FL(400), GM(400), DG(50)
      COMMON /BRNO/A,TRM, DLN,D12,D21, Q1,Q2
      MAXFL = 400
      MAXDG = 50
      MAXGM = 400
c     MAXFL = 200
c     MAXDG = 25
c     MAXGM = 200
C FL(I) = LN( (I-1)FACTORIAL ),                    I=1,2,...,MAXFL
      FL(1) = 0 
      FL(2) = 0 
      FN = 1 
      DO 10 I = 3,MAXFL
      FN = FN + 1 
   10 FL(I) = FL(I-1) + dLOG(FN)
c  10 FL(I) = FL(I-1) + ALOG(FN)
C DG(I) = LN(  2*(I-1)+1  ),                       I = 1,2,...,MAXDG
      DO 20 I = 1,MAXDG
      FN = 2 * I - 1
   20 DG(I) =  dLOG(FN)
c  20 DG(I) = ALOG(FN)
C GM(I) = LN( GAMMA(I+1/2)/GAMMA(3/2) ),           I=1,2,...,MAXGM
      GM(1) = 0 
      FN = 0.5d0
      DO 30 I = 2,MAXGM
      FN = FN + 1.d0
   30 GM(I) = GM(I-1) + dLOG(FN)
c  30 GM(I) = GM(I-1) + ALOG(FN)
C DLN=LN(MASS1/MASS2), D12=LN(1+MASS1/MASS2), D21=LN(1+MASS2/MASS1)
      DLN = dLOG(D)
      D12 = dLOG(1.d0+D)
      D21 = dLOG(1.d0+1.d0/D)
c     DLN = ALOG(D)
c     D12 = ALOG(1.+D)
c     D21 = ALOG(1.+1./D)
      IF (ABS(D-1.).GT..00000001) RETURN
C IMPROVE ACCURACY IF MASS1 = MASS2
      DLN = 0 
c     D12 = 0.6931471805
      d12=dlog(2.d0)
      D21 = D12
      RETURN
      END

