      SUBROUTINE FCTRL                                                  0000000 
      IMPLICIT REAL*8(A-H,O-Z)                                          0000000 
      COMMON /FACTRL/ FACLOG(300)                                       0000000 
      SQRPI=DSQRT(3.14159D0)                                            0000000 
      FACLOG(1)=0.d0                                                    0000000 
      DO 1 N=2,300                                                      0000000 
   1  FACLOG(N)=FACLOG(N-1)+DLOG(DFLOAT(N-1))                           0000000 
      RETURN                                                            0000000 
      END                                                               0000000 
