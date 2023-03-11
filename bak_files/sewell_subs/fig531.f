C/ FIGURE 5.3.1
      RECURSIVE SUBROUTINE DFFT(F,M)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C                              DECLARATIONS FOR ARGUMENTS
      COMPLEX*16 F(2**M)
      INTEGER M
C                              DECLARATIONS FOR LOCAL VARIABLES
      COMPLEX*16 FODD(2**M/2+1),FEVEN(2**M/2+1),D,H,U
C
C  SUBROUTINE DFFT PERFORMS A FAST FOURIER TRANSFORM ON THE COMPLEX
C     VECTOR F, OF LENGTH N=2**M.  THE FOURIER TRANSFORM IS DEFINED BY
C       Y(K) = SUM FROM J=1 TO N OF: EXP[I*2*PI*(K-1)*(J-1)/N]*F(J)
C     WHERE I = SQRT(-1).
C
C  ARGUMENTS
C
C             ON INPUT                          ON OUTPUT
C             --------                          ---------
C
C    F      - THE COMPLEX VECTOR OF LENGTH      THE TRANSFORMED VECTOR
C             2**M TO BE TRANSFORMED.           Y.
C
C    M      - THE LENGTH OF THE VECTOR F
C             IS ASSUMED TO BE 2**M.
C
C-----------------------------------------------------------------------
C                              FOURIER TRANSFORM OF A 1-VECTOR IS
C                              UNCHANGED
      IF (M.EQ.0) RETURN
      N = 2**M
      PI = 3.141592653589793 D0
      H = COS(2*PI/N) + SIN(2*PI/N)*CMPLX(0.0,1.0)
      N2 = N/2
C                              COPY ODD COMPONENTS OF F TO FODD
C                              AND EVEN COMPONENTS TO FEVEN
      DO 5 K=1,N2
         FODD(K) = F(2*K-1)
         FEVEN(K)= F(2*K)
    5 CONTINUE
C                              TRANSFORM N/2-VECTORS FODD AND FEVEN
C                              INTO YODD AND YEVEN
      CALL DFFT(FODD ,M-1)
      CALL DFFT(FEVEN,M-1)
      D = 1.0
C                              Y = (YODD+D*YEVEN , YODD-D*YEVEN)
      DO 10 K=1,N2
         U = D*FEVEN(K)
         F(K)    = FODD(K) + U
         F(N2+K) = FODD(K) - U
         D = D*H
   10 CONTINUE
      RETURN
      END
