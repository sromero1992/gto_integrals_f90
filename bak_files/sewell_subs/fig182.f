C/ FIGURE 1.8.2
C                              GAUSS-SEIDEL METHOD
      PARAMETER (M=10)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION U(0:M,0:M,0:M)
      H = 1.D0/M
C                              SET BOUNDARY KNOWNS TO ZERO PERMANENTLY
C                              AND INTERIOR UNKNOWNS TO ZERO TEMPORARILY
      DO 5 I=0,M
      DO 5 J=0,M
      DO 5 K=0,M
    5 U(I,J,K) = 0.0
C                              BEGIN GAUSS-SEIDEL ITERATION
      DO 20 ITER = 1,10000
C                              UPDATE UNKNOWNS ONLY
         DO 10 I=1,M-1
         DO 10 J=1,M-1
         DO 10 K=1,M-1
            GAUSS = H**2/6.0 + ( U(I+1,J,K) + U(I-1,J,K)
     &                         + U(I,J+1,K) + U(I,J-1,K)
     &                         + U(I,J,K+1) + U(I,J,K-1))/6.0
   10    U(I,J,K) = GAUSS
C                              EVERY 10 ITERATIONS CALCULATE MAXIMUM
C                              RESIDUAL AND CHECK FOR CONVERGENCE
         IF (MOD(ITER,10).NE.0) GO TO 20
         RMAX = 0.0
         DO 15 I=1,M-1
         DO 15 J=1,M-1
         DO 15 K=1,M-1
            RESID = 6*U(I,J,K) - U(I+1,J,K) - U(I-1,J,K)
     &                         - U(I,J+1,K) - U(I,J-1,K)
     &                         - U(I,J,K+1) - U(I,J,K-1) - H**2
            RMAX = MAX(RMAX,ABS(RESID))
   15    CONTINUE
         RMAX = RMAX/H**2
         PRINT *, ITER, RMAX
         IF (RMAX.LE.1.D-10) STOP
   20 CONTINUE
      STOP
      END
