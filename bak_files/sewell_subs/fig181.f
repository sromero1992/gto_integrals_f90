C/ FIGURE 1.8.1
C                              JACOBI METHOD
      PARAMETER (M=10)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION UOLD(0:M,0:M,0:M),UNEW(0:M,0:M,0:M)
      H = 1.D0/M
C                              SET BOUNDARY KNOWNS TO ZERO PERMANENTLY
C                              AND INTERIOR UNKNOWNS TO ZERO TEMPORARILY
      DO 5 I=0,M
      DO 5 J=0,M
      DO 5 K=0,M
    5 UOLD(I,J,K) = 0.0
C                              BEGIN JACOBI ITERATION
      DO 25 ITER = 1,10000
C                              UPDATE UNKNOWNS ONLY
         DO 10 I=1,M-1
         DO 10 J=1,M-1
         DO 10 K=1,M-1
   10    UNEW(I,J,K) = H**2/6.0 + ( UOLD(I+1,J,K) + UOLD(I-1,J,K)
     &                            + UOLD(I,J+1,K) + UOLD(I,J-1,K)
     &                            + UOLD(I,J,K+1) + UOLD(I,J,K-1))/6.0
C                              COPY UNEW ONTO UOLD
         DO 15 I=1,M-1
         DO 15 J=1,M-1
         DO 15 K=1,M-1
   15    UOLD(I,J,K) = UNEW(I,J,K)
C                              EVERY 10 ITERATIONS CALCULATE MAXIMUM
C                              RESIDUAL AND CHECK FOR CONVERGENCE
         IF (MOD(ITER,10).NE.0) GO TO 25
         RMAX = 0.0
         DO 20 I=1,M-1
         DO 20 J=1,M-1
         DO 20 K=1,M-1
            RESID = 6*UOLD(I,J,K) - UOLD(I+1,J,K) - UOLD(I-1,J,K)
     &                            - UOLD(I,J+1,K) - UOLD(I,J-1,K)
     &                            - UOLD(I,J,K+1) - UOLD(I,J,K-1) - H**2
            RMAX = MAX(RMAX,ABS(RESID))
   20    CONTINUE
         RMAX = RMAX/H**2
         PRINT *, ITER, RMAX
         IF (RMAX.LE.1.D-10) STOP
   25 CONTINUE
      STOP
      END
