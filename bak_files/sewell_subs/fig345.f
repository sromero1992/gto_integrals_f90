C/ FIGURE 3.4.5
      SUBROUTINE LR(A,N,ERRLIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C                              DECLARATIONS FOR ARGUMENTS
      DOUBLE PRECISION A(N,N),ERRLIM
      INTEGER N
C                              DECLARATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION SAVE(N)
      IF (N.LE.2) RETURN
C                              USE LR ITERATION TO REDUCE HESSENBERG
C                              MATRIX A TO QUASI-TRIANGULAR FORM
      NITER = 1000*N
      DO 35 ITER=1,NITER
C                              REDUCE A TO UPPER TRIANGULAR FORM USING
C                              GAUSSIAN ELIMINATION (PREMULTIPLY BY
C                              Mij**(-1) MATRICES)
         DO 10 I=1,N-1
            IF (A(I+1,I).EQ.0.0) THEN
               R = 0.0
            ELSE
C                              IF PIVOTING NECESSARY, GIVE UP
               IF (ABS(A(I,I)).LT.ERRLIM) GO TO 40
               R = A(I+1,I)/A(I,I)
            ENDIF
C                              USE SAVE TO SAVE R FOR POST-
C                              MULTIPLICATION PHASE
            SAVE(I) = R
            IF (R.EQ.0.0) GO TO 10
C                              IF MATRIX TRIDIAGONAL, LIMITS ON K
C                              CAN BE:  K = I , I+1
            DO 5 K=I,N
               A(I+1,K) = A(I+1,K) - R*A(I,K)
    5       CONTINUE
   10    CONTINUE
C                              NOW POSTMULTIPLY BY Mij MATRICES
         DO 20 I=1,N-1
            R = SAVE(I)
            IF (R.EQ.0.0) GO TO 20
C                              IF MATRIX TRIDIAGONAL, LIMITS ON K
C                              CAN BE:  K = I , I+1
            DO 15 K=1,I+1
               A(K,I) = A(K,I) + R*A(K,I+1)
   15       CONTINUE
   20    CONTINUE
C                              SET NEARLY ZERO SUBDIAGONALS TO ZERO,
C                              TO AVOID UNDERFLOW.
         DO 25 I=1,N-1
            IF (ABS(A(I+1,I)).LT.ERRLIM) A(I+1,I) = 0.0
   25    CONTINUE
C                              CHECK FOR CONVERGENCE TO "QUASI-
C                              TRIANGULAR" FORM.
         ICONV = 1
         DO 30 I=2,N-1
            IF (A(I,I-1).NE.0.0 .AND. A(I+1,I).NE.0.0) ICONV = 0
   30    CONTINUE
         IF (ICONV.EQ.1) RETURN
   35 CONTINUE
C                              HAS NOT CONVERGED IN NITER ITERATIONS
   40 PRINT 45
   45 FORMAT (' ***** LR ITERATION DOES NOT CONVERGE *****')
      RETURN
      END
