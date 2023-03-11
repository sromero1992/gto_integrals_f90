C/ FIGURE 3.3.4
      SUBROUTINE DEGNON(A,N,EIG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C                              DECLARATIONS FOR ARGUMENTS
      COMPLEX*16 EIG(N)
      DOUBLE PRECISION A(N,N)
      INTEGER N
C
C  SUBROUTINE DEGNON SOLVES THE EIGENVALUE PROBLEM
C
C                A*X = LAMBDA*X
C
C    WHERE A IS A GENERAL REAL MATRIX.
C
C
C  ARGUMENTS
C
C             ON INPUT                          ON OUTPUT
C             --------                          ---------
C
C    A      - THE N BY N MATRIX.                DESTROYED.
C
C    N      - THE SIZE OF MATRIX A.
C
C    EIG    -                                   A COMPLEX N-VECTOR
C                                               CONTAINING THE EIGEN-
C                                               VALUES OF A.
C
C-----------------------------------------------------------------------
C                              EPS = MACHINE FLOATING POINT RELATIVE
C                                    PRECISION
C *****************************
      DATA EPS/2.D-16/
C *****************************
C                              AMAX = MAXIMUM ELEMENT OF A
      AMAX = 0.0
      DO 5 I=1,N
      DO 5 J=1,N
    5 AMAX = MAX(AMAX,ABS(A(I,J)))
      ERRLIM = SQRT(EPS)*AMAX
C                              REDUCTION TO HESSENBERG FORM
      CALL HESSQ(A,N)
C                              REDUCTION TO QUASI-TRIANGULAR FORM
      CALL QR(A,N,ERRLIM)
C                              EXTRACT EIGENVALUES OF QUASI-TRIANGULAR
C                              MATRIX
      I = 1
      DO WHILE (I.LE.N-1)
         IF (A(I+1,I).EQ.0.0) THEN
C                              1 BY 1 BLOCK ON DIAGONAL
            EIG(I) = A(I,I)
            I = I+1
         ELSE
C                              2 BY 2 BLOCK ON DIAGONAL
            DISC = (A(I,I)-A(I+1,I+1))**2 + 4.0*A(I,I+1)*A(I+1,I)
            TERM =  0.5*(A(I,I)+A(I+1,I+1))
            IF (DISC.GE.0.0) THEN
               EIG(I)  = TERM + 0.5*SQRT(DISC)
               EIG(I+1)= TERM - 0.5*SQRT(DISC)
            ELSE
               EIG(I)  = TERM + 0.5*SQRT(-DISC)*CMPLX(0.0,1.0)
               EIG(I+1)= TERM - 0.5*SQRT(-DISC)*CMPLX(0.0,1.0)
            ENDIF
            I = I+2
         ENDIF
      END DO
      IF (I.EQ.N) EIG(N) = A(N,N)
      RETURN
      END
 
      SUBROUTINE HESSQ(A,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C                              DECLARATIONS FOR ARGUMENTS
      DOUBLE PRECISION A(N,N)
      INTEGER N
      IF (N.LE.2) RETURN
C                              USE GIVENS ROTATIONS TO REDUCE A
C                              TO UPPER HESSENBERG FORM
      DO 20 I=2,N-1
         DO 15 J=I+1,N
            IF (A(J,I-1).EQ.0.0) GO TO 15
            DEN = SQRT(A(I,I-1)**2+A(J,I-1)**2)
            C = A(I,I-1)/DEN
            S = A(J,I-1)/DEN
C                              PREMULTIPLY BY Qij**T
            DO 5 K=I-1,N
               PIK = C*A(I,K) + S*A(J,K)
               PJK =-S*A(I,K) + C*A(J,K)
               A(I,K) = PIK
               A(J,K) = PJK
    5       CONTINUE
C                              POSTMULTIPLY BY Qij
            DO 10 K=1,N
               BKI = C*A(K,I) + S*A(K,J)
               BKJ =-S*A(K,I) + C*A(K,J)
               A(K,I) = BKI
               A(K,J) = BKJ
   10       CONTINUE
   15    CONTINUE
   20 CONTINUE
      RETURN
      END
 
      SUBROUTINE QR(A,N,ERRLIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C                              DECLARATIONS FOR ARGUMENTS
      DOUBLE PRECISION A(N,N),ERRLIM
      INTEGER N
C                              DECLARATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION SAVE(2,N)
      IF (N.LE.2) RETURN
C                              USE QR ITERATION TO REDUCE HESSENBERG
C                              MATRIX A TO QUASI-TRIANGULAR FORM
      NITER = 1000*N
      DO 35 ITER=1,NITER
C                              REDUCE A TO UPPER TRIANGULAR FORM USING
C                              ORTHOGONAL REDUCTION (PREMULTIPLY BY
C                              Qij**T MATRICES)
         DO 10 I=1,N-1
            IF (A(I+1,I).EQ.0.0) THEN
               C = 1.0
               S = 0.0
            ELSE
               DEN = SQRT(A(I,I)**2 + A(I+1,I)**2)
               C = A(I,I)/DEN
               S = A(I+1,I)/DEN
            ENDIF
C                              USE SAVE TO SAVE C,S FOR POST-
C                              MULTIPLICATION PHASE
            SAVE(1,I) = C
            SAVE(2,I) = S
            IF (S.EQ.0.0) GO TO 10
C                              IF MATRIX SYMMETRIC, LIMITS ON K
C                              CAN BE:  K = I , MIN(I+2,N)
            DO 5 K=I,N
               PIK = C*A(I,K) + S*A(I+1,K)
               PJK =-S*A(I,K) + C*A(I+1,K)
               A(I,K)   = PIK
               A(I+1,K) = PJK
    5       CONTINUE
   10    CONTINUE
C                              NOW POSTMULTIPLY BY Qij MATRICES
         DO 20 I=1,N-1
            C = SAVE(1,I)
            S = SAVE(2,I)
            IF (S.EQ.0.0) GO TO 20
C                              IF MATRIX SYMMETRIC, LIMITS ON K
C                              CAN BE:  K = MAX(1,I-1) , I+1
            DO 15 K=1,I+1
               BKI = C*A(K,I) + S*A(K,I+1)
               BKJ =-S*A(K,I) + C*A(K,I+1)
               A(K,I)   = BKI
               A(K,I+1) = BKJ
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
      PRINT 40
   40 FORMAT (' ***** QR ITERATION DOES NOT CONVERGE *****')
      RETURN
      END
