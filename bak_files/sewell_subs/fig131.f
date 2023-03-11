C/ FIGURE 1.3.1
      SUBROUTINE DRESLV(A,N,X,C,IPERM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C                              DECLARATIONS FOR ARGUMENTS
      DOUBLE PRECISION A(N,N),X(N),C(N)
      INTEGER N,IPERM(N)
C                              DECLARATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION LJI
C
C  SUBROUTINE DRESLV SOLVES THE LINEAR SYSTEM A*X=C IN O(N**2) TIME,
C    AFTER DLINEQ HAS PRODUCED AN LU DECOMPOSITION OF PA.
C
C  ARGUMENTS
C
C             ON INPUT                          ON OUTPUT
C             --------                          ---------
C
C    A      - THE N BY N COEFFICIENT MATRIX
C             AFTER PROCESSING BY DLINEQ.
C             AS OUTPUT BY DLINEQ, A CONTAINS
C             AN LU DECOMPOSITION OF PA.
C
C    N      - THE SIZE OF MATRIX A.
C
C    X      -                                   AN N-VECTOR CONTAINING
C                                               THE SOLUTION.
C
C    C      - THE RIGHT HAND SIDE N-VECTOR.     DESTROYED.
C
C    IPERM  - THE PERMUTATION VECTOR OF
C             LENGTH N OUTPUT BY DLINEQ.
C
C-----------------------------------------------------------------------
C                              CALCULATE C=P*C, WHERE P IS PERMUTATION
C                              MATRIX DEFINED BY IPERM.
      DO 5 K=1,N
         J = IPERM(K)
         X(K) = C(J)
    5 CONTINUE
      DO 10 K=1,N
         C(K) = X(K)
   10 CONTINUE
C                              BEGIN FORWARD ELIMINATION, TO CALCULATE
C                              C = L**(-1)*C
      DO 20 I=1,N-1
         DO 15 J=I+1,N
C                              RETRIEVE MULTIPLIER SAVED IN A(J,I)
            LJI = A(J,I)
C                              SUBTRACT LJI TIMES C(I) FROM C(J)
            C(J) = C(J) - LJI*C(I)
   15    CONTINUE
   20 CONTINUE
C                              SOLVE U*X = C USING BACK SUBSTITUTION.
      X(N) = C(N)/A(N,N)
      DO 30 I=N-1,1,-1
         SUM = 0.0
         DO 25 J=I+1,N
            SUM = SUM + A(I,J)*X(J)
   25    CONTINUE
         X(I) = (C(I)-SUM)/A(I,I)
   30 CONTINUE
      RETURN
      END
