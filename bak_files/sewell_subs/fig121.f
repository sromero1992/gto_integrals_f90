C/ FIGURE 1.2.1
      SUBROUTINE DLINEQ(A,N,X,B,IPERM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C                              DECLARATIONS FOR ARGUMENTS
      DOUBLE PRECISION A(N,N),X(N),B(N)
      INTEGER N,IPERM(N)
C                              DECLARATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION LJI
C
C  SUBROUTINE DLINEQ SOLVES THE LINEAR SYSTEM A*X=B
C
C  ARGUMENTS
C
C             ON INPUT                          ON OUTPUT
C             --------                          ---------
C
C    A      - THE N BY N COEFFICIENT MATRIX.    THE DIAGONAL AND UPPER
C                                               TRIANGLE OF A CONTAINS U
C                                               AND THE LOWER TRIANGLE
C                                               OF A CONTAINS THE LOWER
C                                               TRIANGLE OF L, WHERE
C                                               PA = LU, P BEING THE
C                                               PERMUTATION MATRIX
C                                               DEFINED BY IPERM.
C
C    N      - THE SIZE OF MATRIX A.
C
C    X      -                                   AN N-VECTOR CONTAINING
C                                               THE SOLUTION.
C
C    B      - THE RIGHT HAND SIDE N-VECTOR.     DESTROYED.
C
C    IPERM  -                                   AN N-VECTOR CONTAINING
C                                               A RECORD OF THE ROW
C                                               INTERCHANGES MADE.  IF
C                                               J = IPERM(K), THEN ROW
C                                               J ENDED UP AS THE K-TH
C                                               ROW.
C
C-----------------------------------------------------------------------
C                              INITIALIZE IPERM = (1,2,3,...,N)
      DO 10 K=1,N
         IPERM(K) = K
   10 CONTINUE
C                              BEGIN FORWARD ELIMINATION
      DO 35 I=1,N-1
C                              SEARCH FROM A(I,I) ON DOWN FOR
C                              LARGEST POTENTIAL PIVOT, A(L,I)
         BIG = ABS(A(I,I))
         L = I
         DO 15 J=I+1,N
            IF (ABS(A(J,I)).GT.BIG) THEN
               BIG = ABS(A(J,I))
               L = J
            ENDIF
   15    CONTINUE
C                              IF LARGEST POTENTIAL PIVOT IS ZERO, 
C                              MATRIX IS SINGULAR
         IF (BIG.EQ.0.0) GO TO 50
C                              SWITCH ROW I WITH ROW L, TO BRING
C                              UP LARGEST PIVOT
         DO 20 K=1,N
            TEMP = A(L,K)
            A(L,K) = A(I,K)
            A(I,K) = TEMP
   20    CONTINUE
C                              SWITCH B(I) AND B(L)
         TEMP = B(L)
         B(L) = B(I)
         B(I) = TEMP
C                              SWITCH IPERM(I) AND IPERM(L)
         ITEMP = IPERM(L)
         IPERM(L) = IPERM(I)
         IPERM(I) = ITEMP
         DO 30 J=I+1,N
C                              CHOOSE MULTIPLIER TO ZERO A(J,I)
            LJI = A(J,I)/A(I,I)
            IF (LJI.NE.0.0) THEN
C                              SUBTRACT LJI TIMES ROW I FROM ROW J
               DO 25 K=I+1,N
                  A(J,K) = A(J,K) - LJI*A(I,K)
   25          CONTINUE
C                              SUBTRACT LJI TIMES B(I) FROM B(J)
               B(J) = B(J) - LJI*B(I)
            ENDIF
C                              SAVE LJI IN A(J,I).  IT IS UNDERSTOOD,
C                              HOWEVER, THAT A(J,I) IS REALLY ZERO.
   30       A(J,I) = LJI
   35 CONTINUE
      IF (A(N,N).EQ.0.0) GO TO 50
C                              SOLVE U*X = B USING BACK SUBSTITUTION.
      X(N) = B(N)/A(N,N)
      DO 45 I=N-1,1,-1
         SUM = 0.0
         DO 40 J=I+1,N
            SUM = SUM + A(I,J)*X(J)
   40    CONTINUE
         X(I) = (B(I)-SUM)/A(I,I)
   45 CONTINUE
      RETURN
C                              MATRIX IS NUMERICALLY SINGULAR.
   50 PRINT 55
   55 FORMAT (' ***** THE MATRIX IS SINGULAR *****')
      RETURN
      END
