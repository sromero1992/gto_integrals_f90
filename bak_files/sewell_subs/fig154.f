C/ FIGURE 1.5.4
      SUBROUTINE DBAND(A,N,NLD,NUD,X,B)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C                              DECLARATIONS FOR ARGUMENTS
      DOUBLE PRECISION A(N,-NLD:NUD+NLD),X(N),B(N)
      INTEGER N,NLD,NUD
C                              DECLARATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION LJI
C
C  SUBROUTINE DBAND SOLVES THE LINEAR SYSTEM A*X=B, WHERE A IS A
C    BAND MATRIX.
C
C  ARGUMENTS
C
C             ON INPUT                          ON OUTPUT
C             --------                          ---------
C
C    A      - THE BAND MATRIX OF SIZE N,        DESTROYED.
C             DIMENSIONED A(N,-NLD:NUD+NLD)
C             IN THE MAIN PROGRAM.  COLUMNS
C             -NLD THROUGH NUD OF A CONTAIN
C             THE NONZERO DIAGONALS OF A.  THE
C             LAST NLD COLUMNS ARE USED AS
C             WORKSPACE (TO HOLD THE FILL-IN
C             IN THE NLD DIAGONALS DIRECTLY
C             ABOVE A).
C
C    N      - THE SIZE OF MATRIX A.
C
C    NLD    - NUMBER OF NONZERO LOWER DIAGONALS
C             IN A, I.E., NUMBER OF DIAGONALS
C             BELOW THE MAIN DIAGONAL.
C
C    NUD    - NUMBER OF NONZERO UPPER DIAGONALS
C             IN A, I.E., NUMBER OF DIAGONALS
C             ABOVE THE MAIN DIAGONAL.
C
C    X      -                                   AN N-VECTOR CONTAINING
C                                               THE SOLUTION.
C
C    B      - THE RIGHT HAND SIDE N-VECTOR.     DESTROYED.
C
C
C-----------------------------------------------------------------------
C                              ZERO TOP NLD DIAGONALS (WORKSPACE)
      DO 10 I=1,N
      DO 10 J=NUD+1,NUD+NLD
         A(I,J) = 0.0
   10 CONTINUE
C                              BEGIN FORWARD ELIMINATION
      DO 35 I=1,N-1
C                              SEARCH FROM AII ON DOWN FOR
C                              LARGEST POTENTIAL PIVOT, ALI
         BIG = ABS(A(I,0))
         L = I
         DO 15 J=I+1,MIN(I+NLD,N)
            IF (ABS(A(J,I-J)).GT.BIG) THEN
               BIG = ABS(A(J,I-J))
               L = J
            ENDIF
   15    CONTINUE
C                              IF LARGEST POTENTIAL PIVOT IS ZERO, 
C                              MATRIX IS SINGULAR
         IF (BIG.EQ.0.0) GO TO 50
C                              SWITCH ROW I WITH ROW L, TO BRING
C                              UP LARGEST PIVOT
         DO 20 K=I,MIN(I+NUD+NLD,N)
            TEMP = A(L,K-L)
            A(L,K-L) = A(I,K-I)
            A(I,K-I) = TEMP
   20    CONTINUE
C                              SWITCH B(I) AND B(L)
         TEMP = B(L)
         B(L) = B(I)
         B(I) = TEMP
         DO 30 J=I+1,MIN(I+NLD,N)
C                              CHOOSE MULTIPLIER TO ZERO AJI
            LJI = A(J,I-J)/A(I,0)
            IF (LJI.NE.0.0) THEN
C                              SUBTRACT LJI TIMES ROW I FROM ROW J
               DO 25 K=I,MIN(I+NUD+NLD,N)
                  A(J,K-J) = A(J,K-J) - LJI*A(I,K-I)
   25          CONTINUE
C                              SUBTRACT LJI TIMES B(I) FROM B(J)
               B(J) = B(J) - LJI*B(I)
            ENDIF
   30    CONTINUE
   35 CONTINUE
      IF (A(N,0).EQ.0.0) GO TO 50
C                              SOLVE U*X = B USING BACK SUBSTITUTION.
      X(N) = B(N)/A(N,0)
      DO 45 I=N-1,1,-1
         SUM = 0.0
         DO 40 J=I+1,MIN(I+NUD+NLD,N)
            SUM = SUM + A(I,J-I)*X(J)
   40    CONTINUE
         X(I) = (B(I)-SUM)/A(I,0)
   45 CONTINUE
      RETURN
C                              MATRIX IS NUMERICALLY SINGULAR.
   50 PRINT 55
   55 FORMAT (' ***** THE MATRIX IS SINGULAR *****')
      RETURN
      END
