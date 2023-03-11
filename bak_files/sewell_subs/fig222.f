C/ FIGURE 2.2.2
      SUBROUTINE DLLSQR(A,M,N,X,B)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C                              DECLARATIONS FOR ARGUMENTS
      DOUBLE PRECISION A(M,N),X(N),B(M)
      INTEGER M,N
C                              DECLARATIONS FOR LOCAL VARIABLES
      INTEGER PIVOT(M)
C
C  SUBROUTINE DLLSQR SOLVES THE LINEAR LEAST SQUARES PROBLEM
C
C          MINIMIZE  2-NORM OF (A*X-B)
C
C
C  ARGUMENTS
C
C             ON INPUT                          ON OUTPUT
C             --------                          ---------
C
C    A      - THE M BY N MATRIX.                DESTROYED.
C
C    M      - THE NUMBER OF ROWS IN A.
C
C    N      - THE NUMBER OF COLUMNS IN A.
C
C    X      -                                   AN N-VECTOR CONTAINING
C                                               THE LEAST SQUARES
C                                               SOLUTION.
C
C    B      - THE RIGHT HAND SIDE M-VECTOR.     DESTROYED.
C
C-----------------------------------------------------------------------
C                              EPS = MACHINE FLOATING POINT RELATIVE
C                                    PRECISION
C *****************************
      DATA EPS/2.D-16/
C *****************************
C                              AMAX = MAXIMUM ELEMENT OF A
      AMAX = 0.0
      DO 5 I=1,M
      DO 5 J=1,N
    5 AMAX = MAX(AMAX,ABS(A(I,J)))
      ERRLIM = 1000*EPS*AMAX
C                              REDUCTION TO ROW ECHELON FORM
      CALL REDQ(A,M,N,B,PIVOT,NPIVOT,ERRLIM)
C                              CAUTION USER IF SOLUTION NOT UNIQUE.
      IF (NPIVOT.NE.N) THEN
         PRINT 10
   10    FORMAT (' NOTE: SOLUTION IS NOT UNIQUE ')
      ENDIF
C                              ASSIGN VALUE OF ZERO TO NON-PIVOT
C                              VARIABLES.
      DO 15 K=1,N
         X(K) = 0.0
   15 CONTINUE
C                              SOLVE FOR PIVOT VARIABLES USING BACK
C                              SUBSTITUTION.
      DO 25 I=NPIVOT,1,-1
         L = PIVOT(I)
         SUM = 0.0
         DO 20 K=L+1,N
            SUM = SUM + A(I,K)*X(K)
   20    CONTINUE
         X(L) = (B(I)-SUM)/A(I,L)
   25 CONTINUE
      RETURN
      END
 
      SUBROUTINE REDQ(A,M,N,B,PIVOT,NPIVOT,ERRLIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C                              DECLARATIONS FOR ARGUMENTS
      DOUBLE PRECISION A(M,N),B(M),ERRLIM
      INTEGER PIVOT(M),M,N,NPIVOT
C                              USE GIVENS ROTATIONS TO REDUCE A
C                              TO ROW ECHELON FORM
      I = 1
      DO 15 L=1,N
C                              USE PIVOT A(I,L) TO KNOCK OUT ELEMENTS
C                              I+1 TO M IN COLUMN L.
         DO 10 J=I+1,M
            IF (A(J,L).EQ.0.0) GO TO 10
            DEN = SQRT(A(I,L)**2+A(J,L)**2)
            C = A(I,L)/DEN
            S = A(J,L)/DEN
C                              PREMULTIPLY A BY Qij**T
            DO 5 K=L,N
               BIK = C*A(I,K) + S*A(J,K)
               BJK =-S*A(I,K) + C*A(J,K)
               A(I,K) = BIK
               A(J,K) = BJK
    5       CONTINUE
C                              PREMULTIPLY B BY Qij**T
            BI = C*B(I) + S*B(J)
            BJ =-S*B(I) + C*B(J)
            B(I) = BI
            B(J) = BJ
   10    CONTINUE
C                              PIVOT A(I,L) IS NONZERO AFTER PROCESSING
C                              COLUMN L--MOVE DOWN TO NEXT ROW, I+1
         IF (ABS(A(I,L)).LE.ERRLIM) A(I,L) = 0.0
         IF (A(I,L).NE.0.0) THEN
            NPIVOT = I
            PIVOT(NPIVOT) = L
            I = I+1
            IF (I.GT.M) RETURN
         ENDIF
   15 CONTINUE
      RETURN
      END
