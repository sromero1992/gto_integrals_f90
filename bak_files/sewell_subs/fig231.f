C/ FIGURE 2.3.1
      SUBROUTINE REDH(A,M,N,B,PIVOT,NPIVOT,ERRLIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C                              DECLARATIONS FOR ARGUMENTS
      DOUBLE PRECISION A(M,N),B(M),ERRLIM
      INTEGER PIVOT(M),M,N,NPIVOT
C                              DECLARATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION W(M)
C                              USE HOUSEHOLDER TRANSFORMATIONS TO
C                              REDUCE A TO ROW ECHELON FORM
      I = 1
      DO 30 L=1,N
C                              USE PIVOT A(I,L) TO KNOCK OUT ELEMENTS
C                              I+1 TO M IN COLUMN L.
         IF (I+1.LE.M) THEN
C                              CHOOSE UNIT M-VECTOR W (WHOSE FIRST
C                              I-1 COMPONENTS ARE ZERO) SUCH THAT WHEN
C                              COLUMN L IS PREMULTIPLIED BY
C                              H = I - 2W*W**T, COMPONENTS I+1 THROUGH
C                              M ARE ZEROED.
            CALL CALW(A(1,L),M,W,I)
C                              PREMULTIPLY A BY H = I - 2W*W**T
            DO 15 K=L,N
               WTA = 0.0
               DO 5 J=I,M
                  WTA = WTA + W(J)*A(J,K)
    5          CONTINUE
               TWOWTA = 2*WTA
               DO 10 J=I,M
                  A(J,K) = A(J,K) - TWOWTA*W(J)
   10          CONTINUE
   15       CONTINUE
C                              PREMULTIPLY B BY H = I - 2W*W**T
            WTA = 0.0
            DO 20 J=I,M
               WTA = WTA + W(J)*B(J)
   20       CONTINUE
            TWOWTA = 2*WTA
            DO 25 J=I,M
               B(J) = B(J) - TWOWTA*W(J)
   25       CONTINUE
         ENDIF
C                              PIVOT A(I,L) IS NONZERO AFTER PROCESSING
C                              COLUMN L--MOVE DOWN TO NEXT ROW, I+1
         IF (ABS(A(I,L)).LE.ERRLIM) A(I,L) = 0.0
         IF (A(I,L).NE.0.0) THEN
            NPIVOT = I
            PIVOT(NPIVOT) = L
            I = I+1
            IF (I.GT.M) RETURN
         ENDIF
   30 CONTINUE
      RETURN
      END
 
      SUBROUTINE CALW(A,M,W,I)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION A(M),W(M)
C                              SUBROUTINE CALW CALCULATES A UNIT
C                              M-VECTOR W (WHOSE FIRST I-1 COMPONENTS
C                              ARE ZERO) SUCH THAT PREMULTIPLYING THE
C                              VECTOR A BY H = I - 2W*W**T ZEROES
C                              COMPONENTS I+1 THROUGH M.
      S = 0.0
      DO 5 J=I,M
         S = S + A(J)**2
         W(J) = A(J)
    5 CONTINUE
      IF (A(I) .GE. 0.0) THEN
         BETA = SQRT(S)
      ELSE
         BETA = -SQRT(S)
      ENDIF
      W(I) = A(I) + BETA
      TWOALP = SQRT(2*BETA*W(I))
C                              TWOALP=0 ONLY IF A(I),...,A(M) ARE ALL
C                              ZERO.  IN THIS CASE, RETURN WITH W=0
      IF (TWOALP.EQ.0.0) RETURN
C                              NORMALIZE W
      DO 10 J=I,M
         W(J) = W(J)/TWOALP
   10 CONTINUE
      RETURN
      END
