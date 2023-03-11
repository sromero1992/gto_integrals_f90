C/ FIGURE 3.4.4
      SUBROUTINE HESSM(A,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C                              DECLARATIONS FOR ARGUMENTS
      DOUBLE PRECISION A(N,N)
      INTEGER N
      IF (N.LE.2) RETURN
C                              USE Mij TRANSFORMATIONS TO REDUCE A
C                              TO UPPER HESSENBERG FORM
      DO 35 I=2,N-1
C                              SEARCH FROM A(I,I-1) ON DOWN FOR
C                              LARGEST POTENTIAL PIVOT, A(L,I-1)
         BIG = ABS(A(I,I-1))
         L = I
         DO 5 J=I+1,N
            IF (ABS(A(J,I-1)).GT.BIG) THEN
               BIG = ABS(A(J,I-1))
               L = J
            ENDIF
    5    CONTINUE
C                              IF ALL SUBDIAGONAL ELEMENTS IN COLUMN
C                              I-1 ALREADY ZERO, GO ON TO NEXT COLUMN
         IF (BIG.EQ.0.0) GO TO 35
C                              PREMULTIPLY BY Pil
C                              (SWITCH ROWS I AND L)
         DO 10 K=I-1,N
            TEMP = A(L,K)
            A(L,K) = A(I,K)
            A(I,K) = TEMP
   10    CONTINUE
C                              POSTMULTIPLY BY Pil**(-1) = Pil
C                              (SWITCH COLUMNS I AND L)
         DO 15 K=1,N
            TEMP = A(K,L)
            A(K,L) = A(K,I)
            A(K,I) = TEMP
   15    CONTINUE
         DO 30 J=I+1,N
            R = A(J,I-1)/A(I,I-1)
            IF (R.EQ.0.0) GO TO 30
C                              PREMULTIPLY BY Mij**(-1)
C                              (SUBTRACT R TIMES ROW I FROM ROW J)
            DO 20 K=I-1,N
               A(J,K) = A(J,K) - R*A(I,K)
   20       CONTINUE
C                              POSTMULTIPLY BY Mij
C                              (ADD R TIMES COLUMN J TO COLUMN I)
            DO 25 K=1,N
               A(K,I) = A(K,I) + R*A(K,J)
   25       CONTINUE
   30    CONTINUE
   35 CONTINUE
      RETURN
      END
