C/ FIGURE 4.6.1
      SUBROUTINE DLPRV(DOTA,B,C,N,M,P,X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C                              DECLARATIONS FOR ARGUMENTS
      DOUBLE PRECISION B(M),C(N),P,X(N),Y(M)
      INTEGER N,M
      EXTERNAL DOTA
C                              DECLARATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION V(M),WK(M),XB(M),D(N),ABINV(M,M)
      INTEGER BASIS(M)
C
C  SUBROUTINE DLPRV USES THE REVISED SIMPLEX METHOD TO SOLVE THE PROBLEM
C
C             MAXIMIZE      P = C(1)*X(1) + ... + C(N)*X(N)
C
C    WITH X(1),...,X(N) NONNEGATIVE, AND
C
C           A(1,1)*X(1) + ... + A(1,N)*X(N)  =  B(1)
C              .                   .              .
C              .                   .              .
C           A(M,1)*X(1) + ... + A(M,N)*X(N)  =  B(M)
C
C    THE LAST M COLUMNS OF A MUST CONTAIN AN IDENTITY MATRIX, AND
C    B(1),...,B(M) MUST BE NONNEGATIVE.
C
C  ARGUMENTS
C
C             ON INPUT                          ON OUTPUT
C             --------                          ---------
C
C    DOTA   - NAME OF A USER-SUPPLIED FUNCTION:
C             DOTA(Z,J) SHOULD RETURN THE
C             DOT PRODUCT OF THE M-VECTOR Z
C             WITH COLUMN J OF A, FOR J .LE. N-M
C             (IDENTITY MATRIX ASSUMED IN LAST 
C             M COLUMNS). 
C
C    B      - A VECTOR OF LENGTH M CONTAINING
C             THE RIGHT HAND SIDES OF THE
C             CONSTRAINTS.  THE COMPONENTS OF
C             B MUST ALL BE NONNEGATIVE.
C
C    C      - A VECTOR OF LENGTH N CONTAINING
C             THE COEFFICIENTS OF THE OBJECTIVE
C             FUNCTION.
C
C    N      - THE NUMBER OF UNKNOWNS (N>M).
C
C    M      - THE NUMBER OF CONSTRAINTS.
C
C    P      -                                   THE MAXIMUM OF THE
C                                               OBJECTIVE FUNCTION.
C
C    X      -                                   A VECTOR OF LENGTH N
C                                               WHICH CONTAINS THE LP
C                                               SOLUTION.
C
C    Y      -                                   A VECTOR OF LENGTH M
C                                               WHICH CONTAINS THE DUAL
C                                               SOLUTION.
C
C-----------------------------------------------------------------------
C                              EPS = MACHINE FLOATING POINT RELATIVE
C                                    PRECISION
C                              BIGNO = A VERY LARGE NUMBER
C *****************************
      DATA EPS,BIGNO/2.D-16,1.D35/
C *****************************
C                              INITIALIZE Ab**(-1) TO IDENTITY
      DO 5 I=1,M
      DO 5 J=1,M
         ABINV(I,J) = 0.0
         IF (I.EQ.J) ABINV(I,J) = 1.0
    5 CONTINUE
C                              BASIS(1),...,BASIS(M) HOLD NUMBERS OF
C                              BASIS VARIABLES
      DO 10 I=1,M
C                              INITIAL BASIS = LAST M VARIABLES
         K = N-M+I
         BASIS(I) = K
C                              INITIALIZE Y TO Ab**(-T)*Cb = Cb
         Y(I) = C(K)
C                              INITIALIZE Xb TO Ab**(-1)*B = B
         XB(I) = B(I)
   10 CONTINUE
C                              THRESH = SMALL NUMBER.  WE ASSUME SCALES
C                              OF A AND C ARE NOT *TOO* DIFFERENT
      THRESH = 0.0
      DO 15 J=1,N
         THRESH = MAX(THRESH,ABS(C(J)))
   15 CONTINUE
      THRESH = 1000*EPS*THRESH
C                              BEGINNING OF SIMPLEX STEP
   20 CONTINUE
C                              D**T = Y**T*A - C**T
         DO 25 J=1,N
            IF (J.LE.N-M) THEN
               D(J) = DOTA(Y,J) - C(J)
            ELSE
               D(J) = Y(J-(N-M)) - C(J)
            ENDIF
   25    CONTINUE
C                              FIND MOST NEGATIVE ENTRY IN D,
C                              IDENTIFYING PIVOT COLUMN JP
         CMIN = -THRESH
         JP = 0
         DO 30 J=1,N
            IF (D(J).LT.CMIN) THEN
               CMIN = D(J)
               JP = J
            ENDIF
   30    CONTINUE
C                              IF ALL ENTRIES NONNEGATIVE (ACTUALLY,
C                              IF GREATER THAN -THRESH) WE ARE THROUGH
         IF (JP.EQ.0) GO TO 80 
C                              V = Ab**(-1)*Ajp (Ajp = COLUMN JP OF A)
         DO 40 I=1,M
            DO 35 J=1,M
   35       WK(J) = ABINV(I,J)
            IF (JP.LE.N-M) THEN
               V(I) = DOTA(WK,JP)
            ELSE
               V(I) = WK(JP-(N-M))
            ENDIF
   40    CONTINUE
C                              FIND SMALLEST POSITIVE RATIO
C                              Xb(I)/V(I), IDENTIFYING PIVOT ROW IP
         RATMIN = BIGNO
         IP = 0
         DO 45 I=1,M
            IF (V(I).GT.THRESH) THEN
               RATIO = XB(I)/V(I)
               IF (RATIO.LT.RATMIN) THEN
                  RATMIN = RATIO
                  IP = I
               ENDIF
            ENDIF
   45    CONTINUE
C                              IF ALL RATIOS NONPOSITIVE, MAXIMUM
C                              IS UNBOUNDED
         IF (IP.EQ.0) THEN
            PRINT 50
   50       FORMAT (' ***** UNBOUNDED MAXIMUM *****')
            RETURN
         ENDIF
C                              ADD X(JP) TO BASIS
         BASIS(IP) = JP
C                              UPDATE Ab**(-1) = E**(-1)*Ab**(-1)
C                                     Xb = E**(-1)*Xb
         DO 55 J=1,M
   55    ABINV(IP,J) = ABINV(IP,J)/V(IP)
         XB(IP) = XB(IP)/V(IP)
         DO 65 I=1,M
            IF (I.EQ.IP) GO TO 65
            DO 60 J=1,M
               ABINV(I,J) = ABINV(I,J) - V(I)*ABINV(IP,J)
   60       CONTINUE
            XB(I) = XB(I) - V(I)*XB(IP)
   65    CONTINUE
C                              CALCULATE Y = Ab**(-T)*Cb
         DO 75 I=1,M
            Y(I) = 0.0
            DO 70 J=1,M
               K = BASIS(J)
               Y(I) = Y(I) + ABINV(J,I)*C(K)
   70       CONTINUE
   75    CONTINUE
      GO TO 20
C                              END OF SIMPLEX STEP
   80 CONTINUE
C                              CALCULATE X
      DO 85 J=1,N
         X(J) = 0.0
   85 CONTINUE
      DO 90 I=1,M
         K = BASIS(I)
         X(K) = XB(I)
   90 CONTINUE
C                              CALCULATE P
      P = 0.0
      DO 95 I=1,N
         P = P + C(I)*X(I)
   95 CONTINUE
      RETURN
      END
