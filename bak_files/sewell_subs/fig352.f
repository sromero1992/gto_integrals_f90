C/ FIGURE 3.5.2
      SUBROUTINE DPOWER(A,N,EIG,V,IUPDAT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C                              DECLARATIONS FOR ARGUMENTS
      COMPLEX*16 EIG,V(N)
      DOUBLE PRECISION A(N,N)
      INTEGER N,IUPDAT
C                              DECLARATIONS FOR LOCAL VARIABLES
      COMPLEX*16 VN(N),VNP1(N),B(N,N),PN,R,RNUM,RDEN
      INTEGER IPERM(N)
C
C  SUBROUTINE DPOWER FINDS ONE EIGENVALUE OF A, AND A CORRESPONDING
C    EIGENVECTOR, USING THE SHIFTED INVERSE POWER METHOD.
C
C  ARGUMENTS
C
C             ON INPUT                          ON OUTPUT
C             --------                          ---------
C
C    A      - THE N BY N MATRIX.
C
C    N      - THE SIZE OF MATRIX A.
C
C    EIG    - A (COMPLEX) INITIAL GUESS AT      AN EIGENVALUE OF A,
C             AN EIGENVALUE.                    NORMALLY THE ONE CLOSEST
C                                               TO THE INITIAL GUESS.
C
C    V      - A (COMPLEX) STARTING VECTOR       AN EIGENVECTOR OF A,
C             FOR THE SHIFTED INVERSE POWER     CORRESPONDING TO THE
C             METHOD.  IF ALL COMPONENTS OF     COMPUTED EIGENVALUE.
C             V ARE ZERO ON INPUT, A RANDOM
C             STARTING VECTOR WILL BE USED.
C
C    IUPDAT - THE NUMBER OF SHIFTED INVERSE
C             POWER ITERATIONS TO BE DONE
C             BETWEEN UPDATES OF P.  IF
C             IUPDAT=1, P WILL BE UPDATED EVERY
C             ITERATION.  IF IUPDAT > 1000,
C             P WILL NEVER BE UPDATED.
C
C-----------------------------------------------------------------------
C                              EPS = MACHINE FLOATING POINT RELATIVE
C                                    PRECISION
C *****************************
      DATA EPS/2.D-16/
C *****************************
      DO 5 I=1,N
         IF (V(I).NE.0.0) GO TO 15
    5 CONTINUE
C                             IF V = 0, GENERATE A RANDOM STARTING
C                             VECTOR
      SEED = N+10000
      DEN = 2.0**31-1.0
      DO 10 I=1,N
         SEED = MOD(7**5*SEED,DEN)
         V(I) = SEED/(DEN+1.0)
   10 CONTINUE
   15 CONTINUE
C                              NORMALIZE V, AND SET VN=V
      VNORM = 0.0
      DO 20 I=1,N
         VNORM = VNORM + ABS(V(I))**2
   20 CONTINUE
      VNORM = SQRT(VNORM)
      DO 25 I=1,N
         V(I) = V(I)/VNORM
         VN(I) = V(I)
   25 CONTINUE
C                              BEGIN SHIFTED INVERSE POWER ITERATION
      NITER = 1000
      DO 60 ITER=0,NITER
         IF (MOD(ITER,IUPDAT).EQ.0) THEN
C                              EVERY IUPDAT ITERATIONS, UPDATE PN
C                              AND SOLVE (A-PN*I)*VNP1 = VN
            PN = EIG
            DO 30 I=1,N
            DO 30 J=1,N
               IF (I.EQ.J) THEN
                  B(I,J) = A(I,J) - PN
               ELSE
                  B(I,J) = A(I,J)
               ENDIF
   30       CONTINUE
            CALL CLINEQ(B,N,VNP1,V,IPERM)
         ELSE
C                              BETWEEN UPDATES, WE CAN USE THE LU
C                              DECOMPOSITION OF B=A-PN*I CALCULATED
C                              EARLIER, TO SOLVE B*VNP1=VN FASTER
            CALL CRESLV(B,N,VNP1,V,IPERM)
         ENDIF
C                              CALCULATE NEW EIGENVALUE ESTIMATE,
C                              PN + (VN*VN)/(VN*VNP1)
         RNUM = 0.0
         RDEN = 0.0
         DO 35 I=1,N
            RNUM = RNUM + VN(I)*VN(I)
            RDEN = RDEN + VN(I)*VNP1(I)
   35    CONTINUE
         R = RNUM/RDEN
         EIG = PN + R
C                              SET V = NORMALIZED VNP1
         VNORM = 0.0
         DO 40 I=1,N
            VNORM = VNORM + ABS(VNP1(I))**2
   40    CONTINUE
         VNORM = SQRT(VNORM)
         DO 45 I=1,N
            V(I) = VNP1(I)/VNORM
   45    CONTINUE
C                              IF R*VNP1 = VN  (R = (VN*VN)/(VN*VNP1) ),
C                              ITERATION HAS CONVERGED.
         ERRMAX = 0.0
         DO 50 I=1,N
            ERRMAX = MAX(ERRMAX,ABS(R*VNP1(I)-VN(I)))
   50    CONTINUE
         IF (ERRMAX.LE.SQRT(EPS)) RETURN
C                              SET VN = V = NORMALIZED VNP1
         DO 55 I=1,N
            VN(I) = V(I)
   55    CONTINUE
   60 CONTINUE
      PRINT 65
   65 FORMAT (' ***** INVERSE POWER METHOD DOES NOT CONVERGE *****')
      RETURN
      END
 
      SUBROUTINE CLINEQ(A,N,X,B,IPERM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C                              DECLARATIONS FOR ARGUMENTS
      COMPLEX*16 A(N,N),X(N),B(N)
      INTEGER N,IPERM(N)
C                              DECLARATIONS FOR LOCAL VARIABLES
      COMPLEX*16 LJI,TEMP,SUM
C
C  SUBROUTINE CLINEQ SOLVES THE COMPLEX LINEAR SYSTEM A*X=B
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
 
      SUBROUTINE CRESLV(A,N,X,C,IPERM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C                              DECLARATIONS FOR ARGUMENTS
      COMPLEX*16 A(N,N),X(N),C(N)
      INTEGER N,IPERM(N)
C                              DECLARATIONS FOR LOCAL VARIABLES
      COMPLEX*16 LJI,SUM
C
C  SUBROUTINE CRESLV SOLVES THE COMPLEX LINEAR SYSTEM A*X=C IN O(N**2)
C    TIME, AFTER CLINEQ HAS PRODUCED AN LU DECOMPOSITION OF PA.
C
C  ARGUMENTS
C
C             ON INPUT                          ON OUTPUT
C             --------                          ---------
C
C    A      - THE N BY N COEFFICIENT MATRIX
C             AFTER PROCESSING BY CLINEQ.
C             AS OUTPUT BY CLINEQ, A CONTAINS
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
C             LENGTH N OUTPUT BY CLINEQ.
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
