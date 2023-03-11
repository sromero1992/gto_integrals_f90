C/ FIGURE 3.2.1
      SUBROUTINE DEGSYM(A,N,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C                              DECLARATIONS FOR ARGUMENTS
      DOUBLE PRECISION A(N,N),X(N,N)
      INTEGER N
C
C  SUBROUTINE DEGSYM SOLVES THE EIGENVALUE PROBLEM
C
C                A*X = LAMBDA*X
C
C    WHERE A IS A SYMMETRIC MATRIX.
C
C
C  ARGUMENTS
C
C             ON INPUT                          ON OUTPUT
C             --------                          ---------
C
C    A      - THE N BY N SYMMETRIC MATRIX.      A DIAGONAL MATRIX,
C                                               WITH THE EIGENVALUES
C                                               OF A ON THE DIAGONAL.
C
C    N      - THE SIZE OF MATRIX A.
C
C    X      -                                   AN N BY N MATRIX WHICH
C                                               CONTAINS THE EIGEN-
C                                               VECTORS OF A IN ITS
C                                               COLUMNS, IN THE SAME
C                                               ORDER AS THE EIGENVALUES
C                                               APPEAR ON THE DIAGONAL.
C
C-----------------------------------------------------------------------
C                              EPS = MACHINE FLOATING POINT RELATIVE
C                                    PRECISION
C *****************************
      DATA EPS/2.D-16/
C *****************************
C                              ANORM = SUM OF ALL SQUARES
C                              X INITIALIZED TO IDENTITY
      ANORM = 0.0
      DO 10 I=1,N
         DO 5 J=1,N
            ANORM = ANORM + A(I,J)**2
            X(I,J) = 0.0
    5    CONTINUE
         X(I,I) = 1.0
   10 CONTINUE
      ERRLIM = 1000*EPS*ANORM
C                              EK = SUM OF OFF-DIAGONAL SQUARES
      EK = 0.0
      DO 15 I=1,N
      DO 15 J=1,N
         IF (I .NE. J) EK = EK + A(I,J)**2
   15 CONTINUE
      IF (EK .LE. ERRLIM) RETURN
      THRESH = 0.5*EK/N/(N-1)
   20 CONTINUE
         DO 40 I=1,N-1
            DO 35 J=I+1,N
C                              IF A(J,I)**2 LESS THAN HALF THE
C                              AVERAGE FOR OFF-DIAGONALS, SKIP IT.
               IF (A(J,I)**2 .LE. THRESH) GO TO 35
C                              KNOCKING OUT A(J,I) WILL DECREASE OFF-
C                              DIAGONAL SUM OF SQUARES BY 2*A(J,I)**2.
               EK = EK - 2*A(J,I)**2
C                              CALCULATE NEW THRESHOLD.
               THRESH = 0.5*EK/N/(N-1)
C                              CALCULATE C,S
               BETA = (A(I,I)-A(J,J))/(2.*A(J,I))
               FRACT = 0.5*BETA/SQRT(1.0+BETA**2)
               S = SQRT(MAX(0.5-FRACT,0.D0))
               C = SQRT(MAX(0.5+FRACT,0.D0))
C                              PREMULTIPLY A BY Qij**T
               DO 25 K=1,N
                  PIK =  C*A(I,K)+S*A(J,K)
                  PJK = -S*A(I,K)+C*A(J,K)
                  A(I,K) = PIK
                  A(J,K) = PJK
   25          CONTINUE
C                              POSTMULTIPLY A AND X BY Qij
               DO 30 K=1,N
                  BKI =  C*A(K,I)+S*A(K,J)
                  BKJ = -S*A(K,I)+C*A(K,J)
                  A(K,I) = BKI
                  A(K,J) = BKJ
                  XKI =  C*X(K,I)+S*X(K,J)
                  XKJ = -S*X(K,I)+C*X(K,J)
                  X(K,I) = XKI
                  X(K,J) = XKJ
   30          CONTINUE
C                              CHECK FOR CONVERGENCE
               IF (EK .LE. ERRLIM) RETURN
   35       CONTINUE
   40    CONTINUE
C                              RETURN TO BEGINNING OF CYCLE
      GO TO 20
      END
