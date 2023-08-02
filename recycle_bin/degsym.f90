SUBROUTINE DEGSYM(A,N,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      !                              DECLARATIONS FOR ARGUMENTS
      DOUBLE PRECISION A(N,N),X(N,N)
      INTEGER N
      !  SUBROUTINE DEGSYM SOLVES THE EIGENVALUE PROBLEM
      !
      !                A*X = LAMBDA*X
      !
      !    WHERE A IS A SYMMETRIC MATRIX.
      !
      !
      !  ARGUMENTS
      !
      !             ON INPUT                          ON OUTPUT
      !             --------                          ---------
      !
      !    A      - THE N BY N SYMMETRIC MATRIX.      A DIAGONAL MATRIX,
      !                                               WITH THE EIGENVALUES
      !                                               OF A ON THE DIAGONAL.
      !
      !    N      - THE SIZE OF MATRIX A.
      !
      !    X      -                                   AN N BY N MATRIX WHICH
      !                                               CONTAINS THE EIGEN-
      !                                               VECTORS OF A IN ITS
      !                                               COLUMNS, IN THE SAME
      !                                               ORDER AS THE EIGENVALUES
      !                                               APPEAR ON THE DIAGONAL.
      !
      !-----------------------------------------------------------------------
      !                              EPS = MACHINE FLOATING POINT RELATIVE
      !                                    PRECISION
      ! *****************************
      DATA EPS/2.D-16/
      ! *****************************
      !                              ANORM = SUM OF ALL SQUARES
      !                              X INITIALIZED TO IDENTITY
      ANORM = 0.0
      DO 10 I=1,N
         DO 5 J=1,N
            ANORM = ANORM + A(I,J)**2
            X(I,J) = 0.0
    5    CONTINUE
         X(I,I) = 1.0
   10 CONTINUE
      ERRLIM = 1000*EPS*ANORM
!                              EK = SUM OF OFF-DIAGONAL SQUARES
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
!                              IF A(J,I)**2 LESS THAN HALF THE
!                              AVERAGE FOR OFF-DIAGONALS, SKIP IT.
               IF (A(J,I)**2 .LE. THRESH) GO TO 35
!                              KNOCKING OUT A(J,I) WILL DECREASE OFF-
!                              DIAGONAL SUM OF SQUARES BY 2*A(J,I)**2.
               EK = EK - 2*A(J,I)**2
!                              CALCULATE NEW THRESHOLD.
               THRESH = 0.5*EK/N/(N-1)
!                              CALCULATE C,S
               BETA = (A(I,I)-A(J,J))/(2.*A(J,I))
               FRACT = 0.5*BETA/SQRT(1.0+BETA**2)
               S = SQRT(MAX(0.5-FRACT,0.D0))
               C = SQRT(MAX(0.5+FRACT,0.D0))
!                              PREMULTIPLY A BY Qij**T
               DO 25 K=1,N
                  PIK =  C*A(I,K)+S*A(J,K)
                  PJK = -S*A(I,K)+C*A(J,K)
                  A(I,K) = PIK
                  A(J,K) = PJK
   25          CONTINUE
!                              POSTMULTIPLY A AND X BY Qij
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
!                              CHECK FOR CONVERGENCE
               IF (EK .LE. ERRLIM) RETURN
   35       CONTINUE
   40    CONTINUE
!                              RETURN TO BEGINNING OF CYCLE
      GO TO 20
END
