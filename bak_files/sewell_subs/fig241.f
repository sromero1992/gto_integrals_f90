C/ FIGURE 2.4.1
      SUBROUTINE DLSQSP(X,N,XD,YD,M,XOUT,YOUT,NOUT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C                              DECLARATIONS FOR ARGUMENTS
      DOUBLE PRECISION X(N),XD(M),YD(M),XOUT(NOUT),YOUT(NOUT)
      INTEGER N,M,NOUT
C                              DECLARATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION Y(N),A(M,N)
C
C  SUBROUTINE DLSQSP CALCULATES A NATURAL CUBIC SPLINE WITH KNOTS AT
C    X(1),...,X(N) WHICH IS THE LEAST SQUARES FIT TO THE DATA POINTS
C    (XD(I),YD(I)), I=1,...,M, AND EVALUATES THIS SPLINE AT THE OUTPUT
C    POINTS XOUT(1),...,XOUT(NOUT).
C
C  ARGUMENTS
C
C             ON INPUT                          ON OUTPUT
C             --------                          ---------
C
C    X      - A VECTOR OF LENGTH N CONTAINING
C             THE SPLINE KNOTS.
C
C    N      - THE NUMBER OF KNOTS.
C             (N.GE.3).
C
C    XD     - A VECTOR OF LENGTH M CONTAINING
C             THE X-COORDINATES OF THE DATA
C             POINTS.
C
C    YD     - A VECTOR OF LENGTH M CONTAINING   DESTROYED.
C             THE Y-COORDINATES OF THE DATA
C             POINTS.
C
C    M      - THE NUMBER OF DATA POINTS.
C
C    XOUT   - A VECTOR OF LENGTH NOUT CONTAINING
C             THE X-COORDINATES AT WHICH THE
C             CUBIC SPLINE IS EVALUATED.  THE
C             ELEMENTS OF XOUT MUST BE IN
C             ASCENDING ORDER.
C
C    YOUT   -                                   A VECTOR OF LENGTH NOUT.
C                                               YOUT(I) CONTAINS THE
C                                               VALUE OF THE SPLINE
C                                               AT XOUT(I).
C
C    NOUT   - THE NUMBER OF OUTPUT POINTS.
C
C-----------------------------------------------------------------------
      ZERO = 0.0D0
      DO 5 J=1,N
         Y(J) = 0.0
    5 CONTINUE
      DO 10 J=1,N
         Y(J) = 1.0
C                             CALCULATE PHI(J,X), NATURAL CUBIC SPLINE
C                             WHICH IS EQUAL TO ONE AT KNOT X(J) AND
C                             ZERO AT OTHER KNOTS.  THEN SET
C                                 A(I,J) = PHI(J,XD(I)), I=1,...,M
         CALL DSPLN(X,Y,N,ZERO,ZERO,XD,A(1,J),M)
         Y(J) = 0.0
   10 CONTINUE
C                             CALL DLLSQR TO MINIMIZE NORM OF A*Y-YD
      CALL DLLSQR(A,M,N,Y,YD)
C                             LEAST SQUARES SPLINE IS
C                              Y(1)*PHI(1,X) + ... + Y(N)*PHI(N,X).
C                             EVALUATE SPLINE AT XOUT(1),...,XOUT(NOUT)
      CALL DSPLN(X,Y,N,ZERO,ZERO,XOUT,YOUT,NOUT)
      RETURN
      END
