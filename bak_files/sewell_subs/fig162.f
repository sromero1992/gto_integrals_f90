C/ FIGURE 1.6.2
      SUBROUTINE DSPLN(X,Y,N,YXX1,YXXN,XOUT,YOUT,NOUT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C                              DECLARATIONS FOR ARGUMENTS
      DOUBLE PRECISION X(N),Y(N),YXX1,YXXN,XOUT(NOUT),YOUT(NOUT)
      INTEGER N,NOUT
C                              DECLARATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION A(N-2,-1:2),COEFF(4,N),SIG(N),R(N)
C
C  SUBROUTINE DSPLN FITS AN INTERPOLATORY CUBIC SPLINE THROUGH THE
C    POINTS (X(I),Y(I)), I=1,...,N, WITH SPECIFIED SECOND DERIVATIVES
C    AT THE END POINTS, AND EVALUATES THIS SPLINE AT THE OUTPUT POINTS
C    XOUT(1),...,XOUT(NOUT).
C
C  ARGUMENTS
C
C             ON INPUT                          ON OUTPUT
C             --------                          ---------
C
C    X      - A VECTOR OF LENGTH N CONTAINING
C             THE X-COORDINATES OF THE DATA
C             POINTS.
C
C    Y      - A VECTOR OF LENGTH N CONTAINING
C             THE Y-COORDINATES OF THE DATA
C             POINTS.
C
C    N      - THE NUMBER OF DATA POINTS
C             (N.GE.3).
C
C    YXX1   - THE SECOND DERIVATIVE OF THE
C             CUBIC SPLINE AT X(1).
C
C    YXXN   - THE SECOND DERIVATIVE OF THE
C             CUBIC SPLINE AT X(N).  (YXX1=0
C             AND YXXN=0 GIVES A NATURAL
C             CUBIC SPLINE)
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
      SIG(1) = YXX1
      SIG(N) = YXXN
C                             SET UP TRIDIAGONAL SYSTEM SATISFIED
C                             BY SECOND DERIVATIVES (SIG(I)=SECOND
C                             DERIVATIVE AT X(I)).
      DO 5 I=1,N-2
         HI   = X(I+1)-X(I)
         HIP1 = X(I+2)-X(I+1)
         R(I) = (Y(I+2)-Y(I+1))/HIP1 - (Y(I+1)-Y(I))/HI
         A(I,-1) = HI/6.0
         A(I, 0) = (HI + HIP1)/3.0
         A(I, 1) = HIP1/6.0
         IF (I.EQ.  1) R(1)   = R(1)   - HI/  6.0*SIG(1)
         IF (I.EQ.N-2) R(N-2) = R(N-2) - HIP1/6.0*SIG(N)
    5 CONTINUE
C                             CALL DBAND TO SOLVE TRIDIAGONAL SYSTEM
      NLD = 1
      NUD = 1
      CALL DBAND(A,N-2,NLD,NUD,SIG(2),R)
C                             CALCULATE COEFFICIENTS OF CUBIC SPLINE
C                             IN EACH SUBINTERVAL
      DO 10 I=1,N-1
         HI = X(I+1)-X(I)
         COEFF(1,I) = Y(I)
         COEFF(2,I) = (Y(I+1)-Y(I))/HI - HI/6.0*(2*SIG(I)+SIG(I+1))
         COEFF(3,I) = SIG(I)/2.0
         COEFF(4,I) = (SIG(I+1)-SIG(I))/(6.0*HI)
   10 CONTINUE
      L = 1
      DO 25 I=1,NOUT
C                             FIND FIRST VALUE OF J FOR WHICH X(J+1) IS
C                             GREATER THAN OR EQUAL TO XOUT(I).  SINCE
C                             ELEMENTS OF XOUT ARE IN ASCENDING ORDER,
C                             WE ONLY NEED CHECK THE KNOTS X(L+1)...X(N)
C                             WHICH ARE GREATER THAN OR EQUAL TO
C                             XOUT(I-1).
         DO 15 J=L,N-1
            JSAVE = J
            IF (X(J+1).GE.XOUT(I)) GO TO 20
   15    CONTINUE
   20    L = JSAVE
C                             EVALUATE CUBIC SPLINE IN INTERVAL
C                             (X(L),X(L+1))
         P = XOUT(I)-X(L)
         YOUT(I) = COEFF(1,L)     + COEFF(2,L)*P
     &           + COEFF(3,L)*P*P + COEFF(4,L)*P*P*P
   25 CONTINUE
      RETURN
      END
