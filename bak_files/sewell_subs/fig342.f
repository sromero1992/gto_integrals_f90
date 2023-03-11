C/ FIGURE 3.4.2
      SUBROUTINE HESSH(A,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C                              DECLARATIONS FOR ARGUMENTS
      DOUBLE PRECISION A(N,N)
      INTEGER N
C                              DECLARATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION W(N)
      IF (N.LE.2) RETURN
C                              USE HOUSEHOLDER TRANSFORMATIONS TO
C                              REDUCE A TO UPPER HESSENBERG FORM
      DO 35 I=2,N-1
C                              CHOOSE UNIT N-VECTOR W (WHOSE FIRST
C                              I-1 COMPONENTS ARE ZERO) SUCH THAT WHEN
C                              COLUMN I-1 IS PREMULTIPLIED BY
C                              H = I - 2W*W**T, COMPONENTS I+1 THROUGH
C                              N ARE ZEROED.
         CALL CALW (A(1,I-1),N,W,I)
C                              PREMULTIPLY A BY H = I - 2W*W**T
         DO 15 K=I-1,N
            WTA = 0.0
            DO 5 J=I,N
               WTA = WTA + W(J)*A(J,K)
    5       CONTINUE
            TWOWTA = 2*WTA
            DO 10 J=I,N
               A(J,K) = A(J,K) - TWOWTA*W(J)
   10       CONTINUE
   15    CONTINUE
C                              POSTMULTIPLY A BY H = I - 2W*W**T
         DO 30 K=1,N
            WTA = 0.0
            DO 20 J=I,N
               WTA = WTA + W(J)*A(K,J)
   20       CONTINUE
            TWOWTA = 2*WTA
            DO 25 J=I,N
               A(K,J) = A(K,J) - TWOWTA*W(J)
   25       CONTINUE
   30    CONTINUE
   35 CONTINUE
      RETURN
      END
