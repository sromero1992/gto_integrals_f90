C/ FIGURE 6.2.4
      SUBROUTINE PCG(A,IROW,JCOL,NZ,X,B,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C                                  DECLARATIONS FOR ARGUMENTS
      DOUBLE PRECISION A(NZ),B(N),X(N)
      INTEGER IROW(NZ),JCOL(NZ)
C                                  DECLARATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION R(N),P(N),API(N),AP(N),LAMBDA
      include 'mpif.h'
C
C  SUBROUTINE PCG SOLVES THE SYMMETRIC LINEAR SYSTEM A*X=B, USING THE
C     CONJUGATE GRADIENT ITERATIVE METHOD.  THE NON-ZEROS OF A ARE STORED
C     IN SPARSE FORMAT.  THE COLUMNS OF A ARE DISTRIBUTED CYCLICALLY OVER
C     THE AVAILABLE PROCESSORS.
C
C  ARGUMENTS
C
C             ON INPUT                          ON OUTPUT
C             --------                          ---------
C
C    A      - A(IZ) IS THE MATRIX ELEMENT IN
C             ROW IROW(IZ), COLUMN JCOL(IZ),
C             FOR IZ=1,...,NZ. ELEMENTS WITH
C             MOD(JCOL(IZ)-1,NPES)=ITASK
C             ARE STORED ON PROCESSOR ITASK.
C
C    IROW   - (SEE A).
C
C    JCOL   - (SEE A).
C
C    NZ     - NUMBER OF NONZEROS STORED ON
C             THE LOCAL PROCESSOR.
C
C    X      -                                   AN N-VECTOR CONTAINING
C                                               THE SOLUTION.
C
C    B      - THE RIGHT HAND SIDE N-VECTOR.
C
C    N      - SIZE OF MATRIX A.
C
C-----------------------------------------------------------------------
C                                  NPES = NUMBER OF PROCESSORS
      CALL MPI_COMM_SIZE (MPI_COMM_WORLD,NPES,IERR)
C                                  ITASK = MY PROCESSOR NUMBER
      CALL MPI_COMM_RANK (MPI_COMM_WORLD,ITASK,IERR)
C                                  X0 = 0
C                                  R0 = B
C                                  P0 = R0
      X(1:N) = 0
      DO 10 I=ITASK+1,N,NPES
         R(I) = B(I)
         P(I) = R(I)
   10 CONTINUE
C                                  NITER = MAX NUMBER OF ITERATIONS
      NITER = 10000
      DO 90 ITER=1,NITER
C                                  AP = A*P
         API(1:N) = 0
         DO 20 IZ=1,NZ
            I = IROW(IZ)
            J = JCOL(IZ)
            API(I) = API(I) + A(IZ)*P(J)
   20    CONTINUE
C                                  MPI_ALLREDUCE COLLECTS THE VECTORS API
C                                  (API = LOCAL(A)*P) FROM ALL PROCESSORS
C                                  AND ADDS THEM TOGETHER, THEN SENDS
C                                  THE RESULT, AP, BACK TO ALL PROCESSORS.
         CALL MPI_ALLREDUCE(API,AP,N,MPI_DOUBLE_PRECISION,
     &    MPI_SUM,MPI_COMM_WORLD,IERR)
C                                  PAP = (P,AP)
         PAPI = 0.0
         DO 30 I=ITASK+1,N,NPES
            PAPI = PAPI + P(I)*AP(I)
   30    CONTINUE
         CALL MPI_ALLREDUCE(PAPI,PAP,1,MPI_DOUBLE_PRECISION,
     &    MPI_SUM,MPI_COMM_WORLD,IERR)
C                                  RP = (R,P)
         RPI = 0.0
         DO 40 I=ITASK+1,N,NPES
            RPI = RPI + R(I)*P(I)
   40    CONTINUE
         CALL MPI_ALLREDUCE(RPI,RP,1,MPI_DOUBLE_PRECISION,
     &    MPI_SUM,MPI_COMM_WORLD,IERR)
C                                  LAMBDA = (R,P)/(P,AP)
         LAMBDA = RP/PAP
C                                  X = X + LAMBDA*P
C                                  R = R - LAMBDA*AP
         DO 50 I=ITASK+1,N,NPES
            X(I) = X(I) + LAMBDA*P(I)
            R(I) = R(I) - LAMBDA*AP(I)
   50    CONTINUE
C                                  RAP = (R,AP)
         RAPI = 0.0
         DO 60 I=ITASK+1,N,NPES
            RAPI = RAPI + R(I)*AP(I)
   60    CONTINUE
         CALL MPI_ALLREDUCE(RAPI,RAP,1,MPI_DOUBLE_PRECISION,
     &    MPI_SUM,MPI_COMM_WORLD,IERR)
C                                  ALPHA = -(R,AP)/(P,AP)
         ALPHA = -RAP/PAP
C                                  P = R + ALPHA*P
         DO 70 I=ITASK+1,N,NPES
            P(I) = R(I) + ALPHA*P(I)
   70    CONTINUE
C                                  RMAX = MAX OF RESIDUAL (R)
         RMAXI = 0
         DO 80 I=ITASK+1,N,NPES
            RMAXI = MAX(RMAXI,ABS(R(I)))
   80    CONTINUE
         CALL MPI_ALLREDUCE(RMAXI,RMAX,1,MPI_DOUBLE_PRECISION,
     &    MPI_MAX,MPI_COMM_WORLD,IERR)
         IF (ITER.EQ.1) THEN
            R0MAX = RMAX
         ELSE
C                                  IF CONVERGED, MERGE PORTIONS OF X
C                                  STORED ON DIFFERENT PROCESSORS
            IF (RMAX.LE.1.D-10*R0MAX) THEN
               IF (ITASK.EQ.0) PRINT *, ' Number of iterations = ',ITER
               CALL MPI_ALLREDUCE(X,R,N,MPI_DOUBLE_PRECISION,
     &          MPI_SUM,MPI_COMM_WORLD,IERR)
               X(1:N) = R(1:N)
               RETURN
            ENDIF
         ENDIF
   90 CONTINUE
C                                  PCG DOES NOT CONVERGE
      IF (ITASK.EQ.0) PRINT 100
  100 FORMAT('***** PCG DOES NOT CONVERGE *****')
      RETURN
      END
