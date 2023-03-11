C/ FIGURE 6.2.5
      PARAMETER (M=20,N=(M-1)**3,NZMAX=7*N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IROW(NZMAX),JCOL(NZMAX),AS(NZMAX),X(N),B(N)
      INCLUDE 'mpif.h'
C                     INITIALIZE MPI
      CALL MPI_INIT (IERR)
C                     NPES = NUMBER OF PROCESSORS
      CALL MPI_COMM_SIZE (MPI_COMM_WORLD,NPES,IERR)
C                     ITASK = MY PROCESSOR NUMBER
      CALL MPI_COMM_RANK (MPI_COMM_WORLD,ITASK,IERR)
C                     SOLVE 1.8.10 USING CONJUGATE GRADIENT ITERATION
      H = 1.d0/M
      L = 0
      NZ = 0
      DO 10 I=1,M-1
      DO 10 J=1,M-1
      DO 10 K=1,M-1
         L = L+1
         IF (MOD(L-1,NPES).EQ.ITASK) THEN
C                     DISTRIBUTE L-TH COLUMN TO PROCESSOR MOD(L-1,NPES)
            NZ = NZ + 1
            IROW(NZ) = L
            JCOL(NZ) = L
            AS(NZ) = 6
            IF (K.NE.1) THEN
               NZ = NZ + 1
               IROW(NZ) = L-1
               JCOL(NZ) = L
               AS(NZ) = -1
            ENDIF
            IF (K.NE.M-1) THEN
               NZ = NZ + 1
               IROW(NZ) = L+1
               JCOL(NZ) = L
               AS(NZ) = -1
            ENDIF
            IF (J.NE.1) THEN
               NZ = NZ + 1
               IROW(NZ) = L-(M-1)
               JCOL(NZ) = L
               AS(NZ) = -1
            ENDIF
            IF (J.NE.M-1) THEN
               NZ = NZ + 1
               IROW(NZ) = L+(M-1)
               JCOL(NZ) = L
               AS(NZ) = -1
            ENDIF
            IF (I.NE.1) THEN
               NZ = NZ + 1
               IROW(NZ) = L-(M-1)**2
               JCOL(NZ) = L
               AS(NZ) = -1
            ENDIF
            IF (I.NE.M-1) THEN
               NZ = NZ + 1
               IROW(NZ) = L+(M-1)**2
               JCOL(NZ) = L
               AS(NZ) = -1
            ENDIF
         ENDIF
         B(L) = H**2
   10 CONTINUE
      CALL PCG(AS,IROW,JCOL,NZ,X,B,N)
C                     SOLUTION AT BOX CENTER SHOULD BE ABOUT 0.056
      IF (ITASK.EQ.0) PRINT *, ' Solution at midpoint = ',X((N+1)/2)
      CALL MPI_FINALIZE(IERR)
      STOP
      END
