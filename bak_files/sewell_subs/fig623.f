C/ FIGURE 6.2.3
      SUBROUTINE PLPRG(A,B,C,N,M,P,X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C                              DECLARATIONS FOR ARGUMENTS
      DOUBLE PRECISION A(M,N),B(M),C(N),P,X(N),Y(M)
      INTEGER N,M
C                              DECLARATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION TAB(M+2,N+M+1),COLJP(M+2),LROWI(N+M),
     &  LROW(N+M)
      INTEGER BASIS(M)
      include 'mpif.h'
C
C  SUBROUTINE PLPRG USES THE SIMPLEX METHOD TO SOLVE THE PROBLEM
C
C             MAXIMIZE      P = C(1)*X(1) + ... + C(N)*X(N)
C
C    WITH X(1),...,X(N) NONNEGATIVE, AND
C
C           A(1,1)*X(1) + ... + A(1,N)*X(N)  =  B(1)
C             .                   .               .
C             .                   .               .
C           A(M,1)*X(1) + ... + A(M,N)*X(N)  =  B(M)
C
C    WHERE B(1),...,B(M) ARE ASSUMED TO BE NONNEGATIVE.
C
C  ARGUMENTS
C
C             ON INPUT                          ON OUTPUT
C             --------                          ---------
C
C    A      - THE M BY N CONSTRAINT COEFFICIENT
C             MATRIX.
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
C    N      - THE NUMBER OF UNKNOWNS.
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
C                              INITIALIZE MPI
      CALL MPI_INIT (IERR)
C                              NPES = NUMBER OF PROCESSORS
      CALL MPI_COMM_SIZE (MPI_COMM_WORLD,NPES,IERR)
C                              ITASK = MY PROCESSOR NUMBER (0,1,...,NPES-1).
C                              I WILL NEVER TOUCH ANY COLUMNS OF TAB EXCEPT
C                              MY COLUMNS, ITASK+1+ K*NPES, K=0,1,2,...
C                              (EXCEPT IN INITIALIZATION STAGE)
      CALL MPI_COMM_RANK (MPI_COMM_WORLD,ITASK,IERR)
C                              BASIS(1),...,BASIS(M) HOLD NUMBERS OF
C                              BASIS VARIABLES.  INITIAL BASIS CONSISTS
C                              OF ARTIFICIAL VARIABLES ONLY
      DO 5 I=1,M
         BASIS(I) = N+I
    5 CONTINUE
C                              INITIALIZE SIMPLEX TABLEAU
      DO 10 I=1,M+2
      DO 10 J=1,N+M+1
         TAB(I,J) = 0.0
   10 CONTINUE
C                              LOAD A INTO UPPER LEFT HAND CORNER
C                              OF TABLEAU
      DO 15 I=1,M
      DO 15 J=1,N
         TAB(I,J) = A(I,J)
   15 CONTINUE
C                              LOAD M BY M IDENTITY TO RIGHT OF A
C                              AND LOAD B INTO LAST COLUMN
      DO 20 I=1,M
         TAB(I,N+I) = 1.0
         TAB(I,N+M+1) = B(I)
   20 CONTINUE
C                              ROW M+1 CONTAINS -C, INITIALLY
      DO 25 J=1,N
         TAB(M+1,J) = -C(J)
   25 CONTINUE
C                              ROW M+2 CONTAINS COEFFICIENTS OF
C                              "ALPHA", WHICH IS TREATED AS +INFINITY
      DO 30 I=1,M
         TAB(M+2,N+I) = 1.0
   30 CONTINUE
C                              CLEAR "ALPHAS" IN LAST ROW
      DO 35 I=1,M
      DO 35 J=1,N+M+1
         TAB(M+2,J) = TAB(M+2,J) - TAB(I,J)
   35 CONTINUE
C                              SIMPLEX METHOD CONSISTS OF TWO PHASES
      DO 90 IPHASE=1,2
         IF (IPHASE.EQ.1) THEN
C                              PHASE I:  ROW M+2 (WITH COEFFICIENTS OF
C                              ALPHA) SEARCHED FOR MOST NEGATIVE ENTRY
            MROW = M+2
            LIM = N+M
         ELSE
C                              PHASE II:  FIRST N ELEMENTS OF ROW M+1
C                              SEARCHED FOR MOST NEGATIVE ENTRY
C                              (COEFFICIENTS OF ALPHA NONNEGATIVE NOW)
            MROW = M+1
            LIM = N
C                              IF ANY ARTIFICIAL VARIABLES LEFT IN
C                              BASIS AT BEGINNING OF PHASE II, THERE
C                              IS NO FEASIBLE SOLUTION
            DO 45 I=1,M
               IF (BASIS(I).GT.N) THEN
                  IF (ITASK.EQ.0) PRINT 40
   40             FORMAT (' ***** NO FEASIBLE SOLUTION *****')
                  RETURN
               ENDIF
   45       CONTINUE
         ENDIF
C                              THRESH = SMALL NUMBER.  WE ASSUME SCALES
C                              OF A AND C ARE NOT *TOO* DIFFERENT
         THRESHI = 0.0
         DO 50 J=ITASK+1,LIM,NPES
            THRESHI = MAX(THRESHI,ABS(TAB(MROW,J)))
   50    CONTINUE
         CALL MPI_ALLREDUCE(THRESHI,THRESH,1,MPI_DOUBLE_PRECISION,
     &   MPI_MAX,MPI_COMM_WORLD,IERR)
         THRESH = 1000*EPS*THRESH
C                              BEGINNING OF SIMPLEX STEP
   55    CONTINUE
C                              COLLECT PORTIONS (LROWI) OF LAST ROW FROM
C                              DIFFERENT PROCESSORS AND MERGE THEM
C                              INTO LROW, USING MPI_ALLREDUCE. 
            DO 56 J=1,LIM
               LROWI(J) = 0
               IF (MOD(J-1,NPES).EQ.ITASK) LROWI(J) = TAB(MROW,J)
   56       CONTINUE
            CALL MPI_ALLREDUCE(LROWI,LROW,LIM,MPI_DOUBLE_PRECISION,
     &      MPI_SUM,MPI_COMM_WORLD,IERR)
C                              FIND MOST NEGATIVE ENTRY IN ROW MROW,
C                              IDENTIFYING PIVOT COLUMN JP.
            CMIN = -THRESH
            JP = 0
            DO 60 J=1,LIM
               IF (LROW(J).LT.CMIN) THEN
                  CMIN = LROW(J)
                  JP = J
               ENDIF
   60       CONTINUE
C                              IF ALL ENTRIES NONNEGATIVE (ACTUALLY,
C                              IF GREATER THAN -THRESH) PHASE ENDS
            IF (JP.EQ.0) GO TO 90
C                              IF I OWN COLUMN JP, SAVE IT IN COLJP
            JTASK = MOD(JP-1,NPES)
            IF (ITASK.EQ.JTASK) THEN
               DO 61 I=1,MROW
                  COLJP(I) = TAB(I,JP)
   61          CONTINUE
            ENDIF
C                              RECEIVE COLJP FROM PROCESSOR THAT OWNS IT 
            CALL MPI_BCAST(COLJP,MROW,MPI_DOUBLE_PRECISION,
     &      JTASK,MPI_COMM_WORLD,IERR)
C                              FIND SMALLEST POSITIVE RATIO
C                              B(*)/TAB(*,JP), IDENTIFYING PIVOT
C                              ROW IP
            RATMIN = BIGNO
            IP = 0
            DO 65 I=1,M
               IF (COLJP(I).GT.THRESH) THEN
                  RATIO = TAB(I,N+M+1)/COLJP(I)
                  IF (RATIO.LT.RATMIN) THEN
                     RATMIN = RATIO
                     IP = I
                  ENDIF
               ENDIF
   65       CONTINUE
C                              IF ALL RATIOS NONPOSITIVE, MAXIMUM
C                              IS UNBOUNDED
            IF (IP.EQ.0) THEN
               IF (ITASK.EQ.0) PRINT 70
   70          FORMAT (' ***** UNBOUNDED MAXIMUM *****')
               RETURN
            ENDIF
C                              ADD X(JP) TO BASIS
            BASIS(IP) = JP
C                              NORMALIZE PIVOT ROW TO MAKE TAB(IP,JP)=1
            AMULT = 1.0/COLJP(IP)
            DO 75 J=ITASK+1,N+M,NPES
               TAB(IP,J) = AMULT*TAB(IP,J)
   75       CONTINUE
            TAB(IP,N+M+1) = AMULT*TAB(IP,N+M+1)
C                              ADD MULTIPLES OF PIVOT ROW TO OTHER
C                              ROWS, TO KNOCK OUT OTHER ELEMENTS IN
C                              PIVOT COLUMN
            DO 85 I=1,MROW
               IF (I.EQ.IP) GO TO 85
               AMULT = COLJP(I)
               DO 80 J=ITASK+1,N+M,NPES
                  TAB(I,J) = TAB(I,J) - AMULT*TAB(IP,J)
   80          CONTINUE
C                              MODIFY B(I) 
               TAB(I,N+M+1) = TAB(I,N+M+1) - AMULT*TAB(IP,N+M+1)
   85       CONTINUE
         GO TO 55
C                              END OF SIMPLEX STEP
   90 CONTINUE
C                              END OF PHASE II; READ X,P,Y FROM
C                              FINAL TABLEAU
      DO 95 J=1,N
         X(J) = 0.0
   95 CONTINUE
      DO 100 I=1,M
         K = BASIS(I)
         X(K) = TAB(I,N+M+1)
  100 CONTINUE
      P = TAB(M+1,N+M+1)
C                              COLLECT PORTIONS (LROWI) OF LAST ROW FROM
C                              DIFFERENT PROCESSORS AND MERGE THEM
C                              INTO LROW, USING MPI_ALLREDUCE. 
      DO 105 J=1,N+M
         LROWI(J) = 0.0
         IF (MOD(J-1,NPES).EQ.ITASK) LROWI(J) = TAB(M+1,J)
  105 CONTINUE
      CALL MPI_ALLREDUCE(LROWI,LROW,N+M,MPI_DOUBLE_PRECISION,
     & MPI_SUM,MPI_COMM_WORLD,IERR)
      DO 110 I=1,M
         Y(I) = LROW(N+I)
  110 CONTINUE
      CALL MPI_FINALIZE(IERR)
      RETURN
      END
