C/ FIGURE 4.6.2
      SUBROUTINE DTRAN(WCAP,SREQ,COST,NW,NS,CMIN,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C                              DECLARATIONS FOR ARGUMENTS
      DOUBLE PRECISION WCAP(NW),SREQ(NS),COST(NW,NS),CMIN,X(NW,NS)
      INTEGER NW,NS
C                              DECLARATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION B(NW+NS),Y(NW+NS),C(NW*NS+NW+NS),
     &  XSOL(NW*NS+NW+NS)
C
C  SUBROUTINE DTRAN SOLVES THE TRANSPORTATION PROBLEM
C
C    MINIMIZE    CMIN = COST(1,1)*X(1,1) + ... + COST(NW,NS)*X(NW,NS)
C
C    WITH X(1,1),...,X(NW,NS) NONNEGATIVE, AND
C
C          X(1,1) + ... + X(1,NS)  .LE. WCAP(1)
C             .              .            .
C             .              .            .
C          X(NW,1)+ ... + X(NW,NS) .LE. WCAP(NW)
C          X(1,1) + ... + X(NW,1)    =  SREQ(1)
C             .              .            .
C             .              .            .
C          X(1,NS)+ ... + X(NW,NS)   =  SREQ(NS)
C
C  ARGUMENTS
C
C             ON INPUT                          ON OUTPUT
C             --------                          ---------
C
C    WCAP   - A VECTOR OF LENGTH NW CONTAINING
C             THE WAREHOUSE CAPACITIES.
C
C    SREQ   - A VECTOR OF LENGTH NS CONTAINING
C             THE STORE REQUIREMENTS.
C
C    COST   - THE NW BY NS COST MATRIX. COST(I,J)
C             IS THE PER UNIT COST TO SHIP FROM
C             WAREHOUSE I TO STORE J.
C
C    NW     - THE NUMBER OF WAREHOUSES.
C
C    NS     - THE NUMBER OF STORES.
C
C    CMIN   -                                   THE TOTAL COST OF THE
C                                               OPTIMAL ROUTING.
C
C    X      -                                   AN NW BY NS MATRIX
C                                               CONTAINING THE OPTIMAL
C                                               ROUTING.  X(I,J) UNITS
C                                               SHOULD BE SHIPPED FROM
C                                               WAREHOUSE I TO STORE J.
C
C-----------------------------------------------------------------------
      COMMON /CTRAN/ NWG,NSG
      EXTERNAL DTRAN2
      NWG = NW
      NSG = NS
      M = NW+NS
      N = NW*NS+NW+NS
C                              LOAD WAREHOUSE CAPACITIES AND STORE
C                              REQUIREMENTS INTO B VECTOR
      DO 5 I=1,NW
    5 B(I) = WCAP(I)
      DO 10 I=1,NS
   10 B(NW+I) = SREQ(I)
C                              FIRST NW*NS ENTRIES IN C ARE -COST(I,J).
C                              NEGATIVE SIGN PRESENT BECAUSE WE WANT
C                              TO MINIMIZE COST
      K = 0
      CNORM = 0.0
      DO 15 I=1,NW
      DO 15 J=1,NS
         K = K+1
         C(K) = -COST(I,J)
         CNORM = MAX(CNORM,ABS(C(K)))
   15 CONTINUE
C                              NEXT NW COSTS ARE ZERO, CORRESPONDING
C                              TO WAREHOUSE CAPACITY SLACK VARIABLES
      DO 20 I=1,NW
         K = K+1
         C(K) = 0.0
   20 CONTINUE
C                              LAST NS COSTS ARE LARGE AND NEGATIVE,
C                              CORRESPONDING TO "ARTIFICIAL WAREHOUSE"
C                              TRANSPORTATION COSTS
      ALPHA = 2*CNORM
      DO 25 I=1,NS
         K = K+1
         C(K) = -ALPHA
   25 CONTINUE
C                              USE REVISED SIMPLEX METHOD TO SOLVE
C                              TRANSPORTATION PROBLEM
      CALL DLPRV(DTRAN2,B,C,N,M,P,XSOL,Y)
C                              IF ANY ARTIFICIAL VARIABLES LEFT, THERE
C                              IS NO FEASIBLE SOLUTION
      DO 35 I=1,NS
         K = N-NS+I
         IF (XSOL(K).NE.0.0) THEN
            PRINT 30
   30       FORMAT (' ***** NO FEASIBLE SOLUTION *****')
            RETURN
         ENDIF
   35 CONTINUE
C                              FORM OPTIMAL ROUTING MATRIX, X
      CMIN = -P
      K = 0
      DO 40 I=1,NW
      DO 40 J=1,NS
         K = K+1
         X(I,J) = XSOL(K)
   40 CONTINUE
      RETURN
      END
 
      FUNCTION DTRAN2(Z,J)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION Z(*)
      COMMON /CTRAN/ NW,NS
C                              DTRAN2=DOT PRODUCT OF Z AND COLUMN J OF A
      JW = (J-1)/NS + 1
      JS = MOD(J-1,NS) + 1
      DTRAN2 = Z(JW) + Z(NW+JS)
      RETURN
      END
