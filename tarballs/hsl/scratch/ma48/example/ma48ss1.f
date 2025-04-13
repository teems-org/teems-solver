C Specification sheet example.  Illustrating rectangular solution and
C     use of scaling routine MC29.
      INTEGER LA,MAXN,MAXNE
      PARAMETER (LA=200,MAXN=20,MAXNE=100)
      REAL             CNTL(10),RINFO(10),ERROR(3),A(LA),W(4*MAXN),
     *                 RHS(MAXN),X(MAXN),R(MAXN),C(MAXN)
      INTEGER JCN(LA),IRN(LA),IW(9*MAXN),KEEP(7*MAXN),ICNTL(20),
     *        INFO(20),JOB,M,N,NE,I,II,J,LP,IFAIL
      LOGICAL TRANS

C     Read in input matrix.
      READ (5,*) M,N,NE
      IF (M.GT.MAXN .OR. N.GT.MAXN .OR. NE.GT.MAXNE) THEN
        WRITE (6,'(A)') ' Error in input data for N and/or NE'
        STOP
      END IF
      READ (5,*) (IRN(I),JCN(I),A(I),I=1,NE)

C     Scale input matrix
      LP = 6
      CALL MC29A(M,N,NE,A,IRN,JCN,R,C,A(NE+1),LP,IFAIL)
      IF (IFAIL.LT.0) THEN
        WRITE (6,'(A)') ' Error in scaling routine'
        STOP
      END IF
      DO 10 I = 1,M
        R(I) = EXP(R(I))
   10 CONTINUE
      DO 20 I = 1,N
        C(I) = EXP(C(I))
   20 CONTINUE
      DO 30 II = 1,NE
        I = IRN(II)
        J = JCN(II)
        A(II) = A(II)*R(I)*C(J)
   30 CONTINUE

C     Set default controls
      CALL MA48I(CNTL,ICNTL)

C     Find pivot order, etc.
      JOB = 1
      CALL MA48A(M,N,NE,JOB,LA,A,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
      IF (INFO(1).LT.0) THEN
        WRITE (6,'(A,I3)') 'Error return from MA48AD with INFO(1) =',
     *                     INFO(1)
        STOP
      END IF

C     Factorize matrix
      CALL MA48B(M,N,NE,JOB,LA,A,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,
     *            RINFO)
      IF (INFO(1).NE.0) THEN
        WRITE (6,'(A,I3/A)') 'Return from MA48B with INFO(1) =',
     *                       INFO(1),'Solution not possible'
        STOP
      END IF
C     Read in right hand side.
      READ (5,*) (RHS(I),I=1,M)
C     Scale the right hand side vector by row weights

      DO 40 I = 1,M
        RHS(I) = RHS(I)*R(I)
   40 CONTINUE

C     Solve linear system.
      TRANS = .FALSE.
      CALL MA48C(M,N,TRANS,JOB,LA,A,IRN,KEEP,CNTL,ICNTL,RHS,X,ERROR,W,
     *            IW,INFO)

C     Scale the right-hand side vector by column weights

      DO 50 I = 1,N
        X(I) = X(I)*C(I)
   50 CONTINUE

C     Print out solution vector.
      WRITE (6,'(/A/(5E13.5))') ' The solution vector is:',
     *                          (X(I),I=1,N)
      END
