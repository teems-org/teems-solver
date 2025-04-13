
      INTEGER  LA, MAXN, MAXNE
      PARAMETER (LA=200, MAXN=20, MAXNE=100)
      REAL             A(LA),CNTL(10),RINFO(10),ERROR(3),RHS(MAXN),
     *         W(4*MAXN),X(MAXN)
      INTEGER  I,ICNTL(20),INFO(20),IRN(LA),IW(9*MAXN),JCN(LA),JOB,
     *         KEEP(7*MAXN),N,NE
      LOGICAL TRANS

C     Read in input matrix.
      READ(5,*) N,NE
      IF (N.GT.MAXN .OR. NE.GT.MAXNE) THEN
        WRITE(6,'(A,2I10)') ' Immediate STOP N or NE too large = ',
     *                   N, NE
        STOP
      END IF

      READ (5,FMT=*) (IRN(I),JCN(I),A(I),I=1,NE)

C     Set default controls
      CALL MA48I(CNTL,ICNTL)

C     Find pivot order, etc.
      JOB = 1
      CALL MA48A(N,N,NE,JOB,LA,A,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
      IF (INFO(1).LT.0) THEN
        WRITE (6,'(A,I3)') 'Error STOP from MA48A/AD with INFO(1) ='
     +    ,INFO(1)
        STOP
      END IF

C     Factorize matrix
      CALL MA48B(N,N,NE,JOB,LA,A,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,
     *            RINFO)
      IF (INFO(1).NE.0) THEN
        WRITE (6,FMT='(A,I3/A)') 'STOP from MA48B/BD with INFO(1) =',
     +    INFO(1),'Solution not possible'
        STOP
      END IF

C     Read in right-hand side.
      READ(5,*) (RHS(I),I=1,N)
C     Solve linear system without iterative refinement.
      TRANS = .FALSE.
      CALL MA48C (N,N,TRANS,JOB,LA,A,IRN,KEEP,CNTL,ICNTL,
     *            RHS,X,ERROR,W,IW,INFO)
C     Print out solution vector.
      WRITE(6,'(/A/(5D13.5))')'The solution vector is',(X(I),I=1,N)

C     Read in and solve system of similar equations.
      READ(5,* ) (A(I),I=1,NE)

C     Use fast factorize option but solve using iterative refinement.
      JOB = 2
      CALL MA48B(N,N,NE,JOB,LA,A,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,
     +            RINFO)

C     Read in right-hand side.
      READ (5,FMT=*) (RHS(I),I=1,N)

C Solve using iterative refinement and error estimation.
      JOB = 4
      CALL MA48C (N,N,TRANS,JOB,LA,A,IRN,KEEP,CNTL,ICNTL,
     *            RHS,X,ERROR,W,IW,INFO)
      WRITE(6,'(/A/(5D13.5))')'The solution vector is',(X(I),I=1,N)
      WRITE(6,'(/A/(3D13.5))')'The ERROR vector is',ERROR
      END

