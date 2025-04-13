! COPYRIGHT (c) 1993 Council for the Central Laboratory
!                    of the Research Councils
!
! Version 2.2.0
! See ChangeLog for version history.
!
      SUBROUTINE MA48AD(M,N,NE,JOB,LA,A,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,
     +                  RINFO)
C Given a sparse matrix, find its block upper triangular form, choose a
C     pivot sequence for each diagonal block, and prepare data
C     structures for actual factorization.
C
C     .. Arguments ..
      INTEGER M,N,NE,JOB,LA
      DOUBLE PRECISION A(LA)
      INTEGER IRN(LA),JCN(LA),KEEP(*)
      DOUBLE PRECISION CNTL(10)
      DOUBLE PRECISION RINFO(10)
      INTEGER ICNTL(20),IW(6*M+3*N),INFO(20)

C
C M must be set by the user to the number of rows.
C      It is not altered by the subroutine.  Restriction:  M > 0.
C N must be set by the user to the number of columns.
C      It is not altered by the subroutine.  Restriction:  N > 0.
C NE must be set by the user to the number of entries in the input
C      matrix. It is not altered by the subroutine. Restriction: NE > 0.
C JOB must be set by the user to 1 for automatic choice of pivots and
C      to 2 if the pivot sequence is specified in KEEP.  If JOB is set
C      to 3, then pivots are chosen automatically from the diagonal as
C      long as this is numerically feasible.
C      It is not altered by the subroutine.  Restriction: 1 <= JOB <= 3.
C LA must be set by the user to the size of A, IRN, and JCN.
C      It is not altered by the subroutine. Restriction LA >= 2*NE.
C      Normally a value of 3*NE will suffice.
C A must have its first NE elements set by the user to hold the matrix
C      entries. They may be in any order. If there is more than one for
C      a particular matrix position, they are accumulated. The first
C      NE entries are not altered by the subroutine. The rest is used
C      as workspace for the active submatrix.
C IRN  is an integer array. Entries 1 to NE must be set to the
C      row indices of the corresponding entries in A.  On return, the
C      leading part holds the row numbers of the permuted matrix, with
C      duplicates excluded. The permuted matrix is block upper
C      triangular with recommended pivots on its diagonal. The entries
C      of the block diagonal part are held contiguously from IRN(1)
C      and the entries of the off-diagonal part are held contiguously
C      backwards from IRN(NE) during the computation.  At the end the
C      off-diagonal blocks are held contiguously immediately after the
C      diagonal blocks.
C JCN  is an integer array. Entries 1 to NE must be set to the column
C      indices of the corresponding entries in A. On return, JCN(k)
C      holds the position in IRN that corresponds to the entry that was
C      input in A(k), k=1,NE.  Entries corresponding to out-of-range
C      indices are first set to 0 then, before exit, to NE+1.
C      If duplicates are found, the first entry of JCN is negated.
C KEEP is an integer array that need not be set on a JOB=1 or JOB=3
C      entry. On entry with JOB=2 and always
C      on a successful exit, KEEP(i) holds the position of row i in the
C      permuted matrix, I=1,M and KEEP(M+j) holds the index of the
C      column that is in position j of the permuted matrix, j=1,N.
C      The rest of the array need not be set on entry. On exit:
C        KEEP(IPTRD+j), IPTRD=M+3*N, holds the position in IRN of the
C          start of the block diagonal part of column j, J=1,N;
C        KEEP(IPTRD+N+1) holds the position that immediately follows
C          the end of the block diagonal part of column N;
C        KEEP(IPTRO+j),IPTRO=IPTRD+N+1, holds the position in IRN of
C          the start of the block off-diagonal part of column j, j=1,N;
C        KEEP(IPTRO+N+1) holds the position that immediately
C          follows the end of the block off-diagonal part of column N;
C          During the computation, the columns of the off-diagonal
C          blocks are held in reverse order and KEEP(IPTRO+N+1) points
C          to the position immediately before the off-diagonal block
C          storage.
C        KEEP(KBLOCK+3), KBLOCK=IPTRO+N+1, holds the number of
C          blocks NB in the block triangular form;
C        KEEP(NBLOCK+3*k), NBLOCK=IPTRO+N-1, holds the number of columns
C          in block k, k=1,NB; and
C        KEEP(MBLOCK+3*k), MBLOCK=IPTRO+N, is negative if block k
C          is triangular or holds the number of rows held in packed
C          storage when processing block k, k=1,NB.
C        KEEP(LBLOCK+k), LBLOCK=KBLOCK+3*NB is accessed only if ICNTL(8)
C          is not equal to 0 and will be set to the number of columns
C          in block k which are to be pivoted on last, k=1,NB.
C CNTL  is a real array of length 10 that must be set by the user
C       as follows and is not altered.
C     CNTL(1)  If this is set to a value less than or equal to one, full
C       matrix processing will be used by MA50A/AD when the density of
C       the reduced matrix reaches CNTL(1).
C     CNTL(2) determines the balance used by MA50A/AD and MA50B/BD
C       between pivoting for sparsity and for stability, values near
C       zero emphasizing sparsity and values near one emphasizing
C       stability.
C     CNTL(3) If this is set to a positive value, any entry whose
C       modulus is less than CNTL(3) will be dropped from the factors
C       calculated by MA50A/AD. The factorization will then require
C       less storage but will be inaccurate.
C     CNTL(4)  If this is set to a positive value, any entry whose
C       modulus is less than CNTL(4) will be regarded as zero from
C       the point of view of rank.
C     CNTL(5:10) are not accessed by this subroutine.
C ICNTL is an integer array of length 20 that must be set by the user
C       as follows and is not altered.
C     ICNTL(1)  must be set to the stream number for error messages.
C       A value less than 1 suppresses output.
C     ICNTL(2) must be set to the stream number for diagnostic output.
C       A value less than 1 suppresses output.
C     ICNTL(3) must be set to control the amount of output:
C       0 None.
C       1 Error messages only.
C       2 Error and warning messages.
C       3 As 2, plus scalar parameters and a few entries of array
C         parameters on entry and exit.
C       4 As 2, plus all parameters on entry and exit.
C     ICNTL(4)  If set to a positive value, the pivot search by MA50A/AD
C       is limited to ICNTL(4) columns. This may result in different
C       fill-in and execution time but could give faster execution.
C     ICNTL(5) The block size to be used for full-matrix processing.
C     ICNTL(6) is the minimum size for a block of the block triangular
C       form.
C     ICNTL(7) If not equal to 0, abort when structurally rank deficient
C       matrix found.
C     ICNTL(8) If set to a value other than zero and JOB = 1 or 3,
C       columns with IW flagged 0 are placed at
C       the end of their respective blocks in the first factorization
C       and the remaining columns are assumed unchanged in subsequent
C       factorizations. If set to a value other than zero and JOB = 2,
C       columns before the first 0 entry of IW are assumed unchanged in
C       subsequent factorizations.
C     ICNTL(9:20) are not accessed by this subroutine.
C IW   is an integer array of length 6*M+3*N that is used as workspace.
C      If ICNTL(8) is not 0, the first N entries of IW must be set on
C      entry so that IW(i), i=1,N is zero for and only for the columns
C      designated to be at the end of the pivot sequence. Since IW is
C      used as workspace, this will be destroyed by the subroutine.
C INFO is an integer array of length 20 that need not be set on entry.
C   INFO(1)  On exit, a negative value  will indicate an error return
C      and a positive value will indicate a warning.
C      Possible nonzero values are:
C      -1  M < 1 or N < 1
C      -2  NE < 1
C      -3  Insufficient space
C      -4  ICNTL(7) not equal to 0 and the matrix is structurally rank
C          deficient.
C      -5  Faulty permutation input on JOB = 2 entry
C      -6  JOB out of range
C      +1  Row or column number out of range (such entries are ignored)
C          and/or more than one entry for the same matrix position.
C      +2  Matrix rank deficient.  Estimated rank in INFO(5).
C          more than one entry for the same matrix position
C      +3  Combination of warnings +1 and +2.
C      +4  Premature switch to full code because of problems in choosing
C          diagonal pivots (JOB = 3 entry).
C      +5  Combination of warnings +1 and +4.
C      +6  Combination of warnings +2 and +4.
C      +7  Combination of warnings +1, +2 and +4.
C    INFO(2) Number of compresses of the files.
C    INFO(3) Minimum storage required to analyse matrix.
C    INFO(4) Minimum storage required to factorize matrix.
C    INFO(5) Upper bound on the rank of the matrix.
C    INFO(6) Number of entries dropped from the data structure.
C    INFO(7) Order of the largest nontriangular block on the diagonal
C            of the block triangular form.
C    INFO(8) Total of the orders of all the nontriangular blocks on the
C            diagonal of the block triangular form.
C    INFO(9) Total no. of entries in all the nontriangular blocks on
C            the diagonal of the block triangular form.
C    INFO(10) Structural rank.
C    INFO(11) Number of multiple entries in the input data.
C    INFO(12) Number of entries with out-of-range indices.
C RINFO is a real array that need not be set on entry. On exit,
C    RINFO(1) holds the number of floating-point operations needed for
C      the factorization.

C     .. Local constants ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.D0)
C
C     .. Local variables ..
      DOUBLE PRECISION CNTL5(10)
      INTEGER EYE,HEADC,I,IB,ICNTL5(20),IDUMMY,INFO5(15),IP,IPTRA,IPTRD,
     +        IPTRO,IPTRP,IQ,ISW,IW13,IW21,IW50,J,JAY,JB,JFIRST,J1,J2,
     +        J3,K,KB,KBLOCK,KD,KK,KL,KO,L,LASTR,LASTC,LBLOCK
      LOGICAL LDUP
      INTEGER LENC,LENP,LENR,LP,MBLOCK,MINBLK,MP,NB,NBLOCK,NC,NDIAG,
     +        NEXTC,NEXTR,NP,NR,NXTEYE,NZA,NZB,NZD,PTRD,PTRO
      DOUBLE PRECISION RINFO5(10),TOL
C CNTL5 passed to MA50A/AD to correspond to dummy argument CNTL.
C EYE   position in arrays.
C HEADC displacement in IW. IW(HEADC+J) holds the position in A of the
C       head of the chain of entries for column J.
C I     row index.
C IB    displacement in IW. IW(IB+K) first holds the index in the
C       permuted matrix of the first column of the k-th block. Later
C       abs(IW(IB+K)) holds the number of columns in block K and is
C       negative for a triangular block.
C ICNTL5 passed to MA50A/AD to correspond to dummy argument ICNTL.
C IDUMMY do loop index not used within the loop.
C INFO5 passed to MA50A/AD to correspond to dummy argument INFO.
C IP    displacement in IW. IW(IP+I), I=1,M holds a row permutation.
C IPTRA displacement in IW. IW(IPTRA+J) holds the position in the
C       reordered matrix of the first entry of column J.
C IPTRD displacement in KEEP. See comment on KEEP.
C IPTRO displacement in KEEP. See comment on KEEP.
C IPTRP displacement in IW. IW(IP+I), I=1,N holds the column starts
C       of a permutation of the matrix, during the block triangular
C       calculation.
C IQ    displacement in IW. IW(IQ+J), J=1,N holds a column permutation.
C ISW   used for swopping two integers
C IW13  displacement in IW. IW(IW13+1) is the first entry of the
C       workarray IW of MC13DD.
C IW21  displacement in KEEP. KEEP(IW21+1) is the first entry of the
C       workarray IW of MC21AD.
C IW50  displacement in IW. IW(IW50+1) is the first entry of the
C       workarray IW of MA50A/AD.
C J     column index.
C JAY   column index.
C JB    block index.
C JFIRST displacement in IW. IW(JFIRST+1) is the first entry of the
C       workarray JFIRST of MA50A/AD.
C J1    first column index of a block.
C J2    last column index of a block.
C J3    last column index of a block ... used in printing
C K     running index for position in matrix.
C KB    block index.
C KBLOCK displacement in KEEP. See comment on KEEP.
C KD    running index for position in matrix.
C KK    running index for position in matrix.
C KL    last position in matrix of current column.
C KO    running index for position in matrix.
C L     length of a block.
C LBLOCK displacement in KEEP. See comment on KEEP.
C LASTR displacement in IW. IW(LASTR+I) holds the position of the last
C       entry encountered in row I.  Used for accumulating duplicates.
C LASTC displacement in IW. IW(LASTC+I) holds the index of the
C       last column encountered that had an entry in row I.
C       LASTC is also used to get workspace in KEEP for MA50AD.
C LDUP  logical flag used to indicate whether duplicates have been
C       found and recorded.
C LENC  displacement in IW. IW(LENC+J) holds the number of entries
C       in column J, excluding duplicates.
C       LENC is also used to get workspace in KEEP for MA50AD.
C LENP  displacement in IW. IW(LENP+J) holds the number of entries
C       in column J of the permuted matrix.
C LENR  displacement in IW. IW(LENR+1) is the first entry of the
C       workarray LENR of MA50A/AD.
C LP Unit for error messages.
C MBLOCK displacement in KEEP. See comment on KEEP.
C MINBLK Minimum size for a block.
C MP Unit for diagnostic messages.
C NB Number of blocks.
C NBLOCK displacement in KEEP. See comment on KEEP.
C NC Number of columns in block.
C NDIAG number entries placed on diagonal by MC21AD.
C NEXTC displacement in IW. IW(NEXTC+1) is the first entry of the
C       workarray NEXTC of MA50A/AD.
C NEXTR displacement in IW. IW(NEXTR+1) is the first entry of the
C       workarray NEXTR of MA50A/AD.
C NP Number of columns in packed storage.
C NR Number of rows in block.
C NXTEYE Next value of EYE.
C NZA   number of entries in the matrix, excluding duplicates.
C NZB   number of entries in the current diagonal block.
C NZD   total number of entries in the diagonal blocks.
C PTRD  displacement in IW. IW(PTRD+J) holds the position of the
C       first entry of the diagonal part of column J.
C PTRO  displacement in IW. IW(PTRO+J) holds the position of the
C       first entry of the off-diagonal part of column J.
C RINFO5 passed to MA50A/AD to correspond to dummy argument RINFO.
C TOL   pivot tolerance.  Entries less than this value are considered
C       as zero.
C
C .. Externals ..
      EXTERNAL MC13DD,MC21AD,MA50AD
C MC13DD    finds permutation for block triangular form
C MC21AD    finds permutation for zero-free diagonal
C MA50A/AD factorizes non-triangular diagonal blocks
      INTRINSIC ABS,MAX
C
C
C Simple data checks
C Set printer streams
      LP = ICNTL(1)
      MP = ICNTL(2)
      IF (M.LE.0 .OR. N.LE.0) THEN
         INFO(1) = -1
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9100) M,N
         GO TO 530
      END IF
      IF (NE.LE.0) THEN
         INFO(1) = -2
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9110) NE
         GO TO 530
      END IF
      IF (LA.LT.2*NE) THEN
         INFO(1) = -3
         INFO(3) = 2*NE
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9120) LA,INFO(3)
         GO TO 530
      END IF
      IF (JOB.LT.1 .OR. JOB.GT.3) THEN
         INFO(1) = -6
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9130) JOB
         GO TO 530
      END IF
C Check input permutations
      IF (JOB.EQ.2) THEN
         DO 10 I = 1,MAX(M,N)
            IW(N+I) = 0
   10    CONTINUE
C Check row permutation
         DO 20 I = 1,M
            J = KEEP(I)
            IF (J.LT.1 .OR. J.GT.M) GO TO 40
            IF (IW(N+J).EQ.1) GO TO 40
            IW(N+J) = 1
   20    CONTINUE
C Check column permutation
         DO 30 I = 1,N
            J = KEEP(M+I)
            IF (J.LT.1 .OR. J.GT.N) GO TO 40
            IF (IW(N+J).EQ.2) GO TO 40
            IW(N+J) = 2
   30    CONTINUE
         GO TO 50
   40    INFO(1) = -5
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9140)
         GO TO 530
      END IF

   50 IF (MP.GT.0 .AND. ICNTL(3).GT.2) THEN
         WRITE (MP,
     +     '(/A/A,I7,A,I6,A,I7,A,I2,A,I7/A,1P,4D12.4/A,4I8/A,3I8)')
     +     ' Entering MA48A/AD with',' M =',M,'     N =',N,'     NE =',
     +     NE,'     JOB =',JOB,'     LA =',LA,' CNTL (1:4) =',
     +     (CNTL(I),I=1,4),' ICNTL(1:4) = ', (ICNTL(I),I=1,4),
     +     ' ICNTL(6:8) = ', (ICNTL(I),I=6,8)
         IF (ICNTL(3).GT.3) THEN
            WRITE (MP,9000) (A(K),IRN(K),JCN(K),K=1,NE)
 9000       FORMAT (' Entries:'/3 (1P,D12.4,2I6))
         ELSE
            WRITE (MP,9000) (A(K),IRN(K),JCN(K),K=1,MIN(9,NE))
         END IF
         IF (JOB.EQ.2) THEN
            WRITE (MP,'(A)') ' Permutations input (JOB=2)'
            IF (ICNTL(3).GT.3) THEN
               WRITE (MP,9010) (KEEP(I),I=1,M)
 9010          FORMAT (' Positions of original rows in the permuted ma',
     +                'trix'/ (10I6))
               WRITE (MP,9020) (KEEP(M+I),I=1,N)
 9020          FORMAT (' Positions of columns of permuted matrix ','in',
     +                ' or','iginal matrix '/ (10I6))
            ELSE
               WRITE (MP,9010) (KEEP(I),I=1,MIN(10,M))
               WRITE (MP,9020) (KEEP(M+I),I=1,MIN(10,N))
            END IF
         END IF
         IF (ICNTL(8).NE.0) THEN
            WRITE (MP,'(A,I6)')
     +        ' Value of IW entries on call with ICNTL(8) =',ICNTL(8)
            IF (ICNTL(3).GT.3) THEN
               WRITE (MP,9030) (IW(I),I=1,N)
 9030          FORMAT (10I6)
            ELSE
               WRITE (MP,9030) (IW(I),I=1,MIN(10,N))
            END IF
         END IF
      END IF
C
C Initializations
      DO 53 I = 1,20
         INFO(I) = 0
   53 CONTINUE
      INFO(3) = NE*2
      INFO(4) = NE
      INFO(10) = MIN(M,N)
      DO 56 I = 1,5
        RINFO(I) = ZERO
   56 CONTINUE
      DO 60 I = 1,4
         CNTL5(I) = CNTL(I)
         ICNTL5(I) = ICNTL(I)
   60 CONTINUE
C Switch off printing from MA50A/AD
      ICNTL5(3) = 0
C Set BLAS control
      ICNTL5(5) = ICNTL(5)
C Default control for restricted pivoting
      ICNTL5(6) = 0
C Set controls if pivot sequence provided
      ICNTL5(7) = 0
      IF (JOB.EQ.2) THEN
         ICNTL5(7) = 2
         ICNTL5(4) = 1
      END IF
C Set control for diagonal pivoting
      IF (JOB.EQ.3) ICNTL5(7) = 1

C Set pivot tolerance
      TOL = MAX(ZERO,CNTL(4))
C Set blcok size for BTF
      MINBLK = MAX(1,ICNTL(6))
C Set flag for scan for duplicates
      LDUP = .FALSE.

C Partition KEEP
      IPTRD = M + 3*N
      IPTRO = IPTRD + N + 1
      NBLOCK = IPTRO + N - 1
      MBLOCK = NBLOCK + 1
      KBLOCK = MBLOCK + 1

C Partition IW
C We assume at this point that there might be data in IW(1:N)
      HEADC = N + 1
      LASTC = HEADC + N
C
C Initialize header array for column links.
      DO 70 J = 1,N
         IW(HEADC+J) = 0
   70 CONTINUE
C Overwrite JCN by column links, checking that rows and column numbers
C     are within range.
      DO 80 K = 1,NE
         I = IRN(K)
         J = JCN(K)
         IF (I.LT.1 .OR. I.GT.M .OR. J.LT.1 .OR. J.GT.N) THEN
            INFO(12) = INFO(12) + 1
            IF (MP.GT.0 .AND. INFO(12).LE.10 .AND.
     +          ICNTL(3).GE.2) WRITE (MP,'(A,I7,A,2I6)')
     +          ' Message from MA48A/AD .. indices for entry ',K,' are',
     +          I,J
            JCN(K) = 0
         ELSE
            JCN(K) = IW(HEADC+J)
            IW(HEADC+J) = K
         END IF
   80 CONTINUE

C
      IF (MINBLK.GE.N .OR. M.NE.N .OR. JOB.GT.1) GO TO 190
C
C Obtain permutations to put matrix in block triangular form.
C
C Note that in this part of the code we know that M is equal to N.

C Get workspace from KEEP for MC21AD and MC13DD
      IW21 = 2*N
      IW13 = IW21
C Partition IW for arrays used by MC21AD and MC13DD
      IPTRA = LASTC + N
      LENC = IPTRA + N
      IB = LENC
      IP = LENC + N
      IPTRP = IP + N
      LENP = IPTRP + N

C Initialize array LASTC
      DO 90 I = 1,N
         IW(LASTC+I) = 0
   90 CONTINUE

C Run through the columns removing duplicates and generating copy of
C     structure by columns.
      LDUP = .TRUE.
      K = 1
      DO 120 J = 1,N
         EYE = IW(HEADC+J)
         IW(IPTRA+J) = K
         DO 100 IDUMMY = 1,NE
            IF (EYE.EQ.0) GO TO 110
            I = IRN(EYE)
C Check for duplicates
            IF (IW(LASTC+I).NE.J) THEN
               IW(LASTC+I) = J
               IRN(NE+K) = I
               K = K + 1
            ELSE
C Record duplicate
               INFO(11) = INFO(11) + 1
               IF (MP.GT.0 .AND. INFO(11).LE.10 .AND.
     +             ICNTL(3).GE.2) WRITE (MP,'(A,I7,A,2I6)')
     +             ' Message from MA48A/AD .. duplicate in position ',K,
     +             ' with indices',I,J
            END IF
            EYE = JCN(EYE)
  100    CONTINUE
  110    IW(LENC+J) = K - IW(IPTRA+J)
  120 CONTINUE

C NZA is number of entries with duplicates and out-of-range indices
C     removed.
      NZA = K - 1

C   Compute permutation for zero-free diagonal
      CALL MC21AD(N,IRN(NE+1),NZA,IW(IPTRA+1),IW(LENC+1),IW(IP+1),NDIAG,
     +           KEEP(IW21+1))
C IW(IP+1) ... col permutation ... new.to.old
      INFO(10) = NDIAG
      IF (NDIAG.LT.N) THEN
         IF (ICNTL(7).NE.0) THEN
            INFO(1) = -4
            IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,
     +          '(A,A/A,I7,A,I7)')
     +        ' Error return from MA48A/AD because matrix structurally '
     +          ,' singular',' order is ',N,' and structural rank',NDIAG
            GO TO 530
         END IF
         GO TO 190
      END IF
C Set pointers and column counts for permuted matrix
CDIR$ IVDEP
      DO 130 J = 1,N
         JAY = IW(IP+J)
         IW(IPTRP+J) = IW(IPTRA+JAY)
         IW(LENP+J) = IW(LENC+JAY)
  130 CONTINUE

C Use Tarjan's algorithm to obtain permutation to block triangular form.
C KEEP(M+1) .. col permutation ...new.to.old
      CALL MC13DD(N,IRN(NE+1),NZA,IW(IPTRP+1),IW(LENP+1),KEEP(M+1),
     +           IW(IB+1),NB,KEEP(IW13+1))
C IW(IB+1) ... col number in permuted matrix of start of blocks

C Merge adjacent 1x1 blocks into triangular blocks
C Change block pointers to block sizes
      DO 140 JB = 2,NB
         IW(IB+JB-1) = IW(IB+JB) - IW(IB+JB-1)
  140 CONTINUE
      IW(IB+NB) = N + 1 - IW(IB+NB)
      IF (IW(IB+1).EQ.1) IW(IB+1) = -1
      KB = 1
      DO 150 JB = 2,NB
         L = IW(IB+JB)
         IF (L.EQ.1 .AND. IW(IB+KB).LE.0) THEN
            IW(IB+KB) = IW(IB+KB) - 1
         ELSE
            KB = KB + 1
            IF (L.EQ.1) THEN
               IW(IB+KB) = -1
            ELSE
               IW(IB+KB) = L
            END IF
         END IF
  150 CONTINUE
      NB = KB
C Merge small blocks
      KB = 1
      DO 160 JB = 2,NB
         IF (ABS(IW(IB+KB)).LT.MINBLK) THEN
            IW(IB+KB) = ABS(IW(IB+KB)) + ABS(IW(IB+JB))
         ELSE
            KB = KB + 1
            IW(IB+KB) = IW(IB+JB)
         END IF
  160 CONTINUE
      NB = KB
C Set size of blocks and triangular block flags
      DO 170 JB = 1,NB
         KEEP(NBLOCK+3*JB) = ABS(IW(IB+JB))
         KEEP(MBLOCK+3*JB) = IW(IB+JB)
  170 CONTINUE
C Record permutations.
C Set KEEP to position of column that is in position j of permuted
C     matrix .... new.to.old
      DO 180 J = 1,N
         KEEP(KEEP(M+J)) = J
         KEEP(M+J) = IW(IP+KEEP(M+J))
  180 CONTINUE
      GO TO 220
C
C Set arrays for the case where block triangular form is not computed
  190 NB = 1
      IF (JOB.EQ.1 .OR. JOB.EQ.3) THEN
         DO 200 I = 1,M
            KEEP(I) = I
  200    CONTINUE
         DO 210 I = 1,N
            KEEP(M+I) = I
  210    CONTINUE
      END IF
C Set number of columns in block
      KEEP(NBLOCK+3) = N
      KEEP(MBLOCK+3) = 0

C Control for forcing columns to end.
  220 IF (ICNTL(8).NE.0) THEN
         LBLOCK = KBLOCK + 3*NB
         IF (JOB.EQ.2) THEN
C Find first column that will be changed in subsequent factorizations.
            DO 230 I = 1,N
               IF (IW(I).EQ.0) GO TO 240
  230       CONTINUE
C Set LBLOCK to number of columns from and including first one flagged
  240       KEEP(LBLOCK+1) = N - I + 1
         ELSE
C Within each block move columns to the end
            J = 1
            DO 270 JB = 1,NB
               KEEP(LBLOCK+JB) = 0
               J2 = J + KEEP(NBLOCK+3*JB) - 1
               J1 = J2
C Jump if triangular block (leave columns in situ)
               IF (KEEP(MBLOCK+3*JB).LT.0) GO TO 260
C Run through columns in block sorting those with IW flag to end
  250          IF (J.EQ.J2) GO TO 260
               IF (IW(KEEP(M+J)).EQ.0) THEN
C Column is put at end
                  KEEP(LBLOCK+JB) = KEEP(LBLOCK+JB) + 1
                  ISW = KEEP(M+J2)
                  KEEP(M+J2) = KEEP(M+J)
                  KEEP(M+J) = ISW
                  J2 = J2 - 1
               ELSE
                  J = J + 1
               END IF
               GO TO 250
  260          J = J1 + 1
  270       CONTINUE
         END IF
      END IF

C Run through the columns to create block-ordered form, removing
C     duplicates, changing row indices, and holding map array in JCN.
C Grab space in IW for LASTR
      LASTR = LASTC + M
C Initialize LASTC
      DO 280 I = 1,M
         IW(LASTC+I) = 0
  280 CONTINUE
      KEEP(KBLOCK+3) = NB
      K = 1
      KK = NE
      J2 = 0
      DO 310 JB = 1,NB
         J1 = J2 + 1
         J2 = J1 + KEEP(NBLOCK+3*JB) - 1
C Run through columns in block JB
C LASTR is used for accumulation of duplicates
C LASTC is used to identify duplicates
         DO 300 JAY = J1,J2
            J = KEEP(M+JAY)
            EYE = IW(HEADC+J)
            KEEP(IPTRD+JAY) = K
            KEEP(IPTRO+JAY) = KK
            IF (KEEP(MBLOCK+3*JB).LT.0) THEN
C Block is triangular.
C Reserve the leading position for the diagonal entry
               IW(LASTC+JAY) = JAY
               A(NE+K) = ZERO
               IRN(NE+K) = JAY
               IW(LASTR+JAY) = K
               K = K + 1
            END IF
            DO 290 IDUMMY = 1,NE
               IF (EYE.EQ.0) GO TO 300
               NXTEYE = JCN(EYE)
               I = KEEP(IRN(EYE))
               IF (IW(LASTC+I).NE.JAY) THEN
C Entry encountered for the first time.
                  IW(LASTC+I) = JAY
                  IF ((I.GE.J1.AND.I.LE.J2) .OR. (M.NE.N)) THEN
C Entry in diagonal block.
                     A(NE+K) = A(EYE)
                     IRN(NE+K) = I
                     IW(LASTR+I) = K
                     JCN(EYE) = K
                     K = K + 1
                  ELSE
C Entry in off-diagonal block.
                     A(NE+KK) = A(EYE)
                     IRN(NE+KK) = I
                     IW(LASTR+I) = KK
                     JCN(EYE) = KK
                     KK = KK - 1
                  END IF
               ELSE
C Entry has already been encountered.
                  IF (.NOT.LDUP) THEN
C Duplicates should be counted here if they were not earlier
                     INFO(11) = INFO(11) + 1
                     IF (MP.GT.0 .AND. INFO(11).LE.10 .AND.
     +                   ICNTL(3).GE.2) WRITE (MP,'(A,I7,A,2I6)')
     +                ' Message from MA48A/AD .. duplicate in position '
     +                   ,EYE,' with indices',IRN(EYE),J
                  END IF
                  KL = IW(LASTR+I)
                  JCN(EYE) = KL
                  A(NE+KL) = A(NE+KL) + A(EYE)
               END IF
               EYE = NXTEYE
  290       CONTINUE
  300    CONTINUE
  310 CONTINUE
      KEEP(IPTRD+N+1) = K
      KEEP(IPTRO+N+1) = KK
      NZD = K - 1
      DO 320 I = 1,K-1
         IRN(I) = IRN(NE+I)
  320 CONTINUE
      DO 325 I = KK+1,NE
         IRN(I) = IRN(NE+I)
  325 CONTINUE
      DO 326 I = K,KK
         IRN(I) = 0
  326 CONTINUE
C
C Repartition IW (none of previous data is needed any more)
      IP = 0
      IQ = M
      PTRD = M + N
      JFIRST = M + 2*N
      LENR = 2*M + 2*N
      LASTR = 3*M + 2*N
      NEXTR = 4*M + 2*N
      IW50 = 5*M + 2*N
      NEXTC = 6*M + 2*N
      PTRO = NEXTC
C Use scratch space in KEEP for work arrays in MA50AD.
      LASTC = M + N
      LENC = LASTC + N
C
C Call MA50A/AD block by block from the back.
      J1 = N + 1
C KB counts number of non-triangular blocks for calculation INFO(4)
      KB = 0
      DO 390 JB = NB,1,-1
C NC is number of columns in block
         NC = KEEP(NBLOCK+3*JB)
         J2 = J1 - 1
         J1 = J2 + 1 - NC
         IF (KEEP(MBLOCK+3*JB).LT.0) THEN
C Action if triangular block
            DO 330 J = J1,J2
C Increment estimate of rank
               IF (ABS(A(NE+KEEP(IPTRD+J))).GT.TOL)
     +             INFO(5) = INFO(5) + 1
               IW(IP+J) = J
               IW(IQ+J) = J
  330       CONTINUE
         ELSE
C Action if non-triangular block
C Copy the indices to IRN, shifting them.
            NZB = KEEP(IPTRD+J2+1) - KEEP(IPTRD+J1)
            DO 340 K = KEEP(IPTRD+J1),KEEP(IPTRD+J2+1) - 1
               IRN(NE+K) = IRN(K) - J1 + 1
  340       CONTINUE
C Copy the column pointers, shifting them.
C K is position immediately before current block in IRN
            K = KEEP(IPTRD+J1) - 1
            DO 350 J = J1,J2
               IW(IQ+J) = KEEP(IPTRD+J) - K
  350       CONTINUE
C NR is number of rows in block
            NR = NC
            IF (NB.EQ.1) NR = M
C Set permutations
C Set to identity because matrix has already been permuted to put
C   pivots on the diagonal.
            IF (JOB.EQ.2) THEN
               DO 360 J = J1,J1 + NR - 1
                  IW(IP+J) = J - J1 + 1
  360          CONTINUE
               DO 370 J = J1,J2
                  IW(PTRD+J) = J - J1 + 1
  370          CONTINUE
            END IF
C Accumulate data on block triangular form
            INFO(7) = MAX(INFO(7),NR)
            INFO(8) = INFO(8) + NC
            INFO(9) = INFO(9) + NZB
C Control for forcing columns to end.
            IF (ICNTL(8).NE.0) ICNTL5(6) = KEEP(LBLOCK+JB)
            CALL MA50AD(NR,NC,NZB,LA-NE-K,A(NE+K+1),IRN(NE+K+1),
     +                  JCN(NE+1),IW(IQ+J1),CNTL5,ICNTL5,IW(IP+J1),
     +                  NP,IW(JFIRST+1),IW(LENR+1),
     +                  IW(LASTR+1),IW(NEXTR+1),IW(IW50+1),IW(PTRD+J1),
     +                  KEEP(LENC+KB+1),KEEP(LASTC+KB+1),IW(NEXTC+1),
     +                  INFO5,RINFO5)
C KEEP(MBLOCK+3*JB) is number rows in packed storage
            KEEP(MBLOCK+3*JB) = NP
C Shift the permutations back
            DO 380 J = J1,J1+NR-1
               IW(IP+J) = IW(IP+J) + J1 - 1
  380       CONTINUE
            DO 385 J = J1,J2
               IW(IQ+J) = IW(IQ+J) + J1 - 1
  385       CONTINUE
C Adjust warning and error flag
            IF (INFO5(1).EQ.1) THEN
               IF (INFO(1).EQ.0 .OR. INFO(1).EQ.4) INFO(1) = INFO(1) + 2
            END IF
            IF (INFO5(1).EQ.2) THEN
               IF (INFO(1).EQ.0 .OR. INFO(1).EQ.2) INFO(1) = INFO(1) + 4
            END IF
            IF (INFO5(1).EQ.3 .AND. INFO(1).GE.0) INFO(1) = 6
            IF (INFO5(1).EQ.-3) INFO(1) = -3
C Accumulate number of garbage collections
            INFO(2) = INFO(2) + INFO5(2)
C Accumulate space for subsequent analyse
            INFO(3) = MAX(INFO(3),NE+K+INFO5(3))

C Accumulate estimate of rank
            INFO(5) = INFO(5) + INFO5(5)
C Accumulate number of dropped entries
            INFO(6) = INFO(6) + INFO5(6)
C Accumulate floating-point count
            RINFO(1) = RINFO(1) + RINFO5(1)
C Calculate data necessary for INFO(4)
C KEEP(LENC+KB) holds the actual storage forecast for MA50BD entry.
            KB = KB + 1
            KEEP(LENC+KB) = INFO5(4) - 2*NR
C KEEP(LASTC+KB) holds the storage plus elbow room forecast
C         for MA50BD entry.
            KEEP(LASTC+KB) = INFO5(4)
         END IF
  390 CONTINUE

C Now calculate storage required for subsequent factorization.
      INFO(4) = NE*2
      K = NE
      DO 400 JB = KB,1,-1
         INFO(4) = MAX(INFO(4),K+KEEP(LASTC+JB))
         K = K + KEEP(LENC+JB)
  400 CONTINUE

      IF (INFO(1).EQ.-3) THEN
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9120) LA,INFO(3)
         GO TO 530
      END IF

C
C Reorder IRN(1:NE) to the revised column permutation, storing in
C     IRN(NE+k) the new position for the entry that was in position k.
C Copy IRN data back.
      DO 410 K = 1,NE
         IRN(NE+K) = IRN(K)
  410 CONTINUE
      DO 420 J = 1,N
         IW(PTRD+J) = KEEP(IPTRD+J)
         IW(PTRO+J) = KEEP(IPTRO+J+1) + 1
  420 CONTINUE
      KD = 1
      KO = NZD + 1
      DO 450 J = 1,N
         KEEP(IPTRD+J) = KD
         JAY = IW(IQ+J)
         KL = NZD
         IF (JAY.NE.N) KL = IW(PTRD+JAY+1) - 1
         DO 430 KK = IW(PTRD+JAY),KL
            IRN(KD) = IW(IP+IRN(NE+KK))
            IRN(NE+KK) = KD
            KD = KD + 1
  430    CONTINUE
         KEEP(IPTRO+J) = KO
         KL = NE
         IF (JAY.NE.1) KL = IW(PTRO+JAY-1) - 1
         DO 440 KK = IW(PTRO+JAY),KL
            IRN(KO) = IW(IP+IRN(NE+KK))
            IRN(NE+KK) = KO
            KO = KO + 1
  440    CONTINUE
  450 CONTINUE
      KEEP(IPTRO+N+1) = KO
C
C Compute the product permutations
      DO 460 I = 1,M
         KEEP(I) = IW(IP+KEEP(I))
  460 CONTINUE
      DO 465 I = 1,N
         IW(IQ+I) = KEEP(M+IW(IQ+I))
  465 CONTINUE
      DO 470 I = 1,N
         KEEP(M+I) = IW(IQ+I)
  470 CONTINUE
C
C Compute (product) map
C Set IRN(NE) to deal with out-of-range indices
      IW(1) = IRN(NE)
      IRN(NE) = NE
      DO 480 K = 1,NE
         JCN(K) = IRN(NE+JCN(K))
  480 CONTINUE
C Reset IRN(NE)
      IRN(NE) = IW(1)

C Warning messages
      IF (INFO(11).GT.0 .OR. INFO(12).GT.0) THEN
         INFO(1) = INFO(1) + 1
         IF (MP.GT.0 .AND. ICNTL(3).GE.2) THEN
            IF (INFO(11).GT.0) WRITE (MP,9150) INFO(11)
            IF (INFO(12).GT.0) WRITE (MP,9160) INFO(12)
         END IF
         IF (INFO(11).GT.0) JCN(1) = -JCN(1)
      END IF

C Set rank estimate to min of structural and INFO(5) estimate.
      IF (INFO(10).LT.INFO(5)) THEN
         INFO(5) = INFO(10)
         IF (INFO(1).NE.2 .AND. INFO(1).NE.3 .AND.
     +       INFO(1).LT.6) INFO(1) = INFO(1) + 2
      END IF

      IF (MP.GT.0 .AND. ICNTL(3).GT.2) THEN
         WRITE (MP,'(/A/A/12I6/A,F12.1)') ' Leaving MA48A/AD with',
     +     ' INFO(1:12)  =',(INFO(I),I=1,12),' RINFO(1) =',RINFO(1)
         WRITE (MP,'(A)') ' Permuted matrix by blocks (IRN)'
         KB = NB
         IF (ICNTL(3).EQ.3) THEN
            WRITE (MP,'(A)')
     +        ' Only first column of up to 10 blocks printed'
            KB = MIN(10,NB)
         END IF
         WRITE (MP,'(A)') ' Diagonal blocks'
         J1 = 1
         DO 500 JB = 1,KB
            J2 = J1 + KEEP(NBLOCK+3*JB) - 1
            IF (J1.LE.J2) WRITE (MP,'(A,I6)') ' Block',JB
            J3 = J2
            IF (ICNTL(3).EQ.3) J3 = J1
            DO 490 J = J1,J3
               WRITE (MP,'(A,I5,(T13,10I6))') ' Column',J,
     +           (IRN(I),I=KEEP(IPTRD+J),KEEP(IPTRD+J+1)-1)
  490       CONTINUE
            J1 = J2 + 1
  500    CONTINUE
         IF (KEEP(IPTRO+N+1).GT.KEEP(IPTRD+N+1)) THEN
            WRITE (MP,'(A)') ' Off-diagonal entries'
            J1 = 1
            DO 520 JB = 1,KB
               J2 = J1 + KEEP(NBLOCK+3*JB) - 1
               J3 = J2
               IF (ICNTL(3).EQ.3) J3 = J1
               DO 510 J = J1,J3
                  IF (KEEP(IPTRO+J+1).GT.KEEP(IPTRO+J)) WRITE (MP,
     +                '(A,I5,(T13,10I6))') ' Column',J,
     +                (IRN(I),I=KEEP(IPTRO+J),KEEP(IPTRO+J+1)-1)
  510          CONTINUE
               J1 = J2 + 1
  520       CONTINUE
         END IF
         IF (ICNTL(3).GT.3) THEN
            WRITE (MP,9040) (JCN(K),K=1,NE)
 9040       FORMAT (' JCN (MAP) ='/ (6X,10I6))
            WRITE (MP,'(A)') ' Contents of KEEP array'
            WRITE (MP,9010) (KEEP(I),I=1,M)
            WRITE (MP,9020) (KEEP(M+I),I=1,N)
            WRITE (MP,'(A)') ' Pointer information from KEEP'
            WRITE (MP,9050) (KEEP(IPTRD+J),J=1,N+1)
 9050       FORMAT (' IPTRD ='/ (8X,10I6))
            WRITE (MP,9060) (KEEP(IPTRO+J),J=1,N+1)
 9060       FORMAT (' IPTRO ='/ (8X,10I6))
            WRITE (MP,'(A)') ' Block structure information from KEEP'
            WRITE (MP,9070) (KEEP(NBLOCK+3*JB),JB=1,NB)
            WRITE (MP,9080) (KEEP(MBLOCK+3*JB),JB=1,NB)
 9070       FORMAT (' NBLOCK (order blocks) ='/ (8X,10I6))
 9080       FORMAT (' MBLOCK (triangular flag and number packed rows) ='
     +             / (8X,10I6))
 9090       FORMAT (' LBLOCK (number of changed columns) ='/ (8X,10I6))
            IF (ICNTL(8).NE.0) WRITE (MP,9090) (KEEP(LBLOCK+JB),JB=1,NB)
         ELSE
            WRITE (MP,9040) (JCN(K),K=1,MIN(10,NE))
            WRITE (MP,'(A)') ' Contents of KEEP array'
            WRITE (MP,9010) (KEEP(I),I=1,MIN(10,M))
            WRITE (MP,9020) (KEEP(M+I),I=1,MIN(10,N))
            WRITE (MP,'(A)') ' Pointer information from KEEP'
            WRITE (MP,9050) (KEEP(IPTRD+J),J=1,MIN(10,N+1))
            WRITE (MP,9060) (KEEP(IPTRO+J),J=1,MIN(10,N+1))
            WRITE (MP,'(A)') ' Block structure information from KEEP'
            WRITE (MP,9070) (KEEP(NBLOCK+3*JB),JB=1,MIN(10,NB))
            WRITE (MP,9080) (KEEP(MBLOCK+3*JB),JB=1,MIN(10,NB))
            IF (ICNTL(8).NE.0) WRITE (MP,9090) (KEEP(LBLOCK+JB),JB=1,
     +          MIN(10,NB))
         END IF
      END IF
  530 RETURN

 9100 FORMAT (' Error return from MA48A/AD because M =',I10,' and N =',
     +       I10)
 9110 FORMAT (' Error return from MA48A/AD because NE =',I10)
 9120 FORMAT (' Error return from MA48A/AD because LA is',I10/' and ',
     +       'must be at least',I10)
 9130 FORMAT (' Error return from MA48A/AD because ','JOB = ',I10)
 9140 FORMAT (' Error return from MA48A/AD because ','faulty permutati',
     +       'ons input when JOB = 2')
 9150 FORMAT (' Message from MA48A/AD ..',I8,' duplicates found')
 9160 FORMAT (' Message from MA48A/AD ..',I8,' out-of-range indices fo',
     +       'und')
      END


      SUBROUTINE MA48BD(M,N,NE,JOB,LA,A,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,
     +                  INFO,RINFO)
C Factorize a sparse matrix, using data provided by MA48A/AD.
C
C     .. Arguments ..
      INTEGER M,N,NE,JOB,LA
      DOUBLE PRECISION A(LA)
      INTEGER IRN(LA),JCN(NE),KEEP(*)
      DOUBLE PRECISION CNTL(10)
      INTEGER ICNTL(20)
      DOUBLE PRECISION W(M)
      INTEGER IW(2*M+2*N),INFO(20)
      DOUBLE PRECISION RINFO(10)
C
C M must be set by the user to the number of rows.
C      It is not altered by the subroutine. Restriction: M > 0.
C N must be set by the user to the number of columns.
C      It is not altered by the subroutine. Restriction: N > 0.
C NE must be set by the user to the number of entries in the input
C      matrix. It is not altered by the subroutine. Restriction: NE > 0.
C JOB  must be set by the user to 1 for a normal entry or to 2 for a
C      fast entry. JOB is not altered by the subroutine.
C      Restriction 1 <= JOB <= 2.
C LA must be set by the user to the size of A and IRN.
C      It is not altered by the subroutine. Restriction: LA > 2*NE.
C A    must be set by the user so that the first NE elements hold the
C      matrix entries in exactly the same order as were input to
C      MA48A/AD. On return, these elements hold the permuted matrix A
C      and the rest contains the factorizations of the block diagonal
C      part (excluding blocks that are triangular).
C IRN  The first KEEP(IPTRO+N+1)-1 (<=NE) entries must be as on return
C      from MA48A/AD and hold
C      the row numbers of the permuted matrix, with duplicates excluded.
C      The permuted matrix is block upper triangular with recommended
C      pivots on its diagonal. The rest of the array need not be set on
C      entry with JOB=1. It is used for the row indices of the
C      factorizations. On entry with JOB=2, it must be unchanged
C      since the entry with JOB=1 and is not altered.
C JCN  must be as on return from MA48A/AD. abs(JCN(k)) holds the
C      position in IRN that corresponds to the entry that was input in
C      A(k),k=1,NE.  It is not altered by the subroutine.
C      If JCN(1) is negative then there are duplicates.
C KEEP must be as on return from MA48A/AD or MA48B/BD.
C      KEEP(i) holds the position of row i in the permuted
C        matrix, I=1,M;
C      KEEP(N+j) holds the index of the  column that is in position j
C        of the permuted matrix, j=1,N;
C      KEEP(IPTRL+k), IPTRL=M+N, k=1,N, need not be set on an
C        entry with JOB=1. It holds pointers, relative to the start
C        of the block, to the columns of the
C        lower-triangular part of the factorization. On entry with
C        JOB=2, it must be unchanged since the entry with JOB=1 and is
C        not altered.
C      KEEP(IPTRU+k), IPTRU=IPTRL+N, k=1,N, need not be set on an
C        entry with JOB=1. It holds pointers, relative to the start
C        of the block, to the columns of the
C        upper-triangular part of the factorization. On an entry with
C        JOB=2, it must be unchanged since the entry with JOB=1 and
C        is not altered.
C      KEEP(IPTRD+j), IPTRD=IPTRU+N, holds the position in IRN of the
C        start of the block diagonal part of column j, J=1,N;
C      KEEP(IPTRD+N+1) holds the position that immediately follows
C        the end of the block diagonal part of column N;
C      KEEP(IPTRO+j),IPTRO=IPTRD+N+1, holds the position in IRN of
C        the start of the block off-diagonal part of column j, j=1,N;
C      KEEP(IPTRO+N+1) holds the position that immediately
C        follows the end of the block off-diagonal part of column N;
C      KEEP(KBLOCK+3), KBLOCK=IPTRO+N+1, holds the number of
C        blocks NB in the block triangular form;
C      KEEP(NBLOCK+3*k), NBLOCK=IPTRO+N-1, holds the number of columns
C        in block k, k=1,NB;
C      KEEP(MBLOCK+3*k), MBLOCK=IPTRO+N, is negative if block k
C        is triangular or holds the number of rows held in packed
C        storage when processing block k, k=1,NB;
C      KEEP(KBLOCK+3*k), k=2,NB need not be set on an entry with
C        JOB=1; on return it holds zero for a triangular block and for
C        a non-triangular block holds the position in A(NEWNE+1) and
C        IRN(NEWNE+1) of the start of the factorization of the block.
C        KEEP(LBLOCK+k), LBLOCK=KBLOCK+3*NB is accessed only if ICNTL(8)
C          is not 0 and holds the number of columns in
C          block k which are to be pivoted on last, k=1,NB.
C      On an entry with JOB=2, KEEP must be unchanged since the entry
C      with JOB=1 and is not altered.
C CNTL  must be set by the user as follows and is not altered.
C     CNTL(2) determines the balance used by MA50A/AD and MA50B/BD
C       between pivoting for sparsity and for stability, values near
C       zero emphasizing sparsity and values near one emphasizing
C       stability.
C     CNTL(3) If this is set to a positive value, any entry whose
C       modulus is less than CNTL(3) will be dropped from the factors
C       calculated by MA50A/AD. The factorization will then require
C       less storage but will be inaccurate.
C     CNTL(4)  If this is set to a positive value, any entry whose
C       modulus is less than CNTL(4) will be regarded as zero from
C       the point of view of rank.
C     CNTL(5:10) are not referenced by the subroutine.
C ICNTL must be set by the user as follows and is not altered.
C     ICNTL(1)  must be set to the stream number for error messages.
C       A value less than 1 suppresses output.
C     ICNTL(2) must be set to the stream number for diagnostic output.
C       A value less than 1 suppresses output.
C     ICNTL(3) must be set to control the amount of output:
C       0 None.
C       1 Error messages only.
C       2 Error and warning messages.
C       3 As 2, plus scalar parameters and a few entries of array
C         parameters on entry and exit.
C       4 As 2, plus all parameters on entry and exit.
C     ICNTL(4)  If set to a positive value, the pivot search by MA50A/AD
C       is limited to ICNTL(4) columns. This may result in different
C       fill-in and execution time.
C     ICNTL(5) The block size to be used for full-matrix processing.
C     ICNTL(6) is not referenced by the subroutine.
C     ICNTL(8) If set to a value other than zero and JOB = 2 or 3,
C       it is assumed that only the last KEEP(LBLOCK+JB), JB = 1,NB
C       columns have changed since the last factorization.
C     ICNTL(9) is not referenced by the subroutine.
C     ICNTL(10) has default value 0. If set to 1, there is an immediate
C       return from MA48B/BD if LA is too small, without continuing the
C       decomposition to compute the size necessary.
C     ICNTL(11) has default value 0. If set to 1 on a JOB=2 call to
C       MA48B/BD and the entries in one of the blocks on the diagonal
C       are unsuitable for the pivot sequence chosen on the previous
C       call, the block is refactorized as on a JOB=1 call.
C     ICNTL(12:20) are not referenced by the subroutine.

C W    is a workarray.
C IW   is a workarray.
C INFO(1) is an integer variable that need not be set on entry. On exit,
C      a nonzero value of INFO(1) will indicate an error return.
C      Possible nonzero values are:
C      -1  M < 1 or N < 1
C      -2  NE < 0
C      -3  Insufficient space
C      -6  JOB out of range or JOB = 2 or 3 after factorization in which
C          entries were dropped.
C      -7  On a call with JOB=2, the matrix entries are numerically
C          unsuitable
C      +2  Matrix is structurally rank deficient. Estimated rank in
C          INFO(5).
C INFO(4) is set to the minimum size required for a further
C      factorization on a matrix with the same pivot sequence.
C INFO(5) is set to computed rank of the matrix.
C INFO(6) is set to number of entries dropped from structure.
C RINFO need not be set on entry. On exit, RINFO(1) holds the number of
C    floating-point operations performed.

C     .. Local constants ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.D0)
C
C     .. Local variables ..
      DOUBLE PRECISION CNTL5(10)
      INTEGER I,ICNTL5(20),INFO5(15),IPTRD,IPTRL,IPTRO,IPTRU,IQB(1),
     +        ITRY,J,JB,JOB5,J1,J2,J3,K,KB,KBLOCK,KK,LBLOCK,LP,
     +        MBLOCK,MP,NB,NBLOCK,NEWNE,NC,NP,NR,NRF,NZB
      DOUBLE PRECISION RINFO5(10),TOL
      LOGICAL TRISNG
C I Temporary DO index.
C ICNTL5 passed to MA50B/BD to correspond to dummy argument ICNTL.
C INFO5 passed to MA50B/BD to correspond to dummy argument INFO.
C IPTRD displacement in KEEP. See comment on KEEP.
C IPTRL displacement in KEEP. See comment on KEEP.
C IPTRO displacement in KEEP. See comment on KEEP.
C IPTRU displacement in KEEP. See comment on KEEP.
C IQB   passed to MA50B/BD to indicate that no column permutation is
C       required.
C ITRY  Loop index for retrying the factorization of a block.
C J     column index.
C JB    block index.
C JOB5  passed to MA50B/BD to correspond to dummy argument JOB
C J1    first column index of a block.
C J2    last column index of a block.
C J3    used in prints to hold last column index of block.
C K     running index for position in matrix.
C KB  number of blocks ... used in printing
C KBLOCK displacement in KEEP. See comment on KEEP.
C KK    running index for position in matrix.
C LBLOCK displacement in KEEP. See comment on KEEP.
C LP Unit for error messages.
C MBLOCK displacement in KEEP. See comment on KEEP.
C MP Unit for diagnostic messages.
C NBLOCK displacement in KEEP. See comment on KEEP.
C NB  number of blocks
C NC    number of columns in block
C NEWNE number of entries in original matrix with duplicates and
C       out-of-range indices removed
C NP Number of columns in packed storage.
C NR    number of rows in block
C NRF   is number of rows in full form for current block.
C NZB   number of entries in current diagonal block
C RINFO5 passed to MA50B/BD to correspond to dummy argument RINFO.
C TOL   pivot tolerance
C TRISNG is flag to indicate singularity in triangular block.

C     .. Externals ..
      EXTERNAL MA50BD
      INTRINSIC MAX

      LP = ICNTL(1)
      MP = ICNTL(2)
C Simple data checks
      IF (M.LE.0 .OR. N.LE.0) THEN
         INFO(1) = -1
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9160) M,N
         GO TO 240
      END IF
      IF (NE.LE.0) THEN
         INFO(1) = -2
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9170) NE
         GO TO 240
      END IF
      IF (LA.LT.2*NE) THEN
         INFO(1) = -3
         INFO(4) = 2*NE
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9180) LA,INFO(4)
         GO TO 240
      END IF
      IF (JOB.LT.1 .OR. JOB.GT.3) THEN
         INFO(1) = -6
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9190) JOB
         GO TO 240
      END IF
      INFO(1) = 0
      DO 10 I = 1,4
         CNTL5(I) = CNTL(I)
         ICNTL5(I) = ICNTL(I)
   10 CONTINUE
C Switch off printing from MA50B/BD
      ICNTL5(3) = 0
C Set BLAS control
      ICNTL5(5) = ICNTL(5)
C Initialize restricted pivoting control
      ICNTL5(6) = 0
      ICNTL5(7) = 0
      ICNTL5(8) = ICNTL(10)

C Partition KEEP and A
      IPTRL = M + N
      IPTRU = IPTRL + N
      IPTRD = IPTRU + N
      IPTRO = IPTRD + N + 1
      NBLOCK = IPTRO + N - 1
      MBLOCK = NBLOCK + 1
      KBLOCK = MBLOCK + 1
      NB = KEEP(KBLOCK+3)
      LBLOCK = KBLOCK + 3*NB
      NEWNE = KEEP(IPTRO+N+1) - 1

      IF (MP.GT.0 .AND. ICNTL(3).GT.2) THEN
         WRITE (MP,'(/A/3(A,I8),A,I2,A,I8/A,1P,3D12.4/A,3I8/A,I8/A,I8)')
     +     ' Entering MA48B/BD with',' M =',M,'     N =',N,'     NE ='
     +     ,NE,'     JOB =',JOB,'     LA =',LA,' CNTL (2:4) =',
     +     (CNTL(I),I=2,4),' ICNTL(1:3) =', (ICNTL(I),I=1,3),
     +     ' ICNTL(5)   =',ICNTL(5),' ICNTL(8)   =',ICNTL(8)
         IF (ICNTL(3).GT.3) THEN
            WRITE (MP,9000) (A(K),K=1,NE)
         ELSE
            WRITE (MP,9000) (A(K),K=1,MIN(10,NE))
         END IF
 9000    FORMAT (' A ='/ (4X,1P,5D12.4))
         WRITE (MP,'(A)') ' Indices for permuted matrix by blocks'
         KB = NB
         IF (ICNTL(3).EQ.3) THEN
            WRITE (MP,'(A)')
     +        ' Only first column of up to 10 blocks printed'
            KB = MIN(10,NB)
         END IF
         WRITE (MP,'(A)') ' Diagonal blocks'
         J1 = 1
         DO 30 JB = 1,KB
            WRITE (MP,'(A,I6)') ' Block',JB
            J2 = J1 + KEEP(NBLOCK+3*JB) - 1
            J3 = J2
            IF(ICNTL(3).EQ.3) J3 = J1
            DO 20 J = J1,J3
               WRITE (MP,'(A,I5,(T13,10I6))') ' Column',J,
     +           (IRN(I),I=KEEP(IPTRD+J),KEEP(IPTRD+J+1)-1)
   20       CONTINUE
            J1 = J2 + 1
   30    CONTINUE
         IF (KEEP(IPTRO+N+1).GT.KEEP(IPTRD+N+1)) THEN
            WRITE (MP,'(A)') ' Off-diagonal entries'
            J1 = 1
            DO 50 JB = 1,KB
               J2 = J1 + KEEP(NBLOCK+3*JB) - 1
               J3 = J2
               IF(ICNTL(3).EQ.3) J3 = J1
               DO 40 J = J1,J3
                  IF (KEEP(IPTRO+J+1).GT.KEEP(IPTRO+J)) WRITE (MP,
     +                '(A,I5,(T13,10I6))') ' Column',J,
     +                (IRN(I),I=KEEP(IPTRO+J),KEEP(IPTRO+J+1)-1)
   40          CONTINUE
               J1 = J2 + 1
   50       CONTINUE
         END IF
         IF (ICNTL(3).GT.3) THEN
            WRITE (MP,9010) (JCN(K),K=1,NE)
 9010       FORMAT (' JCN (MAP) ='/ (6X,10I6))
            WRITE (MP,'(A)') ' Contents of KEEP array'
            WRITE (MP,9020) (KEEP(I),I=1,M)
            WRITE (MP,9030) (KEEP(M+I),I=1,N)
 9020       FORMAT (' Positions of original rows in the permuted matrix'
     +             / (10I6))
 9030       FORMAT (' Positions of columns of permuted matrix ','in ',
     +             'original matrix '/ (10I6))
            WRITE (MP,'(A)') ' Pointer information from KEEP'
            WRITE (MP,9040) (KEEP(IPTRD+J),J=1,N+1)
 9040       FORMAT (' IPTRD ='/ (8X,10I6))
            WRITE (MP,9050) (KEEP(IPTRO+J),J=1,N+1)
 9050       FORMAT (' IPTRO ='/ (8X,10I6))
            IF (JOB.GT.1) THEN
               WRITE (MP,9060) (KEEP(IPTRL+J),J=1,N)
 9060          FORMAT (' IPTRL ='/ (8X,10I6))
               WRITE (MP,9070) (KEEP(IPTRU+J),J=1,N)
 9070          FORMAT (' IPTRU ='/ (8X,10I6))
            END IF
            WRITE (MP,'(A)') ' Block structure information from KEEP'
            WRITE (MP,9080) (KEEP(NBLOCK+3*JB),JB=1,NB)
            WRITE (MP,9090) (KEEP(MBLOCK+3*JB),JB=1,NB)
 9080       FORMAT (' NBLOCK (order blocks) ='/ (8X,10I6))
 9090       FORMAT (' MBLOCK (triangular flag and number packed rows) ='
     +             / (8X,10I6))
 9100       FORMAT (' KBLOCK (position of beginning of block) ='/
     +             (8X,10I6))
 9110       FORMAT (' LBLOCK (number of changed columns) ='/ (8X,10I6))
            IF (JOB.GT.1) THEN
               WRITE (MP,9100) (KEEP(KBLOCK+3*JB),JB=1,NB)
               IF (ICNTL(8).NE.0) WRITE (MP,9110) (KEEP(LBLOCK+JB),JB=1,
     +             NB)
            END IF
         ELSE
            WRITE (MP,9010) (JCN(K),K=1,MIN(10,NE))
            WRITE (MP,'(A)') ' Contents of KEEP array'
            WRITE (MP,9020) (KEEP(I),I=1,MIN(10,M))
            WRITE (MP,9030) (KEEP(M+I),I=1,MIN(10,N))
            WRITE (MP,'(A)') ' Pointer information from KEEP'
            WRITE (MP,9040) (KEEP(IPTRD+J),J=1,MIN(10,N+1))
            WRITE (MP,9050) (KEEP(IPTRO+J),J=1,MIN(10,N+1))
            IF (JOB.GT.1) THEN
               WRITE (MP,9060) (KEEP(IPTRL+J),J=1,MIN(10,N))
               WRITE (MP,9070) (KEEP(IPTRU+J),J=1,MIN(10,N))
            END IF
            WRITE (MP,'(A)') ' Block structure information from KEEP'
            WRITE (MP,9080) (KEEP(NBLOCK+3*JB),JB=1,MIN(10,NB))
            WRITE (MP,9090) (KEEP(MBLOCK+3*JB),JB=1,MIN(10,NB))
            IF (JOB.GT.1) THEN
               WRITE (MP,9100) (KEEP(KBLOCK+3*JB),JB=1,MIN(10,NB))
               IF (ICNTL(8).NE.0) WRITE (MP,9110) (KEEP(LBLOCK+JB),JB=1,
     +             MIN(10,NB))
            END IF
         END IF
      END IF
C
C Initialize INFO and RINFO
      INFO(4) = NE
      INFO(5) = 0
      INFO(6) = 0
      RINFO(1) = ZERO
C Set pivot tolerance
      TOL = MAX(ZERO,CNTL(4))
C
C Use map to sort the new values into A.
C Mapping into first NEWNE locations in array A
      IF (JCN(1).GT.0) THEN
         DO 60 K = 1,NE
            A(NE+K) = A(K)
   60    CONTINUE
CDIR$ IVDEP
         DO 70 K = 1,NE
            A(JCN(K)) = A(NE+K)
   70    CONTINUE
      ELSE
C Duplicates exist
         DO 80 K = 1,NE
            A(NE+K) = A(K)
            A(K) = ZERO
   80    CONTINUE
         A(-JCN(1)) = A(NE+1)
         DO 90 K = 2,NE
            KK = JCN(K)
            A(KK) = A(KK) + A(NE+K)
   90    CONTINUE
      END IF
C
C Call MA50B/BD block by block.
      IQB(1) = 0
      KK = 0
      J2 = 0
      JOB5 = JOB
      DO 150 JB = 1,NB
         NC = KEEP(NBLOCK+3*JB)
         J1 = J2 + 1
         J2 = J1 + NC - 1
         KEEP(KBLOCK+3*JB) = 0
         IF (KEEP(MBLOCK+3*JB).LT.0) THEN
C Action if triangular block
            TRISNG = .FALSE.
            DO 100 J = J1,J2
               IF (ABS(A(KEEP(IPTRD+J))).LE.TOL) TRISNG = .TRUE.
               KEEP(IPTRL+J) = 0
               KEEP(IPTRU+J) = 0
  100       CONTINUE
            IF (.NOT.TRISNG)  THEN
               INFO(5) = INFO(5) + NC
               GO TO 150
            ENDIF
C Block is singular. Treat as non-triangular.
            IF (JOB.EQ.2) THEN
               INFO(1) = -7
               IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,'(A)')
     +              ' Error return from MA48B/BD with JOB=2 because ',
     +              ' the matrix is incompatible with expectations'
               GO TO 240
            ELSE
               KEEP(MBLOCK+3*JB) = NC
            ENDIF
         END IF
C Action if non-triangular block
       DO 145 ITRY = 1,2
C The second iteration of this loop is used only if JOB=2,
C   ICNTL(11)=1, and the call to MA50B with JOB=1 fails.
         NR = NC
         IF (NB.EQ.1) NR = M
         NZB = KEEP(IPTRD+J2+1) - KEEP(IPTRD+J1)
         IF (ICNTL(8).NE.0) ICNTL5(6) = KEEP(LBLOCK+JB)
C Shift the indices and pointers to local values.
         DO 110 K = KEEP(IPTRD+J1),KEEP(IPTRD+J2+1) - 1
            IRN(K) = IRN(K) - J1 + 1
  110    CONTINUE
         K = KEEP(IPTRD+J1) - 1
         DO 115 J = J1,J1+NC-1
            KEEP(IPTRD+J) = KEEP(IPTRD+J) - K
  115    CONTINUE
         DO 120 J = J1,J1+NR-1
            IW(J) = J - J1 + 1
  120    CONTINUE
         NP = KEEP(MBLOCK+3*JB)
         CALL MA50BD(NR,NC,NZB,JOB5,A(K+1),IRN(K+1),KEEP(IPTRD+J1),
     +               CNTL5,ICNTL5,IW(J1),IQB,NP,LA-NEWNE-KK,
     +               A(NEWNE+KK+1),IRN(NEWNE+KK+1),KEEP(IPTRL+J1),
     +               KEEP(IPTRU+J1),W,IW(M+1),INFO5,RINFO5)
C Restore the indices and pointers
         DO 130 J = J1,J2
            KEEP(IPTRD+J) = KEEP(IPTRD+J) + K
  130    CONTINUE
         DO 140 K = KEEP(IPTRD+J1),KEEP(IPTRD+J2+1) - 1
            IRN(K) = IRN(K) + J1 - 1
  140    CONTINUE
C Set warning and error returns
         IF (INFO5(1).EQ.-6) THEN
            INFO(1) = -6
            IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,'(A)')
     +          ' Error return from MA48B/BD with JOB greater than 1',
     +          ' and entries dropped during previous factorization'
            GO TO 240
         END IF
         IF (INFO5(1).LT.-7) THEN
            IF (ICNTL(11).EQ.1 .AND. JOB.EQ.2) THEN
               JOB5 = 1
               IF (LP.GT.0 .AND. ICNTL(3).GE.2) WRITE(LP,'(A,2(A,I4))')
     +          ' Warning from MA48B/BD. Switched from JOB=2 to JOB=1',
     +          ' in block', JB, ' of ',NB
               INFO(1) = INFO(1) + 1
               GO TO 145
            END IF
            INFO(1) = -7
            IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,'(A)')
     +          ' Error return from MA48B/BD with JOB=2 because ',
     +          ' the matrix is incompatible with expectations'
            GO TO 240
         ELSE
            GO TO 147
         END IF
  145  CONTINUE
  147    IF (INFO5(1).EQ.-3) THEN
            INFO(1) = -3
            IF (ICNTL(10).EQ.1) THEN
               KEEP(KBLOCK+3) = NB
               IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE(LP,'(A,2(A,I4))')
     +          ' Error return from MA48B/BD because LA is too small.',
     +          ' In block', JB, ' of ',NB
               GO TO 240
            END IF
         END IF
         IF (INFO(1).EQ.-3) THEN
            INFO(4) = INFO(4) + INFO5(4)
            KK = 0
         ELSE
            INFO(4) = MAX(INFO(4),KK+NEWNE+INFO5(4))
            NRF = IRN(NEWNE+KK+2)
C Set pointer to first entry in block JB in A(NEWNE+1), IRN(NEWNE+1)
            KEEP(KBLOCK+3*JB) = KK + 1
            KK = KK + KEEP(IPTRL+J2) + MAX((NC-KEEP(MBLOCK+3*JB))*
     +           (NRF), (NC-KEEP(MBLOCK+3*JB))+(NRF))
         END IF
C Is matrix rank deficient?
         IF (INFO5(1).EQ.1) THEN
            IF (INFO(1).NE.-3) INFO(1) = MIN(INFO(1)+2,3)
         END IF
C Accumulate stats
         RINFO(1) = RINFO(1) + RINFO5(1)
         INFO(5) = INFO(5) + INFO5(5)
         INFO(6) = INFO(6) + INFO5(6)
  150 CONTINUE

      INFO(4) = MAX(NE*2,INFO(4))
      KEEP(KBLOCK+3) = NB

      IF (INFO(1).EQ.-3) THEN
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9180) LA,INFO(4)
         GO TO 240
      END IF

      IF (MP.GT.0 .AND. ICNTL(3).GT.2) THEN
         WRITE (MP,'(/A/A,I6/A,3I6/A,F12.1)') ' Leaving MA48B/BD with',
     +     ' INFO(1)   = ',INFO(1),' INFO(4:6) = ', (INFO(I),I=4,6),
     +     ' RINFO(1)     =',RINFO(1)
         WRITE (MP,'(A)') ' Permuted matrix by blocks'
         KB = NB
         IF (ICNTL(3).EQ.3) THEN
            WRITE (MP,'(A)')
     +        ' Only first column of up to 10 blocks printed'
            KB = MIN(10,NB)
         END IF
         WRITE (MP,'(A)') ' Diagonal blocks'
         J1 = 1
         DO 170 JB = 1,KB
            J2 = J1 + KEEP(NBLOCK+3*JB) - 1
            IF (J1.LE.J2) WRITE (MP,'(A,I6)') ' Block',JB
            J3 = J2
            IF (ICNTL(3).EQ.3) J3 = MIN(J1,J2)
            DO 160 J = J1,J3
               WRITE (MP,'(A,I5,(T13,3(1PD12.4,I5)))') ' Column',J,
     +           (A(I),IRN(I),I=KEEP(IPTRD+J),KEEP(IPTRD+J+1)-1)
  160       CONTINUE
            J1 = J2 + 1
  170    CONTINUE
         IF (KEEP(IPTRO+N+1).GT.KEEP(IPTRD+N+1)) THEN
            WRITE (MP,'(A)') ' Off-diagonal entries'
            J1 = 1
            DO 190 JB = 1,KB
               J2 = J1 + KEEP(NBLOCK+3*JB) - 1
               J3 = J2
               IF (ICNTL(3).EQ.3) J3 = MIN(J1,J2)
               DO 180 J = J1,J3
                  IF (KEEP(IPTRO+J+1).GT.KEEP(IPTRO+J)) WRITE (MP,
     +                '(A,I5,(T13,3(1P,D12.4,I5)))') ' Column',J,
     +                (A(I),IRN(I),I=KEEP(IPTRO+J),KEEP(IPTRO+J+1)-1)
  180          CONTINUE
               J1 = J2 + 1
  190       CONTINUE
         END IF
         WRITE (MP,'(A)') ' Factorized matrix by blocks'
         J1 = 1
         DO 230 JB = 1,KB
            J2 = J1 + KEEP(NBLOCK+3*JB) - 1
C Jump if triangular block
            IF (KEEP(MBLOCK+3*JB).LT.0) GO TO 220
            NC = J2 - J1 + 1
            NR = NC
            IF (KB.EQ.1) NR = M
            WRITE (MP,'(A,I6)') ' Block',JB
            K = NEWNE
            IF (JB.GT.1) K = KEEP(KBLOCK+3*JB) + NEWNE - 1
            IF (KEEP(IPTRL+J1).GT.KEEP(IPTRU+J1)) WRITE (MP,
     +          '(A,I5,A,(T18,3(1P,D12.4,I5)))') ' Column',J1,' of L',
     +          (A(K+I),IRN(K+I),I=KEEP(IPTRU+J1)+1,KEEP(IPTRL+J1))
            IF (ICNTL(3).EQ.3) GO TO 210
            DO 200 J = J1 + 1,J2
               IF (KEEP(IPTRU+J).GT.KEEP(IPTRL+J-1)) WRITE (MP,
     +             '(A,I5,A,(T18,3(1P,D12.4,I5)))') ' Column',J,' of U',
     +              (A(K+I),IRN(K+I),I=KEEP(IPTRL+J-1)+1,KEEP(IPTRU+J))
               IF (KEEP(IPTRU+J).LT.KEEP(IPTRL+J)) WRITE (MP,
     +             '(A,I5,A,(T18,3(1P,D12.4,I5)))') ' Column',J,' of L',
     +              (A(K+I),IRN(K+I),I=KEEP(IPTRU+J)+1,KEEP(IPTRL+J))
  200       CONTINUE
C Full blocks
  210       WRITE (MP,'(A)') ' Full block'
            WRITE (MP,'(A)') ' Row indices'
            NRF = IRN(K+2)
            K = K + KEEP(IPTRL+J2)
            IF (ICNTL(3).GT.3) THEN
               WRITE (MP,9120) (IRN(K+I),I=1,NRF)
               WRITE (MP,'(A)') ' Column pivoting information'
               WRITE (MP,9120) (IRN(K+I),I=NRF+1,
     +           NRF+NC-KEEP(MBLOCK+3*JB))
               WRITE (MP,'(A)') ' Reals by columns'
               WRITE (MP,9130) (A(K+I),I=1,
     +           (NRF)* (NC-KEEP(MBLOCK+3*JB)))
 9120          FORMAT (10I6)
 9130          FORMAT (1P,5D12.4)
            ELSE
               WRITE (MP,9120) (IRN(K+I),I=1,MIN(10,NRF))
               WRITE (MP,'(A)') ' Column pivoting information'
               WRITE (MP,9120) (IRN(K+I),I=NRF+1,
     +           NRF+MIN(10,NC-KEEP(MBLOCK+3*JB)))
               WRITE (MP,'(A)') ' Reals by columns'
               WRITE (MP,9130) (A(K+I),I=1,
     +           MIN(10, (NRF)* (NC-KEEP(MBLOCK+3*JB))))
            END IF
  220       J1 = J2 + 1
  230    CONTINUE
         IF (JOB.EQ.1 .OR. JOB.EQ.3) THEN
            WRITE (MP,'(A)') ' Contents of KEEP array'
            IF (ICNTL(3).GT.3) THEN
               WRITE (MP,9020) (KEEP(I),I=1,M)
               WRITE (MP,9030) (KEEP(M+I),I=1,N)
               WRITE (MP,'(A)') ' Pointer information from KEEP'
               WRITE (MP,9140) (KEEP(IPTRL+J),J=1,N)
 9140          FORMAT (' IPTRL ='/ (8X,10I6))
               WRITE (MP,9150) (KEEP(IPTRU+J),J=1,N)
 9150          FORMAT (' IPTRU ='/ (8X,10I6))
            ELSE
               WRITE (MP,9020) (KEEP(I),I=1,MIN(10,M))
               WRITE (MP,9030) (KEEP(M+I),I=1,MIN(10,N))
               WRITE (MP,'(A)') ' Pointer information from KEEP'
               WRITE (MP,9140) (KEEP(IPTRL+J),J=1,MIN(10,N))
               WRITE (MP,9150) (KEEP(IPTRU+J),J=1,MIN(10,N))
            END IF
            WRITE (MP,'(A)') ' Block structure information from KEEP'
            WRITE (MP,9100) (KEEP(KBLOCK+3*JB),JB=1,KB)
         END IF
      END IF

  240 RETURN

 9160 FORMAT (' Error return from MA48B/BD because M =',I10,' and N =',
     +       I10)
 9170 FORMAT (' Error return from MA48B/BD because NE =',I10)
 9180 FORMAT (' Error return from MA48B/BD because LA is',I10/' and mu',
     +       'st be at least',I10)
 9190 FORMAT (' Error return from MA48B/BD because ','JOB = ',I10)
      END


      SUBROUTINE MA48CD(M,N,TRANS,JOB,LA,A,IRN,KEEP,CNTL,ICNTL,RHS,X,
     +                  ERROR,W,IW,INFO)

C Solve linear system, using data provided by MA48B/BD.
C
C     .. Arguments ..
      INTEGER M,N
      LOGICAL TRANS
      INTEGER JOB,LA
      DOUBLE PRECISION A(LA)
      INTEGER IRN(LA),KEEP(*)
      DOUBLE PRECISION CNTL(10)
      INTEGER ICNTL(20)
      DOUBLE PRECISION RHS(*),X(*),ERROR(3),W(*)
      INTEGER IW(*),INFO(20)
C
C M must be set by the user to the number of rows in the matrix.
C      It is not altered by the subroutine. Restriction: M > 0.
C N must be set by the user to the number of columns in the matrix.
C      It is not altered by the subroutine. Restriction: N > 0.
C TRANS must be set by the user to indicate whether coefficient matrix
C      is A (TRANS=.FALSE.) or A transpose (TRANS=.TRUE.).
C      It is not altered by the subroutine.
C JOB  must be set by the user to control the solution phase.
C      Restriction: 1 <= JOB <= 4
C      Possible values are:
C       =1 Returns solution only.
C       =2 Returns solution and estimate of backward error.
C       =3 Performs iterative refinement and returns estimate of
C          backward error.
C       =4 As for 3 plus an estimate of the error in the solution.
C      It is not altered by the subroutine.
C LA must be set by the user to the size of arrays A and IRN.
C      It is not altered by the subroutine.
C A    must be set left unchanged from the last call to MA48B/BD.
C      It holds the original matrix in permuted form and the
C      factorization of the block diagonal part (excluding
C      any triangular blocks). It is not altered by the subroutine.
C IRN  must be set left unchanged from the last call to MA48B/BD.
C      It holds the row indices of the factorization in A and
C      the row indices of the original matrix in permuted form.
C      It is not altered by the subroutine.
C KEEP must be as on return from MA48B/BD.
C      It is not altered by the subroutine.
C CNTL  must be set by the user as follows and is not altered.
C     CNTL(I), I=1,4 not accessed by subroutine.
C     CNTL(5) is used in convergence test for termination of iterative
C     refinement.  Iteration stops if successive omegas do not decrease
C     by at least this factor.
C ICNTL must be set by the user as follows and is not altered.
C     ICNTL(1)  must be set to the stream number for error messages.
C       A value less than 1 suppresses output.
C     ICNTL(2) must be set to the stream number for diagnostic output.
C       A value less than 1 suppresses output.
C     ICNTL(3) must be set to control the amount of output:
C       0 None.
C       1 Error messages only.
C       2 Error and warning messages.
C       3 As 2, plus scalar parameters and a few entries of array
C         parameters on entry and exit.
C       4 As 2, plus all parameters on entry and exit.
C     ICNTL(I), I=4,8 not accessed by subroutine
C     ICNTL(9) is maximum number of iterations allowed in iterative
C         refinement.
C RHS  must be set by the user to contain the right-hand side of the
C      equations to be solved. If TRANS is .FALSE., RHS has length M;
C      otherwise, it has length N. It is used as workspace.
C X    need not be set by the user.  On exit, it will contain the
C      solution vector. If TRANS is .FALSE., X has length N; otherwise,
C      it has length M.
C ERROR is a real array of length 3 that need not be set by the user.
C      On exit with JOB >=2 ERROR(1) and ERROR(2) give estimates of
C      backward errors OMEGA1 and OMEGA2.  On exit with JOB = 4,
C      ERROR(3) holds an estimate of the infinity norm in the relative
C      error in the solution.
C W    is a workarray.
C IW   is a workarray.
C INFO need not be set on entry. On exit, it holds the following:
C    INFO(1) A nonzero value will indicate an error return. Possible
C      nonzero values are:
C      -1  M or N < 1
C      -6  JOB out of range
C      -8  Nonconvergence of iterative refinement
C      -9  Failure in MC71A/AD

C     .. Local constants ..
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0.D0,ONE=1.0D0)
C
C     .. Local variables ..
      DOUBLE PRECISION COND(2),CTAU,DXMAX
      INTEGER I,ICNTL5(20),IPTRD,IPTRL,IPTRO,IPTRU,IQB(1),J,JB,JJ,J1,J2,
     +        J3,K,KASE,KB,KBLOCK,KK,KEEP71(5)
      LOGICAL LCOND(2)
      INTEGER LP,MBLOCK,MP,NB,NBLOCK,NC,NE,NEQ,NRF,NVAR
      DOUBLE PRECISION OLDOMG(2),OMEGA(2),OM1,OM2,TAU
C COND  Arioli, Demmel, and Duff condition number
C CTAU is a real variable used to control the splitting of the equations
C     into two categories. This is constant in Arioli, Demmel, and Duff
C     theory.
C DXMAX max-norm of current solution estimate
C I     row index and DO loop variable
C ICNTL5 passed to MA50C/CD to correspond to dummy argument ICNTL.
C IPTRD displacement in KEEP. See comment on KEEP.
C IPTRL displacement in KEEP. See comment on KEEP.
C IPTRO displacement in KEEP. See comment on KEEP.
C IPTRU displacement in KEEP. See comment on KEEP.
C IQB   is passed to MA50C/CD to indicate that no column permutation
C       is used.
C J     column index
C JB    is current block.
C JJ    running index for column.
C J1    is index of first column in block.
C J2    is index of last column in block.
C J3    is index of last column in block used in printing diagnostics.
C K     iteration counter and DO loop variable
C KASE  control for norm calculating routine MC71A/AD
C KB    used in prints to hold number of blocks or less if ICNTL(3)
C       equal to 3.
C KBLOCK displacement in KEEP. See comment on KEEP.
C KEEP71 workspace to preserve MC71A/AD's locals between calls
C KK    DO loop variable
C LCOND LCOND(k) is set to .TRUE. if there are equations in category
C       k, k=1,2
C LP Unit for error messages.
C MBLOCK displacement in KEEP. See comment on KEEP.
C MP Unit for diagnostic messages.
C NB    number of diagonal blocks
C NBLOCK displacement in KEEP. See comment on KEEP.
C NC    is number of columns in packed form for current block.
C NE    is set to the displacement in A/ICN for the beginning of
C       information on the factors.
C NEQ   number of equations in system
C NRF   is number of rows in full form for current block.
C NVAR  number of variables in system
C OLDOMG value of omega from previous iteration. Kept in case it is
C       better.
C OMEGA backward error estimates.
C OM1   value of OMEGA(1)+OMEGA(2) from previous iteration
C OM2   value of OMEGA(1)+OMEGA(2) from current iteration
C TAU   from Arioli, Demmel, and Duff .. used to divide equations
C
C     .. Externals ..
      EXTERNAL MA48DD,MA50CD,MC71AD
      DOUBLE PRECISION EPS
C EPS is the largest real such that 1+EPS is equal to 1.
C MA48D/DD solves block triangular system
C MA50C/CD solves non-blocked system
C MC71A/AD estimates norm of matrix
      INTRINSIC ABS,MAX

C Set streams for errors and warnings
      LP = ICNTL(1)
      MP = ICNTL(2)
C
C Simple data checks
      IF (N.LE.0 .OR. M.LE.0) THEN
         INFO(1) = -1
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9140) M,N
         GO TO 380
      END IF
      IF (JOB.GT.4 .OR. JOB.LT.1) THEN
         INFO(1) = -6
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9150) JOB
         GO TO 380
      END IF
      INFO(1) = 0

C Partition KEEP and A
      IPTRL = M + N
      IPTRU = IPTRL + N
      IPTRD = IPTRU + N
      IPTRO = IPTRD + N + 1
      NBLOCK = IPTRO + N - 1
      MBLOCK = NBLOCK + 1
      KBLOCK = MBLOCK + 1
      NB = KEEP(KBLOCK+3)

      NE = KEEP(IPTRO+N+1) - 1

      OMEGA(1) = ZERO
      OMEGA(2) = ZERO
      ERROR(3) = ZERO
C Initialize EPS
      EPS = EPSILON(EPS)
C CTAU ... 1000 eps (approx)
      CTAU = 1000.*EPS

      IQB(1) = 0

      DO 10 I = 1,7
         ICNTL5(I) = 0
   10 CONTINUE
C Set control for use of BLAS
      ICNTL5(5) = ICNTL(5)

      IF (MP.GT.0 .AND. ICNTL(3).GT.2) THEN
         WRITE (MP,
     +     '(/A/3(A,I8),A,I2/A,L2,A,I7/A,1P,E12.4/A,3I6,2(/A,I6))')
     +     ' Entering MA48C/CD with',' M =',M,'     N =',N,'     LA =',
     +     LA,'      JOB =',JOB,'   TRANS =',TRANS,
     +     '      No. of blocks =',
     +     NB,'   CNTL(5)    = ',CNTL(5),'   ICNTL(1:3) = ',
     +     (ICNTL(I),I=1,3),'   ICNTL(5)   = ',ICNTL(5),
     +     '   ICNTL(9)   = ',ICNTL(9)
         WRITE (MP,'(A)') ' Permuted matrix by blocks'
         KB = NB
         IF (ICNTL(3).EQ.3) THEN
            WRITE (MP,'(A)')
     +        ' Only first column of up to 10 blocks printed'
            KB = MIN(10,NB)
         END IF
         WRITE (MP,'(A)') ' Diagonal blocks'
         J1 = 1
         DO 30 JB = 1,KB
            J2 = J1 + KEEP(NBLOCK+3*JB) - 1
            IF (J1.LE.J2) WRITE (MP,'(A,I6)') ' Block',JB
            J3 = J2
            IF (ICNTL(3).EQ.3) J3 = MIN(J1,J2)
            DO 20 J = J1,J3
               WRITE (MP,'(A,I5,(T13,3(1P,E12.4,I5)))') ' Column',J,
     +           (A(I),IRN(I),I=KEEP(IPTRD+J),KEEP(IPTRD+J+1)-1)
   20       CONTINUE
            J1 = J2 + 1
   30    CONTINUE
         IF (KEEP(IPTRO+N+1).GT.KEEP(IPTRD+N+1)) THEN
            WRITE (MP,'(A)') ' Off-diagonal entries'
            J1 = 1
            DO 50 JB = 1,KB
               J2 = J1 + KEEP(NBLOCK+3*JB) - 1
               J3 = J2
               IF (ICNTL(3).EQ.3) J3 = MIN(J1,J2)
               DO 40 J = J1,J3
                  IF (KEEP(IPTRO+J+1).GT.KEEP(IPTRO+J)) WRITE (MP,
     +                '(A,I5,(T13,3(1P,E12.4,I5)))') ' Column',J,
     +                (A(I),IRN(I),I=KEEP(IPTRO+J),KEEP(IPTRO+J+1)-1)
   40          CONTINUE
               J1 = J2 + 1
   50       CONTINUE
         END IF
         WRITE (MP,'(A)') ' Factorized matrix by blocks'
         J1 = 1
         DO 90 JB = 1,KB
            J2 = J1 + KEEP(NBLOCK+3*JB) - 1
C Jump if triangular block
            IF (KEEP(MBLOCK+3*JB).LT.0) GO TO 80
            NC = J2 - J1 + 1
            WRITE (MP,'(A,I6)') ' Block',JB
            K = NE
            IF (JB.GT.1) K = KEEP(KBLOCK+3*JB) + NE - 1
            IF (KEEP(IPTRL+J1).GT.KEEP(IPTRU+J1)) WRITE (MP,
     +          '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column',J1,' of L',
     +          (A(K+I),IRN(K+I),I=KEEP(IPTRU+J1)+1,KEEP(IPTRL+J1))
            IF (ICNTL(3).EQ.3) GO TO 70
            DO 60 J = J1 + 1,J2
               IF (KEEP(IPTRU+J).GT.KEEP(IPTRL+J-1)) WRITE (MP,
     +             '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column',J,' of U',
     +              (A(K+I),IRN(K+I),I=KEEP(IPTRL+J-1)+1,KEEP(IPTRU+J))
               IF (KEEP(IPTRU+J).LT.KEEP(IPTRL+J)) WRITE (MP,
     +             '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column',J,' of L',
     +              (A(K+I),IRN(K+I),I=KEEP(IPTRU+J)+1,KEEP(IPTRL+J))
   60       CONTINUE
C Full blocks
   70       WRITE (MP,'(A)') ' Full block'
            WRITE (MP,'(A)') ' Row indices'
            NRF = IRN(K+2)
            K = K + KEEP(IPTRL+J2)
            IF (ICNTL(3).GT.3) THEN
               WRITE (MP,9000) (IRN(K+I),I=1,NRF)
               WRITE (MP,'(A)') ' Column pivoting information'
               WRITE (MP,9000) (IRN(K+I),I=NRF+1,
     +           NRF+NC-KEEP(MBLOCK+3*JB))
               WRITE (MP,'(A)') ' Reals by columns'
               WRITE (MP,9010) (A(K+I),I=1,
     +           (NRF)* (NC-KEEP(MBLOCK+3*JB)))
 9000          FORMAT (10I6)
 9010          FORMAT (1P,5D12.4)
            ELSE
               WRITE (MP,9000) (IRN(K+I),I=1,MIN(10,NRF))
               WRITE (MP,'(A)') ' Column pivoting information'
               WRITE (MP,9000) (IRN(K+I),I=NRF+1,
     +           NRF+MIN(10,NC-KEEP(MBLOCK+3*JB)))
               WRITE (MP,'(A)') ' Reals by columns'
               WRITE (MP,9010) (A(K+I),I=1,
     +           MIN(10, (NRF)* (NC-KEEP(MBLOCK+3*JB))))
            END IF
   80       J1 = J2 + 1
   90    CONTINUE
         IF (ICNTL(3).GT.3) THEN
            WRITE (MP,'(A)') ' Contents of KEEP array'
            WRITE (MP,9020) (KEEP(I),I=1,M)
            WRITE (MP,9030) (KEEP(M+I),I=1,N)
 9020       FORMAT (' Positions of original rows in the permuted matrix'
     +             / (10I6))
 9030       FORMAT (' Positions of columns of permuted matrix ','in or',
     +             'ig','inal matrix '/ (10I6))
            WRITE (MP,'(A)') ' Pointer information from KEEP'
            WRITE (MP,9040) (KEEP(IPTRD+J),J=1,N+1)
 9040       FORMAT (' IPTRD ='/ (8X,10I6))
            WRITE (MP,9050) (KEEP(IPTRO+J),J=1,N+1)
 9050       FORMAT (' IPTRO ='/ (8X,10I6))
            WRITE (MP,9060) (KEEP(IPTRL+J),J=1,N)
 9060       FORMAT (' IPTRL ='/ (8X,10I6))
            WRITE (MP,9070) (KEEP(IPTRU+J),J=1,N)
 9070       FORMAT (' IPTRU ='/ (8X,10I6))
            WRITE (MP,'(A)') ' Block structure information from KEEP'
            WRITE (MP,9080) (KEEP(NBLOCK+3*JB),JB=1,NB)
            WRITE (MP,9090) (KEEP(MBLOCK+3*JB),JB=1,NB)
            WRITE (MP,9100) (KEEP(KBLOCK+3*JB),JB=1,NB)
 9080       FORMAT (' NBLOCK (order blocks) ='/ (8X,10I6))
 9090       FORMAT (' MBLOCK (triangular flag and number packed rows) ='
     +             / (8X,10I6))
 9100       FORMAT (' KBLOCK (position of beginning of block) ='/
     +             (8X,10I6))
            IF (TRANS) THEN
               WRITE (MP,9110) (RHS(I),I=1,N)
 9110          FORMAT (' RHS =  ',1P,5D12.4/ (8X,5D12.4))
            ELSE
               WRITE (MP,9110) (RHS(I),I=1,M)
            END IF
         ELSE
            WRITE (MP,'(A)') ' Contents of KEEP array'
            WRITE (MP,9020) (KEEP(I),I=1,MIN(10,M))
            WRITE (MP,9030) (KEEP(M+I),I=1,MIN(10,N))
            WRITE (MP,'(A)') ' Pointer information from KEEP'
            WRITE (MP,9040) (KEEP(IPTRD+K),K=1,MIN(10,N+1))
            WRITE (MP,9050) (KEEP(IPTRO+K),K=1,MIN(10,N+1))
            WRITE (MP,9060) (KEEP(IPTRL+J),J=1,MIN(10,N))
            WRITE (MP,9070) (KEEP(IPTRU+J),J=1,MIN(10,N))
            WRITE (MP,'(A)') ' Block structure information from KEEP'
            WRITE (MP,9080) (KEEP(NBLOCK+3*JB),JB=1,KB)
            WRITE (MP,9090) (KEEP(MBLOCK+3*JB),JB=1,KB)
            WRITE (MP,9100) (KEEP(KBLOCK+3*JB),JB=1,KB)
            IF (TRANS) THEN
               WRITE (MP,9110) (RHS(I),I=1,MIN(10,N))
            ELSE
               WRITE (MP,9110) (RHS(I),I=1,MIN(10,M))
            END IF
         END IF
      END IF

C Apply global permutation to incoming right-hand side
C     and set NEQ and NVAR
      IF (TRANS) THEN
         NEQ = N
         NVAR = M
         DO 100 I = 1,NEQ
            W(I) = RHS(KEEP(M+I))
  100    CONTINUE
      ELSE
         NEQ = M
         NVAR = N
         DO 110 I = 1,NEQ
            W(KEEP(I)) = RHS(I)
  110    CONTINUE
      END IF
C
C
C Straight solution requested with no iterative refinement or error
C     estimate.
      IF (JOB.EQ.1) THEN
C
C Solve system using MA48D/DD or MA50C/CD.
         IF (NB.EQ.1 .AND. KEEP(MBLOCK+3).GE.0) THEN
            CALL MA50CD(M,N,ICNTL5,IQB,KEEP(MBLOCK+3),TRANS,LA-NE,
     +                  A(NE+1),IRN(NE+1),KEEP(IPTRL+1),KEEP(IPTRU+1),W,
     +                  X,W(NEQ+1),INFO)
         ELSE
            CALL MA48DD
     +           (N,NE,LA-NE,A(NE+1),A,IRN(NE+1),IRN,KEEP(IPTRD+1),
     +            KEEP(IPTRO+1),NB,KEEP(NBLOCK+3),KEEP(IPTRL+1),
     +            KEEP(IPTRU+1),W,X,TRANS,ICNTL5,W(NEQ+1))
         END IF
         GO TO 340
      END IF

C Prepare for iterative refinement.
C
C Set initial estimate of solution to zero
      DO 120 I = 1,NVAR
         X(I) = ZERO
  120 CONTINUE
C Save permuted right-hand side for residual calculation
      DO 130 I = 1,NEQ
         RHS(I) = W(I)
  130 CONTINUE

C
C Iterative refinement loop
C Initialize OM1 in case of problems with optimizing compiler
      OM1 = ZERO
      DO 260 K = 1,ICNTL(9)
C
C Solve system using MA48D/DD or MA50C/CD.
         IF (NB.EQ.1 .AND. KEEP(MBLOCK+3).GE.0) THEN
            CALL MA50CD(M,N,ICNTL5,IQB,KEEP(MBLOCK+3),TRANS,LA-NE,
     +                A(NE+1),IRN(NE+1),KEEP(IPTRL+1),KEEP(IPTRU+1),
     +                W,W(NEQ+1),W(M+N+1),INFO)
         ELSE
            CALL MA48DD
     +           (N,NE,LA-NE,A(NE+1),A,IRN(NE+1),IRN,KEEP(IPTRD+1),
     +            KEEP(IPTRO+1),NB,KEEP(NBLOCK+3),KEEP(IPTRL+1),
     +            KEEP(IPTRU+1),W,W(NEQ+1),TRANS,ICNTL5,W(M+N+1))
         END IF
C Update solution
         DO 140 I = 1,NVAR
            X(I) = X(I) + W(NEQ+I)
  140    CONTINUE
C
C Calculate residual using information in A,IRN
C
         DO 150 I = 1,NEQ
C Residual  .. b-Ax
            W(I) = RHS(I)
C |A||x|
            W(NEQ+I) = ZERO
C Sum |a  |, j=1,N (= ||A  ||        )
C       ij               i.  infinity
            W(2*NEQ+I) = ZERO
  150    CONTINUE
         IF (TRANS) THEN
            DO 180 J = 1,N
CDIR$ IVDEP
               DO 160 JJ = KEEP(IPTRD+J),KEEP(IPTRD+J+1) - 1
                  I = IRN(JJ)
                  W(J) = W(J) - A(JJ)*X(I)
                  W(NEQ+J) = W(NEQ+J) + ABS(A(JJ)*X(I))
                  W(2*NEQ+J) = W(2*NEQ+J) + ABS(A(JJ))
  160          CONTINUE
CDIR$ IVDEP
               DO 170 JJ = KEEP(IPTRO+J),KEEP(IPTRO+J+1) - 1
                  I = IRN(JJ)
                  W(J) = W(J) - A(JJ)*X(I)
                  W(NEQ+J) = W(NEQ+J) + ABS(A(JJ)*X(I))
                  W(2*NEQ+J) = W(2*NEQ+J) + ABS(A(JJ))
  170          CONTINUE
  180       CONTINUE
         ELSE
            DO 210 J = 1,N
CDIR$ IVDEP
               DO 190 JJ = KEEP(IPTRD+J),KEEP(IPTRD+J+1) - 1
                  I = IRN(JJ)
                  W(I) = W(I) - A(JJ)*X(J)
                  W(NEQ+I) = W(NEQ+I) + ABS(A(JJ)*X(J))
                  W(2*NEQ+I) = W(2*NEQ+I) + ABS(A(JJ))
  190          CONTINUE
CDIR$ IVDEP
               DO 200 JJ = KEEP(IPTRO+J),KEEP(IPTRO+J+1) - 1
                  I = IRN(JJ)
                  W(I) = W(I) - A(JJ)*X(J)
                  W(NEQ+I) = W(NEQ+I) + ABS(A(JJ)*X(J))
                  W(2*NEQ+I) = W(2*NEQ+I) + ABS(A(JJ))
  200          CONTINUE
  210       CONTINUE
         END IF
C Calculate max-norm of solution
         DXMAX = ZERO
         DO 220 I = 1,NVAR
            DXMAX = MAX(DXMAX,ABS(X(I)))
  220    CONTINUE
C Calculate also omega(1) and omega(2)
C tau is (||A  ||         ||x||   + |b| )*n*1000*epsilon
C            i.  infinity      max     i
         OMEGA(1) = ZERO
         OMEGA(2) = ZERO
         DO 230 I = 1,NEQ
            TAU = (W(2*NEQ+I)*DXMAX+ABS(RHS(I)))*NVAR*CTAU
            IF ((W(NEQ+I)+ABS(RHS(I))).GT.TAU) THEN
C |Ax-b| /(|A||x| + |b|)
C       i               i
               OMEGA(1) = MAX(OMEGA(1),ABS(W(I))/
     +                    (W(NEQ+I)+ABS(RHS(I))))
               IW(I) = 1
            ELSE
C TAU will be zero if all zero row in A, for example
               IF (TAU.GT.ZERO) THEN
C |Ax-b| /(|A||x| + ||A  ||        ||x||   )
C       i        i     i.  infinity     max
                  OMEGA(2) = MAX(OMEGA(2),ABS(W(I))/
     +                       (W(NEQ+I)+W(2*NEQ+I)*DXMAX))
               END IF
               IW(I) = 2
            END IF
  230    CONTINUE
C
C Exit if iterative refinement not being performed
         IF (JOB.EQ.2) GO TO 340
C
C  Stop the calculations if the backward error is small
C
         OM2 = OMEGA(1) + OMEGA(2)
C Jump if converged
C Statement changed because IBM SP held quantities in registers
C        IF ((OM2+ONE).LE.ONE) GO TO 270
         IF (OM2.LE.EPS) GO TO 270
C
C  Check the convergence.
C
         IF (K.GT.1 .AND. OM2.GT.OM1*CNTL(5)) THEN
C  Stop if insufficient decrease in omega.
            IF (OM2.GT.OM1) THEN
C Previous estimate was better ... reinstate it.
               OMEGA(1) = OLDOMG(1)
               OMEGA(2) = OLDOMG(2)
               DO 240 I = 1,NVAR
                  X(I) = W(3*NEQ+I)
  240          CONTINUE
            END IF
            GO TO 270
         END IF
C Hold current estimate in case needed later
         DO 250 I = 1,NVAR
            W(3*NEQ+I) = X(I)
  250    CONTINUE
         OLDOMG(1) = OMEGA(1)
         OLDOMG(2) = OMEGA(2)
         OM1 = OM2
  260 CONTINUE
C End of iterative refinement loop.
      INFO(1) = -8
      IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9170) INFO(1),ICNTL(9)
      GO TO 340

  270 IF (JOB.LE.3) GO TO 340
      IF (M.NE.N) GO TO 340
C
C Calculate condition numbers and estimate of the error.
C
C  Condition numbers obtained through use of norm estimation
C     routine MC71A/AD.
C
C  Initializations
C
      LCOND(1) = .FALSE.
      LCOND(2) = .FALSE.
      DO 280 I = 1,NEQ
         IF (IW(I).EQ.1) THEN
            W(I) = W(NEQ+I) + ABS(RHS(I))
C |A||x| + |b|
            W(NEQ+I) = ZERO
            LCOND(1) = .TRUE.
         ELSE
C |A||x| + ||A  ||        ||x||
C             i.  infinity     max

            W(NEQ+I) = W(NEQ+I) + W(2*NEQ+I)*DXMAX
            W(I) = ZERO
            LCOND(2) = .TRUE.
         END IF
  280 CONTINUE
C
C  Compute the estimate of COND
C
      KASE = 0
      DO 330 K = 1,2
         IF (LCOND(K)) THEN
C MC71A/AD has its own built in limit to the number of iterations
C    allowed. It is this limit that will be used to terminate the
C    following loop.
            DO 310 KK = 1,40
C MC71A/AD calculates norm of matrix
C We are calculating the infinity norm of INV(A).W
               CALL MC71AD(N,KASE,W(3*NEQ+1),COND(K),RHS,IW,KEEP71)
C
C  KASE = 0........ Computation completed
C  KASE = 1........ W * INV(TRANSPOSE(A)) * Y
C  KASE = 2........ INV(A) * W * Y
C                   W is W/W(NEQ+1) .. Y is W(3*NEQ+1)
C
               IF (KASE.EQ.0) GO TO 320
               IF (KASE.EQ.1) THEN
C Solve system using MA48D/DD or MA50C/CD.
                  IF (NB.EQ.1 .AND. KEEP(MBLOCK+3).GE.0) THEN
                     CALL MA50CD(M,N,ICNTL5,IQB,KEEP(MBLOCK+3),
     +                           .NOT.TRANS,LA-NE,A(NE+1),IRN(NE+1),
     +                           KEEP(IPTRL+1),KEEP(IPTRU+1),W(3*NEQ+1),
     +                           W(2*NEQ+1),RHS,INFO)
                  ELSE
                     CALL MA48DD(N,NE,LA-NE,A(NE+1),A,IRN(NE+1),IRN,
     +                           KEEP(IPTRD+1),KEEP(IPTRO+1),NB,
     +                           KEEP(NBLOCK+3),KEEP(IPTRL+1),
     +                           KEEP(IPTRU+1),W(3*NEQ+1),W(2*NEQ+1),
     +                           .NOT.TRANS,ICNTL5,RHS)
                  END IF

                  DO 290 I = 1,M
                     W(3*NEQ+I) = W((K-1)*NEQ+I)*W(2*NEQ+I)
  290             CONTINUE
               END IF
               IF (KASE.EQ.2) THEN
                  DO 300 I = 1,N
                     W(2*NEQ+I) = W((K-1)*NEQ+I)*W(3*NEQ+I)
  300             CONTINUE
C Solve system using MA48D/DD or MA50C/CD.
                  IF (NB.EQ.1 .AND. KEEP(MBLOCK+3).GE.0) THEN
                     CALL MA50CD(M,N,ICNTL5,IQB,KEEP(MBLOCK+3),TRANS,
     +                           LA-NE,A(NE+1),IRN(NE+1),KEEP(IPTRL+1),
     +                           KEEP(IPTRU+1),W(2*NEQ+1),W(3*NEQ+1),
     +                           RHS,INFO)
                  ELSE
                     CALL MA48DD(N,NE,LA-NE,A(NE+1),A,IRN(NE+1),IRN,
     +                           KEEP(IPTRD+1),KEEP(IPTRO+1),NB,
     +                           KEEP(NBLOCK+3),KEEP(IPTRL+1),
     +                           KEEP(IPTRU+1),W(2*NEQ+1),W(3*NEQ+1),
     +                           TRANS,ICNTL5,RHS)
                  END IF
               END IF
  310       CONTINUE
            INFO(1) = -9
            IF (LP.NE.0 .AND. ICNTL(3).GE.1) WRITE (LP,9160)
            GO TO 340
  320       IF (DXMAX.GT.ZERO) COND(K) = COND(K)/DXMAX
            ERROR(3) = ERROR(3) + OMEGA(K)*COND(K)
         END IF
  330 CONTINUE

C Permute solution vector
  340 DO 350 I = 1,NVAR
         W(I) = X(I)
  350 CONTINUE
      IF (.NOT.TRANS) THEN
         DO 360 I = 1,NVAR
            X(KEEP(M+I)) = W(I)
  360    CONTINUE
      ELSE
         DO 370 I = 1,NVAR
            X(I) = W(KEEP(I))
  370    CONTINUE
      END IF
      IF (JOB.GE.2) THEN
         ERROR(1) = OMEGA(1)
         ERROR(2) = OMEGA(2)
      END IF

      IF (MP.GT.0 .AND. ICNTL(3).GT.2) THEN
         WRITE (MP,'(/A,I6)') ' Leaving MA48C/CD with INFO(1) =',INFO(1)
         IF (JOB.GT.1) THEN
            K = 2
            IF (JOB.EQ.4 .AND. INFO(1).NE.-9) K = 3
            WRITE (MP,9120) (ERROR(I),I=1,K)
 9120       FORMAT (' ERROR =',1P,3D12.4)
         END IF
         IF (ICNTL(3).GT.3) THEN
            WRITE (MP,9130) (X(I),I=1,NVAR)
 9130       FORMAT (' X =    ',1P,5D12.4:/ (8X,5D12.4))
         ELSE
            WRITE (MP,9130) (X(I),I=1,MIN(10,NVAR))
         END IF
      END IF
  380 RETURN

 9140 FORMAT (' Error return from MA48C/CD because M =',I10,' and N =',
     +       I10)
 9150 FORMAT (' Error return from MA48C/CD because ','JOB = ',I10)
 9160 FORMAT (' Error return from MA48C/CD because of ','error in MC71',
     +       'A/AD'/' ERROR(3) not calculated')
 9170 FORMAT (' Error return from MA48C/CD because of ','nonconvergenc',
     +       'e of iterative refinement'/' Error INFO(1) = ',I2,'  wit',
     +       'h ICNTL','(9) = ',I10)
      END


      SUBROUTINE MA48DD(N,NE,LA,A,AA,IRN,IRNA,IPTRD,IPTRO,NB,IBLOCK,
     +                  IPTRL,IPTRU,RHS,X,TRANS,ICNTL5,W)
C Solve linear system, using data provided by MA48B/BD.
C
C     .. Arguments ..
      INTEGER N,NE,LA
      DOUBLE PRECISION A(LA),AA(NE)
      INTEGER IRN(LA),IRNA(NE),IPTRD(N+1),IPTRO(N+1),NB
      INTEGER IBLOCK(3,NB),IPTRL(N),IPTRU(N)
      DOUBLE PRECISION RHS(N),X(N)
      LOGICAL TRANS
      INTEGER ICNTL5(20)
      DOUBLE PRECISION W(N)
C
C N must be set by the user to the number of columns in the matrix.
C      It is not altered by the subroutine.
C NE must be set by the user to the size of arrays AA and IRNA.
C      It is not altered by the subroutine.
C LA must be set by the user to the size of arrays A and IRN.
C      It is not altered by the subroutine.
C A    must be set left unchanged from the last call to MA48B/BD.
C      It holds the factorization of the block diagonal part
C      (excluding triangular blocks).
C      It is not altered by the subroutine.
C AA   must be set left unchanged from the last call to MA48B/BD.
C      It holds the original matrix in permuted form.
C      It is not altered by the subroutine.
C IRN  must be set left unchanged from the last call to MA48B/BD.
C      It holds the row indices of the factorization in A.
C      It is not altered by the subroutine.
C IRNA must be set left unchanged from the last call to MA48B/BD.
C      It holds the row indices of the original matrix in permuted form.
C      It is not altered by the subroutine.
C IPTRD must be as on return from MA48B/BD. IPTRD(j) holds the position
C      in IRNA of the start of the block diagonal part of column j.
C      IPTRD(n+1) holds the position immediately after the end of
C      column n.
C      It is not altered by the subroutine.
C IPTRO must be as on return from MA48A/AD. IPTRO(j) holds the position
C      in IRNA of the start of the off block diagonal part of column j.
C      IPTRO(n+1) holds the position immediately after the end of
C      column n.
C      It is not altered by the subroutine.
C NB must be set by the user to the second dimension of IBLOCK.
C      It is not altered by the subroutine.
C IBLOCK must be as on return from MA48B/BD.  IBLOCK(1,k) holds the
C      order of diagonal block k, k=1,2,...  IBLOCK(2,k) holds the
C      number of rows held in packed storage when
C      processing block k, or a negative value for triangular blocks.
C      IBLOCK(3,1) holds the number of blocks and IBLOCK(3,k) holds the
C      position in the factorization of the start of block k, k =
C      2,3,...  IBLOCK is not altered by the subroutine.
C IPTRL must be unchanged since the last call to MA48B/BD.
C      It holds pointers to the
C      columns of the lower-triangular part of the factorization.
C      It is not altered by the subroutine.
C IPTRU must be unchanged since the last call to MA48B/BD.
C      It holds pointers to the
C      columns of the upper-triangular part of the factorization.
C      It is not altered by the subroutine.
C RHS  must be set by the user to contain the right-hand side of
C      the equations to be solved.
C      It is used as workspace.
C X    need not be set by the user.  On exit, it will contain the
C      solution vector.
C TRANS must be set by the user to indicate whether coefficient matrix
C      is A (TRANS=.FALSE.) or A transpose (TRANS=.TRUE.).
C      It is not altered by the subroutine.
C ICNTL5 passed to MA50C/CD to correspond to dummy argument ICNTL.
C W    is a workarray.
C
C
C     .. Local variables ..
      INTEGER I,IFLAG(15),IQB(1),J,JB,JJ,J1,K1,K2,NC,NUMB
C I     row index
C IFLAG error flag from MA50CD .. will never be set
C IQB   is passed to MA50C/CD to indicate that no column permutation
C       is used.
C J     column index
C JB    block index.
C JJ    running index for column.
C J1    position of beginning of block
C K1    index of first column in block
C K2    index of last column in block
C NC    number of columns in block
C NUMB  number of diagonal blocks
C
C     .. Externals ..
      EXTERNAL MA50CD
C
      IQB(1) = 0
      NUMB = IBLOCK(3,1)
      IF (.NOT.TRANS) THEN
C Solve system using block structure information and calls to MA50C
C     for each diagonal block.
C System is block upper triangular.
         K1 = N + 1
         DO 50 JB = NUMB,1,-1
            NC = IBLOCK(1,JB)
            K2 = K1 - 1
            K1 = K1 - NC
            IF (IBLOCK(2,JB).LT.0) THEN
C Process triangular block
               DO 20 J = K2,K1,-1
                  X(J) = RHS(J)/AA(IPTRD(J))
CDIR$ IVDEP
                  DO 10 JJ = IPTRD(J) + 1,IPTRD(J+1) - 1
                     I = IRNA(JJ)
                     RHS(I) = RHS(I) - AA(JJ)*X(J)
   10             CONTINUE
   20          CONTINUE
            ELSE
               J1 = 1
               IF (JB.GT.1) J1 = IBLOCK(3,JB)
               CALL MA50CD(NC,NC,ICNTL5,IQB,IBLOCK(2,JB),TRANS,LA+1-J1,
     +                     A(J1),IRN(J1),IPTRL(K1),IPTRU(K1),RHS(K1),
     +                     X(K1),W,IFLAG)
            END IF
            IF (JB.EQ.1) GO TO 50
C Substitution using off-diagonal block
            DO 40 J = K1,K2
CDIR$ IVDEP
               DO 30 JJ = IPTRO(J),IPTRO(J+1) - 1
                  I = IRNA(JJ)
                  RHS(I) = RHS(I) - AA(JJ)*X(J)
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      ELSE
C Solve system using block structure information and calls to MA50C
C     for each diagonal block.
C System is block lower triangular.
         K2 = 0
         DO 100 JB = 1,NUMB
            NC = IBLOCK(1,JB)
            K1 = K2 + 1
            K2 = K2 + NC
            IF (JB.GT.1) THEN
C Substitution using off-diagonal block
               DO 70 J = K1,K2
                  DO 60 JJ = IPTRO(J),IPTRO(J+1) - 1
                     I = IRNA(JJ)
                     RHS(J) = RHS(J) - AA(JJ)*X(I)
   60             CONTINUE
   70          CONTINUE
            END IF
            IF (IBLOCK(2,JB).LT.0) THEN
C Process triangular block
               DO 90 J = K1,K2
                  DO 80 JJ = IPTRD(J) + 1,IPTRD(J+1) - 1
                     I = IRNA(JJ)
                     RHS(J) = RHS(J) - AA(JJ)*X(I)
   80             CONTINUE
                  X(J) = RHS(J)/AA(IPTRD(J))
   90          CONTINUE
            ELSE
               J1 = 1
               IF (JB.GT.1) J1 = IBLOCK(3,JB)
               CALL MA50CD(NC,NC,ICNTL5,IQB,IBLOCK(2,JB),TRANS,LA+1-J1,
     +                     A(J1),IRN(J1),IPTRL(K1),IPTRU(K1),RHS(K1),
     +                     X(K1),W,IFLAG)
            END IF
  100    CONTINUE
      END IF

      RETURN
      END

      SUBROUTINE MA48ID(CNTL,ICNTL)
C Set default values for the control arrays.

      DOUBLE PRECISION CNTL(10)
      INTEGER ICNTL(20)

C CNTL  is a real array of length 10.
C     CNTL(1)  If this is set to a value less than or equal to one, full
C       matrix processing will be used by MA50A/AD  when the density of
C       the reduced matrix reaches CNTL(1).
C     CNTL(2) determines the balance used by MA50A/AD and MA50B/BD
C       between pivoting for sparsity and for stability, values near
C       zero emphasizing sparsity and values near one emphasizing
C       stability.
C     CNTL(3) If this is set to a positive value, any entry whose
C       modulus is less than CNTL(3) will be dropped from the factors
C       calculated by MA50A/AD. The factorization will then require
C       less storage but will be inaccurate.
C     CNTL(4)  If this is set to a positive value, any entry whose
C       modulus is less than CNTL(4) will be regarded as zero from
C       the point of view of rank.
C     CNTL(5) is used in convergence test for termination of iterative
C       refinement.  Iteration stops if successive omegas do not
C       decrease by at least this factor.
C     CNTL(6:10) are not used.
C ICNTL is an integer array of length 20.
C     ICNTL(1)  must be set to the stream number for error messages.
C       A value less than 1 suppresses output.
C     ICNTL(2) must be set to the stream number for diagnostic output.
C       A value less than 1 suppresses output.
C     ICNTL(3) must be set to control the amount of output:
C       0 None.
C       1 Error messages only.
C       2 Error and warning messages.
C       3 As 2, plus scalar parameters and a few entries of array
C         parameters on entry and exit.
C       4 As 2, plus all parameters on entry and exit.
C     ICNTL(4)  If set to a positive value, the pivot search by MA50A/AD
C       is limited to ICNTL(4) columns. This may result in different
C       fill-in and execution time but could give faster execution.
C     ICNTL(5) The block size to be used for full-matrix processing.
C     ICNTL(6) is the minimum size for a block of the block triangular
C       form.
C     ICNTL(7) If not equal to 0, abort when structurally rank deficient
C       matrix found.
C     ICNTL(8) If set to a value other than zero and JOB = 1 or 3 on
C       entry to MA48A/AD, columns with IW flagged 0 are placed at
C       the end of their respective blocks in the first factorization
C       and the remaining columns are assumed unchanged in subsequent
C       factorizations. If set to a value other than zero and JOB = 2,
C       columns before the first 0 entry of IW are assumed unchanged in
C       subsequent factorizations.  On entry to MA48B/BD, the number
C       of columns that can change in each block is held in array KEEP.
C     ICNTL(9) is the limit on the number of iterative refinements
C        allowed by MA48C/CD
C     ICNTL(10) has default value 0. If set to 1, there is an immediate
C       return from MA48B/BD if LA is too small, without continuing the
C       decomposition to compute the size necessary.
C     ICNTL(11) has default value 0. If set to 1 on a JOB=2 call to
C       MA48B/BD and the entries in one of the blocks on the diagonal
C       are unsuitable for the pivot sequence chosen on the previous
C       call, the block is refactorized as on a JOB=1 call.
C     ICNTL(12:20) are not used.

      INTEGER I

      DO 10 I = 3, 10
         CNTL(I) = 0.0D0
   10 CONTINUE
      CNTL(1) = 0.5D0
      CNTL(2) = 0.1D0
      CNTL(5) = 0.5D0
      ICNTL(1) = 6
      ICNTL(2) = 6
      ICNTL(3) = 2
      ICNTL(4) = 3
      ICNTL(5) = 32
      ICNTL(6) = 1
      ICNTL(7) = 1
      ICNTL(8) = 0
      ICNTL(9) = 10
      DO 20 I = 10, 20
         ICNTL(I) = 0
   20 CONTINUE

      END

