C COPYRIGHT (c) 1993 Council for the Central Laboratory
*                    of the Research Councils
C Original date 20 May 1993
C 24 September 1993. Some IVDEP comments added for speed on the Cray,
C     some minor bugs fixed, defaults changed to BLAS 3 with block size
C     32 and no amalgamation of blocks of the block triangular form.
C 18 October 1993. Minor bug in MA48BD when ICNTL(3)=3 corrected.
C 10 December 1993. Minor bugs in printing corrected.
C 14 June 1994. Minor bug rank calculation corrected.
C 14 June 1994. Minor bugs in printing corrected.
C 12/12/94 Calls of MC13D and MC21A changed to MC13DD and MC21AD
C 4/10/95. Redundant variables ONE removed.
C 1/11/95  IFLAG in MA48DD made of length 7.
C 1/11/95  Temporary variable NP introduced to avoid PFORT warnings.
C 14/2/96  Third argument in calls to MA48DD changed to LA-NE.
C 13/8/96  Bug corrected in print statements at end of MA48BD and
C          start of MA48CD.
C 07/09/96 Resetting of IPTRD and IRN after call to MA50BD moved to
C          immediately after the call so that they will be reset
C          correctly in the event of an error return from MA50BD.
C 12/02/01 Corrections made to the copying of the permutations found
C          by MA50A/AD in rectangular case.
C          Commas added after 1P edit descriptors.
C          Length of ICNTL5 corrected in MA48D/DD
C 09/08/01 use of MC41 changed to MC71
C 6/12/02. The test for convergence of iterative refinement changed to
C avoid any problem with comparisons of numbers held in registers.

C 12th July 2004 Version 1.0.0. Version numbering added.
C 21st February 2005 Version 1.1.0. FD05 dependence changed to FD15.
C 24th May 2006 Version 1.2.0. INFO(4) from MA48A/AD now cannot exceed
C           NE*2.
C 24th May 2006 Version 2.0.0. Sizes of CNTL, ICNTL, INFO, RINFO
C     increased and ICNTL(10) and ICNTL(11) added.
C 13th March 2007 Version 2.1.0. Loops split just ahead of call of
C      MA50B/BD to correct bug.

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
      EXTERNAL FD15AD,MA48DD,MA50CD,MC71AD
      DOUBLE PRECISION EPS,FD15AD
C EPS is the largest real such that 1+EPS is equal to 1.
C FD15A/AD provides machine constants
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
      EPS = FD15AD('E')
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

C COPYRIGHT (c) 1993 Council for the Central Laboratory
*                    of the Research Councils
C Original date 20 May 1993
C 24 September 1993 Some IVDEP comments added for speed on the Cray,
C     some minor bugs fixed, default changed to BLAS 3 with block size
C     32.
C 6 December 1993. Minor bug fixed re threshold test for pivots.
C 4/10/95. IQ in MA50BD made assumed size
C 1/11/95. IQ in MA50CD made assumed size
C 1/11/95. DTRSV not called for zero-sized array.
C 14/2/96. NP initialized to 0.
C 13/11/97 INFO(4) and INFO(6) in MA50AD made to reflect the situation
C          at the point of failure in the case of insufficient storage.
C 17/3/98  In MA50AD, copy a row forward if there is space at its front,
C          rather than put new entry at front. Makes the result
C          repeatable as LA is altered.

C 12th July 2004 Version 1.0.0. Version numbering added.
C 29 November 2006 Version 2.0.0. Sizes of CNTL, ICNTL, INFO, RINFO
C          increased and ICNTL(8) added.

      SUBROUTINE MA50AD(M,N,NE,LA,A,IRN,JCN,IQ,CNTL,ICNTL,IP,NP,JFIRST,
     +                  LENR,LASTR,NEXTR,IW,IFIRST,LENC,LASTC,NEXTC,
     +                  INFO,RINFO)

C MA50A/AD chooses a pivot sequence using a Markowitz criterion with
C     threshold pivoting.

C If  the user requires a more convenient data interface then the MA48
C     package should be used. The MA48 subroutines call the MA50
C     subroutines after checking the user's input data and optionally
C     permute the matrix to block triangular form.

      INTEGER M,N,NE,LA
      DOUBLE PRECISION A(LA)
      DOUBLE PRECISION CNTL(10)
      INTEGER IRN(LA),JCN(LA),IQ(N)
      INTEGER ICNTL(20),IP(M),NP,JFIRST(M),LENR(M),LASTR(M),NEXTR(M),
     +        IW(M),IFIRST(N),LENC(N),LASTC(N),NEXTC(N),INFO(15)
      DOUBLE PRECISION RINFO(10)

C M is an integer variable that must be set to the number of rows.
C      It is not altered by the subroutine.
C N is an integer variable that must be set to the number of columns.
C      It is not altered by the subroutine.
C NE is an integer variable that must be set to the number of entries
C      in the input matrix. It is not altered by the subroutine.
C LA is an integer variable that must be set to the size of A, IRN, and
C      JCN. It is not altered by the subroutine.
C A is an array that holds the input matrix on entry and is used as
C      workspace.
C IRN  is an integer array.  Entries 1 to NE must be set to the
C      row indices of the corresponding entries in A.  IRN is used
C      as workspace and holds the row indices of the reduced matrix.
C JCN  is an integer array that need not be set by the user. It is
C      used to hold the column indices of entries in the reduced
C      matrix.
C IQ is an integer array of length N. On entry, it holds pointers
C      to column starts. During execution, IQ(j) holds the position of
C      the start of column j of the reduced matrix or -IQ(j) holds the
C      column index in the permuted matrix of column j. On exit, IQ(j)
C      holds the index of the column that is in position j of the
C      permuted matrix.
C CNTL must be set by the user as follows and is not altered.
C     CNTL(1)  Full matrix processing will be used if the density of
C       the reduced matrix is MIN(CNTL(1),1.0) or more.
C     CNTL(2) determines the balance between pivoting for sparsity and
C       for stability, values near zero emphasizing sparsity and values
C       near one emphasizing stability. Each pivot must have absolute
C       value at least CNTL(2) times the greatest absolute value in the
C       same column of the reduced matrix.
C     CNTL(3) If this is set to a positive value, any entry of the
C       reduced matrix whose modulus is less than CNTL(3) will be
C       dropped.
C     CNTL(4)  Any entry of the reduced matrix whose modulus is less
C       than or equal to CNTL(4) will be regarded as zero from the
C        point of view of rank.
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
C       4 As 3, plus all parameters on entry and exit.
C     ICNTL(4) If set to a positive value, the pivot search is limited
C       to ICNTL(4) columns (Zlatev strategy). This may result in
C       different fill-in and execution time. If ICNTL(4) is positive,
C       the workspace arrays LASTR and NEXTR are not referenced.
C     ICNTL(5) The block size to be used for full-matrix processing.
C     ICNTL(6) The last ICNTL(6) columns of A must be the last
C       ICNTL(6) columns of the permuted matrix. A value outside the
C       range 1 to N-1 is treated as zero.
C     ICNTL(7) If given the value 1, pivots are limited to
C       the main diagonal, which may lead to a premature switch to full
C       processing if no suitable diagonal entries are available.
C       If given the value 2, IFIRST must be set so that IFIRST(i) is
C       the column in position i of the permuted matrix and IP must
C       be set so that IP(i) < IP(j) if row i is recommended to
C       precede row j in the pivot sequence.
C IP is an integer array of length M that need not be set on entry
C      unless ICNTL(7)=2 (see ICNTL(7) for details of this case).
C      During execution, IP(i) holds the position of the start of row i
C      of the reduced matrix or -IP(i) holds the row index in the
C      permuted matrix of row i. Before exit, IP(i) is made positive.
C NP is an integer variable. It need not be set on entry. On exit,
C     it will be set to the number of columns to be processed in
C     packed storage.
C JFIRST is an integer workarray of length M. JFIRST(i) is the
C      first column of the reduced matrix to have i entries or is
C      zero if no column has i entries.
C LENR is an integer workarray of length M that is used to hold the
C      numbers of entries in the rows of the reduced matrix.
C LASTR is an integer workarray of length M, used only if ICNTL(4) = 0.
C      For rows in the reduced matrix, LASTR(i) indicates the previous
C      row to i with the same number of entries. LASTR(i) is zero if
C      no such row exists.
C NEXTR is an integer workarray of length M, used only if ICNTL(4) = 0
C      or ICNTL(7)=2. If ICNTL(4)=0, for rows in the reduced matrix,
C      NEXTR(i) indicates the next row to i with the same number of
C      entries; and if row i is the last in the chain, NEXTR is
C      equal to zero. If ICNTL(7)=2, NEXTR is a copy of the value of
C      IP on entry.
C IW is an integer array of length M used as workspace and is used to
C     assist the detection of duplicate entries and the sparse SAXPY
C     operations. It is reset to zero each time round the main loop.
C IFIRST is an integer array of length N, used only if ICNTL(4) = 0
C      or ICNTL(7)=2. If ICNTL(4) = 0, it is a workarray; IFIRST(i)
C      points to the first row of the reduced matrix to have i entries
C      or is zero if no row has i entries. If ICNTL(7)=2, IFIRST
C      must be set on entry (see ICNTL(7) for details of this case).
C LENC is an integer workarray of length N that is used to hold
C      the numbers of entries in the columns of the reduced matrix.
C LASTC is an integer workarray of length N.  For columns in the reduced
C      matrix, LASTC(j) indicates the previous column to j with the same
C      number of entries.  If column j is the first in the chain,
C      LASTC(j) is equal to zero.
C NEXTC is an integer workarray of length N.  For columns in the reduced
C      matrix, NEXTC(j) indicates the next column to j with the same
C      number of entries.  If column j is the last column in the chain,
C      NEXTC(j) is zero.
C INFO need not be set on entry. On exit, it holds the following:
C    INFO(1):
C       0  Successful entry.
C      -1  M < 1 or N < 1.
C      -2  NE < 1.
C      -3  Insufficient space.
C      -4  Duplicated entries.
C      -5  Faulty column permutation in IFIRST when ICNTL(7)=2.
C      -6  ICNTL(4) not equal to 1 when ICNTL(7)=2.
C      +1  Rank deficient.
C      +2  Premature switch to full processing because of failure to
C          find a stable diagonal pivot (ICNTL(7)>=1 case only).
C      +3  Both of these warnings.
C    INFO(2) Number of compresses of the arrays.
C    INFO(3) Minimum LA recommended to analyse matrix.
C    INFO(4) Minimum LFACT required to factorize matrix.
C    INFO(5) Upper bound on the rank of the matrix.
C    INFO(6) Number of entries dropped from the data structure.
C    INFO(7) Number of rows processed in full storage.
C RINFO need not be set on entry. On exit, RINFO(1) holds the number of
C    floating-point operations needed for the factorization.

      INTEGER IDAMAX
      EXTERNAL IDAMAX
      EXTERNAL MA50DD
      INTRINSIC ABS,MAX,MIN

      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0D0,ONE=1.0D0)

      DOUBLE PRECISION ALEN,AMULT,ANEW,ASW,AU,COST,CPIV
      INTEGER DISPC,DISPR,EYE,I,IDROP,IDUMMY,IEND,IFILL,IFIR,II,IJ,
     +        IJPOS,IOP,IPIV,IPOS,ISRCH,I1,I2,J,JBEG,JEND,JJ,JLAST,
     +        JMORE,JNEW,JPIV,JPOS,J1,J2,L,LC,LEN,LENPIV,LP,LR
      DOUBLE PRECISION MAXENT
      INTEGER MINC,MORD,MP,MSRCH,NC,NDROP,NEFACT,NEPR,NERED,NE1,NORD,
     +        NORD1,NR,NULLC,NULLI,NULLJ,NULLR,PIVBEG,PIVCOL,PIVEND,
     +        PIVOT
      DOUBLE PRECISION PIVR,PIVRAT,U

C ALEN Real(LEN-1).
C AMULT Temporary variable used to store current multiplier.
C ANEW Temporary variable used to store value of fill-in.
C ASW Temporary variable used when swopping two real quantities.
C AU Temporary variable used in threshold test.
C COST Markowitz cost of current potential pivot.
C CPIV Markowitz cost of best pivot so far found.
C DISPC is the first free location in the column file.
C DISPR is the first free location in the row file.
C EYE Running relative position when processing pivot row.
C I Temporary variable holding row number. Also used as index in DO
C     loops used in initialization of arrays.
C IDROP Temporary variable used to accumulate number of entries dropped.
C IDUMMY DO index not referenced in the loop.
C IEND Position of end of pivot row.
C IFILL is the fill-in to the non-pivot column.
C IFIR Temporary variable holding first entry in chain.
C II Running position for current column.
C IJ Temporary variable holding row/column index.
C IJPOS Position of current pivot in A/IRN.
C IOP holds a running count of the number of rows with entries in both
C     the pivot and the non-pivot column.
C IPIV Row of the pivot.
C IPOS Temporary variable holding position in column file.
C ISRCH Temporary variable holding number of columns searched for pivot.
C I1 Position of the start of the current column.
C I2 Position of the end of the current column.
C J Temporary variable holding column number.
C JBEG Position of beginning of non-pivot column.
C JEND Position of end of non-pivot column.
C JJ Running position for current row.
C JLAST Last column acceptable as pivot.
C JMORE Temporary variable holding number of locations still needed
C     for fill-in in non-pivot column.
C JNEW Position of end of changed non-pivot column.
C JPIV Column of the pivot.
C JPOS Temporary variable holding position in row file.
C J1 Position of the start of the current row.
C J2 Position of the end of the current row.
C L Loop index.
C LC Temporary variable holding previous column in sequence.
C LEN Length of column or row.
C LENPIV Length of pivot column.
C LP Unit for error messages.
C LR Temporary variable holding previous row in sequence.
C MAXENT Temporary variable used to hold value of largest entry in
C    column.
C MINC Minimum number of entries of any row or column of the reduced
C     matrix, or in any column if ICNTL(4) > 0.
C MORD Number of rows ordered, excluding null rows.
C MP Unit for diagnostic messages.
C MSRCH Number of columns to be searched.
C NC Temporary variable holding next column in sequence.
C NDROP Number of entries dropped because of being in a column all of
C   whose entries are smaller than the pivot threshold.
C NEFACT Number of entries in factors.
C NEPR Number of entries in pivot row, excluding the pivot.
C NERED Number of entries in reduced matrix.
C NE1 Temporary variable used to hold number of entries in row/column
C     and to hold temporarily value of MINC.
C NORD Number of columns ordered, excluding null columns beyond JLAST.
C NORD1 Value of NORD at start of step.
C NR Temporary variable holding next row in sequence.
C NULLC Number of structurally zero columns found before any entries
C     dropped for being smaller than CNTL(3).
C NULLR Number of structurally zero rows found before any entries
C     dropped for being smaller than CNTL(3).
C NULLI Number of zero rows found.
C NULLJ Number of zero columns found beyond column JLAST.
C PIVBEG Position of beginning of pivot column.
C PIVCOL Temporary variable holding position in pivot column.
C PIVEND Position of end of pivot column.
C PIVOT Current step in Gaussian elimination.
C PIVR ratio of current pivot candidate to largest in its column.
C PIVRAT ratio of best pivot candidate to largest in its column.
C U Used to hold local copy of CNTL(2), changed if necessary so that it
C    is in range.

      LP = ICNTL(1)
      IF (ICNTL(3).LE.0) LP = 0
      MP = ICNTL(2)
      IF (ICNTL(3).LE.1) MP = 0
      INFO(1) = 0
      INFO(2) = 0
      INFO(3) = NE
      INFO(4) = NE
      INFO(5) = 0
      INFO(6) = 0
      INFO(7) = 0
      RINFO(1) = ZERO
      NP = 0

C Make some simple checks
      IF (M.LT.1 .OR. N.LT.1) GO TO 690
      IF (NE.LT.1) GO TO 700
      IF (LA.LT.NE) THEN
         INFO(3) = NE
         GO TO 710
      END IF

C Initial printing
      IF (MP.GT.0 .AND. ICNTL(3).GT.2) THEN
         WRITE (MP,'(/2(A,I6),A,I8,A,I8/A,1P,4E10.2/A,7I4)')
     +     ' Entering MA50AD with M =',M,' N =',N,' NE =',NE,' LA =',LA,
     +     ' CNTL =',(CNTL(I),I=1,4),' ICNTL =',(ICNTL(I),I=1,7)
         IF (N.EQ.1 .OR. ICNTL(3).GT.3) THEN
            DO 10 J = 1,N - 1
               IF (IQ(J).LT.IQ(J+1)) WRITE (MP,
     +             '(A,I5,(T13,3(1P,E12.4,I5)))') ' Column',J,
     +             (A(II),IRN(II),II=IQ(J),IQ(J+1)-1)
   10       CONTINUE
            IF (IQ(N).LE.NE) WRITE (MP,'(A,I5,(T13,3(1P,E12.4,I5)))')
     +          ' Column',N, (A(II),IRN(II),II=IQ(N),NE)
         ELSE
            IF (IQ(1).LT.IQ(2)) WRITE (MP,
     +          '(A,I5,(T13,3(1P,E12.4,I5)))') ' Column',1,
     +          (A(II),IRN(II),II=IQ(1),IQ(2)-1)
         END IF
         IF (ICNTL(7).EQ.2) THEN
            WRITE (MP,'(A,(T10,10(I7)))') ' IP = ',IP
            WRITE (MP,'(A,(T10,10(I7)))') ' IFIRST = ',IFIRST
         END IF
      END IF

C Initialization of counts etc.
      MINC = 1
      NERED = NE
      U = MIN(CNTL(2),ONE)
      U = MAX(U,ZERO)
      MSRCH = ICNTL(4)
      IF (MSRCH.EQ.0) MSRCH = N
      JLAST = N - ICNTL(6)
      IF (JLAST.LT.1 .OR. JLAST.GT.N) JLAST = N
      NULLI = 0
      NULLJ = 0
      MORD = 0
      NORD = 0
      NDROP = 0
      NEFACT = 0
      DO 20 I = 1,N - 1
         LENC(I) = IQ(I+1) - IQ(I)
   20 CONTINUE
      LENC(N) = NE + 1 - IQ(N)

      IF (CNTL(3).GT.ZERO) THEN
C Drop small entries
         NERED = 0
         DO 40 J = 1,N
            I = IQ(J)
            IQ(J) = NERED + 1
            DO 30 II = I,I + LENC(J) - 1
               IF (ABS(A(II)).GE.CNTL(3)) THEN
                  NERED = NERED + 1
                  A(NERED) = A(II)
                  IRN(NERED) = IRN(II)
               ELSE
                  INFO(6) = INFO(6) + 1
               END IF
   30       CONTINUE
            LENC(J) = NERED + 1 - IQ(J)
   40    CONTINUE
      END IF

      IF (ICNTL(7).EQ.2) THEN
C Column order specified - copy the row ordering array
         DO 50 I = 1,M
            NEXTR(I) = IP(I)
   50    CONTINUE
C Check ICNTL(4)
         IF (ICNTL(4).NE.1) GO TO 740
      END IF

      DISPR = NERED + 1
      DISPC = NERED + 1
C
C Set up row oriented storage.
      DO 60 I = 1,M
         IW(I) = 0
         LENR(I) = 0
         JFIRST(I) = 0
   60 CONTINUE
C Calculate row counts.
      DO 70 II = 1,NERED
         I = IRN(II)
         LENR(I) = LENR(I) + 1
   70 CONTINUE
C Set up row pointers so that IP(i) points to position after end
C     of row i in row file.
      IP(1) = LENR(1) + 1
      DO 80 I = 2,M
         IP(I) = IP(I-1) + LENR(I)
   80 CONTINUE
C Generate row file.
      DO 100 J = 1,N
         I = IQ(J)
         DO 90 II = I,I + LENC(J) - 1
            I = IRN(II)
C Check for duplicate entry.
            IF (IW(I).EQ.J) GO TO 720
            IW(I) = J
            IPOS = IP(I) - 1
            JCN(IPOS) = J
            IP(I) = IPOS
   90    CONTINUE
  100 CONTINUE
      DO 110 I = 1,M
         IW(I) = 0
  110 CONTINUE

C Check for zero rows and (unless ICNTL(4) > 0), compute chains of rows
C    with equal numbers of entries.
      IF (ICNTL(4).LE.0) THEN
         DO 120 I = 1,N
            IFIRST(I) = 0
  120    CONTINUE
         DO 130 I = M,1,-1
            NE1 = LENR(I)
            IF (NE1.GT.0) THEN
               IFIR = IFIRST(NE1)
               IFIRST(NE1) = I
               LASTR(I) = 0
               NEXTR(I) = IFIR
               IF (IFIR.GT.0) LASTR(IFIR) = I
            ELSE
               IP(I) = -M + NULLI
               NULLI = NULLI + 1
            END IF
  130    CONTINUE
      ELSE
         DO 140 I = M,1,-1
            NE1 = LENR(I)
            IF (NE1.EQ.0) THEN
               IP(I) = -M + NULLI
               NULLI = NULLI + 1
            END IF
  140    CONTINUE
      END IF
C Check for zero columns and compute chains of columns with equal
C   numbers of entries.
      DO 150 J = N,1,-1
         NE1 = LENC(J)
         IF (NE1.EQ.0) THEN
            IF (ICNTL(7).NE.2) THEN
               IF (J.LE.JLAST) THEN
                  NORD = NORD + 1
                  IQ(J) = -NORD
                  IF (NORD.EQ.JLAST) THEN
C We have ordered the first N - ICNTL(6) columns.
                     NORD = NORD + NULLJ
                     JLAST = N
                     NULLJ = 0
                  END IF
               ELSE
                  NULLJ = NULLJ + 1
                  IQ(J) = - (JLAST+NULLJ)
               END IF
               LASTC(J) = 0
               NEXTC(J) = 0
            END IF
         ELSE
            IFIR = JFIRST(NE1)
            JFIRST(NE1) = J
            NEXTC(J) = IFIR
            LASTC(J) = 0
            IF (IFIR.GT.0) LASTC(IFIR) = J
         END IF
  150 CONTINUE
      IF (INFO(6).EQ.0) THEN
         NULLC = NORD + NULLJ
         NULLR = NULLI
      END IF

C
C **********************************************
C ****    Start of main elimination loop    ****
C **********************************************
      DO 630 PIVOT = 1,N
C Check to see if reduced matrix should be considered as full.
         IF (NERED.GE. (MIN(CNTL(1),ONE)*(N-NORD))*
     +       (M-MORD)) GO TO 640

         IF (ICNTL(7).EQ.2) THEN
C Column order specified - choose the pivot within the column
            IPIV = 0
            J = IFIRST(PIVOT)
            IF (J.LT.1 .OR. J.GT.N) GO TO 730
            IF (IQ(J).LT.0) GO TO 730
            LEN = LENC(J)
            IF (LEN.LE.0) GO TO 320
            ALEN = LEN - 1
            I1 = IQ(J)
            I2 = I1 + LEN - 1
C Find largest entry in column
            II = IDAMAX(LEN,A(I1),1)
            MAXENT = ABS(A(I1+II-1))
C Is every entry in the column below the pivot threshold?
            IF (MAXENT.LE.CNTL(4)) GO TO 320
            AU = MAX(MAXENT*U,CNTL(4))
C Scan column for pivot
            DO 160 II = I1,I2
               IF (ABS(A(II)).LT.AU) GO TO 160
C Candidate satisfies threshold criterion.
               I = IRN(II)
               IF (IPIV.NE.0) THEN
                  IF (NEXTR(I).GE.NEXTR(IPIV)) GO TO 160
               END IF
               CPIV = ALEN*(LENR(I)-1)
               IJPOS = II
               IPIV = I
               JPIV = J
  160       CONTINUE
            GO TO 330
         END IF

C Find the least number of entries in a row or column (column only if
C   the Zlatev strategy is in use)
         LEN = MINC
         DO 170 MINC = LEN,M - MORD
            IF (JFIRST(MINC).NE.0) GO TO 180
            IF (ICNTL(4).LE.0) THEN
               IF (IFIRST(MINC).NE.0) GO TO 180
            END IF
  170    CONTINUE

C Find the next pivot or a column whose entries are all very small.
C CPIV is the Markowitz cost of the best pivot so far and PIVRAT is the
C      ratio of its absolute value to that of the largest entry in its
C      column.
  180    CPIV = M
         CPIV = CPIV*N
         PIVRAT = ZERO
C Examine columns/rows in order of ascending count.
         ISRCH = 0
         DO 300 LEN = MINC,M - MORD
            ALEN = LEN - 1
C Jump if Markowitz count cannot be bettered.
            IF (CPIV.LE.ALEN**2 .AND. ICNTL(4).LE.0) GO TO 310
            IJ = JFIRST(LEN)
C Scan columns with LEN entries.
            DO 220 IDUMMY = 1,N
C If no more columns with LEN entries, exit loop.
               IF (IJ.LE.0) GO TO 230
               J = IJ
               IJ = NEXTC(J)
               IF (J.GT.JLAST) GO TO 220
C Column J is now examined.
C First calculate multiplier threshold level.
               MAXENT = ZERO
               I1 = IQ(J)
               I2 = I1 + LEN - 1
               II = IDAMAX(LEN,A(I1),1)
               MAXENT = ABS(A(I1+II-1))
C Exit loop if every entry in the column is below the pivot threshold.
               IF (MAXENT.LE.CNTL(4)) GO TO 320
               AU = MAX(MAXENT*U,CNTL(4))
C If diagonal pivoting requested, look for diagonal entry.
               IF (ICNTL(7).EQ.1) THEN
                  DO 190 II = I1,I2
                     IF (IRN(II).EQ.J) GO TO 200
  190             CONTINUE
                  GO TO 220
  200             I1 = II
                  I2 = II
               END IF
C Scan column for possible pivots
               DO 210 II = I1,I2
                  IF (ABS(A(II)).LT.AU) GO TO 210
C Candidate satisfies threshold criterion.
                  I = IRN(II)
                  COST = ALEN*(LENR(I)-1)
                  IF (COST.GT.CPIV) GO TO 210
                  PIVR = ABS(A(II))/MAXENT
                  IF (COST.EQ.CPIV) THEN
                     IF (PIVR.LE.PIVRAT) GO TO 210
                  END IF
C Best pivot so far is found.
                  CPIV = COST
                  IJPOS = II
                  IPIV = I
                  JPIV = J
                  IF (CPIV.LE.ALEN**2 .AND. ICNTL(4).LE.0) GO TO 330
                  PIVRAT = PIVR
  210          CONTINUE
C Increment number of columns searched.
               ISRCH = ISRCH + 1
C Jump if we have searched the number of columns stipulated and found a
C   pivot.
               IF (ISRCH.GE.MSRCH) THEN
                  IF (PIVRAT.GT.ZERO) GO TO 330
               END IF
  220       CONTINUE
C
C Rows with LEN entries now examined.
  230       IF (ICNTL(4).GT.0) GO TO 300
            IF (CPIV.LE.ALEN*(ALEN+1)) GO TO 310
            IF (LEN.GT.N-NORD) GO TO 300
            IJ = IFIRST(LEN)
            DO 290 IDUMMY = 1,M
               IF (IJ.EQ.0) GO TO 300
               I = IJ
               IJ = NEXTR(IJ)
               J1 = IP(I)
               J2 = J1 + LEN - 1
C If diagonal pivoting requested, look for diagonal entry.
               IF (ICNTL(7).EQ.1) THEN
                  DO 240 JJ = J1,J2
                     IF (JCN(JJ).EQ.I) GO TO 250
  240             CONTINUE
                  GO TO 290
  250             J1 = JJ
                  J2 = JJ
               END IF
C Scan row I.
               DO 280 JJ = J1,J2
                  J = JCN(JJ)
                  IF (J.GT.JLAST) GO TO 280
                  COST = ALEN*(LENC(J)-1)
                  IF (COST.GE.CPIV) GO TO 280
C Pivot has best Markowitz count so far. Now check its suitability
C     on numerical grounds by examining other entries in its column.
                  I1 = IQ(J)
                  I2 = I1 + LENC(J) - 1
                  II = IDAMAX(LENC(J),A(I1),1)
                  MAXENT = ABS(A(I1+II-1))
                  DO 260 II = I1,I2 - 1
                     IF (IRN(II).EQ.I) GO TO 270
  260             CONTINUE
  270             JPOS = II
C Exit loop if every entry in the column is below the pivot threshold.
                  IF (MAXENT.LE.CNTL(4)) GO TO 320
                  IF (ABS(A(JPOS)).LT.MAXENT*U) GO TO 280
C Candidate satisfies threshold criterion.
                  CPIV = COST
                  IPIV = I
                  JPIV = J
                  IJPOS = JPOS
                  PIVRAT = ABS(A(JPOS))/MAXENT
                  IF (CPIV.LE.ALEN*(ALEN+1)) GO TO 330
  280          CONTINUE

  290       CONTINUE
C
  300    CONTINUE
  310    IF (PIVRAT.GT.ZERO) GO TO 330
C No pivot found. Switch to full matrix processing.
         INFO(1) = INFO(1) + 2
         IF (MP.GT.0) WRITE (MP,'(A/A)')
     +       ' Warning message from MA50AD: no suitable diagonal pivot',
     +       ' found, so switched to full matrix processing.'
         GO TO 640

C Every entry in the column is below the pivot threshold.
  320    IPIV = 0
         JPIV = J

C The pivot has now been found in position (IPIV,JPIV) in location
C     IJPOS in column file or all entries of column JPIV are very small
C     (IPIV=0).
C Update row and column ordering arrays to correspond with removal
C     of the active part of the matrix. Also update NEFACT.
  330    NEFACT = NEFACT + LENC(JPIV)
         PIVBEG = IQ(JPIV)
         PIVEND = PIVBEG + LENC(JPIV) - 1
         NORD = NORD + 1
         NORD1 = NORD
         IF (NORD.EQ.JLAST) THEN
C We have ordered the first N - ICNTL(6) columns.
            NORD = NORD + NULLJ
            JLAST = N
            NULLJ = 0
         END IF
         IF (ICNTL(4).LE.0) THEN
C Remove active rows from their row ordering chains.
            DO 340 II = PIVBEG,PIVEND
               I = IRN(II)
               LR = LASTR(I)
               NR = NEXTR(I)
               IF (NR.NE.0) LASTR(NR) = LR
               IF (LR.EQ.0) THEN
                  NE1 = LENR(I)
                  IFIRST(NE1) = NR
               ELSE
                  NEXTR(LR) = NR
               END IF
  340       CONTINUE
         END IF
         IF (IPIV.GT.0) THEN
C NEPR is number of entries in strictly U part of pivot row.
            NEPR = LENR(IPIV) - 1
            NEFACT = NEFACT + NEPR
            RINFO(1) = RINFO(1) + CPIV*2 + LENR(IPIV)
            J1 = IP(IPIV)
C Remove active columns from their column ordering chains.
            DO 350 JJ = J1,J1 + NEPR
               J = JCN(JJ)
               LC = LASTC(J)
               NC = NEXTC(J)
               IF (NC.NE.0) LASTC(NC) = LC
               IF (LC.EQ.0) THEN
                  NE1 = LENC(J)
                  JFIRST(NE1) = NC
               ELSE
                  NEXTC(LC) = NC
               END IF
  350       CONTINUE
C Move pivot to beginning of pivot column.
            IF (PIVBEG.NE.IJPOS) THEN
               ASW = A(PIVBEG)
               A(PIVBEG) = A(IJPOS)
               A(IJPOS) = ASW
               IRN(IJPOS) = IRN(PIVBEG)
               IRN(PIVBEG) = IPIV
            END IF
         ELSE
            NEPR = 0
            NE1 = LENC(JPIV)
            IF (CNTL(3).GT.ZERO) NDROP = NDROP + NE1
            IF (NE1.GT.0) THEN
C Remove column of small entries from its column ordering chain.
               LC = LASTC(JPIV)
               NC = NEXTC(JPIV)
               IF (NC.NE.0) LASTC(NC) = LC
               IF (LC.EQ.0) THEN
                  JFIRST(NE1) = NC
               ELSE
                  NEXTC(LC) = NC
               END IF
            END IF
         END IF
C
C Set up IW array so that IW(i) holds the relative position of row i
C    entry from beginning of pivot column.
         DO 360 II = PIVBEG + 1,PIVEND
            I = IRN(II)
            IW(I) = II - PIVBEG
  360    CONTINUE
C LENPIV is length of strictly L part of pivot column.
         LENPIV = PIVEND - PIVBEG
C
C Remove pivot column (including pivot) from row oriented file.
         DO 390 II = PIVBEG,PIVEND
            I = IRN(II)
            LENR(I) = LENR(I) - 1
            J1 = IP(I)
C J2 is last position in old row.
            J2 = J1 + LENR(I)
            DO 370 JJ = J1,J2 - 1
               IF (JCN(JJ).EQ.JPIV) GO TO 380
  370       CONTINUE
  380       JCN(JJ) = JCN(J2)
            JCN(J2) = 0
  390    CONTINUE

C For each active column, add the appropriate multiple of the pivot
C     column to it.
C We loop on the number of entries in the pivot row since the position
C     of this row may change because of compresses.
         DO 600 EYE = 1,NEPR
            J = JCN(IP(IPIV)+EYE-1)
C Search column J for entry to be eliminated, calculate multiplier,
C     and remove it from column file.
C  IDROP is the number of nonzero entries dropped from column J
C        because these fall beneath tolerance level.
            IDROP = 0
            JBEG = IQ(J)
            JEND = JBEG + LENC(J) - 1
            DO 400 II = JBEG,JEND - 1
               IF (IRN(II).EQ.IPIV) GO TO 410
  400       CONTINUE
  410       AMULT = -A(II)/A(IQ(JPIV))
            A(II) = A(JEND)
            IRN(II) = IRN(JEND)
            LENC(J) = LENC(J) - 1
            IRN(JEND) = 0
            JEND = JEND - 1
C Jump if pivot column is a singleton.
            IF (LENPIV.EQ.0) GO TO 600
C Now perform necessary operations on rest of non-pivot column J.
            IOP = 0
C Innermost loop.
CDIR$ IVDEP
            DO 420 II = JBEG,JEND
               I = IRN(II)
               IF (IW(I).GT.0) THEN
C Row i is involved in the pivot column.
                  IOP = IOP + 1
                  PIVCOL = IQ(JPIV) + IW(I)
C Flag IW(I) to show that the operation has been done.
                  IW(I) = -IW(I)
                  A(II) = A(II) + AMULT*A(PIVCOL)
               END IF
  420       CONTINUE

            IF (CNTL(3).GT.ZERO) THEN
C  Run through non-pivot column compressing column so that entries less
C      than CNTL(3) are not stored. All entries less than CNTL(3) are
C      also removed from the row structure.
               JNEW = JBEG
               DO 450 II = JBEG,JEND
                  IF (ABS(A(II)).GE.CNTL(3)) THEN
                     A(JNEW) = A(II)
                     IRN(JNEW) = IRN(II)
                     JNEW = JNEW + 1
                  ELSE
C  Remove non-zero entry from row structure.
                     I = IRN(II)
                     J1 = IP(I)
                     J2 = J1 + LENR(I) - 1
                     DO 430 JJ = J1,J2 - 1
                        IF (JCN(JJ).EQ.J) GO TO 440
  430                CONTINUE
  440                JCN(JJ) = JCN(J2)
                     JCN(J2) = 0
                     LENR(I) = LENR(I) - 1
                  END IF
  450          CONTINUE
               DO 460 II = JNEW,JEND
                  IRN(II) = 0
  460          CONTINUE
               IDROP = JEND + 1 - JNEW
               JEND = JNEW - 1
               LENC(J) = LENC(J) - IDROP
               NERED = NERED - IDROP
               INFO(6) = INFO(6) + IDROP
            END IF

C IFILL is fill-in left to do to non-pivot column J.
            IFILL = LENPIV - IOP
            NERED = NERED + IFILL
            INFO(3) = MAX(INFO(3),NERED+LENC(J))

C Treat no-fill case
            IF (IFILL.EQ.0) THEN
CDIR$ IVDEP
               DO 470 II = PIVBEG + 1,PIVEND
                  I = IRN(II)
                  IW(I) = -IW(I)
  470          CONTINUE
               GO TO 600
            END IF

C See if there is room for fill-in at end of the column.
            DO 480 IPOS = JEND + 1,MIN(JEND+IFILL,DISPC-1)
               IF (IRN(IPOS).NE.0) GO TO 490
  480       CONTINUE
            IF (IPOS.EQ.JEND+IFILL+1) GO TO 540
            IF (JEND+IFILL+1.LE.LA+1) THEN
               DISPC = JEND + IFILL + 1
               GO TO 540
            END IF
            IPOS = LA
            DISPC = LA + 1
C JMORE more spaces for fill-in are required.
  490       JMORE = JEND + IFILL - IPOS + 1
C We now look in front of the column to see if there is space for
C     the rest of the fill-in.
            DO 500 IPOS = JBEG - 1,MAX(JBEG-JMORE,1),-1
               IF (IRN(IPOS).NE.0) GO TO 510
  500       CONTINUE
            IPOS = IPOS + 1
            IF (IPOS.EQ.JBEG-JMORE) GO TO 520
C Column must be moved to the beginning of available storage.
  510       IF (DISPC+LENC(J)+IFILL.GT.LA+1) THEN
               INFO(2) = INFO(2) + 1
               CALL MA50DD(LA,A,IRN,IQ,N,DISPC,.TRUE.)
               JBEG = IQ(J)
               JEND = JBEG + LENC(J) - 1
               PIVBEG = IQ(JPIV)
               PIVEND = PIVBEG + LENC(JPIV) - 1
               IF (DISPC+LENC(J)+IFILL.GT.LA+1) GO TO 705
            END IF
            IPOS = DISPC
            DISPC = DISPC + LENC(J) + IFILL
C Move non-pivot column J.
  520       IQ(J) = IPOS
            DO 530 II = JBEG,JEND
               A(IPOS) = A(II)
               IRN(IPOS) = IRN(II)
               IPOS = IPOS + 1
               IRN(II) = 0
  530       CONTINUE
            JBEG = IQ(J)
            JEND = IPOS - 1
C Innermost fill-in loop which also resets IW.
C We know at this stage that there are IFILL positions free after JEND.
  540       IDROP = 0
            DO 580 II = PIVBEG + 1,PIVEND
               I = IRN(II)
               INFO(3) = MAX(INFO(3),NERED+LENR(I)+1)
               IF (IW(I).LT.0) THEN
                  IW(I) = -IW(I)
                  GO TO 580
               END IF
               ANEW = AMULT*A(II)
               IF (ABS(ANEW).LT.CNTL(3)) THEN
                  IDROP = IDROP + 1
               ELSE
                  JEND = JEND + 1
                  A(JEND) = ANEW
                  IRN(JEND) = I

C Put new entry in row file.
                  IEND = IP(I) + LENR(I)
                  IF (IEND.LT.DISPR) THEN
                     IF (JCN(IEND).EQ.0) GO TO 560
                  ELSE
                     IF (DISPR.LE.LA) THEN
                        DISPR = DISPR + 1
                        GO TO 560
                     END IF
                  END IF
                  IF (IP(I).GT.1) THEN
                     IF (JCN(IP(I)-1).EQ.0) THEN
C Copy row forward
                        IEND = IEND - 1
                        DO 545 JJ = IP(I),IEND
                           JCN(JJ-1) = JCN(JJ)
  545                   CONTINUE
                        IP(I) = IP(I) - 1
                        GO TO 560
                     END IF
                  END IF
                  IF (DISPR+LENR(I).GT.LA) THEN
C Compress.
                     INFO(2) = INFO(2) + 1
                     CALL MA50DD(LA,A,JCN,IP,M,DISPR,.FALSE.)
                     IF (DISPR+LENR(I).GT.LA) GO TO 705
                  END IF
C Copy row to first free position.
                  J1 = IP(I)
                  J2 = IP(I) + LENR(I) - 1
                  IP(I) = DISPR
                  DO 550 JJ = J1,J2
                     JCN(DISPR) = JCN(JJ)
                     JCN(JJ) = 0
                     DISPR = DISPR + 1
  550             CONTINUE
                  IEND = DISPR
                  DISPR = IEND + 1
  560             JCN(IEND) = J
                  LENR(I) = LENR(I) + 1
C End of adjustment to row file.
               END IF
  580       CONTINUE
            INFO(6) = INFO(6) + IDROP
            NERED = NERED - IDROP
            DO 590 II = 1,IDROP
               IRN(JEND+II) = 0
  590       CONTINUE
            LENC(J) = LENC(J) + IFILL - IDROP
C End of scan of pivot row.
  600    CONTINUE


C Remove pivot row from row oriented storage and update column
C     ordering arrays.  Remember that pivot row no longer includes
C     pivot.
         DO 610 EYE = 1,NEPR
            JJ = IP(IPIV) + EYE - 1
            J = JCN(JJ)
            JCN(JJ) = 0
            NE1 = LENC(J)
            LASTC(J) = 0
            IF (NE1.GT.0) THEN
               IFIR = JFIRST(NE1)
               JFIRST(NE1) = J
               NEXTC(J) = IFIR
               IF (IFIR.NE.0) LASTC(IFIR) = J
               MINC = MIN(MINC,NE1)
            ELSE IF (ICNTL(7).NE.2) THEN
               IF (INFO(6).EQ.0) NULLC = NULLC + 1
               IF (J.LE.JLAST) THEN
                  NORD = NORD + 1
                  IQ(J) = -NORD
                  IF (NORD.EQ.JLAST) THEN
C We have ordered the first N - ICNTL(6) columns.
                     NORD = NORD + NULLJ
                     JLAST = N
                     NULLJ = 0
                  END IF
               ELSE
                  NULLJ = NULLJ + 1
                  IQ(J) = - (JLAST+NULLJ)
               END IF
            END IF
  610    CONTINUE
         NERED = NERED - NEPR

C Restore IW and remove pivot column from column file.
C    Record the row permutation in IP(IPIV) and the column
C    permutation in IQ(JPIV), flagging them negative so that they
C    are not confused with real pointers in compress routine.
         IF (IPIV.NE.0) THEN
            LENR(IPIV) = 0
            IW(IPIV) = 0
            IRN(PIVBEG) = 0
            MORD = MORD + 1
            PIVBEG = PIVBEG + 1
            IP(IPIV) = -MORD
         END IF
         NERED = NERED - LENPIV - 1
         DO 620 II = PIVBEG,PIVEND
            I = IRN(II)
            IW(I) = 0
            IRN(II) = 0
            NE1 = LENR(I)
            IF (NE1.EQ.0) THEN
               IF (INFO(6).EQ.0) NULLR = NULLR + 1
               IP(I) = -M + NULLI
               NULLI = NULLI + 1
            ELSE IF (ICNTL(4).LE.0) THEN
C Adjust row ordering arrays.
               IFIR = IFIRST(NE1)
               LASTR(I) = 0
               NEXTR(I) = IFIR
               IFIRST(NE1) = I
               IF (IFIR.NE.0) LASTR(IFIR) = I
               MINC = MIN(MINC,NE1)
            END IF
  620    CONTINUE
         IQ(JPIV) = -NORD1
  630 CONTINUE
C We may drop through this loop with NULLI nonzero.

C ********************************************
C ****    End of main elimination loop    ****
C ********************************************

C Complete the permutation vectors
  640 INFO(5) = MORD + MIN(M-MORD-NULLI,N-NORD-NULLJ)
      DO 650 L = 1,MIN(M-MORD,N-NORD)
         RINFO(1) = RINFO(1) + M - MORD - L + 1 +
     +                   REAL(M-MORD-L)*(N-NORD-L)*2
  650 CONTINUE
      NP = NORD
      INFO(4) = 2 + NEFACT + M*2 + MAX(N-NORD+M-MORD,
     +          (N-NORD)*(M-MORD))
      INFO(6) = INFO(6) + NDROP
      INFO(7) = M - MORD
      DO 660 L = 1,M
         IF (IP(L).LT.0) THEN
            IP(L) = -IP(L)
         ELSE
            MORD = MORD + 1
            IP(L) = MORD
         END IF
  660 CONTINUE
      DO 670 L = 1,N
         IF (IQ(L).LT.0) THEN
            LASTC(L) = -IQ(L)
         ELSE
            IF (NORD.EQ.JLAST) NORD = NORD + NULLJ
            NORD = NORD + 1
            LASTC(L) = NORD
         END IF
  670 CONTINUE
C Store the inverse permutation
      DO 680 L = 1,N
         IQ(LASTC(L)) = L
  680 CONTINUE

C Test for rank deficiency
      IF (INFO(5).LT.MIN(M,N)) INFO(1) = INFO(1) + 1

      IF (MP.GT.0 .AND. ICNTL(3).GT.2) THEN
         WRITE (MP,'(A,I6,A,F12.1/A,7I8)') ' Leaving MA50AD with NP =',
     +     NP,' RINFO(1) =',RINFO(1),' INFO =',(INFO(I),I=1,7)
         IF (ICNTL(3).GT.3) THEN
            WRITE (MP,'(A,(T6,10(I7)))') ' IP = ',IP
            WRITE (MP,'(A,(T6,10(I7)))') ' IQ = ',IQ
         END IF
      END IF

      GO TO 750

C Error conditions.
  690 INFO(1) = -1
      IF (LP.GT.0) WRITE (LP,'(/A/(2(A,I8)))')
     +    ' **** Error return from MA50AD ****',' M =',M,' N =',N
      GO TO 750
  700 INFO(1) = -2
      IF (LP.GT.0) WRITE (LP,'(/A/(A,I10))')
     +    ' **** Error return from MA50AD ****',' NE =',NE
      GO TO 750
  705 INFO(4) =  NEFACT + NERED
      INFO(6) = INFO(6) + NDROP
  710 INFO(1) = -3
      IF (LP.GT.0) WRITE (LP,'(/A/A,I9,A,I9)')
     +    ' **** Error return from MA50AD ****',
     +    ' LA  must be increased from',LA,' to at least',INFO(3)
      GO TO 750
  720 INFO(1) = -4
      IF (LP.GT.0) WRITE (LP,'(/A/(3(A,I9)))')
     +    ' **** Error return from MA50AD ****',' Entry in row',I,
     +    ' and column',J,' duplicated'
      GO TO 750
  730 INFO(1) = -5
      IF (LP.GT.0) WRITE (LP,'(/A/(3(A,I9)))')
     +    ' **** Error return from MA50AD ****',' Fault in component ',
     +    PIVOT,' of column permutation given in IFIRST'
      GO TO 750
  740 INFO(1) = -6
      IF (LP.GT.0) WRITE (LP,'(/A/(3(A,I9)))')
     +    ' **** Error return from MA50AD ****',' ICNTL(4) = ',ICNTL(4),
     +    ' when ICNTL(6) = 2'
  750 END


      SUBROUTINE MA50BD(M,N,NE,JOB,AA,IRNA,IPTRA,CNTL,ICNTL,IP,IQ,NP,
     +                  LFACT,FACT,IRNF,IPTRL,IPTRU,W,IW,INFO,RINFO)
C MA50B/BD factorizes the matrix in AA/IRNA/IPTRA as P L U Q where
C     P and Q are permutations, L is lower triangular, and U is unit
C     upper triangular. The prior information that it uses depends on
C     the value of the parameter JOB.
C
      INTEGER M,N,NE,JOB
      DOUBLE PRECISION AA(NE)
      INTEGER IRNA(NE),IPTRA(N)
      DOUBLE PRECISION CNTL(10)
      INTEGER ICNTL(20),IP(M),IQ(*),NP,LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER IRNF(LFACT),IPTRL(N),IPTRU(N)
      DOUBLE PRECISION W(M)
      INTEGER IW(M+2*N),INFO(15)
      DOUBLE PRECISION RINFO(10)
C
C M is an integer variable that must be set to the number of rows.
C      It is not altered by the subroutine.
C N is an integer variable that must be set to the number of columns.
C      It is not altered by the subroutine.
C NE is an integer variable that must be set to the number of entries
C      in the input matrix.  It is not altered by the subroutine.
C JOB is an integer variable that must be set to the value 1, 2, or 3.
C     If JOB is equal to 1 and any of the first NP recommended pivots
C      fails to satisfy the threshold pivot tolerance, the row is
C      interchanged with the earliest row in the recommended sequence
C      that does satisfy the tolerance. Normal row interchanges are
C      performed in the last N-NP columns.
C     If JOB is equal to 2, then M, N, NE, IRNA, IPTRA, IP, IQ,
C      LFACT, NP, IRNF, IPTRL, and IPTRU must be unchanged since a
C      JOB=1 entry for the same matrix pattern and no interchanges are
C      performed among the first NP pivots; if ICNTL(6) > 0, the first
C      N-ICNTL(6) columns of AA must also be unchanged.
C     If JOB is equal to 3, ICNTL(6) must be in the range 1 to N-1.
C      The effect is as for JOB=2 except that interchanges are
C      performed.
C     JOB is not altered by the subroutine.
C AA is an array that holds the entries of the matrix and
C      is not altered.
C IRNA is an integer array of length NE that must be set to hold the
C      row indices of the corresponding entries in AA. It is not
C      altered.
C IPTRA is an integer array that holds the positions of the starts of
C      the columns of AA. It is not altered by the subroutine.
C CNTL  must be set by the user as follows and is not altered.
C     CNTL(2) determines the balance between pivoting for sparsity and
C       for stability, values near zero emphasizing sparsity and values
C       near one emphasizing stability.
C     CNTL(3) If this is set to a positive value, any entry whose
C       modulus is less than CNTL(3) will be dropped from the factors.
C       The factorization will then require less storage but will be
C       inaccurate.
C     CNTL(4)  Any entry of the reduced matrix whose modulus is less
C       than or equal to CNTL(4) will be regarded as zero from the
C        point of view of rank.
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
C     ICNTL(5) The block size to be used for full-matrix processing.
C       If <=0, the BLAS1 version is used.
C       If =1, the BLAS2 version is used.
C     ICNTL(6) If N > ICNTL(6) > 0, only the columns of A that
C       correspond to the last ICNTL(6) columns of the permuted matrix
C       may change prior to an entry with JOB > 1.
C     ICNTL(8) If this has value 1, there is no attempt to compute a
C       recommended value for LFACT if it is too small.
C IP is an integer array. If JOB=1, it must be set so that IP(I) < IP(J)
C      if row I is recommended to precede row J in the pivot sequence.
C      If JOB>1, it need not be set. If JOB=1 or JOB=3, IP(I) is set
C      to -K when row I is chosen for pivot K and IP is eventually
C      reset to recommend the chosen pivot sequence to a subsequent
C      JOB=1 entry. If JOB=2, IP is not be referenced.
C IQ is an integer array that must be set so that either IQ(J) is the
C      column in position J in the pivot sequence, J=1,2,...,N,
C      or IQ(1)=0 and the columns are taken in natural order.
C      It is not altered by the subroutine.
C NP is an integer variable that holds the number of columns to be
C      processed in packed storage. It is not altered by the subroutine.
C LFACT is an integer variable set to the size of FACT and IRNF.
C      It is not altered by the subroutine.
C FACT is an array that need not be set on a JOB=1 entry and must be
C      unchanged since the previous entry if JOB>1. On return, FACT(1)
C      holds the value of CNTL(3) used, FACT(2) will holds the value
C      of CNTL(4) used, FACT(3:IPTRL(N)) holds the packed part of L/U
C      by columns, and the full part of L/U is held by columns
C      immediately afterwards. U has unit diagonal entries, which are
C      not stored. In each column of the packed part, the entries of
C      U precede the entries of L; also the diagonal entries of L
C      head each column of L and are reciprocated.
C IRNF is an integer array of length LFACT that need not be set on
C      a JOB=1 entry and must be unchanged since the previous entry
C      if JOB>1. On exit, IRNF(1) holds the number of dropped entries,
C      IRNF(2) holds the number of rows MF in full storage,
C      IRNF(3:IPTRL(N)) holds the row numbers of the packed part
C      of L/U, IRNF(IPTRL(N)+1:IPTRL(N)+MF) holds the row indices
C      of the full part of L/U, and IRNF(IPTRL(N)+MF+I), I=1,2,..,N-NP
C      holds the vector IPIV output by MA50GD.
C      If JOB=2, IRNF will be unaltered.
C IPTRL is an integer array that need not be set on a JOB=1 entry and
C     must be unchanged since the previous entry if JOB>1.
C     For J = 1,..., NP, IPTRL(J) holds the position in
C     FACT and IRNF of the end of column J of L.
C     For J = NP+1,..., N, IPTRL(J) is equal to IPTRU(J).
C IPTRU is an integer array that need not be set on a JOB=1 entry and
C     must be unchanged since the previous entry if JOB>1.
C     For J = 1,..., N, IPTRU(J) holds the position in
C     FACT and IRNF of the end of the packed part of column J of U.
C W is an array of length M used as workspace for holding
C      the expanded form of a sparse vector.
C IW is an integer array of length M+2*N used as workspace.
C INFO need not be set on entry. On exit, it holds the following:
C    INFO(1) A negative value will indicate an error return and a
C       positive value a warning. Possible nonzero values are:
C      -1  M < 1 or N < 1.
C      -2  NE < 0.
C      -3  Insufficient space.
C      -4  There are duplicated entries.
C      -5  JOB < 1, 3 when ICNTL(6)=0, or > 3.
C      -6  JOB = 2, but entries were dropped in the corresponding JOB=1
C          entry.
C      -7  NP < 0 or NP > N.
C     -(7+K) Pivot too small in column K when JOB=2.
C      +1  Rank deficient.
C    INFO(4) Minimum storage required to factorize matrix
C            if INFO(1) >= 0. Recommended value for LFACT
C            if ICNTL(8) = 0 and INFO(1) = -3.
C    INFO(5) Computed rank of the matrix.
C    INFO(6) Number of entries dropped from the data structure.
C    INFO(7) Number of rows processed in full storage.
C RINFO need not be set on entry. On exit, RINFO(1) holds the number of
C    floating-point operations performed.

      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0D0,ONE=1.0D0)
      DOUBLE PRECISION AMULT,ASW
      INTEGER BEGCOL
      LOGICAL DROP
      INTEGER ENDCOL,EYE,EYE1,I,IA1,IA2,IF1,IF2,II,IL1,IL2,IPIV,IQPIV,
     +        IU1,IU2,ISW,J,JDUMMY,JJ,JLAST,K,LP
      DOUBLE PRECISION MAXENT
      INTEGER MF,MORD,MP,NEU,NF,NULLC
      DOUBLE PRECISION PIVLIM
      INTEGER RANK
      DOUBLE PRECISION U
C AMULT Temporary variable used to store current multiplier.
C ASW Temporary variable used when swopping two real quantities.
C BEGCOL is pointer to beginning of section of column when pruning.
C DROP True if any entries dropped from current column.
C ENDCOL is pointer to end of section of column when pruning.
C EYE Running position for current column.
C EYE1 Position of the start of second current column.
C I Temporary variable holding row number. Also used as index in DO
C     loops used in initialization of arrays.
C IA1 Position of the start of the current column in AA.
C IA2 Position of the end of the current column in AA.
C IF1 Position of the start of the full submatrix.
C IF2 Position of the end of the full submatrix.
C II Running position for current column.
C IL1 Position of the first entry of the current column of L.
C IL2 Position of the last entry of the current column of L.
C IPIV Position of the pivot in FACT and IRNF.
C IQPIV Recommended position of the pivot row in the pivot sequence.
C IU1 Position of the start of current column of U.
C IU2 Position of the end of the current column of U.
C ISW Temporary variable used when swopping two integer quantities.
C J Temporary variable holding column number.
C JDUMMY DO index not referenced in the loop.
C JJ Running position for current column.
C JLAST The lesser of NP and the last column of A for which no new
C     factorization operations are needed.
C K Temporary variable holding the current pivot step in the elimination
C LP Unit for error messages.
C MAXENT Temporary variable used to hold value of largest entry in
C    column.
C MF Number of rows in full block.
C MORD Number of rows ordered.
C MP Unit for diagnostic messages.
C NEU Number of entries omitted from U and the full block in order to
C    calculate INFO(4) (0 unless INFO(1)=-3).
C NF Number of columns in full block.
C NULLC Number of columns found null before dropping any elements.
C PIVLIM Limit on pivot size.
C RANK Value returned by MA50E/ED or MA50F/FD
C U Used to hold local copy of CNTL(2), changed if necessary so that it
C    is in range.
C
      EXTERNAL MA50ED,MA50FD,MA50GD
      INTRINSIC ABS,MAX,MIN
C LAPACK subroutine for triangular factorization.

      INFO(1) = 0
      INFO(4) = 0
      INFO(5) = 0
      INFO(6) = 0
      INFO(7) = 0
      RINFO(1) = ZERO
      LP = ICNTL(1)
      MP = ICNTL(2)
      IF (ICNTL(3).LE.0) LP = 0
      IF (ICNTL(3).LE.1) MP = 0

C Check input values
      IF (M.LT.1 .OR. N.LT.1) THEN
         INFO(1) = -1
         IF (LP.GT.0) WRITE (LP,'(/A/A,I8,A,I8)')
     +       ' **** Error return from MA50BD ****',' M =',M,' N =',N
         GO TO 550
      END IF
      IF (NE.LE.0) THEN
         INFO(1) = -2
         IF (LP.GT.0) WRITE (LP,'(/A/A,I6)')
     +       ' **** Error return from MA50BD ****',' NE =',NE
         GO TO 550
      END IF
      IF (NP.LT.0 .OR. NP.GT.N) THEN
         INFO(1) = -7
         IF (LP.GT.0) WRITE (LP,'(/A/A,I8,A,I8)')
     +       ' **** Error return from MA50BD ****',' NP =',NP,' N =',N
         GO TO 550
      END IF
      IF (LFACT.LT.MAX(M,NE+2)) THEN
         INFO(4) = MAX(M,NE+2)
         GO TO 520
      END IF
      IF (JOB.EQ.1) THEN
      ELSE IF (JOB.EQ.2 .OR. JOB.EQ.3) THEN
         IF (IRNF(1).NE.0) THEN
            INFO(1) = -6
            IF (LP.GT.0) WRITE (LP,'(/A/A,I1,A)')
     +          ' **** Error return from MA50BD ***',' Call with JOB=',
     +          JOB,' follows JOB=1 call in which entries were dropped'
            GO TO 550
         END IF
      ELSE
         INFO(1) = -5
         IF (LP.GT.0) WRITE (LP,'(/A/A,I2)')
     +       ' **** Error return from MA50BD ****',' JOB =',JOB
         GO TO 550
      END IF

C Print input data
      IF (MP.GT.0) THEN
         IF (ICNTL(3).GT.2) WRITE (MP,
     +       '(/2(A,I6),A,I8,A,I3/A,I8,A,I7/A,1P,4E10.2/A,7I8)')
     +       ' Entering MA50BD with M =',M,' N =',N,' NE =',NE,' JOB =',
     +       JOB,' LFACT =',LFACT,' NP =',NP,' CNTL =',(CNTL(I),I=1,4),
     +       ' ICNTL =',(ICNTL(I),I=1,7)
         IF (ICNTL(3).GT.3) THEN
            WRITE (MP,'(A,(T6,10(I7)))') ' IP = ',IP
            IF (IQ(1).GT.0) THEN
               WRITE (MP,'(A,(T6,10(I7)))') ' IQ = ',(IQ(J),J=1,N)
            ELSE
               WRITE (MP,'(A,(T6,I7))') ' IQ = ',IQ(1)
            END IF
            DO 10 J = 1,N - 1
               IF (IPTRA(J).LT.IPTRA(J+1)) WRITE (MP,
     +             '(A,I5,(T13,3(1P,E12.4,I5)))') ' Column',J,
     +             (AA(II),IRNA(II),II=IPTRA(J),IPTRA(J+1)-1)
   10       CONTINUE
            IF (IPTRA(N).LE.NE) WRITE (MP,
     +          '(A,I5,(T13,3(1P,E12.4,I5)))') ' Column',N,
     +          (AA(II),IRNA(II),II=IPTRA(N),NE)
         END IF
      END IF

C Initializations.
      JLAST = 0
      NULLC = 0
      IF (JOB.GT.1 .AND. ICNTL(6).GT.0 .AND.
     +    ICNTL(6).LT.N) JLAST = MIN(NP,N-ICNTL(6))

      U = MIN(CNTL(2),ONE)
      U = MAX(U,ZERO)
      DO 20 I = 1,M
         IW(I+N) = 0
         W(I) = ZERO
   20 CONTINUE
      MORD = 0
      IF1 = LFACT + 1
      IF2 = 0
      NF = N - NP
      MF = 0
      IL2 = 2
      IF (JLAST.GT.0) IL2 = IPTRL(JLAST)
      NEU = 0

C Jump if JOB is equal to 2.
      IF (JOB.EQ.2) GO TO 370

      IF (JOB.EQ.3) THEN
C Reconstruct IP and set MORD
         DO 30 J = 1,NP
            IA1 = IPTRU(J) + 1
            IF (IA1.GT.IPTRL(J)) GO TO 30
            IF (J.LE.JLAST) THEN
               MORD = MORD + 1
               IP(IRNF(IA1)) = -J
            ELSE
               IP(IRNF(IA1)) = J
            END IF
   30    CONTINUE
         MF = IRNF(2)
         IA1 = IPTRL(N)
         DO 40 J = 1,MF
            IP(IRNF(IA1+J)) = NP + J
   40    CONTINUE
      END IF

C Store copies of column ends ready for pruning
      DO 50 K = 1,JLAST
         IW(M+N+K) = IPTRL(K)
   50 CONTINUE

C Each pass through this main loop processes column K.
      DO 310 K = JLAST + 1,N
         DROP = .FALSE.
         IF (K.EQ.NP+1) THEN
C Set up data structure for full part.
            MF = M - MORD
            IF1 = LFACT + 1 - MF
            II = 0
            DO 60 I = 1,M
               IF (IP(I).GT.0) THEN
                  IW(I+N) = N
                  IRNF(IF1+II) = I
                  II = II + 1
                  IP(I) = NP + II
               END IF
   60       CONTINUE
            IF1 = LFACT + 1 - MAX(MF*NF,MF+NF)
            IF2 = IF1 - 1 + MF*MAX(0,JLAST-NP)
         END IF
         J = K
         IF (IQ(1).GT.0) J = IQ(K)
         IA1 = IPTRA(J)
         IA2 = NE
         IF (J.NE.N) IA2 = IPTRA(J+1) - 1
         IU1 = IL2 + 1
         IU2 = IU1 - 1
         IL1 = IF1 - 1 + IA1 - IA2
         IL2 = IL1 - 1
         INFO(4) = MAX(INFO(4),NEU+LFACT-IL1+IU2+M+1)
         IF (IL1-IU2.LE.M) THEN
            IF (INFO(1).NE.-3) THEN
               INFO(1) = -3
               IF (ICNTL(8).NE.0) GO TO 480
C Get rid of U info.
               NEU = IL2 + LFACT + 1 - MF - IF1
               IF1 = LFACT + 1 - MF
               IF2 = IF1 - 1
               IL2 = 0
               EYE = 0
               DO 80 J = 1,MIN(K-1,NP)
                  IU2 = IPTRU(J)
                  IPTRU(J) = EYE
                  IL2 = IPTRL(J)
                  NEU = NEU + IU2 - IL2
                  DO 70 II = IU2 + 1,IL2
                     EYE = EYE + 1
                     IRNF(EYE) = IRNF(II)
                     FACT(EYE) = FACT(II)
   70             CONTINUE
                  IPTRL(J) = EYE
                  IW(M+N+J) = EYE
   80          CONTINUE
               IU1 = EYE + 1
               IU2 = EYE
               IL1 = IF1 - 1 + IA1 - IA2
               IL2 = IL1 - 1
            END IF
C Quit if LFACT is much too small
            IF (IL1-IU2.LE.M) GO TO 480
         END IF
C Load column K of AA into full vector W and into the back of IRNF.
C Check for duplicates.
         EYE = IL1
         DO 90 II = IA1,IA2
            I = IRNA(II)
            IF (IW(I+N).EQ.-1) GO TO 540
            IW(I+N) = -1
            W(I) = AA(II)
            IRNF(EYE) = I
            EYE = EYE + 1
   90    CONTINUE
C Depth first search to find topological order for triangular solve
C     and structure of column K of L/U
C IW(J) is used to hold a pointer to next entry in column J
C     during the depth-first search at stage K, J = 1,..., N.
C IW(I+N) is set to K when row I has been processed, and to N for rows
C     of the full part once column NP has been passed. It is also
C     used for backtracking, a negative value being used to point to the
C     previous row in the chain.
C IW(M+N+I) is set to the position in FACT and IRNF of the end of the
C     active part of the column after pruning.  It is initially set to
C     IPTRL(I) and is flagged negative when column has been pruned.
C Set IPTRL temporarily for column K so that special code is
C     not required to process this column.
         IPTRL(K) = EYE - 1
         IW(M+N+K) = EYE - 1
C IW(K) is set to beginning of original column K.
         IW(K) = IL1
         J = K
C The outer loop of the depth-first search is executed once for column
C      K and twice for each entry in the upper-triangular part of column
C      K (once to initiate a search in the corresponding column and
C      once when the search in the column is finished).
         DO 120 JDUMMY = 1,2*K
C Look through column J of L (or column K of A). All the entries
C     are entries of the filled-in column K. Store new entries of the
C     lower triangle and continue until reaching an entry of the upper
C     triangle.
            DO 100 II = IW(J),ABS(IW(M+N+J))
               I = IRNF(II)
C Jump if index I already encountered in column K or is in full part.
               IF (IW(I+N).GE.K) GO TO 100
               IF (IP(I).LE.0) GO TO 110
C Entry is in lower triangle. Flag it and store it in L.
               IW(I+N) = K
               IL1 = IL1 - 1
               IRNF(IL1) = I
  100       CONTINUE
            IF (J.EQ.K) GO TO 130
C Flag J, put its row index into U, and backtrack
            IU2 = IU2 + 1
            I = IRNF(IPTRU(J)+1)
            IRNF(IU2) = I
            J = -IW(I+N)
            IW(I+N) = K
            GO TO 120
C Entry in upper triangle.  Move search to corresponding column.
  110       IW(I+N) = -J
            IW(J) = II + 1
            J = -IP(I)
            IW(J) = IPTRU(J) + 2
  120    CONTINUE
C Run through column K of U in the lexicographical order that was just
C     constructed, performing elimination operations.
  130    DO 150 II = IU2,IU1,-1
            I = IRNF(II)
            J = -IP(I)
C Add multiple of column J of L to column K
            EYE1 = IPTRU(J) + 1
            IF (ABS(W(I)).LT.CNTL(3)) GO TO 150
            AMULT = -W(I)*FACT(EYE1)
C Note we are storing negative multipliers
            W(I) = AMULT
            DO 140 EYE = EYE1 + 1,IPTRL(J)
               I = IRNF(EYE)
               W(I) = W(I) + AMULT*FACT(EYE)
  140       CONTINUE
            RINFO(1) = RINFO(1) + ONE + 2*(IPTRL(J)-EYE1)
  150    CONTINUE

C Unload reals of column of U and set pointer
         IF (CNTL(3).GT.ZERO) THEN
            EYE = IU1
            DO 160 II = IU1,IU2
               I = IRNF(II)
               IF (ABS(W(I)).LT.CNTL(3)) THEN
                  INFO(6) = INFO(6) + 1
               ELSE
                  IRNF(EYE) = -IP(I)
                  FACT(EYE) = W(I)
                  EYE = EYE + 1
               END IF
               W(I) = ZERO
  160       CONTINUE
            IU2 = EYE - 1
         ELSE
            DO 170 II = IU1,IU2
               I = IRNF(II)
               IRNF(II) = -IP(I)
               FACT(II) = W(I)
               W(I) = ZERO
  170       CONTINUE
         END IF
         IF (INFO(1).EQ.-3) THEN
            NEU = NEU + IU2 - IU1 + 1
            IU2 = IU1 - 1
         END IF
         IPTRU(K) = IU2
         IF (K.LE.NP) THEN
C Find the largest entry in the column and drop any small entries
            MAXENT = ZERO
            IF (CNTL(3).GT.ZERO) THEN
               EYE = IL1
               DO 180 II = IL1,IL2
                  I = IRNF(II)
                  IF (ABS(W(I)).LT.CNTL(3)) THEN
                     INFO(6) = INFO(6) + 1
                     W(I) = ZERO
                     DROP = .TRUE.
                  ELSE
                     IRNF(EYE) = I
                     EYE = EYE + 1
                     MAXENT = MAX(ABS(W(I)),MAXENT)
                  END IF
  180          CONTINUE
               IL2 = EYE - 1
            ELSE
               DO 190 II = IL1,IL2
                  MAXENT = MAX(ABS(W(IRNF(II))),MAXENT)
  190          CONTINUE
            END IF
C Unload column of L, performing pivoting and moving indexing
C      information.
            PIVLIM = U*MAXENT
            EYE = IU2
            IQPIV = M + N
            IF (IL1.GT.IL2) NULLC = NULLC + 1
            DO 200 II = IL1,IL2
               I = IRNF(II)
               EYE = EYE + 1
               IRNF(EYE) = I
               FACT(EYE) = W(I)
               W(I) = ZERO
C Find position of pivot
               IF (ABS(FACT(EYE)).GE.PIVLIM) THEN
                  IF (ABS(FACT(EYE)).GT.CNTL(4)) THEN
                     IF (IP(I).LT.IQPIV) THEN
                        IQPIV = IP(I)
                        IPIV = EYE
                     END IF
                  END IF
               END IF
  200       CONTINUE
            IL1 = IU2 + 1
            IL2 = EYE
            IF (IL1.LE.IL2) THEN
C Column is not null
               IF (IQPIV.EQ.M+N) THEN
C All entries in the column are too small to be pivotal. Drop them all.
                  IF (CNTL(3).GT.ZERO) INFO(6) = INFO(6) + EYE - IU2
                  IL2 = IU2
               ELSE
                  IF (IL1.NE.IPIV) THEN
C Move pivot to front of L
                     ASW = FACT(IPIV)
                     FACT(IPIV) = FACT(IL1)
                     FACT(IL1) = ASW
                     ISW = IRNF(IL1)
                     IRNF(IL1) = IRNF(IPIV)
                     IRNF(IPIV) = ISW
                  END IF
C Reciprocate pivot
                  INFO(5) = INFO(5) + 1
                  FACT(IL1) = ONE/FACT(IL1)
                  RINFO(1) = RINFO(1) + ONE
C Record pivot row
                  MORD = MORD + 1
                  IP(IRNF(IL1)) = -K
               END IF
            END IF
         ELSE
C Treat column as full
            IL2 = IPTRU(K)
CDIR$ IVDEP
            DO 210 II = LFACT - MF + 1,LFACT
               I = IRNF(II)
               IF2 = IF2 + 1
               FACT(IF2) = W(I)
               W(I) = ZERO
  210       CONTINUE
            IF (INFO(1).EQ.-3) IF2 = IF2 - MF
         END IF
         IW(M+N+K) = IL2
         IPTRL(K) = IL2
         IF (DROP) GO TO 310
C Scan columns involved in update of column K and remove trailing block.
         DO 300 II = IU1,IU2
            I = IRNF(II)
C Jump if column already pruned.
            IF (IW(M+N+I).LT.0) GO TO 300
            BEGCOL = IPTRU(I) + 2
            ENDCOL = IPTRL(I)
C Scan column to see if there is an entry in the current pivot row.
            IF (K.LE.NP) THEN
               DO 220 JJ = BEGCOL,ENDCOL
                  IF (IP(IRNF(JJ)).EQ.-K) GO TO 230
  220          CONTINUE
               GO TO 300
            END IF
C Sort the entries so that those in rows already pivoted (negative IP
C    values) precede the rest.
  230       DO 280 JDUMMY = BEGCOL,ENDCOL
               JJ = BEGCOL
               DO 240 BEGCOL = JJ,ENDCOL
                  IF (IP(IRNF(BEGCOL)).GT.0) GO TO 250
  240          CONTINUE
               GO TO 290
  250          JJ = ENDCOL
               DO 260 ENDCOL = JJ,BEGCOL,-1
                  IF (IP(IRNF(ENDCOL)).LT.0) GO TO 270
  260          CONTINUE
               GO TO 290
  270          ASW = FACT(BEGCOL)
               FACT(BEGCOL) = FACT(ENDCOL)
               FACT(ENDCOL) = ASW
               J = IRNF(BEGCOL)
               IRNF(BEGCOL) = IRNF(ENDCOL)
               IRNF(ENDCOL) = J
               BEGCOL = BEGCOL + 1
               ENDCOL = ENDCOL - 1
  280       CONTINUE
  290       IW(M+N+I) = -ENDCOL
  300    CONTINUE
  310 CONTINUE
      IF (N.EQ.NP) THEN
C Set up data structure for the (null) full part.
         MF = M - MORD
         IF1 = LFACT + 1 - MF
         II = 0
         DO 320 I = 1,M
            IF (IP(I).GT.0) THEN
               IW(I+N) = N
               IRNF(IF1+II) = I
               II = II + 1
               IP(I) = NP + II
            END IF
  320    CONTINUE
         IF1 = LFACT + 1 - MAX(MF*NF,MF+NF)
         IF2 = IF1 - 1 + MF*MAX(0,JLAST-NP)
      END IF
      IF (INFO(5).EQ.MIN(M,N)) THEN
C Restore sign of IP
         DO 330 I = 1,M
            IP(I) = ABS(IP(I))
  330    CONTINUE
      ELSE
C Complete IP
         MORD = NP
         DO 340 I = 1,M
            IF (IP(I).LT.0) THEN
               IP(I) = -IP(I)
            ELSE
               MORD = MORD + 1
               IP(I) = MORD
            END IF
  340    CONTINUE
      END IF
      IRNF(1) = INFO(6)
      IRNF(2) = MF
      INFO(7) = MF
      FACT(1) = CNTL(3)
      FACT(2) = CNTL(4)
      IF (INFO(1).EQ.-3) GO TO 520
C Move full part forward
      IF2 = IF2 - MF*NF
      DO 350 II = 1,MF*NF
         FACT(IL2+II) = FACT(IF1-1+II)
  350 CONTINUE
      DO 360 II = 1,MF
         IRNF(IL2+II) = IRNF(LFACT-MF+II)
  360 CONTINUE
      IF1 = IL2 + 1
      GO TO 440
C
C Fast factor (JOB = 2)
C Each pass through this main loop processes column K.
  370 MF = IRNF(2)
      IF1 = IPTRL(N) + 1
      IF2 = IF1 - 1
      DO 430 K = JLAST + 1,N
         J = K
         IF (IQ(1).GT.0) J = IQ(K)
         IA1 = IPTRA(J)
         IA2 = NE
         IF (J.NE.N) IA2 = IPTRA(J+1) - 1
         IU1 = IL2 + 1
         IU2 = IPTRU(K)
         IL1 = IU2 + 1
         IL2 = IPTRL(K)
C Load column K of A into full vector W
         DO 380 II = IA1,IA2
            W(IRNA(II)) = AA(II)
  380    CONTINUE
C Run through column K of U in lexicographical order, performing
C      elimination operations.
         DO 400 II = IU2,IU1,-1
            J = IRNF(II)
            I = IRNF(IPTRU(J)+1)
C Add multiple of column J of L to column K
            EYE1 = IPTRU(J) + 1
            AMULT = -W(I)*FACT(EYE1)
C Note we are storing negative multipliers
            FACT(II) = AMULT
            W(I) = ZERO
            DO 390 EYE = EYE1 + 1,IPTRL(J)
               I = IRNF(EYE)
               W(I) = W(I) + AMULT*FACT(EYE)
  390       CONTINUE
            RINFO(1) = RINFO(1) + ONE + 2*(IPTRL(J)-EYE1)
  400    CONTINUE
         IF (K.LE.NP) THEN
            IF (IL1.LE.IL2) THEN
C Load column of L.
CDIR$ IVDEP
               DO 410 II = IL1,IL2
                  I = IRNF(II)
                  FACT(II) = W(I)
                  W(I) = ZERO
  410          CONTINUE
C Test pivot. Note that this is the only numerical test when JOB = 2.
               IF (ABS(FACT(IL1)).LE.CNTL(4)) THEN
                  GO TO 530
               ELSE
C Reciprocate pivot
                  INFO(5) = INFO(5) + 1
                  FACT(IL1) = ONE/FACT(IL1)
                  RINFO(1) = RINFO(1) + ONE
               END IF
            END IF
         ELSE
C Treat column as full
            DO 420 II = IF1,IF1 + MF - 1
               I = IRNF(II)
               IF2 = IF2 + 1
               FACT(IF2) = W(I)
               W(I) = ZERO
  420       CONTINUE
         END IF
  430 CONTINUE
      INFO(4) = MAX(IF1+MF+NF-1,IF2)

  440 IF (MF.GT.0 .AND. NF.GT.0) THEN
C Factorize full block
         IF (ICNTL(5).GT.1) CALL MA50GD(MF,NF,FACT(IF1),MF,ICNTL(5),
     +                                  CNTL(4),IRNF(IF1+MF),RANK)
         IF (ICNTL(5).EQ.1) CALL MA50FD(MF,NF,FACT(IF1),MF,CNTL(4),
     +                                  IRNF(IF1+MF),RANK)
         IF (ICNTL(5).LE.0) CALL MA50ED(MF,NF,FACT(IF1),MF,CNTL(4),
     +                                  IRNF(IF1+MF),RANK)
         INFO(5) = INFO(5) + RANK
         DO 450 I = 1,MIN(MF,NF)
            RINFO(1) = RINFO(1) + MF - I + 1 + REAL(MF-I)*(NF-I)*2
  450    CONTINUE
      END IF
      IF (INFO(5).LT.MIN(M,N)) INFO(1) = 1
      IF (MP.GT.0 .AND. ICNTL(3).GT.2) THEN
         WRITE (MP,'(A,I6,A,F12.1/A,I3,A,4I8)')
     +     ' Leaving MA50BD with IRNF(2) =',IRNF(2),
     +     ' RINFO(1) =',RINFO(1),
     +     ' INFO(1) =',INFO(1),' INFO(4:7) =', (INFO(J),J=4,7)
         IF (ICNTL(3).GT.3) THEN
            IF (JOB.NE.2) WRITE (MP,'(A,(T6,10(I7)))') ' IP = ',IP
            DO 460 J = 1,N
               IF (J.GT.1) THEN
                  IF (IPTRL(J-1).LT.IPTRU(J)) WRITE (MP,
     +                '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column',J,
     +                ' of U', (FACT(II),IRNF(II),II=IPTRL(J-1)+1,
     +                IPTRU(J))
               END IF
               IF (IPTRU(J).LT.IPTRL(J)) WRITE (MP,
     +             '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column',J,' of L',
     +              (FACT(II),IRNF(II),II=IPTRU(J)+1,IPTRL(J))
  460       CONTINUE
            WRITE (MP,'(A)') ' Full part'
            WRITE (MP,'((6I12))') (IRNF(IF1+MF+J),J=0,NF-1)
            DO 470 I = 0,MF - 1
               WRITE (MP,'(I4,1P,6E12.4:/(4X,1P,6E12.4))') IRNF(IF1+I),
     +            (FACT(IF1+I+J*MF),J=0,NF-1)
  470       CONTINUE
         END IF
      END IF
      GO TO 550

C Error conditions
C LFACT is much too small or ICNTL(8)/=0. Patch up IP and quit.
  480 DO 490 I = 1,M
         IW(I) = 0
  490 CONTINUE
      DO 500 I = 1,M
         IF (IP(I).GT.0) THEN
            IW(IP(I)) = I
         ELSE
            IP(I) = -IP(I)
         END IF
  500 CONTINUE
      DO 510 I = 1,M
         IF (IW(I).GT.0) THEN
            IP(IW(I)) = K
            K = K + 1
         END IF
  510 CONTINUE
  520 INFO(1) = -3
      IF (LP.GT.0) THEN
         WRITE (LP,'(/A)')' **** Error return from MA50BD **** '
         IF (ICNTL(8).EQ.0) THEN
            WRITE (LP,'(A,I7,A,I7)')' LFACT must be increased from',
     +         LFACT,' to at least',INFO(4)
         ELSE
            WRITE (LP,'(A,I7)')' LFACT must be increased from',LFACT
         END IF
      END IF
      GO TO 550
  530 INFO(1) = - (7+K)
      IF (LP.GT.0) WRITE (LP,'(/A/A,I6,A)')
     +    ' **** Error return from MA50BD **** ',
     +    ' Small pivot found in column',K,' of the permuted matrix.'
      GO TO 550
  540 INFO(1) = -4
      IF (LP.GT.0) WRITE (LP,'(/A/(3(A,I9)))')
     +    ' **** Error return from MA50BD ****',' Entry in row',I,
     +    ' and column',J,' duplicated'
  550 END

      SUBROUTINE MA50CD(M,N,ICNTL,IQ,NP,TRANS,LFACT,FACT,IRNF,IPTRL,
     +                  IPTRU,B,X,W,INFO)
C MA50C/CD uses the factorization produced by
C     MA50B/BD to solve A x = b or (A trans) x = b.
C
      INTEGER M,N,ICNTL(20),IQ(*),NP
      LOGICAL TRANS
      INTEGER LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER IRNF(LFACT),IPTRL(N),IPTRU(N)
      DOUBLE PRECISION B(*),X(*),W(*)
      INTEGER INFO(15)
C
C M  is an integer variable set to the number of rows.
C     It is not altered by the subroutine.
C N  is an integer variable set to the number of columns.
C     It is not altered by the subroutine.
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
C     ICNTL(5) must be set to control the level of BLAS used:
C       0 Level 1 BLAS.
C      >0 Level 2 BLAS.
C IQ is an integer array holding the permutation Q.
C     It is not altered by the subroutine.
C NP is an integer variable that must be unchanged since calling
C     MA50B/BD. It holds the number of rows and columns in packed
C     storage. It is not altered by the subroutine.
C TRANS a logical variable thatmust be set to .TRUE. if (A trans)x = b
C     is to be solved and to .FALSE. if A x = b is to be solved.
C     TRANS is not altered by the subroutine.
C LFACT is an integer variable set to the size of FACT and IRNF.
C     It is not altered by the subroutine.
C FACT is an array that must be unchanged since calling MA50B/BD. It
C     holds the packed part of L/U by columns, and the full part of L/U
C     by columns. U has unit diagonal entries, which are not stored, and
C     the signs of the off-diagonal entries are inverted.  In the packed
C     part, the entries of U precede the entries of L; also the diagonal
C     entries of L head each column of L and are reciprocated.
C     FACT is not altered by the subroutine.
C IRNF is an integer array that must be unchanged since calling
C     MA50B/BD. It holds the row numbers of the packed part of L/U, and
C     the row numbers of the full part of L/U.
C     It is not altered by the subroutine.
C IPTRL is an integer array that must be unchanged since calling
C     MA50B/BD. For J = 1,..., NP, IPTRL(J) holds the position in
C     FACT and IRNF of the end of column J of L.
C     For J = NP+1,..., N, IPTRL(J) is equal to IPTRU(J).
C     It is not altered by the subroutine.
C IPTRU is an integer array that must be unchanged since calling
C     MA50B/BD. For J = 1,..., N, IPTRU(J) holds the position in
C     FACT and IRNF of the end of the packed part of column J of U.
C     It is not altered by the subroutine.
C B is an array that must be set to the vector b.
C     It is not altered.
C X is an array that need not be set on entry. On return, it holds the
C    solution x.
C W is a work array of length max(M,N).
C INFO need not be set on entry. On exit, it holds the following:
C    INFO(1) A nonzero value will indicate an error return. Possible
C      nonzero values are:
C      -1  M < 1 or N < 1

      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0D0)
      INTEGER I,II,IA1,IF1,J,LP,MF,MP,NF
      DOUBLE PRECISION PROD
C I Temporary variable holding row number.
C II Position of the current entry in IRNF.
C IA1 Position of the start of the current row or column.
C IF1 Position of the start of the full part of U.
C J Temporary variable holding column number.
C LP Unit for error messages.
C MF Number of rows held in full format.
C MP Unit for diagnostic messages.
C NF Number of columns held in full format.
C PROD Temporary variable used to accumulate inner products.

      EXTERNAL MA50HD

      LP = ICNTL(1)
      MP = ICNTL(2)
      IF (ICNTL(3).LE.0) LP = 0
      IF (ICNTL(3).LE.1) MP = 0

C Make some simple checks
      IF (M.LT.1 .OR. N.LT.1) GO TO 250

      IF (MP.GT.0 .AND. ICNTL(3).GT.2) WRITE (MP,
     +    '(/2(A,I6),A,I4,A,L2)') ' Entering MA50CD with M=',M,' N =',N,
     +    ' NP =',NP,' TRANS =',TRANS
      IF1 = IPTRL(N) + 1
      MF = IRNF(2)
      NF = N - NP
      IF (MP.GT.0 .AND. ICNTL(3).GT.2) WRITE (MP,
     +    '(A,I5,A,I5)') ' Size of full submatrix',MF,' by',NF
      IF (MP.GT.0 .AND. ICNTL(3).GT.3) THEN
         DO 10 J = 1,N
            IF (J.GT.1) THEN
               IF (IPTRL(J-1).LT.IPTRU(J)) WRITE (MP,
     +             '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column',J,' of U',
     +              (FACT(II),IRNF(II),II=IPTRL(J-1)+1,IPTRU(J))
            END IF
            IF (IPTRU(J).LT.IPTRL(J)) WRITE (MP,
     +          '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column',J,' of L',
     +          (FACT(II),IRNF(II),II=IPTRU(J)+1,IPTRL(J))
   10    CONTINUE
         WRITE (MP,'(A)') ' Full part'
         WRITE (MP,'((6I12))') (IRNF(IF1+MF+J),J=0,NF-1)
         DO 20 I = 0,MF - 1
            WRITE (MP,'(I4,1P,6E12.4:/(4X,1P,6E12.4))') IRNF(IF1+I),
     +        (FACT(IF1+I+J*MF),J=0,NF-1)
   20    CONTINUE
      END IF

      IF (TRANS) THEN
         IF (MP.GT.0 .AND. ICNTL(3).GT.3) WRITE (MP,
     +       '(A4,5F10.4:/(4X,5F10.4))') ' B =', (B(I),I=1,N)
         IF (IQ(1).GT.0) THEN
            DO 30 I = 1,N
               W(I) = B(IQ(I))
   30       CONTINUE
         ELSE
            DO 40 I = 1,N
               W(I) = B(I)
   40       CONTINUE
         END IF
         DO 50 I = 1,M
            X(I) = ZERO
   50    CONTINUE
C Forward substitution through packed part of (U trans).
         DO 70 I = 2,N
            PROD = ZERO
            DO 60 II = IPTRL(I-1) + 1,IPTRU(I)
               PROD = PROD + FACT(II)*W(IRNF(II))
   60       CONTINUE
            W(I) = W(I) + PROD
   70    CONTINUE
C Backsubstitute through the full part of (PL) trans.
         DO 80 I = 1,NF
            X(I) = W(NP+I)
   80    CONTINUE
         IF (MF.GT.0 .AND. NF.GT.0) THEN
            CALL MA50HD(TRANS,MF,NF,FACT(IF1),MF,IRNF(IF1+MF),X,
     +                  ICNTL(5))
         ELSE
            DO 90 I = 1,MF
               X(I) = ZERO
   90       CONTINUE
         END IF
         DO 100 I = MF,1,-1
            J = IRNF(IF1+I-1)
            IF (J.NE.I) X(J) = X(I)
  100    CONTINUE
C Backsubstitute through the packed part of (PL) trans.
         DO 120 I = NP,1,-1
            IA1 = IPTRU(I) + 1
            IF (IA1.GT.IPTRL(I)) GO TO 120
            PROD = ZERO
            DO 110 II = IA1 + 1,IPTRL(I)
               PROD = PROD + FACT(II)*X(IRNF(II))
  110       CONTINUE
            X(IRNF(IA1)) = (W(I)-PROD)*FACT(IA1)
  120    CONTINUE
         IF (MP.GT.0 .AND. ICNTL(3).GT.3) WRITE (MP,
     +       '(A/(4X,5F10.4))') ' Leaving MA50CD with X =', (X(I),I=1,M)
C
      ELSE
         IF (MP.GT.0 .AND. ICNTL(3).GT.3) WRITE (MP,
     +       '(A4,5F10.4:/(4X,5F10.4))') ' B =', (B(I),I=1,M)
C Forward substitution through the packed part of PL
         DO 130 I = 1,M
            W(I) = B(I)
  130    CONTINUE
         DO 150 I = 1,NP
            IA1 = IPTRU(I) + 1
            IF (IA1.LE.IPTRL(I)) THEN
               X(I) = W(IRNF(IA1))*FACT(IA1)
               IF (X(I).NE.ZERO) THEN
CDIR$ IVDEP
                  DO 140 II = IA1 + 1,IPTRL(I)
                     W(IRNF(II)) = W(IRNF(II)) - FACT(II)*X(I)
  140             CONTINUE
               END IF
            END IF
  150    CONTINUE
C Forward substitution through the full part of PL
         IF (MF.GT.0 .AND. NF.GT.0) THEN
            DO 160 I = 1,MF
               W(I) = W(IRNF(IF1+I-1))
  160       CONTINUE
            CALL MA50HD(TRANS,MF,NF,FACT(IF1),MF,IRNF(IF1+MF),W,
     +                  ICNTL(5))
            DO 170 I = 1,NF
               X(NP+I) = W(I)
  170       CONTINUE
         ELSE
            DO 180 I = 1,NF
               X(NP+I) = ZERO
  180       CONTINUE
         END IF
C Back substitution through the packed part of U
         DO 200 J = N,MAX(2,NP+1),-1
            PROD = X(J)
CDIR$ IVDEP
            DO 190 II = IPTRL(J-1) + 1,IPTRU(J)
               X(IRNF(II)) = X(IRNF(II)) + FACT(II)*PROD
  190       CONTINUE
  200    CONTINUE
         DO 220 J = NP,2,-1
            IA1 = IPTRU(J)
            IF (IA1.GE.IPTRL(J)) THEN
               X(J) = ZERO
            ELSE
               PROD = X(J)
CDIR$ IVDEP
               DO 210 II = IPTRL(J-1) + 1,IA1
                  X(IRNF(II)) = X(IRNF(II)) + FACT(II)*PROD
  210          CONTINUE
            END IF
  220    CONTINUE
         IF (NP.GE.1 .AND. IPTRU(1).GE.IPTRL(1)) X(1) = ZERO
         IF (IQ(1).GT.0) THEN
C         Permute X
            DO 230 I = 1,N
               W(I) = X(I)
  230       CONTINUE
            DO 240 I = 1,N
               X(IQ(I)) = W(I)
  240       CONTINUE
         END IF
         IF (MP.GT.0 .AND. ICNTL(3).GT.3) WRITE (MP,
     +       '(A/(4X,5F10.4))') ' Leaving MA50CD with X =', (X(I),I=1,N)
      END IF
      RETURN
C Error condition.
  250 INFO(1) = -1
      IF (LP.GT.0) WRITE (LP,'(/A/2(A,I8))')
     +    ' **** Error return from MA50CD ****',' M =',M,' N =',N
      END

      SUBROUTINE MA50DD(LA,A,IND,IPTR,N,DISP,REALS)
C This subroutine performs garbage collection on the arrays A and IND.
C DISP is the position in arrays A/IND immediately after the data
C     to be compressed.
C     On exit, DISP equals the position of the first entry
C     after the compressed part of A/IND.
C
      INTEGER LA,N,DISP
      DOUBLE PRECISION A(LA)
      INTEGER IPTR(N)
      LOGICAL REALS
      INTEGER IND(LA)
C Local variables.
      INTEGER J,K,KN
C Set the first entry in each row(column) to the negative of the
C     row(column) and hold the column(row) index in the row(column)
C     pointer.  This enables the start of each row(column) to be
C     recognized in a subsequent scan.
      DO 10 J = 1,N
         K = IPTR(J)
         IF (K.GT.0) THEN
            IPTR(J) = IND(K)
            IND(K) = -J
         END IF
   10 CONTINUE
      KN = 0
C Go through arrays compressing to the front so that there are no
C     zeros held in positions 1 to DISP-1 of IND.
C     Reset first entry of each row(column) and the pointer array IPTR.
      DO 20 K = 1,DISP - 1
         IF (IND(K).EQ.0) GO TO 20
         KN = KN + 1
         IF (REALS) A(KN) = A(K)
         IF (IND(K).LE.0) THEN
C First entry of row(column) has been located.
            J = -IND(K)
            IND(K) = IPTR(J)
            IPTR(J) = KN
         END IF
         IND(KN) = IND(K)
   20 CONTINUE
      DISP = KN + 1
      END


      SUBROUTINE MA50ED(M,N,A,LDA,PIVTOL,IPIV,RANK)
**
      INTEGER LDA,M,N,RANK
      DOUBLE PRECISION PIVTOL

      INTEGER IPIV(N)
      DOUBLE PRECISION A(LDA,N)

*
*  Purpose
*  =======
*
*  MA50ED computes an LU factorization of a general m-by-n matrix A.

*  The factorization has the form
*     A = P * L * U * Q
*  where P is a permutation matrix of order m, L is lower triangular
*  of order m with unit diagonal elements, U is upper trapezoidal of
*  order m * n, and Q is a permutation matrix of order n.
*
*  Row interchanges are used to ensure that the entries of L do not
*  exceed 1 in absolute value. Column interchanges are used to
*  ensure that the first r diagonal entries of U exceed PIVTOL in
*  absolute value. If r < m, the last (m-r) rows of U are zero.

*  This is the Level 1 BLAS version.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 1.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 1.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the m by n matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U*Q; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  PIVTOL  (input) DOUBLE PRECISION
*          The pivot tolerance. Any entry with absolute value less
*          than or equal to PIVTOL is regarded as unsuitable to be a
*          pivot.
**
*  IPIV    (output) INTEGER array, dimension (N)
*          The permutations; for 1 <= i <= RANK, row i of the
*          matrix was interchanged with row IPIV(i); for
*          RANK + 1 <= j <= N, column j of the
*          matrix was interchanged with column -IPIV(j).
*
*  RANK    (output) INTEGER
*          The computed rank of the matrix.
*
*  =====================================================================
*
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)

      INTEGER I,J,JP,K
      LOGICAL PIVOT
* I   Row index.
* J   Current column.
* JP  Pivot position.
* K   Main loop index.
* PIVOT True if there is a pivot in current column.

      INTEGER IDAMAX
      EXTERNAL IDAMAX

      EXTERNAL DAXPY,DSCAL,DSWAP
      INTRINSIC ABS

*
      J = 1
      DO 30 K = 1,N

*        Update elements in column J.
         DO 10 I = 1,J - 1
            IF (M.GT.I) CALL DAXPY(M-I,-A(I,J),A(I+1,I),1,A(I+1,J),1)
   10    CONTINUE

*        Find pivot.
         IF (J.LE.M) THEN
            JP = J - 1 + IDAMAX(M-J+1,A(J,J),1)
            IPIV(J) = JP
            PIVOT = ABS(A(JP,J)) .GT. PIVTOL
         ELSE
            PIVOT = .FALSE.
         END IF
         IF (PIVOT) THEN

*           Apply row interchange to columns 1:N+J-K.
            IF (JP.NE.J) CALL DSWAP(N+J-K,A(J,1),LDA,A(JP,1),LDA)
*
*           Compute elements J+1:M of J-th column.
            IF (J.LT.M) CALL DSCAL(M-J,ONE/A(J,J),A(J+1,J),1)

*           Update J
            J = J + 1
*
         ELSE
*
            DO 20 I = J,M
               A(I,J) = ZERO
   20       CONTINUE
*           Apply column interchange and record it.
            IF (K.LT.N) CALL DSWAP(M,A(1,J),1,A(1,N-K+J),1)
            IPIV(N-K+J) = -J
*
         END IF
*
   30 CONTINUE

      RANK = J - 1
*
*     End of MA50ED
*
      END


      SUBROUTINE MA50FD(M,N,A,LDA,PIVTOL,IPIV,RANK)
*
*  -- This is a variant of the LAPACK routine DGETF2 --
*
      INTEGER LDA,M,N,RANK
      DOUBLE PRECISION PIVTOL

      INTEGER IPIV(N)
      DOUBLE PRECISION A(LDA,N)

*
*  Purpose
*  =======
*
*  MA50FD computes an LU factorization of a general m-by-n matrix A.

*  The factorization has the form
*     A = P * L * U * Q
*  where P is a permutation matrix of order m, L is lower triangular
*  of order m with unit diagonal elements, U is upper trapezoidal of
*  order m * n, and Q is a permutation matrix of order n.
*
*  Row interchanges are used to ensure that the entries of L do not
*  exceed 1 in absolute value. Column interchanges are used to
*  ensure that the first r diagonal entries of U exceed PIVTOL in
*  absolute value. If r < m, the last (m-r) rows of U are zero.

*  This is the Level 2 BLAS version.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 1.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 1.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the m by n matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U*Q; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  PIVTOL  (input) DOUBLE PRECISION
*          The pivot tolerance. Any entry with absolute value less
*          than or equal to PIVTOL is regarded as unsuitable to be a
*          pivot.
**
*  IPIV    (output) INTEGER array, dimension (N)
*          The permutations; for 1 <= i <= RANK, row i of the
*          matrix was interchanged with row IPIV(i); for
*          RANK + 1 <= j <= N, column j of the
*          matrix was interchanged with column -IPIV(j).
*
*  RANK    (output) INTEGER
*          The computed rank of the matrix.
*
*  =====================================================================
*
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)

      INTEGER I,J,JP,K
      LOGICAL PIVOT
* I   Row index.
* J   Current column.
* JP  Pivot position.
* K   Main loop index.
* PIVOT True if there is a pivot in current column.

      INTEGER IDAMAX
      EXTERNAL IDAMAX

      EXTERNAL DGEMV,DSCAL,DSWAP

      INTRINSIC ABS

*
      J = 1
      DO 20 K = 1,N

         IF (J.LE.M) THEN
*           Update diagonal and subdiagonal elements in column J.
            CALL DGEMV('No transpose',M-J+1,J-1,-ONE,A(J,1),LDA,A(1,J),
     +                 1,ONE,A(J,J),1)
*          Find pivot.
            JP = J - 1 + IDAMAX(M-J+1,A(J,J),1)
            IPIV(J) = JP
            PIVOT = ABS(A(JP,J)) .GT. PIVTOL
         ELSE
            PIVOT = .FALSE.
         END IF

         IF (PIVOT) THEN

*           Apply row interchange to columns 1:N+J-K.
            IF (JP.NE.J) CALL DSWAP(N+J-K,A(J,1),LDA,A(JP,1),LDA)
*
*           Compute elements J+1:M of J-th column.
            IF (J.LT.M) CALL DSCAL(M-J,ONE/A(J,J),A(J+1,J),1)

            IF (J.LT.N) THEN
*             Compute block row of U.
               CALL DGEMV('Transpose',J-1,N-J,-ONE,A(1,J+1),LDA,A(J,1),
     +                    LDA,ONE,A(J,J+1),LDA)
            END IF

*           Update J
            J = J + 1
*
         ELSE
*
            DO 10 I = J,M
               A(I,J) = ZERO
   10       CONTINUE
*           Apply column interchange and record it.
            IF (K.LT.N) CALL DSWAP(M,A(1,J),1,A(1,N-K+J),1)
            IPIV(N-K+J) = -J
*
         END IF
*
   20 CONTINUE

      RANK = J - 1
*
*     End of MA50FD
*
      END


      SUBROUTINE MA50GD(M,N,A,LDA,NB,PIVTOL,IPIV,RANK)
*
*  -- This is a variant of the LAPACK routine DGETRF --
*
      INTEGER LDA,M,N,NB,RANK
      DOUBLE PRECISION PIVTOL

      INTEGER IPIV(N)
      DOUBLE PRECISION A(LDA,N)

*
*  Purpose
*  =======
*
*  MA50GD computes an LU factorization of a general m-by-n matrix A.
*
*  The factorization has the form
*     A = P * L * U * Q
*  where P is a permutation matrix of order m, L is lower triangular
*  of order m with unit diagonal elements, U is upper trapezoidal of
*  order m * n, and Q is a permutation matrix of order n.
*
*  Row interchanges are used to ensure that the entries of L do not
*  exceed 1 in absolute value. Column interchanges are used to
*  ensure that the first r diagonal entries of U exceed PIVTOL in
*  absolute value. If r < m, the last (m-r) rows of U are zero.
*
*  This is the Level 3 BLAS version.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 1.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 1.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the m by n matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U*Q; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  NB      (input) INTEGER
*          The block size for BLAS3 processing.
*
*  PIVTOL  (input) DOUBLE PRECISION
*          The pivot tolerance. Any entry with absolute value less
*          than or equal to PIVTOL is regarded as unsuitable to be a
*          pivot.
**
*  IPIV    (output) INTEGER array, dimension (N)
*          The permutations; for 1 <= i <= RANK, row i of the
*          matrix was interchanged with row IPIV(i); for
*          RANK + 1 <= j <= N, column j of the
*          matrix was interchanged with column -IPIV(j).
*
*  RANK    (output) INTEGER
*          The computed rank of the matrix.
*
*  =====================================================================
*
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)

      INTEGER I,J,JJ,JP,J1,J2,K
      LOGICAL PIVOT
      DOUBLE PRECISION TEMP

* I   DO index for applying permutations.
* J   Current column.
* JJ  Column in which swaps occur.
* JP  Pivot position.
* J1  Column at start of current block.
* J2  Column at end of current block.
* K   Main loop index.
* PIVOT True if there is a pivot in current column.
* TEMP Temporary variable for swaps.

      EXTERNAL DGEMM,DGEMV,DSWAP,DSCAL,DTRSM,DTRSV

      INTEGER IDAMAX
      EXTERNAL IDAMAX

      INTRINSIC ABS,MIN

*
      J = 1
      J1 = 1
      J2 = MIN(N,NB)
      DO 70 K = 1,N

         IF (J.LE.M) THEN

*          Update diagonal and subdiagonal elements in column J.
            CALL DGEMV('No transpose',M-J+1,J-J1,-ONE,A(J,J1),LDA,
     +                 A(J1,J),1,ONE,A(J,J),1)

*          Find pivot.
            JP = J - 1 + IDAMAX(M-J+1,A(J,J),1)
            IPIV(J) = JP
            PIVOT = ABS(A(JP,J)) .GT. PIVTOL
         ELSE
            PIVOT = .FALSE.
         END IF

         IF (PIVOT) THEN

*           Apply row interchange to columns J1:J2
            IF (JP.NE.J) CALL DSWAP(J2-J1+1,A(J,J1),LDA,A(JP,J1),LDA)
*
*           Compute elements J+1:M of J-th column.
            IF (J.LT.M) CALL DSCAL(M-J,ONE/A(J,J),A(J+1,J),1)

            IF (J+1.LE.J2) THEN
*             Compute row of U within current block
               CALL DGEMV('Transpose',J-J1,J2-J,-ONE,A(J1,J+1),LDA,
     +                    A(J,J1),LDA,ONE,A(J,J+1),LDA)
            END IF

*           Update J
            J = J + 1
*
         ELSE

            DO 10 I = J,M
               A(I,J) = ZERO
   10       CONTINUE
*
*           Record column interchange and revise J2 if necessary
            IPIV(N-K+J) = -J
*           Apply column interchange.
            IF (K.NE.N) CALL DSWAP(M,A(1,J),1,A(1,N-K+J),1)
            IF (N-K+J.GT.J2) THEN
*              Apply operations to new column.
               DO 20 I = J1,J - 1
                  JP = IPIV(I)
                  TEMP = A(I,J)
                  A(I,J) = A(JP,J)
                  A(JP,J) = TEMP
   20          CONTINUE
               IF(J.GT.J1) CALL DTRSV('Lower','No transpose','Unit',
     +                                J-J1,A(J1,J1),LDA,A(J1,J),1)
            ELSE
               J2 = J2 - 1
            END IF
*
         END IF

         IF (J.GT.J2) THEN
*           Apply permutations to columns outside the block
            DO 40 JJ = 1,J1 - 1
               DO 30 I = J1,J2
                  JP = IPIV(I)
                  TEMP = A(I,JJ)
                  A(I,JJ) = A(JP,JJ)
                  A(JP,JJ) = TEMP
   30          CONTINUE
   40       CONTINUE
            DO 60 JJ = J2 + 1,N - K + J - 1
               DO 50 I = J1,J2
                  JP = IPIV(I)
                  TEMP = A(I,JJ)
                  A(I,JJ) = A(JP,JJ)
                  A(JP,JJ) = TEMP
   50          CONTINUE
   60       CONTINUE

            IF (K.NE.N) THEN
*              Update the Schur complement
               CALL DTRSM('Left','Lower','No transpose','Unit',J2-J1+1,
     +                    N-K,ONE,A(J1,J1),LDA,A(J1,J2+1),LDA)
               IF (M.GT.J2) CALL DGEMM('No transpose','No transpose',
     +                                 M-J2,N-K,J2-J1+1,-ONE,A(J2+1,J1),
     +                                 LDA,A(J1,J2+1),LDA,ONE,
     +                                 A(J2+1,J2+1),LDA)
            END IF

            J1 = J2 + 1
            J2 = MIN(J2+NB,N-K+J-1)

         END IF

*
   70 CONTINUE
      RANK = J - 1
*
*     End of MA50GD
*
      END

      SUBROUTINE MA50HD(TRANS,M,N,A,LDA,IPIV,B,ICNTL5)
*
*  -- This is a variant of the LAPACK routine DGETRS --
*     It handles the singular or rectangular case.
*
      LOGICAL TRANS
      INTEGER LDA,M,N,ICNTL5
      INTEGER IPIV(N)
      DOUBLE PRECISION A(LDA,N),B(*)

*
*  Purpose
*  =======
*
*  Solve a system of linear equations
*     A * X = B  or  A' * X = B
*  with a general m by n matrix A using the LU factorization computed
*  by MA50DE, MA50FD, or MA50GD.
*
*  Arguments
*  =========
*
*  TRANS   (input) LOGICAL
*          Specifies the form of the system of equations.
*          = .FALSE. :  A * X = B  (No transpose)
*          = .TRUE.  :  A'* X = B  (Transpose)
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 1.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 1.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The factors L and U from the factorization A = P*L*U
*          as computed by MA50GD.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The permutations; for 1 <= i <= RANK, row i of the
*          matrix was interchanged with row IPIV(i); for
*          RANK + 1 <= j <= N, column j of the
*          matrix was interchanged with column -IPIV(j).

*  B       (input/output) DOUBLE PRECISION array, size max(M,N)
*          On entry, the right hand side vectors B for the system of
*          linear equations.
*          On exit, the solution vectors, X.
*
*  ICNTL5  (input) INTEGER
*          0 for BLAS1 or >0 for BLAS2
*
*  =====================================================================
*
      INTEGER I,K,RANK
C I    Temporary variable.
C K    Temporary variable.
C RANK Rank of matrix.

      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0D0)
      DOUBLE PRECISION TEMP
      INTRINSIC MIN
      EXTERNAL DAXPY,DDOT,DTRSV
      DOUBLE PRECISION DDOT

*   Find the rank
      RANK = 0
      DO 10 RANK = MIN(M,N),1,-1
         IF (IPIV(RANK).GT.0) GO TO 20
   10 CONTINUE

   20 IF (.NOT.TRANS) THEN
*
*        Solve A * X = B.
*
*        Apply row interchanges to the right hand side.
         DO 30 I = 1,RANK
            K = IPIV(I)
            TEMP = B(I)
            B(I) = B(K)
            B(K) = TEMP
   30    CONTINUE
*
*        Solve L*X = B, overwriting B with X.
         IF (ICNTL5.GT.0) THEN
            IF (RANK.GT.0) CALL DTRSV('L','NoTrans','Unit',RANK,A,LDA,B,
     +                                1)
         ELSE
            DO 40 K = 1,RANK - 1
               IF (B(K).NE.ZERO) CALL DAXPY(RANK-K,-B(K),A(K+1,K),1,
     +                                B(K+1),1)
   40       CONTINUE
         END IF

*        Solve U*X = B, overwriting B with X.
         IF (ICNTL5.GT.0) THEN
            IF (RANK.GT.0) CALL DTRSV('U','NoTrans','NonUnit',RANK,A,
     +                                LDA,B,1)
         ELSE
            DO 50 K = RANK,2,-1
               IF (B(K).NE.ZERO) THEN
                  B(K) = B(K)/A(K,K)
                  CALL DAXPY(K-1,-B(K),A(1,K),1,B(1),1)
               END IF
   50       CONTINUE
            IF (RANK.GT.0) B(1) = B(1)/A(1,1)
         END IF

*        Set singular part to zero
         DO 60 K = RANK + 1,N
            B(K) = ZERO
   60    CONTINUE
*
*        Apply column interchanges to the right hand side.
         DO 70 I = RANK + 1,N
            K = -IPIV(I)
            TEMP = B(I)
            B(I) = B(K)
            B(K) = TEMP
   70    CONTINUE

      ELSE
*
*        Solve A' * X = B.

*        Apply column interchanges to the right hand side.
         DO 80 I = N,RANK + 1,-1
            K = -IPIV(I)
            TEMP = B(I)
            B(I) = B(K)
            B(K) = TEMP
   80    CONTINUE
*
*        Solve U'*X = B, overwriting B with X.
*
         IF (ICNTL5.GT.0) THEN
            IF (RANK.GT.0) CALL DTRSV('U','Trans','NonUnit',RANK,A,LDA,
     +                                B,1)
         ELSE
            IF (RANK.GT.0) B(1) = B(1)/A(1,1)
            DO 90 I = 2,RANK
               TEMP = B(I) - DDOT(I-1,A(1,I),1,B(1),1)
               B(I) = TEMP/A(I,I)
   90       CONTINUE
         END IF

*        Solve L'*X = B, overwriting B with X.
         IF (ICNTL5.GT.0) THEN
            IF (RANK.GT.0) CALL DTRSV('L','Trans','Unit',RANK,A,LDA,B,1)
         ELSE
            DO 100 I = RANK - 1,1,-1
               B(I) = B(I) - DDOT(RANK-I,A(I+1,I),1,B(I+1),1)
  100       CONTINUE
         END IF

*        Set singular part to zero
         DO 110 I = RANK + 1,M
            B(I) = ZERO
  110    CONTINUE
*
*        Apply row interchanges to the solution vectors.
         DO 120 I = RANK,1,-1
            K = IPIV(I)
            TEMP = B(I)
            B(I) = B(K)
            B(K) = TEMP
  120    CONTINUE
      END IF

      END

      SUBROUTINE MA50ID(CNTL,ICNTL)
C Set default values for the control arrays.

      DOUBLE PRECISION CNTL(10)
      INTEGER I,ICNTL(20)

      CNTL(1) = 0.5D0
      CNTL(2) = 0.1D0
      DO 10 I = 3,10
         CNTL(I) = 0.0D0
   10 CONTINUE

      ICNTL(1) = 6
      ICNTL(2) = 6
      ICNTL(3) = 1
      ICNTL(4) = 3
      ICNTL(5) = 32
      DO 20 I = 6,20
        ICNTL(I) = 0
   20 CONTINUE

      END
* COPYRIGHT (c) 1988 AEA Technology and
* Council for the Central Laboratory of the Research Councils
C Original date 14 June 2001
C  June 2001: threadsafe version of MC41
C 20/2/02 Cosmetic changes applied to reduce single/double differences

C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC71AD(N,KASE,X,EST,W,IW,KEEP)
C
C      MC71A/AD ESTIMATES THE 1-NORM OF A SQUARE MATRIX A.
C      REVERSE COMMUNICATION IS USED FOR EVALUATING
C      MATRIX-VECTOR PRODUCTS.
C
C
C         N       INTEGER
C                 THE ORDER OF THE MATRIX.  N .GE. 1.
C
C         KASE    INTEGER
C                 SET INITIALLY TO ZERO . IF N .LE. 0 SET TO -1
C                 ON INTERMEDIATE RETURN
C                 = 1 OR 2.
C                ON FINAL RETURN
C                 =  0  ,IF SUCCESS
C                 = -1  ,IF N .LE.0
C
C         X       DOUBLE PRECISION ARRAY OF DIMENSION (N)
C                 IF 1-NORM IS REQUIRED
C                 MUST BE OVERWRITTEN BY
C
C                      A*X,             IF KASE=1,
C                      TRANSPOSE(A)*X,  IF KASE=2,
C
C                 AND MC71 MUST BE RE-CALLED, WITH ALL THE OTHER
C                 PARAMETERS UNCHANGED.
C                 IF INFINITY-NORM IS REQUIRED
C                 MUST BE OVERWRITTEN BY
C
C                      TRANSPOSE(A)*X,  IF KASE=1,
C                      A*X,             IF KASE=2,
C
C                 AND MC71 MUST BE RE-CALLED, WITH ALL THE OTHER
C                 PARAMETERS UNCHANGED.
C
C         EST     DOUBLE PRECISION
C                 CONTAINS AN ESTIMATE (A LOWER BOUND) FOR NORM(A).
C
C         W       DOUBLE PRECISION ARRAY OF DIMENSION (N)
C                 = A*V,   WHERE  EST = NORM(W)/NORM(V)
C                          (V  IS NOT RETURNED).
C         IW      INTEGER(N) USED AS WORKSPACE.
C
C         KEEP    INTEGER ARRAY LENGTH 5 USED TO PRESERVE PRIVATE
C                 DATA, JUMP, ITER, J AND JLAST BETWEEN CALLS,
C                 KEEP(5) IS SPARE.
C
C      REFERENCE
C      N.J. HIGHAM (1987) FORTRAN CODES FOR ESTIMATING
C      THE ONE-NORM OF A
C      REAL OR COMPLEX MATRIX, WITH APPLICATIONS
C      TO CONDITION  ESTIMATION, NUMERICAL ANALYSIS REPORT NO. 135,
C      UNIVERSITY OF MANCHESTER, MANCHESTER M13 9PL, ENGLAND.
C
C      SUBROUTINES AND FUNCTIONS
C
C
C      INTERNAL VARIABLES
C
C     .. Parameters ..
      INTEGER ITMAX
      PARAMETER (ITMAX=5)
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION EST
      INTEGER KASE,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION W(*),X(*)
      INTEGER IW(*),KEEP(5)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ALTSGN,TEMP
      INTEGER I,ITER,J,JLAST,JUMP
C     ..
C     .. External Functions ..
      INTEGER IDAMAX
      EXTERNAL IDAMAX
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,SIGN,NINT,DBLE
C     ..
C     .. Executable Statements ..
C
      IF (N.LE.0) THEN
        KASE = -1
        RETURN

      END IF

      IF (KASE.EQ.0) THEN
        DO 10 I = 1,N
          X(I) = ONE/DBLE(N)
   10   CONTINUE
        KASE = 1
        JUMP = 1
        KEEP(1) = JUMP
        KEEP(2) = 0
        KEEP(3) = 0
        KEEP(4) = 0
        RETURN

      END IF
C
      JUMP  = KEEP(1)
      ITER  = KEEP(2)
      J     = KEEP(3)
      JLAST = KEEP(4)
C
      GO TO (100,200,300,400,500) JUMP
C
C      ................ ENTRY   (JUMP = 1)
C
  100 CONTINUE
      IF (N.EQ.1) THEN
        W(1) = X(1)
        EST = ABS(W(1))
C         ... QUIT
        GO TO 510

      END IF
C
      DO 110 I = 1,N
        X(I) = SIGN(ONE,X(I))
        IW(I) = NINT(X(I))
  110 CONTINUE
      KASE = 2
      JUMP = 2
      GO TO 1010
C
C      ................ ENTRY   (JUMP = 2)
C
  200 CONTINUE
      J = IDAMAX(N,X,1)
      ITER = 2
C
C      MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
C
  220 CONTINUE
      DO 230 I = 1,N
        X(I) = ZERO
  230 CONTINUE
      X(J) = ONE
      KASE = 1
      JUMP = 3
      GO TO 1010
C
C      ................ ENTRY   (JUMP = 3)
C
  300 CONTINUE
C
C      COPY X INTO W
C
      DO 310 I = 1,N
        W(I) = X(I)
  310 CONTINUE
      DO 320 I = 1,N
        IF (NINT(SIGN(ONE,X(I))).NE.IW(I)) GO TO 330
  320 CONTINUE
C
C      REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
      GO TO 410
C
  330 CONTINUE
      DO 340 I = 1,N
        X(I) = SIGN(ONE,X(I))
        IW(I) = NINT(X(I))
  340 CONTINUE
      KASE = 2
      JUMP = 4
      GO TO 1010
C
C      ................ ENTRY   (JUMP = 4)
C
  400 CONTINUE
      JLAST = J
      J = IDAMAX(N,X,1)
      IF ((ABS(X(JLAST)).NE.ABS(X(J))) .AND. (ITER.LT.ITMAX)) THEN
        ITER = ITER + 1
        GO TO 220

      END IF
C
C      ITERATION COMPLETE.  FINAL STAGE.
C
  410 CONTINUE
      EST = ZERO
      DO 420 I = 1,N
        EST = EST + ABS(W(I))
  420 CONTINUE
C
      ALTSGN = ONE
      DO 430 I = 1,N
        X(I) = ALTSGN* (ONE+DBLE(I-1)/DBLE(N-1))
        ALTSGN = -ALTSGN
  430 CONTINUE
      KASE = 1
      JUMP = 5
      GO TO 1010
C
C      ................ ENTRY   (JUMP = 5)
C
  500 CONTINUE
      TEMP = ZERO
      DO 520 I = 1,N
        TEMP = TEMP + ABS(X(I))
  520 CONTINUE
      TEMP = 2.0*TEMP/DBLE(3*N)
      IF (TEMP.GT.EST) THEN
C
C      COPY X INTO W
C
        DO 530 I = 1,N
          W(I) = X(I)
  530   CONTINUE
        EST = TEMP
      END IF
C
  510 KASE = 0
C
 1010 CONTINUE
      KEEP(1) = JUMP
      KEEP(2) = ITER
      KEEP(3) = J
      KEEP(4) = JLAST
      RETURN
C
      END
* COPYRIGHT (c) 1977 AEA Technology
* Original date 8 Oct 1992
C######8/10/92 Toolpack tool decs employed.
C######8/10/92 D version created by name change only.
C 13/3/02 Cosmetic changes applied to reduce single/double differences
C
C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC21AD(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW)
C     .. Scalar Arguments ..
      INTEGER LICN,N,NUMNZ
C     ..
C     .. Array Arguments ..
      INTEGER ICN(LICN),IP(N),IPERM(N),IW(N,4),LENR(N)
C     ..
C     .. External Subroutines ..
      EXTERNAL MC21BD
C     ..
C     .. Executable Statements ..
      CALL MC21BD(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW(1,1),IW(1,2),
     +            IW(1,3),IW(1,4))
      RETURN
C
      END
      SUBROUTINE MC21BD(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,PR,ARP,CV,OUT)
C   PR(I) IS THE PREVIOUS ROW TO I IN THE DEPTH FIRST SEARCH.
C IT IS USED AS A WORK ARRAY IN THE SORTING ALGORITHM.
C   ELEMENTS (IPERM(I),I) I=1, ... N  ARE NON-ZERO AT THE END OF THE
C ALGORITHM UNLESS N ASSIGNMENTS HAVE NOT BEEN MADE.  IN WHICH CASE
C (IPERM(I),I) WILL BE ZERO FOR N-NUMNZ ENTRIES.
C   CV(I) IS THE MOST RECENT ROW EXTENSION AT WHICH COLUMN I
C WAS VISITED.
C   ARP(I) IS ONE LESS THAN THE NUMBER OF NON-ZEROS IN ROW I
C WHICH HAVE NOT BEEN SCANNED WHEN LOOKING FOR A CHEAP ASSIGNMENT.
C   OUT(I) IS ONE LESS THAN THE NUMBER OF NON-ZEROS IN ROW I
C WHICH HAVE NOT BEEN SCANNED DURING ONE PASS THROUGH THE MAIN LOOP.
C
C   INITIALIZATION OF ARRAYS.
C     .. Scalar Arguments ..
      INTEGER LICN,N,NUMNZ
C     ..
C     .. Array Arguments ..
      INTEGER ARP(N),CV(N),ICN(LICN),IP(N),IPERM(N),LENR(N),OUT(N),PR(N)
C     ..
C     .. Local Scalars ..
      INTEGER I,II,IN1,IN2,IOUTK,J,J1,JORD,K,KK
C     ..
C     .. Executable Statements ..
      DO 10 I = 1,N
        ARP(I) = LENR(I) - 1
        CV(I) = 0
        IPERM(I) = 0
   10 CONTINUE
      NUMNZ = 0
C
C
C   MAIN LOOP.
C   EACH PASS ROUND THIS LOOP EITHER RESULTS IN A NEW ASSIGNMENT
C OR GIVES A ROW WITH NO ASSIGNMENT.
      DO 100 JORD = 1,N
        J = JORD
        PR(J) = -1
        DO 70 K = 1,JORD
C LOOK FOR A CHEAP ASSIGNMENT
          IN1 = ARP(J)
          IF (IN1.LT.0) GO TO 30
          IN2 = IP(J) + LENR(J) - 1
          IN1 = IN2 - IN1
          DO 20 II = IN1,IN2
            I = ICN(II)
            IF (IPERM(I).EQ.0) GO TO 80
   20     CONTINUE
C   NO CHEAP ASSIGNMENT IN ROW.
          ARP(J) = -1
C   BEGIN LOOKING FOR ASSIGNMENT CHAIN STARTING WITH ROW J.
   30     CONTINUE
          OUT(J) = LENR(J) - 1
C INNER LOOP.  EXTENDS CHAIN BY ONE OR BACKTRACKS.
          DO 60 KK = 1,JORD
            IN1 = OUT(J)
            IF (IN1.LT.0) GO TO 50
            IN2 = IP(J) + LENR(J) - 1
            IN1 = IN2 - IN1
C FORWARD SCAN.
            DO 40 II = IN1,IN2
              I = ICN(II)
              IF (CV(I).EQ.JORD) GO TO 40
C   COLUMN I HAS NOT YET BEEN ACCESSED DURING THIS PASS.
              J1 = J
              J = IPERM(I)
              CV(I) = JORD
              PR(J) = J1
              OUT(J1) = IN2 - II - 1
              GO TO 70
C
   40       CONTINUE
C
C   BACKTRACKING STEP.
   50       CONTINUE
            J = PR(J)
            IF (J.EQ.-1) GO TO 100
   60     CONTINUE
C
   70   CONTINUE
C
C   NEW ASSIGNMENT IS MADE.
   80   CONTINUE
        IPERM(I) = J
        ARP(J) = IN2 - II - 1
        NUMNZ = NUMNZ + 1
        DO 90 K = 1,JORD
          J = PR(J)
          IF (J.EQ.-1) GO TO 100
          II = IP(J) + LENR(J) - OUT(J) - 2
          I = ICN(II)
          IPERM(I) = J
   90   CONTINUE
C
  100 CONTINUE
C
C   IF MATRIX IS STRUCTURALLY SINGULAR, WE NOW COMPLETE THE
C PERMUTATION IPERM.
      IF (NUMNZ.EQ.N) RETURN
      DO 110 I = 1,N
        ARP(I) = 0
  110 CONTINUE
      K = 0
      DO 130 I = 1,N
        IF (IPERM(I).NE.0) GO TO 120
        K = K + 1
        OUT(K) = I
        GO TO 130
C
  120   CONTINUE
        J = IPERM(I)
        ARP(J) = I
  130 CONTINUE
      K = 0
      DO 140 I = 1,N
        IF (ARP(I).NE.0) GO TO 140
        K = K + 1
        IOUTK = OUT(K)
        IPERM(IOUTK) = I
  140 CONTINUE
      RETURN
C
      END
* COPYRIGHT (c) 1993 AEA Technology and
* Council for the Central Laboratory of the Research Councils
C Original date March 1993
C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC29AD(M,N,NE,A,IRN,ICN,R,C,W,LP,IFAIL)
      INTEGER M,N,NE
      DOUBLE PRECISION A(NE)
      INTEGER IRN(NE),ICN(NE)
      DOUBLE PRECISION R(M),C(N),W(M*2+N*3)
      INTEGER LP,IFAIL
C M is an integer variable that must be set to the number of rows.
C      It is not altered by the subroutine.
C N is an integer variable that must be set to the number of columns.
C      It is not altered by the subroutine.
C NE is an integer variable that must be set to the number of entries.
C      It is not altered by the subroutine.
C A is an array that holds the values of the entries.
C IRN  is an integer array that must be set to the row indices of the
C      entries. It is not altered by the subroutine.
C ICN  is an integer array that must be set to the column indices of the
C      entries. It is not altered by the subroutine.
C R is an array that need not be set on entry. On return, it holds the
C      logarithms of the row scaling factors.
C C is an array that need not be set on entry. On return, it holds the
C      logarithms of the column scaling factors.
C W is a workarray.
C      W(1:M)  holds row non-zero counts (diagonal matrix M).
C      W(M+1:M+N) holds column non-zero counts (diagonal matrix N).
C      W(M+N+J) holds the logarithm of the column I scaling
C         factor during the iteration, J=1,2,...,N.
C      W(M+N*2+J) holds the 2-iteration change in the logarithm
C         of the column J scaling factor, J=1,2,...,N.
C      W(M+N*3+I) is used to save the average logarithm of
C          the entries of row I, I=1,2,...,M.
C LP must be set to the unit number for messages.
C      It is not altered by the subroutine.
C IFAIL need not be set by the user. On return it has one of the
C     following values:
C     0 successful entry.
C     -1 M < 1 or N < 1.
C     -2 NE < 1.

      INTRINSIC LOG,ABS,MIN

C Constants
      INTEGER MAXIT
      PARAMETER (MAXIT=100)
      DOUBLE PRECISION ONE,SMIN,ZERO
      PARAMETER (ONE=1D0,SMIN=0.1,ZERO=0D0)
C MAXIT is the maximal permitted number of iterations.
C SMIN is used in a convergence test on (residual norm)**2

C Local variables
      INTEGER I,I1,I2,I3,I4,I5,ITER,J,K
      DOUBLE PRECISION E,E1,EM,Q,Q1,QM,S,S1,SM,U,V

C Check M, N and NE.
      IFAIL = 0
      IF (M.LT.1 .OR. N.LT.1) THEN
         IFAIL = -1
         GO TO 220
      ELSE IF (NE.LE.0) THEN
         IFAIL = -2
         GO TO 220
      END IF

C     Partition W
      I1 = 0
      I2 = M
      I3 = M + N
      I4 = M + N*2
      I5 = M + N*3

C     Initialise for accumulation of sums and products.
      DO 10 I = 1,M
         R(I) = ZERO
         W(I1+I) = ZERO
   10 CONTINUE
      DO 20 J = 1,N
         C(J) = ZERO
         W(I2+J) = ZERO
         W(I3+J) = ZERO
         W(I4+J) = ZERO
   20 CONTINUE

C     Count non-zeros in the rows, and compute rhs vectors.
      DO 30 K = 1,NE
         U = ABS(A(K))
         IF (U.EQ.ZERO) GO TO 30
         I = IRN(K)
         J = ICN(K)
         IF (MIN(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 30
         U = LOG(U)
         W(I1+I) = W(I1+I) + ONE
         W(I2+J) = W(I2+J) + ONE
         R(I) = R(I) + U
         W(I3+J) = W(I3+J) + U
   30 CONTINUE
C
C     Divide rhs by diag matrices.
      DO 40 I = 1,M
         IF (W(I1+I).EQ.ZERO) W(I1+I) = ONE
         R(I) = R(I)/W(I1+I)
C     Save R(I) for use at end.
         W(I5+I) = R(I)
   40 CONTINUE
      DO 50 J = 1,N
         IF (W(I2+J).EQ.ZERO) W(I2+J) = ONE
         W(I3+J) = W(I3+J)/W(I2+J)
   50 CONTINUE
      SM = SMIN*NE

C     Sweep to compute initial residual vector
      DO 60 K = 1,NE
         IF (A(K).EQ.ZERO) GO TO 60
         I = IRN(K)
         J = ICN(K)
         IF (MIN(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 60
         R(I) = R(I) - W(I3+J)/W(I1+I)
   60 CONTINUE
C
C     Initialise iteration
      E = ZERO
      Q = ONE
      S = ZERO
      DO 70 I = 1,M
         S = S + W(I1+I)*R(I)**2
   70 CONTINUE
      IF (S.LE.SM) GO TO 160

C     Iteration loop
      DO 150 ITER = 1,MAXIT
C    Sweep through matrix to update residual vector
         DO 80 K = 1,NE
            IF (A(K).EQ.ZERO) GO TO 80
            J = ICN(K)
            I = IRN(K)
            IF (MIN(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 80
            C(J) = C(J) + R(I)
   80    CONTINUE
         S1 = S
         S = ZERO
         DO 90 J = 1,N
            V = -C(J)/Q
            C(J) = V/W(I2+J)
            S = S + V*C(J)
   90    CONTINUE
         E1 = E
         E = Q*S/S1
         Q = ONE - E
C      write(*,'(a,i3,a,f12.4)')' Iteration',ITER,' S =',S
         IF (S.LE.SM) E = ZERO
C     Update residual.
         DO 100 I = 1,M
            R(I) = R(I)*E*W(I1+I)
  100    CONTINUE
         IF (S.LE.SM) GO TO 180
         EM = E*E1
C    Sweep through matrix to update residual vector
         DO 110 K = 1,NE
            IF (A(K).EQ.ZERO) GO TO 110
            I = IRN(K)
            J = ICN(K)
            IF (MIN(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 110
            R(I) = R(I) + C(J)
  110    CONTINUE
         S1 = S
         S = ZERO
         DO 120 I = 1,M
            V = -R(I)/Q
            R(I) = V/W(I1+I)
            S = S + V*R(I)
  120    CONTINUE
         E1 = E
         E = Q*S/S1
         Q1 = Q
         Q = ONE - E
C     Special fixup for last iteration.
         IF (S.LE.SM) Q = ONE
C     Update col. scaling powers
         QM = Q*Q1
         DO 130 J = 1,N
            W(I4+J) = (EM*W(I4+J)+C(J))/QM
            W(I3+J) = W(I3+J) + W(I4+J)
  130    CONTINUE
C      write(*,'(a,i3,a,f12.4)')' Iteration',ITER,' S =',S
         IF (S.LE.SM) GO TO 160
C     UPDATE RESIDUAL.
         DO 140 J = 1,N
            C(J) = C(J)*E*W(I2+J)
  140    CONTINUE
  150 CONTINUE
  160 DO 170 I = 1,M
         R(I) = R(I)*W(I1+I)
  170 CONTINUE
C
C     Sweep through matrix to prepare to get row scaling powers
  180 DO 190 K = 1,NE
         IF (A(K).EQ.ZERO) GO TO 190
         I = IRN(K)
         J = ICN(K)
         IF (MIN(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 190
         R(I) = R(I) + W(I3+J)
  190 CONTINUE
C
C     Final conversion to output values.
      DO 200 I = 1,M
         R(I) = R(I)/W(I1+I) - W(I5+I)
  200 CONTINUE
      DO 210 J = 1,N
         C(J) = -W(I3+J)
  210 CONTINUE
      RETURN

C Error returns
  220 IF (LP.GT.0) WRITE (LP,'(/A/A,I3)')
     +    ' **** Error return from MC29AD ****',' IFAIL =',IFAIL

      END
* COPYRIGHT (c) 1976 AEA Technology
* Original date 21 Jan 1993
C       Toolpack tool decs employed.
C	Double version of MC13D (name change only)
C 10 August 2001 DOs terminated with CONTINUE
C 13/3/02 Cosmetic changes applied to reduce single/double differences
C
C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC13DD(N,ICN,LICN,IP,LENR,IOR,IB,NUM,IW)
C     .. Scalar Arguments ..
      INTEGER LICN,N,NUM
C     ..
C     .. Array Arguments ..
      INTEGER IB(N),ICN(LICN),IOR(N),IP(N),IW(N,3),LENR(N)
C     ..
C     .. External Subroutines ..
      EXTERNAL MC13ED
C     ..
C     .. Executable Statements ..
      CALL MC13ED(N,ICN,LICN,IP,LENR,IOR,IB,NUM,IW(1,1),IW(1,2),IW(1,3))
      RETURN

      END
      SUBROUTINE MC13ED(N,ICN,LICN,IP,LENR,ARP,IB,NUM,LOWL,NUMB,PREV)
C
C ARP(I) IS ONE LESS THAN THE NUMBER OF UNSEARCHED EDGES LEAVING
C     NODE I.  AT THE END OF THE ALGORITHM IT IS SET TO A
C     PERMUTATION WHICH PUTS THE MATRIX IN BLOCK LOWER
C     TRIANGULAR FORM.
C IB(I) IS THE POSITION IN THE ORDERING OF THE START OF THE ITH
C     BLOCK.  IB(N+1-I) HOLDS THE NODE NUMBER OF THE ITH NODE
C     ON THE STACK.
C LOWL(I) IS THE SMALLEST STACK POSITION OF ANY NODE TO WHICH A PATH
C     FROM NODE I HAS BEEN FOUND.  IT IS SET TO N+1 WHEN NODE I
C     IS REMOVED FROM THE STACK.
C NUMB(I) IS THE POSITION OF NODE I IN THE STACK IF IT IS ON
C     IT, IS THE PERMUTED ORDER OF NODE I FOR THOSE NODES
C     WHOSE FINAL POSITION HAS BEEN FOUND AND IS OTHERWISE ZERO.
C PREV(I) IS THE NODE AT THE END OF THE PATH WHEN NODE I WAS
C     PLACED ON THE STACK.
C
C
C   ICNT IS THE NUMBER OF NODES WHOSE POSITIONS IN FINAL ORDERING HAVE
C     BEEN FOUND.
C     .. Scalar Arguments ..
      INTEGER LICN,N,NUM
C     ..
C     .. Array Arguments ..
      INTEGER ARP(N),IB(N),ICN(LICN),IP(N),LENR(N),LOWL(N),NUMB(N),
     +        PREV(N)
C     ..
C     .. Local Scalars ..
      INTEGER DUMMY,I,I1,I2,ICNT,II,ISN,IST,IST1,IV,IW,J,K,LCNT,NNM1,STP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MIN
C     ..
C     .. Executable Statements ..
      ICNT = 0
C NUM IS THE NUMBER OF BLOCKS THAT HAVE BEEN FOUND.
      NUM = 0
      NNM1 = N + N - 1
C
C INITIALIZATION OF ARRAYS.
      DO 20 J = 1,N
        NUMB(J) = 0
        ARP(J) = LENR(J) - 1
   20 CONTINUE
C
C
      DO 120 ISN = 1,N
C LOOK FOR A STARTING NODE
        IF (NUMB(ISN).NE.0) GO TO 120
        IV = ISN
C IST IS THE NUMBER OF NODES ON THE STACK ... IT IS THE STACK POINTER.
        IST = 1
C PUT NODE IV AT BEGINNING OF STACK.
        LOWL(IV) = 1
        NUMB(IV) = 1
        IB(N) = IV
C
C THE BODY OF THIS LOOP PUTS A NEW NODE ON THE STACK OR BACKTRACKS.
        DO 110 DUMMY = 1,NNM1
          I1 = ARP(IV)
C HAVE ALL EDGES LEAVING NODE IV BEEN SEARCHED.
          IF (I1.LT.0) GO TO 60
          I2 = IP(IV) + LENR(IV) - 1
          I1 = I2 - I1
C
C LOOK AT EDGES LEAVING NODE IV UNTIL ONE ENTERS A NEW NODE OR
C     ALL EDGES ARE EXHAUSTED.
          DO 50 II = I1,I2
            IW = ICN(II)
C HAS NODE IW BEEN ON STACK ALREADY.
            IF (NUMB(IW).EQ.0) GO TO 100
C UPDATE VALUE OF LOWL(IV) IF NECESSARY.
            LOWL(IV) = MIN(LOWL(IV),LOWL(IW))
   50     CONTINUE
C
C THERE ARE NO MORE EDGES LEAVING NODE IV.
          ARP(IV) = -1
C IS NODE IV THE ROOT OF A BLOCK.
   60     IF (LOWL(IV).LT.NUMB(IV)) GO TO 90
C
C ORDER NODES IN A BLOCK.
          NUM = NUM + 1
          IST1 = N + 1 - IST
          LCNT = ICNT + 1
C PEEL BLOCK OFF THE TOP OF THE STACK STARTING AT THE TOP AND
C     WORKING DOWN TO THE ROOT OF THE BLOCK.
          DO 70 STP = IST1,N
            IW = IB(STP)
            LOWL(IW) = N + 1
            ICNT = ICNT + 1
            NUMB(IW) = ICNT
            IF (IW.EQ.IV) GO TO 80
   70     CONTINUE
   80     IST = N - STP
          IB(NUM) = LCNT
C ARE THERE ANY NODES LEFT ON THE STACK.
          IF (IST.NE.0) GO TO 90
C HAVE ALL THE NODES BEEN ORDERED.
          IF (ICNT.LT.N) GO TO 120
          GO TO 130
C
C BACKTRACK TO PREVIOUS NODE ON PATH.
   90     IW = IV
          IV = PREV(IV)
C UPDATE VALUE OF LOWL(IV) IF NECESSARY.
          LOWL(IV) = MIN(LOWL(IV),LOWL(IW))
          GO TO 110
C
C PUT NEW NODE ON THE STACK.
  100     ARP(IV) = I2 - II - 1
          PREV(IW) = IV
          IV = IW
          IST = IST + 1
          LOWL(IV) = IST
          NUMB(IV) = IST
          K = N + 1 - IST
          IB(K) = IV
  110   CONTINUE
C
  120 CONTINUE
C
C
C PUT PERMUTATION IN THE REQUIRED FORM.
  130 DO 140 I = 1,N
        II = NUMB(I)
        ARP(II) = I
  140 CONTINUE
      RETURN

      END
C COPYRIGHT (c) 1994 Council for the Central Laboratory
*                    of the Research Councils
C Original date 30 November 1995

C JAS   : 21 May 1998
C         LGUARD set to NBF if IFLAG = -6 returned.
C JAS   : 30 November 1998
C         WRITE (LP,FMT=9070) LGUARD,NBF
C  changed to
C         IF (LP.GT.0) WRITE (LP,FMT=9070) LGUARD,NBF
C JAS   : 19 April 1999
C         The workspace in MA52C (IND = 1 or 3) could be too small if
C         MA42B/BD was called with a large number of right-hand sides.
C         Note: bug also fixed in MA42, but because INFO is not passed
C         from MA42 to MA52, we cannot use the actual max frontsizes
C         used by MA42B/BD, but only the bounds NFRONT set by the user.
C         In MA52A/AD changed LGUARD = NBF to IW(1) = NBF (error flag
C         -6) since LGUARD is dim. of IGUARD and cannot be changed.
C JAS   : 14 Oct. 1999
C         Added extra entry MA52F/FD that is the same as MA52B/BD
C         except it saves ROW and COL indices left in front.
C         Thus pivoting need not be restricted to diagonal.
C JAS   : 12 Jan. 2000
C         Can also cope with frontal matrix not square
C         Variables NFAVR and LFVAR replaced by arrays of length 2.
C JAS   : 17 Jan. 2000. No longer initialise X to zero in MA52C/CD
C         if IND=2.
C JAS   : 7 Aug. 2001
C         Added extra entry MA52K/KD that is the same as MA52F/FD
C         except it saves only the nonzero entries left in front
C         (so the frontal matrix is stored in sparse format).
C JAS   : 4 Dec. 2001. Bug fixed when MA52L/LD retunrs IFLAG < 0
C         (also MA52D/DD and MA52G/GD)

C 12th July 2004 Version 1.0.0. Version numbering added.
C 28 February 2008. Version 1.0.1. Comments flowed to column 72.

      SUBROUTINE MA52AD(ICALL,NVAR,IVAR,NELT,NDOMN,IDOMN,NGUARD,LGUARD,
     +                  IGUARD,LIW,IW,LP,IFLAG)
C
C This subroutine generates the guard elements.
C
C ICALL  - must be set by the user to the number of the call
C          to MA52A (must be set to 1 on the first call,
C          to 2 on the second call, and so on).
c          Unchanged on exit.
C NVAR   - must be set by the user to the number of variables
C          in the element. Unchanged on exit. NVAR.ge.1
C IVAR   - integer array of length NVAR. Must be set by the user
C          to hold a list of the variables in the element.
c          Unchanged on exit.
C NELT   - must be set by the user to the number of finite-elements
C          in the whole domain. Unchanged on exit. NELT.gt.2
C IDOMN  - must be set by the user to the index of the subdomain
C          to which the element belongs.
C          Unchanged on exit. IDOMN.ge.1 and IDOMN.le.NDOMN.
C NDOMN  - must be set by the user to the number of  subdomains
C          Unchanged on exit. NDOMN.gt.1
C NGUARD - integer array of length NDOMN. Not set on entry.
C          On final exit, NGUARD(IDOMN) holds the number of variables
C          in the guard element for subdomain IDOMN.
C LGUARD - must be set to the first dimension of the array IGUARD.
C IGUARD - integer array of dimensions LGUARD,NDOMN.
C          Not set on entry. On final exit,  IGUARD(I,IDOMN),
C          I=1,...,NGUARD(IDOMN) is a list of the variables in the
C          guard element for subdomain IDOMN.
C LIW    - must be set by the user to the length of the array IW.
C          LIW must be at least as large as the largest integer
C          used to index a variable in the finite element mesh.
C IW     - integer array of length LIW. Used as workspace.
C LP     - must be set to the stream number for error messages.
C          Unchanged on exit.
C IFLAG  - error flag. Need not be set on entry. A negative value
C          implies a fatal error. +1 is returned if LGUARD is
C          insufficient but allows computation to continue.
C
C     .. Scalar Arguments ..
      INTEGER ICALL,IDOMN,IFLAG,LGUARD,LIW,LP,NDOMN,NELT,NVAR
C     ..
C     .. Array Arguments ..
      INTEGER IGUARD(LGUARD,NDOMN),IVAR(NVAR),IW(LIW),NGUARD(NDOMN)
C     ..
C     .. Local Scalars ..
C
      INTEGER I,J,JDOMN,JSTOP,K,NBF
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX
C     ..
C If ICALL=1 then first call...check input data and
C perform initialisations.
      IF (ICALL.EQ.1) THEN
         IFLAG = 0
         IF (NELT.LE.1) GO TO 90
         IF (NDOMN.LE.1) GO TO 100
         DO 10 I = 1,NDOMN
            NGUARD(I) = 0
   10    CONTINUE
         DO 20 I = 1,LIW
            IW(I) = 0
   20    CONTINUE
      END IF
C
      IF (NVAR.LT.1) GO TO 110
      IF (IDOMN.LT.1 .OR. IDOMN.GT.NDOMN) GO TO 130
C
C Loop over variables in the element
      IF (ICALL.EQ.1) THEN
         DO 30 K = 1,NVAR
            J = IVAR(K)
C Check J.le.LIW
            IF (J.GT.LIW .OR. J.LT.1) GO TO 120
            IW(J) = IDOMN
   30    CONTINUE
      ELSE
         DO 40 K = 1,NVAR
            J = IVAR(K)
C Check J.le.LIW
            IF (J.GT.LIW .OR. J.LT.1) GO TO 120
            JDOMN = IW(J)
            IF (JDOMN.NE.0 .AND. JDOMN.NE.IDOMN) THEN
C Variable J has already appeared in another subdomain and
C so must be on an interface.
               NGUARD(JDOMN) = NGUARD(JDOMN) + 1
               NGUARD(IDOMN) = NGUARD(IDOMN) + 1
               IF (NGUARD(JDOMN).GT.LGUARD .OR.
     +             NGUARD(IDOMN).GT.LGUARD) THEN
                  IF (IFLAG.EQ.0) THEN
                     IFLAG = 1
                     IF (LP.GT.0) WRITE (LP,FMT=9000) ICALL,IFLAG
                     IF (LP.GT.0) WRITE (LP,FMT=9060)
                  END IF
               ELSE
                  IGUARD(NGUARD(JDOMN),JDOMN) = J
                  IGUARD(NGUARD(IDOMN),IDOMN) = J
               END IF
            END IF
            IW(J) = IDOMN
   40    CONTINUE
      END IF
C
      IF (ICALL.EQ.NELT) THEN
         IF (IFLAG.EQ.0) THEN
C We need to remove duplicates from the guard elements.
            DO 50 I = 1,LIW
               IW(I) = 0
   50       CONTINUE
            DO 70 JDOMN = 1,NDOMN
               JSTOP = NGUARD(JDOMN)
               NBF = 0
               DO 60 J = 1,JSTOP
                  I = IGUARD(J,JDOMN)
                  IF (IW(I).LT.JDOMN) THEN
                     IW(I) = JDOMN
                     NBF = NBF + 1
                     IGUARD(NBF,JDOMN) = I
                  END IF
   60          CONTINUE
               NGUARD(JDOMN) = NBF
   70       CONTINUE
         ELSE
            NBF = NGUARD(1)
            DO 80 JDOMN = 2,NDOMN
               NBF = MAX(NBF,NGUARD(JDOMN))
   80       CONTINUE
            IFLAG = -6
            IF (LP.GT.0) WRITE (LP,FMT=9000) ICALL,IFLAG
            IF (LP.GT.0) WRITE (LP,FMT=9070) LGUARD,NBF
            IW(1) = NBF
         END IF
      END IF
      GO TO 140
C
   90 IFLAG = -1
      IF (LP.GT.0) WRITE (LP,FMT=9000) ICALL,IFLAG
      IF (LP.GT.0) WRITE (LP,FMT=9010) NELT
      GO TO 140
  100 IFLAG = -2
      IF (LP.GT.0) WRITE (LP,FMT=9000) ICALL,IFLAG
      IF (LP.GT.0) WRITE (LP,FMT=9020) NDOMN
      GO TO 140
  110 IFLAG = -3
      IF (LP.GT.0) WRITE (LP,FMT=9000) ICALL,IFLAG
      IF (LP.GT.0) WRITE (LP,FMT=9030) NVAR
      GO TO 140
  120 IFLAG = -4
      IF (LP.GT.0) WRITE (LP,FMT=9000) ICALL,IFLAG
      IF (LP.GT.0) WRITE (LP,FMT=9040) K,J
      GO TO 140
  130 IFLAG = -5
      IF (LP.GT.0) WRITE (LP,FMT=9000) ICALL,IFLAG
      IF (LP.GT.0) WRITE (LP,FMT=9050) IDOMN
  140 CONTINUE
      RETURN
 9000 FORMAT (' Error return from MA52A/AD on call ',I6,'. IFLAG = ',I3)
 9010 FORMAT (' Value of NELT out of range. NELT = ',I6)
 9020 FORMAT (' Value of NDOMN out of range. NDOMN = ',I6)
 9030 FORMAT (' Value of NVAR out of range. NVAR = ',I6)
 9040 FORMAT (' Variable  ',I8,' has value  ',I8)
 9050 FORMAT (' Value of IDOMN out of range. IDOMN = ',I6)
 9060 FORMAT (' LGUARD is too small. Complete calling sequence to ',
     +       'determine sufficient value.')
 9070 FORMAT (' Increase LGUARD from ',I8,' to at least ',I8)
      END
C**********************************************************************
      SUBROUTINE MA52BD(LRHS,NRHS,NDF,LAST,LX,X,LW,W,LIW,IW,NFVAR,IFVAR,
     +                  FVAR,FRHS,LFVAR,ISAVE,LP,IFLAG)
C
C *** Copyright Rutherford Appleton Laboratory  November 1993 ***
C *** Any problems contact Jennifer A. Scott at Atlas Centre,
C     Rutherford Appleton Laboratory ***
C
C  This routine dumps the contents of the  buffers into the direct
C  access files used by MA42 and writes what remains in the
C  frontal matrix (and corresponding frontal right-hand sides)
C  into the form of an element matrix and element right-hand side.
C  This routine can only be used after MA42B/BD has been called
C  for at least one element but has not been called for all
C  all the elements in the problem.
C
C     Argument list. * indicates the argument is changed by the routine.
C
C  LRHS   - Integer variable. Must be unchanged since MA42B/BD.
C  NRHS   - Integer variable. Must be unchanged since MA42B/BD.
C  NDF    - Integer variable. Length of array LAST.
C *LAST   - Integer array of length NDF.  On exit LAST is reset
C           to its state on first entry to MA42B/BD.
C  LX     - Integer variable. First dimension of  array X.
C *X      - Real (DP) array of dimensions LX by LRHS. On exit,
C           X(I,I)=0, I=1,2,...,NDF, J=1,2,...,NRHS.
C  LW     - Integer variable. Length of array W.
C  W      - Real (DP) array of length LW.
C  LIW    - Integer variable. Length of array IW.
C  IW     - Integer array of length LIW.
C  LRHS, LAST, NDF, LW, W, LIW, IW must all be unchanged from the
C           last call to MA42B/BD.
C *NFVAR  - Integer variable. Not set on entry.
C           On exit, holds the number of variables
C           left in the front after last call to MA42B/BD
C           (=ISAVE(40)).
C *IFVAR  - Integer array of dimensions (LFVAR,LFVAR).
C           Not set on entry. On exit holds the indices
C           of the variables left in the front after last
C           call to MA42B/B.
C *FVAR   - Real (dp) array of dimensions (LFVAR,LFVAR).
C           Not set on entry. On exit holds values for the
C           variables left in the front after last
C           call to MA42B/B.
C *FRHS   - Real (dp) array of dimensions (LFVAR,LRHS).
C           Not set on entry. On exit holds right-hand sides
C           for the variables left in the front after last
C           call to MA42B/B.
C  LFVAR  - Integer variable. Defines dimension of IFVAR, FVAR,
C           and FRHS. Must be at least ISAVE(40).
C *ISAVE  - Integer array of length 45. Holds parameters which
C           must be preserved between calls to MA42B/BD and between
C           calls to other MA42 routines.
C           ISAVE(1),...,ISAVE(15) were set by MA42I/ID and MA42P/PD.
C           If MA42P/PD was called, ISAVE(6+I) was set to
C           LENBUF(I), I=1,2,3 on the first entry to MA42B/BD.
C           ISAVE(16) used by MA42B/BD if MA42J/JD called.
C           ISAVE(17) is set to 0 by MA42I/ID.
C           On exit from intermediate calls to MA42B/BD, ISAVE(17)
C           holds the number of elements input so far and on exit
C           from the final call to MA42B/BD is set to 0.
C           On return from MA52B/BD,
C           ISAVE(17) = 0 (this allows MA42B/BD
C           to be rerun without rerunning MA42A/AD).
C           ISAVE(18) set by MA42A/AD to number of elements in problem
C           and is unchanged.
C           A non-zero value of ISAVE(19) or ISAVE(20) indicates
C           the occurrence of the non-fatal error 5 or 6
C           respectively, in MA42B/BD.
C           ISAVE(21), ISAVE(22), ISAVE(23) and ISAVE(24)
C           hold the values of variables NMAXE, NRHS,
C           NFRONT(1) and NFRONT(2), resp. as input in the first call
C           to MA42B/BD.
C           ISAVE(25)-ISAVE(29) are used to divide real workspace.
C           ISAVE(30)-set by MA42A/AD to largest integer used to index
C           a variable.
C           ISAVE(31)-flag to indicate if MA42P/PD has been called.
C           ISAVE(32)-ISAVE(37) are used to divide integer workspace.
C           ISAVE(38)- set by MA42B/BD to hold a copy of the
C           error flag INFO(1) if error encountered by MA42B/BD
C           and holds copy of ICNTL(7) from MA42I/ID o.w.
C           On exit, ISAVE(38) set to 3003 to indicate MA52B/BD
C           has been called.
C           ISAVE(39) is set to 0 by MA42I/ID. On exit from
C           intermediate calls to MA42B/BD, ISAVE(39) holds current row
C           front size KFRNT.
C           On return from MA52B/BD, ISAVE(39)=0.
C           ISAVE(40) is set to 0 by MA42I/ID. On exit from
C           intermediate calls to MA42B/BD, ISAVE(40) holds current
C           column front size LFRNT.
C           On return from MA52B/BD, ISAVE(40)=0.
C           ISAVE(41)-ISAVE(45) used to preserve values between calls
C           to MA42F/FD.
C     (ISAVE is altered in the way it would have been if the series
C      of calls to MA42B/BD had been completed, except that the number
C      of records used in the direct access data sets may be different
C      so MKEY and IREC may differ... ISAVE(3+I), ISAVE(9+I)).
C     This array must be unchanged since calls to MA42B/BD.
C   LP    - Integer variable. Set by user to stream number for error
C           messages. Printing suppressed if LP.LE.0
C *IFLAG  - Integer variable. Not set on entry. On successful
C           exit, set to 0. Negative value indicates error.
C
C   Local variables
C
C   FA     - Integer variable. Used in dividing real workspace.
C   FARHS  - Integer variable. Used in dividing real workspace.
C   I      - Integer variable. Do loop variable.
C   J      - Integer variable. Do loop variable.
C   JBUFL  - Integer variable. Holds max(1,ISIZE(2)). Used in dividing
C            real workspace.
C   JBUFR  - Integer variable. Holds 1+SIZE(2). Used in dividing
C            real workspace.
C   JBUFU  - Integer variable. Holds ISIZE(1) (size of U-buffer).
C            Used in dividing real workspace.
C   KPVLNK - Integer variable. Used in dividing integer workspace.
C   LHED   - Integer variable. Used in dividing integer workspace.
C   IFILE  - Integer array of length 3. IFILE(I) is equivalent to
C            ISAVE(I) (I = 1,2,3).
C   IREC   - Integer array of length 3. IREC(I) is equivalent to
C            ISAVE(3+I) (I = 1,2,3).
C   ISIZE  - Integer array of length 3. ISIZE(I) is equivalent to
C            ISAVE(6+I) (I = 1,2,3).
C   MKEY   - Integer array of length 3. MKEY(I) is equivalent to
C            ISAVE(9+I) (I = 1,2,3).
C   NUMBLK - Integer array of length 3. NUMBLK(I) is equivalent to
C            ISAVE(12+I) (I = 1,2,3).
C
C     .. Scalar Arguments ..
      INTEGER IFLAG,LFVAR,LIW,LP,LRHS,LW,LX,NDF,NFVAR,NRHS
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION FRHS(LFVAR,LRHS),FVAR(LFVAR,LFVAR),W(LW),
     +                 X(LX,LRHS)
      INTEGER IFVAR(LFVAR),ISAVE(45),IW(LIW),LAST(NDF)
C     ..
C     .. Local Scalars ..
      INTEGER FA,FARHS,I,J,JBUFL,JBUFR,JBUFU,KPVLNK,LHED,NFRONT
C     ..
C     .. Local Arrays ..
      INTEGER IFILE(3),IREC(3),ISIZE(3),MKEY(3),NUMBLK(3)
C     ..
C     .. External Subroutines ..
      EXTERNAL MA52DD,MA52ED
C     ..
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
      IFLAG = 0
C Immediate return if MA42B/BD has issued an error.
      IF (ISAVE(38).LT.0 .OR. ISAVE(38).EQ.4 .OR. ISAVE(38).EQ.5 .OR.
     +    ISAVE(38).EQ.6) GO TO 50
C Check diagonal pivoting was used by MA42B/BD
      IF (ABS(ISAVE(38)).NE.1001) GO TO 100
C Check that MA42B/BD was not called for all the elements
C in the problem. If it was, then immediate return.
      IF (ISAVE(17).EQ.0) GO TO 60
      DO 10 I = 1,3
         IFILE(I) = ISAVE(I)
         IREC(I) = ISAVE(3+I)
         ISIZE(I) = ISAVE(6+I)
         MKEY(I) = ISAVE(9+I)
         NUMBLK(I) = ISAVE(12+I)
   10 CONTINUE
C Test validity of input parameters.
      NFVAR = ISAVE(40)
      IF (NRHS.NE.ISAVE(22)) GO TO 80
      IF (LFVAR.LT.NFVAR) GO TO 70
      IF (NRHS.GT.0 .AND. LX.LT.NDF) GO TO 90
C
C Perform dump of buffer contents to direct access files
      NFRONT = ISAVE(23)
      JBUFU = ISAVE(25)
      JBUFL = ISAVE(26)
      JBUFR = ISAVE(27)
      FA = ISAVE(28)
      FARHS = ISAVE(29)
      LHED = ISAVE(32)
C
      CALL MA52DD(NFRONT,NRHS,LRHS,W,JBUFL,W(JBUFR),JBUFU,IW,ISIZE(3),
     +            W(FA),W(FARHS),IW(LHED),IFILE,IREC,ISIZE,MKEY,NUMBLK,
     +            IFVAR,FVAR,FRHS,LFVAR,NFVAR,ISAVE,LP,IFLAG)
      IF (IFLAG.LT.0) GO TO 110
C Preserve IREC and MKEY.
      DO 20 I = 1,3
         ISAVE(3+I) = IREC(I)
         ISAVE(9+I) = MKEY(I)
   20 CONTINUE
C
C Reset ISAVE(17), ISAVE(39), ISAVE(40) to 0
C and reset ISAVE(16) to 2 (so MA42A/AD and MA42B/BD
C may be recalled without recalling MA42I/ID).
C
      ISAVE(17) = 0
      ISAVE(39) = 0
      ISAVE(40) = 0
      ISAVE(16) = 2
C Reset LAST
      KPVLNK = ISAVE(35)
      CALL MA52ED(LAST,NDF,NFVAR,IW(LHED),IW(KPVLNK))
C Set X
      DO 40 J = 1,NRHS
         DO 30 I = 1,NDF
            X(I,J) = ZERO
   30    CONTINUE
   40 CONTINUE
C Set ISAVE(38) to indicate call to MA52B/BD complete.
      ISAVE(38) = 3003
      GO TO 110
C
C *** Error returns ***
   50 IFLAG = -1
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) IFLAG
         WRITE (LP,FMT=9040)
      END IF
      GO TO 110
   60 IFLAG = -2
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) IFLAG
         WRITE (LP,FMT=9050)
      END IF
      GO TO 110
   70 IFLAG = -3
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) IFLAG
         WRITE (LP,FMT=9010) LFVAR,NFVAR
      END IF
      GO TO 110
   80 IFLAG = -4
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) IFLAG
         WRITE (LP,FMT=9020)
      END IF
      GO TO 110
   90 IFLAG = -5
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) IFLAG
         WRITE (LP,FMT=9030) NDF
      END IF
      GO TO 110
  100 IFLAG = -8
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) IFLAG
         WRITE (LP,FMT=9060) ISAVE(38)
      END IF
C
  110 RETURN
 9000 FORMAT (' ***** Error return from MA52B/BD *****  IFLAG = ',I3)
 9010 FORMAT (7X,'Length of arrays too small.',/7X,'Increase LFVAR ',
     +       'from ',I8,' to ',I8)
 9020 FORMAT (7X,'NRHS has been changed since the last call to ',
     +       'MA42B/BD.')
 9030 FORMAT (7X,'LX too small. Must be at least ',I8)
 9040 FORMAT (7X,'An error was issued by MA42B/BD.')
 9050 FORMAT (7X,'MA42B/BD was called for all the elements in the ',
     +       'problem.')
 9060 FORMAT (7X,'MA42B/BD was called with ICNTL(7) = ',I3,/7X,'To use',
     +       ' MA52B/BD the user must call MA42B/BD with |ICNTL(7)| =',
     +       ' 1001')
      END
C**********************************************************************
      SUBROUTINE MA52CD(IND,NRHS,LX,B,X,LW,W,LIW,IW,ISAVE,LP,IFLAG)
C
C *** Copyright Rutherford Appleton Laboratory  November 1993 ***
C *** Any problems contact Jennifer A. Scott at Atlas Centre,
C     Rutherford Appleton Laboratory ***
C
C This subroutine is part fo the domain decomposition package
C and is used in solving AX = B
C It calls subroutine MA42D/DD to perform the back-substitution.
C or it calls subroutine MA42E/ED to perform the forward-substitution.
C
C     Argument list. * indicates the argument is changed by the routine.
C
C  IND    - Integer variable.  Indicates back or forward
C           substitution.
C           IND = 1 back sub., solving at same time as factorization.
C           IND = 2 forward sub., further right-hand sides
C           IND = 3 back sub., further right-hand sides
C  NRHS   - Integer variable.  Number of right hand sides.
C  LX     - Integer variable.  Leading dimension of arrays B, X.
C           Must be at least as large as NDF as output from MA42A/AD.
C *B      - Real (DP) array of dimensions LX by NRHS. On entry,
C           if IND.GT.1,  must be set to hold the components
C           of right-hand side for subdomain.
C           On successful exit B(I,J) will hold the Ith component of
C           the partial solution to system J.
C *X      - Real (DP) array of dimensions LX by NRHS. On entry,
C           must be set to hold the components of solution
C           found on combining regions (other entries must be zero).
C           On successful exit X(I,J) will hold the Ith component of
C           the solution to system J.
C  LW     - Integer variable. Length of array W.
C *W      - Real (DP) array of length LW. If direct access
C           files not used, the first LENBUF(1)+LENBUF(2)
C           entries of W must be unchanged since
C           the last call to MA42B/BD on subdomain
C           and these entries are unchanged.
C           Otherwise W is used as workspace.
C  LIW -    Integer variable. Length of array IW.
C *IW     - Integer array of length LIW. If direct access
C           files not used, IW must be unchanged since
C           the last call to MA42B/BD on subdomain and is unchanged.
C           Otherwise IW is used as workspace.
C *LP     - Integer which holds the stream number for error
C           messages.
C  ISAVE  - Integer array of length 45. Holds parameters which
C           must be preserved since call to MA52B/BD on subdomain.
C *IFLAG  - Integer. Error flag.
C
C  Local variables
C
C   IFILE  - Integer array of length 3. IFILE(I) is equivalent to
C            ISAVE(I) (I = 1,2,3).
C   IREC   - Integer array of length 3. IREC(I) is equivalent to
C            ISAVE(3+I) (I = 1,2,3).
C   ISIZE  - Integer array of length 3. ISIZE(I) is equivalent to
C            ISAVE(6+I) (I = 1,2,3).
C   MKEY   - Integer array of length 3. MKEY(I) is equivalent to
C            ISAVE(9+I) (I = 1,2,3).
C   DIMIBF - Integer variable. Used to hold length of integer
C            workspace required by routine.
C   I      - Integer variable. Do loop variable.
C   J      - Integer variable. Do loop variable.
C   L      - Integer variable. Do loop variable.
C L1,L2,L3 _ Integer variables. Equivalent to ISIZE(1),ISIZE(2),
C            ISIZE(3), respectively.
C   LBUFR  - Integer variable. Used to hold length of buffer.
C   LLW    - Integer variable. Used to hold length of real
C            workspace required by routine.
C   NFRONT  - Integer variable. Used to hold max front width
C            (= ISAVE(23)).
C
C     .. Scalar Arguments ..
      INTEGER IFLAG,IND,LIW,LP,LW,LX,NRHS
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION B(LX,NRHS),W(LW),X(LX,NRHS)
      INTEGER ISAVE(45),IW(LIW)
C     ..
C     .. Local Arrays ..
      INTEGER IFILE(3),IREC(3),ISIZE(3),MKEY(3)
C     ..
C     .. Local Scalars ..
      INTEGER DIMIBF,I,IS,JFLAG,L1,L2,L3,LBUFR,LLW,MFRONT,NFRONT,NRHSB
C     ..
C     .. External Subroutines ..
      EXTERNAL MA42DD,MA42ED
C     ..
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C     ..
      IFLAG = 0
      JFLAG = 0
C Recall the direct access files for the subdomain problem
      DO 10 I = 1,3
         IFILE(I) = ISAVE(I)
         IREC(I) = ISAVE(3+I)
         ISIZE(I) = ISAVE(6+I)
         MKEY(I) = ISAVE(9+I)
   10 CONTINUE
C Perform checks
      IF (IND.NE.1 .AND. IND.NE.2 .AND. IND.NE.3) GO TO 110
C Check MA52B/BD was previously called.
      IF (ISAVE(38).NE.3003) GO TO 40
C Check NRHS.
      IF (NRHS.LE.0) GO TO 80
      IF (IND.EQ.1 .AND. NRHS.NE.ISAVE(22)) GO TO 80
      NRHSB = ISAVE(22)
C ISAVE(30) holds copy of NDF for subdomain
C as output from MA42A/AD.
      IF (LX.LT.ISAVE(30)) GO TO 50
C      IF (IND.EQ.2) THEN
C Initialise the solution vector X.
C         DO 30 J = 1,NRHS
C            DO 20 I = 1,ISAVE(30)
C               X(I,J) = ZERO
C   20       CONTINUE
C   30    CONTINUE
C      END IF
C If IND.NE.1 jump to error return if MA42P/PD was
C called with ISTRM(2) = 0 (L-factor was not stored).
C Also jump to error return if MA42B/BD was called with LENBUF(2) = 0.
      IF (IND.NE.1 .AND. ISIZE(2).EQ.0) GO TO 100
      L1 = ISIZE(1)
      L2 = ISIZE(2)
      L3 = ISIZE(3)
C
C Jan. 2000 : we are now not assuming row front = col. front
      MFRONT = ISAVE(23)
      NFRONT = ISAVE(24)
      IF (ISAVE(31).EQ.-1) THEN
C No direct access data sets used.
         DIMIBF = L3
         IF (LIW.LT.DIMIBF) GO TO 70
         IF (IND.EQ.1) THEN
            LLW = NRHS*NFRONT + L1 + L2
            IF (LW.LT.LLW) GO TO 60
            CALL MA42DD(1,NRHS,LX,X,.FALSE.,1,1,B,W(L2+1),L1,
     +                  IW,DIMIBF,W(L2+L1+1),NFRONT,LP,NRHS,
     +                  IFILE,IREC,ISIZE,MKEY,JFLAG)
         ELSE IF (IND.EQ.2) THEN
            LLW = NRHS*MFRONT + L1 + L2
            IF (LW.LT.LLW) GO TO 60
            CALL MA42ED(2,NRHS,LX,B,W,L2,IW,DIMIBF,
     +                  W(L2+L1+1),MFRONT,LP,NRHSB,IFILE,
     +                  ISIZE,MKEY,JFLAG)
         ELSE IF (IND.EQ.3) THEN
            LLW = NRHS*NFRONT + L1 + L2
            IF (LW.LT.LLW) GO TO 60
            CALL MA42DD(1,NRHS,LX,X,.TRUE.,LX,NRHS,B,W(L2+1),
     +                  L1,IW,DIMIBF,W(L2+L1+1),NFRONT,LP,
     +                  NRHSB,IFILE,IREC,ISIZE,MKEY,JFLAG)
         END IF
      ELSE
C Direct access data sets used.
         DIMIBF = L3 + 5 + ISAVE(23) + ISAVE(24)
         IF (LIW.LT.DIMIBF) GO TO 70
         IS = ISAVE(23)*ISAVE(24)
C April 1999: found bug. The workspace could be too small
C if MA42B/BD was called with a large number of right-hand sides
C         LBUFR = L1 + IS
         IF (IND.EQ.1) THEN
            LBUFR = L1 + IS + NRHSB*NFRONT
            LLW = NRHS*NFRONT + LBUFR
            IF (LW.LT.LLW) GO TO 60
            CALL MA42DD(1,NRHS,LX,X,.FALSE.,1,1,B,W,LBUFR,IW,
     +                  DIMIBF,W(LBUFR+1),NFRONT,LP,NRHS,
     +                  IFILE,IREC,ISIZE,MKEY,JFLAG)
         ELSE IF (IND.EQ.2) THEN
            LBUFR = L2 + IS
            LLW = NRHS*MFRONT + LBUFR
            IF (LW.LT.LLW) GO TO 60
            CALL MA42ED(2,NRHS,LX,B,W,LBUFR,IW,DIMIBF,
     +                  W(LBUFR+1),MFRONT,LP,NRHSB,IFILE,
     +                  ISIZE,MKEY,JFLAG)
         ELSE IF (IND.EQ.3) THEN
            LBUFR = L1 + IS + NRHSB*NFRONT
            LLW = NRHS*NFRONT + LBUFR
            IF (LW.LT.LLW) GO TO 60
            CALL MA42DD(1,NRHS,LX,X,.TRUE.,LX,NRHS,B,W,LBUFR,
     +                  IW,DIMIBF,W(LBUFR+1),NFRONT,LP,NRHSB,
     +                  IFILE,IREC,ISIZE,MKEY,JFLAG)
         END IF
      END IF
      IF (JFLAG.LT.0) GO TO 90
      GO TO 120
C
C **** Error returns ****
   40 IFLAG = -1
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) IFLAG
         WRITE (LP,FMT=9060)
      END IF
      GO TO 120
   50 IFLAG = -2
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) IFLAG
         WRITE (LP,FMT=9040) LX,ISAVE(30)
      END IF
      GO TO 120
   60 IFLAG = -3
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) IFLAG
         WRITE (LP,FMT=9030) LW,LLW
      END IF
      GO TO 120
   70 IFLAG = -4
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) IFLAG
         WRITE (LP,FMT=9020) LIW,DIMIBF
      END IF
      GO TO 120
   80 IFLAG = -5
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) IFLAG
         IF (NRHS.LE.0) WRITE (LP,FMT=9050) NRHS
         IF (IND.EQ.1 .AND. NRHS.NE.ISAVE(22)) WRITE (LP,FMT=9010) NRHS
      END IF
      GO TO 120
   90 IFLAG = -6
      IF (LP.GT.0) WRITE (LP,FMT=9000) IFLAG
      GO TO 120
  100 IFLAG = -7
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) IFLAG
         IF (ISAVE(31).EQ.1) WRITE (LP,FMT=9070) IFILE(2)
         IF (ISAVE(31).EQ.-1) WRITE (LP,FMT=9080) ISIZE(2)
      END IF
      GO TO 120
  110 IFLAG = -8
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) IFLAG
         WRITE (LP,FMT=9090) IND
      END IF
C
  120 RETURN
 9000 FORMAT (' ***** Error return from MA52C/CD ***** IFLAG = ',I3)
 9010 FORMAT (7X,'Number of right hand sides (NRHS) is ',I8,/7X,'This ',
     +       'has been changed since call to MA52B/BD (or MA52F/FD).')
 9020 FORMAT (7X,'Length of array IW too small.',/7X,'Increase LIW ',
     +       'from ',I8,' to ',I8)
 9030 FORMAT (7X,'Length of array W too small.',/7X,'Increase LW ',
     +       'from ',I8,' to ',I8)
 9040 FORMAT (7X,'First dimension of array X too small.',/7X,'Increase',
     +       ' LX from ',I8,' to ',I8)
 9050 FORMAT (7X,'NRHS is out of range. NRHS = ',I8)
 9060 FORMAT (7X,'MA52B/BD (or MA52F/FD) was not called previously.')
 9070 FORMAT (7X,'No record of L-factor being stored.',/7X,'MA42P/PD ',
     +       'was called with ISTRM(2) set to ',I3)
 9080 FORMAT (7X,'No record of L-factor being stored.',/7X,'MA42B/BD ',
     +       'was called with LENBUF(2) set to ',I8)
 9090 FORMAT (7X,'Value of IND out of range. IND = ',I6)
      END
C*********************************************************************
      SUBROUTINE MA52DD(NFRONT,NRHS,LRHS,BUFRL,LLB,BUFRU,LUB,IBUFR,
     +                  LIBUFR,FA,FARHS,LHED,IFILE,IREC,ISIZE,MKEY,
     +                  NUMBLK,IFVAR,FVAR,FRHS,LFVAR,NFVAR,ISAVE,LP,
     +                  IFLAG)
C
C     Argument list. * indicates the argument is changed by the routine.
C
C  NFRONT - Integer variable.  Set by user of MA42 to max number
C           of rows/cols in front.
C  NRHS   - Integer variable.  Holds number of right-hand sides.
C  LRHS   - Integer variable. Second dimension of array FARHS.
C  BUFRL  - Real (DP) array of length LLB. In-core buffer holding
C           entries of L.
C  LLB    - Integer variable. Length of array BUFRL.
C  BUFRU  - Real (DP) array of length LUB. In-core buffer holding
C           entries of U.
C  LUB    - Integer variable. Length of array BUFRU.
C  IBUFR  _ Integer array of length LIBUFR. In-core buffer holding
C           integer information on factors L and U.
C  LIBUFR - Integer variable. Length of array IBUFR.
C  FA     - Real (DP) array with dimensions NFRONT,NFRONT.
C           Holds the frontal matrix.
C  FARHS   - Real (DP) array with dimensions NFRONT, LRHS.
C           Holds the right hand sides corresponding to the current
C           frontal matrix.
C  LHED   - Integer array of dimension NFRONT. Used to hold the indices
C           of variables corresponding to columns in the front.
C  IFILE  - Integer array of length 3. Stream numbers for direct
C           access data sets.
C *IREC   - Integer array of length 3. Pointers to first free space
C           in buffers.
C  ISIZE  - Integer array of length 3. Length of buffers.
C *MKEY   - Integer array of length 3. Number of records written
C           to direct access data sets.
C  NUMBLK - Integer array of length 3. Number of records in
C           direct access data sets.
C *IFVAR  - Integer array of dimensions (LFVAR,LFVAR).
C           Not set on entry. On exit holds the indices
C           of the variables left in the front after last
C           call to MA42B/B.
C *FVAR   - Real (dp) array of dimensions (LFVAR,LFVAR).
C           Not set on entry. On exit holds values for the
C           variables left in the front after last
C           call to MA42B/B.
C *FRHS   - Real (dp) array of dimensions (LFVAR,LRHS).
C           Not set on entry. On exit holds right-hand sides
C           for the variables left in the front after last
C           call to MA42B/B.
C  LFVAR  - Integer variable. Defines dimension of IFVAR, FVAR,
C  NFVAR  - Integer variable. Number of variables in front.
C   LP     - Stream for error messages.
C  *IFLAG  - Error flag. Negative values indicate a fatal error.
C            Possible nonzero values are -6,-7.
C
C  Local variables
C
C  I      - Do loop variable.
C  JFLAG  - Error flag when calling MA42L/LD.
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER IFLAG,LFVAR,LIBUFR,LLB,LP,LRHS,LUB,NFRONT,NFVAR,NRHS
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION BUFRL(LLB),BUFRU(LUB),FA(NFRONT,NFRONT),
     +                 FARHS(NFRONT,LRHS),FRHS(LFVAR,LRHS),
     +                 FVAR(LFVAR,LFVAR)
      INTEGER IBUFR(LIBUFR),IFILE(3),IFVAR(LFVAR),IREC(3),ISAVE(31),
     +        ISIZE(3),LHED(NFRONT),MKEY(3),NUMBLK(3)
C     ..
C     .. Local Scalars ..
      INTEGER I,J,JFLAG
C     ..
C     .. External Subroutines ..
      EXTERNAL MA42LD
C     ..
C Jump if direct access not used.
      IF (ISAVE(31).EQ.-1) GO TO 10

      JFLAG = 0
C Write out contents of buffers to direct access data sets.

C Unless IREC(1)=1, BUFRU still contains reals to be output
C to direct access data set.
      IF (IREC(1).NE.1) THEN
C Check there is room to write.
         IF (MKEY(1)+1.GT.NUMBLK(1)) GO TO 70
         MKEY(1) = MKEY(1) + 1
C Fill-in the rest of BUFRU with garbage (so that it is
C not undefined ... can only be undefined if MKEY(1) = 1
C so that the whole buffer has not already been used)
         IF (MKEY(1).EQ.1) THEN
            DO 5 I = IREC(1),LUB
               BUFRU(I) = ZERO
   5        CONTINUE
         END IF
         CALL MA42LD(3,IFILE(1),MKEY(1),BUFRU,LUB,IBUFR,LIBUFR,
     +               LP,JFLAG)
         IF (JFLAG.LT.0) GO TO 80
      END IF

C Unless IREC(2)=1, BUFRL still contains reals to be output
C to direct access data set. If ISIZE(2) = 0 the L factor
C is not being stored.
      IF (ISIZE(2).NE.0 .AND. IREC(2).NE.1) THEN
C Check there is room to write.
         IF (MKEY(2)+1.GT.NUMBLK(2)) GO TO 70
         MKEY(2) = MKEY(2) + 1
C Fill-in the rest of BUFRL with garbage (so that it is
C not undefined ... can only be undefined if MKEY(2) = 1
C so that the whole buffer has not already been used)
         IF (MKEY(2).EQ.1) THEN
            DO 6 I = IREC(2),LLB
               BUFRL(I) = ZERO
   6        CONTINUE
         END IF
         CALL MA42LD(3,IFILE(2),MKEY(2),BUFRL,LLB,IBUFR,LIBUFR,
     +               LP,JFLAG)
         IF (JFLAG.LT.0) GO TO 80
      END IF

C Check there is room to write.
      IF (MKEY(3)+1.GT.NUMBLK(3)) GO TO 70
C Set a flag in IBUFR to indicate end of file (IBUFR(IREC(3))=0).
   10 IBUFR(IREC(3)) = 0
C Fill-in the rest of IBUFR with garbage (so that it is
C not undefined)
      DO 15 I = IREC(3)+1,LIBUFR
         IBUFR(I) = -1
   15 CONTINUE
      IREC(3) = IREC(3) + 1
      MKEY(3) = MKEY(3) + 1
      CALL MA42LD(-3,IFILE(3),MKEY(3),BUFRU,LUB,IBUFR,LIBUFR,
     +            LP,JFLAG)
      IF (JFLAG.LT.0) GO TO 80
C
C Set IFVAR, FVAR and FRHS.
      DO 20 I = 1,NFVAR
         IFVAR(I) = LHED(I)
   20 CONTINUE
      DO 40 I = 1,NFVAR
         DO 30 J = 1,NFVAR
            FVAR(J,I) = FA(J,I)
   30    CONTINUE
C         CALL DCOPY(NFVAR,FA(1,I),1,FVAR(1,I),1)
   40 CONTINUE
      DO 60 I = 1,NRHS
         DO 50 J = 1,NFVAR
            FRHS(J,I) = FARHS(J,I)
   50    CONTINUE
C         CALL DCOPY(NFVAR,FARHS(1,I),1,FRHS(1,I),1)
   60 CONTINUE
      GO TO 90
C
C **** Error returns ****
   70 IFLAG = -6
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) IFLAG
         IF (MKEY(1)+1.GT.NUMBLK(1)) WRITE (LP,FMT=9010) ISIZE(1),
     +       ISIZE(1)* (MKEY(1)+1)
         IF (MKEY(2)+1.GT.NUMBLK(2)) WRITE (LP,FMT=9020) ISIZE(2),
     +       ISIZE(2)* (MKEY(2)+1)
         IF (MKEY(3)+1.GT.NUMBLK(3)) WRITE (LP,FMT=9030) ISIZE(3),
     +       ISIZE(3)* (MKEY(3)+1)
      END IF
      GO TO 90
   80 IFLAG = -7
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) IFLAG
      END IF
   90 RETURN
C
 9000 FORMAT (' ***** Error return from MA52B/BD *****  IFLAG = ',I3)
 9010 FORMAT (7X,'Insufficent storage for U-factor.',/7X,'With present',
     +       ' LENBUF(1) size of ',I8,/7X,'LENFLE(1) must be at least ',
     +       I8)
 9020 FORMAT (7X,'Insufficent storage for L-factor.',/7X,'With present',
     +       ' LENBUF(2) size of ',I8,/7X,'LENFLE(2) must be at least ',
     +       I8)
 9030 FORMAT (7X,'Insufficent storage for integers.',/7X,'With present',
     +       ' LENBUF(3) size of ',I8,/7X,'LENFLE(3) must be at least ',
     +       I8)
      END
C*****************************************************************
      SUBROUTINE MA52ED(LAST,NDF,NFRONT,LHED,KPVLNK)
C
C Reset LAST.
C
C     Argument list. * indicates the argument is changed by the routine.
C
C *LAST   - Integer array of length NDF.  On exit LAST is reset
C           to its state on first entry to MA42B/BD.
C  NDF    - Integer variable. Length of array LAST.
C  NFRONT - Integer variable. Length of arrays LHED, KPVLNK.
C  LHED   - Integer array of dimension NFRONT. Holds the indices
C           of variables corresponding to columns in the front.
C  KPVLNK - Integer array of dimension NFRONT.
C           For each non-fully summed variable in the
C           front, KPVLNK holds the information
C           originally in LAST (i.e. the elt/eqn at which that
C           variable becomes fully summed).
C
C Local variables
C
C  L      - Integer variable. Do loop variable.
C  MFR    - Integer variable.  Holds index of Lth variable in
C           KPVLNK.
C
C     .. Scalar Arguments ..
      INTEGER NDF,NFRONT
C     ..
C     .. Array Arguments ..
      INTEGER KPVLNK(NFRONT),LAST(NDF),LHED(NFRONT)
C     ..
C     .. Local Scalars ..
      INTEGER L,MFR
C     ..
      DO 10 L = 1,NFRONT
         IF (KPVLNK(L).GT.0) GO TO 10
         MFR = LHED(L)
         LAST(MFR) = -KPVLNK(L)
   10 CONTINUE
      RETURN
C
      END

C**********************************************************************
      SUBROUTINE MA52FD(LRHS,NRHS,NDF,LAST,LX,X,LW,W,LIW,IW,NFVAR,IFVAR,
     +                  JFVAR,FVAR,FRHS,LFVAR,ISAVE,LP,IFLAG)
C
C *** Copyright Rutherford Appleton Laboratory  November 1993 ***
C *** Any problems contact Jennifer A. Scott at Atlas Centre,
C     Rutherford Appleton Laboratory ***
C
C  This routine dumps the contents of the  buffers into the direct
C  access files used by MA42 and writes what remains in the
C  frontal matrix (and corresponding frontal right-hand sides)
C  into the form of an element matrix and element right-hand side.
C  This routine can only be used after MA42B/BD has been called
C  for at least one element but has not been called for all
C  all the elements in the problem.
C*** This subroutine is the same as MA52B/BD except
C    row and column indices of variables in front are preserved
C    and we do not assume diagonal pivoting.
C    We are also NOT assuming frontal matrix is square
C    (so that NFVAR and LFVAR are now arrays length 2).
C
C     Argument list. * indicates the argument is changed by the routine.
C
C  LRHS   - Integer variable. Must be unchanged since MA42B/BD.
C  NRHS   - Integer variable. Must be unchanged since MA42B/BD.
C  NDF    - Integer variable. Length of array LAST.
C *LAST   - Integer array of length NDF.  On exit LAST is reset
C           to its state on first entry to MA42B/BD.
C  LX     - Integer variable. First dimension of  array X.
C *X      - Real (DP) array of dimensions LX by LRHS. On exit,
C           X(I,J)=0, I=1,2,...,NDF, J=1,2,...,NRHS.
C  LW     - Integer variable. Length of array W.
C  W      - Real (DP) array of length LW.
C  LIW    - Integer variable. Length of array IW.
C  IW     - Integer array of length LIW.
C  LRHS, LAST, NDF, LW, W, LIW, IW must all be unchanged from the
C           last call to MA42B/BD.
C *NFVAR  - Integer array length 2. Not set on entry.
C           On exit, holds the number of rows and columns
C           left in the front after last call to MA42B/BD
C           (=ISAVE(39),ISAVE(40)).
C *IFVAR  - Integer array of dimension LFVAR(1).
C           Not set on entry. On exit holds the
C           row indices left in the front after last
C           call to MA42B/B.
C *JFVAR  - Integer array of dimensions LFVAR(2).
C           Not set on entry. On exit holds the indices
C           of the variables left in the front after last
C           call to MA42B/B.(corresponds to column indices)
C *FVAR   - Real (dp) array of dimensions (LFVAR(1),LFVAR(2)).
C           Not set on entry. On exit holds values for the
C           variables left in the front after last
C           call to MA42B/B.
C *FRHS   - Real (dp) array of dimensions (LFVAR(1),LRHS).
C           Not set on entry. On exit holds right-hand sides
C           for the variables left in the front after last
C           call to MA42B/B.
C  LFVAR  - Integer array length 2. Defines dimension of IFVAR, JFVAR,
C           FVAR, and FRHS. LFVAR(1) must be at least ISAVE(39)
C           and LFVAR(2) must be at least ISAVE(40).
C *ISAVE  - Integer array of length 45. Holds parameters which
C           must be preserved between calls to MA42B/BD and between
C           calls to other MA42 routines.
C           ISAVE(1),...,ISAVE(15) were set by MA42I/ID and MA42P/PD.
C           If MA42P/PD was called, ISAVE(6+I) was set to
C           LENBUF(I), I=1,2,3 on the first entry to MA42B/BD.
C           ISAVE(16) used by MA42B/BD if MA42J/JD called.
C           ISAVE(17) is set to 0 by MA42I/ID.
C           On exit from intermediate calls to MA42B/BD, ISAVE(17)
C           holds the number of elements input so far and on exit
C           from the final call to MA42B/BD is set to 0.
C           On return from MA52B/BD,
C           ISAVE(17) = 0 (this allows MA42B/BD
C           to be rerun without rerunning MA42A/AD).
C           ISAVE(18) set by MA42A/AD to number of elements in problem
C           and is unchanged.
C           A non-zero value of ISAVE(19) or ISAVE(20) indicates
C           the occurrence of the non-fatal error 5 or 6
C           respectively, in MA42B/BD.
C           ISAVE(21), ISAVE(22), ISAVE(23) and ISAVE(24)
C           hold the values of variables NMAXE, NRHS,
C           NFRONT(1) and NFRONT(2), resp. as input in the first call
C           to MA42B/BD.
C           ISAVE(25)-ISAVE(29) are used to divide real workspace.
C           ISAVE(30)-set by MA42A/AD to largest integer used to index
C           a variable.
C           ISAVE(31)-flag to indicate if MA42P/PD has been called.
C           ISAVE(32)-ISAVE(37) are used to divide integer workspace.
C           ISAVE(38)- set by MA42B/BD to hold a copy of the
C           error flag INFO(1) if error encountered by MA42B/BD
C           and holds copy of ICNTL(7) from MA42I/ID o.w.
C           On exit, ISAVE(38) set to 3003 to indicate MA52B/BD
C           has been called.
C           ISAVE(39) is set to 0 by MA42I/ID. On exit from
C           intermediate calls to MA42B/BD, ISAVE(39) holds current row
C           front size KFRNT.
C           On return from MA52B/BD, ISAVE(39)=0.
C           ISAVE(40) is set to 0 by MA42I/ID. On exit from
C           intermediate calls to MA42B/BD, ISAVE(40) holds current
C           column front size LFRNT.
C           On return from MA52B/BD, ISAVE(40)=0.
C           ISAVE(41)-ISAVE(45) used to preserve values between calls
C           to MA42F/FD.
C     (ISAVE is altered in the way it would have been if the series
C      of calls to MA42B/BD had been completed, except that the number
C      of records used in the direct access data sets may be different
C      so MKEY and IREC may differ... ISAVE(3+I), ISAVE(9+I)).
C     This array must be unchanged since calls to MA42B/BD.
C   LP    - Integer variable. Set by user to stream number for error
C           messages. Printing suppressed if LP.LE.0
C *IFLAG  - Integer variable. Not set on entry. On successful
C           exit, set to 0. Negative value indicates error.
C
C   Local variables
C
C   FA     - Integer variable. Used in dividing real workspace.
C   FARHS  - Integer variable. Used in dividing real workspace.
C   I      - Integer variable. Do loop variable.
C   J      - Integer variable. Do loop variable.
C   JBUFL  - Integer variable. Holds max(1,ISIZE(2)). Used in dividing
C            real workspace.
C   JBUFR  - Integer variable. Holds 1+SIZE(2). Used in dividing
C            real workspace.
C   JBUFU  - Integer variable. Holds ISIZE(1) (size of U-buffer).
C            Used in dividing real workspace.
C   KPVLNK - Integer variable. Used in dividing integer workspace.
C   LHED   - Integer variable. Used in dividing integer workspace.
C   IFILE  - Integer array of length 3. IFILE(I) is equivalent to
C            ISAVE(I) (I = 1,2,3).
C   IREC   - Integer array of length 3. IREC(I) is equivalent to
C            ISAVE(3+I) (I = 1,2,3).
C   ISIZE  - Integer array of length 3. ISIZE(I) is equivalent to
C            ISAVE(6+I) (I = 1,2,3).
C   MKEY   - Integer array of length 3. MKEY(I) is equivalent to
C            ISAVE(9+I) (I = 1,2,3).
C   NUMBLK - Integer array of length 3. NUMBLK(I) is equivalent to
C            ISAVE(12+I) (I = 1,2,3).
C
C     .. Scalar Arguments ..
      INTEGER IFLAG,LIW,LP,LRHS,LW,LX,NDF,NRHS
      INTEGER LFVAR(2),NFVAR(2)
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION FRHS(LFVAR(1),LRHS),FVAR(LFVAR(1),LFVAR(2)),
     +                 W(LW),X(LX,LRHS)
      INTEGER IFVAR(LFVAR(1)),JFVAR(LFVAR(2)),ISAVE(45),IW(LIW),
     +        LAST(NDF)
C     ..
C     .. Local Scalars ..
      INTEGER FA,FARHS,I,J,JBUFL,JBUFR,JBUFU,KPVLNK,KHED,LHED,MFRONT,
     +        NFRONT
C     ..
C     .. Local Arrays ..
      INTEGER IFILE(3),IREC(3),ISIZE(3),MKEY(3),NUMBLK(3)
C     ..
C     .. External Subroutines ..
      EXTERNAL MA52GD,MA52ED
C     ..
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C     ..
      IFLAG = 0
C Immediate return if MA42B/BD has issued an error.
      IF (ISAVE(38).LT.0 .OR. ISAVE(38).EQ.4 .OR. ISAVE(38).EQ.5 .OR.
     +    ISAVE(38).EQ.6) GO TO 50
C Check diagonal pivoting was used by MA42B/BD
C Check that MA42B/BD was not called for all the elements
C in the problem. If it was, then immediate return.
      IF (ISAVE(17).EQ.0) GO TO 60
      DO 10 I = 1,3
         IFILE(I) = ISAVE(I)
         IREC(I) = ISAVE(3+I)
         ISIZE(I) = ISAVE(6+I)
         MKEY(I) = ISAVE(9+I)
         NUMBLK(I) = ISAVE(12+I)
   10 CONTINUE
C Test validity of input parameters.
      NFVAR(1) = ISAVE(39)
      NFVAR(2) = ISAVE(40)
      IF (NRHS.NE.ISAVE(22)) GO TO 80
      IF (LFVAR(1).LT.NFVAR(1)) GO TO 70
      IF (LFVAR(2).LT.NFVAR(2)) GO TO 70
      IF (NRHS.GT.0 .AND. LX.LT.NDF) GO TO 90
C
C Perform dump of buffer contents to direct access files
      MFRONT = ISAVE(23)
      NFRONT = ISAVE(24)
      JBUFU = ISAVE(25)
      JBUFL = ISAVE(26)
      JBUFR = ISAVE(27)
      FA = ISAVE(28)
      FARHS = ISAVE(29)
      LHED = ISAVE(32)
      KHED = ISAVE(33)
C
      CALL MA52GD(MFRONT,NFRONT,NRHS,LRHS,W,JBUFL,W(JBUFR),JBUFU,IW,
     +            ISIZE(3),W(FA),W(FARHS),IW(KHED),IW(LHED),
     +            IFILE,IREC,ISIZE,MKEY,NUMBLK,
     +            IFVAR,JFVAR,FVAR,FRHS,LFVAR,NFVAR,ISAVE,LP,IFLAG)
      IF (IFLAG.LT.0) GO TO 110
C Preserve IREC and MKEY.
      DO 20 I = 1,3
         ISAVE(3+I) = IREC(I)
         ISAVE(9+I) = MKEY(I)
   20 CONTINUE
C
C Reset ISAVE(17), ISAVE(39), ISAVE(40) to 0
C and reset ISAVE(16) to 2 (so MA42A/AD and MA42B/BD
C may be recalled without recalling MA42I/ID).
C
      ISAVE(17) = 0
      ISAVE(39) = 0
      ISAVE(40) = 0
      ISAVE(16) = 2
C Reset LAST
      KPVLNK = ISAVE(35)
      CALL MA52ED(LAST,NDF,NFVAR(2),IW(LHED),IW(KPVLNK))
C Set X
      DO 40 J = 1,NRHS
         DO 30 I = 1,NDF
            X(I,J) = ZERO
   30    CONTINUE
   40 CONTINUE
C Set ISAVE(38) to indicate call to MA52B/BD complete.
      ISAVE(38) = 3003
      GO TO 110
C
C *** Error returns ***
   50 IFLAG = -1
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) IFLAG
         WRITE (LP,FMT=9040)
      END IF
      GO TO 110
   60 IFLAG = -2
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) IFLAG
         WRITE (LP,FMT=9050)
      END IF
      GO TO 110
   70 IFLAG = -3
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) IFLAG
         WRITE (LP,FMT=9010) NFVAR
      END IF
      GO TO 110
   80 IFLAG = -4
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) IFLAG
         WRITE (LP,FMT=9020)
      END IF
      GO TO 110
   90 IFLAG = -5
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) IFLAG
         WRITE (LP,FMT=9030) NDF
      END IF
      GO TO 110
C
  110 RETURN
 9000 FORMAT (' ***** Error return from MA52F/FD *****  IFLAG = ',I3)
 9010 FORMAT (7X,'Length of arrays too small.',/7X,'LFVAR ',
     +       'must have dimensions at least ',2I8)
 9020 FORMAT (7X,'NRHS has been changed since the last call to ',
     +       'MA42B/BD.')
 9030 FORMAT (7X,'LX too small. Must be at least ',I8)
 9040 FORMAT (7X,'An error was issued by MA42B/BD.')
 9050 FORMAT (7X,'MA42B/BD was called for all the elements in the ',
     +       'problem.')
      END
C*********************************************************************
      SUBROUTINE MA52GD(MFRONT,NFRONT,NRHS,LRHS,BUFRL,LLB,BUFRU,
     +                  LUB,IBUFR,LIBUFR,FA,FARHS,KHED,LHED,
     +                  IFILE,IREC,ISIZE,MKEY,
     +                  NUMBLK,IFVAR,JFVAR,FVAR,FRHS,LFVAR,NFVAR,
     +                  ISAVE,LP,IFLAG)
C
C     Argument list. * indicates the argument is changed by the routine.
C
C  MFRONT - Integer variable.  Set by user of MA42 to max number
C           of rows in front.
C  NFRONT - Integer variable.  Set by user of MA42 to max number
C           of cols in front.
C  NRHS   - Integer variable.  Holds number of right-hand sides.
C  LRHS   - Integer variable. Second dimension of array FARHS.
C  BUFRL  - Real (DP) array of length LLB. In-core buffer holding
C           entries of L.
C  LLB    - Integer variable. Length of array BUFRL.
C  BUFRU  - Real (DP) array of length LUB. In-core buffer holding
C           entries of U.
C  LUB    - Integer variable. Length of array BUFRU.
C  IBUFR  _ Integer array of length LIBUFR. In-core buffer holding
C           integer information on factors L and U.
C  LIBUFR - Integer variable. Length of array IBUFR.
C  FA     - Real (DP) array with dimensions MFRONT,NFRONT.
C           Holds the frontal matrix.
C  FARHS   - Real (DP) array with dimensions MFRONT, LRHS.
C           Holds the right hand sides corresponding to the current
C           frontal matrix.
C  KHED   - Integer array of dimension MFRONT. Used to hold the indices
C           of variables corresponding to rows in the front.
C  LHED   - Integer array of dimension NFRONT. Used to hold the indices
C           of variables corresponding to columns in the front.
C  IFILE  - Integer array of length 3. Stream numbers for direct
C           access data sets.
C *IREC   - Integer array of length 3. Pointers to first free space
C           in buffers.
C  ISIZE  - Integer array of length 3. Length of buffers.
C *MKEY   - Integer array of length 3. Number of records written
C           to direct access data sets.
C  NUMBLK - Integer array of length 3. Number of records in
C           direct access data sets.
C *IFVAR  - Integer array of dimension LFVAR(1).
C           Not set on entry. On exit holds the indices
C           of the rows left in the front after last
C           call to MA42B/B.
C *JFVAR  - Integer array of dimension LFVAR(2).
C           Not set on entry. On exit holds the indices
C           of the variables left in the front after last
C           call to MA42B/B.
C *FVAR   - Real (dp) array of dimensions (LFVAR(1),LFVAR(2)).
C           Not set on entry. On exit holds values for the
C           variables left in the front after last
C           call to MA42B/B.
C *FRHS   - Real (dp) array of dimensions (LFVAR(1),LRHS).
C           Not set on entry. On exit holds right-hand sides
C           for the variables left in the front after last
C           call to MA42B/B.
C  LFVAR  - Integer array length 2. Defines dimension of IFVAR,
C           JFVAR, FVAR,
C  NFVAR  - Integer array length 2. Number of rows/cols in front.
C   LP     - Stream for error messages.
C  *IFLAG  - Error flag. Negative values indicate a fatal error.
C            Possible nonzero values are -6,-7.
C
C  Local variables
C
C  I      - Do loop variable.
C  JFLAG  - Error flag when calling MA42L/LD.
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER IFLAG,LIBUFR,LLB,LP,LRHS,LUB,MFRONT,NFRONT,NRHS
      INTEGER LFVAR(2),NFVAR(2)
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION BUFRL(LLB),BUFRU(LUB),FA(MFRONT,NFRONT),
     +                 FARHS(MFRONT,LRHS),FRHS(LFVAR(1),LRHS),
     +                 FVAR(LFVAR(1),LFVAR(2))
      INTEGER IBUFR(LIBUFR),IFILE(3),IFVAR(LFVAR(1)),JFVAR(LFVAR(2)),
     +        IREC(3),ISAVE(31),
     +        ISIZE(3),LHED(NFRONT),MKEY(3),NUMBLK(3),
     +        KHED(MFRONT)
C     ..
C     .. Local Scalars ..
      INTEGER I,J,JFLAG
C     ..
C     .. External Subroutines ..
      EXTERNAL MA42LD
C     ..
C Jump if direct access not used.
      IF (ISAVE(31).EQ.-1) GO TO 10

      JFLAG = 0
C Write out contents of buffers to direct access data sets.

C Unless IREC(1)=1, BUFRU still contains reals to be output
C to direct access data set.
      IF (IREC(1).NE.1) THEN
C Check there is room to write.
         IF (MKEY(1)+1.GT.NUMBLK(1)) GO TO 70
         MKEY(1) = MKEY(1) + 1
C Fill-in the rest of BUFRU with garbage (so that it is
C not undefined ... can only be undefined if MKEY(1) = 1
C so that the whole buffer has not already been used)
         IF (MKEY(1).EQ.1) THEN
            DO 5 I = IREC(1),LUB
               BUFRU(I) = ZERO
   5        CONTINUE
         END IF
         CALL MA42LD(3,IFILE(1),MKEY(1),BUFRU,LUB,IBUFR,LIBUFR,
     +               LP,JFLAG)
         IF (JFLAG.LT.0) GO TO 80
      END IF

C Unless IREC(2)=1, BUFRL still contains reals to be output
C to direct access data set. If ISIZE(2) = 0 the L factor
C is not being stored.
      IF (ISIZE(2).NE.0 .AND. IREC(2).NE.1) THEN
C Check there is room to write.
         IF (MKEY(2)+1.GT.NUMBLK(2)) GO TO 70
         MKEY(2) = MKEY(2) + 1
C Fill-in the rest of BUFRL with garbage (so that it is
C not undefined ... can only be undefined if MKEY(2) = 1
C so that the whole buffer has not already been used)
         IF (MKEY(2).EQ.1) THEN
            DO 6 I = IREC(2),LLB
               BUFRL(I) = ZERO
   6        CONTINUE
         END IF
         CALL MA42LD(3,IFILE(2),MKEY(2),BUFRL,LLB,IBUFR,LIBUFR,
     +               LP,JFLAG)
         IF (JFLAG.LT.0) GO TO 80
      END IF

C Check there is room to write.
      IF (MKEY(3)+1.GT.NUMBLK(3)) GO TO 70
C Set a flag in IBUFR to indicate end of file (IBUFR(IREC(3))=0).
   10 IBUFR(IREC(3)) = 0
C Fill-in the rest of IBUFR with garbage (so that it is
C not undefined)
      DO 15 I = IREC(3)+1,LIBUFR
         IBUFR(I) = -1
   15 CONTINUE
      IREC(3) = IREC(3) + 1
      MKEY(3) = MKEY(3) + 1
      CALL MA42LD(-3,IFILE(3),MKEY(3),BUFRU,LUB,IBUFR,LIBUFR,
     +            LP,JFLAG)
      IF (JFLAG.LT.0) GO TO 80
C
C Set JFVAR, FVAR and FRHS.
      DO 20 I = 1,NFVAR(1)
         IFVAR(I) = KHED(I)
   20 CONTINUE
      DO 25 J = 1,NFVAR(2)
         JFVAR(J) = LHED(J)
   25 CONTINUE
      DO 40 J = 1,NFVAR(2)
         DO 30 I = 1,NFVAR(1)
            FVAR(I,J) = FA(I,J)
   30    CONTINUE
   40 CONTINUE
      DO 60 J = 1,NRHS
         DO 50 I = 1,NFVAR(1)
            FRHS(I,J) = FARHS(I,J)
   50    CONTINUE
   60 CONTINUE
      GO TO 90
C
C **** Error returns ****
   70 IFLAG = -6
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) IFLAG
         IF (MKEY(1)+1.GT.NUMBLK(1)) WRITE (LP,FMT=9010) ISIZE(1),
     +       ISIZE(1)* (MKEY(1)+1)
         IF (MKEY(2)+1.GT.NUMBLK(2)) WRITE (LP,FMT=9020) ISIZE(2),
     +       ISIZE(2)* (MKEY(2)+1)
         IF (MKEY(3)+1.GT.NUMBLK(3)) WRITE (LP,FMT=9030) ISIZE(3),
     +       ISIZE(3)* (MKEY(3)+1)
      END IF
      GO TO 90
   80 IFLAG = -7
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) IFLAG
      END IF
   90 RETURN
C
 9000 FORMAT (' ***** Error return from MA52F/FD *****  IFLAG = ',I3)
 9010 FORMAT (7X,'Insufficent storage for U-factor.',/7X,'With present',
     +       ' LENBUF(1) size of ',I8,/7X,'LENFLE(1) must be at least ',
     +       I8)
 9020 FORMAT (7X,'Insufficent storage for L-factor.',/7X,'With present',
     +       ' LENBUF(2) size of ',I8,/7X,'LENFLE(2) must be at least ',
     +       I8)
 9030 FORMAT (7X,'Insufficent storage for integers.',/7X,'With present',
     +       ' LENBUF(3) size of ',I8,/7X,'LENFLE(3) must be at least ',
     +       I8)
      END
C**********************************************************************
      SUBROUTINE MA52KD(LRHS,NRHS,NDF,LAST,LX,X,LW,W,LIW,IW,NFVAR,IFVAR,
     +                  JFVAR,IP,FVAR,FRHS,LFVAR,ISAVE,LP,IFLAG)
C
C *** Copyright Rutherford Appleton Laboratory  November 1993 ***
C *** Any problems contact Jennifer A. Scott at Atlas Centre,
C     Rutherford Appleton Laboratory ***
C
C  This routine dumps the contents of the  buffers into the direct
C  access files used by MA42 and writes what remains in the
C  frontal matrix (and corresponding frontal right-hand sides)
C  into the form of a sparse matrix and right-hand side vector.
C  This routine can only be used after MA42B/BD has been called
C  for at least one row  but has not been called for all
C  all the rows in the problem.
C
C     Argument list. * indicates the argument is changed by the routine.
C
C  LRHS   - Integer variable. Must be unchanged since MA42B/BD.
C  NRHS   - Integer variable. Must be unchanged since MA42B/BD.
C  NDF    - Integer variable. Length of array LAST.
C *LAST   - Integer array of length NDF.  On exit LAST is reset
C           to its state on first entry to MA42B/BD.
C  LX     - Integer variable. First dimension of  array X.
C *X      - Real (DP) array of dimensions LX by LRHS. On exit,
C           X(I,J)=0, I=1,2,...,NDF, J=1,2,...,NRHS.
C  LW     - Integer variable. Length of array W.
C  W      - Real (DP) array of length LW.
C  LIW    - Integer variable. Length of array IW.
C  IW     - Integer array of length LIW.
C  LRHS, LAST, NDF, LW, W, LIW, IW must all be unchanged from the
C           last call to MA42B/BD.
C *NFVAR  - Integer array length 2. Not set on entry.
C           On exit, NFVAR(1) holds the number of rows
C           left in the front after last call to MA42B/BD
C           (=ISAVE(39)) and
C           NFVAR(2) is set to hold number of nonzero entries
C           left in front.
C *IFVAR  - Integer array of dimension LFVAR(1).
C           Not set on entry. On exit holds the
C           row indices left in the front after last
C           call to MA42B/B.
C *JFVAR  - Integer array of dimensions LFVAR(2).
C           Not set on entry. On exit holds the column indices
C           of the nonzero entries left in the front after last
C           call to MA42B/B, ordered by rows.
C *IP     - Integer array length LFVAR(1)+1. Row pointer array.
C *FVAR   - Real (dp) array of dimensions LFVAR(2).
C           Not set on entry. On exit holds nonzero entries
C           left in the front after last call to MA42B/B.
C *FRHS   - Real (dp) array of dimensions (LFVAR(1),LRHS).
C           Not set on entry. On exit holds right-hand sides
C           for the variables left in the front after last
C           call to MA42B/B.
C  LFVAR  - Integer array length 2. Defines dimension of IFVAR, JFVAR,
C           FVAR, and FRHS. LFVAR(1) must be at least ISAVE(39).
C *ISAVE  - Integer array of length 45. Holds parameters which
C           must be preserved between calls to MA42B/BD and between
C           calls to other MA42 routines.
C           ISAVE(1),...,ISAVE(15) were set by MA42I/ID and MA42P/PD.
C           If MA42P/PD was called, ISAVE(6+I) was set to
C           LENBUF(I), I=1,2,3 on the first entry to MA42B/BD.
C           ISAVE(16) used by MA42B/BD if MA42J/JD called.
C           ISAVE(17) is set to 0 by MA42I/ID.
C           On exit from intermediate calls to MA42B/BD, ISAVE(17)
C           holds the number of elements input so far and on exit
C           from the final call to MA42B/BD is set to 0.
C           On return from MA52B/BD,
C           ISAVE(17) = 0 (this allows MA42B/BD
C           to be rerun without rerunning MA42A/AD).
C           ISAVE(18) set by MA42A/AD to number of elements in problem
C           and is unchanged.
C           A non-zero value of ISAVE(19) or ISAVE(20) indicates
C           the occurrence of the non-fatal error 5 or 6
C           respectively, in MA42B/BD.
C           ISAVE(21), ISAVE(22), ISAVE(23) and ISAVE(24)
C           hold the values of variables NMAXE, NRHS,
C           NFRONT(1) and NFRONT(2), resp. as input in the first call
C           to MA42B/BD.
C           ISAVE(25)-ISAVE(29) are used to divide real workspace.
C           ISAVE(30)-set by MA42A/AD to largest integer used to index
C           a variable.
C           ISAVE(31)-flag to indicate if MA42P/PD has been called.
C           ISAVE(32)-ISAVE(37) are used to divide integer workspace.
C           ISAVE(38)- set by MA42B/BD to hold a copy of the
C           error flag INFO(1) if error encountered by MA42B/BD
C           and holds copy of ICNTL(7) from MA42I/ID o.w.
C           On exit, ISAVE(38) set to 3003 to indicate MA52B/BD
C           has been called.
C           ISAVE(39) is set to 0 by MA42I/ID. On exit from
C           intermediate calls to MA42B/BD, ISAVE(39) holds current row
C           front size KFRNT.
C           On return from MA52B/BD, ISAVE(39)=0.
C           ISAVE(40) is set to 0 by MA42I/ID. On exit from
C           intermediate calls to MA42B/BD, ISAVE(40) holds current
C           column front size LFRNT.
C           On return from MA52B/BD, ISAVE(40)=0.
C           ISAVE(41)-ISAVE(45) used to preserve values between calls
C           to MA42F/FD.
C     (ISAVE is altered in the way it would have been if the series
C      of calls to MA42B/BD had been completed, except that the number
C      of records used in the direct access data sets may be different
C      so MKEY and IREC may differ... ISAVE(3+I), ISAVE(9+I)).
C     This array must be unchanged since calls to MA42B/BD.
C   LP    - Integer variable. Set by user to stream number for error
C           messages. Printing suppressed if LP.LE.0
C *IFLAG  - Integer variable. Not set on entry. On successful
C           exit, set to 0. Negative value indicates error.
C
C   Local variables
C
C   FA     - Integer variable. Used in dividing real workspace.
C   FARHS  - Integer variable. Used in dividing real workspace.
C   I      - Integer variable. Do loop variable.
C   J      - Integer variable. Do loop variable.
C   JBUFL  - Integer variable. Holds max(1,ISIZE(2)). Used in dividing
C            real workspace.
C   JBUFR  - Integer variable. Holds 1+SIZE(2). Used in dividing
C            real workspace.
C   JBUFU  - Integer variable. Holds ISIZE(1) (size of U-buffer).
C            Used in dividing real workspace.
C   KPVLNK - Integer variable. Used in dividing integer workspace.
C   LHED   - Integer variable. Used in dividing integer workspace.
C   IFILE  - Integer array of length 3. IFILE(I) is equivalent to
C            ISAVE(I) (I = 1,2,3).
C   IREC   - Integer array of length 3. IREC(I) is equivalent to
C            ISAVE(3+I) (I = 1,2,3).
C   ISIZE  - Integer array of length 3. ISIZE(I) is equivalent to
C            ISAVE(6+I) (I = 1,2,3).
C   MKEY   - Integer array of length 3. MKEY(I) is equivalent to
C            ISAVE(9+I) (I = 1,2,3).
C   NUMBLK - Integer array of length 3. NUMBLK(I) is equivalent to
C            ISAVE(12+I) (I = 1,2,3).
C
C     .. Scalar Arguments ..
      INTEGER IFLAG,LIW,LP,LRHS,LW,LX,NDF,NRHS
      INTEGER LFVAR(2),NFVAR(2)
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION FRHS(LFVAR(1),LRHS),FVAR(LFVAR(2)),
     +                 W(LW),X(LX,LRHS)
      INTEGER IFVAR(LFVAR(1)),JFVAR(LFVAR(2)),ISAVE(45),IW(LIW),
     +        LAST(NDF),IP(LFVAR(1)+1)
C     ..
C     .. Local Scalars ..
      INTEGER FA,FARHS,I,J,JBUFL,JBUFR,JBUFU,KPVLNK,KHED,LHED,MFRONT,
     +        NFRONT
C     ..
C     .. Local Arrays ..
      INTEGER IFILE(3),IREC(3),ISIZE(3),MKEY(3),NUMBLK(3)
C     ..
C     .. External Subroutines ..
      EXTERNAL MA52LD,MA52ED
C     ..
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C     ..
      IFLAG = 0
C Immediate return if MA42B/BD has issued an error.
      IF (ISAVE(38).LT.0 .OR. ISAVE(38).EQ.4 .OR. ISAVE(38).EQ.5 .OR.
     +    ISAVE(38).EQ.6) GO TO 50
C Check that MA42B/BD was not called for all the elements
C in the problem. If it was, then immediate return.
      IF (ISAVE(17).EQ.0) GO TO 60
      DO 10 I = 1,3
         IFILE(I) = ISAVE(I)
         IREC(I) = ISAVE(3+I)
         ISIZE(I) = ISAVE(6+I)
         MKEY(I) = ISAVE(9+I)
         NUMBLK(I) = ISAVE(12+I)
   10 CONTINUE
C Test validity of input parameters.
      NFVAR(1) = ISAVE(39)
C Temporarily use NFVAR(2) to hold number of cols left in front
C (later will be used to hold number of entries left in front)
      NFVAR(2) = ISAVE(40)
      IF (NRHS.NE.ISAVE(22)) GO TO 80
      IF (LFVAR(1).LT.NFVAR(1)) GO TO 70
      IF (NRHS.GT.0 .AND. LX.LT.NDF) GO TO 90
C
C Perform dump of buffer contents to direct access files
      MFRONT = ISAVE(23)
      NFRONT = ISAVE(24)
      JBUFU = ISAVE(25)
      JBUFL = ISAVE(26)
      JBUFR = ISAVE(27)
      FA = ISAVE(28)
      FARHS = ISAVE(29)
      LHED = ISAVE(32)
      KHED = ISAVE(33)
C
      CALL MA52LD(MFRONT,NFRONT,NRHS,LRHS,W,JBUFL,W(JBUFR),JBUFU,IW,
     +            ISIZE(3),W(FA),W(FARHS),IW(KHED),IW(LHED),
     +            IFILE,IREC,ISIZE,MKEY,NUMBLK,
     +            IFVAR,JFVAR,IP,FVAR,FRHS,LFVAR,NFVAR,ISAVE,LP,IFLAG)
      IF (IFLAG.LT.0) GO TO 110
      IF (LFVAR(2).LT.NFVAR(2)) GO TO 75
C Preserve IREC and MKEY.
      DO 20 I = 1,3
         ISAVE(3+I) = IREC(I)
         ISAVE(9+I) = MKEY(I)
   20 CONTINUE
C
C Reset LAST
      KPVLNK = ISAVE(35)
      CALL MA52ED(LAST,NDF,ISAVE(40),IW(LHED),IW(KPVLNK))

C Reset ISAVE(17), ISAVE(39), ISAVE(40) to 0
C and reset ISAVE(16) to 2 (so MA42A/AD and MA42B/BD
C may be recalled without recalling MA42I/ID).
C
      ISAVE(17) = 0
      ISAVE(39) = 0
      ISAVE(40) = 0
      ISAVE(16) = 2
C Set X
      DO 40 J = 1,NRHS
         DO 30 I = 1,NDF
            X(I,J) = ZERO
   30    CONTINUE
   40 CONTINUE
C Set ISAVE(38) to indicate call to MA52B/BD complete.
      ISAVE(38) = 3003
      GO TO 110
C
C *** Error returns ***
   50 IFLAG = -1
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) IFLAG
         WRITE (LP,FMT=9040)
      END IF
      GO TO 110
   60 IFLAG = -2
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) IFLAG
         WRITE (LP,FMT=9050)
      END IF
      GO TO 110
   70 IFLAG = -3
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) IFLAG
         WRITE (LP,FMT=9010) NFVAR(1)
      END IF
      GO TO 110
   75 IFLAG = -3
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) IFLAG
         WRITE (LP,FMT=9015) NFVAR(2)
      END IF
      GO TO 110
   80 IFLAG = -4
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) IFLAG
         WRITE (LP,FMT=9020)
      END IF
      GO TO 110
   90 IFLAG = -5
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) IFLAG
         WRITE (LP,FMT=9030) NDF
      END IF
      GO TO 110
C
  110 RETURN
 9000 FORMAT (' ***** Error return from MA52K/KD *****  IFLAG = ',I3)
 9010 FORMAT (7X,'LFVAR(1) must be at least ',I8)
 9015 FORMAT (7X,'LFVAR(2) must be at least ',I8)
 9020 FORMAT (7X,'NRHS has been changed since the last call to ',
     +       'MA42B/BD.')
 9030 FORMAT (7X,'LX too small. Must be at least ',I8)
 9040 FORMAT (7X,'An error was issued by MA42B/BD.')
 9050 FORMAT (7X,'MA42B/BD was called for all the elements in the ',
     +       'problem.')
      END
C*********************************************************************
      SUBROUTINE MA52LD(MFRONT,NFRONT,NRHS,LRHS,BUFRL,LLB,BUFRU,
     +                  LUB,IBUFR,LIBUFR,FA,FARHS,KHED,LHED,
     +                  IFILE,IREC,ISIZE,MKEY,
     +                  NUMBLK,IFVAR,JFVAR,IP,FVAR,FRHS,LFVAR,NFVAR,
     +                  ISAVE,LP,IFLAG)
C
C     Argument list. * indicates the argument is changed by the routine.
C
C  MFRONT - Integer variable.  Set by user of MA42 to max number
C           of rows in front.
C  NFRONT - Integer variable.  Set by user of MA42 to max number
C           of cols in front.
C  NRHS   - Integer variable.  Holds number of right-hand sides.
C  LRHS   - Integer variable. Second dimension of array FARHS.
C  BUFRL  - Real (DP) array of length LLB. In-core buffer holding
C           entries of L.
C  LLB    - Integer variable. Length of array BUFRL.
C  BUFRU  - Real (DP) array of length LUB. In-core buffer holding
C           entries of U.
C  LUB    - Integer variable. Length of array BUFRU.
C  IBUFR  _ Integer array of length LIBUFR. In-core buffer holding
C           integer information on factors L and U.
C  LIBUFR - Integer variable. Length of array IBUFR.
C  FA     - Real (DP) array with dimensions MFRONT,NFRONT.
C           Holds the frontal matrix.
C  FARHS   - Real (DP) array with dimensions MFRONT, LRHS.
C           Holds the right hand sides corresponding to the current
C           frontal matrix.
C  KHED   - Integer array of dimension MFRONT. Used to hold the indices
C           of variables corresponding to rows in the front.
C  LHED   - Integer array of dimension NFRONT. Used to hold the indices
C           of variables corresponding to columns in the front.
C  IFILE  - Integer array of length 3. Stream numbers for direct
C           access data sets.
C *IREC   - Integer array of length 3. Pointers to first free space
C           in buffers.
C  ISIZE  - Integer array of length 3. Length of buffers.
C *MKEY   - Integer array of length 3. Number of records written
C           to direct access data sets.
C  NUMBLK - Integer array of length 3. Number of records in
C           direct access data sets.
C *IFVAR  - Integer array of dimension LFVAR(1).
C           Not set on entry. On exit holds the indices
C           of the rows left in the front after last
C           call to MA42B/B.
C *JFVAR  - Integer array of dimension LFVAR(2).
C           Not set on entry. On exit holds the col. indices
C           of the variables left in the front after last
C           call to MA42B/B, ordered by rows.
C *IP     - Integer array length LFVAR(1)+1. Row pointer array.
C *FVAR   - Real (dp) array of dimensions LFVAR(2).
C           Not set on entry. On exit holds values for the
C           variables left in the front after last
C           call to MA42B/B.
C *FRHS   - Real (dp) array of dimensions (LFVAR(1),LRHS).
C           Not set on entry. On exit holds right-hand sides
C           for the variables left in the front after last
C           call to MA42B/B.
C  LFVAR  - Integer array length 2. Defines dimension of IFVAR,
C           JFVAR, FVAR,
C  *NFVAR  - Integer array length 2.
C            On entry holds number of rows and cols in front.
C            On exit holds number of rows and entries in front.
C   LP     - Stream for error messages.
C  *IFLAG  - Error flag. Negative values indicate a fatal error.
C            Possible nonzero values are -6,-7.
C
C  Local variables
C
C  I      - Do loop variable.
C  JFLAG  - Error flag when calling MA42L/LD.
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER IFLAG,LIBUFR,LLB,LP,LRHS,LUB,MFRONT,NFRONT,NRHS
      INTEGER LFVAR(2),NFVAR(2)
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION BUFRL(LLB),BUFRU(LUB),FA(MFRONT,NFRONT),
     +                 FARHS(MFRONT,LRHS),FRHS(LFVAR(1),LRHS),
     +                 FVAR(LFVAR(2))
      INTEGER IBUFR(LIBUFR),IFILE(3),IFVAR(LFVAR(1)),JFVAR(LFVAR(2)),
     +        IREC(3),ISAVE(31),
     +        ISIZE(3),LHED(NFRONT),MKEY(3),NUMBLK(3),
     +        KHED(MFRONT),IP(LFVAR(1)+1)
C     ..
C     .. Local Scalars ..
      INTEGER I,J,JFLAG,K

C     ..
C     .. External Subroutines ..
      EXTERNAL MA42LD
C     ..
C Jump if direct access not used.
      IF (ISAVE(31).EQ.-1) GO TO 10

      JFLAG = 0
C Write out contents of buffers to direct access data sets.

C Unless IREC(1)=1, BUFRU still contains reals to be output
C to direct access data set.
      IF (IREC(1).NE.1) THEN
C Check there is room to write.
         IF (MKEY(1)+1.GT.NUMBLK(1)) GO TO 70
         MKEY(1) = MKEY(1) + 1
C Fill-in the rest of BUFRU with garbage (so that it is
C not undefined ... can only be undefined if MKEY(1) = 1
C so that the whole buffer has not already been used)
         IF (MKEY(1).EQ.1) THEN
            DO 5 I = IREC(1),LUB
               BUFRU(I) = ZERO
   5        CONTINUE
         END IF
         CALL MA42LD(3,IFILE(1),MKEY(1),BUFRU,LUB,IBUFR,LIBUFR,
     +               LP,JFLAG)
         IF (JFLAG.LT.0) GO TO 80
      END IF

C Unless IREC(2)=1, BUFRL still contains reals to be output
C to direct access data set. If ISIZE(2) = 0 the L factor
C is not being stored.
      IF (ISIZE(2).NE.0 .AND. IREC(2).NE.1) THEN
C Check there is room to write.
         IF (MKEY(2)+1.GT.NUMBLK(2)) GO TO 70
         MKEY(2) = MKEY(2) + 1
C Fill-in the rest of BUFRL with garbage (so that it is
C not undefined ... can only be undefined if MKEY(2) = 1
C so that the whole buffer has not already been used)
         IF (MKEY(2).EQ.1) THEN
            DO 6 I = IREC(2),LLB
               BUFRL(I) = ZERO
   6        CONTINUE
         END IF
         CALL MA42LD(3,IFILE(2),MKEY(2),BUFRL,LLB,IBUFR,LIBUFR,
     +               LP,JFLAG)
         IF (JFLAG.LT.0) GO TO 80
      END IF

C Check there is room to write.
      IF (MKEY(3)+1.GT.NUMBLK(3)) GO TO 70
C Set a flag in IBUFR to indicate end of file (IBUFR(IREC(3))=0).
   10 IBUFR(IREC(3)) = 0
C Fill-in the rest of IBUFR with garbage (so that it is
C not undefined)
      DO 15 I = IREC(3)+1,LIBUFR
         IBUFR(I) = -1
   15 CONTINUE
      IREC(3) = IREC(3) + 1
      MKEY(3) = MKEY(3) + 1
      CALL MA42LD(-3,IFILE(3),MKEY(3),BUFRU,LUB,IBUFR,LIBUFR,
     +            LP,JFLAG)
      IF (JFLAG.LT.0) GO TO 80
C
C Set IFVAR, JFVAR, FVAR and FRHS.

C We have to count number of nonzero entries left in FA
C If LFVAR(2) exceeds NFVAR(1)*NFVAR(2) then we know there
C is room in JFVAR and FVAR. If not, we will have to check space
      K = 0
      IP(1) = 1
      IF (LFVAR(2).GE.NFVAR(1)*NFVAR(2)) THEN
         DO 40 I = 1,NFVAR(1)
           DO 30 J = 1,NFVAR(2)
             IF (FA(I,J).NE.ZERO) THEN
               K = K + 1
               FVAR(K) = FA(I,J)
               JFVAR(K) = LHED(J)
             END IF
   30     CONTINUE
          IP(I+1) = K + 1
C Check we do not have a row with no entries. If we
C do we will store the first zero explicitly.
          IF (IP(I+1).EQ.IP(I)) THEN
            K = K + 1
            FVAR(K) = ZERO
            JFVAR(K) = LHED(1)
            IP(I+1) = K+1
          END IF
   40   CONTINUE
      ELSE
         DO 45 I = 1,NFVAR(1)
           IF (K+NFVAR(2).LE.LFVAR(2)) THEN
C We have room for next row, even if all entries are nonzero
             DO 35 J = 1,NFVAR(2)
               IF (FA(I,J).NE.ZERO) THEN
                 K = K + 1
                 FVAR(K) = FA(I,J)
                 JFVAR(K) = LHED(J)
               END IF
   35        CONTINUE
           ELSE
             DO 36 J = 1,NFVAR(2)
               IF (FA(I,J).NE.ZERO) THEN
                 K = K + 1
                 IF (K.LE.LFVAR(2)) THEN
                   FVAR(K) = FA(I,J)
                   JFVAR(K) = LHED(J)
                 END IF
               END IF
   36        CONTINUE
           END IF
           IP(I+1) = K + 1
           IF (IP(I+1).EQ.IP(I)) THEN
             K = K + 1
             FVAR(K) = ZERO
             JFVAR(K) = LHED(1)
             IP(I+1) = K+1
           END IF
   45   CONTINUE
      END IF

C Set NFVAR(2) to number of entries in front
      NFVAR(2) = K

      DO 20 I = 1,NFVAR(1)
         IFVAR(I) = KHED(I)
   20 CONTINUE

      DO 60 J = 1,NRHS
         DO 50 I = 1,NFVAR(1)
            FRHS(I,J) = FARHS(I,J)
   50    CONTINUE
   60 CONTINUE
      GO TO 90
C
C **** Error returns ****
   70 IFLAG = -6
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) IFLAG
         IF (MKEY(1)+1.GT.NUMBLK(1)) WRITE (LP,FMT=9010) ISIZE(1),
     +       ISIZE(1)* (MKEY(1)+1)
         IF (MKEY(2)+1.GT.NUMBLK(2)) WRITE (LP,FMT=9020) ISIZE(2),
     +       ISIZE(2)* (MKEY(2)+1)
         IF (MKEY(3)+1.GT.NUMBLK(3)) WRITE (LP,FMT=9030) ISIZE(3),
     +       ISIZE(3)* (MKEY(3)+1)
      END IF
      GO TO 90
   80 IFLAG = -7
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) IFLAG
      END IF
   90 RETURN
C
 9000 FORMAT (' ***** Error return from MA52K/KD *****  IFLAG = ',I3)
 9010 FORMAT (7X,'Insufficent storage for U-factor.',/7X,'With present',
     +       ' LENBUF(1) size of ',I8,/7X,'LENFLE(1) must be at least ',
     +       I8)
 9020 FORMAT (7X,'Insufficent storage for L-factor.',/7X,'With present',
     +       ' LENBUF(2) size of ',I8,/7X,'LENFLE(2) must be at least ',
     +       I8)
 9030 FORMAT (7X,'Insufficent storage for integers.',/7X,'With present',
     +       ' LENBUF(3) size of ',I8,/7X,'LENFLE(3) must be at least ',
     +       I8)
      END
* COPYRIGHT (c) 1989 AEA Technology
* Original date 19 Oct 1992
C       Toolpack tool decs employed.
C 3/7/98. A changed to be assumed-size.
C 13/3/02 Cosmetic changes applied to reduce single/double differences
C
C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC46AD(NR,NC,LA,YESA,INDEX,A,IPR,JPC,IW,INFORM)
C
C   PURPOSE:
C   --------
C   TO SORT THE PATTERN OF A SPARSE MATRIX A FROM ROW ORDER
C   TO COLUMN ORDER, UNORDERED WITHIN EACH COLUMN (OR VICE VERSA).
C
C
C  *********************************************************************
C
C   ARGUMENT LIST:
C   ---------------
C
C   NR - INTEGER VARIABLE
C      - ON ENTRY MUST BE SET TO THE NUMBER OF ROWS IN THE MATRIX
C      - UNCHANGED ON EXIT
C
C   NC - INTEGER VARIABLE
C      - ON ENTRY MUST BE SET TO THE NUMBER OF COLUMNS IN THE MATRIX
C      - UNCHANGED ON EXIT
C
C   LA  - INTEGER VARIABLE
C       - ON ENTRY, MUST BE SET TO A VALUE AT LEAST AS LARGE AS
C       - THE NUMBER OF NONZEROS IN THE MATRIX
C       - UNCHANGED ON EXIT
C
C   YESA - LOGICAL VARIABLE
C        - IF YESA IS SET TO .FALSE., THE ARRAYS A ARE NOT ACCESSED
C        - UNCHANGED ON EXIT
C
C   A    - REAL (DOUBLE PRECISION) ASSUMED-SIZE ARRAY.
C        - IF YESA IS .FALSE., A NEED NOT BE SET AND IS NOT ALTERED.
C        - IF YESA IS .TRUE.:
C        - ON ENTRY SET TO CONTAIN THE NONZEROS OF THE MATRIX A
C        - STORED BY ROWS. ALL OF THE NONZEROS IN ROW I
C        - MUST APPEAR BEFORE ROW I + 1 (0 < I < NR).
C        - ON EXIT SET TO CONTAIN THE NONZEROS OF THE MATRIX A
C        - STORED BY COLUMN. ALL OF THE NONZEROS IN COLUMN J
C        - MUST APPEAR BEFORE COLUMN J + 1 (0 < J < NC).
C
C   INDEX - INTEGER ARRAY OF LENGTH LA.
C         - ON ENTRY SET TO CONTAIN THE COLUMN INDICES OF THE NONZEROS
C         - OF THE MATRIX A STORED BY ROWS. THE COLUMN INDEX OF THE
C         - NONZERO STORED IN A( K ) MUST APPEAR IN INDEX( K ).
C         - ON EXIT SET TO CONTAIN THE ROW INDICES OF THE NONZEROS
C         - OF THE MATRIX A STORED BY COLUMNS. THE ROW INDEX OF THE
C         - NONZERO STORED IN A( K ) WILL APPEAR IN INDEX( K ).
C
C   IPR - INTEGER ARRAY OF LENGTH NR + 1
C       - ON ENTRY, IPR( I ) MUST CONTAIN THE POSITION IN INDEX
C       - (AND A) OF THE FIRST ENTRY IN ROW I (I=1,...,NR)
C       - IPR( NR + 1 ) MUST BE SET TO 1 + THE NUMBER OF NONZEROS IN A
C       - UNALTERED ON EXIT.
C
C   JPC - INTEGER ARRAY OF LENGTH NC + 1
C       - NEED NOT BE SET ON ENTRY.
C       - ON EXIT, JPC( J ) CONTAINS THE POSITION IN INDEX (AND A)
C       - OF THE FIRST ENTRY IN COLUMN J (J=1,...,NC).
C       - JPC( NC + 1 ) IS SET TO 1 + THE NUMBER OF NONZEROS IN A.
C
C   IW - INTEGER ARRAY OF LENGTH NC.
C      - NOT SET ON ENTRY.
C      - THE ARRAY IS USED AS WORKSPACE.
C
C   INFORM - INTEGER VARIABLE.
C          - NOT SET ON ENTRY.
C          - ON EXIT, INFORM =  0, INDICATES A SUCCESSFUL CALL.
C          -          INFORM = -1, INDICATES NC<1,NR<1 OR LA<NNZ.
C          -          INFORM = -2, INDICATES THAT A COLUMN WITH INDEX
C          -                       SMALLER THAN 1 OR LARGER THAN NC HAS
C          -                       BEEN SPECIFIED.
C
C  ******************************************************************
C
C  NICK GOULD, HARWELL.
C  10TH JANUARY 1989.

C     .. Scalar Arguments ..
      INTEGER INFORM,LA,NC,NR
      LOGICAL YESA
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(*)
      INTEGER INDEX(LA),IPR(NR+1),IW(NC),JPC(NC+1)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AENTRY,ANEW
      INTEGER INEXT,IROW,IROWP1,ISOFAR,ISTART,J,JCOL,JPCJ,K,NEWCOL,
     +        NEWROW,NNZ
C     ..
C
C  CHECK THAT THE INPUT PARAMETERS ARE LARGE ENOUGH.
C
      NNZ = IPR(NR+1) - 1
      IF (NR.LE.0 .OR. NC.LE.0 .OR. LA.LT.NNZ) THEN
         INFORM = -1
         RETURN
C
      END IF
C
C  CLEAR JPC.
C
      DO 10 J = 2,NC + 1
         JPC(J) = 0
   10 CONTINUE
C
C  COUNT THE NUMBER OF ELEMENTS IN EACH COLUMN AND STORE IN JPC.
C  NEGATE THE COLUMN NUMBERS IN INDEX, SO THAT THEY ARE FLAGGED
C  AS INDEX IS GRADUALLY OVERWRITTEN BY ROW NUMBERS.
C
      DO 20 K = 1,NNZ
         J = INDEX(K)
         IF (J.LT.1 .OR. J.GT.NC) THEN
            INFORM = -2
            RETURN
C
         END IF
C
         INDEX(K) = -J
         JPC(J+1) = JPC(J+1) + 1
   20 CONTINUE
C
C  COMPUTE THE STARTING ADDRESSES FOR EACH COLUMN OF A WITHIN
C  INDEX AND A.
C
      JPC(1) = 1
      DO 30 J = 2,NC + 1
         JPC(J) = JPC(J) + JPC(J-1)
   30 CONTINUE
C
C  STORE IN IW THE ROW NUMBER OF THE ENTRY OF A WHOSE COLUMN NUMBER
C  IS THE JPC( J )-TH ENTRY OF INDEX.
C
      ISOFAR = 0
      DO 60 J = 1,NC
         JPCJ = JPC(J)
         DO 40 IROWP1 = ISOFAR + 1,NR
            IF (JPCJ.LT.IPR(IROWP1)) GO TO 50
   40    CONTINUE
         IROWP1 = NR + 1
   50    CONTINUE
         ISOFAR = IROWP1 - 1
         IW(J) = ISOFAR
   60 CONTINUE
C
C  REORDER THE ELEMENTS INTO COLUMN ORDER.
C  FILL IN EACH COLUMN FROM THE FRONT. AS A NEW ENTRY IS PLACED
C  IN COLUMN K INCREASE THE POINTER JPC( K ) BY ONE AND FIND
C  THE NEW ROW, IW( K ), THAT CORRESPONDS TO THE COLUMN NOW
C  POINTED TO BY JPC( K ).
C
      DO 120 J = 1,NC
C
C  DETERMINE THE NEXT UNPLACED ENTRY, ISTART, IN INDEX.
C
   70    CONTINUE
         ISTART = JPC(J)
C
C  SEE IF ALL THE ELEMENTS IN COLUMN J HAVE BEEN ASSIGNED.
C
         IF (ISTART.EQ.JPC(J+1)) GO TO 120
         IF (INDEX(ISTART).GT.0) GO TO 120
C
C  EXTRACT THE ROW AND COLUMN NUMBERS (AND VALUE) OF THE
C  STARTING ELEMENT.
C
         IROW = IW(J)
         JCOL = -INDEX(ISTART)
         IF (YESA) AENTRY = A(ISTART)
C
C  MOVE ELEMENTS IN A CYCLE, ENDING BACK AT COLUMN J.
C
         DO 100 K = ISTART,NNZ
C
C  FIND THE FIRST EMPTY LOCATION IN COLUMN JCOL IN INDEX.
C
            INEXT = JPC(JCOL)
C
C  EXTRACT THE ROW NUMBER OF THE NEXT ELEMENT.
C
            NEWROW = IW(JCOL)
C
C  UPDATE JPC(JCOL), FIND THE NEW ROW NUMBER AND STORE IT IN IW(JCOL).
C
            JPCJ = INEXT + 1
            JPC(JCOL) = JPCJ
            DO 80 IROWP1 = NEWROW + 1,NR
               IF (JPCJ.LT.IPR(IROWP1)) GO TO 90
   80       CONTINUE
            IROWP1 = NR + 1
   90       CONTINUE
            IW(JCOL) = IROWP1 - 1
C
C  IF THE ENTRY BELONGS IN THE J-TH COLUMN, THE CYCLE IS COMPLETE.
C
            IF (JCOL.EQ.J) GO TO 110
C
C  EXTRACT THE COLUMN NUMBER (AND VALUE) OF THE NEXT ELEMENT.
C
            NEWCOL = -INDEX(INEXT)
            IF (YESA) ANEW = A(INEXT)
C
C  STORE THE ROW NUMBER (AND VALUE) OF THE CURRENT ELEMENT.
C
            INDEX(INEXT) = IROW
            IF (YESA) A(INEXT) = AENTRY
C
C  MAKE THE NEXT ELEMENT INTO THE CURRENT ONE.
C
            IROW = NEWROW
            JCOL = NEWCOL
            IF (YESA) AENTRY = ANEW
  100    CONTINUE
C
C  THE CYCLE IS COMPLETE.
C
  110    CONTINUE
C
C  STORE THE ROW NUMBER (AND VALUE) OF THE STARTING ELEMENT.
C
         INDEX(ISTART) = IROW
         IF (YESA) A(ISTART) = AENTRY
         GO TO 70
C
  120 CONTINUE
C
C  REVISE JCP TO POINT TO THE START OF EACH COLUMN OF A.
C
      DO 130 J = NC,2,-1
         JPC(J) = JPC(J-1)
  130 CONTINUE
      JPC(1) = 1
      INFORM = 0
      RETURN
      END
* COPYRIGHT (c) 1992 Council for the Central Laboratory
*                    of the Research Councils
* Original date 2 Apr. 1993
C****************************************************************
C Exploit zeros in front by setting ISAVE(20) > 1
C (HSL 12 version equivalent to ISAVE(20) = 1)
C This version does row and col. swaps to exploit zeros
C in both rows and columns.
C To exploit zeros in front, set ISAVE(20)>1 before first call to
C MA42B/BD.
C Also, PIVBLK = ISAVE(19) in MA42F/FD allows us to wait
C until sufficiently many variables can be eliminated at once
C (as in cache paper). PIVBLK ALSO USED IN MA42J/JD.
C In HSL 12, PIVBLK = 1. If block size other than 1 is wanted,
C user must set ISAVE(19) to min. pivot block before
C calls to MA42J/JD (no action needed for min. block size of 1)
C but in experiments found  16 or 32 better (16 used in MA62)
C This version of code should be used with MP42 (which allows user to
C choose min. pivot block size)
C Note that for equation problems, 16 may be too large.
C USER MUST NOT CHANGE ISAVE(19) BETWEEN CALLS
C TO MA42J/JD AND CALLS TO MA42B/BD (IE USER CAN ONLY RESET
C ISAVE(19) AND ISAVE(20) ONCE) AND ISAVE(20)  MUST NOT BE CHANGED
C BETWEEN CALLS TO MA42B/BD
C****************************************************************
C Minor changes made to code 18 Nov. 1993.
C A bug was found in the case of the error flag INFO(1)=-14 being
C returned (KPVLNK was not properly set in this case and this caused
C unassigned variable error messages to be returned when LAST was
C restored).
C Also, ISAVE(38) set on each return from MA42A/AD, MA42J/JD, MA42P/PD,
C MA42B/BD to hold a copy of INFO(1). On each entry, a check is made
C that an error was not previously encountered.
C
C 4 March 1994 bug found in MA42P/PD. If ISTRM(2)=0 then we need to
C assign NUMBLK(2)=0 (we had NUMBLK(2) unassigned in this case).
C 15 April 1994 :error return -18 not possible from MA42J/JD so
C test for this removed from MA42B/BD. Also, before calling MA42J/JD
C from MA42B/BD for the first time, it is not necessary to move
C data to the start of IW so this has been removed).
C 18 April 1994: argument list for MA42D/DD changed to pass lengths
C of the array Y as LY1 and LY2. Y is of length 1,1 when called from
C MA42B/BD.
C 27 April 1994: change to MA42B/BD so that ISAVE(1)-ISAVE(15) are set
C           on first call (unless MA42P/PD has been called).
C           Want to be able to recall MA42B/BD without recalling
C           MA42P/PD.
C 13 May 1994: -15 cannot be returned by MA42O/OD or by MA42H/HD.
C 20 May 1994: INFO(1)=3 replaced by INFO(1)=2.
C 2 June 1994: NUMPIV made a local variable in MA42O/OD
C         KX made a local variable in MA42N/ND
C         In MA42J/JD, IFSIZE(1) = MAX(INFO(8),KFRNT)
C                      IFSIZE(2) = MAX(INFO(9),LFRNT)
C        (to stop code return a lower bound for the frontsize which is
C         smaller than that provided by user).
C         Changed printed error message when -12 returned.
C 13 June 1994: Immediate return added to MA42B/BD if it is entered with
C          INFO(1)<0.
C 24 June 1994: Bug found in MA42N/ND. We do not have the best pivot on
C          hold if IFORCE>0 and IELIM=1
C 1 Feb 1995: Operation count (RINFO(2)) changed since BLAS
C             cannot take advantage of zeros in frontal matrix.
C 17 March 1995   CLOSE (IFILE(I)) replaced by
C          CLOSE (IFILE(I),STATUS='DELETE')
C 6 August 1995  Ensured INFO(2) = 0 and RINFO(1) = 0.0 if matrix
C          found to be singular (this was not always happening
C          if the computation continued)
C 10 August 1995 Bug in MA42B/BD. IF (ISAVE(31).EQ.31) should read
C          IF (ISAVE(31).EQ.1)
C 14 August 1995 Bug corrected in computation of sign of
C          the determinant.
C 4 October 1995. Filled in buffers on the last time they are
C          written out... June 1996 changed fill in to avoid
C          doing it unnecessarily. We now only do it
C          if the buffer is the first and last one to
C          be written to d.a. file.
C 17 April 1996. Moved statement OFDIAG = 0 from MA42N/ND into
C          MA42F/FD (if the problem has only one element,
C          all variables are static condensation variables
C          and OFDIAG would not be set).
C 4  June 1996.  In MA42B/BD, avoid zeroing the solution vector X
C          if the number of variables in the problem is
C          equal to the largest integer used to index a variable
C          (i.e if NDF = INFO(3)). N.B. Cannot avoid the zeroing
C          in MA42C/CD since we do not pass the number of
C          variables to MA42C/CD in an ISAVE entry (clearly, we should
C          have done!). July 2000 : still have to zero X if
C          problem found to be singular, otherwise components
C          of X may be undefined.
C 11 Nov 1996 After each call to MA42D/DD and MA42E/ED, check
C          the error flag and write error message if appropriate.
C 20 Jan 1997 In MA42G/GD and MA42H/HD, JFLAG should only
C          be checked for an error IF MA42L/LD has been called
C          (o.w. could be undefined). This has been corrected.
C 14 March 1997 In MA42F/FD added a test so that MA42N/ND is only
C          called if some variables have become fully summed since
C          the assembly of the most recent elt/equ.
C 4 Sept. 1997  In MA42J/JD, if called from within MA42B/BD
C          (because frontsize too small), use INFO values to
C          initialise IFSIZE and allow updated lower bounds on the
C          filesizes to be returned to the user.
C 27 October 1998. In MA42J/JD and MA42B/BD changed test on NMAXE
C          so that NMAXE does not have to have the same
C          value on each entry if elements are being used
C          (still check that user has not changed for element
C          to equation entry, or visa versa, by comparing NMAXE with
C          ISAVE(21)
C          Changed  IF (NMAXE.NE.ISAVE(21)) GO TO 170
C          to   IF (ISAVE(21).EQ.1 .AND. NMAXE.GT.1) GO TO 170
C               IF (ISAVE(21).GT.1 .AND. NMAXE.EQ.1) GO TO 170
C          Similar change in MA42J/JD.
C 12 Dec. 1998. INFO(3) is now updated in MA42F/FD (not in MA42H/HD).
C          Change made so that when MA42 used with MA52, INFO(3)
C          contains correct info. on exit (otherwise, as we do not
C          pick all variable as pivots within a subdomain, INFO(3)
C          would not be a count of all variables in subdomain but
C          would hold number of variables eliminated within subdomain)
C 18 December 1998. In MA42F/FD changed
C         INFO(10) = MKEY(1)
C         INFO(11) = MKEY(2)
C         INFO(12) = MKEY(3)
C to
C         INFO(10) = MAX(1,MKEY(1))
C         INFO(11) = MAX(1,MKEY(2)) ... this only if L factored stored
C         INFO(12) = MAX(1,MKEY(3))
C   (since otherwise can get INFO(10)-INFO(12) equal to zero when error
C    -17 return ... INFO values then no use for resetting LENFLE)
C 14 April 1999 In MA42C/CD, bug found in workspace LW if direct
C          access files used.
C 28 April 1999 Error in flop count for static condensation variables
C          corrected. In MA42O/OD
C          OPS = OPS + DBLE(NUMPIV* (NVAR-1)* (MVAR-1)*2)
C          changed to
C          DO 456 J = 1,NUMPIV
C              OPS = OPS + DBLE((NVAR-J)* (MVAR-J)*2)
C    456   CONTINUE
C          NOTE: has no effect for equation entry (since NUMPIV=1)
C 20 May 1999  Error found in MA42D/DD.
C              Replace
C               DO 10 I = 1,IRECD
C                  IBUFR(DIMIBF-IRECD+I) = IBUFR(IR2+I)
C   10          CONTINUE
C              by
C               DO 10 I = IRECD,1,-1
C                  IBUFR(DIMIBF-IRECD+I) = IBUFR(IR2+I)
C   10          CONTINUE
C              and
C               DO 30 I = 1,IPNT
C                  BUFR(DIMBUF-IPNT+I) = BUFR(JR2+I)
C   30          CONTINUE
C              by
C               DO 30 I = IPNT,1,-1
C                  BUFR(DIMBUF-IPNT+I) = BUFR(JR2+I)
C   30          CONTINUE
C 6 August 1999 Error in MA42D/DD. If direct access files not used
C               then JFLAG not set so at end of subroutine jump
C               to return (in MA42E/ED jump already there.)
C 29 November 1999. Error in call to MA42H/HD from MA42N/ND when
C KPRE>0 or LFRE>0.
C Also, changed second dimension of FA in MA42G/GD from NFRONT
C to * (NFRONT causes error when MA42G/GD called from MA42N/ND
C when KPRE>0 or LFRE>0). Similarly, second dimension of FRHS
C changed to *
C 25 March 2002. If INFO(1) = 5 or 6 is returned by MA42K/KD, it
C is necessary to zero the appropraite entries of FA as if the data
C has been output. MA42G/GD, MA42H/HD and MA42K/KD modified.
C
C Changes made to incorporate the option of forcing diagonal pivots
C until the last element has been assembled.
C ICNTL(7)  used to control whether diagonal pivoting is forced.
C ICNTL(7) =  0 (default) Static condensation + off-diag. pivoting
C ICNTL(7) =  1001 Static condensation + diag. pivoting
C ICNTL(7) = -1001 No static condensation + diag. pivoting
C ALL other values of ICNTL(7):
C             No static condensation + off-diag. pivoting.
C Diagonal pivoting only allowed for ELEMENT entry.
C If diagonal pivoting and no warnings issued, on exit
C INFO(1) set to 6 + (Number of off-diagonal pivots)
C ISAVE(38) set to hold abs(ICNTL(7)) UNLESS an error is
C encountered, in which case it holds INFO(1).
C
C MA42N/ND, MA42O/OD altered to cope with choosing diagonal pivots.
C
C Also made a change to pivot search (on or off diagonal).
C In MA42N/ND  added a loop to look only at largest entry in
C column to to see if it can be used as a pivot
C (trying to reduce search time by not looking for largest entry in
C fully summed part unless we have to). This may mean INFO(13),
C INFO(14), and INFO(15) having different values than for previous
C version of code. In particular, some changes to test deck output.
C
C
C *** Reference MA42 suite ***
C *** Any problems contact Iain S. Duff  or Jennifer A. Scott
C     at Atlas Centre, Rutherford Appleton Laboratory ***
C *** Although every effort has been made to ensure robustness and
C     reliability of the subroutines in this MA42 suite,  we
C     disclaim any liability arising through the use or misuse of
C     any of the subroutines in the MA42 suite ***

C 12th July 2004 Version 1.0.0. Version numbering added.
C 18th Jan. 2005 Version 1.1.0. In MA42N/ND changed
C IF (NMAXE.EQ.1 .OR. LAST(MFR).GE.0) THEN
C into two separate tests (problem if NMAXE = 1 and system
C is singular, MFR may be too large)
C 26th April 2007 Version 1.2.0.
C     Statements that extended to column 73 changed.


      SUBROUTINE MA42ID(ICNTL,CNTL,ISAVE)
C This subroutine initializes control parameters.
C
C     Argument list. * indicates the argument is changed by the routine.
C
C  *ICNTL  -Integer array of length 8.
C ICNTL(1) - Stream number for error messages. On exit
C            set to the default value 6.
C ICNTL(2) - Stream number for warnings. On exit
C            set to the default value 6.
C ICNTL(3) - number of bytes for a real record.
C            On exit given default value
C            4 for single precision and 8 for double precision.
C ICNTL(4) - number of bytes for an integer record.
C            On exit given default value 4.
C ICNTL(5) - If ICNTL(5) is positive then, when the number
C            of potential pivot columns is greater than or equal to
C            ICNTL(5), an elimination will be performed even if the
C            best pivot candidate does not satisfy the threshold
C            criterion determined by CNTL(2).
C            On exit set to default value 0.
C ICNTL(6) - If ICNTL(6) and ICNTL(5) are positive
C            then only ICNTL(6) potential pivot
C            columns will be searched for a pivot.
C            The best pivot candidate (largest relative
C            to other nonzeros in its column) from these ICNTL(6)
C            columns will then be used as a pivot.
C            On exit set to default value 0.
C ICNTL(7)  - controls static condensation and diagonal
C            pivoting. On exit set to default value 0.
C ICNTL(7) =  0 (default) Static condensation + off-diag. pivoting
C ICNTL(7) =  1001 Static condensation + diag. pivoting
C ICNTL(7) = -1001 No static condensation + diag. pivoting
C             ALL other values of ICNTL(7):
C             No static condensation + off-diag. pivoting.
C            On exit set to default value 0.
C ICNTL(8) - controls action if singularity detected
C            in matrix. If ICNTL(8) is 0, exit with
C            fatal error immediately singularity detected. If
C            ICNTL(8) is non-zero, a warning
C            is issued when singularity detected but computation
C            proceeds. On exit set to default value 0.
C *CNTL    _ Real (DP) array of length 2. Holds control
C            parameters.
C  CNTL(1) - The matrix will be declared singular if the element of
C            largest absolute value in any column during
C            decomposition has absolute value less than or equal
C            CNTL(1). On exit set to default value 0.0.
C  CNTL(2) - An element of the frontal matrix will only be considered
C            suitable for use as a pivot if it is of absolute value
C            at least as large as CNTL(2) times the element of largest
C            absolute value in its column.
C            On exit set to default value 0.1.
C  *ISAVE  - Integer array of length 45.
C            On successful exit ISAVE(I) = 0 with the following
C            exceptions
C            ISAVE(4) - Initialised to 1.
C            ISAVE(5) - Initialised to 1.
C            ISAVE(6) - Initialised to 1.
C            In subsequent routines, ISAVE(3+I) is equivalent
C            to the local variable IREC(I) (I = 1,2,3).
C            ISAVE(16) - Initialised to 2. ISAVE(16) is used by
C            MA42A/AD to determine whether or not a call is the
C            first call to MA42A/AD. Also used by MA42J/JD.
C            ISAVE(31) - Initialised to -1. If MA42P/PD is called,
C            ISAVE(31) will be changed to +1 (i.e. used as a flag to
C            indicate whether MA42P/PD has been called by the user.
C            This is needed in MA42B/BD and MA42C/CD).
C            ISAVE(19), ISAVE(20) initialised to 1.
C            ISAVE(19) is min. pivot block size
C            ISAVE(20) controls whether zeros exploited in front.
C
C   Local variables
C   I      - Integer variable. Do loop variable.
C
C     .. Array Arguments ..
      DOUBLE PRECISION CNTL(2)
      INTEGER ICNTL(8),ISAVE(45)
C     ..
C     .. Local Scalars ..
      INTEGER I
C     ..
C Set control parameters
      ICNTL(1) = 6
      ICNTL(2) = 6
      ICNTL(3) = 8
      ICNTL(4) = 4
      ICNTL(5) = 0
      ICNTL(6) = 0
      ICNTL(7) = 0
      ICNTL(8) = 0

      CNTL(1) = 0.0D0
      CNTL(2) = 0.1D0

C Initialise ISAVE
      DO 10 I = 1,45
         ISAVE(I) = 0
   10 CONTINUE
      DO 20 I = 1,3
         ISAVE(3+I) = 1
   20 CONTINUE
      ISAVE(16) = 2
      ISAVE(31) = -1
      ISAVE(19) = 1
      ISAVE(20) = 1
C
      RETURN
      END
C*******************************************************************
      SUBROUTINE MA42AD(NVAR,IVAR,NDF,LAST,LENLST,ICNTL,ISAVE,INFO)
C
C *** Reference MA42 suite ***
C *** Copyright Rutherford Appleton Laboratory  November 1992 ***
C *** Any problems contact Iain S. Duff  or Jennifer A. Scott
C     at Atlas Centre, Rutherford Appleton Laboratory ***
C *** Although every effort has been made to ensure robustness and
C     reliability of the subroutines in this MA42 suite,  we
C     disclaim any liability arising through the use or misuse of
C     any of the subroutines in the MA42 suite ***
C
C This subroutine determines, in array LAST (of length LENLST), the
C element or equation in which each variable appears for the last
C time.  This information is required by MA42J/JD and MA42B/BD
C so that they can tell when a variable is fully summed and
C can be used as pivot.
C
C     Argument list. * indicates the argument is changed by the routine.
C
C   NVAR   Integer variable.  Must be set by the user to hold the
C          number of variables in the current elt/eqn.
C   IVAR   Integer array of length NVAR which must be set by the
C          user to hold the global indices of variables in
C          current elt/eqn.
C  *NDF    Integer variable.  On output NDF is set to the maximum index
C          used for a variable.
C  *LAST   Integer array of length LENLST.  On output, LAST(I) holds
C          the last elt/eqn in which variable I appears. If there
C          is no variable with index J (J.LE.NDF), LAST(J)=0.
C   LENLST Integer variable.  Length of array LAST.
C   ICNTL  Integer array of length 8.
C          ICNTL(1) holds stream number for error messages.
C          Other entries not accessed by routine.
C  *ISAVE  Integer array of length 45. Must be unchanged
C          since call to MA42I/ID.
C          On exit only ISAVE(16), ISAVE(18), ISAVE(30),
C          ISAVE(38) have been changed
C          ISAVE(16) - Initialised by MA42I/ID to 2. On entry to
C          MA42A/AD for the first time, ISAVE(16) is set to -2
C          and is unchanged on subsequent calls.
C          ISAVE(18) - Number of elements (or equations) so far.
C          ISAVE(30) - Holds a copy of NDF.
C          ISAVE(38) holds a copy of INFO(1).
C  *INFO   Integer array of length 23.
C          INFO(1) is used for error returns and need not be set
C          by the user before the first call to
C          MA42A/AD for a particular problem. On successful exit,
C          INFO(1)=0.
C          Possible negative values on exit are -1, -2, -3, -4.
C          INFO(I), I = 2,3,..., 25 not accessed by MA42A/AD.
C
C   Local variables
C
C    I     Integer variable. Do loop variable.
C    JVAR  Integer variable. Used to hold the global index
C          of I th variable in the current element/equation.
C    K     Integer variable. Do loop variable.
C    LP    Integer variable. Holds copy of ICNTL(1) (stream number
C          for error messages).
C
C     .. Scalar Arguments ..
      INTEGER LENLST,NDF,NVAR
C     ..
C     .. Array Arguments ..
      INTEGER ICNTL(8),INFO(23),ISAVE(45),IVAR(NVAR),LAST(LENLST)
C     ..
C     .. Local Scalars ..
      INTEGER I,JVAR,K,LP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX
C     ..
      LP = ICNTL(1)
      IF (ISAVE(16).EQ.2) THEN
C Initialization.  Only performed on first entry.
         NDF = 0
         INFO(1) = 0
         ISAVE(38) = ABS(ICNTL(7))
         ISAVE(16) = -2
         ISAVE(18) = 0
         IF (LENLST.LE.0) GO TO 30
         DO 10 I = 1,LENLST
            LAST(I) = 0
   10    CONTINUE
      END IF
C Jump because previous entry has caused error return.
      IF (ISAVE(38).LT.0) GO TO 100
C
C ISAVE(18) indicates elt/equ currently being processed.
      ISAVE(18) = ISAVE(18) + 1
      IF (NVAR.LE.0) GO TO 40
      DO 20 I = 1,NVAR
         JVAR = IVAR(I)
         IF (JVAR.LT.1 .OR. JVAR.GT.LENLST) GO TO 50
         NDF = MAX(JVAR,NDF)
         IF (LAST(JVAR).EQ.ISAVE(18)) GO TO 60
         LAST(JVAR) = ISAVE(18)
   20 CONTINUE
      ISAVE(30) = NDF
      GO TO 90
C
C **** Error returns ****
   30 INFO(1) = -1
      ISAVE(18) = 1
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1),ISAVE(18)
         WRITE (LP,FMT=9010) LENLST
      END IF
      GO TO 90
   40 INFO(1) = -2
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1),ISAVE(18)
         WRITE (LP,FMT=9020) NVAR
      END IF
      GO TO 90
   50 INFO(1) = -3
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1),ISAVE(18)
         WRITE (LP,FMT=9030) I,JVAR
      END IF
      GO TO 90
   60 INFO(1) = -4
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1),ISAVE(18)
         WRITE (LP,FMT=9040) JVAR
         DO 70 K = 1,I
            IF (IVAR(K).EQ.JVAR) GO TO 80
   70    CONTINUE
   80    WRITE (LP,FMT=9050) K,I
      END IF
C
   90 IF (INFO(1).LT.0) ISAVE(38) = INFO(1)

  100 RETURN
 9000 FORMAT (' ***** Error return from MA42A/AD *****  INFO(1) = ',I3,
     +       /' ***** Error discovered after input of elt/eqn *****',I8)
 9010 FORMAT (7X,'Length (LENLST) of array LAST is ',I8)
 9020 FORMAT (7X,'Number of variables in elt/eqn is ',I8)
 9030 FORMAT (7X,'Variable ',I8,' in elt/eqn has value ',I8)
 9040 FORMAT (7X,'More than one occurrence of variable ',I8)
 9050 FORMAT (7X,'Dual occurrences in positions ',I8,' and ',I8)
      END
C**********************************************************************
      SUBROUTINE MA42JD(NVAR,IVAR,NDF,LAST,NMAXE,IFSIZE,ICNTL,ISAVE,
     +                  INFO)
C
C This subroutine performs symbolic assembly and elimination only and
C assumes that variables can be eliminated as soon as they are fully
C summed.  It is called by MA42B/BD if the full assembly and
C elimination has failed because of insufficient space allocated to the
C front matrix.
C
C     Argument list. * indicates the argument is changed by the routine.
C
C  NVAR   - Integer variable. Number of variables in elt/eqn being
C           input.
C  IVAR   - Integer array of length NVAR. Contains indices of variables
C           in elt/eqn being input.
C  NDF    - Integer variable. Total number of variables
C           (number of degrees of freedom).
C *LAST   - Integer array of length NDF.  Must be unchanged
C           since last call to MA42A/AD or since the last call
C           to MA42B/BD (if MA42J/JD is being called by MA42B/BD).
C           For a variable not currently in the front, LAST(I) is the
C           elt/eqn in which variable I appears for the last
C           time and for a (non-fully summed) variable in the front
C           LAST(I) is replaced by -LAST(I). On exit from the final
C           call, LAST will be unchanged since last call to MA42A/AD.
C  NMAXE  - Integer variable. Must exceed 1 for element entry
C           and must be set to 1 for equation entry.
C *IFSIZE - Integer array of length 5. On exit, IFSIZE(1) and IFSIZE(2)
C           hold lower bounds on the frontsizes
C           required for successful symbolic elimination.
C           If the code is called directly by the user,
C           IFSIZE(3), IFSIZE(4), and IFSIZE(5) hold the U, L, and
C           integer file sizes required for successful factorization
C           (assuming a variable may be eliminated as soon as it is
C           fully summed). Note for MA42B/BD the file for U must also
C           allow for right-hand sides so the estimate IFSIZE(3)
C           must be increased by NRHS*NDF. If the code is called from
C           within MA42B/BD, IFSIZE(3), IFSIZE(4), and IFSIZE(5)
C           are not accessed.
C  ICNTL  - Integer array of length 8. Only ICNTL(1) (the stream
C           number for error messages) is accessed.
C *ISAVE  - Integer array of length 45.
C           ISAVE(I), I=16, 17, 39, 40 are accessed
C           and changed by the routine.
C           ISAVE(16) - On entry to first call to MA42J/JD, ISAVE(16)
C           is set to -2 if the call directly follows calls to MA42A/AD
C           and is set to 4 if the call is immediately after INFO(1)=4
C           has been returned from MA42F/FD. On exit from intermediate
C           calls to MA42J/JD, ISAVE(16) is set to 1, and on
C           exit from the final call, ISAVE(16)=2
C           ISAVE(38) updated if there is an error
C *INFO   - Integer array of length 23. Only INFO(1) is
C           accessed and is used as an error flag.
C           Possible nonzero flags are -2, -3, -8, -13.
C
C Local Variables
C
C  LK     - Integer variable. Do loop variable. Loop is over
C           the NVAR variables in the elt/equ.
C  IELL   - Integer variable. Current elt/equ number. (ISAVE(17)).
C  ISTATC - Integer variable. Number of static condensations
C           for current elt/equ.
C  KFRNT  - Integer variable. Number of rows currently in front.
C           (ISAVE(39)).
C  L      - Integer variable. Do loop variable.
C  LFRNT  - Integer variable. Number of cols currently in front.
C           (ISAVE(40)).
C  LP     - Integer variable. Stream number for error messages.
C           (ICNTL(1))
C  MFR    - Integer variable. Used to index variables in the
C           elt/equ.
C
C     .. Scalar Arguments ..
      INTEGER NDF,NMAXE,NVAR
C     ..
C     .. Array Arguments ..
      INTEGER ICNTL(8),IFSIZE(5),INFO(23),ISAVE(45),IVAR(NVAR),LAST(NDF)
C     ..
C     .. Local Scalars ..
      INTEGER IELL,ISTATC,KFRNT,L,LFRNT,LK,LP,MFR,NELL
      INTEGER PIVBLK,KR,J,J1
C We have introduced PIVBLK as smallest pivot block size allowed.
C Note that HSL12 version is equivalent to setting PIVBLK = 1
C (This is default ... user must reset ISAVE(19) if some other value
C is required)
C KR is the number of fully summed variables. We have to save KR
C between calls. We will use ISAVE(41) for this.
C
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX
C     ..
C Immediate return if previous error
      IF (ISAVE(38).LT.0) GO TO 120
      LP = ICNTL(1)
      IELL = ISAVE(17)
      NELL = ISAVE(18)
      IELL = IELL + 1
      KFRNT = ISAVE(39)
      LFRNT = ISAVE(40)

      PIVBLK = ABS(ISAVE(19))

      IF (ISAVE(16).EQ.-2) THEN
C First call to MA42J/JD since MA42A/AD called.
C Initialise IFSIZE and INFO(1).
C IFSIZE(5) is initialised to 1 because of last entry in integer file
C which is an end-of-file flag.
         INFO(1) = 0
         DO 10 L = 1,4
            IFSIZE(L) = 0
   10    CONTINUE
         IFSIZE(5) = 1
C Change flag in ISAVE(16) and store NMAXE.
         ISAVE(16) = -1
         ISAVE(21) = NMAXE
C Jump to error return if NDF is changed from value output from
C MA42A/AD or if NMAXE is less than or equal to 0.
         IF (NDF.NE.ISAVE(30)) GO TO 70
         IF (NMAXE.LE.0) GO TO 50

         ISAVE(41) = 0

      ELSE IF (ISAVE(16).EQ.4) THEN
C First call to MA42J/JD after failure in MA42F/FD due to
C insufficient space. Initialise IFSIZE.
         IFSIZE(1) = MAX(INFO(8),KFRNT)
         IFSIZE(2) = MAX(INFO(9),LFRNT)
         IFSIZE(3) = INFO(4)
         IFSIZE(4) = INFO(6)
         IFSIZE(5) = INFO(7)
C Change flag in ISAVE(16)
         ISAVE(16) = 1

         ISAVE(41) = 0

      END IF
C Check for errors in input data
      IF (NVAR.LE.0) GO TO 30
C      IF (NMAXE.NE.ISAVE(21)) GO TO 50
      IF (ISAVE(21).EQ.1 .AND. NMAXE.GT.1) GO TO 50
      IF (ISAVE(21).GT.1 .AND. NMAXE.EQ.1) GO TO 50
C
C For equ entry, update max number of rows in front.
      IF (NMAXE.EQ.1) KFRNT = KFRNT + 1
C
C Assemble incoming elt/eqn and do symbolic eliminations.
C Loop over variables in incoming elt/eqn.
      KR = ISAVE(41)
      ISTATC = 0
      DO 20 LK = 1,NVAR
         MFR = IVAR(LK)
         IF (MFR.LT.1 .OR. MFR.GT.NDF) GO TO 40
         IF (LAST(MFR).GE.0) THEN
C Variable not yet in front
C Jump to error return if variable is already fully summed.
            IF (LAST(MFR).LT.IELL) GO TO 60
C Test for static condensation.
            IF (LAST(MFR).EQ.IELL) ISTATC = ISTATC + 1
C Perform assembly.
            LFRNT = LFRNT + 1
            IF (LAST(MFR).GT.IELL) LAST(MFR) = -LAST(MFR)
         ELSE
C Variable MFR already in front. Test to see if it can be eliminated.
            IF (-LAST(MFR).EQ.IELL) THEN
C Restore LAST(MFR).
               LAST(MFR) = -LAST(MFR)
               KR = KR + 1
            END IF
         END IF
   20 CONTINUE
C Update max front sizes.
      IF (NMAXE.GT.1) THEN
         IFSIZE(2) = MAX(IFSIZE(2),LFRNT)
         IFSIZE(1) = IFSIZE(2)
         KFRNT = LFRNT
      ELSE IF (NMAXE.EQ.1) THEN
         IFSIZE(1) = MAX(IFSIZE(1),KFRNT)
         IFSIZE(2) = MAX(IFSIZE(2),LFRNT)
      END IF
C Following lines added 7.7.2000
      IF (ISTATC.GT.1 .AND. NMAXE.EQ.1) THEN
         KR = KR + ISTATC
         ISTATC = 0
      END IF
C Update front sizes and file sizes after static condensations.
C Note that for equation entry, we can avoid writing to the
C frontal matrix if a static condensation variable is found.
      IF (ISTATC.GT.0) THEN
         IF (ISTATC.EQ.1 .AND. NMAXE.EQ.1) THEN
            IFSIZE(5) = IFSIZE(5) + 5 + 1 + NVAR
            IFSIZE(3) = IFSIZE(3) + NVAR
            IFSIZE(4) = IFSIZE(4) + 1
            LFRNT = LFRNT - 1
            KFRNT = KFRNT - 1
         ELSE
            IFSIZE(5) = IFSIZE(5) + 5 + LFRNT + KFRNT
            LFRNT = LFRNT - ISTATC
            KFRNT = KFRNT - ISTATC
            IFSIZE(3) = IFSIZE(3) + ISTATC*LFRNT +
     +                  (ISTATC* (ISTATC+1))/2
            IFSIZE(4) = IFSIZE(4) + ISTATC*KFRNT +
     +                  (ISTATC* (ISTATC+1))/2
         END IF
      END IF
C Added this 21.7.00 ... could happen if matrix singular.
      IF (KFRNT.EQ.0 .OR. LFRNT.EQ.0) GO TO 24
C Test to see if there are enough fully summed variables
C to eliminate. Return for next element if there are fewer than
C PIVBLK fully summed variables.
C Update front sizes and file sizes after eliminations.
      IF (KR.LT.PIVBLK .AND. IELL.LT.NELL) GO TO 25
      IF (KR.GT.0) THEN
         IFSIZE(5) = IFSIZE(5) + 5 + LFRNT + KFRNT
         J1 = KR
         DO 23 J = 1,J1
            KFRNT = KFRNT - 1
            LFRNT = LFRNT - 1
            KR = KR - 1
            IF (KFRNT.EQ.0 .OR. LFRNT.EQ.0) THEN
              IFSIZE(3) = IFSIZE(3) + J*LFRNT + (J* (J+1))/2
              IFSIZE(4) = IFSIZE(4) + J*KFRNT + (J* (J+1))/2
              GO TO 24
            END IF
   23    CONTINUE
         IFSIZE(3) = IFSIZE(3) + J1*LFRNT + (J1* (J1+1))/2
         IFSIZE(4) = IFSIZE(4) + J1*KFRNT + (J1* (J1+1))/2
      END IF
  24  CONTINUE
C Update file sizes after eliminations.
C Reset ISAVE(16), ISAVE(17), ISAVE(39), ISAVE(40) before
C final return
      IF (IELL.EQ.NELL) GO TO 100
  25  ISAVE(17) = IELL
      ISAVE(39) = KFRNT
      ISAVE(40) = LFRNT
      ISAVE(41) = KR
      GO TO 110
C
C **** Error returns ****
   30 INFO(1) = -2
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1),IELL
         WRITE (LP,FMT=9010) NVAR
      END IF
      GO TO 80
   40 INFO(1) = -3
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1),IELL
         WRITE (LP,FMT=9020) LK,MFR
      END IF
      GO TO 80
   50 INFO(1) = -8
      IF (LP.GT.0) THEN
         IF (NMAXE.LE.0) THEN
            WRITE (LP,FMT=9070) INFO(1)
            WRITE (LP,FMT=9060) NMAXE
         ELSE
            WRITE (LP,FMT=9000) INFO(1),IELL
            WRITE (LP,FMT=9040) ISAVE(21),NMAXE
         END IF
      END IF
      GO TO 80
   60 INFO(1) = -13
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1),IELL
         WRITE (LP,FMT=9050) MFR
      END IF
      GO TO 80
   70 INFO(1) = -15
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9070) INFO(1)
         WRITE (LP,FMT=9030) ISAVE(30),NDF
      END IF
C Restore LAST after error return.
   80 DO 90 L = 1,NDF
         IF (LAST(L).LT.0) LAST(L) = -LAST(L)
   90 CONTINUE
  100 ISAVE(17) = 0
      ISAVE(39) = 0
      ISAVE(40) = 0
      ISAVE(16) = 2
C
  110 IF (INFO(1).LT.0) ISAVE(38) = INFO(1)
  120 RETURN
 9000 FORMAT (' ***** Error return from MA42J/JD *****  INFO(1) = ',I3,
     +       /' ***** Error discovered after input of elt/eqn *****',I8)
 9010 FORMAT (7X,'Number of variables in elt/eqn is',I8)
 9020 FORMAT (7X,'Variable',I8,' in elt/eqn has value',I8)
 9030 FORMAT (' NDF has been changed from ',I8,' to ',I8)
 9040 FORMAT (7X,'NMAXE has been changed from',I8,' to ',I8)
 9050 FORMAT (7X,'Variable',I8,' is already fully summed')
 9060 FORMAT (7X,'NMAXE is equal to ',I8)
 9070 FORMAT (' ***** Error return from MA42J/JD *****  INFO(1) = ',I3)
      END
C*******************************************************************
      SUBROUTINE MA42PD(ISTRM,LENBUF,LENFLE,ICNTL,ISAVE,INFO)
C
C *** Reference MA42 suite ***
C *** Copyright Rutherford Appleton Laboratory  November 1992 ***
C *** Any problems contact Iain S. Duff  or Jennifer A. Scott
C     at Atlas Centre, Rutherford Appleton Laboratory ***
C *** Although every effort has been made to ensure robustness and
C     reliability of the subroutines in this MA42 suite,  we
C     disclaim any liability arising through the use or misuse of
C     any of the subroutines in the MA42 suite ***
C
C This subroutine is optional. It initializes the direct access data
C sets.
C
C     Argument list. * indicates the argument is changed by the routine.
C
C   ISTRM  - Integer array of length 3.
C            ISTRM(1) - Stream number for U direct access data set
C            ISTRM(2) - Stream number for L direct access data set.
C            ISTRM(3) - Stream number for integer direct access data
C            set.
C   LENBUF - Integer array of length 3.
C            LENBUF(1) - Buffer length for U buffer.
C            LENBUF(2) - Buffer length for L buffer.
C            LENBUF(3) - Buffer length for integer buffer.
C   LENFLE - Integer array of length 3.
C            LENFLE(1) - Length of direct access data set for U.
C            LENFLE(2) - Length of direct access data set for L.
C            LENFLE(3) - Length of direct access data set for
C            integer information.
C   ICNTL  - Integer array of length 8.
C            ICNTL(1) must be set to stream number for error messages.
C            ICNTL(2) must be set to stream number for warnings.
C            ICNTL(3) - must be set to the number
C            of bytes for a real record (in general
C            4 for single precision and 8 for double precision)
C            ICNTL(4) - must be set to the number
C            of bytes for an integer record (in general 4)
C            ICNTL(I), I=5,...,8 not accessed.
C  *ISAVE  - Integer array of length 45.
C            On successful exit the following entries have been changed.
C            ISAVE(1) - Stream number for U direct access data set.
C            ISAVE(2) - Stream number for L direct access data set.
C            ISAVE(3) - Stream number for integer direct access data
C            set.
C            ISAVE(I) is equivalent to the local variable IFILE(I)
C            (I = 1,2,3).
C            ISAVE(7) - Length of U buffer.
C            ISAVE(8) - Length of L buffer.
C            ISAVE(9) - Length of integer buffer.
C            ISAVE(6+I) is equivalent to the local variable ISIZE(I)
C            (I = 1,2,3).
C            ISAVE(13) - Number of records in U direct access data set.
C            ISAVE(14) - Number of records in L direct access data set.
C            ISAVE(15) - Number of records in integer direct access data
C                        set.
C            ISAVE(12+I) is equivalent to the local variable NUMBLK(I)
C            (I = 1,2,3).
C            ISAVE(31) is set to +1 to indicate MA42P/PD has been
C            called.
C            ISAVE(38) updated if there is an error
C  *INFO   - Integer array of length 23.
C            On exit a non-zero value for INFO(1) indicates
C            a fatal error. Possible nonzero values -20,-21,-22,-23,-24.
C            Remainder of array not accessed by routine.
C
C   Local variables
C
C   IFILE  - Integer array of length 3. IFILE(I) is equivalent to
C            ISAVE(I) (I = 1,2,3).
C   ISIZE  - Integer array of length 3. ISIZE(I) is equivalent to
C            ISAVE(6+I) (I = 1,2,3).
C   NUMBLK - Integer array of length 3. NUMBLK(I) is equivalent to
C            ISAVE(12+I) (I = 1,2,3).
C   I      - Integer variable. Do loop variable.
C   ICON   - Integer variable. Used to hold number of bytes in direct
C            access data set record.
C   IOS    - Integer variable. Holds IOSTAT parameter for direct access
C            data sets.
C   LP     - Integer variable. Holds copy of ICNTL(1) (stream number
C            for error messages).
C
C     .. Array Arguments ..
      INTEGER ICNTL(8),INFO(23),ISAVE(45),ISTRM(3),LENBUF(3),LENFLE(3)
C     ..
C     .. Local Arrays ..
      INTEGER IFILE(3),ISIZE(3),NUMBLK(3)
c      CHARACTER*50 FILNAM(3)
C     ..
C     .. Local Scalars ..
      INTEGER I,ICON,IOS,LP,MP
C     ..
C Immediate return if there was a previous error
      IF (ISAVE(38).LT.0) GO TO 110
      INFO(1) = 0
      ISAVE(31) = 1
      LP = ICNTL(1)
      MP = ICNTL(2)
C
      DO 10 I = 1,3
C Check ISTRM(I) is in range and set IFILE(I)
         IF (ISTRM(I).EQ.LP .OR. ISTRM(I).EQ.MP .OR. ISTRM(I).LT.0 .OR.
     +       ISTRM(I).GT.99 .OR. ISTRM(I).EQ.6) GO TO 90
         IF (I.NE.2 .AND. ISTRM(I).EQ.0) GO TO 90
         IFILE(I) = ISTRM(I)
         IF (IFILE(I).GT.0) THEN
C (IFILE(1) and IFILE(3) are positive but IFILE(2) may be 0).
            IF (LENBUF(I).LE.0) GO TO 40
            NUMBLK(I) = LENFLE(I)/LENBUF(I)
            IF (NUMBLK(I).EQ.0) GO TO 50
         ELSE
C (Set LENBUF(2) = 0 if ISTRM(2) was equal to 0).
            LENBUF(I) = 0
            NUMBLK(I) = 0
         END IF
         ISIZE(I) = LENBUF(I)
   10 CONTINUE
C
C Check for duplicate stream numbers
      IF (IFILE(1).EQ.IFILE(2)) GO TO 60
      IF (IFILE(1).EQ.IFILE(3)) GO TO 60
      IF (IFILE(2).EQ.IFILE(3)) GO TO 60
C Open direct access data sets.
      IF (ICNTL(3).LT.0 .OR. ICNTL(4).LT.0) GO TO 80
      DO 20 I = 1,3
c These names are for running problems on parasol
c         if (i.eq.1) filnam(1) = '/tmp/sctu.2'
c         if (i.eq.2) filnam(2) = '/tmp/sctl.2'
c         if (i.eq.3) filnam(3) = '/tmp/scti.2'
c These names are for running problems on rs/6000
c         if (i.eq.1) filnam(1) = '/rutherford/num-ftp/pub/open/ufact'
c         if (i.eq.2) filnam(2) = '/rutherford/num-ftp/pub/open/lfact'
c         if (i.eq.3) filnam(3) = '/rutherford/num-ftp/pub/open/integ'
         ICON = ICNTL(3)
         IF (I.EQ.3) ICON = ICNTL(4)
         IF (IFILE(I).GT.0) THEN
            CLOSE (IFILE(I),STATUS='DELETE')
            OPEN (IFILE(I),ERR=70,ACCESS='DIRECT',
     +           RECL=LENBUF(I)*ICON,IOSTAT=IOS)
c            OPEN (IFILE(I),ERR=70,ACCESS='DIRECT',FILE=FILNAM(I),
c     +           RECL=LENBUF(I)*ICON,IOSTAT=IOS)
c magellan wants status=unknown
c            OPEN (IFILE(I),ERR=70,ACCESS='DIRECT',STATUS='unknown',
c     +           RECL=LENBUF(I)*ICON,IOSTAT=IOS)
         END IF
   20 CONTINUE
C Save IFILE, ISIZE, AND NUMBLK
      DO 30 I = 1,3
         ISAVE(I) = IFILE(I)
         ISAVE(3+I) = 1
         ISAVE(6+I) = ISIZE(I)
         ISAVE(9+I) = 0
         ISAVE(12+I) = NUMBLK(I)
   30 CONTINUE
      ISAVE(16) = 2
      ISAVE(17) = 0
      ISAVE(39) = 0
      ISAVE(40) = 0
      GO TO 100
C
C **** Error returns ****
   40 INFO(1) = -19
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1)
         WRITE (LP,FMT=9010) I,LENBUF(I)
      END IF
      GO TO 100
   50 INFO(1) = -20
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1)
         WRITE (LP,FMT=9020) I,LENFLE(I),I,LENBUF(I)
      END IF
      GO TO 100
   60 INFO(1) = -21
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1)
         IF (IFILE(1).EQ.IFILE(2)) WRITE (LP,FMT=9030) ISTRM(1)
         IF (IFILE(1).EQ.IFILE(3)) WRITE (LP,FMT=9040) ISTRM(1)
         IF (IFILE(2).EQ.IFILE(3)) WRITE (LP,FMT=9050) ISTRM(2)
      END IF
      GO TO 100
   70 INFO(1) = -22
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1)
         WRITE (LP,FMT=9060) IFILE(I),IOS
      END IF
      GO TO 100
   80 INFO(1) = -23
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1)
         WRITE (LP,FMT=9080) ICNTL(3),ICNTL(4)
      END IF
      GO TO 100
   90 INFO(1) = -24
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1)
         WRITE (LP,FMT=9070) I,ISTRM(I)
      END IF
C
  100 IF (INFO(1).LT.0) ISAVE(38) = INFO(1)
  110 RETURN
 9000 FORMAT (' ***** Error return from MA42P/PD ***** INFO(1) = ',I3)
 9010 FORMAT (7X,'Non-positive buffer size (LENBUF(',I1,')) of',I8)
 9020 FORMAT (7X,'File length (LENFLE(',I1,')) ',I8,/7X,'is less than ',
     +       'buffer size (LENBUF(',I1,')) of ',I8)
 9030 FORMAT (7X,'Stream numbers for L and U direct access data sets ',
     +       /7X,'both equal to ',I3)
 9040 FORMAT (7X,'Stream numbers for U and integer direct access data',
     +       ' sets',/7X,'both equal to ',I3)
 9050 FORMAT (7X,'Stream numbers for L and integer direct access data',
     +       ' sets',/7X,'both equal to ',I3)
 9060 FORMAT (7X,'Failure in direct access open',/7X,'file on stream',
     +       I3,5X,'IOSTAT = ',I10)
 9070 FORMAT (7X,'Illegal stream number. ISTRM(',I1,') is set to ',I3)
 9080 FORMAT (7X,'ICNTL(3) and ICNTL(4) are set to ',I3,',',I3)
      END
C**********************************************************************
      SUBROUTINE MA42BD(NVAR,IVAR,NDF,LAST,NMAXE,AVAR,NRHS,RHS,LRHS,LX,
     +                  X,NFRONT,LENBUF,LW,W,LIW,IW,ICNTL,CNTL,ISAVE,
     +                  INFO,RINFO)
C
C *** Reference MA42 suite ***
C *** Copyright Rutherford Appleton Laboratory  November 1992 ***
C *** Any problems contact Iain S. Duff  or Jennifer A. Scott
C     at Atlas Centre, Rutherford Appleton Laboratory ***
C *** Although every effort has been made to ensure robustness and
C     reliability of the subroutines in this MA42 suite,  we
C     disclaim any liability arising through the use or misuse of
C     any of the subroutines in the MA42 suite ***
C
C This subroutine is the user-called driver for the main elimination
C routine of the frontal method (MA42F/FD).  After the last elt/eqn
C has been input, the back-substitution routine (MA42D/DD) is
C automatically called, thus solving for any right hand sides input
C with the matrix.  Most of the checking of the users input data
C is also performed by this subroutine.
C
C     Argument list. * indicates the argument is changed by the routine.
C
C  NVAR   - Integer variable. Number of variables in elt/eqn being
C           input.
C *IVAR   - Integer array of length NVAR. Contains indices of variables
C           in elt/eqn being input. Used locally as workspace.
C  NDF    - Integer variable. Largest integer used to index a variable.
C *LAST   - Integer array of length NDF.  On input LAST(I) must be the
C           elt/eqn in which variable I appears for the last
C           time.  Between calls to MA42B/BD, LAST is used as a
C           work array but on final exit or on detection of a
C           fatal error, LAST is reset and returned unchanged to the
C           user. If variable J (J.LE.NDF) is not used to index
C           a variable, LAST(J) is not accessed.
C  NMAXE  - Integer variable. Leading dimension of element or
C           equation arrays. Must be set to 1 for equation entry.
C  AVAR   - Real (DP) array of dimensions NMAXE by NVAR.
C           Must contain contributions to matrix from elt/eqn
C           being input.
C  NRHS   - Integer variable.  Number of right hand sides.
C  RHS    - Real (DP) array of dimensions NMAXE by LRHS.
C           Must contain right hand sides for elt/eqn being input.
C  LRHS   - Integer variable. Second dimension of RHS array.
C           Must be at least as great as the maximum of 1 and NRHS.
C  LX     - Integer variable.  Leading dimension of array X. If NRHS is
C           at least 1, LX must be at least as large as NDF as output
C           from MA42A/AD.
C *X      - Real (DP) array of dimensions LX by LRHS.
C           On successful exit X(I,J) will hold the Ith component of
C           the solution to system J.
C *NFRONT - Integer array of length 2.  On exit, NFRONT(I),I=1,2 holds
C           respectively the number of rows and columns in the in-core
C           frontal matrix used to hold front during the decomposition.
C           If NMAXE is greater than 1 (element entry), NFRONT(2)
C           is set equal to NFRONT(1).  If there is insufficient
C           space for the decomposition, NFRONT may, on exit, give a
C           lower bound to the space required for a successful
C           decomposition (see call to MA42J/JD).
C  LENBUF - Integer array of length 3. If MA42P/PD was called, LENBUF
C           must be unchanged since that call. Otherwise, LENBUF
C           must hold in-core buffer sizes (no direct access files
C           being used).
C  LW     - Integer variable. Length of array W.
C           must be at least as great as 1+LENBUF(1)+LENBUF(2)+
C           NFRONT(1)*NFRONT(2)+MAX(LRHS*NFRONT(1),NRHS*NFRONT(2).
C *W      - Real (DP) array of length LW. Used as workspace
C           and, on exit from all calls, will hold the solutions.
C  LIW    - Integer variable. Length of array IW.  It must have a
C           value at least as great as LENBUF(3)+2*NFRONT(1)+
C           4*NFRONT(2).
C *IW     - Integer array of length LIW.  Used as workspace.
C  ICNTL  - Integer array of length 8. Holds control parameters.
C           ICNTL(1) must be set by user to stream number for error
C           messages.
C           ICNTL(2) - Stream number for warnings. (default 6)
C           ICNTL(3) - must be set by the user to the number
C           of bytes for a real record (in general
C           4 for single precision and 8 for double precision)
C           ICNTL(4) - must be set by the user to the number
C           of bytes for an integer record (in general 4)
C           If ICNTL(5) is positive then, when the number
C           of potential pivot columns is greater than or equal to
C           ICNTL(6), an elimination will be performed even if the best
C           pivot candidate does not satisfy the threshold criterion
C           determined by CNTL(2). Default is ICNTL(5)=0
C           ICNTL(6) - default 0. If ICNTL(6) and ICNTL(5) are positive
C           then only ICNTL(6) potential pivot
C           columns will be searched for a pivot.
C           The best pivot candidate (largest relative
C           to other nonzeros in its column) from these ICNTL(6)
C           columns will then be used as a pivot.
C           ICNTL(7) - default 0. If
C           ICNTL(7) =  0 Static condensation + off-diag. pivoting
C           ICNTL(7) =  1001 Static condensation + diag. pivoting
C           ICNTL(7) = -1001 No static condensation + diag. pivoting
C           ALL other values of ICNTL(7):
C           No static condensation + off-diag. pivoting.
C           ICNTL(8) - default 0. If ICNTL(8) is 0
C           immediate exit if singularity in matrix detected.
C  CNTL   _ Real (DP) array of length 2. Holds control
C           parameters.
C           CNTL(1) - default 0.0. The matrix will be declared singular
C           if the element of largest absolute value in any column
C           during decomposition is less than or equal CNTL(1).
C           CNTL(2) - default 0.1.
C           An element of the frontal matrix will only be considered
C           suitable for use as a pivot if it has absolute value
C           at least as large as CNTL(2) times the  element of
C           largest absolute value in its column.
C *ISAVE  - Integer array of length 45. Holds parameters which
C           must be preserved between calls to MA42B/BD and between
C           calls to other MA42 routines.
C           ISAVE(1),...,ISAVE(15) were set by MA42I/ID and MA42P/PD.
C           If MA42P/PD was called, only ISAVE(3+I) and ISAVE(9+I),
C           I = 1,2,3 are changed by MA42B/BD.
C           Otherwise, ISAVE(6+I) is set to hold LENBUF(I), I=1,2,3
C           on the first entry to MA42B/BD.
C           ISAVE(16) is only used if MA42J/JD is called because
C           if insufficient space allocated to frontal matrix.
C           ISAVE(16) is used as a flag to indicate when MA42J/JD
C           is called for the first time.
C           ISAVE(17) is set to 0 by MA42I/ID.
C           On exit from intermediate calls to MA42B/BD, ISAVE(17) holds
C           the number of elements input so far. On return from the
C           final call, ISAVE(17) = 0 (this allows MA42B/BD
C           to be rerun without rerunning MA42A/AD).
C           ISAVE(18) set by MA42A/AD to number of elements in problem
C           and is unchanged by MA42B/BD.
C           ISAVE(19) used to hold min. pivot block size.
C           ISAVE(19) and ISAVE(20) are initialised to 1 by MA42I/ID,
C           but user can reset ISAVE(19) to hold min. pivot block size
C           before calls to MA42J/JD.
C           If ISAVE(20) reset to value greater than 2, then
C           zeros in front exploited.
C           A negative value of ISAVE(19) or ISAVE(20) indicates
C           the occurrence of the
C           non-fatal error 5 or 6 respectively.
C           ISAVE(19) or ISAVE(20) restored to original values before
C           final exit.
C           ISAVE(21), ISAVE(22), ISAVE(23) and ISAVE(24)
C           hold the values of variables NMAXE, NRHS,
C           NFRONT(1) and NFRONT(2), resp. as input in the first call.
C           ISAVE(25)-ISAVE(29) are used to divide real workspace.
C           ISAVE(30)-set by MA42A/AD to largest integer used to index
C           a variable.
C           ISAVE(31)-flag to indicate if MA42P/PD has been called.
C           ISAVE(32)-ISAVE(37) are used to divide integer workspace.
C           ISAVE(38)-is set to INFO(1) if an error (or a warning which
C           will lead to an error) is encountered, and to
C           abs(ICNTL(7)) otherwise
C           ISAVE(39) is set to 0 by MA42I/ID. On exit from
C           intermediate calls to MA42B/BD, ISAVE(39) holds current row
C           front size KFRNT. On return from final call, ISAVE(39)=0.
C           ISAVE(40) is set to 0 by MA42I/ID. On exit from
C           intermediate calls to MA42B/BD, ISAVE(40) holds current
C           column front size LFRNT. On return from final call,
C           ISAVE(40)=0.
C           ISAVE(41)-ISAVE(45) used to preserve values between calls
C           to MA42F/FD.
C *INFO   - Integer array of length 23.
C           INFO(1) need not be set by the user.
C           INFO(1) will have value 0 on successful exit. Negative
C           values of INFO(1) indicate a fatal error, while
C           values greater than 0 are non-terminal errors and are
C           associated with a warning or error message.
C           Possible negative values on exit are -2, -3, -5, -6,
C           -7, -8, -9, -10, -11, -12, -13, -14, -15, -16, -17,
C           -18, -25, -26.
C           Possible positive values on exit are 1, 2, 3, 4, 5, 6.
C INFO(2)  - Sign of determinant of matrix
C INFO(3)  - Number of variables.
C INFO(4)  - Number of non-zeros in U.
C INFO(5)  - Total real storage for U (plus right-hand sides).
C INFO(6)  - Number of non-zeros in L. (=Total real storage for L)
C INFO(7)  - Total integer storage.
C INFO(8)  - Maximum front-width (rows).
C INFO(9)  - Maximum front-width (columns).
C INFO(10) - Number of buffers for factors of U.
C INFO(11) - Number of buffers for factors of L.
C INFO(12 )- Number of buffer for integers.
C INFO(13) - Number of columns examined during pivot searches.
C INFO(14) - Number of non-zeros tested for stability as pivots.
C INFO(15) - Number of non-zeros accessed during pivot selection
C            process.
C INFO(16) - Number of pivots chosen which did not satisfy threshold
C            criterion based on the value of CNTL(2).
C INFO(17) - Number of static condensations performed.
C INFO(18) - Potential number of static condensations (may be different
C            from INFO(17) because of stability considerations).
C INFO(19) - Maximum number of buffers required to hold block of pivot
C            rows.
C INFO(20) - Maximum number of buffers required to hold block of pivot
C            cols.
C INFO(21) - Max number of buffers required to hold block of integers.
C INFO(22) - Size of maximum pivot block.
C INFO(23) - Deficiency of matrix (singular case).
C
C *RINFO  _ Real (DP) array of length 2. On final exit,
C           RINFO(1) holds natural log of determinant of the matrix
C           RINFO(2) holds number of operations in innermost loop
C
C   Local variables
C
C   DIMBUF - Integer variable. Used in backsubstitution to
C            hold INFO(19)*JBUFU.
C   DIMIBF - Integer variable. Used in backsubstitution to
C            hold INFO(21)*ISIZE(3).
C   FA     - Integer variable. Used in dividing real workspace.
C            The frontal matrix is held in W(FA+I-1),
C            I=1,2,...,NFRONT(1)*NFRONT(2). (ISAVE(28)).
C   FRHS   - Integer variable. Used in dividing real workspace.
C            The right hand sides corresponding to the current
C            frontal method are held in W(FRHS+I-1),
C            I=1,2,...,NFRONT(1)*LRHS. (ISAVE(29)).
C   I      - Integer variable. Do loop variable.
C   IELL   - Integer variable. IELL is elt/eqn currently being input
C            (ISAVE(17)).
C   JBUFL  - Integer variable. Holds max(1,ISIZE(2)). Used in dividing
C            real workspace.
C   JBUFR  - Integer variable. Holds 1+SIZE(2). Used in dividing
C            real workspace.
C   JBUFU  - Integer variable. Holds ISIZE(1) (size of U-buffer).
C            Used in dividing real workspace.
C   K      - Integer variable. Do loop variable.
C   KDEST  - Integer variable. Used in dividing integer workspace.
C            IW(KDEST+I-1), I=1,2,...,NFRONT(1) are used to hold
C            information on row indices in frontal matrix.
C            (ISAVE(37)).
C   KFRNT  - Integer variable. Holds current row front size
C            (ISAVE(39)).
C   KHED   - Integer variable. Used in dividing integer workspace.
C            The indices corresponding to rows in the front are
C            held in IW(KHED+I-1), I=1,2,...,NFRONT(1).
C            (ISAVE(33)).
C   KPIV   - Integer variable. Used in dividing integer workspace.
C            The indices of the columns of the
C            frontal matrix which are fully summed are held in
C            IW(KPIV+I-1), I=1,2,...,NFRONT(2). (ISAVE(34)).
C   KPVLNK - Integer variable. Used in dividing integer workspace.
C            IW(KPVLNK+I-1), I=1,2,...,NFRONT(2) is used as the
C            inverse array for IW(KPIV+I-1), I=1,2,...,NFRONT(2).
C            (ISAVE(35)).
C   L      - Integer variable. Do loop variable.
C   LDEST  - Integer variable. Used in dividing integer workspace.
C            IW(LDEST+I-1), I=1,2,...,NFRONT(2) are used to hold
C            information on column indices in frontal matrix.
C            (ISAVE(36)).
C   LFL    - Integer variable. Used to hold required file size
C            if INFO(1) = 5 or 6 is returned.
C   LFRNT  - Integer variable. Holds current column front size
C            (ISAVE(40)).
C   LHED   - Integer variable. Used in dividing integer workspace.
C            The indices corresponding to columns in the front are
C            held in IW(LHED+I-1), I=1,2,...,NFRONT(2).(ISAVE(32)).
C   LLIW   - Integer variable. Used to hold amount of integer
C            workspace required by routine (used to check if LIW is
C            sufficient).
C   LLW    - Integer variable. Used to hold amount of real
C            workspace required by routine (used to check if LW is
C            sufficient).
C   LP     - Integer variable. Stream number for error messages.
C   NELL   - Integer variable. Holds  the total number of
C            elts/eqns input to MA42A/AD. (ISAVE(18)).
C   NFR    - Integer variable. Used on first call to MA42J/JD to
C            count the number of non-fully summed variables in the
C            front immediately after failure of MA42F/FD.
C   IFILE  - Integer array of length 3. IFILE(I) is equivalent to
C            ISAVE(I) (I = 1,2,3).
C   IREC   - Integer array of length 3. IREC(I) is equivalent to
C            ISAVE(3+I) (I = 1,2,3).
C   ISIZE  - Integer array of length 3. ISIZE(I) is equivalent to
C            ISAVE(6+I) (I = 1,2,3).
C   MKEY   - Integer array of length 3. MKEY(I) is equivalent to
C            ISAVE(9+I) (I = 1,2,3).
C   NUMBLK - Integer array of length 3. NUMBLK(I) is equivalent to
C            ISAVE(12+I) (I = 1,2,3).
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER LIW,LRHS,LW,LX,NDF,NMAXE,NRHS,NVAR
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION AVAR(NMAXE,NVAR),CNTL(2),RHS(NMAXE,LRHS),
     +                 RINFO(2),W(LW),X(LX,LRHS)
      INTEGER ICNTL(8),INFO(23),ISAVE(45),IVAR(NVAR),IW(LIW),LAST(NDF),
     +        LENBUF(3),NFRONT(2)
C     ..
C     .. Local Scalars ..
      INTEGER DIMBUF,DIMIBF,FA,FRHS,I,IELL,JBUFL,JBUFR,JBUFU,K,KDEST,
     +        KFRNT,KHED,KPIV,KPVLNK,L,LDEST,LFL,LFRNT,LHED,LLIW,LLW,LP,
     +        MFR,NELL,NFR,PIVBLK,BZERO
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION Y(1,1)
      INTEGER IFILE(3),IREC(3),ISIZE(3),MKEY(3),NUMBLK(3)
C     ..
C     .. External Subroutines ..
      EXTERNAL MA42DD,MA42FD,MA42JD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX
C     ..
C IELL is elt/eqn currently being input, NELL is the total number of
C elts/eqns input to MA42A/AD.
      IELL = ISAVE(17)
      NELL = ISAVE(18)
      IELL = IELL + 1
C Return immediately if previous entry has caused a fatal error.
      IF (ISAVE(38).LT.0) GO TO 280
      DO 10 I = 1,3
C Reset ICNTL(1) and/or ICNTL(2) if the user has chosen values
C equal to a stream number.
         IF (ISAVE(31).EQ.1) THEN
            IFILE(I) = ISAVE(I)
            IF (ICNTL(1).EQ.IFILE(I) .AND. IFILE(I).NE.0) ICNTL(1) = 6
            IF (ICNTL(2).EQ.IFILE(I) .AND. IFILE(I).NE.0) ICNTL(2) = 6
         END IF
   10 CONTINUE
      LP = ICNTL(1)
C Initialization phase.
      IF (IELL.EQ.1) THEN
C Reset ICNTL(7) if equation entry and diagonal pivoting requested
         IF (NMAXE.EQ.1 .AND. ABS(ICNTL(7)).EQ.1001) THEN
            IF (ICNTL(7).EQ.1001) ICNTL(7) = 0
            IF (ICNTL(7).EQ.-1001) ICNTL(7) = 1
         END IF
         ISAVE(38) = ABS(ICNTL(7))
         INFO(1) = 0
         INFO(2) = 1
         DO 20 I = 3,23
            INFO(I) = 0
   20    CONTINUE
         DO 30 I = 1,2
            RINFO(I) = ZERO
   30    CONTINUE
         DO 40 I = 1,3
            ISAVE(3+I) = 1
            ISAVE(9+I) = 0
   40    CONTINUE
C If MA42P/PD has not been called, set ISIZE(I) (and ISAVE(6+I))
         IF (ISAVE(31).EQ.-1) THEN
            DO 50 I = 1,3
               ISIZE(I) = LENBUF(I)
               ISAVE(6+I) = ISIZE(I)
               ISAVE(I) = 0
               ISAVE(12+I) = 0
   50       CONTINUE
C Jump to error return if LENBUF(1) or LENBUF(3) is non-positive.
            IF (LENBUF(1).LE.0 .OR. LENBUF(3).LE.0) GO TO 250
C Jump to error return if LENBUF(2) is negative.
            IF (LENBUF(2).LT.0) GO TO 250
         ELSE
C Check LENBUF(I) has not been changed since call to MA42P/PD.
            DO 60 I = 1,3
               IF (LENBUF(I).NE.ISAVE(6+I)) GO TO 250
               ISIZE(I) = ISAVE(6+I)
   60       CONTINUE
         END IF
         IF (NRHS.GE.1 .AND. LX.LT.NDF) GO TO 140
         IF (LRHS.LT.NRHS .OR. LRHS.LT.1) GO TO 190
         IF (NDF.NE.ISAVE(30)) GO TO 230
         IF (NRHS.LT.0) GO TO 240
         IF (NMAXE.LE.0) GO TO 170
C A negative value of ISAVE(19) or ISAVE(20) indicates the
C occurrence of the non-fatal error 5 or 6 respectively.
C Ensure both parameters are positive.
         ISAVE(19) = ABS(ISAVE(19))
         ISAVE(20) = ABS(ISAVE(20))
C Set NFRONT(2) if entry is by elements.
         IF (NMAXE.GT.1) NFRONT(2) = NFRONT(1)
C ISAVE(21), ISAVE(22), ISAVE(23) and ISAVE(24)
C hold the values of variables NMAXE, NRHS,
C NFRONT(1) and NFRONT(2), resp. as input in the first call.
C This is to allow these variables to be checked for consistency
C in future calls.
         ISAVE(21) = NMAXE
         ISAVE(22) = NRHS
         ISAVE(23) = NFRONT(1)
         ISAVE(24) = NFRONT(2)
C Divide up buffer space.
         JBUFU = ISIZE(1)
         JBUFL = MAX(1,ISIZE(2))
         JBUFR = 1 + ISIZE(2)
C Divide up the rest of the real workspace.
         FA = JBUFR + JBUFU
         FRHS = FA + NFRONT(1)*NFRONT(2)
C Jump to error return if insufficient space allocated to real
C workspace.
C MA42F/FD requires workspace FRHS + NFRONT(1)*LRHS
C MA42D/DD requires
C FRHS + NFRONT(2)*NRHS
C Hence sufficient to have FRHS + max(NFRONT(1)*LRHS, NFRONT(2)*NRHS)
         LLW = FRHS + MAX(NFRONT(1)*LRHS,NFRONT(2)*NRHS)
         IF (LW.LT.LLW) GO TO 150
         ISAVE(25) = JBUFU
         ISAVE(26) = JBUFL
         ISAVE(27) = JBUFR
         ISAVE(28) = FA
         ISAVE(29) = FRHS
C Divide up integer workspace.
         LHED = 1 + ISIZE(3)
         KHED = LHED + NFRONT(2)
         KPIV = KHED + NFRONT(1)
         KPVLNK = KPIV + NFRONT(2)
         LDEST = KPVLNK + NFRONT(2)
         KDEST = LDEST + NFRONT(2)
C Jump to error return if insufficient space allocated to integer
C workspace.
         LLIW = KDEST + NFRONT(1) - 1
         IF (LIW.LT.LLIW) GO TO 160
         ISAVE(32) = LHED
         ISAVE(33) = KHED
         ISAVE(34) = KPIV
         ISAVE(35) = KPVLNK
         ISAVE(36) = LDEST
         ISAVE(37) = KDEST
      END IF
C End of initialisation phase.
C
      DO 70 I = 1,3
         IFILE(I) = ISAVE(I)
         IREC(I) = ISAVE(3+I)
         ISIZE(I) = ISAVE(6+I)
         MKEY(I) = ISAVE(9+I)
         NUMBLK(I) = ISAVE(12+I)
   70 CONTINUE
C Jump if we have already exceeded front size.  Symbolic decomposition
C only is now being performed.
      IF (INFO(1).EQ.4) GO TO 120
C Test validity of input parameters.
      IF (NVAR.LE.0) GO TO 130
C      IF (NMAXE.NE.ISAVE(21)) GO TO 170
      IF (ISAVE(21).EQ.1 .AND. NMAXE.GT.1) GO TO 170
      IF (ISAVE(21).GT.1 .AND. NMAXE.EQ.1) GO TO 170
      IF (NMAXE.GT.1 .AND. NVAR.GT.NMAXE) GO TO 130
      IF (NRHS.NE.ISAVE(22)) GO TO 180
      IF (NFRONT(1).NE.ISAVE(23) .OR. NFRONT(1).LE.0) GO TO 200
      IF (NFRONT(2).NE.ISAVE(24) .OR. NFRONT(2).LE.0) GO TO 200
C
C Perform assembly and elimination operations.
C
      IF (IELL.GT.1) THEN
         JBUFU = ISAVE(25)
         JBUFL = ISAVE(26)
         JBUFR = ISAVE(27)
         FA = ISAVE(28)
         FRHS = ISAVE(29)
         LHED = ISAVE(32)
         KHED = ISAVE(33)
         KPIV = ISAVE(34)
         KPVLNK = ISAVE(35)
         LDEST = ISAVE(36)
         KDEST = ISAVE(37)
      END IF
      KFRNT = ISAVE(39)
      LFRNT = ISAVE(40)
      PIVBLK = ABS(ISAVE(19))
      BZERO = ABS(ISAVE(20))
      CALL MA42FD(AVAR,RHS,LRHS,NRHS,IVAR,LAST,NDF,NMAXE,NVAR,NFRONT(1),
     +            NFRONT(2),W,JBUFL,W(JBUFR),JBUFU,IW,ISIZE(3),W(FA),
     +            W(FRHS),IW(LHED),IW(KHED),IW(KPIV),IW(KPVLNK),
     +            IW(LDEST),IW(KDEST),ICNTL,CNTL,IFILE,IREC,ISIZE,MKEY,
     +            NUMBLK,IELL,NELL,KFRNT,LFRNT,ISAVE(41),PIVBLK,
     +            BZERO,INFO,RINFO)
C If the matrix has been found to be singular, return zero determinant.
      IF (INFO(1).EQ.1 .OR. INFO(1).EQ.-14) THEN
         RINFO(1) = ZERO
         INFO(2) = 0
      END IF
C Preserve IREC and MKEY.
      DO 80 I = 1,3
         ISAVE(3+I) = IREC(I)
         ISAVE(9+I) = MKEY(I)
   80 CONTINUE
C  Check for errors from MA42F/FD
      IF (INFO(1).LT.0) GO TO 260
      IF (INFO(1).EQ.4) THEN
C Set a flag in ISAVE(16)
         ISAVE(16) = 4
         LHED = ISAVE(32)
         NFR = 0
C Loop over variables already in the front, adding up the number
C of non-fully summed variables.
         DO 90 L = 1,LFRNT
            MFR = IW(LHED+L-1)
            IF (LAST(MFR).GE.IELL) THEN
               NFR = NFR + 1
               LAST(MFR) = -LAST(MFR)
            END IF
   90    CONTINUE
         KFRNT = KFRNT - (LFRNT-NFR)
         LFRNT = NFR
         ISAVE(39) = KFRNT
         ISAVE(40) = LFRNT
         GO TO 120
      END IF
      ISAVE(17) = IELL
      ISAVE(39) = KFRNT
      ISAVE(40) = LFRNT
      IF (INFO(1).EQ.5 .OR. INFO(1).EQ.6) GO TO 220
      IF (IELL.EQ.NELL) THEN
C Forward elimination pass complete.
C Zero array in which solution is returned.
         IF (NRHS.GT.0) THEN
C Zero array in which solution is returned (only do this
C if NDF.ne.INFO(3), since in this case there are some
C variables I.lt.NDF which never appear and we want X(I)=0.0
C for such variables).
C Also must zero X if problem been found to be singular.
            IF (NDF.NE.INFO(3) .OR. INFO(2).EQ.0) THEN
               DO 110 L = 1,NRHS
                  DO 100 K = 1,NDF
                     X(K,L) = ZERO
  100             CONTINUE
  110          CONTINUE
            END IF
            IF (IFILE(1).EQ.0) THEN
               DIMBUF = JBUFU
            ELSE
               DIMBUF = JBUFU + NFRONT(1)*NFRONT(2)
            END IF
            IF (IFILE(3).EQ.0) THEN
               DIMIBF = ISIZE(3)
            ELSE
               DIMIBF = LIW
            END IF
C Call subroutine to perform back-substitution on right hand sides.
            CALL MA42DD(1,NRHS,LX,X,.FALSE.,1,1,Y,W(JBUFR),DIMBUF,IW,
     +                  DIMIBF,W(FRHS),NFRONT(2),LP,NRHS,IFILE,IREC,
     +                  ISIZE,MKEY,INFO(1))
            IF (INFO(1).LT.0) THEN
               IF (LP.GE.0) WRITE (LP,FMT=9110) INFO(1)
               GO TO 260
            END IF
         END IF
C Reset ISAVE(16), ISAVE(17), ISAVE(39), ISAVE(40) before
C final return
         GO TO 260
      END IF
      GO TO 270
C
C Perform symbolic assembly and elimination only.  This is only
C done if factors have been lost due to front size being
C exceeded at a previous step.
  120 CALL MA42JD(NVAR,IVAR,NDF,LAST,NMAXE,IW,ICNTL,ISAVE,INFO)
C Return if all elt/equ have been considered.
      IF (IELL.EQ.NELL) GO TO 210
      GO TO 270
C
C *** Error returns ***
  130 INFO(1) = -2
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9100) INFO(1),IELL
         WRITE (LP,FMT=9020) NVAR
         IF (NVAR.GT.NMAXE) WRITE (LP,FMT=9220) NMAXE
      END IF
      GO TO 260
  140 INFO(1) = -5
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9110) INFO(1)
         WRITE (LP,FMT=9040) LX,ISAVE(30)
      END IF
      GO TO 260
  150 INFO(1) = -6
      IELL = 1
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9100) INFO(1),IELL
         WRITE (LP,FMT=9030) LW,LLW
      END IF
      GO TO 260
  160 INFO(1) = -7
      IELL = 1
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9100) INFO(1),IELL
         WRITE (LP,FMT=9050) LIW,LLIW
      END IF
      GO TO 260
  170 INFO(1) = -8
      IF (LP.GT.0) THEN
         IF (NMAXE.LE.0) THEN
            WRITE (LP,FMT=9110) INFO(1)
            WRITE (LP,FMT=9060) NMAXE
         ELSE
            WRITE (LP,FMT=9100) INFO(1),IELL
            WRITE (LP,FMT=9000) ISAVE(21),NMAXE
         END IF
      END IF
      GO TO 260
  180 INFO(1) = -9
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9100) INFO(1),IELL
         WRITE (LP,FMT=9070) ISAVE(22),NRHS
      END IF
      GO TO 260
  190 INFO(1) = -10
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9110) INFO(1)
         WRITE (LP,FMT=9080) NRHS,LRHS
      END IF
      GO TO 260
  200 INFO(1) = -11
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9100) INFO(1),IELL
         IF (ISAVE(23).GT.0 .AND. ISAVE(24).GT.0) WRITE (LP,
     +       FMT=9090) ISAVE(23),ISAVE(24),NFRONT
         IF (ISAVE(23).LE.0) WRITE (LP,FMT=9170) ISAVE(23)
         IF (ISAVE(24).LE.0) WRITE (LP,FMT=9180) ISAVE(24)
      END IF
      GO TO 260
  210 INFO(1) = -12
      IF (LP.GT.0) THEN
         IF (NFRONT(1).GE.IW(1) .AND. NFRONT(2).GE.IW(2)) WRITE (LP,
     +       FMT=9270) NFRONT
      END IF
      NFRONT(1) = IW(1)
      NFRONT(2) = IW(2)
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9110) INFO(1)
         WRITE (LP,FMT=9130) NFRONT
      END IF
C We can also use the data in IW(3), IW(4), IW(5) to supply some info.
C on size of files required.
       INFO(4) = MAX(INFO(4),IW(3))
       INFO(5) = MAX(INFO(5),IW(3) + NRHS*NDF)
       INFO(6) = MAX(INFO(6),IW(4))
       INFO(7) = MAX(INFO(7),IW(5))
C In this case, ISAVE(16), ISAVE(17), ISAVE(39), and ISAVE(40)
C have already been reset so jump directly to return.
      GO TO 270
  220 IF (INFO(1).EQ.5 .AND. ISAVE(19).GT.0) ISAVE(19) = -ISAVE(19)
      IF (INFO(1).EQ.6 .AND. ISAVE(20).GT.0) ISAVE(20) = -ISAVE(20)
      IF (IELL.EQ.NELL) THEN
         IF (INFO(1).EQ.5) INFO(1) = -16
         IF (INFO(1).EQ.6) INFO(1) = -17
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9110) INFO(1)
            LFL = INFO(10)*ISIZE(1)
            IF (INFO(10).GT.NUMBLK(1) .AND. ISAVE(20).LT.0) WRITE (LP,
     +          FMT=9140) INFO(5),ISIZE(1),LFL
            IF (ISAVE(19).LT.0 .AND. INFO(5).GT.ISIZE(1)) WRITE (LP,
     +          FMT=9190) INFO(5),ISIZE(1),LFL
            LFL = INFO(11)*ISIZE(2)
            IF (INFO(11).GT.NUMBLK(2) .AND. ISAVE(20).LT.0) WRITE (LP,
     +          FMT=9150) INFO(6),ISIZE(2),LFL
            IF (ISAVE(19).LT.0 .AND. INFO(6).GT.ISIZE(2)) WRITE (LP,
     +          FMT=9200) INFO(6),ISIZE(2),LFL
            LFL = INFO(12)*ISIZE(3)
            IF (INFO(12).GT.NUMBLK(3) .AND. ISAVE(20).LT.0) WRITE (LP,
     +          FMT=9160) INFO(7),ISIZE(3),LFL
            IF (ISAVE(19).LT.0 .AND. INFO(7).GT.ISIZE(3)) WRITE (LP,
     +          FMT=9210) INFO(7),ISIZE(3),LFL
         END IF
         GO TO 260
      ELSE
         ISAVE(17) = IELL
         GO TO 270
      END IF
  230 INFO(1) = -15
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9110) INFO(1)
         WRITE (LP,FMT=9120) ISAVE(30),NDF
      END IF
      GO TO 260
  240 INFO(1) = -18
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9110) INFO(1)
         WRITE (LP,FMT=9010) NRHS
      END IF
      GO TO 260
  250 INFO(1) = -19
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9110) INFO(1)
         IF (ISAVE(31).EQ.1) THEN
            WRITE (LP,FMT=9250) I
         ELSE
            IF (LENBUF(1).LE.0) WRITE (LP,FMT=9230) LENBUF(1)
            IF (LENBUF(3).LE.0) WRITE (LP,FMT=9240) LENBUF(3)
            IF (LENBUF(2).LT.0) WRITE (LP,FMT=9260) LENBUF(2)
         END IF
      END IF
      GO TO 260
C
C Reset ISAVE(17), ISAVE(39), ISAVE(40) to 0
C and reset ISAVE(16) to 2 (so MA42A/AD and MA42B/BD
C may be recalled without recalling MA42I/ID).
C Also restore ISAVE(19) and ISAVE(20)
C
  260 ISAVE(17) = 0
      ISAVE(39) = 0
      ISAVE(40) = 0
      ISAVE(16) = 2
      ISAVE(19) = ABS(ISAVE(19))
      ISAVE(20) = ABS(ISAVE(20))
C
  270 IF (INFO(1).LT.0 .OR. INFO(1).GT.3) ISAVE(38) = INFO(1)
  280 RETURN
 9000 FORMAT (7X,'NMAXE has been changed from ',I8,' to ',I8)
 9010 FORMAT (7X,'Number of right hand sides (NRHS) is',I8)
 9020 FORMAT (7X,'Number of variables in elt/eqn is ',I8)
 9030 FORMAT (7X,'Length of real workspace too small.',
     +       /7X,'Increase LW from ',I8,' to ',I8)
 9040 FORMAT (7X,'First dimension of array X too small.',
     +       /7X,'Increase LX from ',I8,' to ',I8)
 9050 FORMAT (7X,'Length of integer workspace too small.',
     +       /7X,'Increase LIW from ',I8,' to ',I8)
 9060 FORMAT (7X,'NMAXE is equal to ',I8)
 9070 FORMAT (7X,'Change in number of right-hand sides from ',I8,' to ',
     +       I8)
 9080 FORMAT (7X,'You wish to solve for ',I8,' right hand sides',
     +       /7X,'however arrays allow for ',I8,' right hand sides.',
     +       /7X,'Check values of NRHS and LRHS')
 9090 FORMAT (7X,'NFRONT(1) and NFRONT(2) have been changed',
     +       /7X,'from ',I8,' and ',I8,' to ',I8,' and ',I8,
     +       ' respectively.')
 9100 FORMAT (' ***** Error return from MA42B/BD *****  INFO(1) = ',I3,
     +       /' ***** Error discovered after input of elt/eqn *****',I8)
 9110 FORMAT (' ***** Error return from MA42B/BD *****  INFO(1) = ',I3)
 9120 FORMAT (' NDF has been changed from ',I8,' to ',I8)
 9130 FORMAT (7X,'Lower bound on size of NFRONT required is ',2I8)
 9140 FORMAT (7X,'Total storage needed for U-factors is ',I8,
     +       /7X,'At present LENBUF(1) size of ',I8,
     +       /7X,'LENFLE(1) must be at least ',I8)
 9150 FORMAT (7X,'Total storage needed for L-factors is ',I8,
     +       /7X,'At present LENBUF(2) size of ',I8,
     +       /7X,'LENFLE(2) must be at least ',I8)
 9160 FORMAT (7X,'Total storage needed for integers is ',I8,
     +       /7X,'At present LENBUF(3) size of ',I8,
     +       /7X,'LENFLE(3) must be at least ',I8)
 9170 FORMAT (7X,'Front size (NFRONT(1)) non-positive equal to ',I8)
 9180 FORMAT (7X,'Front size (NFRONT(2)) non-positive equal to ',I8)
 9190 FORMAT (7X,'Total storage needed for U-factors is ',I8,
     +       /7X,
     +          'Set LENBUF(1) to this or use direct access data sets.',
     +       /7X,'At present LENBUF(1) size of ',I8,
     +       /7X,'LENFLE(1) must be at least ',I8)
 9200 FORMAT (7X,'Total storage needed for L-factors is ',I8,
     +       /7X,
     +          'Set LENBUF(2) to this or use direct access data sets.',
     +       /7X,'At present LENBUF(2) size of ',I8,
     +       /7X,'LENFLE(2) must be at least',I8)
 9210 FORMAT (7X,'Total storage needed for integers is ',I8,
     +       /7X,
     +          'Set LENBUF(3) to this or use direct access data sets.',
     +       /7X,'At present LENBUF(3) size of ',I8,
     +       /7X,'LENFLE(3) must be at least',I8)
 9220 FORMAT (7X,'But first dimension of element arrays is only ',I8)
 9230 FORMAT (7X,'Non-positive buffer size (LENBUF(1)) of ',I8)
 9240 FORMAT (7X,'Non-positive buffer size (LENBUF(3)) of ',I8)
 9250 FORMAT (7X,'LENBUF(',I1,') has been changed since the call to',
     +          ' MA42P/PD')
 9260 FORMAT (7X,'Negative buffer size (LENBUF(2)) of ',I8)
 9270 FORMAT (7X,'Frontsize of ',2I8,' is too small',
     +       /7X,'but unable',
     +          ' to provide more useful information on a suitable',
     +       /7X,'value for NFRONT.')
      END
C*********************************************************************
      SUBROUTINE MA42CD(TRANS,NRHS,LX,B,X,LW,W,LIW,IW,ICNTL,ISAVE,INFO)
C
C *** Reference MA42 suite ***
C *** Copyright Rutherford Appleton Laboratory  November 1992 ***
C *** Any problems contact Iain S. Duff or Jennifer A. Scott
C     at Atlas Centre, Rutherford Appleton Laboratory ***
C *** Although every effort has been made to ensure robustness and
C     reliability of the subroutines in this MA42 suite,  we
C     disclaim any liability arising through the use or misuse of
C     any of the subroutines in the MA42 suite ***
C
C This subroutine either solves further systems of the form
C      AX = B
C or it solves systems of the form
C      A(T)X = B
C using the factors obtained using MA42B/BD.
C It calls subroutine MA42E/ED to perform the forward elimination
C and then MA42D/DD to perform the back-substitution.
C
C     Argument list. * indicates the argument is changed by the routine.
C
C  TRANS  - Logical variable.  If TRANS = .TRUE. A(T)X = B is to be
C           solved, and if TRANS = .FALSE. AX = B is to be solved.
C  NRHS   - Integer variable.  Number of right hand sides.
C  LX     - Integer variable.  Leading dimension of arrays B and X.
C           Must be at least as large as NDF as output from MA42A/AD.
C *B      - Real (DP) array of dimensions LX by NRHS. Must be
C           set to hold the right hand sides on input.  Is altered by
C           the subroutine.
C *X      - Real (DP) array of dimensions LX by NRHS.
C           On successful exit X(I,J) will hold the Ith component of
C           the solution to system J.
C  LW     - Integer variable. Length of array W.
C           If no direct access data sets used (MA42P/PD not called),
C           LW must be at least LENBUF(1)+LENBUF(2)+
C           NRHS*max(NFRONT(1), NFRONT(2)).  Otherwise,
C           LW must be at least max{LENBUF(1)*INFO(19),
C           LENBUF(2)*INFO(20)}+NRHS*(INFO(8), INFO(9)).
C           (INFO values as output from the final call to MA42B/BD).
C           (These INFO values were copied into ISAVE at end of
C           MA42B/BD).
C *W      - Real (DP) array of length LW. If direct access
C           files not used, the first LENBUF(1)+LENBUF(2)
C           entries of W must be unchanged since
C           the last call to MA42B/BD and these entries are unchanged
C           by MA42C/CD.  Otherwise W is used as workspace.
C  LIW -    Integer variable. Length of array IW. Must be at least
C           as great as LENBUF(3)*INFO(21).
C *IW     - Integer array of length LIW. If direct access
C           files not used, IW must be unchanged since
C           the last call to MA42B/BD and is unchanged
C           by MA42C/CD. Otherwise IW is used as workspace.
C  ICNTL  - Integer array of length 8. Holds control parameters.
C           Only ICNTL(1), which holds the stream number for error
C           messages, is accessed by the routine.
C  ISAVE  - Integer array of length 45. Holds parameters which
C           must be preserved between calls to MA42 routines.
C *INFO   - Integer array of length 23. INFO(1) is an error flag.
C           On exit INFO(1) will equal 0 if run is
C           successful. Negative values indicate an error has occurred.
C           Possible negative values are -5, -6, -7, -18, -25, -27.
C           The rest of the array is not accessed.
C
C  Local variables
C
C   IFILE  - Integer array of length 3. IFILE(I) is equivalent to
C            ISAVE(I) (I = 1,2,3).
C   IREC   - Integer array of length 3. IREC(I) is equivalent to
C            ISAVE(3+I) (I = 1,2,3).
C   ISIZE  - Integer array of length 3. ISIZE(I) is equivalent to
C            ISAVE(6+I) (I = 1,2,3).
C   MKEY   - Integer array of length 3. MKEY(I) is equivalent to
C            ISAVE(9+I) (I = 1,2,3).
C   DIMIBF - Integer variable. Used to hold length of integer
C            workspace required by routine.
C   I      - Integer variable. Do loop variable.
C   L      - Integer variable. Do loop variable.
C L1,L2,L3 _ Integer variables. Equivalent to ISIZE(1),ISIZE(2),
C            ISIZE(3), respectively.
C   LBUFR  - Integer variable. Used to hold length of buffer.
C   LLW    - Integer variable. Used to hold length of real
C            workspace required by routine.
C   LP     - Integer variable. Used to hold stream number for
C            error messages (ICNTL(1)).
C   NFMAX  - Integer variable. Used to hold max(NFRONT(1),NFRONT(2))
C            (= max(ISAVE(23),ISAVE(24)).
C   NRHSB  - Integer variable. Equivalent to ISAVE(22). Holds
C            the number of right-hand sides when MA42B/BD was called.
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER LIW,LW,LX,NRHS
      LOGICAL TRANS
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION B(LX,NRHS),W(LW),X(LX,NRHS)
      INTEGER ICNTL(8),INFO(23),ISAVE(45),IW(LIW)
C     ..
C     .. Local Arrays ..
      INTEGER IFILE(3),IREC(3),ISIZE(3),MKEY(3)
C     ..
C     .. Local Scalars ..
      INTEGER DIMIBF,I,L,L1,L2,L3,LBUFR,LLW,LP,NFMAX,NRHSB
C     ..
C     .. External Subroutines ..
      EXTERNAL MA42DD,MA42ED
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX
C     ..
C Immediate return if a previous call to MA42B/BD caused an error
      IF (ISAVE(38).LT.0) GO TO 100
      INFO(1) = 0
      DO 10 I = 1,3
         IFILE(I) = ISAVE(I)
         IF (ICNTL(1).EQ.IFILE(I) .AND. IFILE(I).NE.0) ICNTL(1) = 6
         IREC(I) = ISAVE(3+I)
         ISIZE(I) = ISAVE(6+I)
         MKEY(I) = ISAVE(9+I)
   10 CONTINUE
      LP = ICNTL(1)
C Perform checks
      IF (NRHS.LT.1) GO TO 70
      NRHSB = ISAVE(22)
C ISAVE(30) holds copy of NDF as output from MA42A/AD.
      IF (LX.LT.ISAVE(30)) GO TO 40
C Initialise solution vector to zero.
      DO 30 L = 1,NRHS
         DO 20 I = 1,ISAVE(30)
            X(I,L) = ZERO
   20    CONTINUE
   30 CONTINUE
C Jump to error return if MA42P/PD was called with ISTRM(2) = 0
C (L-factor was not stored).
C Also jump to error return if MA42B/BD was called with LENBUF(2) = 0.
      IF (ISIZE(2).EQ.0) GO TO 80
      L1 = ISIZE(1)
      L2 = ISIZE(2)
      L3 = ISIZE(3)
C
C      NFMAX = MAX(ISAVE(23),ISAVE(24))
C April 1999: reduced workspace to what is really used.
      NFMAX = MAX(INFO(8),INFO(9))
      IF (ISAVE(31).EQ.-1) THEN
C No direct access data sets used.
         DIMIBF = L3
         IF (LIW.LT.DIMIBF) GO TO 60
         LLW = NRHS*NFMAX + L1 + L2
         IF (LW.LT.LLW) GO TO 50
         IF (TRANS) THEN
C Solve A(T)X = B
            CALL MA42ED(1,NRHS,LX,B,W(L2+1),L1,IW,DIMIBF,W(L2+L1+1),
     +                 NFMAX,LP,NRHSB,IFILE,ISIZE,MKEY,INFO(1))
            IF (INFO(1).LT.0) GO TO 90
            CALL MA42DD(2,NRHS,LX,X,.TRUE.,LX,NRHS,B,W,L2,IW,DIMIBF,
     +                 W(L2+L1+1),NFMAX,LP,NRHSB,IFILE,IREC,ISIZE,MKEY,
     +                 INFO(1))
            IF (INFO(1).LT.0) GO TO 90
         ELSE
C Solve AX = B
            CALL MA42ED(2,NRHS,LX,B,W,L2,IW,DIMIBF,W(L2+L1+1),NFMAX,LP,
     +                 NRHSB,IFILE,ISIZE,MKEY,INFO(1))
            IF (INFO(1).LT.0) GO TO 90
            CALL MA42DD(1,NRHS,LX,X,.TRUE.,LX,NRHS,B,W(L2+1),L1,IW,
     +                 DIMIBF,W(L2+L1+1),NFMAX,LP,NRHSB,IFILE,IREC,
     +                 ISIZE,MKEY,INFO(1))

            IF (INFO(1).LT.0) GO TO 90
         END IF
      ELSE
C Direct access data sets used.
C         DIMIBF = L3 + 5 + ISAVE(23) + ISAVE(24)
         DIMIBF = L3 + 5 + INFO(8) + INFO(9)
         IF (LIW.LT.DIMIBF) GO TO 60
C April 1999: found bug. The workspace could be too small
C if MA42B/BD was called with a large number of right-hand sides
C Note that INFO(22) (largest pivot block size) was not saved in ISAVE
C so we have to use INFO(22) ... user must pass INFO unchanged
C         LBUFR = MAX(L1,L2) + ISAVE(23)*ISAVE(24)
         LBUFR = MAX(L1,L2) + INFO(22)*MAX(INFO(8),INFO(9)+ISAVE(22))
         LLW = NRHS*NFMAX + LBUFR
         IF (LW.LT.LLW) GO TO 50
         IF (TRANS) THEN
C Solve A(T)X = B
            CALL MA42ED(1,NRHS,LX,B,W,LBUFR,IW,DIMIBF,W(LBUFR+1),NFMAX,
     +                 LP,NRHSB,IFILE,ISIZE,MKEY,INFO(1))
            IF (INFO(1).LT.0) GO TO 90
            CALL MA42DD(2,NRHS,LX,X,.TRUE.,LX,NRHS,B,W,LBUFR,IW,DIMIBF,
     +                 W(LBUFR+1),NFMAX,LP,NRHSB,IFILE,IREC,ISIZE,MKEY,
     +                 INFO(1))
            IF (INFO(1).LT.0) GO TO 90
         ELSE
C Solve AX = B
            CALL MA42ED(2,NRHS,LX,B,W,LBUFR,IW,DIMIBF,W(LBUFR+1),NFMAX,
     +                 LP,NRHSB,IFILE,ISIZE,MKEY,INFO(1))
            IF (INFO(1).LT.0) GO TO 90
            CALL MA42DD(1,NRHS,LX,X,.TRUE.,LX,NRHS,B,W,LBUFR,IW,DIMIBF,
     +                 W(LBUFR+1),NFMAX,LP,NRHSB,IFILE,IREC,ISIZE,MKEY,
     +                 INFO(1))
            IF (INFO(1).LT.0) GO TO 90
         END IF
      END IF
      GO TO 100
C
C **** Error returns ****
   40 INFO(1) = -5
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1)
         WRITE (LP,FMT=9060) LX,ISAVE(30)
      END IF
      GO TO 100
   50 INFO(1) = -6
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1)
         WRITE (LP,FMT=9030) LW,LLW
      END IF
      GO TO 100
   60 INFO(1) = -7
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1)
         WRITE (LP,FMT=9020) LIW,DIMIBF
      END IF
      GO TO 100
   70 INFO(1) = -18
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1)
         WRITE (LP,FMT=9010) NRHS
      END IF
      GO TO 100
   80 INFO(1) = -27
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1)
         IF (ISAVE(31).EQ.1) WRITE (LP,FMT=9040) IFILE(2)
         IF (ISAVE(31).EQ.-1) WRITE (LP,FMT=9050) ISIZE(2)
      END IF
      GO TO 100
   90 IF (LP.GT.0) WRITE (LP,FMT=9000) INFO(1)
C
  100 RETURN
 9000 FORMAT (' ***** Error return from MA42C/CD ***** INFO(1) = ',I3)
 9010 FORMAT (7X,'Number of right hand sides (NRHS) is ',I8)
 9020 FORMAT (7X,'Length of integer workspace too small.',
     +       /7X,'Increase LIW from ',I8,' to ',I8)
 9030 FORMAT (7X,'Length of real workspace too small.',
     +       /7X,'Increase LW from ',I8,' to ',I8)
 9040 FORMAT (7X,'No record of L-factor being stored.',
     +       /7X,'MA42P/PD was called with ISTRM(2) set to ',I3)
 9050 FORMAT (7X,'No record of L-factor being stored.',
     +       /7X,'MA42B/BD was called with LENBUF(2) set to ',I8)
 9060 FORMAT (7X,'First dimension of arrays B and X too small.',
     +       /7X,'Increase LX from ',I8,' to ',I8)
      END
C**********************************************************************
      SUBROUTINE MA42DD(IND,NRHS,NDF,X,LYINBF,LY1,LY2,Y,BUFR,DIMBUF,
     +                  IBUFR,DIMIBF,W,LDW,LP,NRHSB,IFILE,IREC,ISIZE,
     +                  MKEY,INFO)
C
C     MA42D/DD is the back-substitution pass of the frontal method
C     equation solver.  Sparse matrix equations of the form
C       (a)       A*X = B   or  (b)   A(T)*X = B
C     where A is a matrix and X and B are vectors, are solved by
C     decomposing A into the product L*U of lower and upper triangular
C     sparse matrices, and then (a) Y=(L**-1)*B and X=(U**-1)*Y
C     or (b) Y=(U(T)**-1)*B and X=(L(T)**-1)*Y  are
C     calculated. MA42D/DD calculates (a) X=(U**-1)*Y or
C     (b) X=(L(T)**-1)*Y.
C
C     Argument list. * indicates the argument is changed by the routine.
C
C  IND    - Integer variable.  Has value 1 or 2.
C           If IND = 1 use U factor (solving A*X = B).
C           If IND = 2 use L factor (solving A(T)*X = B).
C  NRHS   - Integer variable.  Number of right hand sides.
C  NDF    - Integer variable. Total number of variables.
C  LYINBF - Logical control variable.
C *Y      - Real (DP) array of dimensions LY1,LY2.
C           If LYINBF is .TRUE. Y contains the partial solution Y.
C           Otherwise Y is not accessed.
C *X      - Real (DP) array of dimensions NDF,NRHS.
C           On exit, holds the solution vector.
C           Workspace array used to index variables.
C *BUFR   - Real (DP) array of dimension DIMBUF.
C           Workspace used to hold buffers. (Unless IFILE(IND)=0, in
C           which case BUFR holds information from MA42B/BD).
C  DIMBUF - Integer variable. Dimension of BUFR.
C *IBUFR  - Integer array of dimension DIMIBF.
C           Workspace used to hold integer buffers. (Unless IFILE(3)=0,
C           in which case IBUFR holds information from MA42B/BD).
C  DIMIBF - Integer variable. Dimension of IBUFR.
C *W      - Real (DP) array of dimensions LDW,NRHS.
C           Used by the routine as workspace.
C  LDW    - Integer variable. Leading dimension of array W.
C  LP     - Integer variable. Stream number for error messages.
C  NRHSB  - Integer variable. Holds number of right-hand sides
C           when MA42B/BD was called.
C  IFILE  - Integer array of length 3. Stream numbers for direct
C           access data sets.
C  IREC   - Integer array of length 3. Pointers to first free space
C           in buffers.
C  ISIZE  - Integer array of length 3. Length of buffers.
C  MKEY   - Integer array of length 3. Number of records written
C           to direct access data sets.
C  INFO   - Integer variable. Negative value on exit indicates
C           fatal error. Possible negative value -25 caused by
C           failure in direct access read.
C
C Local variables
C
C FIL     - Integer variable. Stream number for real direct
C           access data set (IFILE(IND)).
C I       - Integer variable.  Do loop variable.
C IFIL    - Integer variable. Stream number for integer direct
C           access data set (IFILE(3)).
C IKEY    _ Integer variable. Points to the record
C           in the integer direct access data set being read.
C ILNGTH  - Integer variable. Length of integer record.
C IPNT    - Integer variable. Points to end of record in BUFR.
C IR1     - Integer variable. Points to the first integer in the
C           integer information for the pivot block.
C IR2     - Integer variable. Points to the start of the data
C           which has been read into IBUFR.
C IRECD   - Integer variable. Points to end of record in IBUFR.
C J1,J2   - Integer variables. Used in computing NREAD.
C JFLAG   - Integer variable. Used as error flag for call to MA42L/LD.
C JFRNT1  - Integer variable. Used to hold KFRNT1 or LFRNT1.
C           Use of JFRNT1 simplifies notation (allows
C           same notation for A*X=B and A(T)X=B)
C JR1     - Integer variable. Used to point to entries in pivot
C           block.
C JR2     - Integer variable. Points to the start of the data
C           which has been read into BUFR.
C L       - Integer variable. Do loop variable. Loop over size
C           of pivot block (NPIV).
C K1,K2,K3 - Integer variables. Used to point to entries in pivot
C           block.
C K1PK,K2PK,K3PK - Integer variables. In do loop, K1PK = K1 + K
C           where K is do loop variable. K2PK, K3PK similarly.
C KEY     _ Integer variable. Points to the record
C           in the real direct access data set being read.
C KFRNT   - Integer variable. Number of rows in front.
C KFRNT1  - Integer variable. Holds KFRNT-NPIV.
C L       - Integer variable. Do loop variable. Loop over number
C           of right hand sides (NRHS).
C LBUFR   - Integer variable. Length of real buffer (ISIZE(IND)).
C LCO     - Integer variable. Local column index.
C LENGTH  - Integer variable. Length of real record.
C LFRNT   - Integer variable. Number of columns in front.
C LFRNT1  - Integer variable. Holds LFRNT-NPIV.
C LIBUFR  - Integer variable. Length of integer buffer (ISIZE(3)).
C NLOOP   - Integer variable. Length of do loop. Each pass
C           through loop processes one integer record.
C NPIV    - Integer variable. Size of pivot block.
C NPIV1   - Integer variable. holds (NPIV*(NPIV+1)/2).
C NREAD   - Integer variable. Number of direct access data set records
C           which must be read in to obtain the required reals or
C           integers for a pivot block.
C NREC    - Integer variable. Do loop variable. Ranges from
C           1 to NLOOP.
C
C     .. Parameters ..
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER DIMBUF,DIMIBF,IND,INFO,LDW,LP,LY1,LY2,NDF,NRHS,NRHSB
      LOGICAL LYINBF
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION BUFR(DIMBUF),W(LDW,NRHS),X(NDF,NRHS),Y(LY1,LY2)
      INTEGER IBUFR(DIMIBF),IFILE(3),IREC(3),ISIZE(3),MKEY(3)
C     ..
C     .. External Subroutines ..
      EXTERNAL DGEMM,DGEMV,DTPSV,MA42LD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX
C     ..
C     .. Local Scalars ..
      INTEGER FIL,I,IFIL,IKEY,ILNGTH,IPNT,IR1,IR2,IREAD,IRECD,J1,J2,
     +        JFLAG,JFRNT1,JR1,JR2,K,K1,K1PK,K2,K2PK,K3,K3PK,KEY,KFRNT,
     +        KFRNT1,KRO,L,LBUFR,LCO,LENGTH,LFRNT,LFRNT1,LIBUFR,NLOOP,
     +        NPIV,NPIV1,NREAD,NREC
C     ..
C LBUFR and LIBUFR hold the lengths of the real and integer buffers.
      LBUFR = ISIZE(IND)
      LIBUFR = ISIZE(3)
C MKEY(IND) points to the last record written to the direct access
C data set (IND=1,2,3).
      KEY = MKEY(IND)
      IKEY = MKEY(3)
      FIL = IFILE(IND)
      IFIL = IFILE(3)
      IF (FIL.NE.0 .AND. KEY.NE.0) THEN
C Records are to be read in the reverse order of the output order
C (last out, first in).
C Read last record in direct access data set on stream FIL into last
C space of length LBUFR in BUFR.
         CALL MA42LD(2,FIL,KEY,BUFR(DIMBUF-LBUFR+1),LBUFR,IBUFR,LIBUFR,
     +               LP,JFLAG)
         IF (JFLAG.LT.0) GO TO 150
      END IF
      IF (IFIL.NE.0) THEN
C Read last record in direct access data set on stream IFIL into last
C space of length LIBUFR in IBUFR.
         CALL MA42LD(-2,IFIL,IKEY,BUFR,LBUFR,IBUFR(DIMIBF-LIBUFR+1),
     +               LIBUFR,LP,JFLAG)
         IF (JFLAG.LT.0) GO TO 150
      ELSE
C Set IKEY to 0 to indicate no records to be read on stream IFIL.
         IKEY = 0
      END IF
C Set JR2 and IR2.
      JR2 = DIMBUF - LBUFR
      IR2 = DIMIBF - LIBUFR
C IPNT points to the end of the last record in BUFR.
      IPNT = IREC(IND) - 1
      IF (IREC(IND).EQ.1) IPNT = ISIZE(IND)
C Set IRECD to point to the end of the last record in IBUFR.
C (skip the zero flag in position IREC(3)-1 which is used in
C the forward substitution).
      IRECD = IREC(3) - 2
C Each pass through this loop processes one pivot block.
      NLOOP = MAX(1,MKEY(3))*LIBUFR
      DO 140 NREC = 1,NLOOP
C Return if physical record exhausted.
         IF (IRECD.EQ.0 .AND. IKEY.EQ.0) GO TO 160
C If all records in the integer buffer have been dealt with, read in
C another record from the integer direct access data set.
         IF (IRECD.EQ.0) THEN
C If no room in IBUFR to read in another record, reset IR2
C
            IF (IR2-LIBUFR.LT.0) IR2 = DIMIBF
C Read into last available space of length LIBUFR in IBUFR.
            CALL MA42LD(-2,IFIL,IKEY,BUFR,LBUFR,IBUFR(IR2-LIBUFR+1),
     +                  LIBUFR,LP,JFLAG)
            IF (JFLAG.LT.0) GO TO 150
C Reset IRECD and update IR2.
            IRECD = LIBUFR
            IR2 = IR2 - LIBUFR
         END IF
C ILNGTH is the number of entries in the integer information record
C for the pivot block.
         ILNGTH = IBUFR(IR2+IRECD)
C Calculate how many more (NREAD) direct access records must be
C read in to obtain all the integer information for the pivot block.
         NREAD = 0
         J1 = ABS(ILNGTH) - IRECD
         J2 = J1/LIBUFR
         IF (J1.GT.0) THEN
C Further records must be read from the integer direct access data set.
            NREAD = 1 + J2
            IF (J2*LIBUFR.EQ.J1) NREAD = NREAD - 1
C If there is not room to accommodate the incoming NREAD records
C in the front of IBUFR, copy existing data in IBUFR to the
C end of IBUFR before calling MA42L/LD.
            IF (IR2-NREAD*LIBUFR.LT.0) THEN
               DO 10 I = IRECD,1,-1
c               DO 10 I = 1,IRECD
                  IBUFR(DIMIBF-IRECD+I) = IBUFR(IR2+I)
   10          CONTINUE
C Reset IR2.
               IR2 = DIMIBF - IRECD
            END IF
            DO 20 IREAD = 1,NREAD
               CALL MA42LD(-2,IFIL,IKEY,BUFR,LBUFR,
     +                     IBUFR(IR2-IREAD*LIBUFR+1),LIBUFR,LP,JFLAG)
               IF (JFLAG.LT.0) GO TO 150
   20       CONTINUE
C Reset IRECD.
            IRECD = LIBUFR - (J1-J2*LIBUFR)
            IF (J2*LIBUFR.EQ.J1) IRECD = 0
         ELSE
C Reduce IRECD by length of record.
            IRECD = IRECD - ABS(ILNGTH)
         END IF
C Update IR2.
         IR2 = IR2 - NREAD*LIBUFR
C Jump to next record if record was block associated with singular
C part of matrix.
         IF (ILNGTH.LT.0) GO TO 140
C Required integers now in IBUFR.
C IR1 points to the first integer in the integer information
C on the pivot block.
         IR1 = IR2 + IRECD + 1
         NPIV = IBUFR(IR1+1)
         KFRNT = IBUFR(IR1+2)
         LFRNT = IBUFR(IR1+3)
         LFRNT1 = LFRNT - NPIV
         KFRNT1 = KFRNT - NPIV
         NPIV1 = (NPIV* (NPIV+1))/2
         IF (IND.EQ.1) THEN
            JFRNT1 = LFRNT1
            K1 = IR1 + 3 + KFRNT1
            K2 = IR1 + 3 + KFRNT
            K3 = IR1 + 3 + KFRNT + LFRNT1
            LENGTH = NPIV* (LFRNT1+NRHSB) + NPIV1
         ELSE
            JFRNT1 = KFRNT1
            K1 = IR1 + 3 + KFRNT + LFRNT1
            K2 = IR1 + 3
            K3 = IR1 + 3 + KFRNT1
            LENGTH = NPIV*KFRNT1 + NPIV1
         END IF
C
C Calculate how many more (NREAD) direct access records must be
C read in to obtain all the reals for the pivot block.
         NREAD = 0
         J1 = LENGTH - IPNT
         J2 = J1/LBUFR
         IF (J1.GT.0) THEN
C Further records must be read from the real direct access data set.
            NREAD = 1 + J2
            IF (J2*LBUFR.EQ.J1) NREAD = NREAD - 1
C If there is not room to accommodate the incoming NREAD records
C in the front of BUFR, copy any existing data in BUFR to the
C end of BUFR before calling MA42L/LD.
            IF (JR2-NREAD*LBUFR.LT.0) THEN
               DO 30 I = IPNT,1,-1
c               DO 30 I = 1,IPNT
                  BUFR(DIMBUF-IPNT+I) = BUFR(JR2+I)
   30          CONTINUE
               JR2 = DIMBUF - IPNT
            END IF
            DO 40 IREAD = 1,NREAD
               CALL MA42LD(2,FIL,KEY,BUFR(JR2-IREAD*LBUFR+1),LBUFR,
     +                     IBUFR,LIBUFR,LP,JFLAG)
               IF (JFLAG.LT.0) GO TO 150
   40       CONTINUE
C Reset IPNT
            IPNT = LBUFR - (J1-J2*LBUFR)
            IF (J2*LBUFR.EQ.J1) IPNT = 0
         ELSE
C Reduce IPNT by length of record.
            IPNT = IPNT - LENGTH
         END IF
C Update JR2
         JR2 = JR2 - NREAD*LBUFR
C  Required rows/cols now in BUFR.
C
         IF (LYINBF) THEN
C Copy the partial solution vector Y into the work array W
            DO 60 L = 1,NRHS
               DO 50 K = 1,NPIV
                  K1PK = K1 + K
                  KRO = IBUFR(K1PK)
                  W(K,L) = Y(KRO,L)
   50          CONTINUE
   60       CONTINUE
C
         ELSE IF (IND.EQ.1) THEN
C JR1 points to the first entry in the partial solution corresponding
C to the pivot block.
            JR1 = JR2 + IPNT + 1 + JFRNT1*NPIV + NPIV1
            DO 80 L = 1,NRHSB
               DO 70 K = 1,NPIV
                  W(K,L) = BUFR(JR1)
                  JR1 = JR1 + 1
   70          CONTINUE
   80       CONTINUE
         END IF
C
         IF (JFRNT1.NE.0) THEN
C Put solution vector X into end of work array W.
            DO 100 L = 1,NRHS
               DO 90 K = 1,JFRNT1
                  K2PK = K2 + K
                  LCO = IBUFR(K2PK)
                  W(K+NPIV,L) = X(LCO,L)
   90          CONTINUE
  100       CONTINUE
C
C JR1 points to the first entry in the pivot block
            JR1 = JR2 + IPNT + 1
            IF (NRHS.NE.1) CALL DGEMM('T','N',NPIV,NRHS,JFRNT1,-ONE,
     +                                BUFR(JR1),JFRNT1,W(1+NPIV,1),LDW,
     +                                ONE,W,LDW)
            IF (NRHS.EQ.1) CALL DGEMV('T',JFRNT1,NPIV,-ONE,BUFR(JR1),
     +                                JFRNT1,W(1+NPIV,1),1,ONE,W,1)
         END IF
C
         JR1 = JR2 + IPNT + 1 + JFRNT1*NPIV
         DO 110 L = 1,NRHS
            IF (IND.EQ.1) CALL DTPSV('U','T','U',NPIV,BUFR(JR1),W(1,L),
     +                               1)
            IF (IND.EQ.2) CALL DTPSV('U','T','N',NPIV,BUFR(JR1),W(1,L),
     +                               1)
  110    CONTINUE
C Copy W (first NPIV entries) into the solution vector X
         DO 130 L = 1,NRHS
            DO 120 K = 1,NPIV
               K3PK = K3 + K
               LCO = IBUFR(K3PK)
               X(LCO,L) = W(K,L)
  120       CONTINUE
  130    CONTINUE
C
  140 CONTINUE
      GO TO 160
C Set INFO if MA42L/LD returned an error.
  150 INFO = JFLAG
  160 RETURN
      END
C**********************************************************************
      SUBROUTINE MA42ED(IND,NRHS,NDF,R1,BUFR,DIMBUF,IBUFR,DIMIBF,W,LDW,
     +                  LP,NRHSB,IFILE,ISIZE,MKEY,INFO)
C
C     MA42E/ED is the forward-substitution pass of the frontal method
C     equation solver.  Sparse matrix equations of the form
C       (a)       A*X = B   or  (b)   A(T)*X = B
C     where A is a matrix and X and B are vectors, are solved by
C     decomposing A into the product L*U of lower and upper triangular
C     sparse matrices, and then (a) Y=(L**-1)*B and X=(U**-1)*Y
C     or (b) Y=(U(T)**-1)*B and X=(L(T)**-1)*Y  are
C     calculated. MA42E/ED calculates (a) Y=(L**-1)*B  or
C     (b) Y=(U(T)**-1)*B.
C
C     Argument list. * indicates the argument is changed by the routine.
C
C  IND    - Integer variable.  Has value 1 or 2.
C           If IND = 1 use U factor (solving A(T)*X = B).
C           If IND = 2 use L factor (solving A*X = B).
C  NRHS   - Integer variable.  Number of right hand sides.
C  NDF    - Integer variable. Total number of variables.
C *R1     - Real (DP) array of dimensions NDF, NRHS.
C           On entry holds the right-hand sides (vector B).
C           On exit, holds the partial solution vector (vector Y).
C *BUFR   - Real (DP) array of dimension DIMBUF.
C           Workspace used for buffer. (Unless IFILE(IND)=0, in
C           which case BUFR holds information from MA42B/BD).
C  DIMBUF - Integer variable. Dimension of BUFR.
C *IBUFR  - Integer array of dimension DIMIBF.
C           Workspace used for integer buffer. (Unless IFILE(3)=0, in
C           which case IBUFR holds information from MA42B/BD).
C  DIMIBF - Integer variable. Dimension of IBUFR.
C *W      - Real (DP) array of dimensions LDW,NRHS.
C           Used by the routine as workspace.
C  LDW    - Integer variable. Leading dimension of array W.
C  LP     - Integer variable. Stream number for error messages.
C  NRHSB  - Integer variable. Holds number of right-hand sides
C           when MA42B/BD was called.
C  IFILE  - Integer array of length 3. Stream numbers for direct
C           access data sets.
C  ISIZE  - Integer array of length 3. Length of buffers.
C  MKEY   - Integer array of length 3. Number of records written
C           to direct access data sets.
C  INFO   - Integer variable. Negative value on exit indicates
C           fatal error. Possible negative value -25 caused by
C           failure in direct access read.
C
C Local variables
C
C FIL     - Integer variable. Stream number for real direct
C           access data set (IFILE(IND)).
C I       - Integer variable.  Do loop variable.
C IFIL    - Integer variable. Stream number for integer direct
C           access data set (IFILE(3)).
C IKEY    _ Integer variable. Points to the record
C           in the integer direct access data set being read.
C ILNGTH  - Integer variable. Length of integer record.
C IP1     - Integer variable. Used to point to entries in pivot
C           block.
C IRECD   - Integer variable. Points to start of record in IBUFR.
C IR2     - Integer variable. Points to end of the integer
C           information which has been read into IBUFR.
C ISPACE  - Integer variable. Amount of integer information which has
C           been read into IBUFR but has not yet been used.
C J1,J2   - Integer variables. Used in computing NREAD.
C JFLAG   - Integer variable. Used as error flag for call to MA42L/LD.
C JFRNT   - Integer variable. Holds KFRNT or LFRNT.
C           Use of JFRNT simplifies notation (allows
C           same notation for A*X=B and A(T)X=B).
C JFRNT1  - Integer variable. Used to hold KFRNT1 or LFRNT1.
C           Use of JFRNT1 simplifies notation (allows
C           same notation for A*X=B and A(T)X=B).
C JRECD    - Integer variable. Points to start of record in BUFR.
C JSPACE  - Integer variable. Amount of  information which has
C           been read into BUFR but has not yet been used.
C L       - Integer variable. Do loop variable. Loop over size
C           of pivot block (NPIV).
C K1      - Integer variable. Used to point to entries in pivot
C           block.
C K1PK    - Integer variable. In do loop, K1PK = K1 + K
C           where K is do loop variable.
C KEY     _ Integer variable. Points to the record
C           in the real direct access data set being read.
C KFRNT   - Integer variable. Number of rows in front.
C KFRNT1  - Integer variable. Holds KFRNT-NPIV.
C KRO     - Integer variable. Local row index.
C L       - Integer variable. Do loop variable. Loop over number
C           of right hand sides (NRHS).
C LBUFR   - Integer variable. Length of real buffer (ISIZE(IND)).
C LENGTH  - Integer variable. Length of real record.
C LFRNT   - Integer variable. Number of columns in front.
C LFRNT1  - Integer variable. Holds LFRNT-NPIV.
C LIBUFR  - Integer variable. Length of integer buffer (ISIZE(3)).
C NLOOP   - Integer variable. Length of do loop. Each pass
C           through loop processes one integer record.
C NPIV    - Integer variable. Size of pivot block.
C NPIV1   - Integer variable. holds (NPIV*(NPIV+1)/2).
C NREAD   - Integer variable. Number of direct access data set records
C           which must be read in to obtain the required reals or
C           integers for a pivot block.
C NREC    - Integer variable. Do loop variable. Ranges from
C           1 to NLOOP.
C
C     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D0,ZERO=0.0D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER DIMBUF,DIMIBF,IND,INFO,LDW,LP,NDF,NRHS,NRHSB
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION BUFR(DIMBUF),R1(NDF,NRHS),W(LDW,NRHS)
      INTEGER IBUFR(DIMIBF),IFILE(3),ISIZE(3),MKEY(3)
C     ..
C     .. Local Scalars ..
      INTEGER FIL,I,IFIL,IKEY,ILNGTH,IP1,IR2,IREAD,IRECD,ISPACE,J1,J2,
     +        JFLAG,JFRNT,JFRNT1,JR2,JRECD,JSPACE,K,K1,K1PK,KEY,KFRNT,
     +        KFRNT1,KRO,L,LBUFR,LENGTH,LFRNT,LFRNT1,LIBUFR,NLOOP,NPIV,
     +        NPIV1,NREAD,NREC
C     ..
C     .. External Subroutines ..
      EXTERNAL DGEMM,DGEMV,DTPSV,MA42LD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX
C     ..
C LBUFR and LIBUFR hold the lengths of the real and integer buffers.
      LBUFR = ISIZE(IND)
      LIBUFR = ISIZE(3)
C The records are read in the order they were output.
      KEY = 1
      IKEY = 1
      FIL = IFILE(IND)
      IFIL = IFILE(3)
      IF (FIL.NE.0 .AND. MKEY(IND).NE.0) THEN
         CALL MA42LD(1,FIL,KEY,BUFR,LBUFR,IBUFR,LIBUFR,LP,JFLAG)
         IF (JFLAG.LT.0) GO TO 150
      END IF
      IF (IFIL.NE.0) THEN
         CALL MA42LD(-1,IFIL,IKEY,BUFR,LBUFR,IBUFR,LIBUFR,LP,JFLAG)
C This call to MA42LD reads from direct access data set IFIL into
C buffer IBUFR.  If IFIL is zero, no read is performed and it is
C assumed that the necessary information is already in IBUFR.
         IF (JFLAG.LT.0) GO TO 150
      END IF
C JRECD and IRECD point to the beginning of the first record in
C BUFR and IBUFR.
      JRECD = 1
      IRECD = 1
      JR2 = LBUFR
      JSPACE = LBUFR
      IR2 = LIBUFR
      ISPACE = LIBUFR
C
C Each pass through this loop processes a pivot block.
      NLOOP = MAX(1,MKEY(3))*LIBUFR
      DO 110 NREC = 1,NLOOP
         NREAD = 0
         IF (ISPACE.EQ.0) THEN
C It is necessary to read in an integer record from the integer direct
C access data set to obtain the length (ILGTH) of the integer
C information on the current pivot block.
            IF (IR2+LIBUFR.GT.DIMIBF) THEN
               IR2 = 0
               IRECD = 1
            END IF
            CALL MA42LD(-1,IFIL,IKEY,BUFR,LBUFR,IBUFR(IR2+1),LIBUFR,LP,
     +                  JFLAG)
            IF (JFLAG.LT.0) GO TO 150
C Reset IR2 and ISPACE.
            IR2 = IR2 + LIBUFR
            ISPACE = LIBUFR
         END IF
         ILNGTH = IBUFR(IRECD)
C Return if all records dealt with.
         IF (ILNGTH.EQ.0) GO TO 160
C Calculate how many more (NREAD) direct access records must be
C read in to obtain all the integer information for the pivot block.
         J1 = ABS(ILNGTH) - ISPACE
         J2 = J1/LIBUFR
         IF (J1.GT.0) THEN
C Further records must be read from the integer direct access data set.
            NREAD = 1 + J2
            IF (J2*LIBUFR.EQ.J1) NREAD = NREAD - 1
            IF (IR2+NREAD*LIBUFR.GT.DIMIBF) THEN
C If there is not room to accommodate the incoming NREAD records
C in IBUFR, copy existing data in IBUFR to the beginning of
C IBUFR before calling MA42L/LD.
               DO 10 I = 1,ISPACE
                  IBUFR(I) = IBUFR(IR2-ISPACE+I)
   10          CONTINUE
C Reset IRECD, IR2.
               IRECD = 1
               IR2 = ISPACE
            END IF
            DO 20 IREAD = 1,NREAD
               CALL MA42LD(-1,IFIL,IKEY,BUFR,LBUFR,
     +                     IBUFR(IR2+ (IREAD-1)*LIBUFR+1),LIBUFR,LP,
     +                     JFLAG)
               IF (JFLAG.LT.0) GO TO 150
   20       CONTINUE
         END IF
C Required integers now in IBUFR.
C Update ISPACE
         IF (J1.LE.0) THEN
            ISPACE = ISPACE - ABS(ILNGTH)
         ELSE IF (J2*LIBUFR.EQ.J1) THEN
            ISPACE = 0
         ELSE
            ISPACE = LIBUFR - J1 + J2*LIBUFR
         END IF
C Update IR2
         IR2 = IR2 + NREAD*LIBUFR
C
         NPIV = IBUFR(IRECD+1)
         KFRNT = IBUFR(IRECD+2)
         LFRNT = IBUFR(IRECD+3)
         KFRNT1 = KFRNT - NPIV
         LFRNT1 = LFRNT - NPIV
         NPIV1 = (NPIV* (NPIV+1))/2
         IF (IND.EQ.1) THEN
            JFRNT = LFRNT
            JFRNT1 = LFRNT1
            K1 = IRECD + 3 + KFRNT
            LENGTH = NPIV* (LFRNT1+NRHSB) + NPIV1
         ELSE
            JFRNT = KFRNT
            JFRNT1 = KFRNT1
            K1 = IRECD + 3
            LENGTH = NPIV*KFRNT1 + NPIV1
         END IF
C Jump if the block is the zero block associated with singular
C matrix
         IF (ILNGTH.LT.0) GO TO 120
C
C Calculate how many more (NREAD) direct access records must be
C read in to obtain all the reals for the pivot block.
         NREAD = 0
         J1 = LENGTH - JSPACE
         J2 = J1/LBUFR
         IF (J1.GT.0) THEN
            NREAD = 1 + J2
            IF (J2*LBUFR.EQ.J1) NREAD = NREAD - 1
            IF (JR2+NREAD*LBUFR.GT.DIMBUF) THEN
C Copy to beginning of BUFR to give enough space to read in.
               DO 30 I = 1,JSPACE
                  BUFR(I) = BUFR(JR2-JSPACE+I)
   30          CONTINUE
               JRECD = 1
               JR2 = JSPACE
            END IF
            DO 40 IREAD = 1,NREAD
               CALL MA42LD(1,FIL,KEY,BUFR(JR2+ (IREAD-1)*LBUFR+1),LBUFR,
     +                     IBUFR,LIBUFR,LP,JFLAG)
               IF (JFLAG.LT.0) GO TO 150
   40       CONTINUE
         END IF
C
C Update JSPACE
         IF (J1.LE.0) THEN
            JSPACE = JSPACE - LENGTH
         ELSE IF (J2*LBUFR.EQ.J1) THEN
            JSPACE = 0
         ELSE
            JSPACE = LBUFR - J1 + J2*LBUFR
         END IF
C Update JR2
         JR2 = JR2 + NREAD*LBUFR
C
C Put right hand side vector into work space W.
         DO 60 L = 1,NRHS
            DO 50 K = 1,JFRNT
               K1PK = K1 + K
               KRO = IBUFR(K1PK)
               W(K,L) = R1(KRO,L)
   50       CONTINUE
   60    CONTINUE
C
C Modify right hand side vector.
C IP1 points to the first entry in the right hand side vector
C corresponding to pivot block.
         IP1 = JRECD + JFRNT1*NPIV
         DO 70 L = 1,NRHS
            IF (IND.EQ.1) CALL DTPSV('U','N','U',NPIV,BUFR(IP1),
     +                               W(JFRNT1+1,L),1)

            IF (IND.EQ.2) CALL DTPSV('U','N','N',NPIV,BUFR(IP1),
     +                               W(JFRNT1+1,L),1)
   70    CONTINUE
C
         IF (JFRNT1.NE.0) THEN
            IF (NRHS.NE.1) CALL DGEMM('N','N',JFRNT1,NRHS,NPIV,-ONE,
     +                                BUFR(JRECD),JFRNT1,W(JFRNT1+1,1),
     +                                LDW,ONE,W,LDW)
            IF (NRHS.EQ.1) CALL DGEMV('N',JFRNT1,NPIV,-ONE,BUFR(JRECD),
     +                                JFRNT1,W(JFRNT1+1,1),1,ONE,W,1)
         END IF
C Copy W into right hand side vector R1.
         DO 100 L = 1,NRHS
            DO 80 K = 1,JFRNT1
               K1PK = K1 + K
               KRO = IBUFR(K1PK)
               R1(KRO,L) = W(K,L)
   80       CONTINUE
            DO 90 K = JFRNT1 + 1,JFRNT
               K1PK = K1 + K
               KRO = IBUFR(K1PK)
               R1(KRO,L) = -W(K,L)
   90       CONTINUE
  100    CONTINUE
C
C Update JRECD and IRECD
         JRECD = JRECD + LENGTH
         IRECD = IRECD + ABS(ILNGTH)
C
  110 CONTINUE
  120 DO 140 L = 1,NRHS
         DO 130 K = 1,JFRNT
            K1PK = K1 + K
            KRO = IBUFR(K1PK)
            R1(KRO,L) = ZERO
  130    CONTINUE
  140 CONTINUE
      GO TO 160
C Set INFO if MA42L/LD returned an error.
  150 INFO = JFLAG
  160 RETURN
      END
C**********************************************************************
      SUBROUTINE MA42FD(AVAR,RHS,LRHS,NRHS,IVAR,LAST,NDF,NMAXE,NVAR,
     +                  MFRONT,NFRONT,BUFRL,LLB,BUFRU,LUB,IBUFR,LIBUFR,
     +                  FA,FRHS,LHED,KHED,KPIV,KPVLNK,LDEST,KDEST,ICNTL,
     +                  CNTL,IFILE,IREC,ISIZE,MKEY,NUMBLK,IELL,NELL,
     +                  KFRNT,LFRNT,ICOM,PIVBLK,BZERO,INFO,RINFO)
C
C    MA42F/FD (plus MA42D/DD) solve equations of the form
C                   A*X=B
C    where A is a sparse matrix derived from the element (elt) or
C    equation (eqn) assemblies and X and B are vectors. MA42F/FD
C    decomposes A into the product L*U of lower and upper triangular
C    matrices, and calculates (L**-1)*B=Y. (MA42D/DD calculates
C    X=(U**-1)*Y.) The contributions to matrix A from each elt/eqn
C    in turn are added in to the "frontal" array. Then as many
C    variables as possible are eliminated. A variable is referred to as
C    fully summed after the last elt/eqn which refers to it has been
C    assembled. Any entry corresponding to a fully summed row and
C    column is a candidate for pivot in the elimination, but if
C    possible only pivots which are greater than CNTL(2) times the
C    largest entry in their row are used to reduce roundoff.
C    Factors of L and U are output to buffers through subroutines
C    MA42G/GD and MA42H/HD.

C
C     Argument list. * indicates the argument is changed by the routine.
C
C  AVAR   - Real (DP) array of dimensions NMAXE by NVAR.
C           contains contributions to matrix from elt/eqn being input.
C  RHS    - Real (DP) array of dimensions NMAXE by LRHS.
C           Contains right hand sides for elt/eqn being input.
C  LRHS   - Integer variable. Second dimension of RHS array.
C  NRHS   - Integer variable.  Number of right hand sides.  NRHS<=LRHS.
C  IVAR   - Integer array of length NVAR. Contains indices of variables
C           in elt/eqn being input.  It is used locally as a work array
C           by MA42F/FD to facilitate access to the position of an
C           incoming variable in the front but it is reset and returned
C           unchanged to the user.
C  LAST   - Integer array of length NDF.  On input LAST(I) is the
C           elt/eqn in which variable I appears for the last
C           time.  Between calls to MA42F/FD (MA42B/BD), LAST is used
C           as a work array by MA42F/FD in order to reduce the work
C           required to find whether a variable in an incoming elt/eqn
C           is already in the front.  For a non-fully summed variable in
C           the front, JVAR say, -LAST(JVAR) is the position in the
C           destination vectors KDEST/LDEST of the information for
C           this variable.  On normal final exit or on detection of a
C           fatal error, LAST is reset and returned unchanged to the
C           user.
C  NDF    - Integer variable. Largest integer used to index a variable.
C  NMAXE  - Integer variable. Leading dimension of element or
C           equation arrays. Must be set to 1 for equation entry.
C  NVAR   - Integer variable. Number of variables in elt/eqn being
C           input.
C  MFRONT - Integer variable. Max. number of rows in frontal matrix.
C  NFRONT - Integer variable. Max. number of cols in frontal matrix.
C           For element entry NFRONT=MFRONT
C *BUFRL  - Real (DP) array of length LLB. In-core buffer for
C           entries of L.
C  LLB    - Integer variable. Length of array BUFRL.
C *BUFRU  - Real (DP) array of length LUB. In-core buffer for
C           entries of U.
C  LUB    - Integer variable. Length of array BUFRU.
C *IBUFR  _ Integer array of length LIBUFR. In-core buffer for
C           integer information on factors L and U.
C  LIBUFR - Integer variable. Length of array IBUFR.
C *FA     - Real (DP) array with dimensions MFRONT,NFRONT.
C           Holds the frontal matrix.
C *FRHS   - Real (DP) array with dimensions MFRONT, LRHS.
C           Holds the right hand sides corresponding to the current
C           frontal matrix.
C *LHED   - Integer array of dimension NFRONT. Used to hold the indices
C           of variables corresponding to columns in the front.
C *KHED   _ Integer array of dimension MFRONT. Used to hold the indices
C           of variables corresponding to rows in the front.
C *KPIV   _ Integer array of dimension NFRONT.  Indicates which
C           columns of FA are fully summed. (i.e. potential pivots).
C *KPVLNK - Integer array of dimension NFRONT. For each fully summed
C           column in the front, KPVLNK indicates its position in array
C           KPIV.  That is, it is the inverse array to KPIV.  It is
C           needed to enable KPIV to be updated easily after each
C           elimination.  For each non-fully summed column in the
C           front, KPVLNK holds temporarily the information
C           originally in LAST (i.e. the elt/eqn at which that
C           variable becomes fully summed).
C *LDEST  - Integer array of length NFRONT.  For each non-fully summed
C           variable in the front (global name JVAR say)
C           LDEST(-LAST(JVAR)) is its column position in the frontal
C           matrix.
C *KDEST  - Integer array of length MFRONT.  For each non-fully summed
C           variable in the front (global name JVAR say)
C           KDEST(-LAST(JVAR)) is its row position in the frontal
C           matrix.
C           It is not needed, used, or accessed during equation entries.
C   ICNTL - Integer array of length 8. See MA42B/BD.
C   CNTL  - Real (DP) array of length 2. See MA42B/BD.
C   IFILE - Integer array of length 3. Stream numbers for direct
C           access data sets.
C  *IREC  - Integer array of length 3. Pointers to first free space
C           in buffers.
C   ISIZE - Integer array of length 3. length of buffers.
C  *MKEY  - Integer array of length 3. Number of records written
C           to direct access data sets.
C   NUMBLK- Integer array of length 3. Number of records in
C           direct access data sets.
C   IELL  - Integer variable. Current elt/equ number.
C   NELL  - Integer variable. Total number of elt/equ.
C  *KFRNT - Integer variable. Current number of rows in front.
C  *LFRNT - Integer variable. Current number of cols in front.
C  *ICOM  - Integer array of length 5. Used to hold values which
C           must be preserved between calls to MA42F/FD.
C  PIVBLK - Holds min. pivot block size
C  BZERO  - Controls whether zeros in front exploited (they are
C           exploited if BZERO > 1)
C ICOM(1) - number of fully assembled variables.
C ICOM(2) - points to the last column searched for a pivot.
C ICOM(3) - points to free space.
C ICOM(4) - holds row index of best available pivot.
C ICOM(5) - holds column index of best available pivot.
C  *INFO    Integer array of length 23.
C INFO(1)  - Error flag. Negative values indicate a fatal error.
C INFO(2)  - Sign of determinant of matrix
C INFO(3)  - Number of variables.
C INFO(4)  - Number of non-zeros in U.
C INFO(5)  - Total real storage for U (plus right-hand sides).
C INFO(6)  - Number of non-zeros in L. (=Total real storage for L)
C INFO(7)  - Total integer storage.
C INFO(8)  - Maximum front-width (rows).
C INFO(9)  - Maximum front-width (columns).
C INFO(10) - Number of buffers for factors of U.
C INFO(11) - Number of buffers for factors of L.
C INFO(12) - Number of buffer for integers.
C INFO(13) - Number of columns examined during pivot searches.
C INFO(14) - Number of non-zeros tested for stability as pivots.
C INFO(15) - Number of non-zeros accessed during pivot selection
C            process.
C INFO(16) - Number of pivots chosen which did not satisfy threshold
C            criterion based on the value of CNTL(2).
C INFO(17) - Number of static condensations performed.
C INFO(18) - Potential number of static condensations (may be different
C            from INFO(17) because of stability considerations).
C INFO(19) - Maximum number of buffers required to hold block of pivot
C            rows.
C INFO(20) - Maximum number of buffers required to hold block of pivot
C            cols.
C INFO(21) - Max number of buffers required to hold block of integers.
C INFO(22) - Size of maximum pivot block.
C INFO(23) - Deficiency of matrix (singular case).
C *RINFO   - Real (DP) array of length 2. On final exit,
C            RINFO(1) holds natural log of determinant of the matrix
C            RINFO(2) holds number of operations in innermost loop
C
C  Local variables
C
C  DET    - Real (DP) variable. Equivalent to RINFO(1).
C  OPS    - Real (DP) variable. Equivalent to RINFO(2).
C  I      - Integer variable. Do  loop variable.
C  IFORCE - Integer variable. Number of forced eliminations for
C           current elt/equ.
C  ISTATC - Integer variable. Number of static condensations for
C           current elt/equ.
C  KPIVRX - Integer variable. Row index of best available pivot.
C           Equivalent to ICOM(6).
C  KR     - Integer variable. Number of fully assembled variables.
C           Equivalent to ICOM(1).
C  KS     - Integer variable. Used in cyclic search of pivot columns.
C           Equivalent to ICOM(2).
C  L      - Integer variable.
C  LFREE  - Integer variable. Holds position of first free entry
C           in KPVLNK. Equivalent to ICOM(5).
C  LFRNTT - Integer variable. Size of front after elt/equ has
C           been assembled.
C  LK     - Integer variable. Do loop variable. Ranges from 1 to NVAR.
C  LP     - Integer variable. Stream number for error messages.
C           (ICNTL(1)).
C  LPIVCX - Integer variable. Col. index of best available pivot.
C           Equivalent to ICOM(7).
C  MFR    - Integer variable.  Holds index of LK th variable in
C           current elt/equ (MFR = IVAR(LK)).
C  MP     - Integer variable. Stream number for warnings.
C  MVAR   - Integer variable. Number of rows in incoming elt/equ.
C  NELIM  - Integer variable. Number of fully assembled variables.
C  OFDIAG - Integer variable. On each retrun from MA42N/ND, OFDIAG
C           is the number of off-diagonal pivots selected on that
C           call.
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER IELL,KFRNT,LFRNT,LIBUFR,LLB,LRHS,LUB,MFRONT,NDF,NELL,
     +        NFRONT,NMAXE,NRHS,NVAR,PIVBLK,BZERO
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION AVAR(NMAXE,NVAR),BUFRL(LLB),BUFRU(LUB),CNTL(2),
     +                 FA(MFRONT,NFRONT),FRHS(MFRONT,LRHS),
     +                 RHS(NMAXE,LRHS),RINFO(2)
      INTEGER IBUFR(LIBUFR),ICNTL(8),ICOM(5),IFILE(3),INFO(23),IREC(3),
     +        ISIZE(3),IVAR(NVAR),KDEST(MFRONT),KHED(MFRONT),
     +        KPIV(NFRONT),KPVLNK(NFRONT),LAST(NDF),LDEST(NFRONT),
     +        LHED(NFRONT),MKEY(3),NUMBLK(3)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DET,OPS,PIVOT
      INTEGER I,IFORCE,ISTATC,KPIVRX,KR,KS,L,LFREE,LFRNTT,LK,LP,LPIVCX,
     +        MFR,MP,MVAR,NELIM,OFDIAG
      INTEGER NSIZE,NEW
      LOGICAL LSTAT
C Logical variable LSTAT detemines whether or not we wish to
C exploit any zeros in the front ... we will for equation entry.
C Seems also to be worthwhile for element entry, since there
C there can be a lot of zeros in front if the elements
C are poorly ordered. If elements well-ordered, looking for zeros can
C add an overhead to factorization time but may still give storage
C savings and savings in solve time.
C As a default we will always look for zeros in front ... in future,
C this could be a control parameter?
C     ..
C     .. External Subroutines ..
      EXTERNAL MA42GD,MA42HD,MA42MD,MA42ND,MA42OD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX
C     ..
      LP = ICNTL(1)
      MP = ICNTL(2)
      IF (IFILE(1).EQ.MP .OR. IFILE(2).EQ.MP .OR. IFILE(3).EQ.MP) MP = 0
      IF (IELL.EQ.1) THEN
C
C Initialization.  Only performed on first entry.
         ICOM(1) = 0
         ICOM(2) = 0
         ICOM(3) = 1
         ICOM(4) = 0
         ICOM(5) = 0
C Initialise LDEST.
         DO 10 I = 1,NFRONT
            LDEST(I) = I + 1
   10    CONTINUE
      END IF
C End of initialization phase.

         LSTAT = .FALSE.
         IF (BZERO.GT.1) LSTAT = .TRUE.
c         IF (NMAXE.EQ.1) LSTAT = .TRUE.
c         if (lstat) write (6,*) ' Exploiting zeros in front'

C only perform eliminations if we have at least PIVBLK fully summed
C variables ... (HSL12 version equivalent
C to PIVBLK = 1)

      KR = ICOM(1)
      KS = ICOM(2)
      LFREE = ICOM(3)
      KPIVRX = ICOM(4)
      LPIVCX = ICOM(5)
      DET = RINFO(1)
      OPS = RINFO(2)
      OFDIAG = 0
C
C MVAR is number of rows in incoming elt/eqn.  It is equal 1 for
C an equation entry and is equal to NVAR (number of columns) for
C an element entry.
      MVAR = NVAR
      IF (NMAXE.EQ.1) MVAR = 1
C
C Process columns of input elt/eqn.
C First check to see if incoming elt/eqn will fit in front.
C IFORCE is number of eliminations that must be performed in order
C to create sufficient space in front matrix to assemble incoming
C element or equation.
      IFORCE = 0
C NEW is the number of variables coming into the front for
C the first time
      NEW = 0
C ISTATC is number of static condensations which may be performed on
C incoming elt/eqn.
      ISTATC = 0
      DO 20 LK = 1,NVAR
         MFR = IVAR(LK)
C Check on validity of IVAR indices.
         IF (MFR.GT.NDF .OR. MFR.LE.0) GO TO 50
         IF (LAST(MFR).GE.0) THEN
C Variable not yet in front
C Jump to error return if variable is already fully summed.
            IF (LAST(MFR).LT.IELL) GO TO 70
            NEW = NEW + 1
C Test for static condensation.
            IF (LAST(MFR).EQ.IELL) ISTATC = ISTATC + 1
         END IF
   20 CONTINUE
C
C  If equ entry and ISTATC>1 matrix is structurally singular.
      IF (NMAXE.EQ.1 .AND. ISTATC.GT.1) THEN
         INFO(2) = 0
         DET = ZERO
C Jump to error return if ICNTL(8) is 0.
         IF (ICNTL(8).EQ.0) GO TO 80
C Otherwise set INFO(1) = 1 and issue a warning.
C (Do not overwrite INFO(1) = 2, 5 or 6).
         IF (INFO(1).EQ.0) THEN
            INFO(1) = 1
            IF (MP.GT.0) THEN
               WRITE (MP,FMT=9060) INFO(1),IELL
               WRITE (MP,FMT=9040)
            END IF
         END IF
C Do not perform static condensations.
         ISTATC = 0
      END IF
C
C LFRNTT is number of columns in the front after elt/equ
C  has been assembled (for element entry, LFRNTT
C is also number of rows in front after assembly).
      LFRNTT = LFRNT + NEW
C  INFO(3) holds the number of variables in the problem.
      INFO(3) = INFO(3) + NEW
      IFORCE = MAX(0,LFRNTT-NFRONT)
C For equation entry, check incoming equation (row) can be accommodated.
      IF (NMAXE.EQ.1 .AND. KFRNT.EQ.MFRONT) IFORCE = MAX(IFORCE,1)
C
C Jump to error return if we cannot eliminate sufficient variables
C to leave enough room in front for the incoming elt/eqn.
      IF (IFORCE.GT.KR) GO TO 60
      IF (IFORCE.GT.0 .AND. PIVBLK.GT.1) THEN
C
C There is insufficient room for the incoming elt/eqn in
C the front, but we may still be able to perform some more eliminations
C (because we were waiting for a minimum pivot block to be available).
C These are NOT forced eliminations
         NELIM = KR
         NSIZE = NELIM
         CALL MA42ND(IELL,NDF,LAST,MFRONT,NFRONT,FA,LHED,KHED,LRHS,NRHS,
     +               FRHS,KDEST,LDEST,KPIV,KPVLNK,KPIVRX,LPIVCX,PIVOT,
     +               KFRNT,LFRNT,KR,KS,NMAXE,NELIM,0,LIBUFR,IBUFR,
     +               BUFRL,LLB,BUFRU,LUB,IFILE,IREC,ISIZE,MKEY,NUMBLK,
     +               DET,OPS,NELL,CNTL,ICNTL,INFO,OFDIAG,NSIZE,LSTAT)
         IF (INFO(1).LT.0) GO TO 90
C Update the frontsize after elt/equ has been assembled.
         LFRNTT = LFRNT + NEW
         IFORCE = MAX(0,LFRNTT-NFRONT)
         IF (NMAXE.EQ.1 .AND. KFRNT.EQ.MFRONT) IFORCE = MAX(IFORCE,1)
      END IF

C If there is still insufficient room for the incoming elt/eqn in
C the front and forced eliminations would yield enough room
C then we first do these forced eliminations.
C
      IF (IFORCE.GT.0) THEN
C
         NELIM = IFORCE
         NSIZE = NELIM
         CALL MA42ND(IELL,NDF,LAST,MFRONT,NFRONT,FA,LHED,KHED,LRHS,NRHS,
     +               FRHS,KDEST,LDEST,KPIV,KPVLNK,KPIVRX,LPIVCX,PIVOT,
     +               KFRNT,LFRNT,KR,KS,NMAXE,NELIM,IFORCE,LIBUFR,IBUFR,
     +               BUFRL,LLB,BUFRU,LUB,IFILE,IREC,ISIZE,MKEY,NUMBLK,
     +               DET,OPS,NELL,CNTL,ICNTL,INFO,OFDIAG,NSIZE,LSTAT)
         IF (INFO(1).LT.0) GO TO 90
         IF (INFO(1).EQ.4) GO TO 60
      END IF

C Now perform static condensations
      IF (ICNTL(7).NE.0 .AND. ICNTL(7).NE.1001) ISTATC = 0
      CALL MA42OD(IELL,AVAR,LRHS,NMAXE,RHS,NRHS,IVAR,NDF,LAST,NVAR,
     +            MFRONT,NFRONT,LHED,KHED,KR,KPIV,KPVLNK,LDEST,KDEST,
     +            ISTATC,LFREE,FRHS,BUFRL,LLB,BUFRU,LUB,IBUFR,LIBUFR,
     +            MVAR,KFRNT,LFRNT,FA,IFILE,IREC,ISIZE,MKEY,NUMBLK,DET,
     +            OPS,NELL,CNTL,ICNTL,INFO)
      IF (INFO(1).LT.0) GO TO 90
C
C  Add in elt/eqn information which was waiting for sufficient
C  space to be made by forced eliminations.
      CALL MA42MD(IELL,AVAR,LRHS,NMAXE,RHS,NRHS,IVAR,LAST,NDF,NVAR,
     +            MFRONT,NFRONT,FA,FRHS,LHED,KHED,KPIV,KPVLNK,LDEST,
     +            KDEST,ISTATC,LFREE,KR,KFRNT,LFRNT,INFO)
C
C Now try and eliminate the fully assembled variables.
      IF (KR.EQ.0) GO TO 40
      IF (KR.LT.PIVBLK .AND. IELL.NE.NELL) GO TO 40
C ICOM(1) was the number of fully summed variables on entry
C to MA42F/FD. We don't want to do any searches for pivots if
C there are no new fully summed variables.
      IF (IFORCE.EQ.0 .AND. KR.EQ.ICOM(1) .AND. IELL.NE.NELL) GO TO 40
      NELIM = KR
      IFORCE = 0
C We are adding an extra parameter NSIZE
C The standard code uses NSIZE = KR, and this attempts to
C eliminate as many of the fully summed variables at once as possible.
C If we want to just eliminate one variable at once, set NSIZE = 1.
C Using NSIZE = 1 with PIVBLK = 1
C will ensure no zeros in factors BUT will
C mean only BLAS 2 is used and the integer storage will
C increase (it will be slightly more than the real storage for
C the factors)
      NSIZE = KR
C     NSIZE  = 1
      DO 30 I = 1,KR - NSIZE + 1
         NELIM = NSIZE
         CALL MA42ND(IELL,NDF,LAST,MFRONT,NFRONT,FA,LHED,KHED,LRHS,NRHS,
     +               FRHS,KDEST,LDEST,KPIV,KPVLNK,KPIVRX,LPIVCX,PIVOT,
     +               KFRNT,LFRNT,KR,KS,NMAXE,NELIM,IFORCE,LIBUFR,IBUFR,
     +               BUFRL,LLB,BUFRU,LUB,IFILE,IREC,ISIZE,MKEY,NUMBLK,
     +               DET,OPS,NELL,CNTL,ICNTL,INFO,OFDIAG,NSIZE,LSTAT)
         IF (INFO(1).LT.0) GO TO 90
   30 CONTINUE
   40 CONTINUE
C
      IF (IELL.EQ.NELL) THEN
C Forward elimination complete.
         IF (KR.GT.0) THEN
C Singular matrix. Write to integer buffer.
            CALL MA42HD(-1,KFRNT,LFRNT,KR,IBUFR,LIBUFR,KHED,LHED,MFRONT,
     +                  NFRONT,LP,IFILE(3),IREC(3),ISIZE(3),MKEY(3),
     +                  NUMBLK(3),IELL,NELL,INFO)
            INFO(23) = KR
         END IF
C Output what remains in the buffers to direct access data sets.
         CALL MA42HD(2,1,1,0,IBUFR,LIBUFR,KHED,LHED,MFRONT,NFRONT,LP,
     +               IFILE(3),IREC(3),ISIZE(3),MKEY(3),NUMBLK(3),IELL,
     +               NELL,INFO)
         IF (INFO(1).LT.0) GO TO 90
         CALL MA42GD(1,2,FA,FRHS,NRHS,1,1,1,BUFRU,LUB,MFRONT,
     +               LP,IFILE(1),IREC(1),ISIZE(1),MKEY(1),
     +               NUMBLK(1),IELL,NELL,INFO)
         IF (INFO(1).LT.0) GO TO 90
         CALL MA42GD(2,2,FA,FRHS,NRHS,1,1,1,BUFRL,LLB,MFRONT,
     +               LP,IFILE(2),IREC(2),ISIZE(2),MKEY(2),
     +               NUMBLK(2),IELL,NELL,INFO)
         IF (INFO(1).LT.0) GO TO 90
C Update information for return to user.
         IF (INFO(16).GT.0) THEN
            IF (INFO(1).EQ.0 .OR. INFO(1).EQ.1) INFO(1) = 2
            IF (MP.GT.0) WRITE (MP,FMT=9000) INFO(16)
         END IF
C Transfer information on number of buffers used.
         INFO(10) = MAX(1,MKEY(1))
         INFO(11) = MKEY(2)
         IF (IFILE(2).NE.0) INFO(11) = MAX(1,MKEY(2))
         INFO(12) = MAX(1,MKEY(3))
C If diagonal pivoting was used (abs(ICNTL(7))=1001) then
C return to the user the number of off-diagonal pivots that
C were selected. Off-digaonal pivots can only be selected if
C IELL = NELL
         IF (INFO(1).EQ.0 .AND. ABS(ICNTL(7)).EQ.1001) THEN
            IF (OFDIAG.GT.0) THEN
               INFO(1) = 6 + OFDIAG
               IF (MP.GT.0) WRITE (MP,FMT=9070) OFDIAG
            END IF
         END IF
      END IF
C Obtain next element.
      GO TO 110
C
C **** Error returns ****
   50 INFO(1) = -3
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9050) INFO(1),IELL
         WRITE (LP,FMT=9010) LK,MFR
      END IF
      GO TO 90
   60 INFO(1) = 4
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9050) INFO(1),IELL
         WRITE (LP,FMT=9020) MFRONT,NFRONT
      END IF
      GO TO 90
   70 INFO(1) = -13
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9050) INFO(1),IELL
         WRITE (LP,FMT=9030) MFR
      END IF
      GO TO 90
   80 INFO(1) = -14
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9050) INFO(1),IELL
         WRITE (LP,FMT=9040)
      END IF
C Reset LAST on error return.
   90 DO 100 L = 1,LFRNT
         IF (KPVLNK(L).GT.0) GO TO 100
         MFR = LHED(L)
         LAST(MFR) = -KPVLNK(L)
  100 CONTINUE
C Store ICOM and update RINFO(1), RINFO(2)
  110 ICOM(1) = KR
      ICOM(2) = KS
      ICOM(3) = LFREE
      ICOM(4) = KPIVRX
      ICOM(5) = LPIVCX
      RINFO(1) = DET
      RINFO(2) = OPS
      RETURN
 9000 FORMAT (' ***** Warning from MA42B/BD *****  INFO(1) = 2',
     +       /7X,'Numerical criterion not satisfied by ',I8,' pivots')
 9010 FORMAT (7X,'Variable',I8,' in elt/eqn has value ',I8)
 9020 FORMAT (7X,'NFRONT not large enough, currently equals ',2I8)
 9030 FORMAT (7X,'Variable',I8,' is already fully summed')
 9040 FORMAT (7X,'Matrix found to be singular')
 9050 FORMAT (' ***** Error return from MA42B/BD *****  INFO(1) = ',I3,
     +       /' ***** Error discovered after input of elt/eqn *****',I8)
 9060 FORMAT (' ***** Warning from MA42B/BD *****  INFO(1) = ',I3,
     +       /7X,'after input of elt/equ ',I8)
 9070 FORMAT (' ***** Warning from MA42B/BD *****  ',
     +       /7X,I8,'  off-diagonal pivots were selected.')
      END
C**********************************************************************
      SUBROUTINE MA42GD(IND,JFL,FA,FRHS,NRHS,KFRNT,LFRNT,NPIV,BUFR,
     +                 LBUFR,MFRONT,LP,IFILE,IREC,ISIZE,MKEY,
     +                 NUMBLK,IELL,NELL,INFO)
C
C     MA42G/GD puts the rows of U  or columns of L into the in core
C     buffers. When a buffer is full, the contents are output to
C     direct access data sets. A separate direct access data set
C     holds integer information.
C
C     Note: it is necessary to call MA42G/GD with IND = 1
C           and the MA42G/GD with IND = 2 (in that order) so that
C           FA and FRHS are correctly reset to zero at the and
C           of the call
C
C     Argument list. * indicates the argument is changed by the routine.
C
C  IND   - Integer variable. Used to indicate whether we are writing
C          rows of U or columns of L to main storage buffers.
C          IND = 1 rows of U
C          IND = 2 columns of L
C  JFL   - Integer variable. Indicates nature of call.
C          JFL = 1 except on final call when JFL = 2.
C *FA    - Real (DP) array with dimensions MFRONT,NFRONT.
C          Holds the frontal matrix.
C *FRHS  - Real (DP) array with dimensions MFRONT, LRHS.
C          Holds the right hand sides corresponding to the current
C          frontal matrix.
C  NRHS  - Integer variable.  Number of right hand sides.  NRHS<=LRHS.
C  KFRNT - Integer variable. Number of rows in front.
C  LFRNT - Integer variable. Number of columns in front.
C  NPIV  - Integer variable. Number of rows or columns to be written
C          to storage buffers.
C  BUFR  - Real (DP) array of length LBUFR. The output buffer.
C  LBUFR - Integer variable. Length of array BUFR.
C  MFRONT- Integer variable. Max. number of rows in frontal matrix.
C  NFRONT- Integer variable. Max. number of cols in frontal matrix.
C  LP    - Integer variable. Holds stream number for error messages.
C  IFILE - Integer variable. Holds stream number for direct
C          access data set.
C  ISIZE - Integer variable. Holds length of buffer.
C *MKEY  - Integer variable. Holds number of records written
C          to the direct access data set.
C  NUMBLK- Integer variable. Holds number of records in direct access
C          data set.
C *IREC  - Integer variable.  Points to first available space in
C          buffer.
C  IELL  - Integer variable. Holds current elt/equ number.
C *INFO  - Integer array of length 23.
C          INFO(1) is an error flag. Possible nonzero values
C          on return are -26, 5, 6.
C          INFO(1) = -26 if there is a failure in a direct access write.
C          INFO(1) = 5 if buffer not long enough and no direct access
C                   data set  requested.
C          INFO(1) = 6 if insufficient space allocated to direct access
C                   data set.
C  INFO(4), INFO(5), and INFO(6) are also altered by this routine.
C  INFO(4) holds the number of nonzeros in U.
C  INFO(5) holds the total storage for U together with
C          right-hand sides.
C  INFO(6) holds the number of nonzeros in L (which is also
C          the total storage for L).
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER IELL,IFILE,IND,IREC,ISIZE,JFL,KFRNT,LBUFR,LFRNT,LP,
     +        MFRONT,MKEY,NELL,NPIV,NRHS,NUMBLK
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION BUFR(LBUFR),FA(MFRONT,*),FRHS(MFRONT,*)
      INTEGER INFO(23)
C     ..
C     .. Local Arrays ..
      INTEGER IBUFR(1)
C     ..
C     .. Local Scalars ..
      INTEGER I,I1,IBLOCK,IFRNT1,ISPACE,IUP,J,J1,J2,JBUFR,JFLAG,JFRNT1,
     +        K,K1,KFRNT1,LENGTH,LFRNT1,LIBUFR,NBUFR,NPIV1
C     ..
C     .. External Subroutines ..
      EXTERNAL MA42KD,MA42LD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MIN
C We introduced these variables to allow us to count
C zeros in front.
C      INTEGER NZEROL,NZEROU
C      SAVE NZEROL,NZEROU

C      IF (INFO(6).EQ.0) NZEROL = 0
C      IF (INFO(4).EQ.0) NZEROU = 0
C     ..
      IF (JFL.EQ.1) THEN
C Set local variables KFRNT1, LFRNT1, NPIV1, IFRNT1, JFRNT1, and
C LENGTH and update statistics (INFO(4) and INFO(5) or INFO(6)).
         KFRNT1 = KFRNT - NPIV
         LFRNT1 = LFRNT - NPIV
         NPIV1 = (NPIV* (NPIV+1))/2
         IF (IND.EQ.1) THEN
            IFRNT1 = KFRNT1
            JFRNT1 = LFRNT1
            LENGTH = NPIV* (LFRNT1+NRHS) + NPIV1
            INFO(4) = INFO(4) + NPIV*LFRNT1 + NPIV1
            INFO(5) = INFO(5) + LENGTH
         ELSE
            IFRNT1 = LFRNT1
            JFRNT1 = KFRNT1
            LENGTH = NPIV*KFRNT1 + NPIV1
            INFO(6) = INFO(6) + LENGTH
         END IF
      ELSE
C Last call. Immediate return if nothing left in buffer.
         IF (IREC.EQ.1) THEN
C            IF (IND.EQ.1) WRITE (6,FMT=*)
C     +          'Number of zeros in U factor ',NZEROU
C            IF (IND.EQ.2) WRITE (6,FMT=*)
C     +          'Number of zeros in L factor ',NZEROL
C            IF (IND.EQ.2) WRITE (6,FMT=*)
C     +          'Total Number of zeros in factors ',NZEROU + NZEROL
            GO TO 110
         END IF
      END IF
C If buffer has zero length (which can only be the case
C for the L-buffer) zero entries and return.
      IF (ISIZE.EQ.0) THEN
         DO 20 I = 1,NPIV
            I1 = IFRNT1 + I
            DO 10 J = 1,JFRNT1 + I
               FA(J,I1) = ZERO
   10       CONTINUE
   20    CONTINUE
         GO TO 110
      END IF
C
C MA42K/KD  checks to see if a write to direct access is required
C and, if so, checks that direct access has been requested and
C that the direct access data set has sufficient records.
C MA42L/LD performs write to direct access data set.
C
      CALL MA42KD(JFL,IND,LENGTH,IFILE,ISIZE,MKEY,NUMBLK,IREC,IELL,NELL,
     +           LP,INFO)
      IF (INFO(1).EQ.5 .OR. INFO(1).EQ.6) GO TO 105
C
C Set LIBUFR to 1 for the call to MA42L/LD
      LIBUFR = 1
      IF (JFL.EQ.1) THEN
C Write out pivot rows/columns and right hand sides to buffer.
C Zero entries once they have been output.
         DO 80 IBLOCK = 1,3
C  First time through this loop output rectangular block,
C  second time through output triangular part, and third
C  time through output right-hand sides (if any).
            IF (IBLOCK.EQ.3 .AND. IND.EQ.2) GO TO 80
            IUP = NPIV
            IF (IBLOCK.EQ.3) IUP = NRHS
            DO 70 I = 1,IUP
               I1 = IFRNT1 + I
C ISPACE is space in buffer.
               ISPACE = ISIZE - IREC + 1
C J2 is the number of entries to be output.
               IF (IBLOCK.EQ.1) THEN
                  J = 0
                  J2 = JFRNT1
               ELSE IF (IBLOCK.EQ.2) THEN
                  J = JFRNT1
                  J2 = I
               ELSE
                  J = IFRNT1
                  J2 = NPIV
               END IF
               J1 = J2 - ISPACE
C Compute how many buffers (NBUFR) will be required to output
C the J2 entries.
C If J1 is at most 0, the entries can
C be accommodated in the present buffer and no write to
C the corresponding direct access file is needed.
               NBUFR = 1
               IF (J1.GT.0) NBUFR = 2 + (J1/ISIZE)
C Loop over the number of buffers required.
               DO 60 JBUFR = 1,NBUFR
                  K1 = ISIZE
                  IF (JBUFR.EQ.NBUFR) K1 = J1 - (J1/ISIZE)*ISIZE
                  IF (JBUFR.EQ.1) K1 = MIN(J2,ISPACE)
                  IF (IND.EQ.1) THEN
                     IF (IBLOCK.EQ.1 .OR. IBLOCK.EQ.2) THEN
C Output rows of U to buffer.
                        DO 30 K = 1,K1
C                           IF (FA(I1,J+K).EQ.ZERO) NZEROU = NZEROU + 1
                           BUFR(IREC+K-1) = FA(I1,J+K)
C Do not zero diagonal entries of triangular part of U
                           IF (J+K.NE.LFRNT1+I) FA(I1,J+K) = ZERO
   30                   CONTINUE
                     ELSE
C Output right-hand sides to buffer.
                        DO 40 K = 1,K1
                           BUFR(IREC+K-1) = FRHS(J+K,I)
                           FRHS(J+K,I) = ZERO
   40                   CONTINUE
                     END IF
                  ELSE
C Output columns of L to buffer.
                     DO 50 K = 1,K1
C                        IF (FA(J+K,I1).EQ.ZERO) NZEROL = NZEROL + 1
                        BUFR(IREC+K-1) = FA(J+K,I1)
                        FA(J+K,I1) = ZERO
   50                CONTINUE
                  END IF
                  J = J + K1
                  IF (JBUFR.EQ.NBUFR .AND. IREC+K1.LE.ISIZE) THEN
                     IREC = IREC + K1
                  ELSE
C Write buffer to direct access data set and reset IREC.
                     MKEY = MKEY + 1
                     CALL MA42LD(3,IFILE,MKEY,BUFR,LBUFR,IBUFR,LIBUFR,
     +                          LP,JFLAG)
                     IF (JFLAG.LT.0) GO TO 100
                     IREC = 1
                  END IF
   60          CONTINUE
   70       CONTINUE
   80    CONTINUE

      ELSE
C         IF (IND.EQ.1) WRITE (6,FMT=*) 'Number of zeros in U factor ',
C     +       NZEROU
C         IF (IND.EQ.2) WRITE (6,FMT=*) 'Number of zeros in L factor ',
C     +       NZEROL
C         IF (IND.EQ.2) WRITE (6,FMT=*)
C     +       'Total Number of zeros in factors ',NZEROU + NZEROL
C Final call to MA42G/GD. Unless IREC=1,
C BUFR still contains reals to be output to direct access data set.
         IF (IREC.NE.1) THEN
            MKEY = MKEY + 1
C If we are using d.a. files and we have not yet written the
C buffer to the d.a. file, fill-in the rest of BUFR with garbage
C (so that it is not undefined)
            IF (IFILE.NE.0 .AND. MKEY.EQ.1) THEN
               DO 90 I = IREC,LBUFR
                  BUFR(I) = ZERO
   90          CONTINUE
            END IF
            IF (IFILE.NE.0) THEN
               CALL MA42LD(3,IFILE,MKEY,BUFR,LBUFR,IBUFR,LIBUFR,
     +                     LP,JFLAG)
               IF (JFLAG.LT.0) GO TO 100
            END IF
         END IF
      END IF
      GO TO 110

C **** Fatal error return ****
  100 INFO(1) = JFLAG
      IF (LP.GT.0) WRITE (LP,FMT=9000) INFO(1)
      GO TO 110

  105 CONTINUE
C Zero entries in FA as if they have been output.
      IF (JFL.EQ.1) THEN
C Write out pivot rows/columns and right hand sides to buffer.
C Zero entries once they have been output.
         DO 85 IBLOCK = 1,3
C  First time through this loop output rectangular block,
C  second time through output triangular part, and third
C  time through output right-hand sides (if any).
            IF (IBLOCK.EQ.3 .AND. IND.EQ.2) GO TO 85
            IUP = NPIV
            IF (IBLOCK.EQ.3) IUP = NRHS
            DO 75 I = 1,IUP
               I1 = IFRNT1 + I
C ISPACE is space in buffer.
               ISPACE = ISIZE - IREC + 1
C J2 is the number of entries to be output.
               IF (IBLOCK.EQ.1) THEN
                  J = 0
                  J2 = JFRNT1
               ELSE IF (IBLOCK.EQ.2) THEN
                  J = JFRNT1
                  J2 = I
               ELSE
                  J = IFRNT1
                  J2 = NPIV
               END IF
               J1 = J2 - ISPACE
C Compute how many buffers (NBUFR) will be required to output
C the J2 entries.
C If J1 is at most 0, the entries can
C be accommodated in the present buffer and no write to
C the corresponding direct access file is needed.
               NBUFR = 1
               IF (J1.GT.0) NBUFR = 2 + (J1/ISIZE)
C Loop over the number of buffers required.
               DO 65 JBUFR = 1,NBUFR
                  K1 = ISIZE
                  IF (JBUFR.EQ.NBUFR) K1 = J1 - (J1/ISIZE)*ISIZE
                  IF (JBUFR.EQ.1) K1 = MIN(J2,ISPACE)
                  IF (IND.EQ.1) THEN
                     IF (IBLOCK.EQ.1 .OR. IBLOCK.EQ.2) THEN
                        DO 35 K = 1,K1
                           IF (J+K.NE.LFRNT1+I) FA(I1,J+K) = ZERO
   35                   CONTINUE
                     ELSE
                        DO 45 K = 1,K1
                           FRHS(J+K,I) = ZERO
   45                   CONTINUE
                     END IF
                  ELSE
                     DO 55 K = 1,K1
                        FA(J+K,I1) = ZERO
   55                CONTINUE
                  END IF
                  J = J + K1
                  IF (JBUFR.EQ.NBUFR .AND. IREC+K1.LE.ISIZE) THEN
                     IREC = IREC + K1
                  ELSE
C Update MKEY and IREC.
                     MKEY = MKEY + 1
                     IREC = 1
                  END IF
   65          CONTINUE
   75       CONTINUE
   85    CONTINUE

      ELSE
C Final call to MA42G/GD. Unless IREC=1,
C BUFR still contains reals to be output to direct access data set.
         IF (IREC.NE.1) MKEY = MKEY + 1
      END IF

  110 RETURN
 9000 FORMAT (' ***** Error return from MA42B/BD *****  INFO(1) = ',I3)
      END
C*********************************************************************
      SUBROUTINE MA42HD(JFL,KFRNT,LFRNT,NPIV,IBUFR,LIBUFR,KHED,LHED,
     +                 MFRONT,NFRONT,LP,IFILE,IREC,ISIZE,MKEY,NUMBLK,
     +                 IELL,NELL,INFO)
C
C     MA42H/HD writes the integer information on the rows of U and
C     columns of L to the main (integer) storage buffer.
C     If there is insufficient space in the buffer for the
C     incoming row/col information then the contents of the buffer
C     is output to a direct access data set which is needed because
C     MA42D/DD reads the records in the opposite order from that in
C     which MA42F/FD produces them.
C
C     The integer information which is output is, in order:
C     - The length of the record.
C     - The number of rows and columns in the pivot block (NPIV).
C     - The number of rows in the front (KFRNT).
C     - The number of columns in the front (LFRNT).
C     - The local row indices of the variables in the front (KHED).
C     - The local column indices of the variables in the front (LHED).
C     - The length of the record. (included so records can be read in
C       reverse order).
C
C     Argument list. * indicates the argument is changed by the routine.
C                                                       *
C  JFL   - Integer variable. Indicates nature of call.
C          JFL = -1 when called after singularity detected
C          JFL = 2 on final call
C          otherwise JFL = 1
C  KFRNT - Integer variable. Number of rows in front.
C  LFRNT - Integer variable. Number of columns in front.
C  NPIV  - Integer variable. Number of rows or columns to be written
C          to storage buffers.
C  IBUFR - Integer array of length LIBUFR. The output buffer.
C  ILBUFR- Integer variable. Length of array IBUFR.
C  LHED  - Integer array of dimension NFRONT. Holds the
C          indices of variables corresponding to columns in the
C          front.
C  *KHED - Integer array of dimension MFRONT. Holds the
C          indices of variables corresponding to rows in the
C          front.
C  MFRONT- Integer variable. Max. number of rows in frontal matrix.
C  NFRONT- Integer variable. Max. number of cols in frontal matrix.
C  LP    - Integer variable. Holds stream number for error messages.
C  IFILE - Integer variable. Holds stream number for direct
C          access data set.
C  ISIZE - Integer variable. Holds length of buffer.
C  MKEY  - Integer variable. Holds number of records written
C          to the direct access data set.
C  NUMBLK- Integer variable. Holds number of records in direct access
C          data set.
C *IREC  - Integer variable.  Points to first available space in
C          buffer.
C  IELL  - Integer variable. Holds current elt/equ number.
C *INFO  - Integer array of length 23.
C          INFO(1) is an error flag. Possible nonzero values
C          for INFO(1) on return are -26, 5, 6.
C          INFO(1) = -26 if there is a failure in a direct access write.
C          INFO(1) = 5 if buffer not long enough and no direct access
C                   data set requested.
C          INFO(1) = 6 if insufficient space allocated to direct access
C                   data set.
C  INFO(3) and INFO(7) are also altered by this routine.
C 10.12.1998 edited code so that INFO(3) is NOT updated here
C but is updated in MA42F/FD
C  INFO(3) holds the number of variables in the problem.
C  INFO(7) holds the total amount of storage for the row and
C          column indices).
C  INFO(22) holds largest size of pivot block used.
C
C     .. Scalar Arguments ..
      INTEGER IELL,IFILE,IREC,ISIZE,JFL,KFRNT,LFRNT,LIBUFR,LP,MFRONT,
     +        MKEY,NELL,NFRONT,NPIV,NUMBLK
C     ..
C     .. Array Arguments ..
      INTEGER IBUFR(LIBUFR),INFO(23),KHED(MFRONT),LHED(NFRONT)
C     ..
C     .. Local Scalars ..
      INTEGER I,ISPACE,J,J1,J2,JBUFR,JFLAG,K,K1,LBUFR,LENGTH,NBUFR,
     +        JFIT,JFIT1
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION BUFR(1)
      INTEGER ITEMP(4)
C     ..
C     .. External Subroutines ..
      EXTERNAL MA42KD,MA42LD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,MIN
C     ..
C Set local variable LENGTH to the length of the incoming record.
      IF (JFL.EQ.1) THEN
         LENGTH = 5 + LFRNT + KFRNT
      ELSE IF (JFL.EQ.-1) THEN
C Use a negative flag to indicate block follows detection of
C singularity.
         LENGTH = - (5+LFRNT+KFRNT)
      ELSE
         LENGTH = 1
      END IF
C Update INFO(3) (the total number of variables in the problem).
C Update INFO(7) (the total amount of storage for the row and
C column indices).
C Also update INFO(22) (max size of pivot block used).
      INFO(7)  = INFO(7) + ABS(LENGTH)
      INFO(22) = MAX(INFO(22),NPIV)
C
C MA42K/KD  checks to see if a write to direct access is required
C and, if so, checks that direct access has been requested and
C that the direct access data set has sufficient records.
C MA42L/LD performs write to direct access data set.
C
      CALL MA42KD(JFL,3,ABS(LENGTH),IFILE,ISIZE,MKEY,NUMBLK,IREC,IELL,
     +           NELL,LP,INFO)
      IF (INFO(1).EQ.5 .OR. INFO(1).EQ.6) GO TO 85
C
C Set LBUFR to 1 for the call to MA42L/LD
      LBUFR = 1
      IF (ABS(JFL).EQ.1) THEN
         ITEMP(1) = LENGTH
         ITEMP(2) = NPIV
         ITEMP(3) = KFRNT
         ITEMP(4) = LFRNT
         DO 60 I = 1,4
C On the first pass through the loop ITEMP is output. On the
C second pass the row indices are output and, on the third,
C the column indices are output. On last pass LENGTH is output.
C ISPACE is the space available in the buffer.
            ISPACE = ISIZE - IREC + 1
C J2 is the number of integers to be written to the buffer.
            IF (I.EQ.1) J2 = 4
            IF (I.EQ.2) J2 = KFRNT
            IF (I.EQ.3) J2 = LFRNT
            IF (I.EQ.4) J2 = 1
            J1 = J2 - ISPACE
C Calculate the number (NBUFR) of buffers required to accommodate
C the J2 integers.
C If J1 is at most 0, the entries can
C be accommodated in the present buffer and no write to
C the corresponding direct access file is needed.
            NBUFR = 1
            IF (J1.GT.0) NBUFR = 2 + (J1/ISIZE)
C Initialise J.
            J = 0
            DO 50 JBUFR = 1,NBUFR
               K1 = ISIZE
               IF (JBUFR.EQ.NBUFR) K1 = J1 - (J1/ISIZE)*ISIZE
               IF (JBUFR.EQ.1) K1 = MIN(J2,ISPACE)
               IF (I.EQ.1) THEN
                  DO 10 K = 1,K1
                     IBUFR(IREC+K-1) = ITEMP(J+K)
   10             CONTINUE
               ELSE IF (I.EQ.2) THEN
C Output the row indices.
                  DO 20 K = 1,K1
                     IBUFR(IREC+K-1) = KHED(J+K)
   20             CONTINUE
               ELSE IF (I.EQ.3) THEN
C Output the column indices.
                  DO 30 K = 1,K1
                     IBUFR(IREC+K-1) = LHED(J+K)
   30             CONTINUE
               ELSE
C Output length (LENGTH) of integer record.
                  DO 40 K = 1,K1
                     IBUFR(IREC) = LENGTH
   40             CONTINUE
               END IF
C Update J.
               J = J + K1
               IF (JBUFR.EQ.NBUFR .AND. IREC+K1.LE.ISIZE) THEN
                  IREC = IREC + K1
               ELSE IF (IFILE.EQ.0) THEN
C Not using direct access files but we have run out of space
C ... can happen when IELL = NELL that IREC+K1 = ISIZE+1
C so no room for the flag which is to be set on the last call
                  INFO(1) = 5
                  IF (LP.GT.0) WRITE (LP,FMT=9010) INFO(1),IELL
                  IF (LP.GT.0) WRITE (LP,FMT=9020)
                  MKEY = MKEY + 1
                  IREC = 1
                  INFO(21) = 1
               ELSE
C Output buffer to direct access data set. Reset IREC to point
C to first free location in IBUFR.
                  MKEY = MKEY + 1
                  CALL MA42LD(-3,IFILE,MKEY,BUFR,LBUFR,IBUFR,LIBUFR,LP,
     +                       JFLAG)
                  IF (JFLAG.LT.0) GO TO 80
                  IREC = 1
               END IF
   50       CONTINUE
   60    CONTINUE
C
      ELSE
C Last call to routine.
C Set a flag in IBUFR to indicate end of file (IBUFR(IREC)=0).
C Write out contents of buffer.
         IBUFR(IREC) = 0
C Fill-in the rest of IBUFR with garbage (so that it is
C not undefined)
         MKEY = MKEY + 1
         IF (IFILE.NE.0 .AND. MKEY.EQ.1) THEN
            DO 70 I = IREC + 1,LIBUFR
               IBUFR(I) = -1
   70       CONTINUE
         END IF
         IREC = IREC + 1
         IF (IFILE.NE.0) THEN
            CALL MA42LD(-3,IFILE,MKEY,BUFR,LBUFR,IBUFR,LIBUFR,LP,JFLAG)
            IF (JFLAG.LT.0) GO TO 80
         END IF
      END IF
      GO TO 90

C **** Fatal error return ****
   80 INFO(1) = JFLAG
      IF (LP.GT.0) WRITE (LP,FMT=9000) INFO(1)
      GO TO 90

   85 CONTINUE
      IF (ABS(JFL).EQ.1) THEN
C Calculate number NBUFR of buffers required.
         ISPACE = ISIZE - IREC + 1
         JFIT = LENGTH - ISPACE
         NBUFR = 1
         JFIT1 = JFIT/ISIZE
         IF (JFIT.GT.0) NBUFR = 2 + JFIT1
         IF (JFIT.EQ.0 .AND. IELL.LT.NELL) NBUFR = 2
C Update MKEY and IREC.
         IF (NBUFR.EQ.1) THEN
            IREC = IREC + LENGTH
         ELSE
            MKEY = MKEY + NBUFR - 1
            IREC = 1 + JFIT - JFIT1*ISIZE
         END IF

      ELSE
C Final call
         MKEY = MKEY + 1
      END IF

   90 RETURN
 9000 FORMAT (' ***** Error return from MA42B/BD *****  INFO(1) = ',I3)
 9010 FORMAT (' ***** Error return from MA42B/BD *****  INFO(1) = ',I3,
     +       /' ***** Error discovered after input of elt/eqn *****',I8)
 9020 FORMAT (7X,'LENBUF(3) too small.',
     +        /7X,'Direct access data sets not requested.',
     +        /7X,'Continue computation to find space required.')
      END
C*******************************************************************
      SUBROUTINE MA42KD(JFL,IND,LENGTH,IFILE,ISIZE,MKEY,NUMBLK,IREC,
     +                  IELL,NELL,LP,INFO)
C
C     MA42K/KD  checks to see if a write to direct access is required
C     and, if so, checks that direct access has been requested and
C     that the direct access data set has sufficient records.
C
C     Argument list. * indicates the argument is changed by the routine.
C
C  JFL   - Integer variable. Indicates nature of call.
C          JFL = -1 following detection of singularity.
C          JFL = 2 on final call for a buffer.
C          JFL = 1 otherwise.
C  IND   - Integer variable. Used to indicate whether we are writing
C          rows of U or columns of L or integers.
C          IND = 1 rows of U.
C          IND = 2 columns of L.
C          IND = 3 integers.
C  LENGTH- Integer variable. Holds amount of data to be written.
C  IFILE - Integer variable. Holds stream number for direct
C          access data set.
C  ISIZE - Integer variable. Holds length of buffer.
C *MKEY  - Integer variable. Holds number of records written
C          to the direct access data set.
C  NUMBLK- Integer variable. Holds number of records in direct access
C          data set.
C *IREC  - Integer variable.  Points to first available space in
C          buffer.
C  IELL  - Integer variable. Holds current elt/equ number.
C  LP    - Integer variable. Holds stream number for error messages.
C *INFO  - Integer array of length 23.
C          INFO(1) is an error flag. Possible nonzero values
C          on return are 5, 6.
C          INFO(1) = 5 if buffer not long enough and no direct access
C                   data set requested.
C          INFO(1) = 6 if insufficient space allocated to direct access
C                   data set.
C          INFO(19), INFO(20), INFO(21) also altered by routine.
C
C  Local variables
C
C  INFOIN  - Integer variable. Holds value of INFO(1) on entry.
C  ISPACE  - Integer variable. Holds amount of free space in buffer.
C  JFIT    - Integer variable. Holds the amount of data
C            which will not fit in the current buffer.
C  JFIT1   - Integer variable. HOLDS JFIT/ISIZE.
C  NBUFR   - Integer variable. Holds number of buffers required.
C  NWRITE  - Integer variable. Holds the number of buffers
C            which must be written to the direct access file.
C
C     .. Scalar Arguments ..
      INTEGER IELL,IFILE,IND,IREC,ISIZE,JFL,LENGTH,LP,MKEY,NELL,NUMBLK
C     ..
C     .. Array Arguments ..
      INTEGER INFO(23)
C     ..
C     .. Local Scalars ..
      INTEGER INFOIN,ISPACE,JFIT,JFIT1,NBUFR,NWRITE
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX
C     ..
      INFOIN = INFO(1)
      IF (ABS(JFL).EQ.1) THEN
C Calculate number NBUFR of buffers required.
         ISPACE = ISIZE - IREC + 1
         JFIT = LENGTH - ISPACE
         NBUFR = 1
         JFIT1 = JFIT/ISIZE
         IF (JFIT.GT.0) NBUFR = 2 + JFIT1
         IF (JFIT.EQ.0 .AND. IELL.LT.NELL) NBUFR = 2
C NWRITE is the number of buffers which must be written to
C the direct access file (i.e. the number of full buffers).
         NWRITE = NBUFR
         IF (JFIT.GT.0 .AND. ISIZE*JFIT1.EQ.JFIT) NWRITE = NWRITE - 1
         IF (JFIT.EQ.0 .AND. IELL.LT.NELL) NWRITE = 1
C Update INFO(19), INFO(20), or INFO(21).
         IF (IND.EQ.1) INFO(19) = MAX(INFO(19),NWRITE)
         IF (IND.EQ.2) INFO(20) = MAX(INFO(20),NWRITE)
         IF (IND.EQ.3) INFO(21) = MAX(INFO(21),NWRITE)
         IF (NBUFR.GT.1) THEN
C Write to direct access data set is necessary.
            IF (IFILE.EQ.0) THEN
C No direct access data set requested.
               INFO(1) = 5
            ELSE IF (MKEY+NBUFR-1.GT.NUMBLK) THEN
C Insufficient records allocated to direct access data set.
               INFO(1) = 6
            END IF
         END IF
      END IF

      IF (INFO(1).NE.INFOIN .AND. LP.GT.0) THEN
         WRITE (LP,FMT=9010) INFO(1),IELL
         IF (INFO(1).EQ.5) WRITE (LP,FMT=9020) IND
         IF (INFO(1).EQ.6) WRITE (LP,FMT=9000) IND
      END IF
      RETURN
 9000 FORMAT (7X,'LENFLE(',I1,') too small.',
     +       /7X,'Continue computation to find space required.')
 9010 FORMAT (' ***** Error return from MA42B/BD *****  INFO(1) = ',I3,
     +       /' ***** Error discovered after input of elt/eqn *****',I8)
 9020 FORMAT (7X,'LENBUF(',I1,') too small.',
     +       /7X,'Direct access data sets not requested.',
     +       /7X,'Continue computation to find space required.')
      END
C**********************************************************************
      SUBROUTINE MA42LD(IOPT,IFILE,MKEY,BUFR,LBUFR,IBUFR,LIBUFR,LP,
     +                  JFLAG)
C
C This subroutine handles all the requests for direct access i/o.
C We read or write from direct access data set IFILE to a real
C in-core buffer BUFR of
C length LBUFR or to an integer in-core buffer IBUFR of length LIBUFR.
C MKEY indicates the direct access record to be written or read.
C IOPT indicates whether the real or integer buffer is to be used
C and whether a read or write is to be performed.
C
C     Argument list. * indicates the argument is changed by the routine.
C
C  IOPT   - Integer variable.  Control parameter.
C           IOPT can have six values ...
C           IOPT = +/-1.  Read from d.a. data set and move MKEY
C                        to next record.
C           IOPT = +/-2.  Read from d.a. data set and move MKEY to
C                         previous record.
C           IOPT = +/-3.  Write to direct access data set.
C           (Positive value of IOPT is for reals and
C           negative for integers)
C  IFILE  - Integer variable.  Stream number of data set being
C           read/written.
C *MKEY   - Integer variable.  Pointer to record of direct access data
C           set being written or read. Changed on exit if ABS(IOPT) = 1
C           or ABS(IOPT) = 2.
C  BUFR   - Real (DP) array of length LBUFR.
C           If  IOPT is positive this buffer is
C           being read or written from main storage to the direct
C           access data set.
C  LBUFR  - Integer variable.  Length of array BUFR.
C  IBUFR  - Integer array of length LIBUFR.
C           If  IOPT is negative this buffer is
C           being read or written from main storage to the direct
C           access data set.
C  LIBUFR - Integer variable.  Length of array IBUFR.
C  LP     - Integer variable. Stream number for error messages.
C *JFLAG   - Integer variable. Error flag.
C           If non-zero on exit, error has occurred in
C           direct access read/write. Possible nonzero values -25,-26.
C           JFLAG = -25 indicates failure in direct access read.
C           JFLAG = -26 indicates failure in direct access write.
C
C     .. Scalar Arguments ..
      INTEGER IFILE,IOPT,JFLAG,LBUFR,LIBUFR,LP,MKEY
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION BUFR(LBUFR)
      INTEGER IBUFR(LIBUFR)
C     ..
C     .. Local Scalars ..
      INTEGER IOS
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C Set JFLAG.
      JFLAG = 0
C Immediate return if IFILE is zero
      IF (IFILE.EQ.0) GO TO 30
C
      IF (IOPT.EQ.1 .OR. IOPT.EQ.2) THEN
C Real buffer. Perform direct access read.
         READ (IFILE,REC=MKEY,ERR=10,IOSTAT=IOS) BUFR
      ELSE IF (IOPT.EQ.3) THEN
C Real buffer. Perform direct access write.
         WRITE (IFILE,REC=MKEY,ERR=20,IOSTAT=IOS) BUFR
      ELSE IF (IOPT.EQ.-1 .OR. IOPT.EQ.-2) THEN
C Integer buffer. Perform direct access read.
         READ (IFILE,REC=MKEY,ERR=10,IOSTAT=IOS) IBUFR
      ELSE IF (IOPT.EQ.-3) THEN
C Integer buffer. Perform direct access write.
         WRITE (IFILE,REC=MKEY,ERR=20,IOSTAT=IOS) IBUFR
      END IF
C
C Update MKEY
      IF (ABS(IOPT).EQ.1) MKEY = MKEY + 1
      IF (ABS(IOPT).EQ.2) MKEY = MKEY - 1
C
      GO TO 30
C **** Fatal error returns *****
   10 JFLAG = -25
      IF (LP.GT.0) WRITE (LP,FMT=9000) MKEY,IFILE,IOS
      GO TO 30
   20 JFLAG = -26
      IF (LP.GT.0) WRITE (LP,FMT=9010) MKEY,IFILE,IOS
   30 RETURN
 9000 FORMAT (7X,'Error in direct access read.',/7X,'Record ',I8,' in ',
     +       'file ',I2,' IOSTAT ',I10)
 9010 FORMAT (7X,'Error in direct access write.',/7X,'Record ',I8,' in',
     +       ' file ',I2,' IOSTAT ',I10)
      END
C*********************************************************************
      SUBROUTINE MA42MD(IELL,AVAR,LRHS,NMAXE,RHS,NRHS,IVAR,LAST,NDF,
     +                  NVAR,MFRONT,NFRONT,FA,FRHS,LHED,KHED,KPIV,
     +                  KPVLNK,LDEST,KDEST,ISTATC,LFREE,KR,KFRNT,LFRNT,
     +                  INFO)
C
C Add in contributions from elt/eqn to front matrix and right hand
C sides.
C
C     Argument list. * indicates the argument is changed by the routine.
C
C     IELL      Integer variable. Holds current elt/equ number.
C     AVAR      Real (DP) array of dimensions NMAXE by NVAR.
C               contains contributions from current elt/eqn.
C     LRHS      Integer variable. Second dimension of FRHS array.
C     NMAXE     Integer variable. Leading dimension of element or
C               equation arrays.
C     RHS       Real (DP) array of dimensions NMAXE by LRHS.
C               Contains right hand sides for elt/eqn being input.
C     IVAR      Integer array of length NVAR.
C               Contains indices of variables in current elt/eqn.
C               It is used locally as a work array by the routine
C               to facilitate access to the position of an incoming
C               variable in the front but it is reset and returned
C               unchanged to MA42F/FD.
C     NDF       Integer variable. Largest integer used to index a
C               variable.
C    *LAST      Integer array of length NDF.
C               For a variable not yet in the front, on entry
C               LAST(I) is the elt/eqn in which variable I
C               appears for the last time. For a non-fully summed
C               variable in the front, JVAR say, -LAST(JVAR) is the
C               position in KDEST/LDEST of the information for
C               this variable.
C     NVAR      Integer variable. Number of variables in current
C               elt/equ.
C     MFRONT    Integer variable. Max. number of rows in frontal matrix.
C     NFRONT    Integer variable. Max. number of cols in frontal matrix.
C    *FA        Real (DP) array with dimensions MFRONT,NFRONT.
C               Holds the current frontal matrix
C    *FRHS      Real (DP) array with dimensions MFRONT,LRHS.
C               Holds the right hand sides corresponding to the
C               current frontal matrix
C    *LHED      Integer array of dimension NFRONT. Used to hold the
C               indices of variables corresponding to columns in the
C               front.
C    *KHED      Integer array of dimension MFRONT. Used to hold the
C               indices of variables corresponding to rows in the
C               front.
C    *KPIV      Integer array of dimension NFRONT.  Indicates which
C               columns of FA are fully summed.
C    *KPVLNK    Integer array of dimension NFRONT. For each fully
C               summed column in the front, KPVLNK indicates its
C               position in array KPIV.
C    *LDEST     Integer array of length NFRONT.  For each non-fully
C               summed variable in the front (global name JVAR say)
C               LDEST(-LAST(JVAR)) is its column position in the frontal
C               matrix.
C    *KDEST     Integer array of length MFRONT.  For each non-fully
C               summed variable in the front (global name JVAR say)
C               KDEST(-LAST(JVAR)) is its row position in the frontal
C               matrix.
C     ISTATC    Integer variable. Holds number of variables which are
C               static condensations candidates for current elt/equ.
C    *LFREE     Integer variable. Holds position of first free location
C               in KPVLNK.
C    *KR        Integer variable. On exit holds number of fully
C               assembled variables.
C    *LFRNT     Integer variable. On exit holds number of columns in
C               front. (Only changed if ISTATC=0)
C    *KFRNT     Integer variable. On exit holds number of rows in front.
C               (Only changed if ISTATC=0)
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER IELL,ISTATC,KFRNT,KR,LFREE,LFRNT,LRHS,MFRONT,NDF,NFRONT,
     +        NMAXE,NRHS,NVAR
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION AVAR(NMAXE,NVAR),FA(MFRONT,NFRONT),
     +                 FRHS(MFRONT,LRHS),RHS(NMAXE,LRHS)
      INTEGER INFO(23),IVAR(NVAR),KDEST(MFRONT),KHED(MFRONT),
     +        KPIV(NFRONT),KPVLNK(NFRONT),LAST(NDF),LDEST(NFRONT),
     +        LHED(NFRONT)
C     ..
C     .. Local Scalars ..
      INTEGER IPOS,J,JVAR,K,KC1,KDST,KPOS,KR1,L,LDST,LELL,LFR,LPOS,MFR,
     +        NFREE,NUVAR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX,MIN
C     ..
      IF (ISTATC.EQ.0) THEN
C Place incoming variables not already in front in last rows and
C columns of front matrix.
         IF (NMAXE.EQ.1) THEN
C Equation entry.
            KFRNT = KFRNT + 1
            KHED(KFRNT) = IELL
         END IF
         DO 10 L = 1,NVAR
            MFR = IVAR(L)
            IF (LAST(MFR).GE.0) THEN
C Insert incoming variable in the front.
C Set KDEST/LDEST and KHED/LHED.
               LFRNT = LFRNT + 1
               NFREE = LDEST(LFREE)
               LDEST(LFREE) = LFRNT
               LHED(LFRNT) = MFR
               IF (NMAXE.GT.1) THEN
C Element entry.
                  KFRNT = KFRNT + 1
                  KHED(KFRNT) = MFR
                  KDEST(LFREE) = KFRNT
               END IF
C Save value of LAST and set flag in LAST to indicate variable is
C now in front.
               KPVLNK(LFRNT) = -LAST(MFR)
               LAST(MFR) = -LFREE
               LFREE = NFREE
            END IF
   10    CONTINUE
C Now zero out any so far unused portion of the arrays FA and
C FRHS.  This is to avoid an MFRONT*NFRONT initialization loop.
C KR1 and KC1 hold one more than the max number of rows and
C columns in the front so far.
         KR1 = INFO(8) + 1
         KC1 = INFO(9) + 1
         IF (KR1.LE.KFRNT) THEN
C Update INFO(8).
            INFO(8) = KFRNT
            LFR = MAX(INFO(9),LFRNT)
C Zero out rows unused rows KR1 to KFRNT, columns 1 to LFR
C and corresponding right-hand sides.
            DO 40 K = KR1,KFRNT
               DO 20 L = 1,LFR
                  FA(K,L) = ZERO
   20          CONTINUE
               DO 30 J = 1,NRHS
                  FRHS(K,J) = ZERO
   30          CONTINUE
   40       CONTINUE
         END IF
         IF (KC1.LE.LFRNT) INFO(9) = LFRNT
C Zero out columns KC1 to LFRNT, rows 1 to KR1.
         DO 60 L = KC1,LFRNT
            DO 50 K = 1,MIN(MFRONT,KR1)
               FA(K,L) = ZERO
   50       CONTINUE
   60    CONTINUE
      END IF
C
C Use array IVAR to hold temporarily the column position of the
C variable in the front.  This saves some indirect addressing during
C assembly.
      NUVAR = NVAR - ISTATC
      DO 70 L = 1,NUVAR
         MFR = IVAR(L)
         IPOS = -LAST(MFR)
C LDST is the column position of the variable in the front.
         LDST = LDEST(IPOS)
         IVAR(L) = LDST
   70 CONTINUE
C
C Add into front matrix contributions from incoming elt/eqn
C (contributions from non-static condensation variables).
C Jump if equation entry with one static
C condensation variable (contributions already added).
      IF (NMAXE.EQ.1 .AND. ISTATC.EQ.1) GO TO 130
      IF (NMAXE.EQ.1) THEN
C Equation entry (ISTATC = 0)
         DO 80 K = 1,NVAR
            LDST = IVAR(K)
            FA(KFRNT,LDST) = FA(KFRNT,LDST) + AVAR(1,K)
   80    CONTINUE
         DO 90 J = 1,NRHS
            FRHS(KFRNT,J) = FRHS(KFRNT,J) + RHS(1,J)
   90    CONTINUE
      ELSE
C Element entry
         DO 120 L = 1,NUVAR
C This code is needed to find the row position in the front of the
C variable in position L of the incoming element. LDST is the column
C position, JVAR the global variable associated with this column, KPOS
C the position in the destination vectors KDEST/LDEST, and KDST the
C sought for row of the front matrix.
            LDST = IVAR(L)
            JVAR = LHED(LDST)
            KPOS = -LAST(JVAR)
            KDST = KDEST(KPOS)
            DO 100 K = 1,NUVAR
               LDST = IVAR(K)
               FA(KDST,LDST) = FA(KDST,LDST) + AVAR(L,K)
  100       CONTINUE
            DO 110 J = 1,NRHS
               FRHS(KDST,J) = FRHS(KDST,J) + RHS(L,J)
  110       CONTINUE
  120    CONTINUE
      END IF
C
C Place all variables appearing for the last time in
C KPIV/KPVLNK, restore LAST, and remove from KDEST/LDEST.
  130 DO 140 L = 1,NUVAR
         LDST = IVAR(L)
         JVAR = LHED(LDST)
         LELL = -KPVLNK(LDST)
         IF (LELL.LE.IELL) THEN
            KR = KR + 1
            KPIV(KR) = LDST
            KPVLNK(LDST) = KR
            LPOS = -LAST(JVAR)
            LAST(JVAR) = IELL
            LDEST(LPOS) = LFREE
            LFREE = LPOS
         END IF
C Restore IVAR.
         IVAR(L) = JVAR
  140 CONTINUE
      END
C**********************************************************************
      SUBROUTINE MA42ND(IELL,NDF,LAST,MFRONT,NFRONT,FA,LHED,KHED,LRHS,
     +                  NRHS,FRHS,KDEST,LDEST,KPIV,KPVLNK,KPIVRX,LPIVCX,
     +                  PIVOT,KFRNT,LFRNT,KR,KS,NMAXE,NELIM,IFORCE,
     +                  LIBUFR,IBUFR,BUFRL,LLB,BUFRU,LUB,IFILE,IREC,
     +                  ISIZE,MKEY,NUMBLK,DET,OPS,NELL,CNTL,ICNTL,INFO,
     +                  OFDIAG,NSIZE,LSTAT)
C
C    This routine choses pivots and performs Gauss elimination.
C    Off-diagonal pivoting allowed unless ABS(ICNTL(7))=1001. For
C    off-diagonal pivoting, if there is a potential pivot in
C    a fully summed row and column entry which is greater than
C    CNTL(2) times the largest entry in its row, then use it as
C    pivot. If not, calculate the best potential pivot which does not
C    satisfy the pivot criterion for use in possible
C    later forced eliminations, or if the ICNTL(5) is positive
C    then force elimination immediately.
C    If ABS(ICNTL(7))=1001 pivots can only be chosen from the
C    diagonal, unless IELL=NELL and no more diagonal pivots
C    can be chosen while maintaining stability.
C
C     Argument list. * indicates the argument is changed by the routine.
C
C     IELL      Integer variable. Holds current elt/equ number.
C     NDF       Integer variable. Total number of variables.
C    *LAST      Integer array of length NDF.
C               Unchanged on exit unless error return INFO=-14
C     MFRONT    Integer variable. Max. number of rows in frontal matrix.
C     NFRONT    Integer variable. Max. number of cols in frontal matrix.
C     FA        Real (DP) array with dimensions MFRONT,NFRONT.
C               Holds the frontal matrix.
C    *LHED      Integer array of dimension NFRONT. Holds the indices
C               of variables corresponding to columns in the front.
C               The sign of LHED(L) is used to indicate whether
C               column L of the frontal matrix has zeros within
C               the pivotal rows. Signs are restored before exit.
C    *KHED      Integer array of dimension MFRONT. Holds the indices
C               of variables corresponding to rows in the front.
C               The sign of KHED(K) is used to indicate whether
C               row K of the frontal matrix has zeros within
C               the pivotal columns. Signs are restored before exit.
C     LRHS      Integer variable. Second dimension of FRHS array.
C     NRHS      Integer variable.  Number of right hand sides.
C    *FRHS      Real (DP) array with dimensions MFRONT,LRHS.
C               Holds the right hand sides corresponding to the
C               current frontal matrix (updated by elimination).
C    *LDEST     Integer array of length NFRONT.  For each non-fully
C               summed variable in the front (global name JVAR say)
C               LDEST(-LAST(JVAR)) is its column position in the frontal
C               matrix.
C    *KDEST     Integer array of length MFRONT.  For each non-fully
C               summed variable in the front (global name JVAR say)
C               KDEST(-LAST(JVAR)) is its row position in the frontal
C               matrix.
C    *KPIV      Integer array of dimension NFRONT.  Indicates which
C               columns of FA are fully summed. (i.e. potential pivots).
C    *KPVLINK   Integer array of dimension NFRONT. For each fully summed
C               column in the front, KPVLNK indicates its position
C               in array KPIV.
C    *KPIVRX    Integer variable. Used to hold row index of best
C               available pivot.
C    *LPIVCX    Integer variable. Used to hold column index of best
C               available pivot.
C    *PIVOT     Real (DP) variable. On exit, holds
C               value of pivot.
C     KFRNT     Integer variable. Holds the number of rows in the front.
C     LFRNT     Integer variable. Holds the number of columns in
C               the front.
C    *KR        Integer variable. Holds the number of fully assembled
C               variables.
C    *KS        Integer variable. Used in cyclic search of columns for
C               a pivot column.
C     NMAXE     Integer variable. Leading dimension of element or
C               equation arrays. NMAXE=1 for equation entry.
C     NELIM     Integer variable. Holds number of eliminations to be
C               attempted (NELIM= IFORCE if IFORCE is nonzero and
C               is equal to value of KR on entry otherwise).
C     IFORCE    Integer variable. Holds number of forced eliminations
C               required for current elt/equ.
C     ILBUFR    Integer variable. Length of array IBUFR.
C    *IBUFR     Integer array of length LIBUFR. The output buffer.
C    *BUFRL     Real (DP) array of length LLB. In-core buffer for
C               entries of L.
C     LLB       Integer variable. Length of array BUFRL.
C    *BUFRU     Real (DP) array of length LUB. In-core buffer for
C               entries of U.
C     LUB       Integer variable. Length of array BUFRU.
C     ICNTL    Integer array of length 8. See MA432B/BD for details.
C      CNTL     Real (DP) array of length 2. Matrix singular if
C               largest element in column less than CNTL(1).
C               An element of the frontal matrix will only be considered
C               suitable for use as a pivot if it is at least as large
C               as CNTL(2) times the largest element in its column.
C    *DET       Real (DP) variable.  Used to accumulate log
C               of determinant of matrix.
C    *OPS       Real (DP) variable. Holds number of operations
C               in innermost loop.
C    *INFO      Integer array of length 23. INFO(1) is an error flag.
C               On exit INFO(1) = 4 if eliminations are being forced
C               but only zero pivots are available. (Can only happen
C               if IFORCE > 0 and ICNTL(8) non-zero).
C               On exit INFO(1) = -14 if element matrix is singular.
C               INFO(1) = -26, 5, 6 may be returned from MA42G/GD or
C               MA42H/HD.
C   *OFDIAG     Integer variable. Holds number of off-diag.
C               pivots chosen during the call.

C  Local variables
C
C  GMON    - Real (DP) variable. Stores the ratio of the pivot
C            to the maximum value in its column.
C  GMONX   - Real (DP) variable. Stores the ratio of the pivot
C            to the maximum value in its column
C            for the best pivot found so far.
C  PIVINV  - Real (DP) variable. Holds inverse of pivot.
C  RMAX    - Real (DP) variable. Max. entry in pivot column.
C  SWAP    - Real (DP) variable. Used in swapping pivot row and
C            col with last row and col of front resp.
C  CHFRNT  - Integer variable. Holds change in front size.
C  IAUTO   - Integer variable. Equivalent to ICNTL(5).
C  IELIM   - Integer variable. Do loop variable. Ranges from 1 to
C             NELIM.
C  ISRCH   - Integer variable. Equvalent to ICNTL(6).
C  JSRCH   - Integer variable. Holds number of columns to be
C            searched for a pivot.
C  K       - Integer variable. Do loop variable. Loop over the
C            the rows in the front.
C  KFR     - Integer variable. Holds number of rows in front
C            after elimination.
C  KFRO    - Integer variable. Holds global row index of variable
C            in column KFR.
C  KH      - Integer variable. Used to hold position in fully-summed
C            list of column in which pivot candidate lies.
C  KMAX    - Integer variable. Holds (local) row index of largest
C            entry in pivot column.
C  KPIVRO  - Integer variable. Holds (local) row index of pivot.
C  KX      - Integer variable. Holds position of pivot in KPIV.
C  LPIVCO    Integer variable. Holds (local) column index of pivot.
C  KPIVKR  - Integer variable. Used to adjust KPVLNK to remove pivot
C            column.
C  KPOS    - Integer variable. Holds position in LDEST of variable KFRO.
C  LP      - Integer variable. Stream for error messages. (ICNTL(1)).
C  L       - Integer variable. Do loop variable. Loop over the
C            columns which are searched for a pivot.
C  LCOL    - Integer variable. Holds global column index of fully
C            assembled col.
C  LFR     - Integer variable. Holds number of columns
C            in front after elimination.
C  LFRO    - Integer variable. Holds global column index of variable
C            in column LFR.
C  LK      - Integer variable. Do loop variable.
C  LPOS    - Integer variable. Holds position in LDEST of variable LFRO.
C  LTEMP   - Integer variable. Used in swapping pivot row and
C            col with last row and col of front resp.
C  MFR     - Integer variable. Used to hold global row index of
C            fully assembled rows.
C  MP      - Integer variable. Stream for warnings (ICNTL(2)).
C  PIVX    - Real (DP) variable. Used to hold the best pivot
C            available.
C
C     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D0,ZERO=0.0D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION DET,OPS,PIVOT
      INTEGER IELL,IFORCE,KFRNT,KPIVRX,KR,KS,LFRNT,LIBUFR,LLB,LPIVCX,
     +        LRHS,LUB,MFRONT,NDF,NELIM,NELL,NFRONT,NMAXE,NRHS,NSIZE,
     +        OFDIAG
      LOGICAL LSTAT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION BUFRL(LLB),BUFRU(LUB),CNTL(2),FA(MFRONT,NFRONT),
     +                 FRHS(MFRONT,LRHS)
      INTEGER IBUFR(LIBUFR),ICNTL(8),IFILE(3),INFO(23),IREC(3),ISIZE(3),
     +        KDEST(MFRONT),KHED(MFRONT),KPIV(NFRONT),KPVLNK(NFRONT),
     +        LAST(NDF),LDEST(NFRONT),LHED(NFRONT),MKEY(3),NUMBLK(3)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION GMON,GMONX,PIVINV,PIVX,RMAX,SWAP
      INTEGER CHFRNT,I,IAUTO,IELIM,ISRCH,J,JJ,JSRCH,K,KFR,KFRO,KH,KK,
     +        KMAX,KPIVKR,KPIVRO,KPOS,KPRE,KRTEMP,KX,L,LCOL,LFR,LFRO,LK,
     +        LL,LP,LPIVCO,LPOS,LPRE,LTEMP,MFR,MP,L1
C     ..
C     .. External Functions ..
      INTEGER IDAMAX
      EXTERNAL IDAMAX
C     ..
C     .. External Subroutines ..
      EXTERNAL DAXPY,DGEMM,DGER,DTRSM,MA42GD,MA42HD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DBLE,LOG,MIN
C     ..
C Immediate return if NELIM is 0.
      IF (NELIM.EQ.0) GO TO 420
c      OFDIAG = 0
      LP = ICNTL(1)
      MP = ICNTL(2)
      IAUTO = ICNTL(5)
      ISRCH = ICNTL(6)
      KFR = KFRNT
      LFR = LFRNT

      KRTEMP = KR

      DO 140 IELIM = 1,NELIM
C If KFR or LFR is 0 the matrix is singular.
         IF (KFR.EQ.0 .OR. LFR.EQ.0) THEN
            INFO(2) = 0
            DET = ZERO
C Jump to error return if ICNTL(8) is 0.
            IF (ICNTL(8).EQ.0) GO TO 410
C Otherwise set INFO(1) = 1 and issue a warning.
            IF (INFO(1).EQ.0) THEN
               INFO(1) = 1
               IF (MP.GT.0) THEN
                  WRITE (MP,FMT=9020) INFO(1),IELL
                  WRITE (MP,FMT=9000)
               END IF
            END IF
            GO TO 150
         END IF
C Search for pivot.
         JSRCH = KR
C Restrict pivot search if IAUTO is positive, KR is at
C least IAUTO, and ISRCH is positive
         IF (IAUTO.GT.0 .AND. KR.GE.IAUTO .AND.
     +       ISRCH.GT.0) JSRCH = MIN(KR,ISRCH)
C Search for pivot among fully assembled columns.
C Jump if off-diagonal pivots are allowed.
         IF (ABS(ICNTL(7)).NE.1001) GO TO 20
C Pivots to be chosen from the diagonal.
C GMONX stores the ratio of the pivot to the maximum value in its
C column for the best pivot found so far.
C KX stores the position of this entry in KPIV.
         GMONX = ZERO
         KX = 0
         DO 10 L = 1,JSRCH
C KS controls cyclic search of the pivot columns.
            KS = KS + 1
            IF (KS.GT.KR) KS = 1
            LPIVCO = KPIV(KS)
            KPIVRO = LPIVCO
            PIVOT = FA(KPIVRO,LPIVCO)
C Accumulate statistics on pivoting.
            INFO(13) = INFO(13) + 1
            INFO(15) = INFO(15) + KFR
C KMAX is the row index of the largest entry in column LPIVCO.
            KMAX = IDAMAX(KFR,FA(1,LPIVCO),1)
C RMAX is absolute value of largest entry in pivot column LPIVCO
            RMAX = ABS(FA(KMAX,LPIVCO))
C If all values in column less than tolerance the
C matrix is considered to be singular.
            IF (RMAX.LE.CNTL(1)) THEN
               INFO(2) = 0
               DET = ZERO
C Jump to error return if ICNTL(8) is 0.
               IF (ICNTL(8).EQ.0) GO TO 410
C Otherwise set INFO(1) = 1 and issue a warning.
               IF (INFO(1).EQ.0) THEN
                  INFO(1) = 1
                  IF (MP.GT.0) THEN
                     WRITE (MP,FMT=9020) INFO(1),IELL
                     WRITE (MP,FMT=9000)
                  END IF
               END IF
               GO TO 10
            END IF
            INFO(14) = INFO(14) + 1
C Hold position in fully-summed list of column in which current pivot
C candidate lies.  We must do this in case current column
C fails to yield a better choice.  N.B. The previous choice
C had a multiplier less than CNTL(2).
            KH = KX
            KX = KS
            GMON = ABS(PIVOT)/RMAX
C Jump if GMON > CNTL(2) (pivot found).
            IF (GMON.GT.CNTL(2)) GO TO 60
            KX = KH
            IF (GMON.GT.GMONX) THEN
C Store information on best diagonal pivot candidate so far.
               PIVX = PIVOT
               KX = KS
               GMONX = GMON
               LPIVCX = LPIVCO
               KPIVRX = KPIVRO
            END IF
   10    CONTINUE
C No suitable diagonal pivot.
C If IELL=NELL and pivots are not being forced (i.e. this is
C the very last call to MA42N/ND),
C switch to looking for off-diagonal pivots.
         IF (IELL.EQ.NELL .AND. IFORCE.EQ.0) GO TO 20
C Otherwise, need to assemble another elt/eqn into the front
C or use best available diagonal pivot if we need to
C force eliminations.
         IF (ABS(PIVOT).LE.CNTL(1)) THEN
C Best available diagonal pivot is too small to use.
            IF (IFORCE.EQ.0) GO TO 150
C Eliminations were being forced, so there is not room to assemble
C another elt/equ. Give warning INFO(1)= 4 and return.
            INFO(1) = 4
            GO TO 420
         ELSE IF (IFORCE.EQ.0 .AND. (IAUTO.EQ.0.OR.KR.LT.IAUTO)) THEN
C Do not need to force eliminations  and either ICNTL(5)=0 or
C the number of fully assembled variables is less than ICNTL(5).
C Perform eliminations with pivots already chosen and then
C return and assemble another elt/eqn into the front.
            GO TO 150
         ELSE
C Choose best available diagonal pivot if we need to force
C eliminations to create space in front or IAUTO is positive and the
C number of fully assembled columns KR is large enough
            INFO(16) = INFO(16) + 1
            PIVOT = PIVX
            LPIVCO = LPIVCX
            KPIVRO = KPIVRX
            GO TO 60
         END IF
C
C Look for pivots (on or off the diagonal)
   20    CONTINUE
C GMONX stores the ratio of the pivot to the maximum value in its
C     column for the best pivot found so far.
C KX stores the position of this entry in KPIV.
         GMONX = ZERO
         KX = 0
         PIVOT = ZERO
C In this loop we are looking only at largest entry in column to
C to see if can be used as a pivot (trying to reduce search time
C by not looking for largest entry in fully summed part unless
C we have to)
         DO 30 L = 1,JSRCH
C KS controls cyclic search of the pivot columns.
            KS = KS + 1
            IF (KS.GT.KR) KS = 1
            LPIVCO = KPIV(KS)
C Accumulate statistics on pivoting.
            INFO(13) = INFO(13) + 1
            INFO(15) = INFO(15) + KFR
            KMAX = IDAMAX(KFR,FA(1,LPIVCO),1)
C RMAX is absolute value of largest entry in pivot column LPIVCO
            RMAX = ABS(FA(KMAX,LPIVCO))
C If all values in column less than tolerance the
C matrix is considered to be singular.
            IF (RMAX.LE.CNTL(1)) THEN
               INFO(2) = 0
               DET = ZERO
C Jump to error return if ICNTL(8) is 0.
               IF (ICNTL(8).EQ.0) GO TO 410
C Otherwise set INFO(1) = 1 and issue a warning.
               IF (INFO(1).EQ.0) THEN
                  INFO(1) = 1
                  IF (MP.GT.0) THEN
                     WRITE (MP,FMT=9020) INFO(1),IELL
                     WRITE (MP,FMT=9000)
                  END IF
               END IF
               GO TO 30
            END IF
C Hold position in fully-summed list of column in which current pivot
C candidate lies.  We must do this in case current column
C fails to yield a better choice.  N.B. The previous choice
C had a multiplier less than CNTL(2).
            KH = KX
            KX = KS
            MFR = KHED(KMAX)
C Jump if equation entry or largest entry is fully summed.
            IF (NMAXE.EQ.1) THEN
               KPIVRO = KMAX
               PIVOT = FA(KPIVRO,LPIVCO)
               GO TO 60
            ELSE IF (LAST(MFR).GE.0) THEN
               KPIVRO = KMAX
               PIVOT = FA(KPIVRO,LPIVCO)
               GO TO 60
            END IF
   30    CONTINUE
         DO 50 L = 1,JSRCH
C KS controls cyclic search of the pivot columns.
            KS = KS + 1
            IF (KS.GT.KR) KS = 1
            LPIVCO = KPIV(KS)
C Accumulate statistics on pivoting.
            INFO(13) = INFO(13) + 1
            INFO(15) = INFO(15) + KFR
            KMAX = IDAMAX(KFR,FA(1,LPIVCO),1)
C RMAX is absolute value of largest entry in pivot column LPIVCO
            RMAX = ABS(FA(KMAX,LPIVCO))
C If all values in column less than tolerance the
C matrix is considered to be singular.
            IF (RMAX.LE.CNTL(1)) THEN
               INFO(2) = 0
               DET = ZERO
C Jump to error return if ICNTL(8) is 0.
               IF (ICNTL(8).EQ.0) GO TO 410
C Otherwise set INFO(1) = 1 and issue a warning.
               IF (INFO(1).EQ.0) THEN
                  INFO(1) = 1
                  IF (MP.GT.0) THEN
                     WRITE (MP,FMT=9020) INFO(1),IELL
                     WRITE (MP,FMT=9000)
                  END IF
               END IF
               GO TO 50
            END IF
C Hold position in fully-summed list of column in which current pivot
C candidate lies.  We must do this in case current column
C fails to yield a better choice.  N.B. The previous choice
C had a multiplier less than CNTL(2).
            KH = KX
            KX = KS
            MFR = KHED(KMAX)
C Jump if equation entry or largest entry is fully summed.
            IF (NMAXE.EQ.1) THEN
               KPIVRO = KMAX
               PIVOT = FA(KPIVRO,LPIVCO)
               GO TO 60
            ELSE IF (LAST(MFR).GE.0) THEN
               KPIVRO = KMAX
               PIVOT = FA(KPIVRO,LPIVCO)
               GO TO 60
            END IF
C Largest entry was not fully-summed.
C Find the largest entry in the fully-summed part.
            PIVOT = ZERO
C Run through column choosing as potential pivot the largest
C value in the fully assembled rows. Update INFO(14).
            DO 40 K = 1,KFR
               MFR = KHED(K)
               IF (LAST(MFR).GE.0) THEN
C Fully-summed variable and hence available for use as pivot.
                  INFO(14) = INFO(14) + 1
                  IF (ABS(FA(K,LPIVCO)).GE.ABS(PIVOT)) THEN
                     KPIVRO = K
                     PIVOT = FA(KPIVRO,LPIVCO)
                  END IF
               END IF
   40       CONTINUE
            GMON = ABS(PIVOT)/RMAX
C Jump if GMON > CNTL(2) (pivot found).
            IF (GMON.GT.CNTL(2)) GO TO 60
            KX = KH
            IF (GMON.GT.GMONX) THEN
C Store information on best pivot candidate so far.
               PIVX = PIVOT
               KX = KS
               GMONX = GMON
               LPIVCX = LPIVCO
               KPIVRX = KPIVRO
            END IF
   50    CONTINUE
         IF (ABS(PIVOT).LE.CNTL(1)) THEN
C All columns are zero columns. Need to
C assemble another elt/eqn into the front. If eliminations were
C being forced there is not room to do this give warning INFO(1)=4.
            IF (IFORCE.EQ.0) GO TO 150
            INFO(1) = 4
            GO TO 420
         ELSE IF (IFORCE.EQ.0 .AND. (IAUTO.EQ.0.OR.KR.LT.IAUTO)) THEN
C Do not need to force eliminations  and either ICNTL(5)=0 or
C the number of fully assembled variables is less than ICNTL(5).
C Perform eliminations with pivots already chosen and then
C return and assemble another elt/eqn into the front.
            GO TO 150
         ELSE
C Choose best available pivot if we need to force eliminations to
C create space in front or IAUTO is positive and the number
C of fully assembled columns KR is large enough
            INFO(16) = INFO(16) + 1
            PIVOT = PIVX
            LPIVCO = LPIVCX
            KPIVRO = KPIVRX
         END IF
   60    CONTINUE
C Pivot chosen.
         IF (PIVOT.LT.ZERO) INFO(2) = -INFO(2)
C Update OFDIAG
         IF (KPIVRO.NE.LPIVCO) OFDIAG = OFDIAG + 1

C Adjust KPIV and KPVLNK to remove pivot column from
C fully assembled list.
         IF (KPVLNK(LFR).GT.0) KX = KPVLNK(LFR)
         KPIVKR = KPIV(KR)
         KPIV(KX) = KPIVKR
         KPVLNK(KPIVKR) = KX
         IF (KPVLNK(LFR).LT.0) KPVLNK(LPIVCO) = KPVLNK(LFR)
C
C Change KDEST/LDEST according to pivot permutation.
         LFRO = LHED(LFR)
         LPOS = -LAST(LFRO)
         IF (LPOS.GT.0) LDEST(LPOS) = LPIVCO
         IF (NMAXE.GT.1) THEN
            KFRO = KHED(KFR)
            KPOS = -LAST(KFRO)
            IF (KPOS.GT.0) KDEST(KPOS) = KPIVRO
         END IF
C Swap pivotal row and col with row KFR and col LFR of FA resp.
         IF (KFR.NE.KPIVRO) THEN
            LTEMP = KHED(KFR)
            KHED(KFR) = KHED(KPIVRO)
            KHED(KPIVRO) = LTEMP

            DO 90 L = 1,LFRNT
               SWAP = FA(KFR,L)
               FA(KFR,L) = FA(KPIVRO,L)
               FA(KPIVRO,L) = SWAP
   90       CONTINUE
            INFO(2) = -INFO(2)
         END IF
         IF (LFR.NE.LPIVCO) THEN
            LTEMP = LHED(LFR)
            LHED(LFR) = LHED(LPIVCO)
            LHED(LPIVCO) = LTEMP

            DO 100 K = 1,KFRNT
               SWAP = FA(K,LFR)
               FA(K,LFR) = FA(K,LPIVCO)
               FA(K,LPIVCO) = SWAP
  100       CONTINUE
            INFO(2) = -INFO(2)
         END IF
C Change sign of pivotal column.
         DO 110 K = 1,KFR - 1
            FA(K,LFR) = -FA(K,LFR)
  110    CONTINUE
C Store pivot in FA(KFR,LFR)
         FA(KFR,LFR) = -PIVOT
C Update columns in FA corresponding to fully assembled variables
C (we only do this if we are going to look for more pivots)
         KR = KR - 1

C This is if we want to force one elimination at a time.
         IF (NSIZE.EQ.1) THEN
            KRTEMP = KR
            KR = 0
         END IF

         PIVINV = ONE/FA(KFR,LFR)
         DO 120 L = 1,KR
            LCOL = KPIV(L)
C Scale entries in pivotal row by pivot.
            FA(KFR,LCOL) = -FA(KFR,LCOL)*PIVINV
            CALL DAXPY(KFR-1,FA(KFR,LCOL),FA(1,LFR),1,FA(1,LCOL),1)
  120    CONTINUE
C
         IF (KFR.NE.KPIVRO) THEN
C Swap right-hand sides in rows KFR and KPIVRO.
            DO 130 L = 1,NRHS
               SWAP = FRHS(KFR,L)
               FRHS(KFR,L) = FRHS(KPIVRO,L)
               FRHS(KPIVRO,L) = SWAP
  130       CONTINUE
         END IF
C
C Update DET
         DET = DET + LOG(ABS(PIVOT))
C
C Reduce front width
         KFR = KFR - 1
         LFR = LFR - 1
C
  140 CONTINUE

C All pivots now chosen.
C
  150 CHFRNT = LFRNT - LFR
C KR is now the number of fully assembled columns not used as pivotal
C columns. If KR is nonzero, need to do some column interchanges to
C get the fully assembled columns to the last columns of the
C frontal matrix.
      DO 200 L = 1,KR
         LCOL = KPIV(L)
C If LCOL already lies in the last block of columns, a swap is
C unnecessary.
         IF (LCOL.GE.LFR-KR+1) GO TO 200
C Look for a non-fully summed column to swap.
         DO 190 LK = LFR,LFR - KR + 1,-1
            IF (KPVLNK(LK).LT.0) THEN
C Column LK is not fully summed and we may swap.
               KPVLNK(LCOL) = KPVLNK(LK)
               KPIV(L) = LK
               KPVLNK(LK) = L
               IF (LPIVCX.EQ.LCOL) LPIVCX = LK
C Change LDEST according to permutation.
               LFRO = LHED(LK)
               LPOS = -LAST(LFRO)
               LDEST(LPOS) = LCOL
C Swap columns LCOL and LK
               LTEMP = LHED(LK)
               LHED(LK) = LHED(LCOL)
               LHED(LCOL) = LTEMP

               DO 160 K = 1,KFRNT
                  SWAP = FA(K,LK)
                  FA(K,LK) = FA(K,LCOL)
                  FA(K,LCOL) = SWAP
  160          CONTINUE
               INFO(2) = -INFO(2)
C If diagonal pivoting is being used, row swaps
C must also be performed
               IF (ABS(ICNTL(7)).EQ.1001 .AND. KPIVRO.EQ.LPIVCO) THEN
                  IF (NMAXE.GT.1) THEN
                     KFRO = KHED(LK)
                     KPOS = -LAST(KFRO)
                     KDEST(KPOS) = LCOL
                  END IF
C Swap rows LCOL and LK
                  LTEMP = KHED(LK)
                  KHED(LK) = KHED(LCOL)
                  KHED(LCOL) = LTEMP
                  DO 170 K = 1,LFRNT
                     SWAP = FA(LK,K)
                     FA(LK,K) = FA(LCOL,K)
                     FA(LCOL,K) = SWAP
  170             CONTINUE
                  INFO(2) = -INFO(2)
C Swap right-hand sides in rows LCOL and LK.
                  DO 180 K = 1,NRHS
                     SWAP = FRHS(LK,K)
                     FRHS(LK,K) = FRHS(LCOL,K)
                     FRHS(LCOL,K) = SWAP
  180             CONTINUE
               END IF
               GO TO 200
            END IF
  190    CONTINUE
  200 CONTINUE
C
C Now have the pivot rows/columns in the last rows/columns
C of the frontal matrix.
C If the number of pivot rows/columns
C (that is, the change in the frontsize)
C is odd AND the remaining frontsize is odd,
C need to change the sign of the determinant
      IF ((CHFRNT/2)*2.NE.CHFRNT) THEN
         IF ((KFR/2)*2.NE.KFR) INFO(2) = -INFO(2)
         IF ((LFR/2)*2.NE.LFR) INFO(2) = -INFO(2)
      END IF

      KPRE = 0
      LPRE = 0
      IF (.NOT.LSTAT) GO TO 370

C We want to exploit zeros in the front.

C We want to swap columns which have zero entries
C in each of the pivotal rows to the front of the matrix

C LPRE is number of cols with zero entry in each of the pivot rows.
C The number of columns which are not fully summed is LFR - KR
C Loop over these columns, in reverse order.
      DO 240 L = LFR - KR,1,-1
          IF (L.EQ.LPRE) GO TO 250
C Is column L of the frontal matrix to be swapped?
          DO 210 LK = 1,CHFRNT
             KPIVRO = KFR + LK
C Jump to next col if col L has at least one nonzero entry in
C the pivot rows
             IF (FA(KPIVRO,L).NE.ZERO) GO TO 240
  210     CONTINUE
C All entries in col L in pivotal rows are zero
          LPRE = LPRE + 1
C Swap col L with the first  col
C which does NOT have a zero in all the pivot rows.
            L1 = LPRE
            DO 220 LL = L1,L - 1
               DO 215 LK = 1,CHFRNT
                  KPIVRO = KFR + LK
C Swap col L with col LL if col LL has at least one nonzero entry
                  IF (FA(KPIVRO,LL).NE.ZERO) GO TO 225
  215          CONTINUE
C All entries in col LL are zero but no swap needed (already in place)
               LPRE = LPRE + 1
 220       CONTINUE
C All remaining columns have zero entries
C in each of the pivotal rows so no more swapping
            GO TO 250
  225       CONTINUE

C Swap column LL with column L
             J = LHED(L)
             JJ = LHED(LL)
C Here we are allowing all possible combinations of
C column L,LL being fully/non-fully summed
C (only needed when we are doing one elimination at a time)
            IF (KPVLNK(L).LT.0) THEN
C Column L is not fully summed
               IF (KPVLNK(LL).LT.0) THEN
C Column LL is not fully summed either ... swap
                  I = KPVLNK(L)
                  KPVLNK(L) = KPVLNK(LL)
                  KPVLNK(LL) = I
C Change LDEST according to permutation.
                  LDEST(-LAST(J)) = LL
                  LDEST(-LAST(JJ)) = L
               ELSE
C Col. LL is fully summed
                  I = KPVLNK(LL)
                  KPVLNK(LL) = KPVLNK(L)
                  KPVLNK(L) = I
                  KPIV(I) = L
                  LDEST(-LAST(J)) = LL
               END IF
            ELSE
C Col L is fully summed
               IF (KPVLNK(LL).LT.0) THEN
C Column LL is not fully summed
                  I = KPVLNK(L)
                  KPVLNK(L) = KPVLNK(LL)
                  KPVLNK(LL) = I
                  KPIV(I) = LL
                  LDEST(-LAST(JJ)) = L
               ELSE
C Col. LL is also fully summed
                  I = KPVLNK(L)
                  LCOL = KPVLNK(LL)
                  KPVLNK(L) = LCOL
                  KPIV(LCOL) = L
                  KPVLNK(LL) = I
                  KPIV(I) = LL
               END IF
            END IF

C This is all we need if we are working only with columns
C which are not fully summed
C           I = KPVLNK(L)
C           KPVLNK(L) = KPVLNK(LL)
C           KPVLNK(LL) = I
C Change LDEST according to permutation.
C           LDEST(-LAST(J)) = LL
C           LDEST(-LAST(JJ)) = L

            LHED(L) = JJ
            LHED(LL) = J

            INFO(2) = -INFO(2)

            DO 230 K = 1,KFRNT
               SWAP = FA(K,L)
               FA(K,L) = FA(K,LL)
               FA(K,LL) = SWAP
  230       CONTINUE

  240 CONTINUE
  250 CONTINUE

C We now want to swap rows which have zero entries
C in each of the pivotal columns to the front of the matrix

C The number of rows which are not fully summed is KFR-KR
C Loop over these rows, in reverse order
      DO 310 K = KFR - KR,1,-1
         IF (K.EQ.KPRE) GO TO 320
C Is row K to be swapped?
         DO 260 LK = 1,CHFRNT
             LPIVCO = LFR + LK
C Jump to next row if row K has at least one nonzero entry
             IF (FA(K,LPIVCO).NE.ZERO) GO TO 310
  260    CONTINUE
C All entries in row K in pivotal columns are zero
         KPRE = KPRE + 1
C Swap row K with the first  row
C which does NOT have a zero in all the pivot columns.
         L1 = KPRE
         DO 280 KK = L1,K - 1
            DO 270 LK = 1,CHFRNT
               LPIVCO = LFR + LK
C Swap row K with row KK if row KK has at least one nonzero entry
               IF (FA(KK,LPIVCO).NE.ZERO) GO TO 290
  270       CONTINUE
C All entries in row KK are zero but no swap needed (already in place)
            KPRE = KPRE + 1
  280    CONTINUE
C All remaining rows have zero entries
C in each of the pivotal columns so no more swapping
         GO TO 320
  290    CONTINUE

C Swap  row KK with row K
         J = KHED(K)
         JJ = KHED(KK)
C Change KDEST according to permutation.
         IF (NMAXE.GT.1) THEN
            IF (LAST(J).LT.0) KDEST(-LAST(J)) = KK
            IF (LAST(JJ).LT.0) KDEST(-LAST(JJ)) = K
         END IF

         KHED(K) = JJ
         KHED(KK) = J

         INFO(2) = -INFO(2)

         DO 300 L = 1,LFRNT
            SWAP = FA(K,L)
            FA(K,L) = FA(KK,L)
            FA(KK,L) = SWAP
  300    CONTINUE

C Swap right-hand sides in rows K and KK.
         DO 305 LK = 1,NRHS
            SWAP = FRHS(K,LK)
            FRHS(K,LK) = FRHS(KK,LK)
            FRHS(KK,LK) = SWAP
  305    CONTINUE

  310 CONTINUE
  320 CONTINUE

  370 CONTINUE

      IF (CHFRNT.GT.0) THEN
C Update FRHS and FA.
         IF (CHFRNT.EQ.1) THEN
            PIVINV = -ONE/FA(KFRNT,LFRNT)
C Scale by pivot and then update rest of FRHS.
            DO 380 L = 1,NRHS
               FRHS(KFRNT,L) = FRHS(KFRNT,L)*PIVINV
  380       CONTINUE

            IF (KFR.NE.KPRE .AND. NRHS.NE.0) CALL DGER(KFR-KPRE,NRHS,
     +          ONE,FA(KPRE+1,LFRNT),1,FRHS(KFRNT,1),MFRONT,
     +          FRHS(KPRE+1,1),MFRONT)

            IF (LFR.NE.KR+LPRE) THEN
C Scale pivot row by pivot (KFR+1 = KFRNT)
               DO 390 K = 1,LFR - KR - LPRE
                  FA(KFRNT,LPRE+K) = FA(KFRNT,LPRE+K)*PIVINV
  390          CONTINUE
C Use Level 2 to update rest of FA.
               IF (KFR.NE.KPRE) CALL DGER(KFR-KPRE,LFR-KR-LPRE,ONE,
     +                                    FA(KPRE+1,LFR+1),1,
     +                                    FA(KFR+1,LPRE+1),MFRONT,
     +                                    FA(KPRE+1,LPRE+1),MFRONT)
            END IF

         ELSE
C Use Level 3 BLAS
            CALL DTRSM('L','U','N','N',CHFRNT,NRHS,-ONE,FA(KFR+1,LFR+1),
     +                 MFRONT,FRHS(KFR+1,1),MFRONT)
            IF (KFR.NE.KPRE .AND. NRHS.NE.0) CALL DGEMM('N','N',
     +          KFR-KPRE,NRHS,CHFRNT,ONE,FA(KPRE+1,LFR+1),MFRONT,
     +          FRHS(KFR+1,1),MFRONT,ONE,FRHS(1+KPRE,1),MFRONT)
            IF (LFR.NE.KR+LPRE) THEN
C Update rows KFR+1, KFR+2,...,KFRNT, columns LPRE+1,...,LFR-KR of FA.
               CALL DTRSM('L','U','N','N',CHFRNT,LFR-KR-LPRE,-ONE,
     +                    FA(KFR+1,LFR+1),MFRONT,FA(KFR+1,LPRE+1),
     +                    MFRONT)
C Update rest of FA.
               IF (KFR.NE.KPRE) CALL DGEMM('N','N',KFR-KPRE,LFR-KR-LPRE,
     +                                     CHFRNT,ONE,FA(KPRE+1,LFR+1),
     +                                     MFRONT,FA(KFR+1,LPRE+1),
     +                                     MFRONT,ONE,FA(KPRE+1,LPRE+1),
     +                                     MFRONT)
            END IF
         END IF
C Write to integer buffer.
         CALL MA42HD(1,KFRNT-KPRE,LFRNT-LPRE,CHFRNT,IBUFR,LIBUFR,
     +               KHED(1+KPRE),LHED(1+LPRE),MFRONT-KPRE,NFRONT-LPRE,
     +               LP,IFILE(3),IREC(3),ISIZE(3),MKEY(3),NUMBLK(3),
     +               IELL,NELL,INFO)
         IF (INFO(1).LT.0) GO TO 420
C Write to U buffer.
         CALL MA42GD(1,1,FA(KPRE+1,LPRE+1),FRHS(1+KPRE,1),NRHS,
     +               KFRNT-KPRE,LFRNT-LPRE,CHFRNT,BUFRU,LUB,MFRONT,
     +               LP,IFILE(1),IREC(1),ISIZE(1),MKEY(1),
     +               NUMBLK(1),IELL,NELL,INFO)
         IF (INFO(1).LT.0) GO TO 420
C Write to L buffer.
         CALL MA42GD(2,1,FA(KPRE+1,LPRE+1),FRHS(1+KPRE,1),NRHS,
     +               KFRNT-KPRE,LFRNT-LPRE,CHFRNT,BUFRL,LLB,MFRONT,
     +               LP,IFILE(2),IREC(2),ISIZE(2),MKEY(2),
     +               NUMBLK(2),IELL,NELL,INFO)
         IF (INFO(1).LT.0) GO TO 420
C Update OPS
         DO 400 J = 1,CHFRNT
            OPS = OPS + DBLE((KFRNT-J-KPRE)* (LFRNT-J-LPRE)*2)
  400    CONTINUE
C Update front sizes.
         LFRNT = LFR
         KFRNT = KFR
      END IF

      IF (NSIZE.EQ.1) KR = KRTEMP

      GO TO 420
C
C Error return
  410 INFO(1) = -14
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9010) INFO(1),IELL
         WRITE (LP,FMT=9000)
      END IF
C
  420 RETURN
 9000 FORMAT (7X,'Matrix found to be singular')
 9010 FORMAT (' ***** Error return from MA42B/BD *****  INFO(1) = ',I3,
     +       /' ***** Error discovered after input of elt/eqn *****',I8)
 9020 FORMAT (' ***** Warning from MA42B/BD *****  INFO(1) = ',I3,
     +       /7X,'after input of elt/equ ',I8)
      END
C**********************************************************************
      SUBROUTINE MA42OD(IELL,AVAR,LRHS,NMAXE,RHS,NRHS,IVAR,NDF,LAST,
     +                  NVAR,MFRONT,NFRONT,LHED,KHED,KR,KPIV,KPVLNK,
     +                  LDEST,KDEST,ISTATC,LFREE,FRHS,BUFRL,LLB,BUFRU,
     +                  LUB,IBUFR,LIBUFR,MVAR,KFRNT,LFRNT,FA,IFILE,IREC,
     +                  ISIZE,MKEY,NUMBLK,DET,OPS,NELL,CNTL,ICNTL,INFO)
C
C  This routine performs static condensations
C
C     Argument list. * indicates the argument is changed by the routine.
C
C     IELL      Integer variable. Holds current elt/equ number.
C     AVAR      Real (DP) array of dimensions NMAXE by NVAR.
C               contains contributions from current elt/eqn.
C     LRHS      Integer variable. Second dimension of FRHS array.
C     NMAXE     Integer variable. Leading dimension of element or
C               equation arrays. NMAXE=1 for equation entry.
C     RHS       Real (DP) array of dimensions NMAXE by LRHS.
C               Contains right hand sides for elt/eqn being input.
C     NRHS      Integer variable.  Number of right hand sides.
C     IVAR      Integer array of length NVAR.
C               Contains indices of variables in current elt/eqn.
C               It is used locally as a work array by the routine
C               to facilitate access to the position of an incoming
C               variable in the front but it is reset and returned
C               unchanged to MA42F/FD.
C     NDF       Integer variable. Total number of variables.
C    *LAST      Integer array of length NDF.
C               For a variable not yet in the front, on entry
C               LAST(I) is the elt/eqn in which variable I
C               appears for the last time. For a non-fully summed
C               variable in the front, JVAR say, -LAST(JVAR) is the
C               position in KDEST/LDEST of the information for
C               this variable.
C     NVAR      Integer variable. Number of variables in current
C               elt/equ.
C     MFRONT    Integer variable. Max. number of rows in frontal matrix.
C     NFRONT    Integer variable. Max. number of cols in frontal matrix.
C    *LHED      Integer array of dimension NFRONT. Used to hold the
C               indices of variables corresponding to columns in the
C               front.
C    *KHED      Integer array of dimension MFRONT. Used to hold the
C               indices of variables corresponding to rows in the
C               front.
C    *KPIV      Integer array of dimension NFRONT.  Indicates which
C               columns of FA are fully summed.
C    *KR        Integer variable. On exit holds number of fully
C               assembled variables. (Only changed if the number
C               of static condensations performed NUMPIV is less
C               than ISTATC).
C    *KPVLNK    Integer array of dimension NFRONT. For each fully
C               summed column in the front, KPVLNK indicates its
C               position in array KPIV.
C    *LDEST     Integer array of length NFRONT.  For each non-fully
C               summed variable in the front (global name JVAR say)
C               LDEST(-LAST(JVAR)) is its column position in the
C               frontal matrix.
C    *KDEST     Integer array of length MFRONT.  For each non-fully
C               summed variable in the front (global name JVAR say)
C               KDEST(-LAST(JVAR)) is its row position in the
C               frontal matrix.
C     ISTATC    Integer variable. Number of static condensation
C               candidates in current elt/equ
C    *LFREE     Integer variable. Points to first free location in
C               LDEST, KDEST
C    *DET       Real (DP) variable. Accumulates natural
C               logarithm of determinant of matrix.
C    *TMPRHS    Real (DP) array of length LRHS. Used as a working
C               vector to hold right hand sides when forward elimination
C               is being performed on them.
C    *BUFRL     Real (DP) array of length LLB. In-core buffer
C               for factors of L.
C     LLB       Integer variable. Length of array BUFRL.
C    *BUFRU     Real (DP) array of length LUB. In-core buffer
C               for factors of U.
C     LUB       Integer variable. Length of array BUFRU.
C    *IBUFR     Integer array of length LIBUFR. In-core buffer
C               for integer information on factors.
C     LIBUFR    Integer variable. Length of array IBUFR.
C     MVAR      Integer variable. Set to NVAR for elt entry and to 1
C               for equ entry.
C    *LFRNT     Integer variable. On exit holds number of columns in
C               front after eliminations.
C    *KFRNT     Integer variable. On exit holds number of rows in front
C               after eliminations.
C     FA        Real (DP) array with dimensions MFRONT,NFRONT.
C               Holds the current frontal matrix. Used by the routine
C               but unchanged on exit.
C     IFILE     Integer array of length 3. Stream numbers for direct
C               access data sets.
C    *IREC      Integer array of length 3. Pointers to first free space
C               in buffers.
C     ISIZE     Integer array of length 3. Length of buffers.
C    *MKEY      Integer array of length 3. Number of records written
C               to direct access data sets.
C     NUMBLK    Integer array of length 3. Number of records in
C               direct access data sets.
C      CNTL     Real (DP) array of length 2. Matrix singular if
C               largest element in column less than CNTL(1).
C               An element of the frontal matrix will only be considered
C               suitable for use as a pivot if it is at least as large
C               as CNTL(2) times the largest element in its column.
C     ICNTL     Integer array of length 8.
C    *INFO      Integer array of length 23. INFO(1) is error
C               flag. Possible fatal errors
C               returned by this routine are INFO(1) = -14,-26.
C               INFO = -14 if element matrix is singular
C               INFO = -26 if there is a failure in direct access write.
C
C     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D0,ZERO=0.0D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION DET,OPS
      INTEGER IELL,ISTATC,KFRNT,KR,LFREE,LFRNT,LIBUFR,LLB,LRHS,LUB,
     +        MFRONT,MVAR,NDF,NELL,NFRONT,NMAXE,NRHS,NVAR
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION AVAR(NMAXE,NVAR),BUFRL(LLB),BUFRU(LUB),CNTL(2),
     +                 FA(MFRONT,NFRONT),FRHS(MFRONT,LRHS),
     +                 RHS(NMAXE,LRHS)
      INTEGER IBUFR(LIBUFR),ICNTL(8),IFILE(3),INFO(23),IREC(3),ISIZE(3),
     +        IVAR(NVAR),KDEST(MFRONT),KHED(MFRONT),KPIV(NFRONT),
     +        KPVLNK(NFRONT),LAST(NDF),LDEST(NFRONT),LHED(NFRONT),
     +        MKEY(3),NUMBLK(3)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION PIVINV,PIVOT,RMAX,SMAX,SWAP
      INTEGER I,IMAX,IPIV,IPOS,IRHS,ISWAP,J,J1,JCOL,JJ,JSTATC,JVAR,K,
     +        KDST,KFR,KFRNEW,KFROLD,KPOS,KR1,L,LC1,LDST,LFR,LFRNEW,
     +        LFROLD,LK,LP,MFR,MOVE,MP,NFREE,NUMPIV,NUVAR,NUVRP1
      LOGICAL LDET
C     ..
C     .. External Functions ..
      INTEGER IDAMAX
      EXTERNAL IDAMAX
C     ..
C     .. External Subroutines ..
      EXTERNAL DGER,MA42GD,MA42HD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,LOG,MAX,MIN,DBLE
C     ..
C Immediate return if ISTATC=0.
      IF (ISTATC.EQ.0) GO TO 490
      LP = ICNTL(1)
      MP = ICNTL(2)
C NUVAR is number of non-fully summed rows/columns in
C incoming elt/equ.
      NUVAR = NVAR - ISTATC
      JCOL = NVAR
      J1 = 1
C Sort static condensation variable(s) to end of elt/equ.
      DO 60 MOVE = 1,ISTATC
C JCOL is column currently being tested for static condensation.
         MFR = IVAR(JCOL)
         IF (LAST(MFR).NE.IELL) THEN
C Not in static condensation set.
            LDET = LAST(MFR) .GT. 0
C Look for a column to swap.
            DO 10 J = J1,JCOL - 1
               MFR = IVAR(J)
               IF (LAST(MFR).EQ.IELL) GO TO 20
   10       CONTINUE
C Swap columns J and JCOL.
   20       ISWAP = IVAR(J)
            IVAR(J) = IVAR(JCOL)
            IVAR(JCOL) = ISWAP
            J1 = J + 1
            IF (LDET) INFO(2) = -INFO(2)
            DO 30 I = 1,MVAR
               SWAP = AVAR(I,J)
               AVAR(I,J) = AVAR(I,JCOL)
               AVAR(I,JCOL) = SWAP
   30       CONTINUE
            IF (NMAXE.GT.1) THEN
C Element entry. Swap rows J and JCOL.
               IF (LDET) INFO(2) = -INFO(2)
               DO 40 I = 1,NVAR
                  SWAP = AVAR(J,I)
                  AVAR(J,I) = AVAR(JCOL,I)
                  AVAR(JCOL,I) = SWAP
   40          CONTINUE
C Swap right hand sides.
               DO 50 IRHS = 1,NRHS
                  SWAP = RHS(J,IRHS)
                  RHS(J,IRHS) = RHS(JCOL,IRHS)
                  RHS(JCOL,IRHS) = SWAP
   50          CONTINUE
            END IF
         END IF
         JCOL = JCOL - 1
   60 CONTINUE
C
C Place incoming variables  which are not static condensation
C candidates and are not already in front in last rows and
C columns of front matrix.
      DO 70 LK = 1,NUVAR
         MFR = IVAR(LK)
         IF (LAST(MFR).GE.0) THEN
C Insert incoming variable in the front.
C Set KDEST/LDEST and KHED/LHED.
            LFRNT = LFRNT + 1
            NFREE = LDEST(LFREE)
            LDEST(LFREE) = LFRNT
            LHED(LFRNT) = MFR
            IF (NMAXE.GT.1) THEN
               KFRNT = KFRNT + 1
               KHED(KFRNT) = MFR
               KDEST(LFREE) = KFRNT
            END IF
C Save value of LAST and set flag in LAST to indicate variable is
C now in front.
            KPVLNK(LFRNT) = -LAST(MFR)
            LAST(MFR) = -LFREE
            LFREE = NFREE
         END IF
   70 CONTINUE
C
C NUMPIV is the actual number of static condensations which are
C performed on incoming elt/eqn (NUMPIV may be less than ISTATC
C because of stability considerations)
      IF (NMAXE.EQ.1) THEN
C Equation entry. Must have ISTATC = 1.
C Set KHED/LHED for static condensation variable.
         MFR = IVAR(NVAR)
         KFRNT = KFRNT + 1
         LFRNT = LFRNT + 1
         KHED(KFRNT) = IELL
         LHED(LFRNT) = MFR
C Perform elimination within equ using entry (1,NVAR) as pivot.
C (Divide entries in equation and corresponding right-hand side
C by pivot). Check pivot is greater than CNTL(1) (otherwise matrix
C is considered to be singular).
         PIVOT = AVAR(1,NVAR)
         IF (ABS(PIVOT).LE.CNTL(1)) THEN
C Matrix is singular.
            NUMPIV = 0
            INFO(2) = 0
            DET = ZERO
C Jump to error return if ICNTL(8) is 0.
            IF (ICNTL(8).EQ.0) GO TO 460
C Otherwise set INFO(1) = 1 and issue a warning.
            IF (INFO(1).EQ.0) THEN
               INFO(1) = 1
               IF (MP.GT.0) THEN
                  WRITE (MP,FMT=9020) INFO(1),IELL
                  WRITE (MP,FMT=9010)
               END IF
            END IF
            GO TO 300
         END IF
         DET = DET + LOG(ABS(PIVOT))
         IF (PIVOT.LT.ZERO) INFO(2) = -INFO(2)
         AVAR(1,NVAR) = -PIVOT
         PIVINV = ONE/PIVOT
         DO 80 IRHS = 1,NRHS
            RHS(1,IRHS) = RHS(1,IRHS)*PIVINV
   80    CONTINUE
         DO 90 J = 1,NUVAR
            AVAR(1,J) = AVAR(1,J)*PIVINV
   90    CONTINUE
C
      ELSE
C Element entry.
C Set KHED/LHED for static condensation rows/columns.
         NUVRP1 = NUVAR + 1
         DO 100 I = NUVRP1,NVAR
            MFR = IVAR(I)
            LFRNT = LFRNT + 1
            LHED(LFRNT) = MFR
            KFRNT = KFRNT + 1
            KHED(KFRNT) = MFR
  100    CONTINUE
C Perform static condensation within element matrix itself.
         JCOL = NVAR
C Jump if we are only allowing diagonal pivots
         IF (ABS(ICNTL(7)).EQ.1001) GO TO 200
         DO 190 JSTATC = 1,ISTATC
C Try to find pivot.  We scan columns from the back forwards to
C minimize swap operations.
            DO 110 J = JCOL,NUVRP1,-1
               INFO(13) = INFO(13) + 1
               INFO(15) = INFO(15) + JCOL
C Find max entry in column J.
               IMAX = IDAMAX(JCOL,AVAR(1,J),1)
               RMAX = ABS(AVAR(IMAX,J))
               IF (RMAX.LE.CNTL(1)) THEN
C Matrix is structurally singular.
                  INFO(2) = 0
                  DET = ZERO
C Jump to error return if ICNTL(8) is 0.
                  IF (ICNTL(8).EQ.0) THEN
                     NUMPIV = JSTATC - 1
                     GO TO 460
                  ELSE IF (INFO(1).EQ.0) THEN
C Otherwise set INFO(1) = 1 and issue a warning.
                     INFO(1) = 1
                     IF (MP.GT.0) THEN
                        WRITE (MP,FMT=9020) INFO(1),IELL
                        WRITE (MP,FMT=9010)
                     END IF
                  END IF
                  IF (J.EQ.NUVRP1) THEN
C No more static condensations can be performed without violating
C threshold criterion.
                     NUMPIV = JSTATC - 1
                     GO TO 300
                  END IF
                  GO TO 110
               END IF
               INFO(14) = INFO(14) + 1
C Jump if row JCOL is suitable for use as pivot. (Pivot is in
C row JCOL, column J)
               IF (ABS(AVAR(JCOL,J)).GE.CNTL(2)*RMAX) GO TO 150
C Jump if maximum entry is fully summed (pivot found).
C (Pivot is in row IMAX, column J)
               IF (IMAX.GE.NUVRP1) GO TO 120
C Find max entry in fully summed part.
               INFO(14) = INFO(14) + JCOL - NUVAR
               IMAX = IDAMAX(JCOL-NUVAR,AVAR(NUVRP1,J),1) + NUVAR
               SMAX = ABS(AVAR(IMAX,J))
C Jump if maximum fully-summed entry is big enough.
C (Pivot is in row IMAX, column J)
               IF (SMAX.GE.CNTL(2)*RMAX) GO TO 120
C Otherwise, if J = NUVRP1 no more static condensations
C can be performed without violating threshold criterion.
               IF (J.EQ.NUVRP1) THEN
                  NUMPIV = JSTATC - 1
                  GO TO 300
               END IF
  110       CONTINUE
C Permute row IMAX to the end.
  120       DO 130 JJ = 1,NVAR
               SWAP = AVAR(JCOL,JJ)
               AVAR(JCOL,JJ) = AVAR(IMAX,JJ)
               AVAR(IMAX,JJ) = SWAP
  130       CONTINUE
            INFO(2) = -INFO(2)
C Adjust KHED array.  Pivot row will be in position
C KFRNT+JCOL-NVAR of frontal matrix and has just been moved from
C position KFRNT+IMAX-NVAR.
            KFRNEW = KFRNT + JCOL - NVAR
            KFROLD = KFRNT + IMAX - NVAR
            ISWAP = KHED(KFRNEW)
            KHED(KFRNEW) = KHED(KFROLD)
            KHED(KFROLD) = ISWAP
C Permute right hand side array.
            DO 140 IRHS = 1,NRHS
               SWAP = RHS(JCOL,IRHS)
               RHS(JCOL,IRHS) = RHS(IMAX,IRHS)
               RHS(IMAX,IRHS) = SWAP
  140       CONTINUE
  150       IF (J.NE.JCOL) THEN
C Swap columns J and JCOL.
               DO 160 I = 1,NVAR
                  SWAP = AVAR(I,J)
                  AVAR(I,J) = AVAR(I,JCOL)
                  AVAR(I,JCOL) = SWAP
  160          CONTINUE
               INFO(2) = -INFO(2)
               LFRNEW = LFRNT + JCOL - NVAR
               LFROLD = LFRNT + J - NVAR
               ISWAP = LHED(LFRNEW)
               LHED(LFRNEW) = LHED(LFROLD)
               LHED(LFROLD) = ISWAP
            END IF
C Perform elimination within element.
            PIVOT = AVAR(JCOL,JCOL)
            DET = DET + LOG(ABS(PIVOT))
            IF (PIVOT.LT.ZERO) INFO(2) = -INFO(2)
            AVAR(JCOL,JCOL) = -PIVOT
C Change sign of pivotal column and scale pivotal row by pivot..
            PIVINV = ONE/PIVOT
            DO 170 I = 1,JCOL - 1
               AVAR(I,JCOL) = -AVAR(I,JCOL)
               AVAR(JCOL,I) = AVAR(JCOL,I)*PIVINV
  170       CONTINUE
C Scale right hand sides by pivot.
            DO 180 IRHS = 1,NRHS
               RHS(JCOL,IRHS) = RHS(JCOL,IRHS)*PIVINV
  180       CONTINUE
C Update rest of AVAR and RHS using BLAS
            IF (JCOL.NE.1) CALL DGER(JCOL-1,JCOL-1,ONE,AVAR(1,JCOL),1,
     +                               AVAR(JCOL,1),NMAXE,AVAR(1,1),NMAXE)
            IF (JCOL.NE.1 .AND. NRHS.NE.0) CALL DGER(JCOL-1,NRHS,ONE,
     +          AVAR(1,JCOL),1,RHS(JCOL,1),NMAXE,RHS(1,1),NMAXE)
            JCOL = JCOL - 1
  190    CONTINUE
         GO TO 290

  200    CONTINUE
C Search for diagonal pivots only
         DO 280 JSTATC = 1,ISTATC
            DO 210 J = JCOL,NUVRP1,-1
               INFO(13) = INFO(13) + 1
               INFO(15) = INFO(15) + JCOL
C Find max entry in column J.
               IMAX = IDAMAX(JCOL,AVAR(1,J),1)
               RMAX = ABS(AVAR(IMAX,J))
               IF (RMAX.LE.CNTL(1)) THEN
C Matrix is structurally singular.
                  INFO(2) = 0
                  DET = ZERO
C Jump to error return if ICNTL(8) is 0.
                  IF (ICNTL(8).EQ.0) THEN
                     NUMPIV = JSTATC - 1
                     GO TO 460
                  ELSE IF (INFO(1).EQ.0) THEN
C Otherwise set INFO(1) = 1 and issue a warning.
                     INFO(1) = 1
                     IF (MP.GT.0) THEN
                        WRITE (MP,FMT=9020) INFO(1),IELL
                        WRITE (MP,FMT=9010)
                     END IF
                  END IF
                  IF (J.EQ.NUVRP1) THEN
C No more static condensations can be performed without violating
C threshold criterion.
                     NUMPIV = JSTATC - 1
                     GO TO 300
                  END IF
                  GO TO 210
               END IF
               INFO(14) = INFO(14) + 1
C Check to see if the diagonal entry (J,J) can be used as a pivot.
               IF (ABS(AVAR(J,J)).GE.CNTL(2)*RMAX) GO TO 220
C Otherwise, if J = NUVRP1 no more static condensations
C can be performed without violating threshold criterion.
               IF (J.EQ.NUVRP1) THEN
                  NUMPIV = JSTATC - 1
                  GO TO 300
               END IF
  210       CONTINUE
C Permute row J to the end.
  220       IF (J.NE.JCOL) THEN
               DO 230 JJ = 1,NVAR
                  SWAP = AVAR(JCOL,JJ)
                  AVAR(JCOL,JJ) = AVAR(J,JJ)
                  AVAR(J,JJ) = SWAP
  230          CONTINUE
C Adjust KHED array.  Pivot row will be in position
C KFRNT+JCOL-NVAR of frontal matrix and has just been moved from
C position KFRNT+J-NVAR.
               KFRNEW = KFRNT + JCOL - NVAR
               KFROLD = KFRNT + J - NVAR
               ISWAP = KHED(KFRNEW)
               KHED(KFRNEW) = KHED(KFROLD)
               KHED(KFROLD) = ISWAP
C Permute right hand side array.
               DO 240 IRHS = 1,NRHS
                  SWAP = RHS(JCOL,IRHS)
                  RHS(JCOL,IRHS) = RHS(J,IRHS)
                  RHS(J,IRHS) = SWAP
  240          CONTINUE
C Swap columns J and JCOL.
               DO 250 I = 1,NVAR
                  SWAP = AVAR(I,J)
                  AVAR(I,J) = AVAR(I,JCOL)
                  AVAR(I,JCOL) = SWAP
  250          CONTINUE
               LFRNEW = LFRNT + JCOL - NVAR
               LFROLD = LFRNT + J - NVAR
               ISWAP = LHED(LFRNEW)
               LHED(LFRNEW) = LHED(LFROLD)
               LHED(LFROLD) = ISWAP
            END IF
C Perform elimination within element.
            PIVOT = AVAR(JCOL,JCOL)
            DET = DET + LOG(ABS(PIVOT))
            IF (PIVOT.LT.ZERO) INFO(2) = -INFO(2)
            AVAR(JCOL,JCOL) = -PIVOT
C Change sign of pivotal column and scale pivotal row by pivot.
            PIVINV = ONE/PIVOT
            DO 260 I = 1,JCOL - 1
               AVAR(I,JCOL) = -AVAR(I,JCOL)
               AVAR(JCOL,I) = AVAR(JCOL,I)*PIVINV
  260       CONTINUE
C Scale right hand sides by pivot.
            DO 270 IRHS = 1,NRHS
               RHS(JCOL,IRHS) = RHS(JCOL,IRHS)*PIVINV
  270       CONTINUE
C Update rest of AVAR and RHS using BLAS
            IF (JCOL.NE.1) CALL DGER(JCOL-1,JCOL-1,ONE,AVAR(1,JCOL),1,
     +                               AVAR(JCOL,1),NMAXE,AVAR(1,1),NMAXE)
            IF (JCOL.NE.1 .AND. NRHS.NE.0) CALL DGER(JCOL-1,NRHS,ONE,
     +          AVAR(1,JCOL),1,RHS(JCOL,1),NMAXE,RHS(1,1),NMAXE)
            JCOL = JCOL - 1
  280    CONTINUE

  290    CONTINUE

      END IF
C
      NUMPIV = ISTATC
C Static condensations have been performed.
C Update INFO(17) and INFO(18).
  300 INFO(17) = INFO(17) + NUMPIV
      INFO(18) = INFO(18) + ISTATC
C Put static condensation rows/cols in FA (and FRHS) to be written
C out to the buffers.
C First zero out any so far unused portion of the arrays.
      KR1 = INFO(8) + 1
      LC1 = INFO(9) + 1
      IF (KR1.LE.KFRNT) THEN
         INFO(8) = KFRNT
         LFR = MAX(INFO(9),LFRNT)
         DO 330 K = KR1,KFRNT
            DO 310 L = 1,LFR
               FA(K,L) = ZERO
  310       CONTINUE
            DO 320 J = 1,NRHS
               FRHS(K,J) = ZERO
  320       CONTINUE
  330    CONTINUE
      END IF
      IF (LC1.LE.LFRNT) INFO(9) = LFRNT
      DO 350 L = LC1,LFRNT
         DO 340 K = 1,MIN(MFRONT,KR1)
            FA(K,L) = ZERO
  340    CONTINUE
  350 CONTINUE

      KFR = KFRNT
      LFR = LFRNT

C For equation entry only need to write to the front if NUMPIV=0.
      IF (NMAXE.EQ.1 .AND. NUMPIV.EQ.1) GO TO 455
C
C Use IVAR to temporarily hold the column position of the variable
C in the front
      DO 360 LK = 1,NUVAR
         MFR = IVAR(LK)
         IPOS = -LAST(MFR)
C LDST is the column position of the variable in the front.
         LDST = LDEST(IPOS)
         IVAR(LK) = LDST
  360 CONTINUE
C
C For equation entry only need to write to the front if NUMPIV=0.
      IF (NMAXE.EQ.1) THEN
         IF (NUMPIV.EQ.0) THEN
C Set pivot row in FA and FRHS.
          DO 370 J = 1,NVAR - 1
            LDST = IVAR(J)
            FA(KFR,LDST) = AVAR(1,J)
  370     CONTINUE
          FA(KFR,LFR) = AVAR(1,NVAR)
          DO 380 J = 1,NRHS
            FRHS(KFR,J) = RHS(1,J)
  380     CONTINUE
         END IF
      ELSE
C Loop over the static condensation variables.
C Place static condensation rows and columns in last rows
C and columns of FA.
         JCOL = NVAR
         DO 440 IPIV = 1,ISTATC
C Set row KFR
            DO 390 J = 1,NUVAR
               LDST = IVAR(J)
               FA(KFR,LDST) = AVAR(JCOL,J)
  390       CONTINUE
            DO 400 J = NUVAR + 1,JCOL
               FA(KFR,LFR+J-JCOL) = AVAR(JCOL,J)
  400       CONTINUE
C Set column LFR
            DO 410 J = 1,NUVAR
               LDST = IVAR(J)
               JVAR = LHED(LDST)
               KPOS = -LAST(JVAR)
               KDST = KDEST(KPOS)
               FA(KDST,LFR) = AVAR(J,JCOL)
  410       CONTINUE
            DO 420 J = NUVAR + 1,JCOL
               FA(KFR+J-JCOL,LFR) = AVAR(J,JCOL)
  420       CONTINUE
C Place right hand side in row KFR of FRHS.
            DO 430 J = 1,NRHS
               FRHS(KFR,J) = RHS(JCOL,J)
  430       CONTINUE
            KFR = KFR - 1
            LFR = LFR - 1
            JCOL = JCOL - 1
  440    CONTINUE
      END IF

C Restore IVAR
      DO 450 LK = 1,NUVAR
         LDST = ABS(IVAR(LK))
         JVAR = LHED(LDST)
         IVAR(LK) = JVAR
  450 CONTINUE

  455 CONTINUE
C Now have the pivot rows/columns in the last rows/columns
C of the frontal matrix. If the number NUMPIV of pivot rows/columns
C is odd AND the remaining frontsize is odd,
C need to change the sign of the determinant
      IF ((NUMPIV/2)*2.NE.NUMPIV) THEN
         KFR = KFRNT - NUMPIV
         LFR = LFRNT - NUMPIV
         IF ((KFR/2)*2.NE.KFR) INFO(2) = -INFO(2)
         IF ((LFR/2)*2.NE.LFR) INFO(2) = -INFO(2)
      END IF

C Write to buffers rows and columns corresponding to static
C condensations which have been performed (NUMPIV rows and cols).

C For equation entry, with ISTATC=NUMPIV=1, then we do not need
C to write the static condenstion pivotal row to the frontal
C matrix ... we can write AVAR and IVAR directly to the buffers
C (if we do write to the frontal matrx, we will end up
C storing a lot of zero entries)

C For element entry, if we have used off-diag. pivoting, then
C we need global row and column indices, so we will need both KHED
C and LHED.

      IF (NUMPIV.GT.0) THEN
          IF (NMAXE.EQ.1) THEN
             CALL MA42HD(1,1,NVAR,NUMPIV,IBUFR,LIBUFR,KHED(KFRNT),IVAR,
     +                   1,NVAR,LP,IFILE(3),IREC(3),ISIZE(3),MKEY(3),
     +                   NUMBLK(3),IELL,NELL,INFO)
             IF (INFO(1).EQ.-26) GO TO 490

             CALL MA42GD(1,1,AVAR,RHS,NRHS,MVAR,NVAR,NUMPIV,BUFRU,
     +                   LUB,NMAXE,LP,IFILE(1),IREC(1),ISIZE(1),
     +                   MKEY(1),NUMBLK(1),IELL,NELL,INFO)
             IF (INFO(1).EQ.-26) GO TO 490
             CALL MA42GD(2,1,AVAR,RHS,NRHS,MVAR,NVAR,NUMPIV,BUFRL,
     +                   LLB,NMAXE,LP,IFILE(2),IREC(2),ISIZE(2),
     +                   MKEY(2),NUMBLK(2),IELL,NELL,INFO)
             IF (INFO(1).EQ.-26) GO TO 490
          ELSE
             CALL MA42HD(1,KFRNT,LFRNT,NUMPIV,IBUFR,LIBUFR,KHED,LHED,
     +                   MFRONT,NFRONT,LP,IFILE(3),IREC(3),ISIZE(3),
     +                   MKEY(3),NUMBLK(3),IELL,NELL,INFO)
             IF (INFO(1).EQ.-26) GO TO 490
             CALL MA42GD(1,1,FA,FRHS,NRHS,KFRNT,LFRNT,NUMPIV,BUFRU,
     +                   LUB,MFRONT,LP,IFILE(1),IREC(1),ISIZE(1),
     +                   MKEY(1),NUMBLK(1),IELL,NELL,INFO)
             IF (INFO(1).EQ.-26) GO TO 490
             CALL MA42GD(2,1,FA,FRHS,NRHS,KFRNT,LFRNT,NUMPIV,BUFRL,
     +                   LLB,MFRONT,LP,IFILE(2),IREC(2),ISIZE(2),
     +                   MKEY(2),NUMBLK(2),IELL,NELL,INFO)
             IF (INFO(1).EQ.-26) GO TO 490
C
         END IF
C Update OPS, KFRNT,LFRNT
         DO 456 J = 1,NUMPIV
            OPS = OPS + DBLE((NVAR-J)* (MVAR-J)*2)
  456    CONTINUE
         KFRNT = KFRNT - NUMPIV
         LFRNT = LFRNT - NUMPIV
      END IF

      GO TO 470
  460 INFO(1) = -14
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1),IELL
         WRITE (LP,FMT=9010)
      END IF
      NUMPIV = 0
C
C Put static condenstaion variables which are not being eliminated
C into KPIV/KPVLNK.
  470 DO 480 J = 1,ISTATC - NUMPIV
         KR = KR + 1
         LDST = LFRNT - J + 1
         KPIV(KR) = LDST
         KPVLNK(LDST) = KR
  480 CONTINUE
  490 RETURN
 9000 FORMAT (' ***** Error return from MA42B/BD *****  INFO(1) = ',I3,
     +       /' ***** Error discovered after input of elt/eqn *****',I8)
 9010 FORMAT (7X,'Matrix found to be singular')
 9020 FORMAT (' ***** Warning from MA42B/BD *****  INFO(1) = ',I3,
     +       /7X,'after input of elt/equ ',I8)
      END
! COPYRIGHT (c) 2008 Science and Technology Facilities Council
! 31 July 2008 Version 1.0.0. Original version.

! Written by:  Jonathan Hogg
!
! This module exists to allow the inclusion of mpif.h without requiring
! a fixed source form within any other MPI package than this one.
!
! 31st July 2008 Version 1.0.0 Initial version
      module hsl_mp01
         implicit none
         include 'mpif.h'
      end module
* COPYRIGHT (c) 1988 AEA Technology
* Original date 17 Feb 2005

C 17th February 2005 Version 1.0.0. Replacement for FD05.

      DOUBLE PRECISION FUNCTION FD15AD(T)
C----------------------------------------------------------------
C  Fortran 77 implementation of the Fortran 90 intrinsic
C    functions: EPSILON, TINY, HUGE and RADIX.  Note that
C    the RADIX result is returned as DOUBLE PRECISION.
C
C  The CHARACTER argument specifies the type of result:
C
C   'E'  smallest positive real number: 1.0 + DC(1) > 1.0, i.e.
C          EPSILON(DOUBLE PRECISION)
C   'T'  smallest full precision positive real number, i.e.
C          TINY(DOUBLE PRECISION)
C   'H'  largest finite positive real number, i.e.
C          HUGE(DOUBLE PRECISION)
C   'R'  the base of the floating point arithematic, i.e.
C          RADIX(DOUBLE PRECISION)
C
C    any other value gives a result of zero.
C----------------------------------------------------------------
      CHARACTER T

      IF ( T.EQ.'E' ) THEN
         FD15AD = EPSILON(1.0D0)
      ELSE IF ( T.EQ.'T' ) THEN
         FD15AD = TINY(1.0D0)
      ELSE IF ( T.EQ.'H' ) THEN
         FD15AD = HUGE(1.0D0)
      ELSE IF ( T.EQ.'R' ) THEN
         FD15AD = DBLE(RADIX(1.0D0))
      ELSE
         FD15AD = 0.0D0
      ENDIF
      RETURN
      END
