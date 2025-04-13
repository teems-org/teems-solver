!We have add MYIW1 and MYLMX.
!The code lines related to IW1 and LMX also have been changed accordingly
!We have change the dimension of CGLOB to the maximum number of columns in
!a single block (including border part) times NBLOCK

! COPYRIGHT (c) 2003 Council for the Central Laboratory
!                    of the Research Councils
! Original date 13 March 2003

! 4 feb 2004 : Incorrect solves time sent to root.

! Note: Formerly, we used some pointer assignments (see TEMP, ITEMP etc).
! This is because on DEC Alpha, complier produces very inefficient
! code when pointer components of a derived data type are
! used in a subroutine call and when array sections are used.

! To derive single precision version from double precision version
! (a) Change HSL_MP48_DOUBLE to HSL_MP48_SINGLE
! and HSL_MP48_DATA_DOUBLE to HSL_MP48_DATA_SINGLE
! (b) Change WP = KIND(0.0D0) to WP = KIND(0.0)
!     BUT Be careful: MPI_WTIME is double precision and so
!     need DP = KIND(0.0D0)
! (c) Change routine names from double version to single
!     (eg MP48_MA60AD changed to MP48_MA60A). Also MP48AD to MP48A.
!     Change BLAS routine names.
! (d) Change MPI_DOUBLE_PRECISION to MPI_REAL

! 12th July 2004 Version 1.0.0. Version numbering added.
! 18th Jan. 2005 Version 1.1.0. Bug fix. VALNAM allocatable
!             (not allocatable).
! 28 Aug 2007 Version 1.2.0. Changed lengths of arrays used by MA48
!             to be consistent with HSL 2007 version of MA48.
! 17 July 2008 Version 2.0.0. Bug fixes to allow number of processors
!             to exceed number of submatrices. Corrections also to spec.
!             INCLUDE replaced by USE HSL_MP01 and source form changed.
!             Pointer array components replaced by allocatables.
!             Private communicator added.

      MODULE HSL_MP48_DATA_DOUBLE

      USE HSL_MP01
      IMPLICIT NONE

      INTEGER, PARAMETER :: WP = KIND(0.0D0)
      INTEGER, PARAMETER :: DP = KIND(0.0D0)

      TYPE MP48_DATA_PRIVATE

      INTEGER ERROR
      INTEGER FJOB  ! Initialsed to 0 on JOB =3 call.
!                     Incremented by 1 on each JOB = 4. We
!                     do not want to recall MA48A/AD if more
!                     than one matrix with the same sparsity pattern
!                     is to be factorized ... only call MA48A/AD
!                     if FJOB = 1
      INTEGER LGUARD
      INTEGER LORDER
      INTEGER NBLOCK
      INTEGER NEQ
      INTEGER OLDJOB
      INTEGER, DIMENSION (:), allocatable :: IDUP
      INTEGER, DIMENSION (:), allocatable :: IGLOB
      INTEGER, DIMENSION (:), allocatable :: JPTR
      INTEGER, DIMENSION (:), allocatable :: LPTR
      INTEGER, DIMENSION (:), allocatable :: IPTR
      INTEGER, DIMENSION (:), allocatable :: VPTR
      INTEGER, DIMENSION (:), allocatable :: LFACT
      INTEGER, DIMENSION (:), allocatable :: LIRNF
      INTEGER, DIMENSION (:), allocatable :: NCOLA
      INTEGER, DIMENSION (:), allocatable :: NF
      INTEGER, DIMENSION (:), allocatable :: NP
      INTEGER, DIMENSION (:), allocatable :: NEQSUB
      INTEGER, DIMENSION (:), allocatable :: SUBPTR
      INTEGER, DIMENSION (:), allocatable :: INULL
      INTEGER, DIMENSION (:,:), allocatable :: IMAP
      INTEGER, DIMENSION (:,:), allocatable :: JFVAR
      INTEGER, DIMENSION (:,:), allocatable :: GLOBAL
      INTEGER, DIMENSION (:,:), allocatable :: CGLOB
      INTEGER, DIMENSION (:,:), allocatable :: IP
      INTEGER, DIMENSION (:,:), allocatable :: IQ
      INTEGER, DIMENSION (:,:), allocatable :: IQ_COPY
      INTEGER, DIMENSION (:,:), allocatable :: IQM_COPY

! MA48 variables/arrays
      INTEGER LA48,LKEEP48
      INTEGER,  DIMENSION (:), allocatable :: IRN
      INTEGER,  DIMENSION (:), allocatable :: JCN
      REAL(WP), DIMENSION (:), allocatable :: A
      INTEGER :: ICNTL_48(20)
      REAL(WP) :: CNTL_48(10)
      INTEGER, DIMENSION (:), allocatable :: KEEP48

! MP48_MA60 arrays
      INTEGER,  DIMENSION (:,:), allocatable :: ICNTL_60
      REAL(WP), DIMENSION (:,:), allocatable :: CNTL_60
      INTEGER,  DIMENSION(:,:),  allocatable :: IPTRL
      INTEGER,  DIMENSION(:,:),  allocatable :: IPTRU

! For MPI
      INTEGER :: COMM  = MPI_COMM_NULL  ! Communicator

      END TYPE MP48_DATA_PRIVATE

      TYPE MP48_DATA

! For MPI
      INTEGER COMM   ! Communicator
      INTEGER RANK   ! Rank

      INTEGER JOB              ! Controls action.
! JOB = 1 : initialisation call
! JOB = 2 : pre analyse
! JOB = 3 : analyse using MP48_MA60A/AD
! JOB = 4 : factorize submatrices in parallel
! JOB = 5 : solve for right-hand sides
! JOB = 6 : deallocate all arrays (final call)
! We can combine calls using 23 (2+3), 24 (2+3+4) or 25 (2+3+4+5)
      INTEGER FACT_JOB       ! Controls factorization action.
      INTEGER IOSTAT         ! IOSTAT parameter.
      INTEGER BORDER         ! Number of border columns.
      INTEGER NBLOCK         ! Number of submatrices.
      INTEGER NEQ            ! Number of rows/columns in whole problem.
      INTEGER N              ! N = NEQ.
      INTEGER NE             ! Length of row variable lists.
      INTEGER NE_INTER       ! Number of entries in interface prob.
      INTEGER NULL           ! number of null cols.
      INTEGER NINTER         ! No. of cols. in interface problem
      INTEGER RINTER         ! No. of rows in interface problem
      INTEGER MAXSBCOLS      ! Maximum No. of cols in single block (including border part)
      INTEGER ICNTL(20)      ! Control parameters.
! ICNTL(1)     Error output unit (default 6)
! ICNTL(2)     Warning output unit (default 6)
! ICNTL(3)     Diagnostics output unit (default 6)
! ICNTL(4)     Diagnostics format (default 1)
!              ICNTL(4) <= 0 : suppress diagnostics
!              1 Errors and warnings only.
!              2 As 1 plus additioanl diagnostic printing.
!              3 As 2, but timings of parts of the code (wall clock
!                times in seconds) are also printed.
! ICNTL(5)     is used by MP48 to control the
!              full-matrix processing. It has default value 32.
!              Possible values are:
!              0 Level 1 BLAS used.
!              1 Level 2 BLAS used.
!              2 Level 3 BLAS used by MP48
!              with data%JOB = 4, with block column size ICNTL(5).
!              Level 2 BLAS used by MP48 with data%JOB = 5
!              (solve phase).
! ICNTL(6)     has default value 2. When a real or an integer array
!              that holds data for the factors is too small,
!              it is reallocated with its size
!              increased by the factor ICNTL(6).
! ICNTL(7)     Controls how the user wishes to supply the submatrix
!              data.
! ICNTL(8)     Default value 3. If ICNTL(8) has a positive value,
!              each pivot search in MP48 is limited to a maximum of
!              ICNTL(8) columns.  If ICNTL(8) is set to the value 0,
!              full Markowitz is used to find the best pivot.
!              This requires extra workspace and is usually only a
!              little slower, but can occasionally be
!              very slow.  It may result in reduced fill-in.
! ICNTL(9)     controls the action if duplicate entries within a row are
!              found. If ICNTL(9) = 0 (the default) and duplicates are
!              found, the values of duplicates are summed.
!              If ICNTL(9) nonzero and duplicates are found, the program
!              terminates with an error flag -3.
! ICNTL(10)    is used to control whether or not the user wishes to
!              supply a list of the submatrices that are to be
!              factorized by each processor.  If ICNTL(10) = 0 (the
!              default), the list is generated automatically.
!              Otherwise, the user  must supply a list using INV_LIST.
! ICNTL(11)    controls whether the user wants to use sequential files.
!              If ICNTL(11) = 0 (default), sequential files are NOT
!              used. Otherwise, factors held sequential files.
! ICNTL(12)    ICNTL(12) controls whether or not the user
!              wants the sequential files used to hold the matrix
!              factors to be deleted at the end of the computation.
!              If ICNTL(12) = 0 (the default), when  the final
!              call is made to MP48A/AD (JOB = 6), the sequential
!              files are deleted. Not deleted otherwise.
! ICNTL(13)    controls whether the user allocates the solution
!              vector data%X. If nonzero, user must allocate
!              data%X. Defualt value 0.
! ICNTL(14) to ICNTL(20) are not currently used but are set to zero


      REAL(WP) CNTL(10)    ! Control parameters.
!  CNTL(1)   has default value 0.5. MP48 switches to full matrix
!            processing if the ratio of number of entries
!            in the reduced matrix to the number that it would have
!            as a full matrix is CNTL(1) or more.
!            A value greater than 1.0 is treated as 1.0
!  CNTL(2)   has default value 0.01 and determines the balance in
!            MP48 between pivoting for sparsity and for
!            stability, values near zero emphasizing sparsity and
!            values near one emphasizing stability.  If CNTL(2) < 0.0,
!            it is regarded as having the value 0.0;
!            if CNTL(2) > 1.0, it is regarded as having the value 1.0.
!  CNTL(3)   has default value zero.  If it is set to a positive value,
!            MP48 will drop from the factors any entry whose
!            modulus is less than CNTL(3).  The factorization will
!            then require less storage but will be inaccurate.
!  CNTL(4)   has default value zero. If MP48 finds a
!            column of the reduced matrix with entries all of
!            modulus less than or equal to CNTL(4), all such entries
!            are dropped from the factorization (and contribute
!            to the count in data%DROP). They also require every pivot
!            to have absolute value greater than CNTL(4).



      INTEGER ERROR             ! Error/Warning flag.

! -1         MPI has not been initialised by the user. Immediate return.
!            An error message is printed on the default output unit.
!  Error and warning diagnostics for JOB = 2
! -2         Either NBLOCK.le. 1 or NBLOCK.gt.NEQ.
! -4         One or more variables in EQVAR is out of range entries
!            (less than 1 or greater than NEQ).
!            This error is also returned if EQVAR is not allocated
!            or if allocated with size less than NE or if
!            NE is less than EQPTR(NEQ+1)-1.
! -5         Error detected in EQPTR.
! -6         ICNTL(7) out of range. Possible values 1,2,3,4,5.
! -7         Either the array NEQSB is not allocated or is
!            of size less than
!            NBLOCK or NEQSB(JBLOCK).lt.1 for one or more
!            submatrix  JBLOCK.
! -9         Either the array INV_LIST is not allocated or is
!            of size less than data%NBLOCK or an
!            entry in INV_LIST out of range (ICNTL(10) nonzero).
! -11        Error in Fortran ALLOCATE statement. The STAT parameter
!            is returned. If the user is not
!            using files (ICNTL(11)=0),
!            it may be possible to avoid this error by rerunning
!            with ICNTL(11).ne.0.
! -13        JOB does not have the same value on all processes
!            or invalid value.
! -19        The call follows a call with JOB = 6.
!  Error and warning diagnostics for JOB = 3
!  -3        Duplicate entries have been found
!            (data%ICNTL(9) nonzero). The number of duplicate
!            entries is data%IDUP.
! -10        If ICNTL(7) = 1 or 2, VALNAM is not allocated
!            or is allocated with the length
!            is than NBLOCK. If ICNTL(7) = 3 or 4, VALUES is not
!            allocated or is allocated but length of VALUES is too
!            small. If ICNTL(7) = 5, RVAL is not allocated
!            or is allocated but length of RVAL is too small.
! -11        Error in Fortran ALLOCATE statement. The STAT parameter
!            is returned.
! -13        JOB does not have the same value on all processes
!            or invalid value.
! -14        Error in Fortran INQUIRE statement. IOSTAT parameter
!            is returned.
! -17        Error in Fortran OPEN statement.
! -19        An error was returned on a previous call or the call
!            follows a call with JOB = 1 (no JOB = 2 call) or
!            follows a call with JOB = 6.
! -20        Failed to find a unit to which a file could be connected.
! +1         Duplicate entries have been found
!            (data%ICNTL(9) = 0). The values of such entries are
!            summed. The number of duplicate entries is data%IDUP.
!  Error and warning diagnostics for JOB = 4
! -10        If ICNTL(7) = 1 or 2, VALNAM is not allocated
!            or is allocated with the length
!            is than NBLOCK. If ICNTL(7) = 3 or 4, VALUES is not
!            allocated or is allocated but length of VALUES is too
!            small. If ICNTL(7) = 5, RVAL is not allocated
!            or is allocated but length of RVAL is too small.
! -11        Error in Fortran ALLOCATE statement. The STAT parameter
!            is returned.
! -12        FACT_JOB has invalid value.
! -13        JOB does not have the same value on all processes
!            or has aninvalid value.
! -14        Error in Fortran INQUIRE statement. IOSTAT parameter
!            is returned.
! -15        On a call with data%FACT_JOB = 3, the matrix entries are
!            unsuitable for the pivot sequence chosen on the
!            previous data%FACT_JOB = 1 or 2 call, or either the order
!            or the number of nonzero entries in the interface matrix
!            differs from the  earlier call.
! -17        Error in Fortran OPEN statement.
! -19        An error was returned on a previous call or the call
!            follows a call with JOB = 1 (no JOB = 2 call) or
!            follows a call with JOB = 6.
! -20        Failed to find a unit to which a file could be connected.
! -21        Interface matrix is structurally rank deficient.
! -25        data%FILES is either not allocated or is allocated but
!            with incorrect size (ICNTL(11)<0).
! +2         The order data%NINTER of the interface matrix is greater
!            than the number data%BORDER of columns in the border.
!            The matrix is singular.
! Possible warnings are:
!  Error and warning diagnostics for JOB = 5
! -11        Error in Fortran ALLOCATE statement. The STAT parameter
!            is returned.
! -13        JOB does not have the same value on all processes
!            or invalid value.
! -14        Error in Fortran INQUIRE statement. IOSTAT parameter
!            is returned.
! -17        Error in Fortran OPEN statement.
! -19        An error was returned on a previous call
!            or the call was not  preceded by a
!            call with JOB = 4, or follows a call with JOB = 6.
! -20        Failed to find a unit to which a file could be connected.
! -22        data%B is either not allocated or
!            is allocated but with incorrect size.
! -23        data%X is either not allocated or is allocated
!            but with incorrect size (ICNTL(13) non zero only).
!  Error and warning diagnostics for JOB = 6
! -13        JOB does not have the same value on all processes
!            or invalid value.
! -14        Error in Fortran INQUIRE statement. IOSTAT parameter
!            is returned.
! -20        Failed to find a unit to which a file could be connected.

      LOGICAL LINTERF      !  Set to TRUE if we get as far as solving
!                             the interface problem
      INTEGER ERROR_SUBBLOCK
! Information returned to user

      INTEGER IDUP         ! Holds the number of duplicate variables
!                            in the row variable lists
      INTEGER STAT         ! STAT parameter.
      INTEGER NPROC        ! Number of processes.

! Number of variables in guard rows
      INTEGER, DIMENSION (:), allocatable :: NGUARD

! Number of submatrices dealt with by each processor
      INTEGER, DIMENSION (:), allocatable :: ICOUNT

! List of submatrices dealt with by each processor
      INTEGER, DIMENSION (:), allocatable :: LIST
! On each processor, IBLOCK holds a list of the submatrices
! assigned to it
      INTEGER, DIMENSION (:), allocatable :: IBLOCK
! ENTRIES holds the length of the variable list for each submatrix
      INTEGER, DIMENSION (:), allocatable :: ENTRIES

! allocatable for LIST
      INTEGER, DIMENSION (:), allocatable :: IPLIST

      INTEGER, DIMENSION (:), allocatable :: INV_LIST
! List of interface columns
      INTEGER, DIMENSION (:), allocatable :: INTER_COL
! Arrays for factors
      REAL(WP), DIMENSION(:), allocatable  :: FACT
      INTEGER,  DIMENSION(:), allocatable  :: IRNF

! Arrays for MP48_MA60 info. on submatrices
      INTEGER,  DIMENSION(:,:), allocatable :: INFO_MA60
      REAL(WP), DIMENSION(:), allocatable   :: OPS
      REAL(WP), DIMENSION(:), allocatable   :: OPSA
! holds flop counts on submatrices
      REAL(WP), DIMENSION(:,:), allocatable  :: STORAGE
! holds real and integer storage for partial factorization of
! submatrices

! Arrays for MA48 info. on interface
      INTEGER INFO_MA48(20)
      REAL(WP)  :: RINFO_MA48(10)

      INTEGER  DROP
      INTEGER  EST_RANK
      REAL(WP) FLOPS
      REAL(WP) NZ
      REAL(WP) STORINT

! Scalars/Arrays for timings
      REAL(DP) :: TIMEI
      REAL(DP) :: TIME_MA48A,TIME_MA48F,TIME_MA48S
      REAL(DP), DIMENSION (:), allocatable :: TIMEA
      REAL(DP), DIMENSION (:), allocatable :: TIMEF
      REAL(DP), DIMENSION (:), allocatable :: TIMES

! NOTE: ALL THESE ARRAYS ARE ONLY SIGNIFICANT ON THE HOST (RANK=0)
! Arrays which must be set up by user.
      INTEGER,DIMENSION(:), allocatable :: EQPTR
      INTEGER,DIMENSION(:), allocatable :: EQVAR
      INTEGER,DIMENSION(:), allocatable :: NEQSB

      REAL(WP),DIMENSION(:), allocatable :: VALUES
      REAL(WP),DIMENSION(:), allocatable :: RVAL

      CHARACTER,DIMENSION(:), allocatable :: VALNAM*128
      CHARACTER,DIMENSION(:,:), allocatable :: FILES*128

      REAL(WP),DIMENSION(:), allocatable :: B
! Solution vector
      REAL(WP),DIMENSION(:), allocatable :: X

      TYPE (MP48_DATA_PRIVATE) private

      END TYPE MP48_DATA

      END MODULE HSL_MP48_DATA_DOUBLE

      MODULE HSL_MP48_DOUBLE
        USE HSL_MP48_DATA_DOUBLE
        PUBLIC MP48AD

        INTERFACE
          SUBROUTINE MP48AD(data)
            USE HSL_MP48_DATA_DOUBLE
            IMPLICIT NONE
            TYPE (MP48_DATA) data
          END SUBROUTINE MP48AD

          SUBROUTINE MP48BD(data)
            USE HSL_MP48_DATA_DOUBLE
            IMPLICIT NONE
            TYPE (MP48_DATA) data
          END SUBROUTINE MP48BD

          SUBROUTINE MP48CD(data)
            USE HSL_MP48_DATA_DOUBLE
            IMPLICIT NONE
            TYPE (MP48_DATA) data
          END SUBROUTINE MP48CD

          SUBROUTINE MP48DD(data)
            USE HSL_MP48_DATA_DOUBLE
            IMPLICIT NONE
            TYPE (MP48_DATA) data
          END SUBROUTINE MP48DD

          SUBROUTINE MP48ED(M,N,KEEP,LKEEP,INFO)
            USE HSL_MP48_DATA_DOUBLE
            IMPLICIT NONE
            INTEGER M,N,LKEEP
            INTEGER KEEP(LKEEP),INFO(5)
          END SUBROUTINE MP48ED

          SUBROUTINE MP48LD(data,ST)
            USE HSL_MP48_DATA_DOUBLE
            IMPLICIT NONE
            TYPE (MP48_DATA) data
            INTEGER ST
          END SUBROUTINE MP48LD

          SUBROUTINE MP48MD(data,NPROC,IFLAG,ST)
            USE HSL_MP48_DATA_DOUBLE
            IMPLICIT NONE
            TYPE (MP48_DATA) data
            INTEGER NPROC,ST
            INTEGER IFLAG(NPROC)
          END SUBROUTINE MP48MD

          SUBROUTINE MP48ND(data,NPROC,IFLAG,FLAG)
            USE HSL_MP48_DATA_DOUBLE
            IMPLICIT NONE
            TYPE (MP48_DATA) data
            INTEGER NPROC,FLAG
            INTEGER IFLAG(NPROC)
          END SUBROUTINE MP48ND

          SUBROUTINE MP48_MA60AD(M,N,NE,LA,A,IRN,JCN,IQ,CNTL, &
                                 ICNTL,IP,NP,JFIRST,LENR,LASTR,NEXTR, &
                                 IW,IFIRST,LENC,LASTC,NEXTC, &
                                 INFO,RINFO)
            USE HSL_MP48_DATA_DOUBLE
            IMPLICIT NONE
            INTEGER M,N,NE,LA,NP
            REAL(WP) A(LA)
            INTEGER ICNTL(10),IRN(LA),JCN(LA),IQ(N), &
                    IP(M),JFIRST(M),LENR(M),LASTR(M),NEXTR(M), &
                    IW(M),IFIRST(N),LENC(N),LASTC(N),NEXTC(N),INFO(20)
            REAL(WP) CNTL(4),RINFO
          END SUBROUTINE MP48_MA60AD

          SUBROUTINE MP48_MA60BD(M,N,NE,JOB,AA,IRNA,IPTRA,CNTL,ICNTL, &
                                 IP,IQ,NP,RANK,NNA,LFACT,FACT,LIRNF, &
                                 IRNF,IPTRL,IPTRU,W,IW,INFO,RINFO)
            USE HSL_MP48_DATA_DOUBLE
            IMPLICIT NONE
            INTEGER M,N,NE,JOB,RANK,NNA,NP,LFACT,LIRNF
            REAL(WP) CNTL(4),W(M),RINFO,AA(NE),FACT(LFACT)
            INTEGER IRNA(NE),IPTRA(N),IRNF(LIRNF),IPTRL(N),IPTRU(N), &
                    ICNTL(10),IP(M),IQ(*),IW(M+2*N),INFO(20)
          END SUBROUTINE MP48_MA60BD

          SUBROUTINE MP48_MA60CD(M,N,ICNTL,IQ,NP,LORU,LFACT,FACT,LIRNF, &
                            IRNF,IPTRL,IPTRU,B,X,W,INFO)
            USE HSL_MP48_DATA_DOUBLE
            IMPLICIT NONE
            INTEGER LFACT,LIRNF,M,N,NP,LORU
            REAL(WP) B(MAX(M,N)),FACT(LFACT),X(MAX(M,N)),W(MAX(M,N))
            INTEGER INFO(2),ICNTL(7),IQ(N),IRNF(LIRNF),IPTRL(N),IPTRU(N)
          END SUBROUTINE MP48_MA60CD

          SUBROUTINE MP48_MA60DD(LA,A,IND,IPTR,N,DISP,REALS)
            USE HSL_MP48_DATA_DOUBLE
            IMPLICIT NONE
            INTEGER LA,N,DISP
            REAL(WP) A(LA)
            LOGICAL REALS
            INTEGER IPTR(N),IND(LA)
          END SUBROUTINE MP48_MA60DD

          SUBROUTINE MP48_MA60ED(M,N,A,LDA,PIVTOL,IPIV,RANK,NNO)
            USE HSL_MP48_DATA_DOUBLE
            IMPLICIT NONE
            INTEGER LDA,M,N,RANK,NNO
            REAL(WP) PIVTOL
            INTEGER IPIV(N)
            REAL(WP) A(LDA,N)
          END SUBROUTINE MP48_MA60ED

          SUBROUTINE MP48_MA60FD(M,N,A,LDA,PIVTOL,IPIV,RANK,NNO)
            USE HSL_MP48_DATA_DOUBLE
            IMPLICIT NONE
            INTEGER LDA,M,N,RANK,NNO
            REAL(WP) PIVTOL
            INTEGER IPIV(N)
            REAL(WP) A(LDA,N)
          END SUBROUTINE MP48_MA60FD

          SUBROUTINE MP48_MA60GD(M,N,A,LDA,NB,PIVTOL,IPIV,RANK,NNO)
            USE HSL_MP48_DATA_DOUBLE
            IMPLICIT NONE
            INTEGER LDA,M,N,NB,RANK,NNO
            REAL(WP) PIVTOL
            INTEGER IPIV(N)
            REAL(WP) A(LDA,N)
          END SUBROUTINE MP48_MA60GD

          SUBROUTINE MP48_MA60ID(CNTL,ICNTL)
            USE HSL_MP48_DATA_DOUBLE
            IMPLICIT NONE
            REAL(WP) CNTL(4)
            INTEGER ICNTL(10)
          END SUBROUTINE MP48_MA60ID

          SUBROUTINE MP48_MA60HD(LORU,M,N,A,LDA,IPIV,B,ICNTL5,NNO)
            USE HSL_MP48_DATA_DOUBLE
            IMPLICIT NONE
            INTEGER LORU,LDA,M,N,ICNTL5,NNO
            INTEGER IPIV(N)
            REAL(WP) A(LDA,N),B(*)
          END SUBROUTINE MP48_MA60HD


        END INTERFACE
      END MODULE HSL_MP48_DOUBLE
!********************************************************
      SUBROUTINE MP48AD(data)

      USE HSL_MP48_DATA_DOUBLE
      IMPLICIT NONE

      TYPE (MP48_DATA) data

!     .. Parameters ..
      REAL(WP), PARAMETER :: ZERO = 0.0_WP

!     .. Local Scalars ..
      INTEGER DGN,ERCODE,ERR,FLAG,FLAGMIN,FLAGMAX, &
              I,ICALL,IPROC,IOS,ISTRT,ISTOP, &
              J,J1,J2,JBLOCK,JROW, &
              K,KL, &
              L,L1,LCNT,LDGN,LGUARD, &
              LLPTR,LORDER,LP, &
              N,NBLOCK,NEQ,NIN,NPROC,NVAR,NCOL,NL,NNULL, &
              RANK,ST,STRM,WRN
      LOGICAL EX,OPEN,LERR,LFLAG,LWRN

!     .. Local Arrays ..
      REAL(WP), DIMENSION (:), ALLOCATABLE :: LOAD
      REAL(WP), DIMENSION (:), ALLOCATABLE :: W

      REAL(DP) T1,T2

      INTEGER IBUF(10),ICONTROL(30),STAT(MPI_STATUS_SIZE,20)
      INTEGER, DIMENSION (:), ALLOCATABLE :: IFLAG
      INTEGER, DIMENSION (:), ALLOCATABLE :: IW
      INTEGER, DIMENSION (:), ALLOCATABLE :: IWORK
      INTEGER, DIMENSION (:,:), ALLOCATABLE :: IGUARD

!     .. External Subroutines ..
!     EXTERNAL MA52AD

!     .. MPI routines
      EXTERNAL MPI_BARRIER,MPI_BCAST,MPI_COMM_SIZE, &
               MPI_INITIALIZED,MPI_GATHER,MPI_ALLREDUCE

!     .. Intrinsic Functions ..
      INTRINSIC MAX

! Check mpi has been initialized
      CALL MPI_INITIALIZED(LFLAG,ERCODE)
! LFLAG is true if MPI_INIT has been called, FALSE otherwise
      IF (.NOT. LFLAG) THEN
         data%ERROR = -1
         WRITE (*,FMT=9020) data%ERROR
         WRITE (*,*) 'MPI is not running'
         RETURN
      END IF

      IF (data%JOB == 1) THEN
         CALL MPI_COMM_DUP(data%COMM,data%private%COMM,ERCODE)
         CALL MPI_COMM_RANK(data%private%COMM,data%RANK,ERCODE)
         CALL MPI_COMM_SIZE(data%private%COMM,data%NPROC,ERCODE)
         RANK = data%RANK
         NPROC = data%NPROC
! ======================
! Set control parameters
! ======================

! Error output unit
         data%ICNTL(1) = 6
! Warning output unit
         data%ICNTL(2) = 6
! Diagnostics output unit
         data%ICNTL(3) = 6
! Level of diagnostic printing
         data%ICNTL(4) = 1

! data%ICNTL(5) controls full-matrix processing.
         data%ICNTL(5) = 32
! data%ICNTL(6) controls increasing factor storage
         data%ICNTL(6) = 2
! data%ICNTL(7) controls input of data.
         data%ICNTL(7) = 1
! data%ICNTL(8) controls number of cols. searched for pivots.
         data%ICNTL(8) = 3
! data%ICNTL(9) controls action if duplicates found.
         data%ICNTL(9) = 0
! data%ICNTL(10) controls whether lists of submatrices for each
! processor supplied or not.
         data%ICNTL(10) = 0
! data%ICNTL(11) controls whether files used for factor data.
         data%ICNTL(11) = 0
! data%ICNTL(12) controls whether files are to be deleted.
         data%ICNTL(12) = 0
! data%ICNTL(13) controls whether user allocates solution vector.
         data%ICNTL(13) = 0

         data%ICNTL(14:20) = 0

         data%CNTL(1) = 0.5_WP
         data%CNTL(2) = 0.01_WP
         data%CNTL(3:10) = ZERO

         data%ERROR   = 0
         data%DROP    = 0
         data%FLOPS   = ZERO
         data%IOSTAT  = 0
         data%NZ      = ZERO
         data%STORINT = ZERO

         data%private%OLDJOB = 1
         data%private%ERROR = 0

         CALL MPI_BARRIER(data%private%COMM,ERCODE)
         RETURN
      END IF

! ===============
! Initializations
! ===============
! This is where we enter if JOB.ne.1

! Ensure communicator is defined to allow normal error return.
      IF(data%private%COMM == MPI_COMM_NULL) data%private%COMM = data%COMM
! Make sure all processes are sychronized.
      CALL MPI_BARRIER(data%private%COMM,ERCODE)
      RANK = data%RANK
      NPROC = data%NPROC
! Start by broadcasting the control parameters
! We are doing this on each call (does mean control parameters
! can have been changed by user beteween calls ... could cause problems)
      CALL MPI_BCAST(data%ICNTL,30,MPI_INTEGER,0,data%private%COMM, &
                     ERCODE)
      CALL MPI_BCAST(data%CNTL,10,MPI_DOUBLE_PRECISION, &
                     0,data%private%COMM,ERCODE)

! Set diagnostics controls
      DGN = data%ICNTL(3)
      LDGN = data%ICNTL(4)
      IF (DGN < 0) LDGN = 0
! Set error controls
      LERR = .FALSE.
      IF (data%ICNTL(1).GE.0 .AND. data%ICNTL(4).GE.1) LERR = .TRUE.
      ERR = data%ICNTL(1)
! Set warning controls
      LWRN = .FALSE.
      IF (data%ICNTL(2).GE.0 .AND. data%ICNTL(4).GE.1) LWRN = .TRUE.
      WRN = data%ICNTL(2)

! ICONTROL(7) controls input of matrices.
      ICONTROL(7) =  data%ICNTL(7)
! ICONTROL(10) controls whether lists of submatrices for each processor
! supplied or not.
      ICONTROL(10) =  data%ICNTL(10)
! ICONTROL(11) controls whether files used to reduce storage
! requirements.
      ICONTROL(11) =  data%ICNTL(11)
! ICONTROL(12) controls whether files are to be deleted.
      ICONTROL(12) =  data%ICNTL(12)

! JOB = 1 : initialisation call
! JOB = 2 : pre analyse
! JOB = 3 : analyse
! JOB = 4 : factorize submatrices in parallel
! JOB = 5 : solve for right-hand sides
! JOB = 6 : deallocate all arrays (final call)
! Also combination calls allowed
! 23 = 2 + 3
! 24 = 2 + 3 + 4
! 25 = 2 + 3 + 4 + 5

      CALL MPI_BARRIER(data%private%COMM,ERCODE)
! Check that value of JOB is the same on all processes
      CALL MPI_ALLREDUCE(data%JOB,FLAGMAX,1,MPI_INTEGER,MPI_MAX, &
                         data%private%COMM,ERCODE)
      CALL MPI_ALLREDUCE(data%JOB,FLAGMIN,1,MPI_INTEGER,MPI_MIN, &
                         data%private%COMM,ERCODE)

      IF ( FLAGMAX  /=  FLAGMIN ) THEN
         data%ERROR = -13
         IF (RANK == 0 .AND. LERR) THEN
            WRITE (ERR,FMT=9020) data%ERROR
            WRITE (ERR,'(/A)') ' Different values of data%JOB'
         END IF
         GO TO 890
      END IF
      IF ( FLAGMAX  >  6 .OR. FLAGMIN < 2 ) THEN
        IF (FLAGMAX /= 23 .AND. FLAGMAX /= 24 .AND. &
             FLAGMAX /= 25) THEN
          data%ERROR = -13
          IF (RANK == 0 .AND. LERR) THEN
            WRITE (ERR,FMT=9020) data%ERROR
            WRITE (ERR,'(/A)') ' Value of data%JOB not valid'
          END IF
          GO TO 890
        END IF
      END IF

      IF (RANK == 0 .AND. LDGN > 1) THEN
         IF (data%JOB /= 4) WRITE (DGN,FMT='(/A,I5)') &
           ' Entering MP48 with data%JOB  = ',data%JOB
         IF (data%JOB == 4) WRITE (DGN,FMT='(/A,I2,A,I2)') &
           ' Entering MP48 with data%JOB  = ',data%JOB, &
           ' and data%FACT_JOB = ',data%FACT_JOB
      END IF

! If JOB = 24 or 25, we will ignore whatever user has in FACT_JOB and
! set it equal to 1 (first call)
      IF (data%JOB == 24 .OR. data%JOB == 25) THEN
         data%FACT_JOB = 1
         data%private%FJOB = 1
      END IF

! Check calling sequence
      IF (data%JOB == 2 .OR. data%JOB == 23 .OR. &
          data%JOB == 24 .OR. data%JOB == 25) THEN
        IF (data%private%OLDJOB == 6) THEN
          data%ERROR = -19
          IF (RANK == 0 .AND. LERR) THEN
            WRITE (ERR,FMT=9000) data%JOB,data%ERROR
            WRITE (ERR,FMT='(/A)') &
          ' Call does not follow a call with data%JOB = 1'
          END IF
          GO TO 890
        END IF
      END IF
! Set data%OLDJOB to 2 to show call with JOB = 2 has been made.
      IF (data%JOB == 2) data%private%OLDJOB = 2

      IF (data%JOB == 3) THEN
         IF (data%private%ERROR < 0 .OR. data%private%OLDJOB == 1 &
           .OR. data%private%OLDJOB == 6) THEN
           data%ERROR = -19
           IF (RANK == 0 .AND. LERR) THEN
               WRITE (ERR,FMT=9000) data%JOB,data%ERROR
               IF (data%private%OLDJOB == 1 &
                  .OR. data%private%OLDJOB == 6) THEN
                  WRITE (ERR,FMT=9130) data%private%OLDJOB
               ELSE
                  WRITE (ERR,FMT=9120)
               END IF
            END IF
            GO TO 890
         END IF
      END IF
      IF (data%JOB == 3 .OR. data%JOB == 23) THEN
! Set OLDJOB to 3 to show call with JOB = 3 has been made.
         data%private%OLDJOB = 3
! Initialise data%private%FJOB
         data%private%FJOB = 0
      END IF
      IF (data%JOB == 3) GO TO 200

      IF (data%JOB == 4) THEN
         IF (data%private%ERROR < 0 .OR. data%private%OLDJOB == 1 &
           .OR. data%private%OLDJOB == 2 .OR. &
           data%private%OLDJOB == 6) THEN
           data%ERROR = -19
           IF (RANK == 0 .AND. LERR) THEN
               WRITE (ERR,FMT=9000) data%JOB,data%ERROR
               IF (data%private%OLDJOB == 1 &
                  .OR. data%private%OLDJOB == 6) THEN
                  WRITE (ERR,FMT=9130) data%private%OLDJOB
               ELSE
                  WRITE (ERR,FMT=9120)
               END IF
            END IF
            GO TO 890
         END IF
         CALL MPI_BCAST(data%FACT_JOB,1,MPI_INTEGER,0,data%private%COMM, &
                        ERCODE)
         IF (data%FACT_JOB < 1 .OR. data%FACT_JOB > 3) THEN
           data%ERROR = -12
         ELSE IF (data%private%FJOB == 0 .AND. &
                  data%FACT_JOB == 3) THEN
           data%ERROR = -12
         ELSE IF (data%FACT_JOB == 3) THEN
           IF (data%CNTL(3) /= ZERO .OR. data%CNTL(4) /= ZERO) &
           data%ERROR = -12
         END IF
         IF (data%ERROR == -12 .AND. RANK == 0 .AND. LERR) THEN
            WRITE (ERR,FMT=9000) data%JOB,data%ERROR
            WRITE (ERR,FMT=9080) &
           'data%FACT_JOB','data%FACT_JOB',data%FACT_JOB
            IF (data%private%FJOB == 0 .AND. data%FACT_JOB == 3) &
            WRITE (ERR,FMT=9160)
         END IF
         IF (data%ERROR == -12) GO TO 890
! Increment data%private%FJOB by 1.
! data%private%FJOB is the number of calls with JOB = 4
! made following a JOB = 3 call
         data%private%FJOB = data%private%FJOB + 1
         CALL MPI_BCAST(data%private%FJOB,1,MPI_INTEGER,0, &
                        data%private%COMM,ERCODE)
      END IF

      IF (data%JOB == 4 .OR. data%JOB == 24) THEN
! Set data%OLDJOB to 4 to show call with JOB = 4 has been made.
         data%private%OLDJOB = 4
      END IF
      IF (data%JOB == 4) GO TO 300

      IF (data%JOB == 5) THEN
         IF (data%private%OLDJOB == 1 &
           .OR. data%private%OLDJOB == 2 &
           .OR. data%private%OLDJOB == 3 &
           .OR. data%private%OLDJOB == 6) THEN
            data%ERROR = -19
            IF (RANK == 0 .AND. LERR) THEN
               WRITE (ERR,FMT=9000) data%JOB,data%ERROR
               WRITE (ERR,FMT=9130) data%private%OLDJOB
            END IF
            GO TO 890
         END IF
      END IF
      IF (data%JOB == 5 .OR. data%JOB == 25) THEN
! Set data%OLDJOB to 5 to show call with JOB = 5 has been made.
         data%private%OLDJOB = 5
      END IF
      IF (data%JOB == 5) GO TO 400

      IF (data%JOB == 6) THEN
! Set data%OLDJOB to 6 to show call with JOB = 6 has been made.
         data%private%OLDJOB = 6
         GO TO 900
      END IF

! At this point, JOB = 2 or 23 or 24 or 25
      data%ERROR = 0
      FLAG = 0

! Initialise STAT parameter
      data%STAT = 0
      data%LINTERF = .FALSE.

! ================================
! Print some information (by root)
! ================================
      IF (RANK == 0 .AND. LDGN > 1) THEN

         WRITE (DGN,FMT='(A,I8)') &
         ' Number of processes is    = ',NPROC
         WRITE (DGN,FMT='(A,I8)') &
         ' Number of submatrices is  = ',data%NBLOCK
         WRITE (DGN,FMT='(A,I8)') &
         ' Order of matrix NEQ       = ',data%NEQ
         WRITE (DGN,9140) (data%ICNTL(I),I=1,13)
         WRITE (DGN,9150) (data%CNTL(I),I=1,4)
      END IF

      IF (RANK == 0) THEN
         NBLOCK = data%NBLOCK
         NEQ = data%NEQ

! Initialise parameters so that they are not undefined
         data%IDUP = 0
         data%ERROR_SUBBLOCK = 0

! Perform error checks on the data
! There must be at least one row in each submatrix
! so the total number of rows must be at least NBLOCK
         IF (NBLOCK.LE.1 .OR. NBLOCK > NEQ ) THEN
            data%ERROR = -2
            FLAG = -2
            IF (LERR) THEN
               WRITE (ERR,FMT=9000) data%JOB,data%ERROR
               IF (NBLOCK.LE.1) THEN
                  WRITE (ERR,FMT=9080) ' data%NBLOCK', &
                 'data%NBLOCK',NBLOCK
               ELSE IF (NBLOCK > NEQ) THEN
                  WRITE (ERR,FMT=9080)' data%NEQ','data%NEQ',NEQ
               END IF
            END IF
            GO TO 20
         END IF

! Check EQPTR allocated and check its size
         IF (.NOT. allocated(data%EQPTR)) THEN
            FLAG = -5
         ELSE
            IF (SIZE(data%EQPTR) < data%NEQ+1) FLAG = -5
         END IF
         IF (FLAG == -5) THEN
            data%ERROR = -5
            IF (LERR) THEN
               WRITE (ERR,FMT=9000) data%JOB,data%ERROR
               IF (.NOT. allocated(data%EQPTR)) THEN
                  WRITE (ERR,FMT=9010) ' data%EQPTR'
               ELSE IF (SIZE(data%EQPTR) < data%NEQ+1) THEN
                  WRITE (ERR,FMT=9100) ' data%EQPTR',SIZE(data%EQPTR), &
                  data%NEQ+1
               END IF
            END IF
            GO TO 20
         END IF

         IF (ICONTROL(7) < 1 .OR. ICONTROL(7) > 5) THEN
            data%ERROR = -6
            FLAG = -6
            IF (LERR) THEN
               WRITE (ERR,FMT=9000) data%JOB,data%ERROR
               WRITE (ERR,FMT=9080) &
              'data%ICNTL(7)','data%ICNTL(7)',ICONTROL(7)
            END IF
            GO TO 20
         END IF

! Check NEQSB has been allocated and check its size
         IF (.NOT. allocated(data%NEQSB)) THEN
            FLAG = -7
         ELSE
            IF (SIZE(data%NEQSB) < NBLOCK) FLAG = -7
         END IF
         IF (FLAG == -7) THEN
            data%ERROR = -7
            IF (LERR) THEN
               WRITE (ERR,FMT=9000) data%JOB,data%ERROR
               IF (.NOT. allocated(data%NEQSB)) THEN
                  WRITE (ERR,FMT=9010) ' data%NEQSB'
               ELSE IF (SIZE(data%NEQSB) < NBLOCK) THEN
                  WRITE (ERR,FMT=9100) ' data%NEQSB',SIZE(data%NEQSB), &
                     NBLOCK
               END IF
            END IF
            GO TO 20
         END IF
         DO 8 JBLOCK = 1,NBLOCK
            IF (data%NEQSB(JBLOCK).LE.0) THEN
               data%ERROR = -7
               data%ERROR_SUBBLOCK = JBLOCK
               IF (LERR) THEN
                  IF (FLAG == 0) THEN
                     WRITE (ERR,FMT=9000) data%JOB,data%ERROR
                     WRITE (ERR,FMT='(A)') ' Check data%NEQSB'
                  END IF
                  WRITE (ERR,FMT=9090) &
                  ' data%NEQSB',data%NEQSB(JBLOCK),JBLOCK
               END IF
               FLAG = -7
            END IF
    8    CONTINUE
         IF (FLAG == -7) GO TO 20

! Check size of data%EQVAR
! First check it has been allocated
         IF (.NOT. allocated(data%EQVAR)) THEN
            FLAG = -4
         ELSE
            IF (SIZE(data%EQVAR) < data%EQPTR(NEQ+1)-1) FLAG = -4
         END IF
         IF (FLAG == -4) THEN
            data%ERROR = -4
            IF (LERR) THEN
               WRITE (ERR,FMT=9000) data%JOB,data%ERROR
               IF (.NOT. allocated(data%EQVAR)) THEN
                  WRITE (ERR,FMT=9010) ' data%EQVAR'
               ELSE IF (SIZE(data%EQVAR) < data%EQPTR(NEQ+1)-1) THEN
                  WRITE (ERR,FMT=9100) 'data%EQVAR', &
                  SIZE(data%EQVAR), data%EQPTR(NEQ+1)-1
               END IF
            END IF
            GO TO 20
         END IF
         data%NE = data%EQPTR(NEQ+1)-1

         IBUF(1) = NBLOCK
         IBUF(2) = NEQ
         IBUF(3) = data%NE

      END IF

! Broadcast error status from root
   20 CALL MPI_BARRIER(data%private%COMM,ERCODE)
      CALL MPI_BCAST(data%ERROR,1,MPI_INTEGER,0,data%private%COMM,ERCODE)
! Terminate if error occurred
      IF (data%ERROR < 0) GO TO 890

! Broadcast data from root
      CALL MPI_BCAST(IBUF,3,MPI_INTEGER,0,data%private%COMM,ERCODE)
      NBLOCK  = IBUF(1)
      NEQ     = IBUF(2)
      N       = NEQ
      data%NE = IBUF(3)
      data%NBLOCK = NBLOCK
      data%NEQ = NEQ
      data%N   = N

! Allocate arrays
      DEALLOCATE (IFLAG,STAT=ST)
      ALLOCATE (IFLAG(NPROC),STAT=ST)
! Initialise the error flag (so that it is always well-defined)
      IF (ST == 0) IFLAG(1:NPROC) = 0

      DEALLOCATE (data%ICOUNT, &
                  data%NGUARD, &
                  data%private%NCOLA,data%ENTRIES, &
                  STAT=ST)
      ALLOCATE (data%ICOUNT(0:NPROC-1), &
                data%NGUARD(1:NBLOCK), &
                data%private%NCOLA(1:NBLOCK),data%ENTRIES(1:NBLOCK), &
                STAT=ST)
      IF (ST /= 0) GO TO 100
      data%ICOUNT(0:NPROC-1) = 0

      IF (RANK == 0) THEN
         DEALLOCATE (data%private%SUBPTR,STAT=ST)
         ALLOCATE (data%private%SUBPTR(NBLOCK+1),STAT=ST)
      ELSE
         DEALLOCATE (data%NEQSB,STAT=ST)
         ALLOCATE (data%NEQSB(1:NBLOCK),STAT=ST)
      END IF
      IF (ST /= 0) GO TO 100

      FLAG = 0
      IF (ICONTROL(10) == 0 .OR. ICONTROL(10) == 1) THEN
         DEALLOCATE (data%INV_LIST,STAT=ST)
         ALLOCATE (data%INV_LIST(NBLOCK),STAT=ST)
         IF (ST /= 0) GO TO 100

      ELSE
! data%INV_LIST has been set by user. Check array allocated
! and with correct size.
         IF (RANK == 0) THEN
            IF (.NOT. allocated(data%INV_LIST)) THEN
               FLAG = -9
            ELSE
               IF (SIZE(data%INV_LIST) < NBLOCK) FLAG = -9
            END IF
            IF (FLAG == 0) THEN
               DO JBLOCK = 1,NBLOCK
                  IPROC = data%INV_LIST(JBLOCK)
                  IF (IPROC < 0 .OR. IPROC > NPROC-1) FLAG = -9
               END DO
            END IF
         ELSE
            DEALLOCATE (data%INV_LIST,STAT=ST)
            ALLOCATE (data%INV_LIST(NBLOCK),STAT=ST)
            IF (ST /= 0) GO TO 100
         END IF
      END IF

      DEALLOCATE (data%IPLIST,data%LIST, &
                  data%IBLOCK,STAT=ST)
      ALLOCATE (data%IPLIST(0:NPROC), &
                data%LIST(1:NBLOCK), &
                data%IBLOCK(1:NBLOCK), &
                STAT=ST)

  100 CALL MPI_BARRIER(data%private%COMM,ERCODE)
      CALL MP48MD(data,NPROC,IFLAG,ST)
! Broadcast error status from root
      CALL MPI_BCAST(data%ERROR,1,MPI_INTEGER,0,data%private%COMM,ERCODE)
! Terminate if error occurred
      IF (data%ERROR < 0) GO TO 890

      IF (ICONTROL(10) /= 0 .AND. ICONTROL(10) /= 1) THEN
         IF (RANK == 0 .AND. FLAG == -9) THEN
            data%ERROR = -9
            WRITE (ERR,FMT=9000) data%JOB,data%ERROR
            IF (.NOT. allocated(data%INV_LIST)) THEN
               WRITE (ERR,FMT=9010) 'data%INV_LIST'
            ELSE IF (SIZE(data%INV_LIST) < NBLOCK) THEN
               WRITE (ERR,FMT='(A)') &
             ' Size of data%INV_LIST less than data%NBLOCK'
            ELSE
               WRITE (ERR,FMT='(A)')' Error in data%INV_LIST'
            END IF
         END IF
         CALL MPI_BCAST(data%ERROR,1,MPI_INTEGER,0,data%private%COMM,ERCODE)
! Terminate if error occurred
         IF (data%ERROR < 0) GO TO 890

         CALL MPI_BCAST(data%INV_LIST,NBLOCK, &
                        MPI_INTEGER,0,data%private%COMM,ERCODE)
      END IF

      IF (RANK == 0) THEN

! data%ENTRIES holds the length of variable list for each submatrix
! LORDER is max. number of rows in a block
! (and is first dimension of GLOBAL)
         LORDER =1
         LLPTR = 1
         data%ENTRIES(1:NBLOCK)  = 0
         DO 120 JBLOCK = 1,NBLOCK
            LORDER = MAX(LORDER,data%NEQSB(JBLOCK)+1)
! Loop over rows in submatrix
            DO L = 1,data%NEQSB(JBLOCK)
               NVAR =  data%EQPTR(LLPTR+1) - data%EQPTR(LLPTR)
               IF (NVAR < 1) THEN
                  !print *, "nvar ",NVAR,"llpt",LLPTR
                  FLAG = -5
                  data%ERROR = -5
                  IF (LERR) THEN
                     !print *, "here "
                     WRITE (ERR,FMT=9000) data%JOB,data%ERROR
                     WRITE (ERR,FMT='(A)') ' Error in data%EQPTR.'
                     WRITE (ERR,FMT='(A)') ' One or more rows are null.'
                  END IF
                  GO TO 180
               END IF
               LLPTR = LLPTR + 1
! Subdomain allocatables
               data%ENTRIES(JBLOCK) = &
               data%ENTRIES(JBLOCK) + NVAR
            END DO
  120    CONTINUE

! Set SUBPTR(K) to point to the beginning of the list of rows
! in submatrix K. SUBPTR(NBLOCK+1) points to one after data%NEQSB.
         data%private%SUBPTR(1) = 1
         DO JBLOCK = 2,NBLOCK + 1
            data%private%SUBPTR(JBLOCK) = data%NEQSB(JBLOCK-1) + &
                                          data%private%SUBPTR(JBLOCK-1)
         END DO

         T1 = MPI_WTIME()

! Generate interface variable lists using MA52A/AD (done on root).
         LGUARD = data%ENTRIES(1)
         DO JBLOCK = 2,NBLOCK
            LGUARD = MAX(LGUARD,data%ENTRIES(JBLOCK))
         END DO
         LGUARD = LGUARD/2

         DEALLOCATE (IW,IWORK,STAT=ST)
         ALLOCATE (IW(NEQ), &
                   IWORK(2*NEQ),STAT=ST)
         IF (ST /= 0) GO TO 180
         DEALLOCATE (data%private%CGLOB, &
                     data%private%GLOBAL, &
                     data%private%INULL, &
                     STAT=ST)
         ALLOCATE (data%private%CGLOB(data%MAXSBCOLS,NBLOCK), & !NEQ
                   data%private%GLOBAL(LORDER,NBLOCK), &
                   data%private%INULL(NEQ), &
                   STAT=ST)
         IF (ST /= 0) GO TO 180

  130    DEALLOCATE (IGUARD,STAT=ST)
         ALLOCATE (IGUARD(1:LGUARD,1:NBLOCK),STAT=ST)
         IF (ST /= 0) GO TO 180

         ICALL = 0
         data%private%INULL(1:NEQ) = 0
         NNULL = NEQ
         L1 = 1
         DO 160 JBLOCK = 1,NBLOCK
            JROW = 1
            DO 150 L = 1,data%NEQSB(JBLOCK)
               J1 = data%EQPTR(L1)
               J2 = data%EQPTR(L1+1)-1
               NVAR = J2 - J1 + 1
               ICALL = ICALL + 1
! Store global row number of row JROW
               data%private%GLOBAL(JROW,JBLOCK) = ICALL
! If not wanted, switch off MA52 printing.
! Note that MA52 checks for indices that are out-of-range
               LP = -1
               CALL MA52AD(ICALL,NVAR,data%EQVAR(J1:J2), &
                           NEQ,NBLOCK,JBLOCK,data%NGUARD,LGUARD, &
                           IGUARD,NEQ,IWORK,LP,FLAG)
! Check for out of range entries
               IF (FLAG == -4) THEN
                  data%ERROR = -4
                  IF (LERR) THEN
                    WRITE (ERR,FMT=9000) data%JOB,data%ERROR
                    WRITE (ERR,FMT='(A)')' Check row variable lists'
                    WRITE (ERR,FMT='(A)') &
                  ' One or more variables is out of range '
                  END IF
                  GO TO 180
! If failed because insufficient space, allocate more and recall
               ELSE IF (FLAG == -6) THEN
                  LGUARD = IWORK(1)
                  GO TO 130
               END IF
               JROW = JROW + 1
               L1 = L1 + 1
! Check if we have any null cols (do this after we checked for
! out of range entries in EQVAR)
               DO J = J1,J2
                 K = data%EQVAR(J)
                 IF (data%private%INULL(K) == 0) THEN
                   data%private%INULL(K) = 1
                   NNULL = NNULL - 1
                 END IF
              END DO

  150       CONTINUE
  160    CONTINUE
         data%NULL = NNULL
!        write (6,*) 'null',NNULL
         IF (data%NULL == 0) DEALLOCATE (data%private%INULL,STAT=ST)

      END IF

  180 CALL MPI_BARRIER(data%private%COMM,ERCODE)
      IF (RANK == 0 .AND. ST /= 0) CALL MP48LD(data,ST)
! Broadcast error status from root
      CALL MPI_BCAST(data%ERROR,1,MPI_INTEGER,0,data%private%COMM,ERCODE)
! Terminate if error occurred
      IF (data%ERROR < 0) GO TO 890

      IF (RANK == 0) THEN
! Add up total number of cols in border variables
         IWORK(1:NEQ) = 0
         NIN = 0
         DO JBLOCK = 1,NBLOCK
           DO L = 1,data%NGUARD(JBLOCK)
             I = IGUARD(L,JBLOCK)
             IF (IWORK(I) == 0) THEN
! Interface variable ... set flag
               NIN = NIN + 1
               IWORK(I) = 1
             END IF
           END DO
         END DO
! Store number of border columns
         data%BORDER = NIN

         IWORK(1:2*NEQ) = 0
         DO 190 JBLOCK = 1,NBLOCK
! Flag columns that appear within the block and add up
! number of columns in block
            NCOL = 0
            L1 = data%private%SUBPTR(JBLOCK)
            DO L = 1,data%NEQSB(JBLOCK)
               J1 = data%EQPTR(L1)
               J2 = data%EQPTR(L1+1)-1
               DO J = J1,J2
                  K = data%EQVAR(J)
                  IF (IWORK(K) == 0) THEN
                     IWORK(K) = 1
                     NCOL = NCOL + 1
                  END IF
               END DO
               L1 = L1 + 1
            END DO
! Flag columns that are in the border for the block
            DO L = 1,data%NGUARD(JBLOCK)
               K = IGUARD(L,JBLOCK)
! Flag cols. in border
               IWORK(K) = 2
            END DO

! We need to locally renumber the columns in the block
! Note: the same column can have a different local number
!       in each block.

            NL = 0
            KL = NCOL - data%NGUARD(JBLOCK)
            L1 = data%private%SUBPTR(JBLOCK)
            DO L = 1,data%NEQSB(JBLOCK)
               J1 = data%EQPTR(L1)
               J2 = data%EQPTR(L1+1)-1
               DO J = J1,J2
                  K = data%EQVAR(J)
                  !if(L1.EQ.7071819)print *, "m ", J, "var ",data%EQVAR(J),"iwk ",IWORK(K),"iwk +k ",IWORK(NEQ+K)
                  IF (IWORK(K) == 1) THEN
! First appearance of col. index K
                     NL = NL + 1
! Overwrite EQVAR with LOCAL column indices
                     data%EQVAR(J) = NL
                     IWORK(NEQ+K) = -NL
! Restore IWORK ready for next block
                     IWORK(K) = 0
! CGLOB holds global col. indices
                     data%private%CGLOB(NL,JBLOCK) = K
!                    write (6,*) 'jblock,l,j,k,nl',jblock,l,j,k,nl
                  ELSE IF (IWORK(K) == 2) THEN
! First appearance of col. index K which lies in border
                     KL = KL + 1
                     data%EQVAR(J) = KL
                     IWORK(NEQ+K) = -KL
                     IWORK(K) = 0
                     data%private%CGLOB(KL,JBLOCK) = K
!                    write (6,*) 'jblock,l,j,k,kl',jblock,l,j,k,kl
                  ELSE
                     data%EQVAR(J) = -IWORK(NEQ+K)
                  END IF
                  !if(L1.EQ.7071819)print *, "m ", J, "var ",data%EQVAR(J)
               END DO
               L1 = L1 + 1
            END DO

! Store the number of columns NL not in the border
            data%private%NCOLA(JBLOCK) = NL
! Store total number of cols. in submatrix
            data%private%NCOLA(JBLOCK) = &
               data%private%NCOLA(JBLOCK)+data%NGUARD(JBLOCK)

  190    CONTINUE

         IF (LDGN > 1) THEN
            WRITE (DGN,FMT='(A,/8(I10))') &
              ' The number of rows in the submatrices is:', &
              (data%NEQSB(JBLOCK),JBLOCK=1,NBLOCK)
            WRITE (DGN,FMT='(A,/8(I10))') &
              ' The number of columns in the submatrices is:', &
              (data%private%NCOLA(JBLOCK),JBLOCK=1,NBLOCK)
            WRITE (DGN,FMT='(A,/8(I10))') &
              ' The number of border columns for the submatrices is:', &
              (data%NGUARD(JBLOCK),JBLOCK=1,NBLOCK)
            WRITE (DGN,FMT='(A,/8(I10))') &
              ' The number of entries in the submatrices is:', &
                (data%ENTRIES(JBLOCK),JBLOCK=1,NBLOCK)
         END IF
         T2 = MPI_WTIME()
         data%TIMEI = T2 - T1
         IF (LDGN.GE.3) &
            WRITE (DGN,FMT='(A,F9.3,A)') ' Host took ',data%TIMEI, &
          ' seconds to generate the list of border columns'

         DEALLOCATE (IGUARD,STAT=ST)

! data%EQVAR now holds local column indices, with the cols in the border
! of each block ordered last.
      END IF

! Broadcast NGUARD from root
      CALL MPI_BCAST(data%NGUARD,NBLOCK,MPI_INTEGER,0,data%private%COMM,ERCODE)
! Broadcast LORDER from root
      CALL MPI_BCAST(LORDER,1,MPI_INTEGER,0,data%private%COMM,ERCODE)
      data%private%LORDER = LORDER
! Broadcast NEQSB, NCOLA,data%ENTRIES from root
      CALL MPI_BCAST(data%NEQSB,NBLOCK,MPI_INTEGER,0,data%private%COMM,ERCODE)
      CALL MPI_BCAST(data%private%NCOLA,NBLOCK,MPI_INTEGER,0, &
                     data%private%COMM,ERCODE)
      CALL MPI_BCAST(data%ENTRIES,NBLOCK,MPI_INTEGER,0,data%private%COMM, &
                     ERCODE)

! Prepare for analyse for each submatrix.
! We will do this in parallel and just share the submatrices
! out to the processes in order.
! Reset ICONTROL(10) if user has single process and asked for
! host NOT to be involved in submatrix factorizations (not possible!)
      IF (ICONTROL(10) == 1 .AND. NPROC == 1) ICONTROL(10) = 0

      IF (ICONTROL(10) == 0) THEN
! Host is involved in submatrix factorizations
        DO JBLOCK = 1,NBLOCK
          IPROC = JBLOCK - (JBLOCK/NPROC)*NPROC
! Add up number of submatrices process IPROC deals with
          data%ICOUNT(IPROC) = data%ICOUNT(IPROC) + 1
          data%INV_LIST(JBLOCK) = IPROC
        END DO
      ELSE IF (ICONTROL(10) == 1) THEN
! User does not want the host (process 0) to be involved in
! the subdomain factorizations
        IPROC = 0
        DO JBLOCK = 1,NBLOCK
          IF (IPROC.EQ.NPROC-1) IPROC = 0
          IPROC = IPROC + 1
          data%ICOUNT(IPROC) = data%ICOUNT(IPROC) + 1
          data%INV_LIST(JBLOCK) = IPROC
        END DO
      ELSE
! User has supplied INV_LIST
        DO  JBLOCK = 1,NBLOCK
           IPROC = data%INV_LIST(JBLOCK)
           data%ICOUNT(IPROC) = data%ICOUNT(IPROC) + 1
        END DO
      END IF

! Fill in data%LIST from the back.
      data%IPLIST(0) = data%ICOUNT(0) + 1
      DO IPROC = 1,NPROC-1
         data%IPLIST(IPROC) = data%IPLIST(IPROC-1) + &
                              data%ICOUNT(IPROC)
      END DO
      DO JBLOCK = 1,NBLOCK
         IPROC = data%INV_LIST(JBLOCK)
         data%IPLIST(IPROC) = data%IPLIST(IPROC) - 1
         J = data%IPLIST(IPROC)
         data%LIST(J) = JBLOCK
      END DO
      data%IPLIST(NPROC) = NBLOCK + 1
      DO IPROC = 0,NPROC-1
         IF (data%ICOUNT(IPROC) == 0) data%IPLIST(IPROC) = NBLOCK + 1
      END DO

! Initialise IBLOCK to hold list of submatrices for each process
      ISTRT = data%IPLIST(RANK)
      ISTOP = data%IPLIST(RANK+1) - 1
      LCNT = 1
      DO I = ISTRT,ISTOP
         JBLOCK = data%LIST(I)
         data%IBLOCK(LCNT) = JBLOCK
         LCNT = LCNT + 1
      END DO

      CALL MPI_BARRIER(data%private%COMM,ERCODE)

      data%private%ERROR = data%ERROR
      data%private%NEQ = NEQ
      data%private%NBLOCK = NBLOCK
      IF (RANK == 0) data%private%LGUARD = LGUARD
      IF (data%JOB == 2) GO TO 950

! ==================================
! Ready for analyse
! ==================================

  200 CONTINUE
! At this point JOB = 3 or 23 or 24 or 25
      CALL MP48BD(data)
      IF (data%JOB == 3 .OR. data%JOB == 23) GO TO 950

! ==================================
! Ready to prepare for factorization
! ==================================

  300 CONTINUE
! At this point JOB = 4 or 24 or 25
      CALL MP48CD(data)
      IF (data%JOB == 4 .OR. data%JOB == 24) GO TO 950

! ==================================
! Solve for right-hand sides
! ==================================
! Right-hand sides are stored in data%B(1:NEQ)

  400 CONTINUE
! At this point JOB = 5 or 25
      CALL MP48DD(data)
      GO TO 950

! ==================================
! General error return
! ==================================
  890 CONTINUE
      CALL MPI_BCAST(data%ERROR,1,MPI_INTEGER,0,data%private%COMM,ERCODE)
      data%private%ERROR = data%ERROR

! ==========================
! Delete files
! ==========================
! Delete files for submatrix JBLOCK
! (we can set ICONTROL(12) to be nonzero if files to be saved
! but beware that files local to processor)

  900 CONTINUE
      FLAG = 0
      IF (ICONTROL(11) == 0) GO TO 940
! Will not try to delete files if an error was found with them
      IF (data%ERROR == -25) GO TO 940
      IF (ICONTROL(12) == 0) THEN
        IF (allocated(data%ICOUNT) .AND. allocated(data%FILES)) THEN
! We must delete files used
          DEALLOCATE (IFLAG,STAT=ST)
          ALLOCATE (IFLAG(NPROC),STAT=ST)
! Find a unit number
          DO STRM = 8,99
            IF (STRM == DGN .OR. STRM == ERR .OR. STRM == WRN) &
                CYCLE
            INQUIRE (UNIT=STRM,IOSTAT=IOS,ERR=910,EXIST=EX, &
                     OPENED=OPEN)
            IF (EX .AND. .NOT.OPEN) GO TO 910
          END DO
! No unit found. Jump to error return.
          FLAG = -20
! Gather error status from on root
  910     CALL MPI_BARRIER(data%private%COMM,ERCODE)
          IF (IOS /= 0) FLAG = -14
          data%STAT = ST
          CALL MP48ND(data,NPROC,IFLAG,FLAG)
          data%IOSTAT = IOS
! Broadcast error status from root
          CALL MPI_BCAST(data%ERROR,1,MPI_INTEGER,0,data%private%COMM,ERCODE)
! Terminate if error occurred
          IF (data%ERROR == -14 .OR. data%ERROR == -20) GO TO 940

          DO LCNT = 1,data%ICOUNT(RANK)
            JBLOCK = data%IBLOCK(LCNT)
            CLOSE (UNIT=STRM,STATUS='DELETE')
            OPEN (UNIT=STRM,IOSTAT=IOS,STATUS='REPLACE',ERR=300, &
                  FILE=data%FILES(1,JBLOCK),FORM='UNFORMATTED')
            CLOSE (UNIT=STRM,STATUS='DELETE')
            OPEN (UNIT=STRM,IOSTAT=IOS,STATUS='REPLACE',ERR=300, &
                  FILE=data%FILES(2,JBLOCK),FORM='UNFORMATTED')
            CLOSE (UNIT=STRM,STATUS='DELETE')
          END DO

! On root, must also delete interface files
          IF (data%RANK == 0) THEN
            CLOSE (UNIT=STRM,STATUS='DELETE')
            OPEN (UNIT=STRM,IOSTAT=IOS,STATUS='REPLACE',ERR=300, &
                  FILE=data%FILES(1,data%NBLOCK+1),FORM='UNFORMATTED')
            CLOSE (UNIT=STRM,STATUS='DELETE')
            OPEN (UNIT=STRM,IOSTAT=IOS,STATUS='REPLACE',ERR=300, &
                  FILE=data%FILES(2,data%NBLOCK+1),FORM='UNFORMATTED')
            CLOSE (UNIT=STRM,STATUS='DELETE')
          END IF
        END IF

      END IF
  940 CONTINUE

      IF (data%JOB /= 6) GO TO 950

      DEALLOCATE (data%NGUARD,STAT=ST)
      DEALLOCATE (data%ICOUNT,STAT=ST)
      DEALLOCATE (data%LIST,STAT=ST)
      DEALLOCATE (data%IBLOCK,STAT=ST)
      DEALLOCATE (data%IPLIST,STAT=ST)
      IF (data%ICNTL(10) == 0) DEALLOCATE (data%INV_LIST,STAT=ST)
      DEALLOCATE (data%ENTRIES,STAT=ST)
      DEALLOCATE (data%FACT,STAT=ST)
      DEALLOCATE (data%IRNF,STAT=ST)
      DEALLOCATE (data%OPS,STAT=ST)
      DEALLOCATE (data%OPSA,STAT=ST)
      DEALLOCATE (data%STORAGE,STAT=ST)
      DEALLOCATE (data%INFO_MA60,STAT=ST)
! 23.01.03  Bug pointed out by Iain Strachen. data%NEQSB is
! allocated by the user on the host so we do not deallocate it here.
!     DEALLOCATE (data%NEQSB,STAT=ST)
      IF (data%RANK /= 0) DEALLOCATE (data%NEQSB,STAT=ST)
      IF (data%RANK /= 0) DEALLOCATE (data%X,STAT=ST)
      IF (data%RANK == 0 .AND. data%ICNTL(13) == 0) &
       DEALLOCATE (data%X,STAT=ST)
      DEALLOCATE (data%INTER_COL,STAT=ST)
      DEALLOCATE (data%TIMEA,STAT=ST)
      DEALLOCATE (data%TIMEF,STAT=ST)
      DEALLOCATE (data%TIMES,STAT=ST)
      DEALLOCATE (data%private%IDUP,STAT=ST)
      IF (data%RANK == 0) DEALLOCATE (data%private%IGLOB,STAT=ST)
      DEALLOCATE (data%private%JPTR,STAT=ST)
      DEALLOCATE (data%private%VPTR,STAT=ST)
      DEALLOCATE (data%private%IPTR,STAT=ST)
      DEALLOCATE (data%private%LPTR,STAT=ST)
      DEALLOCATE (data%private%LFACT,STAT=ST)
      DEALLOCATE (data%private%LIRNF,STAT=ST)
      DEALLOCATE (data%private%NCOLA,STAT=ST)
      DEALLOCATE (data%private%NF,STAT=ST)
      DEALLOCATE (data%private%NP,STAT=ST)
      DEALLOCATE (data%private%NEQSUB,STAT=ST)
      DEALLOCATE (data%private%SUBPTR,STAT=ST)
      DEALLOCATE (data%private%IMAP,STAT=ST)
      DEALLOCATE (data%private%JFVAR,STAT=ST)
      DEALLOCATE (data%private%GLOBAL,STAT=ST)
      DEALLOCATE (data%private%CGLOB,STAT=ST)
      DEALLOCATE (data%private%IP,STAT=ST)
      DEALLOCATE (data%private%IQ,STAT=ST)
      DEALLOCATE (data%private%IQ_COPY,STAT=ST)
      DEALLOCATE (data%private%IQM_COPY,STAT=ST)
      DEALLOCATE (data%private%INULL,STAT=ST)
      DEALLOCATE (data%private%IRN,STAT=ST)
      DEALLOCATE (data%private%JCN,STAT=ST)
      DEALLOCATE (data%private%A,STAT=ST)
      DEALLOCATE (data%private%KEEP48,STAT=ST)
      DEALLOCATE (data%private%ICNTL_60,STAT=ST)
      DEALLOCATE (data%private%CNTL_60,STAT=ST)
      DEALLOCATE (data%private%IPTRL,STAT=ST)
      DEALLOCATE (data%private%IPTRU,STAT=ST)

      IF(data%private%COMM /= data%COMM) &
         CALL MPI_COMM_FREE(data%private%COMM,ERCODE)

  950 CONTINUE
      DEALLOCATE (W,STAT=ST)

      DEALLOCATE (IGUARD,STAT=ST)
      DEALLOCATE (IFLAG,STAT=ST)
      DEALLOCATE (IW,STAT=ST)
      DEALLOCATE (IWORK,STAT=ST)
      DEALLOCATE (LOAD,STAT=ST)

      RETURN

! =========================
! Error and warning formats
! =========================
 9000 FORMAT (/' Error return with data%JOB = ',I2,'. ', &
                'data%ERROR = ',I4)
 9010 FORMAT (' ',A,' has not been allocated')
 9020 FORMAT (/' Error return with data%ERROR = ',I4)
 9080 FORMAT (' Value of ', A,' out of range. ', A, ' = ',I6)
 9090 FORMAT (' ',A,' = ',I10,' for submatrix ',I3)
 9100 FORMAT (' ',A,' is of size ',I10,' Increase to at least ',I10)
 9120 FORMAT (' Error returned by previous call')
 9130 FORMAT (' Call follows a call with data%JOB = ',I3)
 9140 FORMAT (' ICNTL    : ', 5(I10,1X),/, &
              '            ', 5(I10,1X),/, &
              '            ', 3(I10,1X))
 9150 FORMAT (' CNTL     : ', 4(1PES12.5,1X))
 9160 FORMAT (' On first call with data%JOB = 4', &
        ' must have data%FACT_JOB = 1 or 2')
      END SUBROUTINE MP48AD
!*************************************************

      SUBROUTINE MP48BD(data)

! This subroutine performs analyse phase using MP48_MA60A/AD (JOB=3)
      USE HSL_MP48_DATA_DOUBLE
      IMPLICIT NONE

      TYPE (MP48_DATA) data

!     .. Parameters ..
      REAL(WP), PARAMETER :: ONE = 1.0_WP
      REAL(WP), PARAMETER :: ZERO = 0.0_WP
      REAL(DP), PARAMETER :: ZEROD = 0.0_DP

!     .. Local Scalars ..

      REAL(DP) T1,T2
      REAL(WP) RINFO
      INTEGER BLANK,DEST,DGN,ERCODE,ERR, &
              FLAG, &
              I,IFLAG46,IPROC,ISTRT,ISTOP,ICOUNT_RANK,IDUP,IDUP1, &
              J,J1,J2,J1V,J2V,JBLOCK,JJ, &
              K,KSTART,KSTOP, &
              L,L1,L2,LA,LCNT,LDGN,LORDER,LMX, &
              M,MMX, &
              NBLOCK,NEQ,NPROC,NC,NZBLK,NCMX,NZMX,NVAR,NP, &
              RANK,STRM,SRC,VALCRD,WRN

      INTEGER(8) MYLMX
      INTEGER IOS        !  Holds IOSTAT parameter for files.
      INTEGER ST         !  Holds STAT parameter (ALLOCATE statement)
      INTEGER STD        !  Holds STAT parameter (DEALLOCATE statement)
      LOGICAL LERR,LWRN
      LOGICAL EX         !  Used by INQUIRE. Set to .TRUE. if
                         !  unit exists.
      LOGICAL OPEN       !  Used by INQUIRE. Set to .TRUE. if
                         !  file open.
!     .. Local Arrays ..

      REAL(WP) CONTROL(30)
      REAL(WP), DIMENSION (:), ALLOCATABLE :: A
      REAL(WP), DIMENSION (:), ALLOCATABLE :: VALUES

      CHARACTER,DIMENSION(:), ALLOCATABLE :: VALNAM*128
!     CHARACTER,DIMENSION(:), allocatable :: VALNAM*128
      CHARACTER (LEN = 128) :: TEMP_VALNAM

      INTEGER ICONTROL(30),LEN(2),STAT(MPI_STATUS_SIZE,20)

      INTEGER, DIMENSION (:), ALLOCATABLE :: COLPTR
      INTEGER, DIMENSION (:), ALLOCATABLE :: IFLAG
      INTEGER, DIMENSION (:), ALLOCATABLE :: IFIRST
      INTEGER, DIMENSION (:), ALLOCATABLE :: IW
      INTEGER, DIMENSION (:), ALLOCATABLE :: IW1
      INTEGER(8), DIMENSION (:), ALLOCATABLE :: MYIW1
      INTEGER, DIMENSION (:), ALLOCATABLE :: IRN
      INTEGER, DIMENSION (:), ALLOCATABLE :: JCN
      INTEGER, DIMENSION (:), ALLOCATABLE :: JFIRST
      INTEGER, DIMENSION (:), ALLOCATABLE :: JPTR
      INTEGER, DIMENSION (:), ALLOCATABLE :: LPTR
      INTEGER, DIMENSION (:), ALLOCATABLE :: LASTC
      INTEGER, DIMENSION (:), ALLOCATABLE :: LASTR
      INTEGER, DIMENSION (:), ALLOCATABLE :: LENC
      INTEGER, DIMENSION (:), ALLOCATABLE :: LENR
      INTEGER, DIMENSION (:), ALLOCATABLE :: NEQSUB
      INTEGER, DIMENSION (:), ALLOCATABLE :: NVARSB
      INTEGER, DIMENSION (:), ALLOCATABLE :: NEXTC
      INTEGER, DIMENSION (:), ALLOCATABLE :: NEXTR
      INTEGER, DIMENSION (:), ALLOCATABLE :: ROWPTR
      INTEGER, DIMENSION (:), ALLOCATABLE :: VPTR
      INTEGER, DIMENSION (:,:), ALLOCATABLE :: SBUF

      LOGICAL, DIMENSION (:), ALLOCATABLE :: LMP48_MA60

!     .. External Subroutines ..
      EXTERNAL MC46AD

!     .. MPI routines
      EXTERNAL MPI_BARRIER,MPI_BCAST,MPI_COMM_RANK, &
               MPI_GATHER,MPI_SEND,MPI_RECV,MPI_REDUCE

!     .. Intrinsic Functions ..
      INTRINSIC MAX

      FLAG = 0
      RANK = data%RANK
      NPROC = data%NPROC

! Set diagnostics controls
      DGN = data%ICNTL(3)
      LDGN = data%ICNTL(4)
      IF (DGN < 0) LDGN = 0
! Set error controls
      LERR = .FALSE.
      IF (data%ICNTL(1).GE.0 .AND. data%ICNTL(4).GE.1) LERR = .TRUE.
      ERR = data%ICNTL(1)
! Set warning controls
      LWRN = .FALSE.
      IF (data%ICNTL(2).GE.0 .AND. data%ICNTL(4).GE.1) LWRN = .TRUE.
      WRN = data%ICNTL(2)

! ICONTROL(5) controls full matrix processing.
      ICONTROL(5) =  MAX(0,data%ICNTL(5))
! ICONTROL(6) controls increasing space when factorization space too
! small.
      ICONTROL(6) =  data%ICNTL(6)
      IF (data%ICNTL(6).LE.0) ICONTROL(6) = 2
! ICONTROL(7) controls input of matrices.
      ICONTROL(7) =  data%ICNTL(7)
! ICONTROL(8) controls number of cols. searched for pivots
      ICONTROL(8) =  data%ICNTL(8)
! ICONTROL(9) controls action if duplicates found
      ICONTROL(9) =  data%ICNTL(9)

      CONTROL(1) =  MAX(ZERO,data%CNTL(1))
      CONTROL(1) =  MIN(ONE,CONTROL(1))
      CONTROL(2) =  MAX(ZERO,data%CNTL(2))
      CONTROL(2) =  MIN(ONE,CONTROL(2))
      CONTROL(3) =  MAX(ZERO,data%CNTL(3))
      CONTROL(4) =  MAX(ZERO,data%CNTL(4))

      IF (RANK == 0) THEN
         NBLOCK = data%private%NBLOCK
         NEQ    = data%private%NEQ

! Check data for errors. Check required arrays have been allocated.
         IF (ICONTROL(7) == 1 .OR. ICONTROL(7) == 2) THEN
! Check size of array VALNAM
! First check it has been allocated
            IF (.NOT. allocated(data%VALNAM)) THEN
               FLAG = -10
            ELSE
               IF (SIZE(data%VALNAM) < data%NBLOCK) FLAG = -10
            END IF
            IF (FLAG == -10) THEN
               data%ERROR = -10
               IF (LERR) THEN
                  WRITE (ERR,FMT=9000) data%JOB,data%ERROR
                  IF (.NOT. allocated(data%VALNAM)) THEN
                     WRITE (ERR,FMT=9020) ' data%VALNAM'
                  ELSE IF (SIZE(data%VALNAM) < data%NBLOCK) THEN
                     WRITE (ERR,FMT=9100) ' data%VALNAM', &
                     SIZE(data%VALNAM),data%NBLOCK
                  END IF
               END IF
               GO TO 20
            END IF
         ELSE IF (ICONTROL(7) == 3) THEN
! Check size of array VALUES
            IF (.NOT. allocated(data%VALUES)) THEN
               FLAG = -10
            ELSE
               VALCRD = SIZE(data%VALUES)
               IF (VALCRD < data%NE) FLAG = -10
            END IF
            IF (FLAG == -10) THEN
               data%ERROR = -10
               IF (LERR) THEN
                  WRITE (ERR,FMT=9000) data%JOB,data%ERROR
                  IF (.NOT. allocated(data%VALUES)) THEN
                     WRITE (ERR,FMT=9020) ' data%VALUES'
                  ELSE IF (VALCRD < data%NE) THEN
                     WRITE (ERR,FMT=9100) &
                  ' data%VALUES',VALCRD,data%NE
                  END IF
               END IF
               GO TO 20
            END IF
         END IF

      END IF
! Broadcast error status from root
   20 CALL MPI_BCAST(data%ERROR,1,MPI_INTEGER,0,data%private%COMM,ERCODE)
! Terminate if error occurred
      IF (data%ERROR < 0) GO TO 400

      NBLOCK = data%NBLOCK
      NEQ = data%NEQ

      DEALLOCATE (IFLAG,STAT=ST)
      ALLOCATE (IFLAG(NPROC),STAT=ST)
! Initialise the error flag (so that it is always well-defined)
      IF (ST == 0) IFLAG(1:NPROC) = 0

      FLAG = 0
      IF (ICONTROL(7) == 4) THEN
! On each process, check array VALUES
         IF (.NOT. allocated(data%VALUES)) THEN
            FLAG = -10
         ELSE
            VALCRD = SIZE(data%VALUES)
            IF (VALCRD < data%NE) FLAG = -11
         END IF
      ELSE IF (ICONTROL(7) == 5) THEN
         IF (.NOT. allocated(data%RVAL)) THEN
            FLAG = -10
         ELSE
            VALCRD = SIZE(data%RVAL)
! Check length of RVAL. Loop over submatrices assigned
! to process
            L = 0
            DO LCNT = 1,data%ICOUNT(RANK)
               JBLOCK = data%IBLOCK(LCNT)
               L = L + data%ENTRIES(JBLOCK)
            END DO
            IF (VALCRD < L) FLAG = -10
         END IF
      END IF
      IF (ICONTROL(7) == 4 .OR. ICONTROL(7) == 5) THEN
         CALL MPI_GATHER(FLAG,1,MPI_INTEGER,IFLAG,1,MPI_INTEGER,0, &
                      data%private%COMM,ERCODE)
! Check for error on root
         IF (RANK == 0) THEN
            DO 21 IPROC = 1,NPROC
               IF (IFLAG(IPROC) /= 0) THEN
                  data%ERROR = -10
                  IF (LERR) THEN
                     WRITE (ERR,FMT=9000) data%JOB,data%ERROR
                     IF (ICONTROL(7) == 4) THEN
                       IF (IFLAG(IPROC) == -10) THEN
                         WRITE (ERR,FMT=9020) 'data%VALUES'
                       ELSE
                         WRITE (ERR,FMT=9110) 'data%VALUES'
                       END IF
                     ELSE
                       IF (IFLAG(IPROC) == -10) THEN
                         WRITE (ERR,FMT=9020) 'data%RVAL'
                       ELSE
                         WRITE (ERR,FMT=9110) 'data%RVAL'
                       END IF
                     END IF
                  END IF
                  GO TO 22
               END IF
   21       CONTINUE
         END IF
      END IF
! Broadcast error status from root
   22 CALL MPI_BCAST(data%ERROR,1,MPI_INTEGER,0,data%private%COMM,ERCODE)
! Terminate if error occurred
      IF (data%ERROR < 0) GO TO 400

      ST = 0
      IOS = 0
      ICOUNT_RANK = data%ICOUNT(RANK)
      DEALLOCATE (JPTR,LPTR,LMP48_MA60,NEQSUB,VALNAM,VPTR, &
                  NVARSB,SBUF,STAT=ST)
      DEALLOCATE (data%private%IPTR,data%private%JPTR, &
                  data%private%LPTR,data%private%VPTR, &
                  data%private%NEQSUB,STAT=ST)
      ALLOCATE (JPTR(1:NBLOCK), &
                LPTR(1:NBLOCK), &
                LMP48_MA60(1:NBLOCK), &
                NEQSUB(1:NBLOCK), &
                VALNAM(1:ICOUNT_RANK), &
                VPTR(1:NBLOCK), &
                data%private%IPTR(1:NBLOCK), &
                data%private%JPTR(1:ICOUNT_RANK), &
                data%private%LPTR(1:ICOUNT_RANK), &
                data%private%VPTR(1:ICOUNT_RANK), &
                data%private%NEQSUB(1:ICOUNT_RANK), &
                NVARSB(1:ICOUNT_RANK), &
                SBUF(1:30,1:NBLOCK), &
                STAT=ST)
      IF (ST /= 0) GO TO 30
! Ensure VPTR always defined
      VPTR(1:NBLOCK) = 1

      LORDER = data%private%LORDER
      IF (RANK /= 0) THEN
         DEALLOCATE (data%private%GLOBAL,STAT=ST)
         ALLOCATE (data%private%GLOBAL(1:LORDER,1:NBLOCK), &
                   STAT=ST)
         IF (ST /= 0) GO TO 30
      END IF

! Initialise LMP48_MA60 to FALSE (will be set to TRUE when MP48_MA60A/AD
! has been called successfully)
      LMP48_MA60(1:NBLOCK) = .FALSE.

! Find suitable unit number if we are to read data from file
! Ensure STRM is defined
      IOS = 0
      FLAG = 0
      STRM = 1
      IF (ICONTROL(7).GE.3) GO TO 30
      IF (ICONTROL(7) == 2 .AND. RANK /= 0) GO TO 30
      DO STRM = 8,99
         IF (STRM == DGN .OR. STRM == ERR .OR. STRM == WRN) &
             CYCLE
         INQUIRE (UNIT=STRM,IOSTAT=IOS,ERR=30,EXIST=EX, &
                  OPENED=OPEN)
         IF (EX .AND. .NOT.OPEN) GO TO 25
      END DO
! No unit found. Jump to error return.
      FLAG = -20
      GO TO 30
   25 CONTINUE

! Gather error status from on root
   30 CALL MPI_BARRIER(data%private%COMM,ERCODE)
      IF (IOS /= 0) FLAG = -14
      data%STAT = ST
      CALL MP48ND(data,NPROC,IFLAG,FLAG)
      data%IOSTAT = IOS
! Broadcast error status from root
      CALL MPI_BCAST(data%ERROR,1,MPI_INTEGER,0,data%private%COMM,ERCODE)
! Terminate if error occurred
      IF (data%ERROR < 0) GO TO 400

! IPTR points to beginning of variable list for submatrix JBLOCK.
      data%private%IPTR(1) = 1
      DO JBLOCK = 1,NBLOCK-1
         data%private%IPTR(JBLOCK+1) = data%private%IPTR(JBLOCK) + &
                                       data%ENTRIES(JBLOCK)
      END DO

      FLAG = 0
      NZMX = 0
      NCMX = 0
      MMX  = 0
! Send data from root to other processes.
      IF (RANK == 0) THEN

        IF (ICONTROL(7) == 2) THEN
! Find length needed for VALUES
           L = 1
           DO JBLOCK = 1,NBLOCK
              L = MAX(L,data%ENTRIES(JBLOCK))
           END DO
           DEALLOCATE (VALUES,STAT=ST)
           ALLOCATE (VALUES(1:L),STAT=ST)
           IF (ST /= 0) GO TO 190
         END IF

         DO 80 IPROC = 1,NPROC-1
            DEST = IPROC
            ISTRT = data%IPLIST(IPROC)
            ISTOP = data%IPLIST(IPROC+1) - 1
            LCNT = 1
            DO 70 I = ISTRT,ISTOP
               JBLOCK = data%LIST(I)
! Send data for rows in submatrix JBLOCK.

               CALL MPI_SEND(data%EQPTR(data%private%SUBPTR(JBLOCK)), &
                             data%NEQSB(JBLOCK)+1,MPI_INTEGER,DEST,3, &
                             data%private%COMM,ERCODE)
               CALL MPI_SEND(data%EQVAR(data%private%IPTR(JBLOCK)), &
                             data%ENTRIES(JBLOCK),MPI_INTEGER,DEST,4, &
                             data%private%COMM,ERCODE)
               CALL MPI_SEND(data%private%GLOBAL(1,JBLOCK), &
                             data%NEQSB(JBLOCK), &
                             MPI_INTEGER,DEST,5,data%private%COMM,ERCODE)

               IF (ICONTROL(7) == 3) THEN
                  CALL MPI_SEND(data%VALUES(data%private%IPTR(JBLOCK)), &
                                data%ENTRIES(JBLOCK), &
                                MPI_DOUBLE_PRECISION, &
                                DEST,6,data%private%COMM,ERCODE)

               ELSE IF (ICONTROL(7) == 2) THEN
                  OPEN (UNIT=STRM,IOSTAT=IOS,ERR=190, &
                        FILE=data%VALNAM(JBLOCK),FORM='UNFORMATTED')
                  READ (STRM) VALUES(1:data%ENTRIES(JBLOCK))
                  CLOSE (UNIT=STRM,STATUS='KEEP')

                  CALL MPI_SEND(VALUES,data%ENTRIES(JBLOCK), &
                                MPI_DOUBLE_PRECISION, &
                                DEST,6,data%private%COMM,ERCODE)

               ELSE IF (ICONTROL(7) == 1) THEN
                  TEMP_VALNAM = ADJUSTL(data%VALNAM(JBLOCK))
                  BLANK = LEN_TRIM(TEMP_VALNAM)
                  CALL MPI_SEND(BLANK,1,MPI_INTEGER,DEST,7, &
                                data%private%COMM,ERCODE)
                  CALL MPI_SEND(TEMP_VALNAM,BLANK, &
                                MPI_CHARACTER,DEST,8,data%private%COMM,ERCODE)

               END IF
               LCNT = LCNT + 1

   70       CONTINUE
   80    CONTINUE
         IF (ICONTROL(7) == 2) DEALLOCATE (VALUES,STAT=STD)


! Set up parameters for submatrices that root process will deal with
         DO LCNT = 1,ICOUNT_RANK
            JBLOCK = data%IBLOCK(LCNT)
            NEQSUB(LCNT) = data%NEQSB(JBLOCK)
            NVARSB(LCNT) = data%ENTRIES(JBLOCK)
            LPTR(LCNT)   = data%private%SUBPTR(JBLOCK)
            JPTR(LCNT)   = data%private%IPTR(JBLOCK)
            M = NEQSUB(LCNT)
            NC = data%private%NCOLA(JBLOCK)
            NZMX = MAX(NZMX,NVARSB(LCNT))
            NCMX = MAX(NCMX,NC)
            MMX  = MAX(MMX,M)
         END DO

         IF (ICONTROL(7) == 1) THEN
            DO LCNT = 1,ICOUNT_RANK
               JBLOCK = data%IBLOCK(LCNT)
               TEMP_VALNAM = ADJUSTL(data%VALNAM(JBLOCK))
               BLANK = LEN_TRIM(TEMP_VALNAM)
               VALNAM(LCNT) = TEMP_VALNAM(1:BLANK)
            END DO

         ELSE IF (ICONTROL(7) == 3 .OR. ICONTROL(7) == 4) THEN
            DO LCNT = 1,ICOUNT_RANK
               JBLOCK = data%IBLOCK(LCNT)
! VPTR points to beginning of value list for submatrix JBLOCK.
               VPTR(LCNT) = data%private%IPTR(JBLOCK)
            END DO

         ELSE IF (ICONTROL(7) == 5) THEN
            DO LCNT = 1,ICOUNT_RANK
               JBLOCK = data%IBLOCK(LCNT)
               IF (LCNT == 1) THEN
                  VPTR(LCNT) = 1
               ELSE
                  VPTR(LCNT) = VPTR(LCNT-1) + NVARSB(LCNT-1)
               END IF
            END DO

         END IF

! If the root is NOT dealing with submatrix 1, then the data
! in data%EQPTR, data%EQVAR, data%VALUES will not be at the start

      ELSE

         LEN(1:2) = 0
         DO LCNT = 1,ICOUNT_RANK
            JBLOCK = data%IBLOCK(LCNT)
            NEQSUB(LCNT) = data%NEQSB(JBLOCK)
            NVARSB(LCNT) = data%ENTRIES(JBLOCK)
! Add up how much space needed
            LEN(1) = LEN(1) + NEQSUB(LCNT) + 1
            LEN(2) = LEN(2) + NVARSB(LCNT)
            IF (LCNT == 1) THEN
               LPTR(LCNT) = 1
               JPTR(LCNT) = 1
            ELSE
               LPTR(LCNT) = LPTR(LCNT-1) + NEQSUB(LCNT-1) + 1
               JPTR(LCNT) = JPTR(LCNT-1) + NVARSB(LCNT-1)
            END IF
            VPTR(LCNT) = JPTR(LCNT)
            IF (ICONTROL(7) == 4) VPTR(LCNT) = data%private%IPTR(JBLOCK)
         END DO


         DEALLOCATE (data%EQPTR,data%EQVAR,STAT=ST)
         ALLOCATE (data%EQPTR(1:LEN(1)), &
                   data%EQVAR(1:LEN(2)),STAT=ST)
         IF (ST /= 0) GO TO 190

         IF (ICONTROL(7) == 2 .OR. ICONTROL(7) == 3) THEN
! We need space to receive values sent from root
            DEALLOCATE (data%VALUES,STAT=ST)
            ALLOCATE (data%VALUES(1:LEN(2)),STAT=ST)
            IF (ST /= 0) GO TO 190
         END IF

         DO 100 LCNT = 1,ICOUNT_RANK

            JBLOCK = data%IBLOCK(LCNT)
! Receive data from the root
            CALL MPI_RECV(data%EQPTR(LPTR(LCNT)),NEQSUB(LCNT)+1, &
                          MPI_INTEGER,0,3,data%private%COMM,STAT,ERCODE)

            CALL MPI_RECV(data%EQVAR(JPTR(LCNT)),NVARSB(LCNT), &
                          MPI_INTEGER,0,4,data%private%COMM,STAT,ERCODE)


            CALL MPI_RECV(data%private%GLOBAL(1,JBLOCK),NEQSUB(LCNT), &
                          MPI_INTEGER,0,5,data%private%COMM,STAT,ERCODE)

            IF (ICONTROL(7) == 2 .OR. ICONTROL(7) == 3) THEN
               CALL MPI_RECV(data%VALUES(VPTR(LCNT)),NVARSB(LCNT), &
                             MPI_DOUBLE_PRECISION, &
                             0,6,data%private%COMM,STAT,ERCODE)

            ELSE IF (ICONTROL(7) == 1) THEN
               CALL MPI_RECV(BLANK,1,MPI_INTEGER,0,7, &
                             data%private%COMM,STAT,ERCODE)
               CALL MPI_RECV(VALNAM(LCNT),BLANK, &
                             MPI_CHARACTER,0,8,data%private%COMM,STAT,ERCODE)
               VALNAM(LCNT) = VALNAM(LCNT)(1:BLANK)

            END IF
            M = NEQSUB(LCNT)
            NC = data%private%NCOLA(JBLOCK)
            NZMX = MAX(NZMX,NVARSB(LCNT))
            NCMX = MAX(NCMX,NC)
            MMX  = MAX(MMX,M)

  100    CONTINUE

      END IF
      DEALLOCATE (data%private%SUBPTR,STAT=STD)


      DEALLOCATE (COLPTR,IFIRST,IW,IW1,JFIRST,LASTC,LASTR,LENC,LENR, &
                  NEXTC,NEXTR,ROWPTR, &
                  STAT=ST)
      ALLOCATE (COLPTR(1:NCMX+1), &
                IFIRST(1:NCMX), &
                IW(1:NEQ), &
                IW1(1:NEQ), &
                JFIRST(1:MMX), &
                LASTC(1:NCMX), &
                LENR(1:MMX), &
                LENC(1:NCMX), &
                NEXTC(1:NCMX), &
                ROWPTR(1:LORDER+1), &
                STAT=ST)
      DEALLOCATE (MYIW1,STAT=ST)
      ALLOCATE (MYIW1(1:NEQ),STAT=ST)
      IF (ST /= 0) GO TO 190
      IF (ICONTROL(8) == 0) THEN
                ALLOCATE (LASTR(1:MMX), &
                NEXTR(1:MMX), &
                STAT=ST)
      ELSE
                ALLOCATE (LASTR(1), &
                NEXTR(1), &
                STAT=ST)
      END IF
      IF (ST /= 0) GO TO 190
      DEALLOCATE (data%OPS,data%OPSA, &
                  data%INFO_MA60, &
                  data%TIMEA, &
                  data%TIMEF, &
                  data%private%CNTL_60, &
                  data%private%ICNTL_60, &
                  data%private%IDUP, &
                  data%private%IMAP, &
                  data%private%IQ, &
                  data%private%IP, &
                  data%private%NF, &
                  data%private%NP, &
                  STAT=ST)
      ALLOCATE (data%OPS(1:NBLOCK),data%OPSA(1:NBLOCK), &
                data%INFO_MA60(1:20,1:NBLOCK), &
                data%TIMEA(1:NBLOCK), &
                data%TIMEF(1:NBLOCK), &
                data%private%CNTL_60(1:4,ICOUNT_RANK), &
                data%private%ICNTL_60(1:10,ICOUNT_RANK), &
                data%private%IDUP(1:ICOUNT_RANK), &
                data%private%IMAP(1:NZMX,1:ICOUNT_RANK), &
                data%private%IQ(1:NCMX+1,1:ICOUNT_RANK), &
                data%private%IP(1:LORDER,1:ICOUNT_RANK), &
                data%private%NF(1:NBLOCK), &
                data%private%NP(1:NBLOCK), &
                STAT=ST)
      IF (ST /= 0) GO TO 190

      data%INFO_MA60(1:20,1:NBLOCK) = 0
      data%OPS(1:NBLOCK) = ZERO
      data%OPSA(1:NBLOCK) = ZERO
      data%TIMEA(1:NBLOCK) = ZEROD
      data%TIMEF(1:NBLOCK) = ZEROD

! Gather error status from allocations
  190 CALL MPI_BARRIER(data%private%COMM,ERCODE)
      data%STAT = ST
      data%IOSTAT = IOS
      CALL MP48ND(data,NPROC,IFLAG,FLAG)
! Broadcast error status from root
      CALL MPI_BCAST(data%ERROR,1,MPI_INTEGER,0,data%private%COMM,ERCODE)
! Terminate if error occurred
      IF (data%ERROR < 0) GO TO 400

! RECALL: EQVAR holds local column indices, with the cols in the border
! of each block ordered last.
! For MP48_MA60A/AD we need row indices and column allocatables for
! each block. COLPTR will hold col. allocatables. MC46 will be used
! for this, but first squeeze out any duplicate entries

      FLAG = 0
      IDUP = 0
      MYIW1(1:NEQ) = 0
      MYLMX = 0
! Look for duplicates and set up mapping array
      DO 240 LCNT = 1,ICOUNT_RANK
         IDUP1 = 0
         MYIW1(1:NEQ) = 0
         data%private%IMAP(1:NVARSB(LCNT),LCNT) = 0
         JBLOCK = data%IBLOCK(LCNT)
! M is number of rows in the block
         M = NEQSUB(LCNT)
! L1 points to beginning of col. indices for block JBLOCK
         L1 = LPTR(LCNT)
! J1 points to first entry in EQVAR for JBLOCK
         J1 = JPTR(LCNT)
         KSTART = 1
         J = 0
         JJ = 0
! Loop over the rows in the block
         DO 230 L = 1,M
           NVAR = data%EQPTR(L1+L) - data%EQPTR(L1+L-1)
           KSTOP = KSTART + NVAR
           DO 220 K = KSTART,KSTOP - 1
             I = data%EQVAR(J1+K-1)
             !if(L1+L.EQ.7071819)print *, "L ", L1+L,"I ",I,"indx",J1+K-1,"iw ",MYIW1(I),"llmx ",L+MYLMX
             J = J + 1
             IF (MYIW1(I) < L+MYLMX) THEN
               JJ = JJ + 1
               MYIW1(I) = L + MYLMX
             ELSE
! We have a duplicate in row L
               IDUP1 = IDUP1 + 1
               !if(L1+L.EQ.7071819) print *, "dub in ", L1+L,"I ",I
             END IF
! What was Jth entry in the submatrix maps to JJ-th entry
! (if we have duplicates more than one entry maps to same entry;
! if no duplicates J = JJ)
             data%private%IMAP(J,LCNT) = JJ
  220      CONTINUE
           KSTART = KSTOP
  230    CONTINUE
         MYLMX = MYLMX + NEQ
         data%private%IDUP(LCNT) = IDUP1
         IDUP = IDUP + IDUP1
  240 CONTINUE

      CALL MPI_BARRIER(data%private%COMM,ERCODE)
! On each process, IDUP is now the number of duplicate entries found.
! We must sum these to get total number of duplicates.
      CALL MPI_REDUCE(IDUP,data%IDUP,1,MPI_INTEGER,MPI_SUM,0, &
                      data%private%COMM,ERCODE)
      CALL MPI_BCAST(data%IDUP,1,MPI_INTEGER,0,data%private%COMM,ERCODE)
! Issue warning if duplicates found
      IF (data%RANK == 0 .AND. data%IDUP > 0) THEN
         IF (ICONTROL(9) == 0) THEN
            data%ERROR = data%ERROR + 1
            IF (LWRN) THEN
               WRITE (WRN,FMT=9010) data%JOB,data%ERROR
               WRITE (WRN,FMT=9060) data%IDUP
            END IF
         ELSE
            data%ERROR = -3
            IF (LERR) THEN
               WRITE (ERR,FMT=9000) data%JOB,data%ERROR
               WRITE (ERR,FMT=9060) data%IDUP
            END IF
         END IF
      END IF

! Broadcast error status from root
      CALL MPI_BCAST(data%ERROR,1,MPI_INTEGER,0,data%private%COMM,ERCODE)
! Terminate if error occurred
      IF (data%ERROR < 0) GO TO 400

      FLAG = 0
! Do not know how long the arrays A,IRN,JCN need to be ... choose
! eight times max. number of entries in submatrix
      LA = 8*NZMX
! this gives prob. in test deck
!! fixed!!      LA = NZMX+3
! !!! we need to monitor this for real examples ... is it
! !!! too much (look at size of info(2), the no. of garbage collections)
! Note: we could get different results when running on a different
! number of processes since NZMX will depend on number of processes


  250 DEALLOCATE (A,IRN,JCN,STAT=ST)
      ALLOCATE (A(1:LA),IRN(1:LA),JCN(1:LA),STAT=ST)
      IF (ST /= 0 .AND. FLAG == -3) THEN
! Try reducing space required (if FLAG = -3, we have had to increase
! LA so no point in now reducing it below min. for MP48_MA60A/AD)
         LA = data%INFO_MA60(3,JBLOCK)
         DEALLOCATE (A,IRN,JCN,STAT=ST)
         ALLOCATE (A(1:LA),IRN(1:LA),JCN(1:LA),STAT=ST)
         IF (ST /= 0) GO TO 300
      END IF
  251 IF (ST /= 0 .AND. FLAG == 0) THEN
         LA = LA/2
         DEALLOCATE (A,IRN,JCN,STAT=ST)
         ALLOCATE (A(1:LA),IRN(1:LA),JCN(1:LA),STAT=ST)
         IF (ST /= 0) THEN
            IF (LA > 2*NZMX) GO TO 251
            GO TO 300
         END IF
      END IF

! Call MP48_MA60A/AD to perform analyse.
! Loop over submatrices.
      FLAG = 0
      IF (data%IDUP > 0) MYIW1(1:NEQ) = 0
      MYLMX = 0
      IF (ICONTROL(7) == 1 .OR. &
        (ICONTROL(7) == 2 .AND. RANK == 0)) THEN
! Find length required for VALUES.
         L = 1
         DO LCNT = 1,ICOUNT_RANK
            L = MAX(L,NVARSB(LCNT))
         END DO
         DEALLOCATE (data%VALUES,STAT=ST)
         ALLOCATE (data%VALUES(1:L),STAT=ST)
         IF (ST /= 0) GO TO 300
      END IF

      T1 = 0.0
      DO 280 LCNT = 1,ICOUNT_RANK
         T1 = MPI_WTIME()
         JBLOCK = data%IBLOCK(LCNT)

! If this is the second time we are going through this analyse
! loop because we have had to increase LA, we only
! want to recall MP48_MA60A/AD if LA was actually too small for JBLOCK
! (LA may only have been too small for only one of the submatrices
! dealt with by the current processor).
! If MP48_MA60A/AD has been called successfully for submatrix, jump
! to end of loop)
         IF (LMP48_MA60(JBLOCK)) GO TO 280

         M = NEQSUB(LCNT)
         NC = data%private%NCOLA(JBLOCK)
! L1 points to beginning of col. indices for block JBLOCK
         L1 = LPTR(LCNT)
         L2 = L1 + NEQSUB(LCNT)
! J1 points to first entry in EQVAR for JBLOCK
! and J2 points to the last entry
         J1 = JPTR(LCNT)
         J2 = J1 + NVARSB(LCNT) - 1
         J1V = VPTR(LCNT)
         J2V = J1V + NVARSB(LCNT) - 1
! Set NZBLK to number of entries in block
         NZBLK = NVARSB(LCNT) - data%private%IDUP(LCNT)

         IF (ICONTROL(7) == 1 .OR. &
            (ICONTROL(7) == 2 .AND. RANK == 0)) THEN
! Read in VALUES from sequential file.
! We only hold values for one block at a time.
            IF (ICONTROL(7) == 1) OPEN (UNIT=STRM,IOSTAT=IOS,ERR=300, &
                  FILE=VALNAM(LCNT),FORM='UNFORMATTED')
            IF (ICONTROL(7) == 2) OPEN (UNIT=STRM,IOSTAT=IOS,ERR=300, &
                  FILE=data%VALNAM(JBLOCK),FORM='UNFORMATTED')
            READ (STRM) data%VALUES(1:NVARSB(LCNT))
            CLOSE (UNIT=STRM)

            J1V = 1
            J2V = J1V + NVARSB(LCNT) - 1
         END IF

! We have to squeeze out any duplicates before we call MC46
         IF (data%private%IDUP(LCNT) == 0) THEN
! No duplicates in this block so do copy of data
           IF (ICONTROL(7) == 5) THEN
              A(1:NZBLK) = data%RVAL(J1V:J2V)
           ELSE
              A(1:NZBLK) = data%VALUES(J1V:J2V)
           END IF
           IRN(1:NZBLK) = data%EQVAR(J1:J2)
           ROWPTR(1) = 1
           DO L = 1,M
             ROWPTR(L+1) = ROWPTR(L) + &
                           data%EQPTR(L1+L) - data%EQPTR(L1+L-1)
           END DO
         ELSE
! Duplicates to be got rid of.
           A(1:NZBLK) = ZERO
           ROWPTR(1) = 1
           J = 0
           KSTART = 1
! Loop over the rows
           IF (ICONTROL(7) == 5) THEN
             DO 270 L = 1,M
               NVAR = data%EQPTR(L1+L) - data%EQPTR(L1+L-1)
               KSTOP = KSTART + NVAR
               ROWPTR(L+1) = ROWPTR(L)
               DO 260 K = KSTART,KSTOP-1
                 I = data%EQVAR(J1+K-1)
                 J = J + 1
                 JJ = data%private%IMAP(J,LCNT)
                 IF (MYIW1(I) < L+MYLMX) THEN
                    MYIW1(I) = L + MYLMX
                    ROWPTR(L+1) = ROWPTR(L+1) + 1
! Set col. index in IRN
                    IRN(JJ) = I
                 END IF
                 A(JJ) = A(JJ) + data%RVAL(J1V+K-1)
  260          CONTINUE
               KSTART = KSTOP
  270        CONTINUE
             MYLMX = MYLMX + NEQ
           ELSE
             DO 271 L = 1,M
               NVAR = data%EQPTR(L1+L) - data%EQPTR(L1+L-1)
               KSTOP = KSTART + NVAR
               ROWPTR(L+1) = ROWPTR(L)
               DO 261 K = KSTART,KSTOP-1
                 I = data%EQVAR(J1+K-1)
                 J = J + 1
                 JJ = data%private%IMAP(J,LCNT)
                 IF (MYIW1(I) < L+MYLMX) THEN
                    MYIW1(I) = L + MYLMX
                    ROWPTR(L+1) = ROWPTR(L+1) + 1
! Set col. index in IRN
                    IRN(JJ) = I
                 END IF
                 A(JJ) = A(JJ) + data%VALUES(J1V+K-1)
  261          CONTINUE
               KSTART = KSTOP
  271        CONTINUE
             MYLMX = MYLMX + NEQ
           END IF
        END IF

        CALL MC46AD(M,NC,NZBLK,.TRUE.,IRN,A,ROWPTR,COLPTR,IW,IFLAG46)
! There should not be any errors from MC46
        IF (IFLAG46 < 0) THEN
          FLAG = -101
          GO TO 300
        END IF

! IRN now holds row indices, stored by cols.
! COLPTR holds col.allocatables.

! Set control parameters for MP48_MA60A/AD
         CALL MP48_MA60ID(data%private%CNTL_60(1:4,LCNT), &
                          data%private%ICNTL_60(1:10,LCNT))

! ICNTL_60(4,LCNT): the pivot search is limited to ICNTL_60(4,LCNT)
! cols (Zlatev strategy). Default is 3.
         data%private%ICNTL_60(4,LCNT) = ICONTROL(8)

         data%private%ICNTL_60(5,LCNT) = ICONTROL(5)
! Set ICNTL(6) equal to the number of columns in the border.
! These cols are then forced to end ... MP48_MA60A/AD terminates
! after considering the first N - ICNTL(6) columns
         data%private%ICNTL_60(6,LCNT) = data%NGUARD(JBLOCK)

! New parameter: ICNTL_60(8,LCNT) is the
! number of compresses allowed in the analysis phase before
! an error return is invoked.  Default value is 10.

! Switch off printing within MP48_MA60
         data%private%ICNTL_60(1,LCNT) = -1
         data%private%ICNTL_60(2,LCNT) = -1
         data%private%ICNTL_60(3,LCNT) = -1
!        data%private%ICNTL_60(3,LCNT) =  4
!        data%private%ICNTL_60(3,LCNT) =  3

! Switch from sparse to full matrix processing
         data%private%CNTL_60(1,LCNT) = CONTROL(1)
! Threshold pivoting parameter
         data%private%CNTL_60(2,LCNT) = CONTROL(2)
! Drop tolerances
         data%private%CNTL_60(3:4,LCNT) = CONTROL(3:4)

         data%private%IQ(1:NC,LCNT) = COLPTR(1:NC)

! Ready for MP48_MA60AD.
         CALL MP48_MA60AD(M,NC,NZBLK,LA,A,IRN,JCN,data%private%IQ(1:NC,LCNT), &
               data%private%CNTL_60(1:4,LCNT), &
               data%private%ICNTL_60(1:10,LCNT),data%private%IP(1:M,LCNT), &
               NP,JFIRST,LENR,LASTR,NEXTR,IW,IFIRST,LENC, &
               LASTC,NEXTC,data%INFO_MA60(1:20,JBLOCK),RINFO)

         data%OPSA(JBLOCK) = RINFO

! Monitor space used with this write statement.
! if we have done several compressions or error is -7,
! better to set LA larger
!         write (6,'(a,i3/i3,2i8,i4)')
!     &   'on exit ma60a/ad for submatrix ',
!     &   jblock,data%INFO_MA60(1,JBLOCK),la,data%INFO_MA60(3,JBLOCK),
!     &   data%INFO_MA60(2,JBLOCK)
!         write (6,'(a,2i4)') 'Submatrix  info(7) ',
!     &   jblock,data%INFO_MA60(7,JBLOCK)

! Check for errors from MP48_MA60A/AD.
! Only possible errors are that there was insufficient space
! or too many compressions of the data. In both cases,
! LA should be increased to at least data%INFO_MA60(3,JBLOCK)
         IF (data%INFO_MA60(1,JBLOCK) < 0) THEN
            FLAG = data%INFO_MA60(1,JBLOCK)
            IF (FLAG == -3 .OR. FLAG == -7) THEN
               FLAG = -3
               LA = MAX(ICONTROL(6)*LA,data%INFO_MA60(3,JBLOCK))
            ELSE
               FLAG = -102
            END IF
            GO TO 290
         ELSE
            LMP48_MA60(JBLOCK) = .TRUE.
         END IF

         data%private%NP(JBLOCK) = NP
! Set NF to be NC-NP
! NP is no cols processed in packed (ie sparse) storage so that
! NF is no. cols processed in full storage
         data%private%NF(JBLOCK) = NC - NP
         T2 = MPI_WTIME()
         data%TIMEA(JBLOCK) = data%TIMEA(JBLOCK) + T2 - T1

!         write (6,*) 'no. of garbage collections for submatrix ',
!     &   jblock,data%INFO_MA60(2,JBLOCK)

  280 CONTINUE
  290 IF (FLAG == -3) GO TO 250

! Synchronize
  300 CALL MPI_BARRIER(data%private%COMM,ERCODE)
      data%STAT = ST
      data%IOSTAT = IOS
      CALL MP48ND(data,NPROC,IFLAG,FLAG)
 ! Broadcast error status from root
      CALL MPI_BCAST(data%ERROR,1,MPI_INTEGER,0,data%private%COMM,ERCODE)
! Terminate if error occurred
      IF (data%ERROR < 0) GO TO 400

! Collect information arrays on root.
      IF (RANK == 0) THEN
         data%DROP = ZERO
         DO IPROC = 1,NPROC-1
            SRC = IPROC
            ISTRT = data%IPLIST(IPROC)
            ISTOP = data%IPLIST(IPROC+1) - 1
            DO I = ISTRT,ISTOP
               JBLOCK = data%LIST(I)
               CALL MPI_RECV(data%INFO_MA60(1,JBLOCK),20,MPI_INTEGER, &
                             SRC,1,data%private%COMM,STAT,ERCODE)
               CALL MPI_RECV(data%OPSA(JBLOCK),1, &
                             MPI_DOUBLE_PRECISION, &
                             SRC,2,data%private%COMM,STAT,ERCODE)
               CALL MPI_RECV(data%private%NF(JBLOCK),1,MPI_INTEGER, &
                             SRC,3,data%private%COMM,STAT,ERCODE)
               CALL MPI_RECV(data%private%NP(JBLOCK),1,MPI_INTEGER, &
                             SRC,4,data%private%COMM,STAT,ERCODE)
               CALL MPI_RECV(data%TIMEA(JBLOCK),1,MPI_DOUBLE_PRECISION, &
                             SRC,5,data%private%COMM,STAT,ERCODE)
               data%DROP   = data%DROP + data%INFO_MA60(6,JBLOCK)
            END DO
         END DO
         DO LCNT = 1,ICOUNT_RANK
            JBLOCK = data%IBLOCK(LCNT)
! Accumulate statistics
            data%DROP   = data%DROP + data%INFO_MA60(6,JBLOCK)
         END DO

      ELSE
         DO LCNT = 1,ICOUNT_RANK
            JBLOCK = data%IBLOCK(LCNT)
            CALL MPI_SEND(data%INFO_MA60(1,JBLOCK),20,MPI_INTEGER, &
                          0,1,data%private%COMM,ERCODE)
            CALL MPI_SEND(data%OPSA(JBLOCK),1, &
                          MPI_DOUBLE_PRECISION, &
                          0,2,data%private%COMM,ERCODE)
            CALL MPI_SEND(data%private%NF(JBLOCK),1,MPI_INTEGER, &
                          0,3,data%private%COMM,ERCODE)
            CALL MPI_SEND(data%private%NP(JBLOCK),1,MPI_INTEGER, &
                          0,4,data%private%COMM,ERCODE)
            CALL MPI_SEND(data%TIMEA(JBLOCK),1,MPI_DOUBLE_PRECISION, &
                          0,5,data%private%COMM,ERCODE)
         END DO
      END IF

  400 CALL MPI_BARRIER(data%private%COMM,ERCODE)
      IF (RANK == 0 .AND. data%ERROR.GE.0) THEN

         IF (LDGN.GE.2) WRITE (DGN,FMT='(/A/6(ES12.5))') &
       ' The predicted number of flops needed for each submatrix is:', &
         data%OPSA(1:NBLOCK)

           IF (LDGN.GE.3) WRITE (DGN,FMT='(A/8(F9.3))') &
      ' The analyse wall clock time (secs) for each submatrix is:', &
           (data%TIMEA(JBLOCK),JBLOCK=1,NBLOCK)

      END IF

      CALL MPI_BCAST(data%ERROR,1,MPI_INTEGER,0,data%private%COMM,ERCODE)
      data%private%ERROR = data%ERROR

      IF (data%ERROR.GE.0) THEN
        DO LCNT = 1,ICOUNT_RANK
          data%private%NEQSUB(LCNT) = NEQSUB(LCNT)
          data%private%JPTR(LCNT) = JPTR(LCNT)
          data%private%LPTR(LCNT) = LPTR(LCNT)
          data%private%VPTR(LCNT) = VPTR(LCNT)
        END DO
      END IF

      DEALLOCATE (A,STAT=ST)
      DEALLOCATE (VALUES,STAT=ST)

      DEALLOCATE (COLPTR,STAT=ST)
      DEALLOCATE (IFLAG,STAT=ST)
      DEALLOCATE (IFIRST,STAT=ST)
      DEALLOCATE (IW,STAT=ST)
      DEALLOCATE (IW1,STAT=ST)
      DEALLOCATE (MYIW1,STAT=ST)
      DEALLOCATE (IRN,STAT=ST)
      DEALLOCATE (JCN,STAT=ST)
      DEALLOCATE (JFIRST,STAT=ST)
      DEALLOCATE (JPTR,STAT=ST)
      DEALLOCATE (LPTR,STAT=ST)
      DEALLOCATE (LENC,STAT=ST)
      DEALLOCATE (LENR,STAT=ST)
      DEALLOCATE (LASTC,STAT=ST)
      DEALLOCATE (LASTR,STAT=ST)
      DEALLOCATE (NEXTC,STAT=ST)
      DEALLOCATE (NEXTR,STAT=ST)
      DEALLOCATE (NVARSB,STAT=ST)
      DEALLOCATE (NEQSUB,STAT=ST)
      DEALLOCATE (ROWPTR,STAT=ST)
      DEALLOCATE (VPTR,STAT=ST)
      DEALLOCATE (SBUF,STAT=ST)
      DEALLOCATE (LMP48_MA60,STAT=ST)

 9000 FORMAT (/' Error return with data%JOB = ',I2,'. ', &
                'data%ERROR = ',I4)
 9010 FORMAT (/' Warning return with data%JOB = ',I2,'. ', &
                'data%ERROR = ',I4)
 9020 FORMAT (' ',A,' has not been allocated')
 9060 FORMAT (' Number of duplicates found  = ',I6)
 9100 FORMAT (' ',A,' is of size ',I10,' Increase to at least ',I10)
 9110 FORMAT (' Size of ',A,' is too small on one or more process')
      END SUBROUTINE MP48BD

!*************************************************

      SUBROUTINE MP48CD(data)

! This subroutine performs factorize phase using MP48_MA60B/BD
! JOB = 4 call.

      USE HSL_MP48_DATA_DOUBLE
      IMPLICIT NONE

      TYPE (MP48_DATA) data

!     .. Parameters ..
      REAL(WP), PARAMETER :: ZERO = 0.0_WP
!     REAL(DP), PARAMETER :: ZEROD = 0.0_DP

!     .. Local Scalars ..

      REAL(WP) RINFO
!     REAL(DP) TIM60
      REAL(DP) T1,T2
      INTEGER BLANK,DEST,DGN,ERCODE,ERR, &
              FLAG, &
              I,I1,I2,I1F,I2F,I3,I4,II1,II2,INFO4,INFO10,II,IDUM, &
              IPROC,ISTRT,ISTOP,ICOUNT_RANK,INTERF,IFLAG46,IF1I,IOUT, &
              J,J1,J2,J1V,J2V,JBLOCK,JOB50,JOB48,JJ,JPOS, &
              K,KCOPY,KL,KL1,KSTART,KSTOP, &
              L,L1,L2,LA,LFACT,LCNT,LDGN,LMX, &
              LIRNF,LICOPY,LFCOPY,LKEEP48,LF,LENGTH, &
              M,MXINT,MF,MNF,MNFMX, &
              NBLOCK,NC,NE,NEQ,NF,NPROC,NU,N10,N100,N1000, &
              NZBLK,NCMX,NZMX,NINTER,RINTER,NVAR,NNA,NP, &
              RANK,STRM,SRC,VALCRD,WRN,STD,TOTNE_SCHUR

      INTEGER IOS
!  Holds IOSTAT parameter for files.
      INTEGER ST
!  Holds STAT parameter (ALLOCATE statement)
      LOGICAL LERR,LWRN
      LOGICAL EX
!  Used by INQUIRE. Set to .TRUE. if  unit exists.
      LOGICAL OPEN
!  Used by INQUIRE. Set to .TRUE. if file open.
      LOGICAL LTEMP

!     .. Local Arrays ..

      CHARACTER (LEN = 10) :: NUM
      REAL(WP) :: CONTROL(30)

      REAL(WP), DIMENSION (:), ALLOCATABLE :: A
      REAL(WP), DIMENSION (:), ALLOCATABLE :: W
      REAL(WP), DIMENSION (:), ALLOCATABLE :: VALUES
      REAL(WP), DIMENSION (:), ALLOCATABLE :: FACT_COPY
      REAL(WP), DIMENSION (:), ALLOCATABLE :: FF_COPY

!     CHARACTER,DIMENSION(:), allocatable :: VALNAM*128
      CHARACTER,DIMENSION(:), ALLOCATABLE :: VALNAM*128
      CHARACTER (LEN = 128) :: TEMP_VALNAM

      INTEGER ICONTROL(20),STAT(MPI_STATUS_SIZE,20),JNFO(5)

      INTEGER, DIMENSION (:), ALLOCATABLE :: COLPTR
      INTEGER, DIMENSION (:), ALLOCATABLE :: IRN
      INTEGER, DIMENSION (:), ALLOCATABLE :: IBLOCK
      INTEGER, DIMENSION (:), ALLOCATABLE :: IFLAG
      INTEGER, DIMENSION (:), ALLOCATABLE :: IW
      INTEGER, DIMENSION (:), ALLOCATABLE :: IW1
      INTEGER, DIMENSION (:), ALLOCATABLE :: KPCOPY
      INTEGER, DIMENSION (:), ALLOCATABLE :: NEQSUB
      INTEGER, DIMENSION (:), ALLOCATABLE :: NVARSB
      INTEGER, DIMENSION (:), ALLOCATABLE :: ROWPTR
      INTEGER, DIMENSION (:), ALLOCATABLE :: VPTR
      INTEGER, DIMENSION (:), ALLOCATABLE :: NE_SCHUR
      INTEGER, DIMENSION (:,:), ALLOCATABLE :: IQM

!     .. External Subroutines ..
      EXTERNAL MA48ID,MA48AD,MA48BD,MC46AD

!     .. MPI routines
      EXTERNAL MPI_BARRIER,MPI_BCAST,MPI_COMM_RANK, &
               MPI_GATHER,MPI_SEND,MPI_RECV

!     .. Intrinsic Functions ..
      INTRINSIC MAX

      FLAG = 0
      RANK = data%RANK
      NPROC = data%NPROC

! Set diagnostics controls
      DGN = data%ICNTL(3)
      LDGN = data%ICNTL(4)
      IF (DGN < 0) LDGN = 0
! Set error controls
      LERR = .FALSE.
      IF (data%ICNTL(1).GE.0 .AND. data%ICNTL(4).GE.1) LERR = .TRUE.
      ERR = data%ICNTL(1)
! Set warning controls
      LWRN = .FALSE.
      IF (data%ICNTL(2).GE.0 .AND. data%ICNTL(4).GE.1) LWRN = .TRUE.
      WRN = data%ICNTL(2)

      ICONTROL(1:20) =  data%ICNTL(1:20)
      CONTROL(1:4)   =  data%CNTL(1:4)

! Deallocate arrays that might still be allocated
      DEALLOCATE (data%private%A,STAT=ST)

      NBLOCK = data%private%NBLOCK
      NEQ    = data%private%NEQ
      IF (RANK == 0) THEN

! Check data for errors. Check required arrays have been allocated.
        IF (data%FACT_JOB /= 1) THEN
! Arrays can have changed since the call with JOB = 3 so perform
! tests again.
         IF (ICONTROL(7) == 1 .OR. ICONTROL(7) == 2) THEN
! Check size of array VALNAM
! First check it has been allocated
            IF (.NOT. allocated(data%VALNAM)) THEN
               FLAG = -10
            ELSE
               IF (SIZE(data%VALNAM) < data%NBLOCK) FLAG = -10
            END IF
            IF (FLAG == -10) THEN
               data%ERROR = -10
               IF (LERR) THEN
                  WRITE (ERR,FMT=9000) data%JOB,data%ERROR
                  IF (.NOT. allocated(data%VALNAM)) THEN
                     WRITE (ERR,FMT=9020) ' data%VALNAM'
                  ELSE IF (SIZE(data%VALNAM) < data%NBLOCK) THEN
                     WRITE (ERR,FMT=9100) ' data%VALNAM', &
                     SIZE(data%VALNAM),data%NBLOCK
                  END IF
               END IF
               GO TO 20
            END IF
         ELSE IF (ICONTROL(7) == 3) THEN
! Check size of array VALUES
            IF (.NOT. allocated(data%VALUES)) THEN
               FLAG = -10
            ELSE
               VALCRD = SIZE(data%VALUES)
               IF (VALCRD < data%NE) FLAG = -10
            END IF
            IF (FLAG == -10) THEN
               data%ERROR = -10
               IF (LERR) THEN
                  WRITE (ERR,FMT=9000) data%JOB,data%ERROR
                  IF (.NOT. allocated(data%VALUES)) THEN
                     WRITE (ERR,FMT=9020) ' data%VALUES'
                  ELSE IF (VALCRD < data%NE) THEN
                     WRITE (ERR,FMT=9100) &
                  ' data%VALUES',VALCRD,data%NE
                  END IF
               END IF
               GO TO 20
            END IF
         END IF
        END IF
        IF (ICONTROL(11) < 0) THEN
! Check size of array FILES
! First check it has been allocated
            IF (.NOT. allocated(data%FILES)) THEN
               FLAG = -25
            ELSE IF (SIZE(data%FILES,1) < 2) THEN
               FLAG = -25
            ELSE
               IF (SIZE(data%FILES,2) < data%NBLOCK+1) FLAG = -25
            END IF
            IF (FLAG == -25) THEN
               data%ERROR = -25
               IF (LERR) THEN
                  WRITE (ERR,FMT=9000) data%JOB,data%ERROR
                  IF (.NOT. allocated(data%FILES)) THEN
                     WRITE (ERR,FMT=9020) ' data%FILES'
                  ELSE IF (SIZE(data%FILES,1) < 2) THEN
                     WRITE (ERR,FMT=9120) 'data%FILES', &
                     SIZE(data%FILES,1),2
                  ELSE IF (SIZE(data%FILES,2) < data%NBLOCK+1) THEN
                     WRITE (ERR,FMT=9130) 'data%FILES', &
                     SIZE(data%FILES,2),data%NBLOCK+1
                  END IF
               END IF
               GO TO 20
            END IF
        END IF

      END IF
! Broadcast error status from root
   20 CALL MPI_BCAST(FLAG,1,MPI_INTEGER,0,data%private%COMM,ERCODE)
! Terminate if error occurred
      IF (FLAG < 0) GO TO 890

      DEALLOCATE (IFLAG,STAT=ST)
      ALLOCATE (IFLAG(NPROC),STAT=ST)
! Initialise the error flag (so that it is always well-defined)
      IF (ST == 0) IFLAG(1:NPROC) = 0

      FLAG = 0
      IF (data%FACT_JOB /= 1) THEN
! If ICNTL(7) = 4 or 5, check VALUES/RVAL
      IF (ICONTROL(7) == 4) THEN
! On each process, check array VALUES
         IF (.NOT. allocated(data%VALUES)) THEN
            FLAG = -10
         ELSE
            VALCRD = SIZE(data%VALUES)
            IF (VALCRD < data%NE) FLAG = -11
         END IF
      ELSE IF (ICONTROL(7) == 5) THEN
         IF (.NOT. allocated(data%RVAL)) THEN
            FLAG = -10
         ELSE
            VALCRD = SIZE(data%RVAL)
! Check length of RVAL. Loop over submatrices assigned
! to process
            L = 0
            DO LCNT = 1,data%ICOUNT(RANK)
               JBLOCK = data%IBLOCK(LCNT)
               L = L + data%ENTRIES(JBLOCK)
            END DO
            IF (VALCRD < L) FLAG = -11
         END IF
      END IF
      IF (ICONTROL(7) == 4 .OR. ICONTROL(7) == 5) THEN
         CALL MPI_GATHER(FLAG,1,MPI_INTEGER,IFLAG,1,MPI_INTEGER,0, &
                      data%private%COMM,ERCODE)
! Check for error on root
         IF (RANK == 0) THEN
            DO 21 IPROC = 1,NPROC
               IF (IFLAG(IPROC) /= 0) THEN
                  data%ERROR = -10
                  IF (LERR) THEN
                     WRITE (ERR,FMT=9000) data%JOB,data%ERROR
                     IF (ICONTROL(7) == 4) THEN
                       IF (IFLAG(IPROC) == -10) THEN
                         WRITE (ERR,FMT=9020) 'data%VALUES'
                       ELSE
                         WRITE (ERR,FMT=9110) 'data%VALUES'
                       END IF
                     ELSE
                       IF (IFLAG(IPROC) == -10) THEN
                         WRITE (ERR,FMT=9020) 'data%RVAL'
                       ELSE
                         WRITE (ERR,FMT=9110) 'data%RVAL'
                       END IF
                     END IF
                  END IF
                  GO TO 22
               END IF
   21       CONTINUE
         END IF
      END IF
! Broadcast error status from root
   22 CALL MPI_BCAST(data%ERROR,1,MPI_INTEGER,0,data%private%COMM,ERCODE)
! Terminate if error occurred
      IF (data%ERROR < 0) GO TO 890
      END IF

      ST = 0
      IOS = 0
      ICOUNT_RANK = data%ICOUNT(RANK)

      DEALLOCATE (IBLOCK,NEQSUB,NVARSB,VALNAM,VPTR,STAT=ST)
      ALLOCATE (IBLOCK(1:ICOUNT_RANK), &
                NEQSUB(1:ICOUNT_RANK), &
                NVARSB(1:ICOUNT_RANK), &
                VALNAM(1:ICOUNT_RANK), &
                VPTR(1:ICOUNT_RANK),STAT=ST)
      IF (ST /= 0) GO TO 30

      IF (data%private%FJOB == 1) THEN
         DEALLOCATE (data%private%LFACT,data%private%LIRNF,STAT=ST)
         ALLOCATE (data%private%LFACT(1:ICOUNT_RANK), &
                   data%private%LIRNF(1:ICOUNT_RANK),STAT=ST)
         IF (ST /= 0) GO TO 30
      END IF

! Set data from analyse phase
      DO LCNT = 1,ICOUNT_RANK
         IBLOCK(LCNT) = data%IBLOCK(LCNT)
         NEQSUB(LCNT) = data%private%NEQSUB(LCNT)
         JBLOCK = IBLOCK(LCNT)
         NVARSB(LCNT) = data%ENTRIES(JBLOCK)
         VPTR(LCNT) = data%private%VPTR(LCNT)
      END DO

      INTERF = NBLOCK + 1
      IF (ICONTROL(11).GE.0) THEN
         DEALLOCATE (data%FILES,STAT=ST)
         ALLOCATE (data%FILES(1:2,1:INTERF),STAT=ST)
         IF (ST /= 0) GO TO 30
      ELSE IF (ICONTROL(11) < 0 .AND. RANK /= 0) THEN
! User has chosen file names on root ... we must broadcast them to
! other processes
         DEALLOCATE (data%FILES,STAT=ST)
         ALLOCATE (data%FILES(1:2,1:INTERF),STAT=ST)
         IF (ST /= 0) GO TO 30
      END IF

! Find suitable unit number if we are to read data from file
! Also needed if ICONTROL(11) nonzero
! Ensure STRM is defined
      FLAG = 0
      STRM = 1
      IF (ICONTROL(11) == 0) THEN
         IF (ICONTROL(7).GE.3) GO TO 30
         IF (ICONTROL(7) == 2 .AND. RANK /= 0) GO TO 30
      END IF

      DO STRM = 8,99
         IF (STRM == DGN .OR. STRM == ERR .OR. STRM == WRN) &
             CYCLE
         INQUIRE (UNIT=STRM,IOSTAT=IOS,ERR=30,EXIST=EX, &
                  OPENED=OPEN)
         IF (EX .AND. .NOT.OPEN) GO TO 25
      END DO
! No unit found. Jump to error return.
      FLAG = -20
      GO TO 30
   25 CONTINUE

   30 CALL MPI_BARRIER(data%private%COMM,ERCODE)
      IF (IOS /= 0) FLAG = -14
      data%STAT = ST
      CALL MP48ND(data,NPROC,IFLAG,FLAG)
      data%IOSTAT = IOS
! Broadcast error status from root
      CALL MPI_BCAST(data%ERROR,1,MPI_INTEGER,0,data%private%COMM,ERCODE)
! Terminate if error occurred
      IF (data%ERROR < 0) GO TO 890

      IF (ICONTROL(11) > 0) THEN
! Choose names for files to hold factors
         DO 60 JBLOCK = 1,NBLOCK
             NUM(1:10) = '0123456789'
             N1000 = JBLOCK/1000
             N100 = (JBLOCK - N1000*1000)/100
             N10 = (JBLOCK - N1000*1000 - N100*100)/10
             NU = JBLOCK - N1000*1000 - N100*100 - N10*10
             N1000 = N1000 + 1
             N100 = N100 + 1
             N10 = N10 + 1
             NU = NU + 1
             data%FILES(1,JBLOCK) = &
             'fact.'//NUM(N1000:N1000)//NUM(N100:N100)// &
              NUM(N10:N10)//NUM(NU:NU)
             data%FILES(2,JBLOCK) = &
             'integ.'//NUM(N1000:N1000)//NUM(N100:N100)// &
              NUM(N10:N10)//NUM(NU:NU)
   60    CONTINUE
         data%FILES(1,INTERF) = 'fact_interf'
         data%FILES(2,INTERF) = 'integ_interf'
      ELSE IF (ICONTROL(11) < 0) THEN
! Broadcast user-chosen file names from the root (don't
! need to send interface names)
         CALL MPI_BCAST(data%FILES,NBLOCK*2*128,MPI_CHARACTER,0, &
                        data%private%COMM,ERCODE)

      END IF

      FLAG = 0
! Send name of data files root to other processes.
! If data%FACT_JOB = 1 then this is first call and values
! have not been changed since analyse, so if CONTROL(7) = 3 the
! array VALUES has correct data.
! If data%FACT_JOB .ne. 1, we must send the required data from the host.
      IF (RANK == 0) THEN

         IF (data%FACT_JOB /= 1 .AND. ICONTROL(7) == 2) THEN
! Check how large array VALUES must be.
            L = 1
            DO JBLOCK = 1,NBLOCK
               L = MAX(L,data%ENTRIES(JBLOCK))
            END DO
            DEALLOCATE (VALUES,STAT=ST)
            ALLOCATE (VALUES(1:L),STAT=ST)
            IF (ST /= 0) GO TO 300
         END IF

         DO 80 IPROC = 1,NPROC-1
            DEST = IPROC
            ISTRT = data%IPLIST(IPROC)
            ISTOP = data%IPLIST(IPROC+1) - 1
            LCNT = 1
            DO 70 I = ISTRT,ISTOP
               JBLOCK = data%LIST(I)

               IF (data%FACT_JOB /= 1 .AND. ICONTROL(7) == 3) THEN
                  CALL MPI_SEND(data%VALUES(data%private%IPTR(JBLOCK)), &
                                data%ENTRIES(JBLOCK), &
                                MPI_DOUBLE_PRECISION, &
                                DEST,6,data%private%COMM,ERCODE)

               ELSE IF (data%FACT_JOB /= 1 .AND. ICONTROL(7) == 2) THEN
! Read values from host and send out to DEST
                  OPEN (UNIT=STRM,IOSTAT=IOS,ERR=300, &
                        FILE=data%VALNAM(JBLOCK),FORM='UNFORMATTED')
                  READ (STRM) VALUES(1:data%ENTRIES(JBLOCK))
                  CLOSE (UNIT=STRM,STATUS='KEEP')

                  CALL MPI_SEND(VALUES,data%ENTRIES(JBLOCK), &
                                MPI_DOUBLE_PRECISION, &
                                DEST,6,data%private%COMM,ERCODE)


               ELSE IF (ICONTROL(7) == 1) THEN
! Each process will read in values ... send out names of files
! where matrix entries are.
                  TEMP_VALNAM = ADJUSTL(data%VALNAM(JBLOCK))
                  BLANK = LEN_TRIM(TEMP_VALNAM)
                  CALL MPI_SEND(BLANK,1,MPI_INTEGER,DEST,7, &
                                data%private%COMM,ERCODE)
                  CALL MPI_SEND(TEMP_VALNAM,BLANK, &
                                MPI_CHARACTER,DEST,8,data%private%COMM,ERCODE)

               END IF
               LCNT = LCNT + 1

   70       CONTINUE
   80    CONTINUE
         IF (data%FACT_JOB /= 1 .AND. ICONTROL(7) == 2) &
            DEALLOCATE (VALUES,STAT=STD)

         IF (ICONTROL(7) == 1) THEN
            DO LCNT = 1,ICOUNT_RANK
               JBLOCK = data%IBLOCK(LCNT)
               TEMP_VALNAM = ADJUSTL(data%VALNAM(JBLOCK))
               BLANK = LEN_TRIM(TEMP_VALNAM)
               VALNAM(LCNT) = TEMP_VALNAM(1:BLANK)
            END DO
        END IF

      ELSE

! Add up how much space needed
         L = 0
         DO LCNT = 1,ICOUNT_RANK
            JBLOCK = data%IBLOCK(LCNT)
            L = L + NVARSB(LCNT)
         END DO

         IF (data%FACT_JOB /= 1) THEN
            IF (ICONTROL(7) == 2 .OR. ICONTROL(7) == 3) THEN
               DEALLOCATE (data%VALUES,STAT=ST)
               ALLOCATE (data%VALUES(1:L),STAT=ST)
               IF (ST /= 0) GO TO 300
            END IF
         END IF

! Receive data from the root
         DO 100 LCNT = 1,ICOUNT_RANK
            IF (data%FACT_JOB /= 1) THEN
              IF (ICONTROL(7) == 2 .OR. ICONTROL(7) == 3) &
                CALL MPI_RECV(data%VALUES(VPTR(LCNT)), &
                              NVARSB(LCNT),MPI_DOUBLE_PRECISION, &
                              0,6,data%private%COMM,STAT,ERCODE)
            END IF
            IF (ICONTROL(7) == 1) THEN
               CALL MPI_RECV(BLANK,1,MPI_INTEGER,0,7, &
                             data%private%COMM,STAT,ERCODE)
               CALL MPI_RECV(VALNAM(LCNT),BLANK, &
                             MPI_CHARACTER,0,8,data%private%COMM,STAT,ERCODE)
               VALNAM(LCNT) = VALNAM(LCNT)(1:BLANK)

            END IF
  100    CONTINUE

      END IF

      NZMX  = 0
      NCMX  = 0
      DO LCNT = 1,ICOUNT_RANK
         JBLOCK = IBLOCK(LCNT)
         NC = data%private%NCOLA(JBLOCK)
         NZMX = MAX(NZMX,NVARSB(LCNT))
         NCMX = MAX(NCMX,NC)
      END DO

      IF (data%private%FJOB == 1) THEN
        DO LCNT = 1,ICOUNT_RANK
           JBLOCK = IBLOCK(LCNT)
! INFO_MA60(4,JBLOCK) is length real space required for
! factorization on JBLOCK.
! INFO_MA60(10,JBLOCK) is length of integer space.
           data%private%LFACT(LCNT) = data%INFO_MA60(4,JBLOCK)
           data%private%LIRNF(LCNT) = data%INFO_MA60(10,JBLOCK)
        END DO

! We only allocate IPTRL and IPTRU if this is the first call
! with JOB = 4 (ie if data%private%FJOB = 1)
        DEALLOCATE (data%private%IPTRL,data%private%IPTRU, &
                    STAT=ST)
        ALLOCATE (data%private%IPTRL(1:NCMX,1:ICOUNT_RANK), &
                  data%private%IPTRU(1:NCMX,1:ICOUNT_RANK), &
                  STAT=ST)
        IF (ST /= 0) GO TO 300
      END IF

      DEALLOCATE (A,STAT=ST)
      DEALLOCATE (COLPTR,STAT=ST)
      DEALLOCATE (IQM,STAT=ST)
      DEALLOCATE (IRN,STAT=ST)
      DEALLOCATE (ROWPTR,STAT=ST)
      ALLOCATE (A(1:NZMX), &
                COLPTR(1:NCMX+1), &
                IQM(1:NCMX+1,1:ICOUNT_RANK), &
                IRN(1:NZMX), &
                ROWPTR(1:data%private%LORDER+1), &
                STAT=ST)
      IF (ST /= 0) GO TO 300

      DEALLOCATE (FACT_COPY,NE_SCHUR,STAT=ST)
!      IF (RANK == 0) THEN
!         MNFMX  = 0
!         DO JBLOCK = 1,NBLOCK
!            MF = data%INFO_MA60(7,JBLOCK)
!            NF = data%private%NF(JBLOCK)
!            MNF = MF*NF
!            MNFMX = MAX(MNFMX,MNF)
!         END DO
!         ALLOCATE
!     &    (FACT_COPY(1:MNFMX,1:NBLOCK),NE_SCHUR(1:NBLOCK),STAT=ST)
!      ELSE
         LICOPY  = 0
         LFCOPY  = 0
         DO LCNT = 1,data%ICOUNT(RANK)
            JBLOCK = IBLOCK(LCNT)
            MF = data%INFO_MA60(7,JBLOCK)
            NF = data%private%NF(JBLOCK)
            LICOPY = LICOPY + MF
            LFCOPY = LFCOPY + MF*NF
         END DO
         ALLOCATE (FACT_COPY(1:LFCOPY),NE_SCHUR(1:NBLOCK),STAT=ST)
!      END IF
      IF (ST /= 0) GO TO 300
      FLAG = 0
      LTEMP = .TRUE.

! Set job parameter for MP48_MA60BD.

! If data%FACT_JOB = 1 and data%private%FJOB = 1 then JOB50 = 1
! (first call since analyse)

! If data%FACT_JOB = 2  and data%private%FJOB = 1 then JOB50 = 1
! (first call since analyse ... the difference
! between this and the above case, is that with data%FACT_JOB = 1
! the numerical values must not have changed between the call to
! analyse phase and the factorise phase. On this call, the numerical
! values may have changed, so we have had to resend data
! from host).

! If data%FACT_JOB = 2 and data%private%FJOB > 1 then JOB50 = 1
! (more than one call since analyse ... values changed so have to send
! data from host). The difference between this call and one
! with data%private%FJOB = 1 is that we save later by not having
! to redo the calls to MA48A/AD (which can be expensive).

! If data%FACT_JOB = 3 then JOB50 = 2 (fast factorize).

! Note:
! if data%FACT_JOB = 3 then we must also have data%private%FJOB > 1
! (as we must have made a call with data%FACT_JOB = 1 or 2)
! and hence data%private%FJOB = 1 indicates a first call

      JOB50 = 1
      IF (data%FACT_JOB == 3) JOB50 = 2

      FLAG = 0
  250 CONTINUE
      IF (JOB50 == 1 .AND. ICONTROL(11) == 0) THEN
! First run. Not using files to hold factor data
         LFACT = 0
         LIRNF = 0
         DO LCNT = 1,ICOUNT_RANK
            JBLOCK = IBLOCK(LCNT)
            LFACT = LFACT + data%private%LFACT(LCNT)
            LIRNF = LIRNF + data%private%LIRNF(LCNT)
         END DO
         DEALLOCATE (data%FACT,data%IRNF,STAT=ST)
         ALLOCATE (data%FACT(1:LFACT),data%IRNF(1:LIRNF),STAT=ST)
         IF (ST /= 0) GO TO 300
      END IF

! Call MP48_MA60B/BD to perform factorize.
! Loop over submatrices.
      IF (data%IDUP > 0) THEN
         DEALLOCATE (IW1,STAT=ST)
         ALLOCATE (IW1(1:NEQ),STAT=ST)
         IW1(1:NEQ) = 0
      END IF
      LMX = 0
      IF (ST /= 0) GO TO 300
      FLAG = 0
      I2 = 0
      I2F = 0
      I4 = 1
      IF (ICONTROL(7) == 1 .OR. &
        (ICONTROL(7) == 2 .AND. RANK == 0)) THEN
! Find length required for VALUES.
         L = 1
         DO LCNT = 1,ICOUNT_RANK
            L = MAX(L,NVARSB(LCNT))
         END DO
         DEALLOCATE (data%VALUES,STAT=ST)
         ALLOCATE (data%VALUES(1:L),STAT=ST)
         IF (ST /= 0) GO TO 300
      END IF

      IF (ICONTROL(11) /= 0) THEN
         L = 1
         L1 = 1
         DO LCNT = 1,ICOUNT_RANK
            JBLOCK = IBLOCK(LCNT)
            IF (JOB50 == 1) THEN
               LFACT = data%private%LFACT(LCNT)
               LIRNF = data%private%LIRNF(LCNT)
               L = MAX(L,LFACT)
               L1 = MAX(L1,LIRNF)
            ELSE
              INFO4 = data%INFO_MA60(4,JBLOCK)
              INFO10 = data%INFO_MA60(10,JBLOCK)
              L = MAX(L,INFO4)
              L1 = MAX(L1,INFO10)
            END IF
         END DO
         DEALLOCATE (data%FACT,data%IRNF,STAT=ST)
         ALLOCATE (data%FACT(1:L),data%IRNF(1:L1),STAT=ST)
         IF (ST /= 0) GO TO 300
      END IF

      T1 = 0.0
      DO 280 LCNT = 1,ICOUNT_RANK
         JBLOCK = IBLOCK(LCNT)
         T1 = MPI_WTIME()
         LFACT = data%private%LFACT(LCNT)
         LIRNF = data%private%LIRNF(LCNT)
         M = NEQSUB(LCNT)
         NC = data%private%NCOLA(JBLOCK)
         IF (ICONTROL(11) /= 0) THEN
! Files to be used for factor data
            IF (JOB50 /= 1) THEN
! JOB50 = 2
! We are doing a fast factorize, so FACT and IRNF will have to be
! read back in using the stored lengths
              INFO4 = data%INFO_MA60(4,JBLOCK)
              INFO10 = data%INFO_MA60(10,JBLOCK)
              LFACT = INFO4
              LIRNF = INFO10
              CLOSE (UNIT=STRM)
              OPEN (UNIT=STRM,IOSTAT=IOS,ERR=300, &
                    FILE=data%FILES(1,JBLOCK),FORM='UNFORMATTED')
              READ (STRM) data%FACT(1:LFACT)
              CLOSE (UNIT=STRM)
              OPEN (UNIT=STRM,IOSTAT=IOS,ERR=300, &
                    FILE=data%FILES(2,JBLOCK),FORM='UNFORMATTED')
              READ (STRM) data%IRNF(1:LIRNF)
              CLOSE (UNIT=STRM)
            END IF
            I1 = 1
            I1F = 1
            I2 = LIRNF
            I2F = LFACT

         ELSE
! Files not used ... data for each submatrix written
! consecutively to FACT and IRNF
            I1 = I2 + 1
            I2 = I1 + LIRNF - 1
            I1F = I2F + 1
            I2F = I1F + LFACT - 1
         END IF

! L1 points to beginning of col. indices for block JBLOCK
         L1 = data%private%LPTR(LCNT)
! J1 points to first entry in EQVAR for JBLOCK
! and J2 points to the last entry
         J1 = data%private%JPTR(LCNT)
         J2 = J1 + NVARSB(LCNT) - 1
         J1V = VPTR(LCNT)
         J2V = J1V  + NVARSB(LCNT) - 1
! Set NZBLK to number of entries in block
         NZBLK = NVARSB(LCNT) - data%private%IDUP(LCNT)

         IF (ICONTROL(7) == 1 .OR. &
            (ICONTROL(7) == 2 .AND. RANK == 0)) THEN
! Read in VALUES from sequential file.
            IF (ICONTROL(7) == 1) OPEN (UNIT=STRM,IOSTAT=IOS,ERR=300, &
                  FILE=VALNAM(LCNT),FORM='UNFORMATTED')
            IF (ICONTROL(7) == 2) OPEN (UNIT=STRM,IOSTAT=IOS,ERR=300, &
                  FILE=data%VALNAM(JBLOCK),FORM='UNFORMATTED')
            READ (STRM) data%VALUES(1:NVARSB(LCNT))
            CLOSE (UNIT=STRM)
            J1V = 1
            J2V = J1V + NVARSB(LCNT) - 1
          END IF

! We have to squeeze out any duplicates before we call MC46
         IF (data%private%IDUP(LCNT) == 0) THEN
! No duplicates in this block so take a copy of data
           IF (ICONTROL(7) == 5) THEN
              A(1:NZBLK) = data%RVAL(J1V:J2V)
           ELSE
              A(1:NZBLK) = data%VALUES(J1V:J2V)
           END IF
           IRN(1:NZBLK) = data%EQVAR(J1:J2)
           ROWPTR(1) = 1
           DO L = 1,M
             ROWPTR(L+1) = ROWPTR(L) + &
                           data%EQPTR(L1+L) - data%EQPTR(L1+L-1)
           END DO
         ELSE
! Duplicates to be got rid of.
           A(1:NZBLK) = ZERO
           ROWPTR(1) = 1
           KSTART = 1
           J = 0
! Loop over the rows
           IF (ICONTROL(7) == 5) THEN
             DO 266 L = 1,M
               NVAR = data%EQPTR(L1+L) - data%EQPTR(L1+L-1)
               KSTOP = KSTART + NVAR
               ROWPTR(L+1) = ROWPTR(L)
               DO 265 K = KSTART,KSTOP - 1
                 I = data%EQVAR(J1+K-1)
                 J = J + 1
                 JJ = data%private%IMAP(J,LCNT)
                 IF (IW1(I) < L+LMX) THEN
                    IW1(I) = L + LMX
                    ROWPTR(L+1) = ROWPTR(L+1) + 1
! Set col. index in IRN
                    IRN(JJ) = I
                 END IF
                 A(JJ) = A(JJ) + data%RVAL(J1V+K-1)
  265          CONTINUE
               KSTART = KSTOP
  266        CONTINUE
             LMX = LMX + NEQ
           ELSE
             DO 268 L = 1,M
               NVAR = data%EQPTR(L1+L) - data%EQPTR(L1+L-1)
               KSTOP = KSTART + NVAR
               ROWPTR(L+1) = ROWPTR(L)
               DO 267 K = KSTART,KSTOP - 1
                 I = data%EQVAR(J1+K-1)
                 J = J + 1
                 JJ = data%private%IMAP(J,LCNT)
                 IF (IW1(I) < L+LMX) THEN
                    IW1(I) = L + LMX
                    ROWPTR(L+1) = ROWPTR(L+1) + 1
! Set col. index in IRN
                    IRN(JJ) = I
                 END IF
                 A(JJ) = A(JJ) + data%VALUES(J1V+K-1)
  267          CONTINUE
               KSTART = KSTOP
  268        CONTINUE
             LMX = LMX + NEQ
           END IF
         END IF

! Allocate IW so it will also be sufficient for MP48_MA60B/BD
         DEALLOCATE (IW,STAT=ST)
         ALLOCATE (IW(1:M+2*NC),STAT=ST)
         IF (ST /= 0) GO TO 300
         DEALLOCATE (W,STAT=ST)
         ALLOCATE (W(1:M),STAT=ST)
         IF (ST /= 0) GO TO 300

        CALL MC46AD(M,NC,NZBLK,.TRUE.,IRN,A,ROWPTR,COLPTR,IW,IFLAG46)
! There should not be any errors from MC46
        IF (IFLAG46 < 0) THEN
          FLAG = -101
          GO TO 290
        END IF
! IRN now holds row indices, stored by cols.
! COLPTR holds col.allocatables.

! Ready for MP48_MA60B/BD

         CALL MP48_MA60BD(M,NC,NZBLK,JOB50,A,IRN,COLPTR, &
               data%private%CNTL_60(1:4,LCNT), &
               data%private%ICNTL_60(1:10,LCNT),data%private%IP(1:M,LCNT), &
               data%private%IQ(1:NC,LCNT),data%private%NP(JBLOCK),LF,NNA, &
               LFACT,data%FACT(I1F:I2F),LIRNF,data%IRNF(I1:I2), &
               data%private%IPTRL(1:NC,LCNT),data%private%IPTRU(1:NC,LCNT), &
               W,IW,data%INFO_MA60(1:20,JBLOCK),RINFO)

         data%OPS(JBLOCK) = RINFO

! We will want some of the info. returned by MP48_MA60B/BD
!    INFO_MA60(4) Minimum real storage required to factorize matrix.
!    INFO_MA60(6) Number of entries dropped from the data structure.
!    INFO_MA60(7) Number of rows processed in full storage.
!    INFO_MA60(8) Number of reals in factors.
!    INFO_MA60(9) Number of integers in factors.
!    INFO_MA60(10) Minimum integer storage required to factorize matrix.

         T2 = MPI_WTIME()
         data%TIMEF(JBLOCK) = data%TIMEF(JBLOCK) + T2 - T1

! Check for errors.
         FLAG = 0
! Only possible values are -7+K if JOB50 = 2
         IF (data%INFO_MA60(1,JBLOCK) < 0) THEN
           IF (data%INFO_MA60(1,JBLOCK) == -3) THEN
! Insufficent space
              LFACT = MAX(ICONTROL(6)*LFACT,data%INFO_MA60(4,JBLOCK))
              LIRNF = MAX(ICONTROL(6)*LIRNF,data%INFO_MA60(10,JBLOCK))
              data%private%LFACT(LCNT) = LFACT
              data%private%LIRNF(LCNT) = LIRNF
              FLAG = -3
           ELSE IF (JOB50 == 2 .AND. &
                     data%INFO_MA60(1,JBLOCK) < -7) THEN
              FLAG = -8
           ELSE
!             write (6,*) 'error',jblock,data%INFO_MA60(1,JBLOCK)
              FLAG = -103
           END IF
           GO TO 290
         END IF

! LF is the number of pivots chosen from the full part
! We will store this in data%INFO_MA60(11,JBLOCK)
         data%INFO_MA60(11,JBLOCK) = LF

!*******
! The following block of code is only for writing out the row
! number of the rows within the block that are not pivoted on.
! these rows form the "bottom" border ... information needed
! by Gene and Felix
!
!         MF = data%INFO_MA60(7,JBLOCK)
! MF is number of rows in full part
!         LF = data%INFO_MA60(11,JBLOCK)
! LF is number of pivots selected in full part
! We will only have rows left as part of border if MF - LF > 0
!         IF (data%FACT_JOB == 1) THEN
!          IF (MF - LF > 0) THEN
!            IF1I = data%private%IPTRL(NC,LCNT) + I1
!!            write (6,*) ' Rows in full part'
!!            write (6,'(8i4)')
!!     &      data%IRNF(IF1I:IF1I+MF-1)
!!            write (6,'(a,12i4)')
!!     &    ' IPIV. NC,NP,LF',NC,data%private%NP(JBLOCK),LF
!!            write (6,'(8i4)')
!!     &      data%IRNF(IF1I+MF:IF1I+MF+NC-data%private%NP(JBLOCK)-1)
!            DO I = 1,M
!               IW(I) = I
!            END DO
!            DO I = 1,LF
!               J = data%IRNF(IF1I+MF+I-1)
!               K = IW(I)
!               IW(I) = IW(J)
!               IW(J) = K
!            END DO
!            IOUT = 80 + JBLOCK
!!            write (6,'(12i4)') iw(1:m)
! Rows that have not been chosen as pivots in full part
! are rows iw(lf+1:mf) ... these are local indices to
! full part, but we want global indices
!!            write (6,'(a,12i4)') 'rows',data%IRNF(IF1I-1+IW(LF+1:MF))
!!            write (IOUT,'(a,i3,a)')
!!     &     'For block ',JBLOCK, ' the rows in the border are'
!            write (IOUT,'(12i6)') JBLOCK,MF-LF
!            write (IOUT,'(12i6)') data%IRNF(IF1I-1+IW(LF+1:MF))
!          ELSE
!!            write (6,'(a)') 'No rows left for border'
!            write (IOUT,'(12i6)') JBLOCK,MF-LF
!
!          END IF
!         END IF
! End of block of code
!*******

! LF holds the number of pivots chosen in the dense matrix.
! NP is the number of columns processed in packed storage.

! Optionally write out data held in FACT and IRNF to sequential
! files to reduce in-core storage requirements
         IF (ICONTROL(11) /= 0) THEN
            INFO4 = data%INFO_MA60(4,JBLOCK)
            INFO10 = data%INFO_MA60(10,JBLOCK)

            CLOSE (UNIT=STRM)
            OPEN (UNIT=STRM,IOSTAT=IOS,ERR=300, &
                  FILE=data%FILES(1,JBLOCK),STATUS='REPLACE', &
                  ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
            WRITE (STRM) data%FACT(1:INFO4)
            CLOSE (UNIT=STRM)
            OPEN (UNIT=STRM,IOSTAT=IOS,ERR=300, &
                  FILE=data%FILES(2,JBLOCK),STATUS='REPLACE', &
                  ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
            WRITE (STRM) data%IRNF(1:INFO10)
            CLOSE (UNIT=STRM)
         END IF

! Store data left after submatrix factorization complete.
! We do this before sending to host
         MF = data%INFO_MA60(7,JBLOCK)
         LF = data%INFO_MA60(11,JBLOCK)
         NF = data%private%NF(JBLOCK)
         II1 = data%private%IPTRL(NC,LCNT)
         L = I4
         MNF = (MF - LF)*(NF - LF)

! If L+MNF > SIZE(FACT_COPY) then we do not have enough room
! to save the Schur complement.
         IF (L + MNF > SIZE(FACT_COPY)) THEN
! Increase space for FACT_COPY
             DO STRM = 8,99
               IF (STRM == DGN .OR. STRM == ERR .OR. STRM == WRN) &
                CYCLE
               INQUIRE (UNIT=STRM,IOSTAT=IOS,ERR=300,EXIST=EX, &
                      OPENED=OPEN)
               IF (EX .AND. .NOT.OPEN) GO TO 272
             END DO
! No unit found. Jump to error return.
             FLAG = -20
             GO TO 300
 272        CONTINUE

             IF (SIZE(FACT_COPY) > 0) THEN
               INQUIRE (IOLENGTH=LENGTH) FACT_COPY
               OPEN (UNIT=STRM,IOSTAT=IOS,ERR=300,STATUS='SCRATCH', &
                     RECL=LENGTH,FORM='UNFORMATTED', &
                     ACTION='READWRITE')
               WRITE (UNIT=STRM) FACT_COPY
               REWIND (UNIT=STRM)
             END IF
             DEALLOCATE (FACT_COPY,STAT=ST)
             ALLOCATE (FACT_COPY(L+MNF),STAT=ST)
             IF (ST /= 0) GO TO 300
             IF (L > 1) THEN
               READ (UNIT=STRM) FACT_COPY(1:L-1)
               CLOSE(UNIT=STRM)
             END IF

         END IF

         DO 273 I = LF+1,MF
            DO J = LF+1,NF
              II2 = II1 + I1F + (J-1)*MF + I - 1
              FACT_COPY(L) = data%FACT(II2)
              L = L + 1
            END DO
  273    CONTINUE
         NE_SCHUR(JBLOCK) = L - I4
         I4 = I4 + MNF

! Generate modified IQ (ONLY needed if interface has more entries
! than just border cols.)
!        write (6,*) 'jblock,nf,lf,data%NGUARD(JBLOCK),mf',
!     &  jblock,nf,lf,data%NGUARD(JBLOCK),mf
         IF (NF-LF == data%NGUARD(JBLOCK) .OR. MF == 0) GO TO 280
         DO I = 1,NC
           IQM(I,LCNT) = I
         END DO
! Loop backwards over the border cols.
         DO 275 I = NC-data%private%NP(JBLOCK),1,-1
            KL = data%IRNF(II1+I1+MF+I-1)
            IF (KL > 0 .OR. KL == -I) GO TO 275
            II = IQM(I,LCNT)
            IQM(I,LCNT) = IQM(-KL,LCNT)
            IQM(-KL,LCNT) = II
  275    CONTINUE
!        write (6,*) 'iqm',iqm(1:nc,lcnt)
  280 CONTINUE

  290 CONTINUE
! 13 Dec 2002 added following
      IF (ICONTROL(11) /= 0) DEALLOCATE (data%FACT,data%IRNF,STAT=ST)

! NOTE: if files NOT used and LFACT and LIRNF have been increased,
! and not enough space to take temporary
! copy of factors generated so far, process will have to start again
      IF (FLAG == -3) GO TO 250

      DEALLOCATE (IW1,IW,W,STAT=ST)
      DEALLOCATE (A,COLPTR,ROWPTR,IRN,STAT=ST)

! On host, allocate IQ_COPY to which we will send the IQ for submatrices
      IF (RANK == 0) THEN
         NCMX  = 0
         DO JBLOCK = 1,NBLOCK
            NC = data%private%NCOLA(JBLOCK)
            NCMX = MAX(NCMX,NC)
         END DO
         DEALLOCATE (IW,STAT=ST)
         ALLOCATE (IW(NEQ),STAT=ST)
         IF (ST /= 0) GO TO 300
         DEALLOCATE (data%private%IQ_COPY, &
                     data%private%IQM_COPY, STAT=ST)
         ALLOCATE (data%private%IQ_COPY(1:NCMX+1,1:NBLOCK), &
                   data%private%IQM_COPY(1:NCMX+1,1:NBLOCK),STAT=ST)
         IF (ST /= 0) GO TO 300
         DEALLOCATE (data%STORAGE,STAT=ST)
         ALLOCATE (data%STORAGE(1:2,1:NBLOCK),STAT=ST)
         IF (ST /= 0) GO TO 300
      END IF

! Synchronize
  300 CALL MPI_BARRIER(data%private%COMM,ERCODE)
      IF (IOS /= 0) FLAG = -17
      data%IOSTAT = IOS
      data%STAT = ST
      CALL MP48ND(data,NPROC,IFLAG,FLAG)
 ! Broadcast error status from root
      CALL MPI_BCAST(data%ERROR,1,MPI_INTEGER,0,data%private%COMM,ERCODE)
! Terminate if error occurred
      IF (data%ERROR < 0) GO TO 890

! Collect information arrays and IQ on root.
      IF (RANK == 0) THEN
         MXINT = 0
         data%FLOPS = ZERO
         data%NZ = ZERO
         data%STORINT = ZERO
         data%DROP = 0
         TOTNE_SCHUR = 0
         DO IPROC = 1,NPROC-1
            SRC = IPROC
            ISTRT = data%IPLIST(IPROC)
            ISTOP = data%IPLIST(IPROC+1) - 1
            DO I = ISTRT,ISTOP
               JBLOCK = data%LIST(I)
               NC = data%private%NCOLA(JBLOCK)
               CALL MPI_RECV(data%INFO_MA60(1,JBLOCK),20,MPI_INTEGER, &
                             SRC,1,data%private%COMM,STAT,ERCODE)
               CALL MPI_RECV(data%OPS(JBLOCK),1, &
                             MPI_DOUBLE_PRECISION, &
                             SRC,2,data%private%COMM,STAT,ERCODE)
               CALL MPI_RECV(data%private%IQ_COPY(1,JBLOCK),NC, &
                             MPI_INTEGER,SRC,3,data%private%COMM,STAT,ERCODE)

! Accumulate statistics
               data%FLOPS = data%FLOPS + data%OPS(JBLOCK)
               data%DROP  = data%DROP  + data%INFO_MA60(6,JBLOCK)
               data%NZ    = data%NZ + REAL(data%INFO_MA60(8,JBLOCK),WP)
               data%STORINT = data%STORINT + &
                              REAL(data%INFO_MA60(9,JBLOCK),WP)
               data%STORAGE(1,JBLOCK) = &
                  REAL(data%INFO_MA60(8,JBLOCK),WP)
               data%STORAGE(2,JBLOCK) = &
                  REAL(data%INFO_MA60(9,JBLOCK),WP)

               MF = data%INFO_MA60(7,JBLOCK)
               LF = data%INFO_MA60(11,JBLOCK)
               NF = data%private%NF(JBLOCK)
               MNF = (MF - LF)*(NF - LF)
               CALL MPI_RECV(NE_SCHUR(JBLOCK),1,MPI_INTEGER, &
                             SRC,7,data%private%COMM,STAT,ERCODE)
               TOTNE_SCHUR = TOTNE_SCHUR + NE_SCHUR(JBLOCK)
               MXINT = MAX(MXINT,NF-LF)
               IF (NF-LF /= data%NGUARD(JBLOCK)) &
               CALL MPI_RECV(data%private%IQM_COPY(1,JBLOCK),NC, &
                             MPI_INTEGER,SRC,8,data%private%COMM,STAT,ERCODE)

               CALL MPI_RECV(data%TIMEF(JBLOCK),1,MPI_DOUBLE_PRECISION, &
                             SRC,9,data%private%COMM,STAT,ERCODE)
            END DO
         END DO
         DO LCNT = 1,ICOUNT_RANK
            JBLOCK = IBLOCK(LCNT)
! Accumulate statistics
            TOTNE_SCHUR = TOTNE_SCHUR + NE_SCHUR(JBLOCK)
            data%FLOPS = data%FLOPS + data%OPS(JBLOCK)
            data%DROP = data%DROP + data%INFO_MA60(6,JBLOCK)
            data%NZ = data%NZ + REAL(data%INFO_MA60(8,JBLOCK),WP)
            data%STORINT = data%STORINT + &
                              REAL(data%INFO_MA60(9,JBLOCK),WP)
            data%STORAGE(1,JBLOCK) = &
                  REAL(data%INFO_MA60(8,JBLOCK),WP)
            data%STORAGE(2,JBLOCK) = &
                  REAL(data%INFO_MA60(9,JBLOCK),WP)

            NC = data%private%NCOLA(JBLOCK)
            data%private%IQ_COPY(1:NC,JBLOCK) = &
                  data%private%IQ(1:NC,LCNT)
            LF = data%INFO_MA60(11,JBLOCK)
            NF = data%private%NF(JBLOCK)
            MXINT = MAX(MXINT,NF-LF)
            MF = data%INFO_MA60(7,JBLOCK)
            IF (NF-LF == data%NGUARD(JBLOCK) .OR. MF == 0) CYCLE
            data%private%IQM_COPY(1:NC,JBLOCK) = IQM(1:NC,LCNT)
         END DO

      ELSE

         I4 = 1
         DO LCNT = 1,ICOUNT_RANK
            JBLOCK = IBLOCK(LCNT)
            NC = data%private%NCOLA(JBLOCK)
            CALL MPI_SEND(data%INFO_MA60(1,JBLOCK),20,MPI_INTEGER, &
                          0,1,data%private%COMM,ERCODE)
            CALL MPI_SEND(data%OPS(JBLOCK),1, &
                          MPI_DOUBLE_PRECISION, &
                          0,2,data%private%COMM,ERCODE)
            CALL MPI_SEND(data%private%IQ(1,LCNT),NC,MPI_INTEGER, &
                          0,3,data%private%COMM,ERCODE)

            MF = data%INFO_MA60(7,JBLOCK)
            LF = data%INFO_MA60(11,JBLOCK)
            NF = data%private%NF(JBLOCK)
            MNF = (MF - LF)*(NF - LF)
            CALL MPI_SEND(NE_SCHUR(JBLOCK),1,MPI_INTEGER, &
                          0,7,data%private%COMM,ERCODE)
            I4 = I4 + MNF
            IF (NF-LF /= data%NGUARD(JBLOCK)) &
            CALL MPI_SEND(IQM(1,LCNT),NC,MPI_INTEGER, &
                          0,8,data%private%COMM,ERCODE)
            CALL MPI_SEND(data%TIMEF(JBLOCK),1,MPI_DOUBLE_PRECISION, &
                          0,9,data%private%COMM,ERCODE)
         END DO
      END IF
      DEALLOCATE(IQM,STAT=ST)

! Synchronize
      CALL MPI_BARRIER(data%private%COMM,ERCODE)
      IF (RANK == 0) THEN
         DEALLOCATE (FF_COPY,STAT=ST)
         ALLOCATE (FF_COPY(TOTNE_SCHUR),STAT=ST)
         IF (ST /= 0) CALL MP48LD(data,ST)
      END IF
      CALL MPI_BCAST(data%ERROR,1,MPI_INTEGER,0,data%private%COMM,ERCODE)
! Terminate if error occurred
      IF (data%ERROR < 0) GO TO 890
      IF (RANK == 0) THEN
         DO IPROC = 1,NPROC-1
            SRC = IPROC
            ISTRT = data%IPLIST(IPROC)
            ISTOP = data%IPLIST(IPROC+1) - 1
            DO I = ISTRT,ISTOP
               JBLOCK = data%LIST(I)
               JPOS = 1
               DO J = 1,JBLOCK-1
                  JPOS = JPOS + NE_SCHUR(J)
               END DO
               IF (NE_SCHUR(JBLOCK) /= 0) &
               CALL MPI_RECV(FF_COPY(JPOS),NE_SCHUR(JBLOCK), &
                             MPI_DOUBLE_PRECISION, &
                             SRC,1,data%private%COMM,STAT,ERCODE)
            END DO
         END DO
         I4 = 1
         DO LCNT = 1,ICOUNT_RANK
            JBLOCK = IBLOCK(LCNT)
            JPOS = 1
            DO J = 1,JBLOCK-1
               JPOS = JPOS + NE_SCHUR(J)
            END DO
            DO J = 1,NE_SCHUR(JBLOCK)
               FF_COPY(JPOS+J-1) = FACT_COPY(I4+J-1)
            END DO
            MF = data%INFO_MA60(7,JBLOCK)
            LF = data%INFO_MA60(11,JBLOCK)
            NF = data%private%NF(JBLOCK)
            MNF = (MF - LF)*(NF - LF)
            I4 = I4 + MNF
         END DO
      ELSE
         I4 = 1
         DO LCNT = 1,ICOUNT_RANK
            JBLOCK = IBLOCK(LCNT)
            MF = data%INFO_MA60(7,JBLOCK)
            LF = data%INFO_MA60(11,JBLOCK)
            NF = data%private%NF(JBLOCK)
            MNF = (MF - LF)*(NF - LF)
            IF (NE_SCHUR(JBLOCK) /= 0) &
            CALL MPI_SEND(FACT_COPY(I4),NE_SCHUR(JBLOCK), &
                          MPI_DOUBLE_PRECISION,0,1,data%private%COMM,ERCODE)
            I4 = I4 + MNF
         END DO
      END IF
      DEALLOCATE(FACT_COPY,STAT=ST)

! Synchronize
      CALL MPI_BARRIER(data%private%COMM,ERCODE)
      IF (RANK == 0) THEN
         IF (LDGN.GE.2) THEN
            WRITE (DGN,FMT='(/A/6(ES12.5))') &
          ' The flops performed for each submatrix and total is:', &
            data%OPS(1:NBLOCK),data%FLOPS
            WRITE (DGN,FMT='(A/6(ES12.5))') &
          ' The real factor storage for each submatrix is:', &
            data%STORAGE(1,1:NBLOCK)
            WRITE (DGN,FMT='(A/6(ES12.5))') &
          ' The integer factor storage for each submatrix is:', &
            data%STORAGE(2,1:NBLOCK)
            IF (data%DROP > 0) WRITE (DGN,FMT='(A/6(I8))') &
          ' The number of dropped entries for each submatrix is:', &
            data%INFO_MA60(6,1:NBLOCK)
            IF (LDGN.GE.3) WRITE (DGN,FMT='(A/8(F9.3))') &
      ' The factor wall clock time (secs) for each submatrix is:', &
            (data%TIMEF(JBLOCK),JBLOCK=1,NBLOCK)
         END IF
!!! Following commented out since we only wish to allow host write out
! We only need the following if the number of submatrices is
! not equal to the number of processes
!!         IF (NBLOCK /= NPROC) THEN
!!          DO IPROC = 0,NPROC-1
!!            TIM60 = ZEROD
!!            RINFO = ZERO
!!            ISTRT = data%IPLIST(IPROC)
!!            ISTOP = data%IPLIST(IPROC+1) - 1
!!            DO I = ISTRT,ISTOP
!!               JBLOCK = data%LIST(I)
!!               TIM60 = TIM60 + data%TIMEF(JBLOCK)
!!               RINFO = RINFO + data%OPS(JBLOCK)
!!            END DO
!!            IF (LDGN.GE.3) THEN
!!               WRITE (DGN,FMT='(A,I2,A,ES12.5,2X,F9.3)')
!!     &        ' Processor ', IPROC,
!!     &  ' flop count and wall clock time (secs) for factorizations:',
!!     &          RINFO,TIM60
!!            ELSE IF (LDGN.GE.2) THEN
!!               WRITE (DGN,FMT='(A,I3,A,ES12.5)')
!!     &       ' For processor ', IPROC,
!!     &       ' the flop count for factorizations:',RINFO
!!            END IF
!!          END DO
!!         END IF
!!!
      END IF

! Broadcast MXINT (max border size) from root
      CALL MPI_BCAST(MXINT,1,MPI_INTEGER,0,data%private%COMM,ERCODE)

      DEALLOCATE (data%private%JFVAR,STAT=ST)
      ALLOCATE (data%private%JFVAR(1:MXINT,1:NBLOCK),STAT=ST)
      IF (ST /= 0) GO TO 330

      IF (RANK == 0) THEN
         IW(1:NEQ) = 0
         L1 = 0
         DO 320 JBLOCK = 1,NBLOCK
            LF = data%INFO_MA60(11,JBLOCK)
            NF = data%private%NF(JBLOCK)
            NP = data%private%NP(JBLOCK)
            MF = data%INFO_MA60(7,JBLOCK)
            IF (NF-LF==data%NGUARD(JBLOCK) .OR. MF==0) THEN
             DO 310 L = 1,NF-LF
! Cols that are left are the last cols. in the border
! KL is local index
               J  = NP + LF
               KL = data%private%IQ_COPY(J+L,JBLOCK)
! Global index corresponding to local index KL is held in CGLOB(KL)
               K  = data%private%CGLOB(KL,JBLOCK)
               IF (IW(K) == 0) THEN
! First appearance of col. K
                  L1 = L1 + 1
! Give border col. with global index K the local index L1
                  IW(K) = L1
               END IF
  310       CONTINUE
           ELSE
            DO 315 L = 1,NF-LF
! Cols that are left are the last cols. in the border
! KL is local index
               J  = NP + LF
               KL = data%private%IQ_COPY(J+L,JBLOCK)
! Global index corresponding to local index KL is held in CGLOB(KL)
               KL = data%private%IQM_COPY(KL,JBLOCK)
               K  = data%private%CGLOB(KL,JBLOCK)
               IF (IW(K) == 0) THEN
! First appearance of col. K
                  L1 = L1 + 1
! Give border col. with global index K the local index L1
                  IW(K) = L1
               END IF
  315       CONTINUE
            END IF
  320    CONTINUE

! No. of cols. in interface problem is L1+NULL
! If this is a JOB50 = 2 call, we must check that this
! is the same as on the earlier call; otherwise we cannot proceed.
         IF (JOB50 == 2 .AND. L1+data%NULL /= data%NINTER) THEN
            data%ERROR = -15
            IF (LERR) THEN
               WRITE (ERR,FMT=9000) data%JOB,data%ERROR
               WRITE (ERR,FMT=9045)
            END IF
            GO TO 330
         END IF

         data%NINTER = L1 + data%NULL

! Check size of interface compared to size of border
         IF (data%BORDER /= data%NINTER) THEN
            IF (data%ERROR < 2) THEN
               data%ERROR = data%ERROR + 2
               IF (LWRN) THEN
                  WRITE (WRN,FMT=9010) data%JOB,data%ERROR
                  WRITE (WRN,FMT=9080)
               END IF
            END IF
         END IF

! List of columns belonging to interface is to be held in INTER_COL
         DEALLOCATE (data%INTER_COL,STAT=ST)
         ALLOCATE (data%INTER_COL(data%NINTER),STAT=ST)
         IF (ST /= 0) GO TO 330
      END IF

  330 CALL MPI_BARRIER(data%private%COMM,ERCODE)
      CALL MP48MD(data,NPROC,IFLAG,ST)
! Broadcast error status from root
      CALL MPI_BCAST(data%ERROR,1,MPI_INTEGER,0,data%private%COMM,ERCODE)
! Terminate if error occurred
      IF (data%ERROR < 0) GO TO 890

! Set JFVAR(:,JBLOCK) to hold a list of the columns left for block
! JBLOCK, renumbered so that border columns are 1,2,...,NINTER

      IF (RANK == 0) THEN
         IW(1:NEQ) = 0
         L1 = 0
         DO 340 JBLOCK = 1,NBLOCK
            LF = data%INFO_MA60(11,JBLOCK)
            NF = data%private%NF(JBLOCK)
            NP = data%private%NP(JBLOCK)
            MF = data%INFO_MA60(7,JBLOCK)
            IF (NF-LF==data%NGUARD(JBLOCK) .OR. MF==0) THEN
              DO 331 L = 1,NF-LF
! Cols that are left are the last cols. in the border
! KL is local index
                J  = NP + LF
                KL = data%private%IQ_COPY(J+L,JBLOCK)
! Global index corresponding to local index KL is held in CGLOB(KL)
                K  = data%private%CGLOB(KL,JBLOCK)
!                write (6,*)'lf,j,jblock,k',lf,j,jblock,k
                IF (IW(K) == 0) THEN
! First appearance of col. K
                  L1 = L1 + 1
! Give border col. with global index K the local index L1
                  IW(K) = L1
                  data%private%JFVAR(L,JBLOCK) = L1
                  data%INTER_COL(L1) = K
!               write (6,*)'jblock,l,l1,k',jblock,l,l1,k
                ELSE
! Col. K has already appeared ... its local index is held in IW(K)
                  data%private%JFVAR(L,JBLOCK) = IW(K)
                END IF
  331        CONTINUE
           ELSE
! This is the when we have extra cols left (ie original border
! cols plus cols. not eliminated for stability reasons). In this
! case we have to use data%private%IQM_COPY
              DO 332 L = 1,NF-LF
! Cols that are left are the last cols. in the border
! KL is local index
                J  = NP + LF
                KL = data%private%IQ_COPY(J+L,JBLOCK)
!                  write (6,*) 'l,j+l,kl',l,j+l,kl
                KL = data%private%IQM_COPY(KL,JBLOCK)
! Global index corresponding to local index KL is held in CGLOB(KL)
                K  = data%private%CGLOB(KL,JBLOCK)
!               write (6,'(a,8i4)')
!     &        'iqm,l,lf,nf,j,jblock,kl,k',l,lf,nf,j,jblock,kl,k
               IF (IW(K) == 0) THEN
! First appearance of col. K
                  L1 = L1 + 1
! Give border col. with global index K the local index L1
                  IW(K) = L1
                  data%private%JFVAR(L,JBLOCK) = L1
                  data%INTER_COL(L1) = K
!                 write (6,*)'iqm,jblock,l,l1,k',jblock,l,l1,k
                ELSE
! Col. K has already appeared ... its local index is held in IW(K)
                  data%private%JFVAR(L,JBLOCK) = IW(K)
                END IF
  332        CONTINUE
           END IF
  340    CONTINUE
         IF (data%NULL /= 0) THEN
           DO I = 1,NEQ
             IF (data%private%INULL(I) == 0) THEN
! Col I is null
               L1 = L1 + 1
               data%INTER_COL(L1) = I
             END IF
           END DO
         END IF
         DEALLOCATE (IW,STAT=ST)
         L1 = 0
         DO JBLOCK = 1,NBLOCK
            MF = data%INFO_MA60(7,JBLOCK)
            LF = data%INFO_MA60(11,JBLOCK)
            L1 = L1 + (MF-LF)
         END DO

! No. rows in interface
         data%RINTER = L1
      END IF


! Broadcast JFVAR and data%NINTER from root
      CALL MPI_BCAST(data%private%JFVAR,MXINT*NBLOCK,MPI_INTEGER,0, &
                     data%private%COMM,ERCODE)
      CALL MPI_BCAST(data%NINTER,1,MPI_INTEGER,0,data%private%COMM,ERCODE)
      CALL MPI_BCAST(data%RINTER,1,MPI_INTEGER,0,data%private%COMM,ERCODE)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! At this point factorization of submatrices is complete.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      data%EST_RANK = data%NEQ
! data%EST_RANK holds estimate of rank of A. It is equal to NEQ
! unless interface problem has rank less then NINTER.

! Solve interface problem using MA48 on the host only.
! If the subproblems are totally disjoint then
! there will be no interface problem to solve so we
! should not try to call MA48.
      IF (data%NINTER == 0) GO TO 890
      IF (RANK == 0) THEN

       data%TIME_MA48A = ZERO
       data%TIME_MA48F = ZERO
       T1 = MPI_WTIME()
! We only want to call initialise and analyse if
! this is the first call with a matrix having given sparsity
! pattern (data%private%FJOB = 1)
       IF (data%private%FJOB == 1) THEN
! First call.
         CALL MA48ID(data%private%CNTL_48,data%private%ICNTL_48)

! Switch off printing from MA48
         data%private%ICNTL_48(1) = -1
         data%private%ICNTL_48(2) = -1
         data%private%ICNTL_48(3) = -1
!        data%private%ICNTL_48(3) =  4

! ICNTL_48(4): the pivot search is limited to ICNTL_60(4,LCNT)
! cols (Zlatev strategy). Default is 3.
         data%private%ICNTL_48(4) = ICONTROL(8)

         data%private%ICNTL_48(5) = ICONTROL(5)
! Handle rank deficient
         data%private%ICNTL_48(7) = 0

! Switch control from sparse to full matrix processing (default 0.5)
! Smaller value will speed up MA48AD but means higher flop count
! and so MA48BD may be slower
         data%private%CNTL_48(1) = CONTROL(1)
! Set threshold pivoting parameter
         data%private%CNTL_48(2) = CONTROL(2)

         data%private%CNTL_48(3:4) = CONTROL(3:4)

! Find number of entries in interface matrix
         NE = 0
         DO JBLOCK = 1,NBLOCK
           MF = data%INFO_MA60(7,JBLOCK)
           LF = data%INFO_MA60(11,JBLOCK)
           NF = data%private%NF(JBLOCK)
           NE = NE + (MF - LF)*(NF - LF)
         END DO
! Arrays of length LA ... may choose LA too small
         LA = MAX(8*NE,100)
       END IF

  350  NINTER = data%NINTER
       RINTER = data%RINTER

       IF (data%private%FJOB == 1) THEN
         LKEEP48 = RINTER+5*NINTER + 4*NINTER/data%private%ICNTL_48(6)+7
         DEALLOCATE (IW,STAT=ST)
         DEALLOCATE (data%private%A,STAT=ST)
         DEALLOCATE (data%private%IRN,STAT=ST)
         DEALLOCATE (data%private%JCN,STAT=ST)
         DEALLOCATE (data%private%KEEP48,STAT=ST)

         ALLOCATE (IW(1:6*RINTER+3*NINTER), &
                   data%private%A(1:LA), &
                   data%private%IRN(1:LA), &
                   data%private%JCN(1:LA), &
                   data%private%KEEP48(1:LKEEP48), &
                   STAT=ST)
         IF (ST /= 0) GO TO 420
! Ensure KEEP48 is defined (otherwise I get an error later when
! I try to take a copy)
           data%private%KEEP48(1:LKEEP48) = ZERO

         ELSE
! For data%private%FJOB > 1 call, IRN and JCN must be unchanged since
! the call to MA48A/AD, but A is reset.
! Note: we do not call analyse phase MA48A/AD.

           LA = data%private%LA48
           LKEEP48 = data%private%LKEEP48
           DEALLOCATE (data%private%A,STAT=ST)
           ALLOCATE (data%private%A(1:LA),STAT=ST)
           IF (ST /= 0) GO TO 420
         END IF

! Copy entries to IRN, JCN and A
         L1 = 0
         L2 = 0
         I2 = 1
         DO 370 JBLOCK = 1,NBLOCK
! Need to loop over interface rows.
! IRNF holds row indices
            MF = data%INFO_MA60(7,JBLOCK)
            LF = data%INFO_MA60(11,JBLOCK)
            NF = data%private%NF(JBLOCK)
!           write (6,*) 'mf,lf,nf',jblock,mf,lf,nf
            IF (data%private%FJOB == 1) THEN
              DO I = 1,MF-LF
                L1 = L1 + 1
                DO J = 1,NF-LF
                   IF (FF_COPY(I2) /= ZERO) THEN
                     L2 = L2 + 1
                     data%private%IRN(L2) = L1
                     data%private%JCN(L2) = data%private%JFVAR(J,JBLOCK)
                     data%private%A(L2)   = FF_COPY(I2)
!                    write (6,'(a,5i4,g11.4)') 'l2',l2,l1,j,jblock,
!     &              data%private%JFVAR(J,JBLOCK),FACT_COPY(I2,JBLOCK)
                   END IF
                   I2 = I2 + 1
                END DO
              END DO
            ELSE
! Subsequent factorize. Must set entries of A
              DO I = 1,MF-LF
                DO J = 1,NF-LF
                   IF (FF_COPY(I2) /= ZERO) THEN
                     L2 = L2 + 1
                     data%private%A(L2) = FF_COPY(I2)
                   END IF
                   I2 = I2 + 1
                END DO
              END DO
            END IF
  370    CONTINUE

         IF (data%private%FJOB == 1) THEN
            data%NE_INTER = L2
         ELSE IF (L2 /= data%NE_INTER) THEN
! Fast factor ... if L2.ne.data%NE_INTER then different
! number of nonzeros left in the interface and so we can't proceed
! (user will have to restart and recall JOB = 3)
            data%ERROR = -15
            IF (LERR) THEN
               WRITE (ERR,FMT=9000) data%JOB,data%ERROR
               WRITE (ERR,FMT=9040)
            END IF
            GO TO 420
         END IF
         NE = data%NE_INTER

! If NE = 0 (ie no nonzero entries in interface matrix) then MA48 will
! fail ... could happen in singular case.
! Of course, not likely in a real problem.
         IF (NE == 0) GO TO 420

! Only call analyse phase if this is a first call
! (data%private%FJOB = 1)
         IF (data%private%FJOB == 1) THEN
           JOB48 = 1
           CALL MA48AD(RINTER,NINTER,NE,JOB48,LA,data%private%A(1:LA), &
                data%private%IRN(1:LA),data%private%JCN(1:LA), &
                data%private%KEEP48(1:LKEEP48), &
                data%private%CNTL_48,data%private%ICNTL_48, &
                IW,data%INFO_MA48,data%RINFO_MA48)

! Check for errors
           IF (data%INFO_MA48(1) < 0) THEN
             IF (data%INFO_MA48(1) == -3) THEN
! Insufficient space ... reset LA. INFO(4) is suggested value
               LA = MAX(ICONTROL(6)*LA,data%INFO_MA48(4))
               GO TO 350
             ELSE IF (data%INFO_MA48(1) == -4) THEN
! Rank deficient
               data%ERROR = -21
               IF (LERR) THEN
                  WRITE (ERR,FMT=9000) data%JOB,data%ERROR
                  WRITE (ERR,FMT=9150)
               END IF
               GO TO 420
             ELSE
               data%ERROR = -102
               IF (LERR) THEN
                  WRITE (ERR,FMT=9000) data%JOB,data%ERROR
                  WRITE (ERR,FMT=9050)
                  WRITE (ERR,FMT=9060) data%INFO_MA48(1)
               END IF
               GO TO 420
             END IF
           END IF

! Take copy of KEEP as it must be unchanged between calling MA48A/AD
! and MA48B/BD (in case we have to recall MA48B/BD because
! of insufficient space). Can only happen if JOB50=1.
           IF (JOB50 == 1) THEN
              DEALLOCATE (KPCOPY,STAT=ST)
              ALLOCATE (KPCOPY(1:LKEEP48),STAT=ST )
              IF (ST /= 0) GO TO 420

              KPCOPY(1:LKEEP48) = data%private%KEEP48(1:LKEEP48)

! Check that LA is at least INFO(4).
! If not, we must increase length of A and IRN
               IF (LA < data%INFO_MA48(4)) THEN
                  LA = data%INFO_MA48(4)
                  DEALLOCATE (data%private%A,data%private%IRN, &
                              STAT=ST)
                  ALLOCATE (data%private%A(1:LA), &
                            data%private%IRN(1:LA),STAT=ST )
                  IF (ST /= 0) GO TO 420

! Reset first NE entries of A and IRN
                  L1 = 0
                  L2 = 0
                  I2 = 1
                  DO 375 JBLOCK = 1,NBLOCK
                    MF = data%INFO_MA60(7,JBLOCK)
                    LF = data%INFO_MA60(11,JBLOCK)
                    NF = data%private%NF(JBLOCK)
                    DO I = 1,MF-LF
                      L1 = L1 + 1
                      DO J = 1,NF-LF
                         IF (FF_COPY(I2) /= ZERO) THEN
                           L2 = L2 + 1
                           data%private%IRN(L2) = L1
                           data%private%A(L2) = FF_COPY(I2)
                         END IF
                         I2 = I2 + 1
                      END DO
                    END DO
  375            CONTINUE

               END IF
            END IF

           T2 = MPI_WTIME()
           data%TIME_MA48A = T2 - T1

         END IF

         T1 = MPI_WTIME()

         DEALLOCATE (W,STAT=ST)
         ALLOCATE (W(1:RINTER),STAT=ST )
         DEALLOCATE (IW,STAT=ST)
         ALLOCATE (IW(1:2*RINTER+2*NINTER),STAT=ST )
         IF (ST /= 0) GO TO 420

  380      CALL MA48BD(RINTER,NINTER,NE,JOB50,LA,data%private%A(1:LA), &
                data%private%IRN(1:LA),data%private%JCN(1:NE), &
                data%private%KEEP48(1:LKEEP48),data%private%CNTL_48, &
                data%private%ICNTL_48,W,IW, &
                data%INFO_MA48,data%RINFO_MA48)

         IF (data%INFO_MA48(1) < 0) THEN
            IF (JOB50 == 1 .AND. data%INFO_MA48(1) == -3) THEN
! Insufficient space ... reset LA. INFO(4) is suggested value
! This can only happen on JOB50 = 1 call.
               LA = MAX(ICONTROL(6)*LA,data%INFO_MA48(4))
               DEALLOCATE (data%private%A,data%private%IRN, &
                           STAT=ST)
               ALLOCATE (data%private%A(1:LA), &
                         data%private%IRN(1:LA),STAT=ST )
               IF (ST /= 0) GO TO 420

               data%private%KEEP48(1:LKEEP48) = KPCOPY(1:LKEEP48)
! Reset first NE entries of A and IRN
                L1 = 0
                L2 = 0
                I2 = 1
                DO 385 JBLOCK = 1,NBLOCK
                  MF = data%INFO_MA60(7,JBLOCK)
                  LF = data%INFO_MA60(11,JBLOCK)
                  NF = data%private%NF(JBLOCK)
                  DO I = 1,MF-LF
                    L1 = L1 + 1
                    DO J = 1,NF-LF
                       IF (FF_COPY(I2) /= ZERO) THEN
                         L2 = L2 + 1
                         data%private%IRN(L2) = L1
                         data%private%A(L2) = FF_COPY(I2)
                       END IF
                       I2 = I2 + 1
                    END DO
                  END DO
  385          CONTINUE
               GO TO 380
            ELSE IF (JOB50 == 2 .AND. data%INFO_MA48(1) < -7) THEN
               data%ERROR = -15
               IF (LERR) THEN
                  WRITE (ERR,FMT=9000) data%JOB,data%ERROR
                  WRITE (ERR,FMT=9030)
               END IF
               GO TO 420
            ELSE
               data%ERROR = -102
               IF (LERR) THEN
                  WRITE (ERR,FMT=9000) data%JOB,data%ERROR
                  WRITE (ERR,FMT=9055)
                  WRITE (ERR,FMT=9060) data%INFO_MA48(1)
               END IF
               GO TO 420
            END IF

         ELSE IF (data%INFO_MA48(1) == 2) THEN
! Singular with estimated rank INFO_MA48(5)
! Compute estimated rank of A
            data%EST_RANK = data%NEQ - (NINTER - data%INFO_MA48(5))
            IF (data%ERROR < 2) THEN
               data%ERROR = data%ERROR + 2
               IF (LWRN) THEN
                  WRITE (WRN,FMT=9010) data%JOB,data%ERROR
                  WRITE (WRN,FMT=9080)
               END IF
            END IF
         END IF

         data%private%LA48 = LA
         data%private%LKEEP48 = LKEEP48

! MA48B/BD has been successful. If we are using files for data
! (ICONTROL(11) /= 0) then we need to write out data%private%A
! and data%private%IRN
         IF (ICONTROL(11) /= 0) THEN

            CLOSE (UNIT=STRM)
            OPEN (UNIT=STRM,IOSTAT=IOS,ERR=415, &
                  FILE=data%FILES(1,NBLOCK+1),STATUS='REPLACE', &
                  ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
            WRITE (STRM) data%private%A
            CLOSE (UNIT=STRM)
            OPEN (UNIT=STRM,IOSTAT=IOS,ERR=415, &
                  FILE=data%FILES(2,NBLOCK+1),STATUS='REPLACE', &
                  ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
            WRITE (STRM) data%private%IRN
            CLOSE (UNIT=STRM)
            GO TO 420
  415       data%ERROR = -17
            IF (LERR) THEN
               WRITE (ERR,FMT=9000) data%JOB,data%ERROR
               WRITE (ERR,FMT=9140)
            END IF
            DEALLOCATE (data%private%A,data%private%IRN,STAT=ST)
         END IF

      END IF
! Broadcast error status from root
  420 CALL MPI_BARRIER(data%private%COMM,ERCODE)
      IF (RANK == 0 .AND. ST /= 0) CALL MP48LD(data,ST)
      CALL MPI_BCAST(data%ERROR,1,MPI_INTEGER,0,data%private%COMM,ERCODE)
! Terminate if error occurred
      IF (data%ERROR < 0) GO TO 890

      T2 = MPI_WTIME()
      data%TIME_MA48F = T2 - T1

! Interface problem now factorized.
! Compute (and optionally print) statistics for interface factorization.

      IF (RANK == 0) THEN
        IF (NE /= 0) THEN
         CALL MP48ED(RINTER,NINTER,data%private%KEEP48(1:LKEEP48),LKEEP48,JNFO)
! Real storage for factors is held in JNFO(5)
! Integer storage = real storage
         data%NZ      =  data%NZ + REAL(JNFO(5),WP)
         data%STORINT =  data%STORINT + REAL(JNFO(5),WP)
         data%FLOPS    = data%FLOPS + data%RINFO_MA48(1)
         data%DROP     = data%DROP + data%INFO_MA48(6)
        ELSE
! MA48 not run (no entries in interface)
         JNFO(5) = 0
         data%RINFO_MA48(1) = ZERO
        END IF
        IF (LDGN.GE.3) THEN
          IF (data%private%FJOB == 1) THEN
            WRITE (DGN,'(/A/A,2I8/A,I8,3(/A,ES12.5)/A,F9.3/A,F9.3)') &
          ' Interface problem: ', &
          ' Number of rows/columns            = ', &
            data%RINTER,data%NINTER, &
          ' Number of entries                 = ', data%NE_INTER, &
          ' Flops                             = ', data%RINFO_MA48(1), &
          ' Factor  storage                   = ', REAL(JNFO(5)), &
          ' Integer storage                   = ', REAL(JNFO(5)), &
          ' Analyse   wall clock time (secs)  = ', data%TIME_MA48A, &
          ' Factorize wall clock time (secs)  = ', data%TIME_MA48F
          ELSE
            WRITE (DGN,'(/A/A,2I8/A,I8,3(/A,ES12.5)/A,F9.3)') &
          ' Interface problem: ', &
          ' Number of rows/columns            = ', &
            data%RINTER,data%NINTER, &
          ' Number of entries                 = ', data%NE_INTER, &
          ' Flops                             = ', data%RINFO_MA48(1), &
          ' Factor  storage                   = ', REAL(JNFO(5)), &
          ' Integer storage                   = ', REAL(JNFO(5)), &
          ' Factorize wall clock time (secs)  = ', data%TIME_MA48F
          END IF
        ELSE IF (LDGN.GE.2) THEN
            WRITE(DGN,'(/A/A,2I8/A,I8,3(/A,ES12.5))') &
          ' Interface problem: ', &
          ' Number of rows/columns  = ',   data%RINTER,data%NINTER, &
          ' Number of entries       = ',   data%NE_INTER, &
          ' Flops                   = ',   data%RINFO_MA48(1), &
          ' Factor storage          = ',   REAL(JNFO(5)), &
          ' Integer storage         = ',   REAL(JNFO(5))
        END IF

        IF (LDGN.GE.2) THEN
           WRITE (DGN,'(/A,3(/A,ES12.5),/A,I8 )') &
          ' Global statistics: ', &
          ' Total flops                           = ', &
            data%FLOPS, &
          ' Total real factor storage             = ', &
            data%NZ, &
          ' Total integer factor storage          = ', &
            data%STORINT, &
          ' Estimated rank of matrix              = ', &
            data%EST_RANK
        END IF
      END IF

      CALL MPI_BCAST(data%NE_INTER,1,MPI_INTEGER,0,data%private%COMM,ERCODE)
! Return to user for right hand sides
  890 CALL MPI_BARRIER(data%private%COMM,ERCODE)
      CALL MPI_BCAST(data%ERROR,1,MPI_INTEGER,0,data%private%COMM,ERCODE)
      data%private%ERROR = data%ERROR

      DEALLOCATE (A,STAT=ST)
      DEALLOCATE (W,STAT=ST)
      DEALLOCATE (VALUES,STAT=ST)
      DEALLOCATE (FACT_COPY,STAT=ST)

      DEALLOCATE (COLPTR,STAT=ST)
      DEALLOCATE (IRN,STAT=ST)
      DEALLOCATE (IBLOCK,STAT=ST)
      DEALLOCATE (IFLAG,STAT=ST)
      DEALLOCATE (IW,STAT=ST)
      DEALLOCATE (IW1,STAT=ST)
      DEALLOCATE (KPCOPY,STAT=ST)
      DEALLOCATE (NEQSUB,STAT=ST)
      DEALLOCATE (NVARSB,STAT=ST)
      DEALLOCATE (ROWPTR,STAT=ST)
      DEALLOCATE (VPTR,STAT=ST)
      DEALLOCATE (NE_SCHUR,STAT=ST)
      DEALLOCATE (IQM,STAT=ST)

      IF (ICONTROL(7) == 1 .OR. ICONTROL(7) == 2) &
         DEALLOCATE (data%VALUES, STAT = ST)
      IF (ICONTROL(7) == 3 .AND. RANK /= 0) &
         DEALLOCATE (data%VALUES, STAT = ST)

 9000 FORMAT (/' Error return with data%JOB = ',I2,'. ', &
                'data%ERROR = ',I4)
 9010 FORMAT (/' Warning return with data%JOB = ',I2,'. ', &
                'data%ERROR = ',I4)
 9020 FORMAT (' ',A,' has not been allocated')
 9030 FORMAT (' Small pivot found during factorization of interface')
 9040 FORMAT (' Number of entries in interface matrix differs from the', &
              ' earlier call with data%JOB = 4')
 9045 FORMAT (' Order of interface matrix differs from the', &
              ' earlier call with data%JOB = .')
 9050 FORMAT (' Unexpected error return from MA48A/AD')
 9055 FORMAT (' Unexpected error return from MA48B/BD')
 9060 FORMAT (' MA48 error flag = ',I3)
 9080 FORMAT (' The matrix is singular')
 9100 FORMAT (' ',A,' is of size ',I10,' Increase to at least ',I10)
 9110 FORMAT (' Size of ',A,' is too small on one or more process')
 9120 FORMAT (' First dimension of ',A,' is of extent ',I10,/ &
              ' Increase to at least ',I10)
 9130 FORMAT (' Second dimension of ',A,' is of extent ',I10,/ &
              ' Increase to at least ',I10)
 9140 FORMAT (' Failure in OPEN statement.')
 9150 FORMAT (' Interface problem is structurally rank deficient.')
      END SUBROUTINE MP48CD
!*************************************************

      SUBROUTINE MP48DD(data)

! This subroutine performs solve phase using MP48_MA60CD. JOB = 5 call.
      USE HSL_MP48_DATA_DOUBLE
      IMPLICIT NONE

      REAL(WP), PARAMETER :: ZERO = 0.0_WP
      REAL(DP), PARAMETER :: ZEROD = 0.0_DP

      TYPE (MP48_DATA) data

!     .. Local Scalars ..

      REAL(DP) T1,T2,T3,T4
      INTEGER DGN,ERCODE,ERR,FLAG, &
              I,I1,I1F,I2,I2F,IPROC,ISTRT,ISTOP,ICOUNT_RANK, &
              J,JBLOCK,JOB48, &
              K, &
              L,L1,LA,LDGN,LFACT,LIRNF,LCNT,LORDER, &
              LORU,LX,LKEEP48,LF, &
              M,MMX,MF, &
              NBLOCK,NEQ,NPROC,NP,NC,NINTER,RINTER,NF, &
              RANK,STRM,SRC,WRN

      INTEGER IOS
!  Holds IOSTAT parameter for files.
      INTEGER ST
!  Holds STAT parameter (ALLOCATE statement)
      LOGICAL LERR,LWRN
      LOGICAL EX
!  Used by INQUIRE. Set to .TRUE. if unit exists.
      LOGICAL OPEN
!  Used by INQUIRE. Set to .TRUE. if file open.

!     .. Local Arrays ..

      REAL(WP) CONTROL(30),ERR48(3)

      REAL(WP), DIMENSION (:), ALLOCATABLE :: W
      REAL(WP), DIMENSION (:), ALLOCATABLE :: WLOCAL
      REAL(WP), DIMENSION (:), ALLOCATABLE :: X1
      REAL(WP), DIMENSION (:,:), ALLOCATABLE :: B1
      REAL(WP), DIMENSION (:,:), ALLOCATABLE :: XLOCAL

      INTEGER ICONTROL(20),STAT(MPI_STATUS_SIZE,20)

      INTEGER, DIMENSION (:), ALLOCATABLE :: IBLOCK
      INTEGER, DIMENSION (:), ALLOCATABLE :: IFLAG
      INTEGER, DIMENSION (:), ALLOCATABLE :: IW
      INTEGER, DIMENSION (:), ALLOCATABLE :: NEQSUB
      INTEGER, DIMENSION (:), ALLOCATABLE :: NVARSB

!     .. External Subroutines ..
      EXTERNAL MA48CD

!     .. MPI routines
      EXTERNAL MPI_BARRIER,MPI_BCAST,MPI_COMM_RANK, &
               MPI_GATHER,MPI_SEND,MPI_RECV

!     .. Intrinsic Functions ..
      INTRINSIC MAX

      FLAG = 0
      RANK = data%RANK
      NPROC = data%NPROC

! Set diagnostics controls
      DGN = data%ICNTL(3)
      LDGN = data%ICNTL(4)
      IF (DGN < 0) LDGN = 0
! Set error controls
      LERR = .FALSE.
      IF (data%ICNTL(1).GE.0 .AND. data%ICNTL(4).GE.1) LERR = .TRUE.
      ERR = data%ICNTL(1)
! Set warning controls
      LWRN = .FALSE.
      IF (data%ICNTL(2).GE.0 .AND. data%ICNTL(4).GE.1) LWRN = .TRUE.
      WRN = data%ICNTL(2)

      ICONTROL(1:20) =  data%ICNTL(1:20)
      CONTROL(1:4)   =  data%CNTL(1:4)

      NBLOCK = data%private%NBLOCK
      NEQ    = data%private%NEQ
      LORDER = data%private%LORDER

      IF (RANK == 0) THEN
! Check data for errors. Check required arrays have been allocated.
! Check size of array B
! First check it has been allocated
         IF (.NOT. allocated(data%B)) THEN
            FLAG = -22
         ELSE IF (SIZE(data%B) < NEQ) THEN
            FLAG = -22
         END IF
         IF (FLAG == -22) THEN
            data%ERROR = -22
            IF (LERR) THEN
               WRITE (ERR,FMT=9000) data%JOB,data%ERROR
               IF (.NOT. allocated(data%B)) THEN
                  WRITE (ERR,FMT=9020) ' data%B'
               ELSE IF (SIZE(data%B) < NEQ) THEN
                  WRITE (ERR,FMT=9030) 'data%B',SIZE(data%B),NEQ
               END IF
            END IF
            GO TO 20
         END IF
! Check size of array X (ICNTL(13) nonzer only)
! First check it has been allocated
         IF (data%ICNTL(13) /= 0) THEN
            IF (.NOT. allocated(data%X)) THEN
               FLAG = -23
            ELSE IF (SIZE(data%X) < NEQ) THEN
               FLAG = -23
            END IF
            IF (FLAG == -23) THEN
               data%ERROR = -23
               IF (LERR) THEN
                  WRITE (ERR,FMT=9000) data%JOB,data%ERROR
                  IF (.NOT. allocated(data%X)) THEN
                     WRITE (ERR,FMT=9020) ' data%X'
                  ELSE IF (SIZE(data%X) < NEQ) THEN
                     WRITE (ERR,FMT=9030) 'data%X',SIZE(data%X),NEQ
                  END IF
               END IF
               GO TO 20
            END IF
         END IF
      END IF
! Broadcast error status from root
   20 CALL MPI_BCAST(data%ERROR,1,MPI_INTEGER,0,data%private%COMM,ERCODE)
! Terminate if error occurred
      IF (data%ERROR < 0) GO TO 890

      ST = 0
      IOS = 0
      ICOUNT_RANK = data%ICOUNT(RANK)

      DEALLOCATE (IFLAG,STAT=ST)
      ALLOCATE (IFLAG(NPROC),STAT=ST)
! Initialise the error flag (so that it is always well-defined)
      IF (ST == 0) IFLAG(1:NPROC) = 0

      DEALLOCATE (IBLOCK,NEQSUB,NVARSB,STAT=ST)
      ALLOCATE (IBLOCK(1:ICOUNT_RANK),NEQSUB(1:ICOUNT_RANK), &
                NVARSB(1:ICOUNT_RANK),STAT=ST)
      IF (ST /= 0) GO TO 30

! Set data from analyse phase
      DO LCNT = 1,ICOUNT_RANK
         IBLOCK(LCNT) = data%IBLOCK(LCNT)
         NEQSUB(LCNT) = data%private%NEQSUB(LCNT)
         JBLOCK = IBLOCK(LCNT)
         NVARSB(LCNT) = data%ENTRIES(JBLOCK)
      END DO

      LX = 0
      MMX = 0
      DO LCNT = 1,ICOUNT_RANK
         JBLOCK = IBLOCK(LCNT)
         M   = NEQSUB(LCNT)
         NC  = data%private%NCOLA(JBLOCK)
         LX  = MAX(LX,MAX(M,NC))
         LF = data%INFO_MA60(11,JBLOCK)
         MF = data%INFO_MA60(7,JBLOCK)
         NP = data%private%NP(JBLOCK)
         LX  = MAX(LX,NP+MF)
         MMX = MAX(MMX,M)
      END DO
      IF (RANK == 0) THEN
         DO IPROC = 1,NPROC-1
            SRC = IPROC
            ISTRT = data%IPLIST(IPROC)
            ISTOP = data%IPLIST(IPROC+1) - 1
            DO I = ISTRT,ISTOP
               JBLOCK = data%LIST(I)
               M = data%NEQSB(JBLOCK)
               NC  = data%private%NCOLA(JBLOCK)
               LX  = MAX(LX,MAX(M,NC))
               LF = data%INFO_MA60(11,JBLOCK)
               MF = data%INFO_MA60(7,JBLOCK)
               NP = data%private%NP(JBLOCK)
               LX  = MAX(LX,NP+MF)
            END DO
         END DO
         DEALLOCATE (XLOCAL,STAT=ST)
         ALLOCATE (XLOCAL(1:LX, 1:NBLOCK),STAT=ST )
      ELSE
         DEALLOCATE (XLOCAL,STAT=ST)
         ALLOCATE (XLOCAL(1:LX, 1:ICOUNT_RANK),STAT=ST )
      END IF
      IF (ST /= 0) GO TO 30

      DEALLOCATE (B1,WLOCAL,STAT=ST)
!!! CURRENT VERSION OF MP48_MA60C/CD REQUIRES B TO HAVE DIM(M,N)
      ALLOCATE (B1(1:LX, 1:ICOUNT_RANK),WLOCAL(1:LX),STAT=ST )
      IF (ST /= 0) GO TO 30

      IF (RANK /= 0) THEN
         DEALLOCATE (data%B,STAT=ST)
         ALLOCATE (data%B(1:NEQ),STAT=ST)
         IF (ST /= 0) GO TO 30
      END IF

! Allocate array for interface solution
      NINTER = data%NINTER
      RINTER = data%RINTER
      DEALLOCATE (X1,STAT=ST)
      ALLOCATE (X1(1:NINTER),STAT=ST )
      IF (ST /= 0) GO TO 30
! If the interface problem had no nonzero entries (possible in
! singular case) ensure X1 is defined
      IF (data%NE_INTER == 0) X1(1:NINTER) = ZERO

      DEALLOCATE (data%TIMES,STAT=ST)
      ALLOCATE (data%TIMES(1:NBLOCK),STAT=ST)
      IF (ST /= 0) GO TO 30

! Find suitable unit number if factor data held in files
! Ensure STRM is defined
      FLAG = 0
      STRM = 1
      IF (ICONTROL(11) == 0) GO TO 30

      DO STRM = 8,99
         IF (STRM == DGN .OR. STRM == ERR .OR. STRM == WRN) &
             CYCLE
         INQUIRE (UNIT=STRM,IOSTAT=IOS,ERR=30,EXIST=EX, &
                  OPENED=OPEN)
         IF (EX .AND. .NOT.OPEN) GO TO 25
      END DO
! No unit found. Jump to error return.
      FLAG = -20
      GO TO 30
   25 CONTINUE

! Gather error status from allocations
   30 CALL MPI_BARRIER(data%private%COMM,ERCODE)
      IF (IOS /= 0) FLAG = -14
      data%STAT = ST
      CALL MP48ND(data,NPROC,IFLAG,FLAG)
      data%IOSTAT = IOS
! Broadcast error status from root
      CALL MPI_BCAST(data%ERROR,1,MPI_INTEGER,0,data%private%COMM,ERCODE)
! Terminate if error occurred
      IF (data%ERROR < 0) GO TO 890

      CALL MPI_BCAST(data%B,NEQ,MPI_DOUBLE_PRECISION, &
                     0,data%private%COMM,ERCODE)

      IF (ICONTROL(11) /= 0) THEN
         L = 1
         L1 = 1
         DO LCNT = 1,ICOUNT_RANK
            JBLOCK = IBLOCK(LCNT)
            LFACT = data%INFO_MA60(4,JBLOCK)
            LIRNF = data%INFO_MA60(10,JBLOCK)
            L = MAX(L,LFACT)
            L1 = MAX(L1,LIRNF)
         END DO
         DEALLOCATE (data%FACT,data%IRNF,STAT=ST)
         ALLOCATE (data%FACT(1:L),data%IRNF(1:L1),STAT=ST)
         IF (ST /= 0) GO TO 120
      END IF

      FLAG = 0
      I2 = 0
      I2F = 0
      DO 100 LCNT = 1,ICOUNT_RANK
         T1 = MPI_WTIME()
         JBLOCK = IBLOCK(LCNT)
         M = NEQSUB(LCNT)
         NC = data%private%NCOLA(JBLOCK)
         DO K = 1,M
! J is the global index.
            J = data%private%GLOBAL(K,JBLOCK)
            B1(K,LCNT) = data%B(J)
         END DO
         IF (ICONTROL(11) /= 0) THEN
! Files used for factor data
            LFACT = data%INFO_MA60(4,JBLOCK)
            LIRNF = data%INFO_MA60(10,JBLOCK)
            I1 = 1
            I2 = LIRNF
            I1F = 1
            I2F = LFACT
            CLOSE (UNIT=STRM)
            OPEN (UNIT=STRM,IOSTAT=IOS,ERR=120, &
                  FILE=data%FILES(1,JBLOCK),FORM='UNFORMATTED')
            READ (STRM) data%FACT(1:LFACT)
            CLOSE (UNIT=STRM)
            OPEN (UNIT=STRM,IOSTAT=IOS,ERR=120, &
                  FILE=data%FILES(2,JBLOCK),FORM='UNFORMATTED')
            READ (STRM) data%IRNF(1:LIRNF)
            CLOSE (UNIT=STRM)
         ELSE
! Files not used ... data for each submatrix written
! consecutively to FACT and IRNF
            LFACT = data%private%LFACT(LCNT)
            LIRNF = data%private%LIRNF(LCNT)
            I1 = I2 + 1
            I2 = I1 + LIRNF - 1
            I1F = I2F + 1
            I2F = I1F + LFACT - 1
         END IF
!         IF (data%private%ICNTL_60(3,LCNT).GE.3 .AND.
!     &       data%private%ICNTL_60(2,LCNT) > 0)
!     &   WRITE (data%private%ICNTL_60(2,LCNT),'(/A,I3)')
!     & ' Calling MP48_MA60C/CD for forward elimination for submatrix',
!     &   JBLOCK

         LORU =  1
         IF (RANK == 0) THEN
! Use XLOCAL(1,JBLOCK)
            CALL MP48_MA60CD(M,NC,data%private%ICNTL_60(1:7,LCNT), &
                 data%private%IQ(1:NC,LCNT),data%private%NP(JBLOCK),LORU, &
                 LFACT,data%FACT(I1F:I2F),LIRNF,data%IRNF(I1:I2), &
                 data%private%IPTRL(1:NC,LCNT),data%private%IPTRU(1:NC,LCNT), &
                 B1(1,LCNT),XLOCAL(1,JBLOCK),WLOCAL,data%INFO_MA60(1:2,JBLOCK))
         ELSE
! Use XLOCAL(1,LCNT)
            CALL MP48_MA60CD(M,NC,data%private%ICNTL_60(1:7,LCNT), &
                 data%private%IQ(1:NC,LCNT),data%private%NP(JBLOCK),LORU, &
                 LFACT,data%FACT(I1F:I2F),LIRNF,data%IRNF(I1:I2), &
                 data%private%IPTRL(1:NC,LCNT),data%private%IPTRU(1:NC,LCNT), &
                 B1(1,LCNT),XLOCAL(1,LCNT),WLOCAL,data%INFO_MA60(1:2,JBLOCK))
         END IF

! No errors possible
         T2 = MPI_WTIME()
         data%TIMES(JBLOCK) = T2 - T1
  100 CONTINUE

  120 CALL MPI_BARRIER(data%private%COMM,ERCODE)
      DEALLOCATE (WLOCAL,STAT=ST)
      data%STAT = ST
      data%IOSTAT = IOS
      CALL MP48ND(data,NPROC,IFLAG,FLAG)
! Broadcast error status from root
      CALL MPI_BCAST(data%ERROR,1,MPI_INTEGER,0,data%private%COMM,ERCODE)
! Terminate if error occurred
      IF (data%ERROR < 0) GO TO 890

! At this point, XLOCAL holds partial solution vector.
! Send to host (no need if no interface variables
! or no nonzero entries in interface matrix).
      IF (data%NINTER /= 0 .AND. data%NE_INTER /= 0) THEN
       IF (RANK == 0) THEN
         DO IPROC = 1,NPROC-1
            SRC = IPROC
            ISTRT = data%IPLIST(IPROC)
            ISTOP = data%IPLIST(IPROC+1) - 1
            DO I = ISTRT,ISTOP
               JBLOCK = data%LIST(I)
!               M = data%NEQSB(JBLOCK)
!               CALL MPI_RECV(XLOCAL(1,JBLOCK),M,MPI_DOUBLE_PRECISION,
!     &                       SRC,1,data%private%COMM,STAT,ERCODE)
                LF = data%INFO_MA60(11,JBLOCK)
                MF = data%INFO_MA60(7,JBLOCK)
                NP = data%private%NP(JBLOCK)
                IF (MF > LF) &
                  CALL MPI_RECV(XLOCAL(NP+LF+1,JBLOCK),MF-LF, &
                                MPI_DOUBLE_PRECISION, &
                                SRC,1,data%private%COMM,STAT,ERCODE)
            END DO
         END DO

       ELSE
         DO LCNT = 1,ICOUNT_RANK
            JBLOCK = IBLOCK(LCNT)
! Send data to the root
!            M = NEQSUB(LCNT)
!            CALL MPI_SEND(XLOCAL(1,LCNT),M,MPI_DOUBLE_PRECISION,
!     &                    0,1,data%private%COMM,ERCODE)
             LF = data%INFO_MA60(11,JBLOCK)
             MF = data%INFO_MA60(7,JBLOCK)
             NP = data%private%NP(JBLOCK)
             IF (MF > LF) &
               CALL MPI_SEND(XLOCAL(NP+LF+1,LCNT),MF-LF, &
                          MPI_DOUBLE_PRECISION, &
                          0,1,data%private%COMM,ERCODE)
         END DO
       END IF
      END IF

! For interface work only on root
      IF (RANK == 0) THEN
         T3 = MPI_WTIME()
! Jump if no interface variables
! Also jump if no nonzero entries in interface matrix
         IF (data%NINTER == 0 .OR. data%NE_INTER == 0) GO TO 215

! Copy required part of right-hand sides B into an array B1
! that has rows 1 to NINTER
         DEALLOCATE (B1,IW,W,STAT=ST)
         ALLOCATE (B1(1:RINTER,1:1), &
                   IW(1:RINTER),W(1:3*(RINTER+NINTER)),STAT=ST )
         IF (ST /= 0) GO TO 220

         LA = data%private%LA48
! If files used, read data back in
         IF (ICONTROL(11) /= 0) THEN
            DEALLOCATE (data%private%A,data%private%IRN,STAT=ST)
            ALLOCATE (data%private%A(1:LA), &
                      data%private%IRN(1:LA),STAT=ST)
            IF (ST /= 0) GO TO 220

            CLOSE (UNIT=STRM)
            OPEN (UNIT=STRM,IOSTAT=IOS,ERR=220, &
                  FILE=data%FILES(1,NBLOCK+1),FORM='UNFORMATTED')
            READ (STRM) data%private%A
            CLOSE (UNIT=STRM)
            OPEN (UNIT=STRM,IOSTAT=IOS,ERR=220, &
                  FILE=data%FILES(2,NBLOCK+1),FORM='UNFORMATTED')
            READ (STRM) data%private%IRN
            CLOSE (UNIT=STRM)
         END IF

! We need to ensure B1 is defined (could miss some
! components in following do loop if interface problem has
! some rows with all entries zero)
         B1(1:RINTER,1) = ZERO
         K = 1
         DO JBLOCK = 1,NBLOCK
            LF = data%INFO_MA60(11,JBLOCK)
            MF = data%INFO_MA60(7,JBLOCK)
            NP = data%private%NP(JBLOCK)
            M = data%NEQSB(JBLOCK)
            DO I = LF+1,MF
!             write (6,*) 'jblock,lf,mf,np,k,i',jblock,lf,mf,np,k,i
              B1(K,1) = XLOCAL(NP+I,JBLOCK)
!             write (6,*) b1(k,1)
              K = K + 1
            END DO
         END DO

! Assume no iterative refinement
         JOB48 = 1
         LKEEP48 =  data%private%LKEEP48

         CALL MA48CD(RINTER,NINTER,.FALSE.,JOB48,LA,data%private%A(1:LA), &
              data%private%IRN(1:LA), &
                     data%private%KEEP48(1:LKEEP48),data%private%CNTL_48, &
                     data%private%ICNTL_48,B1,X1,ERR48,W,IW, &
                     data%INFO_MA48)

         DEALLOCATE (IW,W,STAT=ST)

  215    IF (data%ICNTL(13) == 0) THEN
            DEALLOCATE (data%X,STAT=ST)
            ALLOCATE (data%X(1:NEQ),STAT=ST )
            IF (ST /= 0) GO TO 220
         END IF

         DO K = 1,NINTER
            J = data%INTER_COL(K)
            data%X(J) = X1(K)
!           write (6,*) 'k,j,x(j)',k,j,data%x(j)
         END DO

         T4 = MPI_WTIME()
         data%TIME_MA48S = T4 - T3
      END IF

! Broadcast error status from root
  220 CALL MPI_BARRIER(data%private%COMM,ERCODE)
      IF (RANK == 0 .AND. ST /= 0) CALL MP48LD(data,ST)
      IF (RANK == 0 .AND. IOS /= 0) THEN
         data%ERROR = -17
         IF (LERR) THEN
            WRITE (ERR,FMT=9000) data%JOB,data%ERROR
            WRITE (ERR,FMT=9140)
         END IF
      END IF
      CALL MPI_BCAST(data%ERROR,1,MPI_INTEGER,0,data%private%COMM,ERCODE)
! Terminate if error occurred
      IF (data%ERROR < 0) GO TO 890
      CALL MPI_BCAST(X1,NINTER,MPI_DOUBLE_PRECISION,0,data%private%COMM,ERCODE)

! For non-interface variables, XLOCAL holds partial solution vector.
! For interface variables, X1 holds solution.
! For back sub., we need partial solution vector to be in B1
!!! current version of MP48_MA60C/CD requires B to have dim. maxn(m,n)
      DEALLOCATE (B1,STAT=ST)
      ALLOCATE (B1(1:LX, 1:ICOUNT_RANK),STAT=ST )
      IF (ST /= 0) GO TO 260

      IF (ICONTROL(11) /= 0) THEN
         L = 1
         L1 = 1
         DO LCNT = 1,ICOUNT_RANK
            JBLOCK = IBLOCK(LCNT)
            LFACT = data%INFO_MA60(4,JBLOCK)
            LIRNF = data%INFO_MA60(10,JBLOCK)
            L = MAX(L,LFACT)
            L1 = MAX(L1,LIRNF)
         END DO
         DEALLOCATE (data%FACT,data%IRNF,STAT=ST)
         ALLOCATE (data%FACT(1:L),data%IRNF(1:L1),STAT=ST)
         IF (ST /= 0) GO TO 260
      END IF

      FLAG = 0
      I2 = 0
      I2F = 0
      DO 250 LCNT = 1,ICOUNT_RANK

         T1 = MPI_WTIME()
         JBLOCK = IBLOCK(LCNT)
         M = NEQSUB(LCNT)
         NC = data%private%NCOLA(JBLOCK)
         NP = data%private%NP(JBLOCK)
         LF = data%INFO_MA60(11,JBLOCK)
         NF = data%private%NF(JBLOCK)

         DEALLOCATE (W,STAT=ST)
         ALLOCATE (W(1:MAX(M,NC)),STAT=ST)
         IF (ST /= 0) GO TO 260

         IF (RANK == 0) THEN
            B1(1:NP+LF,LCNT) = XLOCAL(1:NP+LF,JBLOCK)
         ELSE
            B1(1:NP+LF,LCNT) = XLOCAL(1:NP+LF,LCNT)
         END IF

         DO L = 1,NF-LF
            K = L + NP + LF
            L1 = data%private%JFVAR(L,JBLOCK)
!           write (6,'(a,4i5,g11.4)')
!     &    'jblock,l,k,l1,x1(l1)',jblock,l,k,l1,x1(l1)
            B1(K,LCNT) = X1(L1)
         END DO
! Make sure B1(1:max(M,NC)) defined ... needed by MP48_MA60C/CD
! (at least in current version)
         B1(NF+NP+1:MAX(M,NC),LCNT) = ZERO

         IF (ICONTROL(11) /= 0) THEN
! Files used for factor data so read data back in
            LFACT = data%INFO_MA60(4,JBLOCK)
            LIRNF = data%INFO_MA60(10,JBLOCK)
            I1 = 1
            I2 = LIRNF
            I1F = 1
            I2F = LFACT
            CLOSE (UNIT=STRM)
            OPEN (UNIT=STRM,IOSTAT=IOS,ERR=260, &
                  FILE=data%FILES(1,JBLOCK),FORM='UNFORMATTED')
            READ (STRM) data%FACT(1:LFACT)
            CLOSE (UNIT=STRM)
            OPEN (UNIT=STRM,IOSTAT=IOS,ERR=260, &
                  FILE=data%FILES(2,JBLOCK),FORM='UNFORMATTED')
            READ (STRM) data%IRNF(1:LIRNF)
            CLOSE (UNIT=STRM)
         ELSE
! Files not used ... data for each submatrix written
! consecutively to FACT and IRNF
            LFACT = data%private%LFACT(LCNT)
            LIRNF = data%private%LIRNF(LCNT)
            I1 = I2 + 1
            I2 = I1 + LIRNF - 1
            I1F = I2F + 1
            I2F = I1F + LFACT - 1
         END IF

         LORU = 0
!         IF (data%private%ICNTL_60(3,LCNT).GE.3 .AND.
!     &       data%private%ICNTL_60(2,LCNT) > 0)
!     &   WRITE (data%private%ICNTL_60(2,LCNT),'(/A,I3)')
!     & ' Calling MP48_MA60CD for backward substitution for submatrix',
!     &   JBLOCK
         IF (RANK == 0) THEN
! Use XLOCAL(1,JBLOCK)
            CALL MP48_MA60CD(M,NC,data%private%ICNTL_60(1:7,LCNT), &
                 data%private%IQ(1:NC,LCNT),data%private%NP(JBLOCK),LORU, &
                 LFACT,data%FACT(I1F:I2F),LIRNF,data%IRNF(I1:I2), &
                 data%private%IPTRL(1:NC,LCNT),data%private%IPTRU(1:NC,LCNT), &
                 B1(1,LCNT),XLOCAL(1,JBLOCK),W, &
                 data%INFO_MA60(1:2,JBLOCK))
         ELSE
! Use XLOCAL(1,LCNT)
            CALL MP48_MA60CD(M,NC,data%private%ICNTL_60(1:7,LCNT), &
                 data%private%IQ(1:NC,LCNT),data%private%NP(JBLOCK),LORU, &
                 LFACT,data%FACT(I1F:I2F),LIRNF,data%IRNF(I1:I2), &
                 data%private%IPTRL(1:NC,LCNT),data%private%IPTRU(1:NC,LCNT), &
                             B1(1,LCNT),XLOCAL(1,LCNT),W, &
                             data%INFO_MA60(1:2,JBLOCK))
         END IF
         DEALLOCATE (W,STAT=ST)
         T2 = MPI_WTIME()
         data%TIMES(JBLOCK) = data%TIMES(JBLOCK) + T2 - T1

  250 CONTINUE

! Synchronize
  260 CALL MPI_BARRIER(data%private%COMM,ERCODE)
      data%STAT = ST
      data%IOSTAT = IOS
      CALL MP48ND(data,NPROC,IFLAG,FLAG)
 ! Broadcast error status from root
      CALL MPI_BCAST(data%ERROR,1,MPI_INTEGER,0,data%private%COMM,ERCODE)
! Terminate if error occurred
      IF (data%ERROR < 0) GO TO 890

! At this point, XLOCAL holds solution vector.
! Send to host.
      IF (RANK == 0) THEN
        DO IPROC = 1,NPROC-1
          SRC = IPROC
          ISTRT = data%IPLIST(IPROC)
          ISTOP = data%IPLIST(IPROC+1) - 1
          DO I = ISTRT,ISTOP
            JBLOCK = data%LIST(I)
            NC = data%private%NCOLA(JBLOCK)
            CALL MPI_RECV(XLOCAL(1,JBLOCK),NC,MPI_DOUBLE_PRECISION, &
                          SRC,1,data%private%COMM,STAT,ERCODE)
            CALL MPI_RECV(data%TIMES(JBLOCK),1,MPI_DOUBLE_PRECISION, &
                          SRC,2,data%private%COMM,STAT,ERCODE)
          END DO
        END DO

      ELSE
        DO LCNT = 1,ICOUNT_RANK
          JBLOCK = IBLOCK(LCNT)
! Send data to the root
          NC = data%private%NCOLA(JBLOCK)
          CALL MPI_SEND(XLOCAL(1,LCNT),NC,MPI_DOUBLE_PRECISION, &
                        0,1,data%private%COMM,ERCODE)
          CALL MPI_SEND(data%TIMES(JBLOCK),1,MPI_DOUBLE_PRECISION, &
                        0,2,data%private%COMM,ERCODE)
        END DO
      END IF

      IF (RANK == 0) THEN
         IF (LDGN.GE.3) WRITE (DGN,FMT='(A,F9.3)') &
      ' The solve wall clock time (secs) for interface problem is:', &
        data%TIME_MA48S
         IF (LDGN.GE.3) WRITE (DGN,FMT='(A/8(F9.3))') &
      ' The solve wall clock time (secs) for each submatrix is:', &
         (data%TIMES(JBLOCK),JBLOCK=1,NBLOCK)
! Copy XLOCAL into appropriate part of X
         DO JBLOCK = 1,NBLOCK
            LF = data%INFO_MA60(11,JBLOCK)
            NP = data%private%NP(JBLOCK)
            NF = data%private%NF(JBLOCK)
            NC = data%private%NCOLA(JBLOCK)
            MF = data%INFO_MA60(7,JBLOCK)
            IF (NF-LF == data%NGUARD(JBLOCK) .OR. MF == 0) THEN
              DO K = 1,NC
!              DO K = 1,NP+LF
                J = data%private%CGLOB(K,JBLOCK)
                data%X(J) = XLOCAL(K,JBLOCK)
!               write (6,'(a,3i5,e16.8)')
!     &        'jblock,k,j,x(j)',jblock,k,j,data%x(j)
              END DO
            ELSE
              DO K = 1,NC
!              DO K = 1,NP+LF
                L = data%private%IQM_COPY(K,JBLOCK)
                J = data%private%CGLOB(L,JBLOCK)
                data%X(J) = XLOCAL(L,JBLOCK)
!               write (6,'(a,4i5,e16.8)')
!     &        'jblock,k,j,x(j)',jblock,k,l,j,data%x(j)
              END DO
            END IF
         END DO
      END IF

  890 CALL MPI_BARRIER(data%private%COMM,ERCODE)
      CALL MPI_BCAST(data%ERROR,1,MPI_INTEGER,0,data%private%COMM,ERCODE)
      data%private%ERROR = data%ERROR

      DEALLOCATE (B1,STAT=ST)
      DEALLOCATE (X1,STAT=ST)
      DEALLOCATE (W,STAT=ST)
      DEALLOCATE (XLOCAL,STAT=ST)
      DEALLOCATE (WLOCAL,STAT=ST)

      DEALLOCATE (IBLOCK,STAT=ST)
      DEALLOCATE (IFLAG,STAT=ST)
      DEALLOCATE (IW,STAT=ST)
      DEALLOCATE (NEQSUB,STAT=ST)
      DEALLOCATE (NVARSB,STAT=ST)

 9000 FORMAT (/' Error return with data%JOB = ',I2,'. ', &
                'data%ERROR = ',I4)
 9020 FORMAT (' ',A,' has not been allocated')
 9030 FORMAT (' ',A,' is of size ',I10,' Increase to at least ',I10)
 9140 FORMAT (' Failure in OPEN statement.')
      END SUBROUTINE MP48DD

!*************************************************

      SUBROUTINE MP48ED(M,N,KEEP,LKEEP,INFO)
! This subroutine computes statistics for MA48
!     MP48E/ED should be called after MA48BD/BD ..
! Changed on 11/03/02 to avoid bug when matrix was triangular
! M,N,KEEP .. the same as for that call ... LKEEP the dimension of KEEP.
! INFO(1)  Number of original entries excluding duplicates and
!          out-of-range
! INFO(2)  Number of entries in off-diagonal blocks of block triangular
!          form
! INFO(3)  Number of entries in triangular blocks on diagonal
! INFO(4)  Number of entries in L/U decomp of non-triangular blocks on
!          diagonal
! INFO(5)  Set to INFO(2)+INFO(3)+INFO(4) on exit.
!          This is the real storage required for factors.
!          The same storage is required for integers.
! INFO(3) and INFO(4) are separate because they are stored separately
!          .. the INFO(3) entries are also used to calculate residuals.

      USE HSL_MP48_DATA_DOUBLE
      IMPLICIT NONE
      INTEGER M,N,LKEEP
      INTEGER KEEP(LKEEP),INFO(5)
! Declare internal variables
      INTEGER IPTRL,IPTRU,IPTRD,IPTRO,NBLOCK,MBLOCK,KBLOCK
      INTEGER NB,KB,JB,J1,J2,NC,NR,LC
      INTRINSIC MAX

! allocatables for KEEP array
      IPTRL = M + N
      IPTRU = IPTRL + N
      IPTRD = IPTRU + N
      IPTRO = IPTRD + N + 1
      NBLOCK = IPTRO + N - 1
      MBLOCK = NBLOCK + 1
      KBLOCK = MBLOCK + 1

      NB = KEEP(KBLOCK+3)

! Number of entries in original matrix with dups and o-o-r omitted
      INFO(1) = KEEP(IPTRO+N+1) - 1

! Number entries in off-diagonal blocks
      INFO(2) = KEEP(IPTRO+N+1) - KEEP(IPTRO+1)

! Number of entries in triangular diagonal blocks
      INFO(3) = 0
      J2 = 0
      LC = 0
      DO 10 JB = 1,NB
         NC = KEEP(NBLOCK+3*JB)
         J1 = J2 + 1
         J2 = J1 + NC - 1
         IF (KEEP(MBLOCK+3*JB) < 0) THEN
            INFO(3) = INFO(3) + KEEP(IPTRD+J2+1) - KEEP(IPTRD+J1)
         ELSE
            KB = JB
            NR = NC
            LC = J2
         END IF
   10 CONTINUE

! Number of entries in factors of diagonal blocks
      INFO(4) = 0
      IF (LC /= 0) THEN
! At least one non-triangular block on diagonal
! Set NC to number of columns in last non-triangular block
        NC = NR
        IF (NB == 1) THEN
          NR = M
        ELSE
          INFO(4) = KEEP(KBLOCK+3*KB) - 1
        END IF
        INFO(4) = INFO(4) + KEEP(IPTRL+LC) + &
                 MAX(((NC-KEEP(MBLOCK+3*KB))+(NR-KEEP(MBLOCK+3*KB))), &
                      ((NC-KEEP(MBLOCK+3*KB))*(NR-KEEP(MBLOCK+3*KB))))
      END IF
      INFO(5) = INFO(2) + INFO(3) + INFO(4)
      END SUBROUTINE MP48ED

!********************************************************
      SUBROUTINE MP48LD(data,ST)

      USE HSL_MP48_DATA_DOUBLE
      IMPLICIT NONE
      TYPE (MP48_DATA) data
      INTEGER ST
      INTEGER ERR
      LOGICAL LERR

      data%ERROR = -11
      data%STAT = ST
      LERR = .FALSE.
      IF (data%ICNTL(1).GE.0 .AND. data%ICNTL(4).GE.1) LERR = .TRUE.
      ERR = data%ICNTL(1)
      IF (LERR) THEN
         WRITE (ERR,FMT=9000) data%JOB,data%ERROR
         WRITE (ERR,FMT=9050) ST
      END IF
 9000 FORMAT (/' Error return with data%JOB = ',I2,'.', &
               ' data%ERROR = ',I4)
 9050 FORMAT (' Failure in ALLOCATE statement. STAT = ',I5)

      END SUBROUTINE MP48LD

!********************************************************
      SUBROUTINE MP48MD(data,NPROC,IFLAG,ST)

      USE HSL_MP48_DATA_DOUBLE
      IMPLICIT NONE
      TYPE (MP48_DATA) data
      INTEGER NPROC,ST
      INTEGER IFLAG(NPROC)

      INTEGER FLAG,ERCODE,IPROC
      EXTERNAL MPI_GATHER

      FLAG = 0
      CALL MPI_GATHER(ST,1,MPI_INTEGER,IFLAG,1,MPI_INTEGER,0, &
                      data%private%COMM,ERCODE)
! Check for error on root
      IF (data%RANK == 0) THEN
         DO IPROC = 1,data%NPROC
            IF (IFLAG(IPROC) /= 0) THEN
               ST = IFLAG(IPROC)
               CALL MP48LD(data,ST)
               RETURN
            END IF
         END DO
      END IF
      data%STAT = ST

      END SUBROUTINE MP48MD

!********************************************************
      SUBROUTINE MP48ND(data,NPROC,IFLAG,FLAG)

      USE HSL_MP48_DATA_DOUBLE
      IMPLICIT NONE
      TYPE (MP48_DATA) data
      INTEGER NPROC,FLAG
      INTEGER IFLAG(NPROC)

      INTEGER ERCODE,ERR,IPROC
      LOGICAL LERR
      EXTERNAL MPI_GATHER

      ERR = data%ICNTL(1)
      LERR = .FALSE.
      IF (data%ICNTL(1).GE.0 .AND. data%ICNTL(4).GE.1) LERR = .TRUE.
      IF (data%STAT /= 0) FLAG = -11
      IF (data%IOSTAT /= 0) FLAG = -17

      CALL MPI_GATHER(FLAG,1,MPI_INTEGER,IFLAG,1,MPI_INTEGER,0, &
                      data%private%COMM,ERCODE)
! Check for error on root
      IF (data%RANK == 0) THEN
         DO IPROC = 1,data%NPROC
            IF (IFLAG(IPROC).GE.0) CYCLE
            IF (IFLAG(IPROC) == -17) THEN
               data%ERROR = -17
               FLAG = -17
               IF (LERR) THEN
                  WRITE (ERR,FMT=9000) data%JOB,data%ERROR
                  WRITE (ERR,FMT=9140)
               END IF
               RETURN
            ELSE IF (IFLAG(IPROC) == -14) THEN
               data%ERROR = -14
               FLAG = -14
               IF (LERR) THEN
                  WRITE (ERR,FMT=9000) data%JOB,data%ERROR
                  WRITE (ERR,FMT=9070) IFLAG(IPROC)
                END IF
                RETURN
            ELSE IF (IFLAG(IPROC) == -11) THEN
               data%ERROR = -11
               FLAG = -11
               IF (LERR) THEN
                  WRITE (ERR,FMT=9000) data%JOB,data%ERROR
                  WRITE (ERR,FMT=9050)
                END IF
                RETURN
            ELSE IF (IFLAG(IPROC) == -20) THEN
               data%ERROR = -20
               FLAG = -20
               IF (LERR) THEN
                  WRITE (ERR,FMT=9000) data%JOB,data%ERROR
                  WRITE (ERR,FMT=9180)
               END IF
               RETURN
            ELSE IF (IFLAG(IPROC) == -8) THEN
               data%ERROR = -15
               FLAG = -15
               IF (LERR) THEN
                  WRITE (ERR,FMT=9000) data%JOB,data%ERROR
                  WRITE (ERR,FMT=9035)
               END IF
               RETURN
            ELSE IF (IFLAG(IPROC) == -101) THEN
               data%ERROR = -101
               FLAG = -101
               IF (LERR) THEN
                  WRITE (ERR,FMT=9000) data%JOB,data%ERROR
                  WRITE (ERR,FMT=9030)
               END IF
               RETURN
            ELSE IF (IFLAG(IPROC) == -102) THEN
               data%ERROR = -102
               FLAG = -102
               IF (LERR) THEN
                  WRITE (ERR,FMT=9000) data%JOB,data%ERROR
                  WRITE (ERR,FMT=9040)
               END IF
               RETURN
            ELSE IF (IFLAG(IPROC) == -103) THEN
               data%ERROR = -103
               FLAG = -103
               IF (LERR) THEN
                  WRITE (ERR,FMT=9000) data%JOB,data%ERROR
                  WRITE (ERR,FMT=9045)
               END IF
               RETURN
            ELSE IF (IFLAG(IPROC) < 0) THEN
                IF (LERR) WRITE (ERR,FMT=9160)
                FLAG = -101
                data%ERROR = -101
                RETURN
            END IF
         END DO
      END IF
 9000 FORMAT (/' Error return with data%JOB = ',I2,'.', &
               ' data%ERROR = ',I4)
 9030 FORMAT (' Unexpected error return from MC46A/AD')
 9035 FORMAT (' Small pivot found during submatrix factorization')
 9040 FORMAT (' Unexpected error return from submatrix analyse')
 9045 FORMAT (' Unexpected error return from submatrix factorize')
 9050 FORMAT (' Failure in ALLOCATE statement. ')
 9070 FORMAT (' Failure in INQUIRE statement. IOSTAT = ',I5)
 9140 FORMAT (' Failure in OPEN statement.')
 9160 FORMAT (' Unexpected error return')
 9180 FORMAT (' Failed to find a unit to which a file', &
              ' may be connected.')

      END SUBROUTINE MP48ND
!********************************************************

!!! The following routine named MP48_MA60 are a version
!!! of MA50 developed for use by MP48.
!!! Note that not all the lines of the code will
!!! be executed since most of the error returns
!!! should not be possible. Also we do not offer the
!!! user the option of printing messages etc.
      SUBROUTINE MP48_MA60AD(M,N,NE,LA,A,IRN,JCN,IQ,CNTL,ICNTL,IP,NP, &
                             JFIRST, LENR,LASTR,NEXTR,IW,IFIRST, &
                             LENC,LASTC,NEXTC,INFO,RINFO)

! MP48_MA60A/AD chooses a pivot sequence using a Markowitz criterion
!     with threshold pivoting.

      USE HSL_MP48_DATA_DOUBLE
      IMPLICIT NONE
      REAL(WP), PARAMETER :: ZERO = 0.0_WP
      REAL(WP), PARAMETER :: ONE = 1.0_WP

      INTEGER M,N,NE,LA,NP
      REAL(WP) RINFO
      REAL(WP) A(LA),CNTL(4)
      INTEGER IRN(LA),JCN(LA),IQ(N)
      INTEGER ICNTL(10),IP(M),JFIRST(M),LENR(M),LASTR(M),NEXTR(M), &
              IW(M),IFIRST(N),LENC(N),LASTC(N),NEXTC(N),INFO(20)

! M is an integer variable that must be set to the number of rows.
!      It is not altered by the subroutine.
! N is an integer variable that must be set to the number of columns.
!      It is not altered by the subroutine.
! NE is an integer variable that must be set to the number of entries
!      in the input matrix. It is not altered by the subroutine.
! LA is an integer variable that must be set to the size of A, IRN, and
!      JCN. It is not altered by the subroutine.
! A is an array that holds the input matrix on entry and is used as
!      workspace.
! IRN  is an integer array.  Entries 1 to NE must be set to the
!      row indices of the corresponding entries in A.  IRN is used
!      as workspace and holds the row indices of the reduced matrix.
! JCN  is an integer array that need not be set by the user. It is
!      used to hold the column indices of entries in the reduced
!      matrix.
! IQ is an integer array of length N. On entry, it holds allocatables
!      to column starts. During execution, IQ(j) holds the position of
!      the start of column j of the reduced matrix or -IQ(j) holds the
!      column index in the permuted matrix of column j. On exit, IQ(j)
!      holds the index of the column that is in position j of the
!      permuted matrix.
! CNTL must be set by the user as follows and is not altered.
!     CNTL(1)  Full matrix processing will be used if the density of
!       the reduced matrix is MIN(CNTL(1),1.0) or more.
!     CNTL(2) determines the balance between pivoting for sparsity and
!       for stability, values near zero emphasizing sparsity and values
!       near one emphasizing stability. Each pivot must have absolute
!       value at least CNTL(2) times the greatest absolute value in the
!       same column of the reduced matrix.
!     CNTL(3) If this is set to a positive value, any entry of the
!       reduced matrix whose modulus is less than CNTL(3) will be
!       dropped.
!     CNTL(4)  Any entry of the reduced matrix whose modulus is less
!       than or equal to CNTL(4) will be regarded as zero from the
!        point of view of rank.
! ICNTL must be set by the user as follows and is not altered.
!     ICNTL(1)  must be set to the stream number for error messages.
!       A value less than 1 suppresses output.
!     ICNTL(2) must be set to the stream number for diagnostic output.
!       A value less than 1 suppresses output.
!     ICNTL(3) must be set to control the amount of output:
!       0 None.
!       1 Error messages only.
!       2 Error and warning messages.
!       3 As 2, plus scalar parameters and a few entries of array
!         parameters on entry and exit.
!       4 As 3, plus all parameters on entry and exit.
!     ICNTL(4) If set to a positive value, the pivot search is limited
!       to ICNTL(4) columns (Zlatev strategy). This may result in
!       different fill-in and execution time. If ICNTL(4) is positive,
!       the workspace arrays LASTR and NEXTR are not referenced.
!     ICNTL(5) The block size to be used for full-matrix processing.
!     ICNTL(6) The last ICNTL(6) columns of A must be the last
!       ICNTL(6) columns of the permuted matrix. A value outside the
!       range 1 to N-1 is treated as zero.
!     ICNTL(7) If given the value 1, pivots are limited to
!       the main diagonal, which may lead to a premature switch to full
!       processing if no suitable diagonal entries are available.
!       If given the value 2, IFIRST must be set so that IFIRST(i) is
!       the column in position i of the permuted matrix and IP must
!       be set so that IP(i) < IP(j) if row i is recommended to
!       precede row j in the pivot sequence.
!     ICNTL(8) Is the number of compresses allowed in the analysis
!       phase before an error return is invoked.  Default value is 10.
!     ICNTL(9) and ICNTL(10) are not currently used.
! IP is an integer array of length M that need not be set on entry
!      unless ICNTL(7)=2 (see ICNTL(7) for details of this case).
!      During execution, IP(i) holds the position of the start of row i
!      of the reduced matrix or -IP(i) holds the row index in the
!      permuted matrix of row i. Before exit, IP(i) is made positive.
! NP is an integer variable. It need not be set on entry. On exit,
!     it will be set to the number of columns processed in
!     packed storage.
! JFIRST is an integer workarray of length M. JFIRST(i) is the
!      first column of the reduced matrix to have i entries or is
!      zero if no column has i entries.
! LENR is an integer workarray of length M that is used to hold the
!      numbers of entries in the rows of the reduced matrix.
! LASTR is an integer workarray of length M, used only if ICNTL(4) = 0.
!      For rows in the reduced matrix, LASTR(i) indicates the previous
!      row to i with the same number of entries. LASTR(i) is zero if
!      no such row exists.
! NEXTR is an integer workarray of length M, used only if ICNTL(4) = 0
!      or ICNTL(7)=2. If ICNTL(4)=0, for rows in the reduced matrix,
!      NEXTR(i) indicates the next row to i with the same number of
!      entries; and if row i is the last in the chain, NEXTR is
!      equal to zero. If ICNTL(7)=2, NEXTR is a copy of the value of
!      IP on entry.
! IW is an integer array of length M used as workspace and is used to
!     assist the detection of duplicate entries and the sparse SAXPY
!     operations. It is reset to zero each time round the main loop.
! IFIRST is an integer array of length N, used only if ICNTL(4) = 0
!      or ICNTL(7)=2. If ICNTL(4) = 0, it is a workarray; IFIRST(i)
!      points to the first row of the reduced matrix to have i entries
!      or is zero if no row has i entries. If ICNTL(7)=2, IFIRST
!      must be set on entry (see ICNTL(7) for details of this case).
! LENC is an integer workarray of length N that is used to hold
!       the numbers of entries in the columns of the reduced matrix.
! LASTC is an integer workarray of length N.  For columns in the reduced
!      matrix, LASTC(j) indicates the previous column to j with the same
!      number of entries.  If column j is the first in the chain,
!      LASTC(j) is equal to zero.
! NEXTC is an integer workarray of length N.  For columns in the reduced
!      matrix, NEXTC(j) indicates the next column to j with the same
!      number of entries.  If column j is the last column in the chain,
!      NEXTC(j) is zero.
! INFO need not be set on entry. On exit, it holds the following:
!    INFO(1):
!       0  Successful entry.
!      -1  M < 1 or N < 1.
!      -2  NE < 1.
!      -3  Insufficient space.
!      -4  Duplicated entries.
!      -5  Faulty column permutation in IFIRST when ICNTL(7)=2.
!      -6  ICNTL(4) not equal to 1 when ICNTL(7)=2.
!      -7  More than ICNTL(xxx) compresses required
!      &1  Rank deficient.
!      &2  Premature switch to full processing because of failure to
!          find a stable diagonal pivot (ICNTL(7)>=1 case only).
!      &3  Both of these warnings.
!    INFO(2) Number of compresses of the arrays.
!    INFO(3) Minimum LA recommended to analyse matrix.
!    INFO(4) Minimum LFACT required to factorize matrix.
!    INFO(5) Upper bound on the rank of the matrix.
!    INFO(6) Number of entries dropped from the data structure.
!    INFO(7) Number of rows processed in full storage.
!    INFO(10) Minimum LIRNF required to factorize matrix.
!    INFO(11) Step at which failure because of too many compresses.
! RINFO need not be set on entry. On exit, it holds the number of
!    floating-point operations needed for the factorization.

      INTEGER IDAMAX
      EXTERNAL IDAMAX
      EXTERNAL MP48_MA60DD
      INTRINSIC ABS,MAX,MIN

      REAL(WP) ALEN,AMULT,ANEW,ASW,AU,COST,CPIV
      INTEGER DISPC,DISPR,EYE,I,IDROP,IDUMMY,IEND,IFILL,IFIR,II,IJ, &
              IJPOS,IOP,IPIV,IPOS,ISRCH,I1,I2,J,JBEG,JEND,JJ,JLAST, &
              JMORE,JNEW,JPIV,JPOS,J1,J2,L,LC,LEN,LENPIV,LP,LR
      REAL(WP) MAXENT
      INTEGER MINC,MORD,MP,MSRCH,NC,NDROP,NEFACT,NEPR,NERED,NE1,NORD, &
              NR,NULLC,NULLI,NULLJ1,NULLJ2,NULLR,PIVBEG,PIVCOL,PIVEND, &
              PIVOT,FULLCD,ISW
      REAL(WP) PIVR,PIVRAT,U

! ALEN Real(LEN-1).
! AMULT Temporary variable used to store current multiplier.
! ANEW Temporary variable used to store value of fill-in.
! ASW Temporary variable used when swopping two real quantities.
! AU Temporary variable used in threshold test.
! COST Markowitz cost of current potential pivot.
! CPIV Markowitz cost of best pivot so far found.
! DISPC is the first free location in the column file.
! DISPR is the first free location in the row file.
! EYE Running relative position when processing pivot row.
! I Temporary variable holding row number. Also used as index in DO
!     loops used in initialization of arrays.
! IDROP Temporary variable used to accumulate number of entries dropped.
! IDUMMY DO index not referenced in the loop.
! IEND Position of end of pivot row.
! IFILL is the fill-in to the non-pivot column.
! IFIR Temporary variable holding first entry in chain.
! II Running position for current column.
! IJ Temporary variable holding row/column index.
! IJPOS Position of current pivot in A/IRN.
! IOP holds a running count of the number of rows with entries in both
!     the pivot and the non-pivot column.
! IPIV Row of the pivot.
! IPOS Temporary variable holding position in column file.
! ISRCH Temporary variable holding number of columns searched for pivot.
! I1 Position of the start of the current column.
! I2 Position of the end of the current column.
! J Temporary variable holding column number.
! JBEG Position of beginning of non-pivot column.
! JEND Position of end of non-pivot column.
! JJ Running position for current row.
! JLAST Last column acceptable as pivot.
! JMORE Temporary variable holding number of locations still needed
!     for fill-in in non-pivot column.
! JNEW Position of end of changed non-pivot column.
! JPIV Column of the pivot.
! JPOS Temporary variable holding position in row file.
! J1 Position of the start of the current row.
! J2 Position of the end of the current row.
! L Loop index.
! LC Temporary variable holding previous column in sequence.
! LEN Length of column or row.
! LENPIV Length of pivot column.
! LP Unit for error messages.
! LR Temporary variable holding previous row in sequence.
! MAXENT Temporary variable used to hold value of largest entry in
!    column.
! MINC Minimum number of entries of any row or column of the reduced
!     matrix, or in any column if ICNTL(4) > 0.
! MORD Number of rows ordered, excluding null rows.
! MP Unit for diagnostic messages.
! MSRCH Number of columns to be searched.
! NC Temporary variable holding next column in sequence.
! NDROP Number of entries dropped because of being in a column all of
!   whose entries are smaller than the pivot threshold.
! NEFACT Number of entries in factors.
! NEPR Number of entries in pivot row, excluding the pivot.
! NERED Number of entries in reduced matrix.
! NE1 Temporary variable used to hold number of entries in row/column
!     and to hold temporarily value of MINC.
! NORD Number of columns ordered.
! NR Temporary variable holding next row in sequence.
! NULLC Number of structurally zero columns found before any entries
!     dropped for being smaller than CNTL(3).
! NULLR Number of structurally zero rows found before any entries
!     dropped for being smaller than CNTL(3).
! NULLI Number of zero rows found.
! NULLJ1 Number of zero columns found in non-border columns.
! NULLJ2 Number of zero columns found in border columns.
! PIVBEG Position of beginning of pivot column.
! PIVCOL Temporary variable holding position in pivot column.
! PIVEND Position of end of pivot column.
! PIVOT Current step in Gaussian elimination.
! PIVR ratio of current pivot candidate to largest in its column.
! PIVRAT ratio of best pivot candidate to largest in its column.
! U Used to hold local copy of CNTL(2), changed if necessary so that it
!    is in range.

      LP = ICNTL(1)
      IF (ICNTL(3).LE.0) LP = 0
      MP = ICNTL(2)
      IF (ICNTL(3).LE.1) MP = 0
      INFO(1) = 0
      INFO(2) = 0
      INFO(3) = NE
      INFO(4) = NE+2
      INFO(5) = 0
      INFO(6) = 0
      INFO(7) = 0
      INFO(10) = NE+2
      RINFO = ZERO
      NP = 0

! Make some simple checks
      IF (M.LT.1 .OR. N.LT.1) GO TO 690
      IF (NE.LT.1) GO TO 700
      IF (LA.LT.NE) THEN
         INFO(3) = NE
         GO TO 710
      END IF

! Initial printing
      IF (MP.GT.0 .AND. ICNTL(3).GT.2) THEN
         WRITE (MP,'(/2(A,I6),A,I8,A,I8/A,1P,4E10.2/A,8I4)') &
           ' Entering Analyse with M =',M,' N =',N,' NE =',NE, &
           ' LA =',LA,' CNTL =',CNTL,' ICNTL =',(ICNTL(J),J=1,8)
         IF (N.EQ.1 .OR. ICNTL(3).GT.3) THEN
            DO 10 J = 1,N - 1
               IF (IQ(J).LT.IQ(J+1)) WRITE (MP, &
                   '(A,I5,(T13,3(1P,E12.4,I5)))') ' Column',J, &
                   (A(II),IRN(II),II=IQ(J),IQ(J+1)-1)
   10       CONTINUE
            IF (IQ(N).LE.NE) WRITE (MP,'(A,I5,(T13,3(1P,E12.4,I5)))') &
                ' Column',N, (A(II),IRN(II),II=IQ(N),NE)
         ELSE
            IF (IQ(1).LT.IQ(2)) WRITE (MP, &
                '(A,I5,(T13,3(1P,E12.4,I5)))') ' Column',1, &
                (A(II),IRN(II),II=IQ(1),IQ(2)-1)
         END IF
         IF (ICNTL(7).EQ.2) THEN
            WRITE (MP,'(A,(T10,10(I7)))') ' IP = ',IP
            WRITE (MP,'(A,(T10,10(I7)))') ' IFIRST = ',IFIRST
         END IF
      END IF

! Initialization of counts etc.
      FULLCD = 0
      MINC = 1
      NERED = NE
      U = MIN(CNTL(2),ONE)
      U = MAX(U,ZERO)
      MSRCH = ICNTL(4)
      IF (MSRCH.EQ.0) MSRCH = N
      JLAST = N - ICNTL(6)
      IF (JLAST.LT.1 .OR. JLAST.GT.N) JLAST = N
      NULLI = 0
      NULLJ1 = 0
      NULLJ2 = 0
      MORD = 0
      NORD = 0
      NDROP = 0
      NEFACT = 0
      DO 20 I = 1,N - 1
         LENC(I) = IQ(I+1) - IQ(I)
   20 CONTINUE
      LENC(N) = NE + 1 - IQ(N)

!     write (6,'(A/(20I8))') 'LENC',(LENC(I),I=1,N)

      IF (CNTL(3).GT.ZERO) THEN
! Drop small entries
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
!                 write (6,*) 'Dropping entry .. J,I,A(II)',J,I,A(II)
                  INFO(6) = INFO(6) + 1
               END IF
   30       CONTINUE
            LENC(J) = NERED + 1 - IQ(J)
   40    CONTINUE
      END IF

      IF (ICNTL(7).EQ.2) THEN
! Column order specified - copy the row ordering array
         DO 50 I = 1,M
            NEXTR(I) = IP(I)
   50    CONTINUE
! Check ICNTL(4)
         IF (ICNTL(4).NE.1) GO TO 740
      END IF

      DISPR = NERED + 1
      DISPC = NERED + 1

! Set up row oriented storage.
      DO 60 I = 1,M
         IW(I) = 0
         LENR(I) = 0
         JFIRST(I) = 0
   60 CONTINUE
! Calculate row counts.
      DO 70 II = 1,NERED
         I = IRN(II)
         LENR(I) = LENR(I) + 1
   70 CONTINUE
! Set up row allocatables so that IP(i) points to position after end
! of row i in row file.
      IP(1) = LENR(1) + 1
      DO 80 I = 2,M
         IP(I) = IP(I-1) + LENR(I)
   80 CONTINUE
! Generate row file.
      DO 100 J = 1,N
         I = IQ(J)
         DO 90 II = I,I + LENC(J) - 1
            I = IRN(II)
! Check for duplicate entry.
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

! Check for zero rows and (unless ICNTL(4) > 0), compute chains of rows
! with equal numbers of entries.
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
! Check for zero columns and compute chains of columns with equal
!  numbers of entries.
      DO 150 J = N,1,-1
         NE1 = LENC(J)
         IF (NE1.EQ.0) THEN
            IF (ICNTL(7).NE.2) THEN
               IF (J.LE.JLAST) THEN
!       write (6,*) 'Oops .. problem for now, Nord,j,jlast =',
!     &            NORD,J,JLAST
                  NULLJ1 = NULLJ1 + 1
                  IQ(J) = -(N+1)
                  IF (NORD+NULLJ1.EQ.JLAST) THEN
! We have ordered the first N - ICNTL(6) columns.
                     GO TO 640
                  END IF
               ELSE
!       write (6,*) 'Oops .. nullj incremented in loop 150',
!     &    NULLJ1,NULLJ2,J,JLAST
                  NULLJ2 = NULLJ2 + 1
                  IQ(J) = -(N+2)
               END IF
               LASTC(J) = 0
               NEXTC(J) = 0
            END IF
         ELSE
! Column is non-null
            IFIR = JFIRST(NE1)
            JFIRST(NE1) = J
            NEXTC(J) = IFIR
            LASTC(J) = 0
            IF (IFIR.GT.0) LASTC(IFIR) = J
         END IF
  150 CONTINUE
      IF (INFO(6).EQ.0) THEN
         NULLC = NORD + NULLJ1 + NULLJ2
         NULLR = NULLI
      END IF

! **********************************************
! ****    Start of main elimination loop    ****
! **********************************************
      DO 630 PIVOT = 1,N
! Check to see if reduced matrix should be considered as full.

         IF (NERED.GE. (MIN(CNTL(1),ONE)*(N-NORD))* &
             (M-MORD)) GO TO 640

! Jump if we have chosen as many pivots as are allowed.  Rest of matrix
! columns are in the border.  Action is same as switch to full
! code.
         IF (FULLCD.EQ.1) GO TO 640

         IF (ICNTL(7).EQ.2) THEN
! Column order specified - choose the pivot within the column
            IPIV = 0
            J = IFIRST(PIVOT)
            IF (J.LT.1 .OR. J.GT.N) GO TO 730
            IF (IQ(J).LT.0) GO TO 730
            LEN = LENC(J)
            IF (LEN.LE.0) GO TO 320
            ALEN = LEN - 1
            I1 = IQ(J)
            I2 = I1 + LEN - 1
! Find largest entry in column
            II = IDAMAX(LEN,A(I1),1)
            MAXENT = ABS(A(I1+II-1))
! Is every entry in the column below the pivot threshold?
!           if (maxent.le.cntl(4)) write (6,*) 'maxent',maxent
            IF (MAXENT.LE.CNTL(4)) GO TO 320
            AU = MAX(MAXENT*U,CNTL(4))
!            write (6,*) 'AU',AU
! Scan column for pivot
            DO 160 II = I1,I2
               IF (ABS(A(II)).LT.AU) GO TO 160
! Candidate satisfies threshold criterion.
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

! Find the least number of entries in a row or column (column only if
! the Zlatev strategy is in use)
         LEN = MAX(1,MINC)
         DO 170 MINC = LEN,M - MORD
            IF (JFIRST(MINC).NE.0) GO TO 180
            IF (ICNTL(4).LE.0) THEN
               IF (IFIRST(MINC).NE.0) GO TO 180
            END IF
  170    CONTINUE

! Find the next pivot or a column whose entries are all very small.
! CPIV is the Markowitz cost of the best pivot so far and PIVRAT is the
! ratio of its absolute value to that of the largest entry in its
! column.
  180    CPIV = M
         CPIV = CPIV*N
         PIVRAT = ZERO
! Examine columns/rows in order of ascending count.
         ISRCH = 0
         DO 300 LEN = MINC,M - MORD
            ALEN = LEN - 1
! Jump if Markowitz count cannot be bettered.
            IF (CPIV.LE.ALEN**2 .AND. ICNTL(4).LE.0) GO TO 310
            IJ = JFIRST(LEN)
! Scan columns with LEN entries.
            DO 220 IDUMMY = 1,N
! If no more columns with LEN entries, exit loop.
               IF (IJ.LE.0) GO TO 230
               J = IJ
!      write (6,*) 'Testing column ',J,' at stage ',PIVOT
               IJ = NEXTC(J)
               IF (J.GT.N - ICNTL(6)) GO TO 220
! Column J is now examined.
! First calculate multiplier threshold level.
               MAXENT = ZERO
               I1 = IQ(J)
               I2 = I1 + LEN - 1
               II = IDAMAX(LEN,A(I1),1)
               MAXENT = ABS(A(I1+II-1))
! Exit loop if every entry in the column is below the pivot threshold.
!              if (maxent.le.cntl(4)) write (6,*) 'maxent..',maxent
               IF (MAXENT.LE.CNTL(4)) GO TO 320
               AU = MAX(MAXENT*U,CNTL(4))
!              write (6,*) 'AU..',AU
! If diagonal pivoting requested, look for diagonal entry.
               IF (ICNTL(7).EQ.1) THEN
                  DO 190 II = I1,I2
                     IF (IRN(II).EQ.J) GO TO 200
  190             CONTINUE
                  GO TO 220
  200             I1 = II
                  I2 = II
               END IF
! Scan column for possible pivots
               DO 210 II = I1,I2
                  IF (ABS(A(II)).LT.AU) GO TO 210
! Candidate satisfies threshold criterion.
                  I = IRN(II)
                  COST = ALEN*(LENR(I)-1)
                  IF (COST.GT.CPIV) GO TO 210
                  PIVR = ABS(A(II))/MAXENT
                  IF (COST.EQ.CPIV) THEN
                     IF (PIVR.LE.PIVRAT) GO TO 210
                  END IF
! Best pivot so far is found.
                  CPIV = COST
                  IJPOS = II
                  IPIV = I
                  JPIV = J
                  IF (CPIV.LE.ALEN**2 .AND. ICNTL(4).LE.0) GO TO 330
                  PIVRAT = PIVR
  210          CONTINUE
! Increment number of columns searched.
               ISRCH = ISRCH + 1
! Jump if we have searched the number of columns stipulated and found a
! pivot.
               IF (ISRCH.GE.MSRCH) THEN
                  IF (PIVRAT.GT.ZERO) GO TO 330
               END IF
  220       CONTINUE

! Rows with LEN entries now examined.
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
! If diagonal pivoting requested, look for diagonal entry.
               IF (ICNTL(7).EQ.1) THEN
                  DO 240 JJ = J1,J2
                     IF (JCN(JJ).EQ.I) GO TO 250
  240             CONTINUE
                  GO TO 290
  250             J1 = JJ
                  J2 = JJ
               END IF
! Scan row I.
               DO 280 JJ = J1,J2
                  J = JCN(JJ)
                  IF (J.GT.N - ICNTL(6)) GO TO 280
                  COST = ALEN*(LENC(J)-1)
                  IF (COST.GE.CPIV) GO TO 280
! Pivot has best Markowitz count so far. Now check its suitability
!     on numerical grounds by examining other entries in its column.
                  I1 = IQ(J)
                  I2 = I1 + LENC(J) - 1
                  II = IDAMAX(LENC(J),A(I1),1)
                  MAXENT = ABS(A(I1+II-1))
                  DO 260 II = I1,I2 - 1
                     IF (IRN(II).EQ.I) GO TO 270
  260             CONTINUE
  270             JPOS = II
! Exit loop if every entry in the column is below the pivot threshold.
!               if (maxent.le.cntl(4)) write (6,*) 'maxent...',maxent
!               if (maxent.le.cntl(4)) write (6,*) 'u...',u
                  IF (MAXENT.LE.CNTL(4)) GO TO 320
                  IF (ABS(A(JPOS)).LT.MAXENT*U) GO TO 280
! Candidate satisfies threshold criterion.
                  CPIV = COST
                  IPIV = I
                  JPIV = J
                  IJPOS = JPOS
                  PIVRAT = ABS(A(JPOS))/MAXENT
                  IF (CPIV.LE.ALEN*(ALEN+1)) GO TO 330
  280          CONTINUE

  290       CONTINUE

  300    CONTINUE
  310    IF (PIVRAT.GT.ZERO) GO TO 330
! No pivot found. Switch to full matrix processing.
         INFO(1) = INFO(1) + 2
         IF (MP.GT.0) WRITE (MP,'(A/A)') &
             ' Warning message from Analyse: no suitable diagonal', &
             ' pivot found, so switched to full matrix processing.'
         GO TO 640

! Every entry in the column is below the pivot threshold.
  320    IPIV = 0
         JPIV = J

! The pivot has now been found in position (IPIV,JPIV) in location
! IJPOS in column file or all entries of column JPIV are very small
! (IPIV=0).
! Update row and column ordering arrays to correspond with removal
! of the active part of the matrix. Also update NEFACT.
  330    NEFACT = NEFACT + LENC(JPIV)
!        write (6,*) 'Pivot found in position ',IPIV,JPIV,
!     &            ' at stage ',PIVOT
         PIVBEG = IQ(JPIV)
         PIVEND = PIVBEG + LENC(JPIV) - 1

         IF (IPIV.GT.0) THEN
           NORD = NORD + 1
         ELSE
           NULLJ1 = NULLJ1 + 1
           IQ(J) = -(N+1)
         END IF
         IF (NORD+NULLJ1.EQ.JLAST) THEN
!        write (6,'(A,A/(20I8))')
!     &   'We have ordered all the columns we can',' NORD, NULLJ,
!     &    N, ICNTL(6), JLAST',NORD,NULLJ1,N,ICNTL(6),JLAST
! We have ordered the first N - ICNTL(6) columns.
! In this code we now switch to the full code
            FULLCD = 1
         END IF
         IF (ICNTL(4).LE.0) THEN
! Remove active rows from their row ordering chains.
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
! NEPR is number of entries in strictly U part of pivot row.
            NEPR = LENR(IPIV) - 1
            NEFACT = NEFACT + NEPR
            RINFO = RINFO + CPIV*2 + LENR(IPIV)
            J1 = IP(IPIV)
! Remove active columns from their column ordering chains.
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
! Move pivot to beginning of pivot column.
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
!           if (cntl(3).gt.zero) write (6,*) 'ndrop,ne1',ndrop,ne1
            IF (NE1.GT.0) THEN
! Remove column of small entries from its column ordering chain.
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

! Set up IW array so that IW(i) holds the relative position of row i
! entry from beginning of pivot column.
         DO 360 II = PIVBEG + 1,PIVEND
            I = IRN(II)
            IW(I) = II - PIVBEG
  360    CONTINUE
! LENPIV is length of strictly L part of pivot column.
         LENPIV = PIVEND - PIVBEG

! Remove pivot column (including pivot) from row oriented file.
         DO 390 II = PIVBEG,PIVEND
            I = IRN(II)
            LENR(I) = LENR(I) - 1
            J1 = IP(I)
! J2 is last position in old row.
            J2 = J1 + LENR(I)
            DO 370 JJ = J1,J2 - 1
               IF (JCN(JJ).EQ.JPIV) GO TO 380
  370       CONTINUE
  380       JCN(JJ) = JCN(J2)
            JCN(J2) = 0
  390    CONTINUE

! For each active column, add the appropriate multiple of the pivot
! column to it.
! We loop on the number of entries in the pivot row since the position
! of this row may change because of compresses.
         DO 600 EYE = 1,NEPR
            J = JCN(IP(IPIV)+EYE-1)
! Search column J for entry to be eliminated, calculate multiplier,
! and remove it from column file.
! IDROP is the number of nonzero entries dropped from column J
! because these fall beneath tolerance level.
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
!      write (6,*) 'Setting entry',JEND,' of IRN to zero'
            IRN(JEND) = 0
            JEND = JEND - 1
! Jump if pivot column is a singleton.
            IF (LENPIV.EQ.0) GO TO 600
! Now perform necessary operations on rest of non-pivot column J.
            IOP = 0
! Innermost loop.
            DO 420 II = JBEG,JEND
               I = IRN(II)
               IF (IW(I).GT.0) THEN
! Row i is involved in the pivot column.
                  IOP = IOP + 1
                  PIVCOL = IQ(JPIV) + IW(I)
! Flag IW(I) to show that the operation has been done.
                  IW(I) = -IW(I)
                  A(II) = A(II) + AMULT*A(PIVCOL)
               END IF
  420       CONTINUE

            IF (CNTL(3).GT.ZERO) THEN
!  Run through non-pivot column compressing column so that entries less
!  than CNTL(3) are not stored. All entries less than CNTL(3) are
!  also removed from the row structure.
               JNEW = JBEG
               DO 450 II = JBEG,JEND
                  IF (ABS(A(II)).GE.CNTL(3)) THEN
                     A(JNEW) = A(II)
                     IRN(JNEW) = IRN(II)
                     JNEW = JNEW + 1
                  ELSE
!  Remove non-zero entry from row structure.
                     I = IRN(II)
!                    write (6,*) 'Remove entry in row',I
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
!                write (6,*) 'Setting entry',II,' of IRN to zero in 460'
                 IRN(II) = 0
  460          CONTINUE
               IDROP = JEND + 1 - JNEW
               JEND = JNEW - 1
               LENC(J) = LENC(J) - IDROP
               NERED = NERED - IDROP
               INFO(6) = INFO(6) + IDROP
            END IF

! IFILL is fill-in left to do to non-pivot column J.
            IFILL = LENPIV - IOP
            NERED = NERED + IFILL
            INFO(3) = MAX(INFO(3),NEFACT+NERED+LENC(J))

! Treat no-fill case
            IF (IFILL.EQ.0) THEN
               DO 470 II = PIVBEG + 1,PIVEND
                  I = IRN(II)
                  IW(I) = -IW(I)
  470          CONTINUE
               GO TO 600
            END IF

!            write (6,*) 'DISPC,JEND,IFILL',DISPC,JEND,IFILL
!            write (6,*) 'IRN ...',(IRN(IPOS),IPOS=1,DISPC-1)
!            write (6,*) 'A   ...',(A(IPOS),IPOS=1,DISPC-1)

! See if there is room for fill-in at end of the column.
            DO 480 IPOS = JEND + 1,MIN(JEND+IFILL,DISPC-1)
               IF (IRN(IPOS).NE.0) GO TO 490
  480       CONTINUE
            IF (IPOS.EQ.JEND+IFILL+1) GO TO 540
            IF (JEND+IFILL+1.LE.LA+1) THEN
               DISPC = JEND + IFILL + 1
               GO TO 540
            END IF
            IPOS = LA
!           write (6,*) 'Setting DISPC to LA+1'
! JMORE more spaces for fill-in are required.
  490       JMORE = JEND + IFILL - IPOS + 1
! We now look in front of the column to see if there is space for
! the rest of the fill-in.
            DO 500 IPOS = JBEG - 1,MAX(JBEG-JMORE,1),-1
               IF (IRN(IPOS).NE.0) GO TO 510
  500       CONTINUE
            IPOS = IPOS + 1
            IF (IPOS.EQ.JBEG-JMORE) GO TO 520
! Column must be moved to the beginning of available storage.
  510       IF (DISPC+LENC(J)+IFILL.GT.LA+1) THEN
               IF (INFO(2).EQ.ICNTL(8)) THEN
                 INFO(1) = -7
                 INFO(3) = MAX(INFO(3),NEFACT+NERED)
                 INFO(11) = PIVOT
                 GO TO 750
               END IF
               INFO(2) = INFO(2) + 1
               CALL MP48_MA60DD(LA,A,IRN,IQ,N,DISPC,.TRUE.)
               JBEG = IQ(J)
               JEND = JBEG + LENC(J) - 1
               PIVBEG = IQ(JPIV)
               PIVEND = PIVBEG + LENC(JPIV) - 1
               IF (DISPC+LENC(J)+IFILL.GT.LA+1) GO TO 705
            END IF
            IPOS = DISPC
            DISPC = DISPC + LENC(J) + IFILL
! Move non-pivot column J.
  520       IQ(J) = IPOS
!           write (6,*) 'After 520, J,IPOS',J,IPOS
            DO 530 II = JBEG,JEND
               A(IPOS) = A(II)
               IRN(IPOS) = IRN(II)
               IPOS = IPOS + 1
!              write (6,*) 'Setting',II,' to zero'
               IRN(II) = 0
  530       CONTINUE
            JBEG = IQ(J)
            JEND = IPOS - 1
! Innermost fill-in loop which also resets IW.
! We know at this stage that there are IFILL positions free after JEND.
  540       IDROP = 0
            DO 580 II = PIVBEG + 1,PIVEND
               I = IRN(II)
               INFO(3) = MAX(INFO(3),NEFACT+NERED+LENR(I)+1)
               IF (IW(I).LT.0) THEN
                  IW(I) = -IW(I)
                  GO TO 580
               END IF
               ANEW = AMULT*A(II)
               IF (ABS(ANEW).LT.CNTL(3)) THEN
!                 write (6,*) 'IDROP,ANEW,AMULT',IDROP,ANEW,AMULT
                  IDROP = IDROP + 1
               ELSE
                  JEND = JEND + 1
                  A(JEND) = ANEW
                  IRN(JEND) = I

! Put new entry in row file.
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
! Copy row forward
                        IEND = IEND - 1
                        DO 545 JJ = IP(I),IEND
                           JCN(JJ-1) = JCN(JJ)
  545                   CONTINUE
                        IP(I) = IP(I) - 1
                        GO TO 560
                     END IF
                  END IF
                  IF (DISPR+LENR(I).GT.LA) THEN
! Compress.
                    IF (INFO(2).EQ.ICNTL(8)) THEN
                      INFO(1) = -7
                      INFO(3) = MAX(INFO(3),NEFACT+NERED)
                      INFO(11) = PIVOT
                      GO TO 750
                    END IF
                    INFO(2) = INFO(2) + 1
                    CALL MP48_MA60DD(LA,A,JCN,IP,M,DISPR,.FALSE.)
                    IF (DISPR+LENR(I).GT.LA) GO TO 705
                  END IF
! Copy row to first free position.
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
! End of adjustment to row file.
               END IF
  580       CONTINUE
            INFO(6) = INFO(6) + IDROP
            NERED = NERED - IDROP
            DO 590 II = 1,IDROP
!              write (6,*) 'Drop loop .. setting',JEND+II,' to zero'
               IRN(JEND+II) = 0
  590       CONTINUE
            LENC(J) = LENC(J) + IFILL - IDROP
! End of scan of pivot row.
  600    CONTINUE


! Remove pivot row from row oriented storage and update column
! ordering arrays.  Remember that pivot row no longer includes pivot.
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

! In this case matrix is singular and we will take other action
! By moving singular part to trailing columns
                  IQ(J) = -(N+1)
!                 write (6,*) 'recording null column ... j,jlast,nullj',
!     &           J,JLAST,NULLJ1,NULLJ2
                  NULLJ1 = NULLJ1 + 1
                  IF (NORD+NULLJ1.EQ.JLAST) THEN
! We have ordered the first N - ICNTL(6) columns.

!                  write (6,*)
!     &            'A second source of problems .. nord,nullj',
!     &            NORD,NULLJ1,NULLJ2
                     FULLCD = 1
                  END IF
               ELSE
!                  write (6,*)
!     &           'Another source of problems .. j,nullj,jlast',
!     &            J,NULLJ1,NULLJ2,JLAST
                  NULLJ2 = NULLJ2 + 1
                  IQ(J) = -(N+2)
               END IF
            END IF
  610    CONTINUE
         NERED = NERED - NEPR

! Restore IW and remove pivot column from column file.
! Record the row permutation in IP(IPIV) and the column
! permutation in IQ(JPIV), flagging them negative so that they
! are not confused with real allocatables in compress routine.
         IF (IPIV.NE.0) THEN
            LENR(IPIV) = 0
            IW(IPIV) = 0
!           write (6,*) 'Pivot removal',PIVBEG,' to zero'
            IRN(PIVBEG) = 0
            MORD = MORD + 1
            PIVBEG = PIVBEG + 1
            IP(IPIV) = -MORD
         END IF
         NERED = NERED - LENPIV - 1
         DO 620 II = PIVBEG,PIVEND
            I = IRN(II)
            IW(I) = 0
!           write (6,*) 'Pivott removal',II,' to zero'
            IRN(II) = 0
            NE1 = LENR(I)
            IF (NE1.EQ.0) THEN
               IF (INFO(6).EQ.0) NULLR = NULLR + 1
               IP(I) = -M + NULLI
               NULLI = NULLI + 1
            ELSE IF (ICNTL(4).LE.0) THEN
! Adjust row ordering arrays.
               IFIR = IFIRST(NE1)
               LASTR(I) = 0
               NEXTR(I) = IFIR
               IFIRST(NE1) = I
               IF (IFIR.NE.0) LASTR(IFIR) = I
               MINC = MIN(MINC,NE1)
            END IF
  620    CONTINUE
!        write (6,*) 'IPIV,JPIV,NORD',IPIV,JPIV,NORD
         IF (IPIV.GT.0) IQ(JPIV) = -NORD
  630 CONTINUE
! We may drop through this loop with NULLI nonzero.

! ********************************************
! ****    End of main elimination loop    ****
! ********************************************

! Complete the permutation vectors
  640 INFO(5) = MORD + MIN(M-MORD-NULLI,N-NORD-NULLJ1-NULLJ2)
      DO 650 L = 1,MIN(M-MORD,N-NORD)
         RINFO = RINFO + M - MORD - L + 1 + REAL(M-MORD-L)*(N-NORD-L)*2
  650 CONTINUE
      NP = NORD
      INFO(4) = MAX(INFO(4),2 + NEFACT + M*2 + (N-NORD)*(M-MORD))
      INFO(6) = INFO(6) + NDROP
      INFO(7) = M - MORD
      INFO(10) = MAX(INFO(10),2 + NEFACT + M*2 + (N-NORD+M-MORD))
      DO 660 L = 1,M
         IF (IP(L).LT.0) THEN
            IP(L) = -IP(L)
         ELSE
            MORD = MORD + 1
            IP(L) = MORD
         END IF
  660 CONTINUE
!     write (6,'(A/(20I8))') 'N, IQ and LASTC after 660',
!     &                       N,(IQ(I),I=1,N),(LASTC(I),I=1,N)
!     write (6,*) 'IQ',(IQ(L),L=1,N)
      DO 670 L = 1,N
         IF (IQ(L).LT.0) THEN
           IF (IQ(L).EQ.-(N+1)) THEN
             NORD = NORD + 1
             LASTC(L) = NORD
           ELSE
            LASTC(L) = -IQ(L)
           END IF
         END IF
  670 CONTINUE
!     write (6,*) 'Issue in loop 670'
!     write (6,'(A/(20I8))') 'L,NORD, NULLJ1,JLAST',L,NORD,NULLJ1,JLAST
      DO 675 L = 1,N
         IF (IQ(L).GT.0 .OR. IQ(L).EQ.-(N+2)) THEN
            NORD = NORD + 1
            LASTC(L) = NORD
         END IF
  675 CONTINUE
!     write (6,*) 'LASTC',(LASTC(L),L=1,N)

!     write (6,'(A/(20I8))') 'N, IQ and LASTC after 670', &
!                            N,(IQ(I),I=1,N),(LASTC(I),I=1,N)
! Store the inverse permutation
      DO 680 L = 1,N
         IQ(LASTC(L)) = L
  680 CONTINUE

! Test for rank deficiency
      IF (INFO(5).LT.MIN(M,N)) INFO(1) = INFO(1) + 1

      IF (MP.GT.0 .AND. ICNTL(3).GT.2) THEN
         WRITE (MP,'(A,I6,A,F12.1/A,7I8/A,7I8)') &
         ' Leaving Analyse with NP =', &
           NP,' RINFO =',RINFO, &
         ' INFO(1:7)  =',INFO(1:7), &
         ' INFO(10)   =',INFO(10)
         IF (ICNTL(3).GT.3) THEN
            WRITE (MP,'(A,(T6,10(I7)))') ' IP = ',IP
            WRITE (MP,'(A,(T6,10(I7)))') ' IQ = ',IQ
         END IF
      END IF

      GO TO 750

! Error conditions.
  690 INFO(1) = -1
      IF (LP.GT.0) WRITE (LP,'(/A/(2(A,I8)))') &
          ' **** Error return from Analyse ****',' M =',M,' N =',N
      GO TO 750
  700 INFO(1) = -2
      IF (LP.GT.0) WRITE (LP,'(/A/(A,I10))') &
          ' **** Error return from Analyse ****',' NE =',NE
      GO TO 750
  705 INFO(4) =  MAX(NE+2,NEFACT + NERED)
      INFO(10) =  MAX(NE+2,NEFACT + NERED)
      INFO(6) = INFO(6) + NDROP
  710 INFO(1) = -3
      IF (LP.GT.0) WRITE (LP,'(/A/A,I9,A,I9)') &
          ' **** Error return from Analyse ****', &
          ' LA  must be increased from',LA,' to at least',INFO(3)
      GO TO 750
  720 INFO(1) = -4
      IF (LP.GT.0) WRITE (LP,'(/A/(3(A,I9)))') &
          ' **** Error return from Analyse ****',' Entry in row',I, &
          ' and column',J,' duplicated'
      GO TO 750
  730 INFO(1) = -5
      IF (LP.GT.0) WRITE (LP,'(/A/(3(A,I9)))') &
          ' **** Error return from Analyse ****',' Fault in component ', &
          PIVOT,' of column permutation given in IFIRST'
      GO TO 750
  740 INFO(1) = -6
      IF (LP.GT.0) WRITE (LP,'(/A/(3(A,I9)))') &
          ' **** Error return from Analyse ****',' ICNTL(4) = ', &
          ICNTL(4),' when ICNTL(6) = 2'
  750 CONTINUE

      END SUBROUTINE MP48_MA60AD

!********************************************************
      SUBROUTINE MP48_MA60BD(M,N,NE,JOB,AA,IRNA,IPTRA,CNTL,ICNTL,IP, &
                             IQ,NP,RANK,NNA,LFACT,FACT,LIRNF,IRNF, &
                             IPTRL,IPTRU,W,IW,INFO,RINFO)

! MP48_MA60B/BD factorizes the matrix in AA/IRNA/IPTRA as P L U Q where
!     P and Q are permutations, L is lower triangular, and U is unit
!     upper triangular. The prior information that it uses depends on
!     the value of the parameter JOB.

      USE HSL_MP48_DATA_DOUBLE
      IMPLICIT NONE

      INTEGER M,N,NE,JOB,NNA,NP,LFACT,LIRNF,RANK
      REAL(WP) AA(NE),CNTL(4)
      INTEGER IRNA(NE),IPTRA(N),IW(M+2*N),INFO(20)
      INTEGER ICNTL(10),IP(M),IQ(*),IRNF(LIRNF),IPTRL(N),IPTRU(N)
      REAL(WP) RINFO,FACT(LFACT),W(M)

! M is an integer variable that must be set to the number of rows.
!      It is not altered by the subroutine.
! N is an integer variable that must be set to the number of columns.
!      It is not altered by the subroutine.
! NE is an integer variable that must be set to the number of entries
!      in the input matrix.  It is not altered by the subroutine.
! JOB is an integer variable that must be set to the value 1, 2, or 3.
!     If JOB is equal to 1 and any of the first NP recommended pivots
!      fails to satisfy the threshold pivot tolerance, the row is
!      interchanged with the earliest row in the recommended sequence
!      that does satisfy the tolerance. Normal row interchanges are
!      performed in the last N-NP columns.
!     If JOB is equal to 2, then M, N, NE, IRNA, IPTRA, IP, IQ,
!      IRNF, IPTRL, and IPTRU must be unchanged since a
!      JOB=1 entry for the same matrix pattern and no interchanges are
!      performed among the first NP pivots.
!     JOB is not altered by the subroutine.
! AA is an array that holds the entries of the matrix and
!      is not altered.
! IRNA is an integer array of length NE that must be set to hold the
!      row indices of the corresponding entries in AA. It is not
!      altered.
! IPTRA is an integer array that holds the positions of the starts of
!      the columns of AA. It is not altered by the subroutine.
! CNTL  must be set by the user as follows and is not altered.
!     CNTL(2) determines the balance between pivoting for sparsity and
!       for stability, values near zero emphasizing sparsity and values
!       near one emphasizing stability.
!     CNTL(3) If this is set to a positive value, any entry whose
!       modulus is less than CNTL(3) will be dropped from the factors.
!       The factorization will then require less storage but will be
!       inaccurate.
!     CNTL(4)  Any entry of the reduced matrix whose modulus is less
!       than or equal to CNTL(4) will be regarded as zero from the
!        point of view of rank.
! ICNTL must be set by the user as follows and is not altered.
!     ICNTL(1)  must be set to the stream number for error messages.
!       A value less than 1 suppresses output.
!     ICNTL(2) must be set to the stream number for diagnostic output.
!       A value less than 1 suppresses output.
!     ICNTL(3) must be set to control the amount of output:
!       0 None.
!       1 Error messages only.
!       2 Error and warning messages.
!       3 As 2, plus scalar parameters and a few entries of array
!         parameters on entry and exit.
!       4 As 2, plus all parameters on entry and exit.
!     ICNTL(5) The block size to be used for full-matrix processing.
!       If <=0, the BLAS1 version is used.
!       If =1, the BLAS2 version is used.
!     ICNTL(6) This should be set by the user to the number of columns
!       in the border.
! IP is an integer array. If JOB=1, it must be set so that IP(I) < IP(J)
!      if row I is recommended to precede row J in the pivot sequence.
!      If JOB>1, it need not be set. If JOB=1 or JOB=3, IP(I) is set
!      to -K when row I is chosen for pivot K and IP is eventually
!      reset to recommend the chosen pivot sequence to a subsequent
!      JOB=1 entry. If JOB=2, IP is not be referenced.
! IQ is an integer array that must be set so that either IQ(J) is the
!      column in position J in the pivot sequence, J=1,2,...,N,
!      or IQ(1)=0 and the columns are taken in natural order.
!      It is not altered by the subroutine.
! NNA is an integer variable that holds the number of columns
!      that are not pivoted on.  It will be greater or equal to
!      ICNTL(6).
! RANK Need not be set on entry.  Is set by subroutine to the number
!      This value is returned by MP48_MA60E/ED or MP48_MA60F/FD
!      of pivots chosen in the dense matrix.
! NP is an integer variable that holds the number of columns
!      processed in packed storage. It is not altered by the subroutine.
! LFACT is an integer variable set to the size of FACT.
!      It is not altered by the subroutine.
! FACT is an array that need not be set on a JOB=1 entry and must be
!      unchanged since the previous entry if JOB>1. On return, FACT(1)
!      holds the value of CNTL(3) used, FACT(2) will holds the value
!      of CNTL(4) used, FACT(3:IPTRL(N)) holds the packed part of L/U
!      by columns, and the full part of L/U is held by columns
!      immediately afterwards. U has unit diagonal entries, which are
!      not stored. In each column of the packed part, the entries of
!      U precede the entries of L; also the diagonal entries of L
!      head each column of L and are reciprocated.
! LIRNF is an integer variable set to the size of IRNF.
!      It is not altered by the subroutine.
! IRNF is an integer array of length LIRNF that need not be set on
!      a JOB=1 entry and must be unchanged since the previous entry
!      if JOB>1. On exit, IRNF(1) holds the number of dropped entries,
!      IRNF(2) holds the number of rows MF in full storage,
!      IRNF(3:IPTRL(N)) holds the row numbers of the packed part
!      of L/U, IRNF(IPTRL(N)+1:IPTRL(N)+MF) holds the row indices
!      of the full part of L/U, and IRNF(IPTRL(N)+MF+I), I=1,2,..,N-NP
!      holds the vector IPIV output by MP48_MA60GD.
!      If JOB=2, IRNF will be unaltered.
! IPTRL is an integer array that need not be set on a JOB=1 entry and
!     must be unchanged since the previous entry if JOB>1.
!     For J = 1,..., NP, IPTRL(J) holds the position in
!     FACT and IRNF of the end of column J of L.
!     For J = NP+1,..., N, IPTRL(J) is equal to IPTRU(J).
! IPTRU is an integer array that need not be set on a JOB=1 entry and
!     must be unchanged since the previous entry if JOB>1.
!     For J = 1,..., N, IPTRU(J) holds the position in
!     FACT and IRNF of the end of the packed part of column J of U.
! W is an array of length M used as workspace for holding
!      the expanded form of a sparse vector.
! IW is an integer array of length M+2*N used as workspace.
! INFO need not be set on entry. On exit, it holds the following:
!    INFO(1) A negative value will indicate an error return and a
!       positive value a warning. Possible nonzero values are:
!      -1  M < 1 or N < 1.
!      -2  NE < 0.
!      -3  Insufficient space.
!      -4  There are duplicated entries.
!      -5  JOB < 1 or > 2.
!      -6  JOB = 2, but entries were dropped in the corresponding JOB=1
!          entry.
!      -7  NP < 0 or NP > N.
!     -(7+K) Pivot too small in column K when JOB=2.
!      &1  Rank deficient.
!    INFO(4) Minimum real storage required to factorize matrix (or
!            recommended value for LFACT if INFO(1) = -3).  After
!            JOB = 1 entry or for subsequent JOB=2 entry.
!    INFO(5) Computed rank of the matrix.
!    INFO(6) Number of entries dropped from the data structure.
!    INFO(7) Number of rows processed in full storage.
!    INFO(8) Number of reals in factors.
!    INFO(9) Number of integers in factors.
!    INFO(10)Minimum integer storage required to factorize matrix (or
!            recommended value for LIRNF if INFO(1) = -3).  After
!            JOB = 1 entry or for subsequent JOB=2 entry.
! RINFO need not be set on entry. On exit, it holds the number of
!    floating-point operations performed.

      REAL(WP), PARAMETER :: ZERO = 0.0_WP
      REAL(WP), PARAMETER :: ONE = 1.0_WP

      REAL(WP) AMULT,ASW
      INTEGER BEGCOL
      LOGICAL DROP
      INTEGER ENDCOL,EYE,EYE1,I,IA1,IA2,IF1I,IF2I,IF1R,IF2R,II,IL1, &
              IL2,IPIV,IQPIV, &
              IU1,IU2,ISW,J,JDUMMY,JJ,K,LP
      REAL(WP) MAXENT
      INTEGER MF,MORD,MP,NEU,NF,NNO,NULLC
      REAL(WP) PIVLIM
      REAL(WP) U
! AMULT Temporary variable used to store current multiplier.
! ASW Temporary variable used when swopping two real quantities.
! BEGCOL is allocatable to beginning of section of column when pruning.
! DROP True if any entries dropped from current column.
! ENDCOL is allocatable to end of section of column when pruning.
! EYE Running position for current column.
! EYE1 Position of the start of second current column.
! I Temporary variable holding row number. Also used as index in DO
!     loops used in initialization of arrays.
! IA1 Position of the start of the current column in AA.
! IA2 Position of the end of the current column in AA.
! IF1I Position of the start of the full submatrix (integer storage).
! IF2I Position of the end of the full submatrix (integer storage).
! IF1R Position of the start of the full submatrix (real storage).
! IF2R Position of the end of the full submatrix (real storage).
! II Running position for current column.
! IL1 Position of the first entry of the current column of L.
! IL2 Position of the last entry of the current column of L.
! IPIV Position of the pivot in FACT and IRNF.
! IQPIV Recommended position of the pivot row in the pivot sequence.
! IU1 Position of the start of current column of U.
! IU2 Position of the end of the current column of U.
! ISW Temporary variable used when swopping two integer quantities.
! J Temporary variable holding column number.
! JDUMMY DO index not referenced in the loop.
! JJ Running position for current column.
! K Temporary variable holding the current pivot step in the elimination
! LP Unit for error messages.
! MAXENT Temporary variable used to hold value of largest entry in
!    column.
! MF Number of rows in full block.
! MORD Number of rows ordered.
! MP Unit for diagnostic messages.
! NEU Number of entries omitted from U and the full block in order to
!    calculate INFO(4) (0 unless INFO(1)=-3).
! NF Number of columns in full block.
! NNO is an integer variable that holds the number of columns
!      that are ineligible for pivoting because they are in the border.
! NULLC Number of columns found null before dropping any elements.
! PIVLIM Limit on pivot size.
! U Used to hold local copy of CNTL(2), changed if necessary so that it
!    is in range.

      INTRINSIC ABS,MAX,MIN
! LAPACK subroutine for triangular factorization.

      INFO(1) = 0
      INFO(5) = 0
      INFO(6) = 0
      INFO(7) = 0
! Initialize real and integer space required for factorization
      INFO(4) = MAX(M,NE+2)
      INFO(10) = MAX(M,NE+2)
      RINFO = ZERO
      LP = ICNTL(1)
      MP = ICNTL(2)
      IF (ICNTL(3).LE.0) LP = 0
      IF (ICNTL(3).LE.1) MP = 0

      NNO = ICNTL(6)
      NNA = NNO

! Check input values
      IF (M.LT.1 .OR. N.LT.1) THEN
         INFO(1) = -1
         IF (LP.GT.0) WRITE (LP,'(/A/A,I8,A,I8)') &
             ' **** Error return from Factorise ****',' M =',M,' N =',N
         GO TO 550
      END IF
      IF (NE.LE.0) THEN
         INFO(1) = -2
         IF (LP.GT.0) WRITE (LP,'(/A/A,I6)') &
             ' **** Error return from Factorise ****',' NE =',NE
         GO TO 550
      END IF
      IF (NP.LT.0 .OR. NP.GT.N) THEN
         INFO(1) = -7
         IF (LP.GT.0) WRITE (LP,'(/A/A,I8,A,I8)') &
             ' **** Error return from Factorise ****', &
             ' NP =',NP,' N =',N
         GO TO 550
      END IF
      IF (LFACT.LT.MAX(M,NE+2) .OR. LIRNF.LT.MAX(M,NE+2)) THEN
         INFO(4) = MAX(M,NE+2)
         INFO(10) = MAX(M,NE+2)
         GO TO 520
      END IF
      IF (JOB.EQ.1) THEN
      ELSE IF (JOB.EQ.2 .OR. JOB.EQ.3) THEN
         IF (IRNF(1).NE.0) THEN
            INFO(1) = -6
            IF (LP.GT.0) WRITE (LP,'(/A/A,I1,A)') &
                ' **** Error return from Factorise ***', &
                ' Call with JOB=', &
                JOB,' follows JOB=1 call in which entries were dropped'
            GO TO 550
         END IF
      ELSE
         INFO(1) = -5
         IF (LP.GT.0) WRITE (LP,'(/A/A,I2)') &
             ' **** Error return from Factorise ****',' JOB =',JOB
         GO TO 550
      END IF

! Print input data
      IF (MP.GT.0) THEN
         IF (ICNTL(3).GT.2) WRITE (MP, &
             '(/2(A,I6),A,I8,A,I3/2(A,I8),A,I7/A,1P,4E10.2/A,8I4)') &
             ' Entering Factorise with M =',M,' N =',N,' NE =',NE, &
             ' JOB =',JOB,' LFACT =',LFACT,' LIRNF =',LIRNF, &
           ' NP =',NP,' CNTL =',CNTL,' ICNTL =',(ICNTL(J),J=1,8)
         IF (ICNTL(3).GT.3) THEN
            WRITE (MP,'(A,(T6,10(I7)))') ' IP = ',IP
            IF (IQ(1).GT.0) THEN
               WRITE (MP,'(A,(T6,10(I7)))') ' IQ = ',(IQ(J),J=1,N)
            ELSE
               WRITE (MP,'(A,(T6,I7))') ' IQ = ',IQ(1)
            END IF
            DO 10 J = 1,N - 1
               IF (IPTRA(J).LT.IPTRA(J+1)) WRITE (MP, &
                   '(A,I5,(T13,3(1P,E12.4,I5)))') ' Column',J, &
                   (AA(II),IRNA(II),II=IPTRA(J),IPTRA(J+1)-1)
   10       CONTINUE
            IF (IPTRA(N).LE.NE) WRITE (MP, &
                '(A,I5,(T13,3(1P,E12.4,I5)))') ' Column',N, &
                (AA(II),IRNA(II),II=IPTRA(N),NE)
         END IF
      END IF

! Initializations.
      NULLC = 0

      U = MIN(CNTL(2),ONE)
      U = MAX(U,ZERO)
      DO 20 I = 1,M
         IW(I+N) = 0
         W(I) = ZERO
   20 CONTINUE
      MORD = 0
      IF1I = LIRNF + 1
      IF2I = 0
      IF1R = LFACT + 1
      NF = N - NP
      MF = 0
      IL2 = 2
      NEU = 0

! Jump if JOB is equal to 2.
      IF (JOB.EQ.2) GO TO 370

      IF (JOB.EQ.3) THEN
! Reconstruct IP and set MORD
         DO 30 J = 1,NP
            IA1 = IPTRU(J) + 1
            IF (IA1.GT.IPTRL(J)) GO TO 30
            IP(IRNF(IA1)) = J
   30    CONTINUE
         MF = IRNF(2)
         IA1 = IPTRL(N)
         DO 40 J = 1,MF
            IP(IRNF(IA1+J)) = NP + J
   40    CONTINUE
      END IF

! Store copies of column ends ready for pruning
!     DO 50 K = 1,JLAST
!        IW(M+N+K) = IPTRL(K)
!  50 CONTINUE

! Each pass through this main loop processes column K.
      DO 310 K = 1,N
         DROP = .FALSE.
         IF (K.EQ.NP+1) THEN
! Set up data structure for full part.
            MF = M - MORD
            IF1I = LIRNF + 1 - MF
            IF1R = LFACT + 1 - MF*NF
            II = 0
            DO 60 I = 1,M
               IF (IP(I).GT.0) THEN
                  IW(I+N) = N
                  IRNF(IF1I+II) = I
                  II = II + 1
                  IP(I) = NP + II
               END IF
   60       CONTINUE
            IF1I = LIRNF + 1 - (MF+NF)
            IF2I = IF1I - 1
! IF2R points to position before first entry of full block
            IF2R = LFACT - MF*NF
         END IF
         J = K
         IF (IQ(1).GT.0) J = IQ(K)
         IA1 = IPTRA(J)
         IA2 = NE
         IF (J.NE.N) IA2 = IPTRA(J+1) - 1
!        write (6,*) 'K,J,IF1I,IF1R,IA1,IA2,NP', &
!                    K,J,IF1I,IF1R,IA1,IA2,NP
         IU1 = IL2 + 1
         IU2 = IU1 - 1
         IL1 = IF1I - 1 + IA1 - IA2
         IL2 = IL1 - 1
!        write (6,*) 'K,INFO(4),NEU,LFACT,IL1,IL2,M', &
!                    K,INFO(4),NEU,LFACT,IL1,IL2,M
         INFO(4) = MAX(INFO(4),NEU+LFACT-IL1+IU2+M+1)
!        write (6,*) 'INFO(10),NEU,LIRNF,IL1,IU2,M', &
!                    INFO(10),NEU,LIRNF,IL1,IU2,M
         INFO(10) = MAX(INFO(10),NEU+LIRNF-IL1+IU2+M+1)
         IF (IL1-IU2.LE.M) THEN
            IF (INFO(1).NE.-3) THEN
! Get rid of U info.
               INFO(1) = -3
               NEU = IL2 + LIRNF + 1 - MF - IF1I
               IF1I = LIRNF + 1 - MF
               IF2I = IF1I - 1
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
               IL1 = IF1I - 1 + IA1 - IA2
               IL2 = IL1 - 1
            END IF
! Quit if LIRNF is much too small
            IF (IL1-IU2.LE.M) GO TO 480
         END IF
! Load column K of AA into full vector W and into the back of IRNF.
! Check for duplicates.
         EYE = IL1
         DO 90 II = IA1,IA2
            I = IRNA(II)
            IF (IW(I+N).EQ.-1) GO TO 540
            IW(I+N) = -1
            W(I) = AA(II)
            IRNF(EYE) = I
            EYE = EYE + 1
   90    CONTINUE
! Depth first search to find topological order for triangular solve
!     and structure of column K of L/U
! IW(J) is used to hold a allocatable to next entry in column J
!     during the depth-first search at stage K, J = 1,..., N.
! IW(I+N) is set to K when row I has been processed, and to N for rows
!     of the full part once column NP has been passed. It is also
!     used for backtracking, a negative value being used to point to the
!     previous row in the chain.
! IW(M+N+I) is set to the position in FACT and IRNF of the end of the
!     active part of the column after pruning.  It is initially set to
!     IPTRL(I) and is flagged negative when column has been pruned.
! Set IPTRL temporarily for column K so that special code is
!     not required to process this column.
         IPTRL(K) = EYE - 1
         IW(M+N+K) = EYE - 1
! IW(K) is set to beginning of original column K.
         IW(K) = IL1
         J = K
! The outer loop of the depth-first search is executed once for column
!      K and twice for each entry in the upper-triangular part of column
!      K (once to initiate a search in the corresponding column and
!      once when the search in the column is finished).
         DO 120 JDUMMY = 1,2*K
! Look through column J of L (or column K of A). All the entries
!     are entries of the filled-in column K. Store new entries of the
!     lower triangle and continue until reaching an entry of the upper
!     triangle.
            DO 100 II = IW(J),ABS(IW(M+N+J))
               I = IRNF(II)
! Jump if index I already encountered in column K or is in full part.
               IF (IW(I+N).GE.K) GO TO 100
               IF (IP(I).LE.0) GO TO 110
! Entry is in lower triangle. Flag it and store it in L.
               IW(I+N) = K
               IL1 = IL1 - 1
               IRNF(IL1) = I
  100       CONTINUE
            IF (J.EQ.K) GO TO 130
! Flag J, put its row index into U, and backtrack
            IU2 = IU2 + 1
            I = IRNF(IPTRU(J)+1)
            IRNF(IU2) = I
            J = -IW(I+N)
            IW(I+N) = K
            GO TO 120
! Entry in upper triangle.  Move search to corresponding column.
  110       IW(I+N) = -J
            IW(J) = II + 1
            J = -IP(I)
            IW(J) = IPTRU(J) + 2
  120    CONTINUE
! Run through column K of U in the lexicographical order that was just
!     constructed, performing elimination operations.
  130    DO 150 II = IU2,IU1,-1
            I = IRNF(II)
            J = -IP(I)
! Add multiple of column J of L to column K
            EYE1 = IPTRU(J) + 1
            IF (ABS(W(I)).LT.CNTL(3)) GO TO 150
            AMULT = -W(I)*FACT(EYE1)
! Note we are storing negative multipliers
            W(I) = AMULT
            DO 140 EYE = EYE1 + 1,IPTRL(J)
               I = IRNF(EYE)
               W(I) = W(I) + AMULT*FACT(EYE)
  140       CONTINUE
            RINFO = RINFO + ONE + 2*(IPTRL(J)-EYE1)
  150    CONTINUE

! Unload reals of column of U and set allocatable
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
! Find the largest entry in the column and drop any small entries
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
! Unload column of L, performing pivoting and moving indexing
!      information.
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
! Find position of pivot
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
! Column is not null
               IF (IQPIV.EQ.M+N) THEN
! All entries in the column are too small to be pivotal. Drop them all.
                  IF (CNTL(3).GT.ZERO) INFO(6) = INFO(6) + EYE - IU2
                  IL2 = IU2
               ELSE
                  IF (IL1.NE.IPIV) THEN
! Move pivot to front of L
                     ASW = FACT(IPIV)
                     FACT(IPIV) = FACT(IL1)
                     FACT(IL1) = ASW
                     ISW = IRNF(IL1)
                     IRNF(IL1) = IRNF(IPIV)
                     IRNF(IPIV) = ISW
                  END IF
! Reciprocate pivot
                  INFO(5) = INFO(5) + 1
                  FACT(IL1) = ONE/FACT(IL1)
                  RINFO = RINFO + ONE
! Record pivot row
                  MORD = MORD + 1
                  IP(IRNF(IL1)) = -K
               END IF
            END IF
         ELSE
! Treat column as full
            IL2 = IPTRU(K)
            DO 210 II = LIRNF - MF + 1,LIRNF
               I = IRNF(II)
               IF2R = IF2R + 1
               FACT(IF2R) = W(I)
               W(I) = ZERO
  210       CONTINUE
            IF (INFO(1).EQ.-3) IF2R = IF2R - MF
         END IF
         IW(M+N+K) = IL2
         IPTRL(K) = IL2
         IF (DROP) GO TO 310
! Scan columns involved in update of column K and remove trailing block.
         DO 300 II = IU1,IU2
            I = IRNF(II)
! Jump if column already pruned.
            IF (IW(M+N+I).LT.0) GO TO 300
            BEGCOL = IPTRU(I) + 2
            ENDCOL = IPTRL(I)
! Scan column to see if there is an entry in the current pivot row.
            IF (K.LE.NP) THEN
               DO 220 JJ = BEGCOL,ENDCOL
                  IF (IP(IRNF(JJ)).EQ.-K) GO TO 230
  220          CONTINUE
               GO TO 300
            END IF
! Sort the entries so that those in rows already pivoted (negative IP
!    values) precede the rest.
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
! Set up data structure for the (null) full part.
         MF = M - MORD
         IF1I = LIRNF + 1 - MF
         II = 0
         DO 320 I = 1,M
            IF (IP(I).GT.0) THEN
               IW(I+N) = N
               IRNF(IF1I+II) = I
               II = II + 1
               IP(I) = NP + II
            END IF
  320    CONTINUE
         IF1I = LIRNF + 1 - (MF+NF)
         IF1R = LFACT - (MF*NF)
         IF2I = IF1I - 1
      END IF
      IF (INFO(5).EQ.MIN(M,N)) THEN
! Restore sign of IP
         DO 330 I = 1,M
            IP(I) = ABS(IP(I))
  330    CONTINUE
      ELSE
! Complete IP
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
! Move full part forward

! For this code we now finish leaving Schur complement as dense block
! in last positions of arrays FACT and IRNF.
      IF2R = IF2R - MF*NF
      IF2I = IF2I - (MF+NF)
      DO 350 II = 1,MF*NF
         FACT(IL2+II) = FACT(IF1R-1+II)
  350 CONTINUE
      DO 360 II = 1,MF
         IRNF(IL2+II) = IRNF(LIRNF-MF+II)
  360 CONTINUE
      IF1R = IL2 + 1
      IF1I = IL2 + 1
      GO TO 440

! Fast factor (JOB = 2)
! Each pass through this main loop processes column K.
  370 MF = IRNF(2)
      INFO(7) = MF
      IF1R = IPTRL(N) + 1
      IF1I = IPTRL(N) + 1
      IF2I = IF1I - 1
      IF2R = IF1R - 1
      DO 430 K = 1,N
         J = K
         IF (IQ(1).GT.0) J = IQ(K)
         IA1 = IPTRA(J)
         IA2 = NE
         IF (J.NE.N) IA2 = IPTRA(J+1) - 1
         IU1 = IL2 + 1
         IU2 = IPTRU(K)
         IL1 = IU2 + 1
         IL2 = IPTRL(K)
! Load column K of A into full vector W
         DO 380 II = IA1,IA2
            W(IRNA(II)) = AA(II)
  380    CONTINUE
! Run through column K of U in lexicographical order, performing
!      elimination operations.
         DO 400 II = IU2,IU1,-1
            J = IRNF(II)
            I = IRNF(IPTRU(J)+1)
! Add multiple of column J of L to column K
            EYE1 = IPTRU(J) + 1
            AMULT = -W(I)*FACT(EYE1)
! Note we are storing negative multipliers
            FACT(II) = AMULT
            W(I) = ZERO
            DO 390 EYE = EYE1 + 1,IPTRL(J)
               I = IRNF(EYE)
               W(I) = W(I) + AMULT*FACT(EYE)
  390       CONTINUE
            RINFO = RINFO + ONE + 2*(IPTRL(J)-EYE1)
  400    CONTINUE
         IF (K.LE.NP) THEN
            IF (IL1.LE.IL2) THEN
! Load column of L.
               DO 410 II = IL1,IL2
                  I = IRNF(II)
                  FACT(II) = W(I)
                  W(I) = ZERO
  410          CONTINUE
! Test pivot. Note that this is the only numerical test when JOB = 2.
               IF (ABS(FACT(IL1)).LE.CNTL(4)) THEN
                  GO TO 530
               ELSE
! Reciprocate pivot
                  INFO(5) = INFO(5) + 1
                  FACT(IL1) = ONE/FACT(IL1)
                  RINFO = RINFO + ONE
               END IF
            END IF
         ELSE
! Treat column as full
            DO 420 II = IF1I,IF1I + MF - 1
               I = IRNF(II)
               IF2R = IF2R + 1
               FACT(IF2R) = W(I)
               W(I) = ZERO
  420       CONTINUE
         END IF
  430 CONTINUE
! Set real and integer space for subsequent factorization
!      write (6,*) 'INFO(4),IF2R',INFO(4),IF2R
      INFO(4)  = MAX(INFO(4),IF2R)
!      write (6,*) 'INFO(10),IF1I+MF+NF-1',INFO(4),IF1I+MF+NF-1
      INFO(10) = MAX(INFO(10),IF1I+MF+NF-1)

  440 IF (MF.GT.0 .AND. NF.GT.0) THEN
! Factorize full block
         IF (ICNTL(5).GT.1) &
            CALL MP48_MA60GD(MF,NF,FACT(IF1R),MF,ICNTL(5), &
                             CNTL(4),IRNF(IF1I+MF),RANK,NNO)
         IF (ICNTL(5).EQ.1) &
            CALL MP48_MA60FD(MF,NF,FACT(IF1R),MF,CNTL(4), &
                             IRNF(IF1I+MF),RANK,NNO)
         IF (ICNTL(5).LE.0) &
            CALL MP48_MA60ED(MF,NF,FACT(IF1R),MF,CNTL(4), &
                             IRNF(IF1I+MF),RANK,NNO)
         NNA = NF - RANK
         INFO(4)  = MAX(INFO(4), MF*NF+2)
         INFO(10) = MAX(INFO(10), MF+NF+2)
         INFO(5)  = INFO(5) + RANK
         DO 450 I = 1,MIN(MF,NF)
            RINFO = RINFO + MF - I + 1 + REAL(MF-I)*(NF-I)*2
  450    CONTINUE
      ELSE
        RANK = 0
      END IF
! Total number reals
      INFO(8) = IF1R+(MF*NF)-1
      INFO(4) = MAX(INFO(8),M,NE+2)
      INFO(4) = MAX(INFO(8),INFO(4))
! Total number integers
      INFO(9) = IF1I+MF+NF-1-2
!     INFO(10)= MAX(INFO(9)+2,M,NE+2)
      IF (INFO(5).LT.MIN(M,N)) INFO(1) = 1
      IF (MP.GT.0 .AND. ICNTL(3).GT.2) THEN
         WRITE (MP,'(A,I6,A,F12.1/A,I3/A,7I8)') &
           ' Leaving Factorise with IRNF(2) =',IRNF(2), &
           ' RINFO =',RINFO, &
           ' INFO(1) =',INFO(1),' INFO(4:10) =', (INFO(J),J=4,10)
         IF (ICNTL(3).GT.3) THEN
            IF (JOB.NE.2) WRITE (MP,'(A,(T6,10(I7)))') ' IP = ',IP
            DO 460 J = 1,N
               IF (J.GT.1) THEN
                  IF (IPTRL(J-1).LT.IPTRU(J)) WRITE (MP, &
                      '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column',J, &
                      ' of U', (FACT(II),IRNF(II),II=IPTRL(J-1)+1, &
                      IPTRU(J))
               END IF
               IF (IPTRU(J).LT.IPTRL(J)) WRITE (MP, &
                   '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column',J,' of L', &
                    (FACT(II),IRNF(II),II=IPTRU(J)+1,IPTRL(J))
  460       CONTINUE
            IF (MF.GT.0) THEN
              WRITE (MP,'(A)') ' Full part'
              WRITE (MP,'((6I12))') (IRNF(IF1I+MF+J),J=0,NF-1)
              DO 470 I = 0,MF - 1
                 WRITE (MP,'(I4,1P,6E12.4:/(4X,1P,6E12.4))') &
                  IRNF(IF1I+I),(FACT(IF1R+I+J*MF),J=0,NF-1)
!                 do j = 0,nf-1
!                   write (6,'(6i4,d12.4)')
!     &             IPTRL(n),if1,j,mf,nf,IF1I+I+J*MF,FACT(IF1I+I+J*MF)
!                 end do
  470          CONTINUE
            ELSE
               WRITE (MP,'(A)') ' No Full part .. zero rows'
            END IF
         END IF
      END IF
      GO TO 550

! Error conditions
! LIRNF is much too small. Patch up IP and quit.
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
         WRITE (LP,'(/A/A,I7,A,I7/A,I7,A,I7)') &
           ' **** Error return from Factorise **** ', &
           ' either LFACT must be increased from ',LFACT, &
           ' to at least',INFO(4), &
           ' or     LIRNF must be increased from ', &
             LIRNF,' to at least',INFO(10)
      END IF
      GO TO 550
  530 INFO(1) = - (7+K)
      IF (LP.GT.0) WRITE (LP,'(/A/A,I6,A)') &
          ' **** Error return from Factorise **** ', &
          ' Small pivot found in column',K,' of the permuted matrix.'
      GO TO 550
  540 INFO(1) = -4
      IF (LP.GT.0) WRITE (LP,'(/A/(3(A,I9)))') &
          ' **** Error return from Factorise ****',' Entry in row', &
          I,' and column',J,' duplicated'
  550 CONTINUE

      END SUBROUTINE MP48_MA60BD

!****************************************************
      SUBROUTINE MP48_MA60DD(LA,A,IND,IPTR,N,DISP,REALS)
! This subroutine performs garbage collection on the arrays A and IND.
! DISP is the position in arrays A/IND immediately after the data
!     to be compressed.
!     On exit, DISP equals the position of the first entry
!     after the compressed part of A/IND.

      USE HSL_MP48_DATA_DOUBLE
      IMPLICIT NONE

      INTEGER LA,N,DISP
      REAL(WP) A(LA)
      INTEGER IPTR(N)
      LOGICAL REALS
      INTEGER IND(LA)
! Local variables.
      INTEGER ISW,J,K,KN
! Set the first entry in each row(column) to the negative of the
! row(column) and hold the column(row) index in the row(column)
! allocatable.  This enables the start of each row(column) to be
! recognized in a subsequent scan.
      ISW = 0
      DO 10 J = 1,N
         K = IPTR(J)
         IF (K.GT.0) THEN
           IF (IND(K).GE.0) THEN
             IPTR(J) = IND(K)
             IND(K) = -J
           ELSE
             ISW = 1
             IPTR(J) = -(N+1)
           END IF
         END IF
   10 CONTINUE
      KN = 0
! Go through arrays compressing to the front so that there are no
! zeros held in positions 1 to DISP-1 of IND.
! Reset first entry of each row(column) and the allocatable array IPTR.
      DO 20 K = 1,DISP - 1
         IF (IND(K).EQ.0) GO TO 20
         KN = KN + 1
         IF (REALS) A(KN) = A(K)
         IF (IND(K).LE.0) THEN
! First entry of row(column) has been located.
            J = -IND(K)
            IND(K) = IPTR(J)
            IPTR(J) = KN
         END IF
         IND(KN) = IND(K)
   20 CONTINUE
      DISP = KN + 1
      IF (ISW.EQ.1) THEN
        DO 30 J=2,N
          IF (IPTR(J).EQ.-(N+1)) THEN
            IPTR(J) = IPTR(J-1)
          END IF
   30   CONTINUE
      END IF

      END SUBROUTINE MP48_MA60DD

!********************************************************
      SUBROUTINE MP48_MA60ED(M,N,A,LDA,PIVTOL,IPIV,RANK,NNO)

      USE HSL_MP48_DATA_DOUBLE
      IMPLICIT NONE

      INTEGER LDA,M,N,RANK,NNO
      REAL(WP) PIVTOL

      INTEGER IPIV(N)
      REAL(WP) A(LDA,N)

!  Purpose
!  MP48_MA60ED computes an LU factorization of a general m-by-n
!  matrix A.

!  The factorization has the form
!     A = P * L * U * Q
!  where P is a permutation matrix of order m, L is lower triangular
!  of order m with unit diagonal elements, U is upper trapezoidal of
!  order m * n, and Q is a permutation matrix of order n.
!
!  Row interchanges are used to ensure that the entries of L do not
!  exceed 1 in absolute value. Column interchanges are used to
!  ensure that the first r diagonal entries of U exceed PIVTOL in
!  absolute value. If r < m, the last (m-r) rows of U are zero.

!  This is the Level 1 BLAS version.

!  Arguments

!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 1.
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 1.
!  A       (input/output) REAL(WP) array, dimension (LDA,N)
!          On entry, the m by n matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U*Q; the unit diagonal elements of L are not stored.
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!  PIVTOL  (input) REAL(WP)
!          The pivot tolerance. Any entry with absolute value less
!          than or equal to PIVTOL is regarded as unsuitable to be a
!          pivot.
!  IPIV    (output) INTEGER array, dimension (N)
!          The permutations; for 1 <= i <= RANK, row i of the
!          matrix was interchanged with row IPIV(i); for
!          RANK + 1 <= j <= N, column j of the
!          matrix was interchanged with column -IPIV(j).
!  RANK    (output) INTEGER
!          The computed rank of the matrix.
!  NNO     Number of columns (at end) from which pivots cannot be chosen

      REAL(WP), PARAMETER :: ZERO = 0.0_WP
      REAL(WP), PARAMETER :: ONE = 1.0_WP

      INTEGER I,JP,K,NCOL
      LOGICAL PIVOT
! I   Row index.
! JP  Pivot position.
! K   Main loop index.
! NCOL Number of columns from which pivots can be chosen
! PIVOT True if there is a pivot in current column.

      INTEGER IDAMAX
      EXTERNAL IDAMAX

      EXTERNAL DAXPY,DSCAL,DSWAP
      INTRINSIC ABS

      NCOL = N - NNO

      RANK = 0
      DO 30 K = 1,N-NNO

! Update elements in column RANK+1.
         DO 10 I = 1, RANK
            IF (I.LE.MIN(M-1,NCOL)) &
            CALL DAXPY(M-I,-A(I,RANK+1),A(I+1,I),1,A(I+1,RANK+1),1)
   10    CONTINUE

! Find pivot.
         JP = RANK + IDAMAX(M-RANK,A(RANK+1,RANK+1),1)
         PIVOT = ABS(A(JP,RANK+1)) .GT. PIVTOL

         IF (PIVOT) THEN

! Update RANK
            RANK = RANK + 1

            IPIV(RANK) = JP
! Apply row interchange to columns 1:N
            IF (JP.NE.RANK) &
              CALL DSWAP(N,A(RANK,1),LDA,A(JP,1),LDA)

! Compute elements RANK+1:M of RANK-th column.
            IF (RANK.LT.M) &
              CALL DSCAL(M-RANK,ONE/A(RANK,RANK),A(RANK+1,RANK),1)

         ELSE

            DO 20 I = RANK+1,M
               A(I,RANK+1) = ZERO
   20       CONTINUE
! Apply column interchange and record it.
            IF (RANK+1.LT.NCOL) CALL DSWAP(M,A(1,RANK+1),1,A(1,NCOL),1)
            IPIV(NCOL) = -(RANK+1)
            NCOL = NCOL - 1

         END IF

         IF (RANK.GE.MIN(NCOL,M)) GO TO 40

   30 CONTINUE

! Update border columns and set IPIV for them
   40 DO 50 K = N-NNO+1,N
         DO 45 I = 1, RANK
            IF (I.LT.M) &
            CALL DAXPY(M-I,-A(I,K),A(I+1,I),1,A(I+1,K),1)
   45    CONTINUE
         IPIV(K) = -K
   50 CONTINUE

! Set IPIV for remaining columns
      DO 55 I = RANK+1,NCOL
        IPIV(I) = -I
   55 CONTINUE

      END SUBROUTINE MP48_MA60ED

!********************************************************
      SUBROUTINE MP48_MA60FD(M,N,A,LDA,PIVTOL,IPIV,RANK,NNO)

!  -- This is a variant of the LAPACK routine DGETF2 --

      USE HSL_MP48_DATA_DOUBLE
      IMPLICIT NONE

      INTEGER LDA,M,N,RANK,NNO
      REAL(WP) PIVTOL

      INTEGER IPIV(N)
      REAL(WP) A(LDA,N)


!   Purpose
!  MP48_MA60FD computes an LU factorization of a general m-by-n
!  matrix A.

!  The factorization has the form
!     A = P * L * U * Q
!  where P is a permutation matrix of order m, L is lower triangular
!  of order m with unit diagonal elements, U is upper trapezoidal of
!  order m * n, and Q is a permutation matrix of order n.

!  Row interchanges are used to ensure that the entries of L do not
!  exceed 1 in absolute value. Column interchanges are used to
!  ensure that the first r diagonal entries of U exceed PIVTOL in
!  absolute value. If r < m, the last (m-r) rows of U are zero.

!  This is the Level 2 BLAS version.

!  Arguments
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 1.
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 1.
!  A       (input/output) REAL(WP) array, dimension (LDA,N)
!          On entry, the m by n matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U*Q; the unit diagonal elements of L are not stored.
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!  PIVTOL  (input) REAL(WP)
!          The pivot tolerance. Any entry with absolute value less
!          than or equal to PIVTOL is regarded as unsuitable to be a
!          pivot.
!  IPIV    (output) INTEGER array, dimension (N)
!          The permutations; for 1 <= i <= RANK, row i of the
!          matrix was interchanged with row IPIV(i); for
!          RANK + 1 <= j <= N, column j of the
!          matrix was interchanged with column -IPIV(j).
!  RANK    (output) INTEGER
!          The computed rank of the matrix.
!  NNO     Number of columns (at end) from which pivots cannot be chosen


      REAL(WP), PARAMETER :: ZERO = 0.0_WP
      REAL(WP), PARAMETER :: ONE = 1.0_WP

      INTEGER I,JP,K,NCOL
      LOGICAL PIVOT
! I   Row index.
! JP  Pivot position.
! K   Main loop index.
! PIVOT True if there is a pivot in current column.

      INTEGER IDAMAX
      EXTERNAL IDAMAX

      EXTERNAL DGEMV,DSCAL,DSWAP

      INTRINSIC ABS

      NCOL = N - NNO
      RANK = 0

      DO 20 K = 1,N-NNO

         IF (RANK.LE.M-1) THEN
! Update diagonal and subdiagonal elements in column RANK+1.
            CALL DGEMV('No transpose',M-RANK,RANK,-ONE,A(RANK+1,1),LDA, &
                       A(1,RANK+1), &
                       1,ONE,A(RANK+1,RANK+1),1)
! Find pivot.
            JP = RANK + IDAMAX(M-RANK,A(RANK+1,RANK+1),1)
            PIVOT = ABS(A(JP,RANK+1)) .GT. PIVTOL
         END IF

         IF (PIVOT) THEN

! Update RANK
            RANK = RANK + 1

            IPIV(RANK) = JP

! Apply row interchange to columns 1:N.
            IF (JP.NE.RANK) CALL DSWAP(N,A(RANK,1),LDA,A(JP,1),LDA)

! Compute elements RANK+1:M of RANK-th column.
            IF (RANK.LT.M) &
              CALL DSCAL(M-RANK,ONE/A(RANK,RANK),A(RANK+1,RANK),1)

            IF (RANK.LT.N) THEN
! Compute row of U.
               CALL DGEMV('Transpose',RANK-1,N-RANK,-ONE, &
                          A(1,RANK+1),LDA,A(RANK,1), &
                          LDA,ONE,A(RANK,RANK+1),LDA)
            END IF

         ELSE

            DO 10 I = RANK+1,M
               A(I,RANK+1) = ZERO
   10       CONTINUE
! Apply column interchange and record it.
            IF (RANK+1.LT.NCOL) CALL DSWAP(M,A(1,RANK+1),1,A(1,NCOL),1)
            IPIV(NCOL) = -(RANK+1)
            NCOL = NCOL - 1

         END IF

         IF (RANK.GE.MIN(NCOL,M)) GO TO 40

   20 CONTINUE

! Update remaining block and set IPIV for the border columns
   40 IF (RANK.LT.M .AND. NNO.GT.0) &
         CALL DGEMM('No transpose','No transpose',M-RANK,NNO,RANK, &
                 -1.0D0,A(RANK+1,1),LDA,A(1,N-NNO+1),LDA, &
                 1.0D0,A(RANK+1,N-NNO+1),LDA)
      DO 50 I = N-NNO+1,N
        IPIV(I) = -I
   50 CONTINUE
! Set IPIV for remaining columns
      DO 60 I = RANK+1,NCOL
        IPIV(I) = -I
   60 CONTINUE

      END SUBROUTINE MP48_MA60FD

!********************************************************

      SUBROUTINE MP48_MA60GD(M,N,A,LDA,NB,PIVTOL,IPIV,RANK,NNO)

!  -- This is a variant of the LAPACK routine DGETRF --

      USE HSL_MP48_DATA_DOUBLE
      IMPLICIT NONE

      INTEGER LDA,M,N,NB,RANK,NNO
      REAL(WP) PIVTOL

      INTEGER IPIV(N)
      REAL(WP) A(LDA,N)

!  Purpose
!  =======
!  MP48_MA60GD computes an LU factorization of a general m-by-n
!  matrix A. The factorization has the form
!     A = P * L * U * Q
!  where P is a permutation matrix of order m, L is lower triangular
!  of order m with unit diagonal elements, U is upper trapezoidal of
!  order m * n, and Q is a permutation matrix of order n.
!  Row interchanges are used to ensure that the entries of L do not
!  exceed 1 in absolute value. Column interchanges are used to
!  ensure that the first r diagonal entries of U exceed PIVTOL in
!  absolute value. If r < m, the last (m-r) rows of U are zero.

!  This is the Level 3 BLAS version.

!  Arguments
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 1.
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 1.
!  A       (input/output) REAL(WP) array, dimension (LDA,N)
!          On entry, the m by n matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U*Q; the unit diagonal elements of L are not stored.
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!  NB      (input) INTEGER
!          The block size for BLAS3 processing.
!  PIVTOL  (input) REAL(WP)
!          The pivot tolerance. Any entry with absolute value less
!          than or equal to PIVTOL is regarded as unsuitable to be a
!          pivot.
!  IPIV    (output) INTEGER array, dimension (N)
!          The permutations; for 1 <= i <= RANK, row i of the
!          matrix was interchanged with row IPIV(i); for
!          RANK + 1 <= j <= N, column j of the
!          matrix was interchanged with column -IPIV(j).
!  RANK    (output) INTEGER
!          The computed rank of the matrix.
!  NNO     Number of columns (at end) from which pivots cannot be chosen


      REAL(WP), PARAMETER :: ZERO = 0.0_WP
      REAL(WP), PARAMETER :: ONE = 1.0_WP

      INTEGER I,JJ,JP,J1,J2,K,NCOL
      LOGICAL PIVOT
      REAL(WP) TEMP

! I   DO index for applying permutations.
! JJ  Column in which swaps occur.
! JP  Pivot position.
! J1  Column at start of current block.
! J2  Column at end of current block.
! K   Main loop index.
! PIVOT True if there is a pivot in current column.
! TEMP Temporary variable for swaps.

      EXTERNAL DGEMM,DGEMV,DSWAP,DSCAL,DTRSM,DTRSV

      INTEGER IDAMAX
      EXTERNAL IDAMAX

      INTRINSIC ABS,MIN

      NCOL = N - NNO
      RANK = 0

      J1 = 1
      J2 = MIN(NCOL,NB)

!     write (6,*) 'M,N,NCOL,NNO,RANK,J2,NB',M,N,NCOL,NNO,RANK,J2,NB

      DO 70 K = 1,N-NNO

!        write (6,*) 'Start of loop .. K =',K

         IF (RANK.LE.M-1) THEN

!          Update diagonal and subdiagonal elements in column RANK+1.
            CALL DGEMV('No transpose',M-RANK,RANK+1-J1,-ONE, &
                       A(RANK+1,J1),LDA, &
                       A(J1,RANK+1),1,ONE,A(RANK+1,RANK+1),1)

!          Find pivot.
            JP = RANK + IDAMAX(M-RANK,A(RANK+1,RANK+1),1)
            PIVOT = ABS(A(JP,RANK+1)) .GT. PIVTOL
!           write (6,*) 'jp,rank,pivot',jp,rank,pivot
         END IF

         IF (PIVOT) THEN
! Update RANK
            RANK = RANK+1
            IPIV(RANK) = JP

!           Apply row interchange to columns J1:J2
            IF (JP.NE.RANK) &
               CALL DSWAP(J2-J1+1,A(RANK,J1),LDA,A(JP,J1),LDA)

!           Compute elements RANK+1:M of RANK-th column.
            IF (RANK.LT.M) &
              CALL DSCAL(M-RANK,ONE/A(RANK,RANK),A(RANK+1,RANK),1)

            IF (RANK.LT.J2) THEN
!             Compute row of U within current block
               CALL DGEMV('Transpose',RANK-J1,J2-RANK,-ONE, &
                          A(J1,RANK+1),LDA, &
                          A(RANK,J1),LDA,ONE,A(RANK,RANK+1),LDA)
            END IF

         ELSE

            DO 10 I = RANK+1,M
               A(I,RANK+1) = ZERO
   10       CONTINUE

! Record column interchange and revise J2 if necessary
            IPIV(NCOL) = -(RANK+1)
!           write (6,*) 'NCOL',NCOL
! Apply column interchange.
            IF (RANK+1.LT.NCOL) &
                 CALL DSWAP(M,A(1,RANK+1),1,A(1,NCOL),1)
            IF (NCOL.GT.J2) THEN
!              write (6,*) 'J1,RANK',J1,RANK
! Perform row exchanges on new column
               DO 20 I = J1,RANK
                  JP = IPIV(I)
                  TEMP = A(I,RANK+1)
                  A(I,RANK+1) = A(JP,RANK+1)
                  A(JP,RANK+1) = TEMP
   20          CONTINUE
!              Update new column according to earlier pivots in block.
               IF(RANK+1.GT.J1) &
                  CALL DTRSV('Lower','No transpose','Unit', &
                  RANK+1-J1,A(J1,J1),LDA,A(J1,RANK+1),1)
            ELSE
               J2 = J2 - 1
!              write (6,*) 'Changing J2 to',J2
            END IF
            NCOL = NCOL - 1
!           write (6,*) 'Changing NCOL to',NCOL

         END IF

         IF (RANK.EQ.J2 .OR. &
             (RANK.GE.MIN(NCOL,M))) THEN
! Apply permutations to columns outside the block
! Columns from earlier blocks
!           write (6,*) 'RANK,J1,J2',RANK,J1,J2
            DO 40 JJ = 1,J1 - 1
               DO 30 I = J1,RANK
                  JP = IPIV(I)
                  TEMP = A(I,JJ)
                  A(I,JJ) = A(JP,JJ)
                  A(JP,JJ) = TEMP
   30          CONTINUE
   40       CONTINUE
!  Rest of still active columns
            DO 60 JJ = J2 + 1,NCOL
               DO 50 I = J1,RANK
                  JP = IPIV(I)
                  TEMP = A(I,JJ)
                  A(I,JJ) = A(JP,JJ)
                  A(JP,JJ) = TEMP
   50          CONTINUE
   60       CONTINUE
!  Border columns
            DO 65 JJ = N-NNO+1,N
               DO 55 I = J1,RANK
                  JP = IPIV(I)
                  TEMP = A(I,JJ)
                  A(I,JJ) = A(JP,JJ)
                  A(JP,JJ) = TEMP
   55          CONTINUE
   65       CONTINUE

            IF (K.NE.N) THEN
! Update the Schur complement
! Update future active columns
               IF (J2.LT.NCOL) &
               CALL DTRSM('Left','Lower','No transpose','Unit', &
                          RANK-J1+1, &
                          NCOL-J2,ONE,A(J1,J1),LDA,A(J1,J2+1),LDA)
! Update border columns
               IF (NNO.GT.0) &
               CALL DTRSM('Left','Lower','No transpose','Unit', &
                          RANK-J1+1, &
                          NNO,ONE,A(J1,J1),LDA,A(J1,N-NNO+1),LDA)
! Update Schur in future active columns
               IF (M.GT.J2) THEN
                  CALL DGEMM('No transpose','No transpose', &
                             M-J2,NCOL-J2,J2-J1+1,-ONE,A(J2+1,J1), &
                             LDA,A(J1,J2+1),LDA,ONE, &
                             A(J2+1,J2+1),LDA)
! Update Schur in border columns
               IF (NNO.GT.0) &
                  CALL DGEMM('No transpose','No transpose', &
                             M-J2,NNO,J2-J1+1,-ONE,A(J2+1,J1), &
                             LDA,A(J1,N-NNO+1),LDA,ONE, &
                             A(J2+1,N-NNO+1),LDA)
               END IF
            END IF

            J1 = J2 + 1
            J2 = MIN(J2+NB,NCOL)

         END IF

         IF (RANK.GE.MIN(NCOL,M)) GO TO 80

   70 CONTINUE

! Set IPIV for the border columns
   80 DO 90 I = N-NNO+1,N
        IPIV(I) = -I
   90 CONTINUE

! Set IPIV for remaining columns
      DO 95 I = RANK+1,NCOL
        IPIV(I) = -I
   95 CONTINUE

      END SUBROUTINE MP48_MA60GD

!********************************************************

      SUBROUTINE MP48_MA60ID(CNTL,ICNTL)
! Set default values for the control arrays.

      USE HSL_MP48_DATA_DOUBLE
      IMPLICIT NONE

      REAL(WP) CNTL(4)
      INTEGER ICNTL(10)

      CNTL(1) = 0.5D0
      CNTL(2) = 0.1D0
      CNTL(3) = 0.0D0
      CNTL(4) = 0.0D0
      ICNTL(1) = 6
      ICNTL(2) = 6
      ICNTL(3) = 1
      ICNTL(4) = 3
      ICNTL(5) = 32
      ICNTL(6) = 0
      ICNTL(7) = 0
      ICNTL(8) = 10

      END SUBROUTINE MP48_MA60ID

!********************************************************
      SUBROUTINE MP48_MA60CD(M,N,ICNTL,IQ,NP,LORU,LFACT,FACT,LIRNF, &
                             IRNF,IPTRL,IPTRU,B,X,W,INFO)

! MP48_MA60C/CD uses the LU factorization produced by
!     MP48_MA60B/BD to solve L x = b or  U x = b

      USE HSL_MP48_DATA_DOUBLE
      IMPLICIT NONE

      INTEGER M,N,ICNTL(7),IQ(N),NP,LORU
      INTEGER LFACT,LIRNF
      REAL(WP) FACT(LFACT)
      INTEGER IRNF(LIRNF),IPTRL(N),IPTRU(N)
      REAL(WP) B(MAX(M,N)),X(MAX(M,N)),W(MAX(M,N))
      INTEGER INFO(2)

! M  is an integer variable set to the number of rows.
!     It is not altered by the subroutine.
! N  is an integer variable set to the number of columns.
!     It is not altered by the subroutine.
! ICNTL must be set by the user as follows and is not altered.
!     ICNTL(1)  must be set to the stream number for error messages.
!       A value less than 1 suppresses output.
!     ICNTL(2) must be set to the stream number for diagnostic output.
!       A value less than 1 suppresses output.
!     ICNTL(3) must be set to control the amount of output:
!       0 None.
!       1 Error messages only.
!       2 Error and warning messages.
!       3 As 2, plus scalar parameters and a few entries of array
!         parameters on entry and exit.
!       4 As 2, plus all parameters on entry and exit.
!     ICNTL(5) must be set to control the level of BLAS used:
!       0 Level 1 BLAS.
!      >0 Level 2 BLAS.
! IQ is an integer array holding the permutation Q.
!     It is not altered by the subroutine.
! NP is an integer variable that must be unchanged since calling
!    MP48_MA60B/BD. It holds the number of rows and columns in packed
!    storage. It is not altered by the subroutine.
! LORU is an integer variable set to indicate whether we are solving for
!    a lower triangular (LORU=1) or an upper triangular (LORU \= 1)
!     matrix.
! LFACT is an integer variable set to the size of FACT.
!    It is not altered by the subroutine.
! FACT is an array that must be unchanged since calling MP48_MA60BD. It
!     holds the packed part of L/U by columns, and the full part of L/U
!     by columns. U has unit diagonal entries, which are not stored, and
!     the signs of the off-diagonal entries are inverted.  In the packed
!     part, the entries of U precede the entries of L; also the diagonal
!     entries of L head each column of L and are reciprocated.
!     FACT is not altered by the subroutine.
! LIRNF is an integer variable set to the size of IRNF.
!     It is not altered by the subroutine.
! IRNF is an integer array that must be unchanged since calling
!     MP48_MA60B/BD. It holds the row numbers of the packed part of
!     L/U, and the row numbers of the full part of L/U.
!     It is not altered by the subroutine.
! IPTRL is an integer array that must be unchanged since calling
!     MP48_MA60B/BD. For J = 1,..., NP, IPTRL(J) holds the position in
!     FACT and IRNF of the end of column J of L.
!     For J = NP+1,..., N, IPTRL(J) is equal to IPTRU(J).
!     It is not altered by the subroutine.
! IPTRU is an integer array that must be unchanged since calling
!     MP48_MA60B/BD. For J = 1,..., N, IPTRU(J) holds the position in
!     FACT and IRNF of the end of the packed part of column J of U.
!     It is not altered by the subroutine.
! B is an array that must be set to the vector b.
!     It is not altered.
! X is an array that need not be set on entry. On return, it holds the
!    solution x.
! W is a work array of length max(M,N).
! INFO need not be set on entry. On exit, it holds the following:
!    INFO(1) A nonzero value will indicate an error return. Possible
!      nonzero values are:
!      -1  M < 1 or N < 1
!      In this case INFO(2) will be set to one of the transgressing
!      values.

      REAL(WP), PARAMETER :: ZERO = 0.0_WP

      INTEGER I,II,IA1,IF1,J,LP,MF,MP,NF,NNO
      REAL(WP) PROD
! I Temporary variable holding row number.
! II Position of the current entry in IRNF.
! IA1 Position of the start of the current row or column.
! IF1 Position of the start of the full part of U.
! J Temporary variable holding column number.
! LP Unit for error messages.
! MF Number of rows held in full format.
! MP Unit for diagnostic messages.
! NF Number of columns held in full format.
! NNO holds the number of columns that are in the border.
! PROD Temporary variable used to accumulate inner products.

      LP = ICNTL(1)
      MP = ICNTL(2)
      IF (ICNTL(3).LE.0) LP = 0
      IF (ICNTL(3).LE.1) MP = 0

      NNO = ICNTL(6)
! Make some simple checks
      IF (M < 1 .OR. N < 1) GO TO 250

      IF (MP > 0 .AND. ICNTL(3) > 2) WRITE (MP, &
          '(/A/2(A,I6),A,I4,A,I4,A,I6)') &
          ' Entering Solve with','M =',M,' N =',N, &
          ' NP =',NP,' LORU =',LORU, &
          ' NNO =',NNO

      IF1 = IPTRL(N) + 1
      MF  = IRNF(2)
      NF  = N - NP

      IF (MP > 0 .AND. ICNTL(3) > 2) WRITE (MP, &
          '(A,I5,A,I5)') ' Size of full submatrix',MF,' by',NF
      IF (MP > 0 .AND. ICNTL(3) > 3) THEN
         DO 10 J = 1,N
            IF (J > 1) THEN
               IF (IPTRL(J-1) < IPTRU(J)) WRITE (MP, &
                   '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column',J,' of U', &
                    (FACT(II),IRNF(II),II=IPTRL(J-1)+1,IPTRU(J))
            END IF
            IF (IPTRU(J) < IPTRL(J)) WRITE (MP, &
                '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column',J,' of L', &
                (FACT(II),IRNF(II),II=IPTRU(J)+1,IPTRL(J))
   10    CONTINUE
         IF (MF > 0) THEN
           WRITE (MP,'(A)') ' Full part'
           WRITE (MP,'((6I12))') (IRNF(IF1+MF+J),J=0,NF-1)
           DO 20 I = 0,MF - 1
             WRITE (MP,'(I4,1P,6E12.4:/(4X,1P,6E12.4))') IRNF(IF1+I), &
               (FACT(IF1+I+J*MF),J=0,NF-1)
   20      CONTINUE
         END IF
      END IF

      IF (LORU == 1) THEN
! Solve for L
        IF (MP > 0 .AND. ICNTL(3) > 4) WRITE (MP, &
           '(A4,5ES12.4:/(4X,5ES12.4))') ' B =', (B(I),I=1,M)

! Forward substitution through the packed part of PL
        W(1:M) = B(1:M)
! We have to ensure all of X is defined
        X(1:M) = ZERO
        DO 150 I = 1,NP
          IA1 = IPTRU(I) + 1
          IF (IA1.LE.IPTRL(I)) THEN
             X(I) = W(IRNF(IA1))*FACT(IA1)
             IF (X(I) /= ZERO) THEN
                DO II = IA1 + 1,IPTRL(I)
                   W(IRNF(II)) = W(IRNF(II)) - FACT(II)*X(I)
                END DO
             END IF
             W(IRNF(IA1)) = X(I)
          END IF
  150   CONTINUE

        IF (MP > 0 .AND. ICNTL(3) > 4) WRITE (MP, &
           '(A4,5ES12.4:/(4X,5ES12.4))') ' W =', (W(I),I=1,M)
! Forward substitution through the full part of PL
        IF (MF > 0 .AND. NF > 0) THEN
          DO I = 1,MF
             W(I) = W(IRNF(IF1+I-1))
          END DO
          IF (MP > 0 .AND. ICNTL(3) > 4) WRITE (MP, &
             '(A4,5ES12.4:/(4X,5ES12.4))') ' W =', (W(I),I=1,MF)
          CALL MP48_MA60HD(LORU,MF,NF,FACT(IF1),MF,IRNF(IF1+MF),W, &
                           ICNTL(5),NNO)

          X(NP+1:NP+MF) = W(1:MF)
        ELSE
          X(NP+1:NP+MF) = ZERO
        END IF

      ELSE

! Solve for U.
! Back substitution through the full part of U
         DO I = 1,MAX(M,N)
           X(I) = B(I)
         END DO
         IF (MP > 0 .AND. ICNTL(3) > 3) &
            WRITE (MP,'(A/(4X,5ES12.4))') &
          ' Entering Solve with X =', (X(I),I=1,MAX(M,N))
         IF (MF > 0 .AND. NF > 0) THEN
         IF (MP > 0 .AND. ICNTL(3) > 4) WRITE (MP, &
             '(A4,5ES12.4:/(4X,5ES12.4))') ' X =', (X(I),I=1,MAX(M,N))
             CALL MP48_MA60HD(LORU,MF,NF,FACT(IF1),MF,IRNF(IF1+MF), &
                              X(NP+1),ICNTL(5),NNO)
          END IF
          IF (MP > 0 .AND. ICNTL(3) > 4) WRITE (MP, &
              '(A4,5ES12.4:/(4X,5ES12.4))') ' X       =', (X(I),I=1,M)
! Back substitution through the packed part of U
         DO 200 J = N,MAX(2,NP+1),-1
            PROD = X(J)
            DO II = IPTRL(J-1) + 1,IPTRU(J)
               X(IRNF(II)) = X(IRNF(II)) + FACT(II)*PROD
            END DO
  200    CONTINUE
         DO 220 J = NP,2,-1
            IA1 = IPTRU(J)
            IF (IA1.GE.IPTRL(J)) THEN
               X(J) = ZERO
            ELSE
               PROD = X(J)
               DO II = IPTRL(J-1) + 1,IA1
                  X(IRNF(II)) = X(IRNF(II)) + FACT(II)*PROD
               END DO
            END IF
  220    CONTINUE
         IF (NP.GE.1 .AND. IPTRU(1).GE.IPTRL(1)) X(1) = ZERO
         IF (IQ(1) > 0) THEN
! Permute X
            W(1:N) = X(1:N)
            DO I = 1,N
               X(IQ(I)) = W(I)
            END DO
         END IF
      END IF
      IF (MP > 0 .AND. ICNTL(3) > 3) THEN
         IF (LORU == 1) WRITE (MP, &
        '(A/(4X,5ES12.4))') ' Leaving Solve with X =', X(1:M)
         IF (LORU /= 1) WRITE (MP, &
        '(A/(4X,5ES12.4))') ' Leaving Solve with X =', X(1:N)
      END IF
      RETURN
! Error condition.
  250 INFO(1) = -1
      IF (LP > 0) WRITE (LP,'(/A/2(A,I8))') &
          ' **** Error return from Solve ****',' M =',M,' N =',N

      END SUBROUTINE MP48_MA60CD
!********************************************************
      SUBROUTINE MP48_MA60HD(LORU,M,N,A,LDA,IPIV,B,ICNTL5,NNO)

      USE HSL_MP48_DATA_DOUBLE
      IMPLICIT NONE

      INTEGER LORU,LDA,M,N,ICNTL5,NNO
      INTEGER IPIV(N)
      REAL(WP) A(LDA,N),B(*)

!  Solve a system of linear equations
!     L * X = B  or  U * X = B
!  with a triangular matrix A using the LU factorization computed
!  by MP48_MA60E/ED, MP48_MA60F/FD, or MP48_MA60G/GD.

!  LORU    (input) INTEGER
!          Specifies the form of the system of equations.
!          L * X = B  (LORU=1)
!          U * X = B  (LORU=2)
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 1.
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 1.
!  A       (input) REAL(WP) array, dimension (LDA,N)
!          The factors L and U from the factorization A = P*L*U
!          as computed by MP48_MA60G/GD.
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!  IPIV    (input) INTEGER array, dimension (N)
!          The permutations; for 1 <= i <= RANK, row i of the
!          matrix was interchanged with row IPIV(i); for
!          RANK + 1 <= j <= N, column j of the
!          matrix was interchanged with column -IPIV(j).
!  B       (input/output) REAL(WP) array, size max(M,N)
!          On entry, the right hand side vectors B for the system of
!          linear equations.
!          On exit, the solution vectors, X.
!  ICNTL5  (input) INTEGER
!          0 for BLAS1 or >0 for BLAS2
!  NNO     (input) INTEGER
!          Holds number of columns in border

      INTEGER I,J,K,NCOL,RANK
! I    Temporary variable.
! K    Temporary variable.
! NCOL Number of columns from which pivots can be chosen.
! RANK Rank of matrix.

      REAL(WP), PARAMETER :: ZERO = 0.0_WP
      REAL(WP), PARAMETER :: ONE = 1.0_WP

      REAL(WP) TEMP
      INTRINSIC MIN
      EXTERNAL DAXPY,DGEMV,DTRSV

      NCOL = N - NNO
! Find the rank
      RANK = 0
      DO RANK = MIN(M,NCOL),1,-1
         IF (IPIV(RANK) > 0) EXIT
      END DO

      IF (LORU == 1) THEN
! Solve L * X = B.
! Apply row interchanges to the right hand side.
         DO 30 I = 1,RANK
            K = IPIV(I)
            TEMP = B(I)
            B(I) = B(K)
            B(K) = TEMP
   30    CONTINUE

! Solve L*X = B, overwriting B with X.
         IF (ICNTL5 > 0) THEN
!           write (6,*) 'icntl5',icntl5
           IF (RANK > 0) &
              CALL DTRSV('L','NoTrans','Unit',RANK,A,LDA,B,1)
           IF (RANK > 0 .AND. RANK < LDA) &
              CALL DGEMV('NoTrans',M-RANK,RANK,-ONE, &
                          A(RANK+1,1),LDA,B(1),1,ONE,B(RANK+1),1)
         ELSE
             DO K = 1,RANK
                IF (B(K) /= ZERO .AND. M > K) &
                   CALL DAXPY(M-K,-B(K),A(K+1,K),1,B(K+1),1)
             END DO
         END IF
      ELSE

! Solve U*X = B, overwriting B with X.
         IF (ICNTL5 > 0) THEN
            IF (NCOL < N .AND. NNO > 0) &
               CALL DGEMV('NoTrans',RANK,NNO,-ONE, &
                           A(1,NCOL+1),LDA,B(NCOL+1),1,ONE,B,1)
            IF (RANK > 0) &
               CALL DTRSV('U','NoTrans','NonUnit',RANK,A,LDA,B,1)
         ELSE
            DO 50 K = RANK,1,-1
              DO J = NCOL+1,N
                B(K) = B(K) - B(J)*A(K,J)
              END DO
              B(K) = B(K)/A(K,K)
              IF (K > 1 .AND. B(K) /= ZERO) &
                 CALL DAXPY(K-1,-B(K),A(1,K),1,B(1),1)
   50       CONTINUE
        END IF

! Apply column interchanges to the right hand side.
        DO 70 I = RANK+1,N
           K = -IPIV(I)
           TEMP = B(I)
           B(I) = B(K)
           B(K) = TEMP
   70   CONTINUE
      END IF

      END SUBROUTINE MP48_MA60HD

