module constants
  implicit none
  INTEGER, PARAMETER :: DPC = KIND(0.0D0)
  !INTEGER, PARAMETER :: DPC = KIND(0.0)
  INTEGER, PARAMETER :: FSORD = 1
  !switch to single by setting FSORD = 0, DPC and hsl_mp48d.f90 accordingly
end module constants

SUBROUTINE SPEC48_SINGLE(indata,irn1,jcn1,b1,values1,x,neleperrow,ai1,fcomm)
  USE HSL_MP48_DOUBLE
  use HSL_mc66_double
  use constants
  IMPLICIT NONE

  !INTEGER, PARAMETER :: myreal = KIND( 1.0D+0 )!SELECTED_REAL_KIND(8)!kind(real*8)!KIND( 0.0D0 )
  !INTEGER, PARAMETER :: myint = KIND( 1)!selected_int_kind(8)!kind(interger*4)!KIND( 1)
  !INTEGER, PARAMETER :: myint64 = selected_int_kind(16)!kind(interger*4)!KIND( 1)
  !type  datawrp
  !integer (8) :: nz
  !integer (8) :: m
  !integer (8) :: mpisize
  !integer (8) :: nblock
  !integer (8) :: nsbbdblocks
  !integer , pointer :: irn(:)
  !integer , pointer :: jcn(:)
  !real , pointer :: b(:)
  !real , pointer :: values(:)
  !end type datawrp
  !type(datawrp) indata
  integer (8) indata(*)
  integer (4) fcomm(*),fcomm1
  integer(4) irn1(*),jcn1(*),ai1(*),neleperrow(*)
  real (kind=DPC) b1(*),values1(*),x(*)
  ! mc66 controller:
  type (mc66_control) :: control

  ! random number seed
  type (fa14_seed) :: seed

  !nz: number of nonzeros
  !n: number of rows.
  !m: number of columns
  !irn: row indices of the matrix
  !jcn: column indices of the matrix
  !integer (kind=myint) :: nz
  !integer (kind=myint) :: n,m
  integer nz,n,m
  !integer, pointer :: irn(:)
  !integer, pointer :: jcn(:)
  !integer (kind=myint), pointer :: neqsb(:)

  ! rowptr: rowptr[1:nblocks] is the starting row index for
  !         each diagonal block. rowptr[nblocks+1]=m
  ! colptr: colptr[1:nblocks] is the starting column index
  !         for each diagonal block. colptr[nblocks+1]
  !         is the starting column index of the border.
  !integer (kind = myint), allocatable, dimension (:) :: rowptr
  !integer (kind = myint), allocatable, dimension (:) :: colptr
  integer, pointer :: rowptr(:)
  integer, pointer :: colptr(:)

  ! column order: column_order(i) is the original column index
  !   of the i-column of the reordered matrix
  ! row ordering: row_order(i) is the original row index
  !   of the i-row of the reordered matrix
  ! row_position: row_position(i) is the new row index
  !   of the original row i.
  ! column_position: column_position(i) is the new column index
  !   of the original column i.
  !integer (kind = myint), allocatable, dimension (:) :: row_order,&
  !  column_order,column_position !,row_position
  integer, pointer :: row_order(:),column_order(:),column_position(:) !,row_position
  ! nblocks: number of blocks
  ! netcut: column size of size of border
  ! info: info tag
  ! kblocks: the actual number of blocks in the SBBD form
  ! rowdiff: row dimension imbalance in percentage term
  integer nblocks,netcut,info,kblocks
  !integer (kind = myint) :: nblocks
  !integer (kind = myint) :: netcut
  !integer (kind = myint) :: info
  !integer (kind = myint) :: kblocks
  !real (kind = myreal) :: rowdiff
  DOUBLE PRECISION rowdiff
  integer i,j,k,h,o,p
  TYPE (MP48_DATA) data
  INTEGER ERCODE,ST
  LOGICAL FLAG
  integer myid, numprocs
  fcomm1=fcomm(1)

  ! Program to illustrate use of MP48.
!  CALL MPI_INIT(ERCODE)
  ! Define a communicator for the package
  data%COMM = fcomm1 !MPI_COMM_WORLD
  ! Initialize package
  data%JOB = 1
  CALL MP48AD(data)
  ! Reset control parameters (if required)
  ! Read all values on host
  data%ICNTL(7) = 3
  IF (data%RANK.EQ.0) THEN

    !write(*,"(a,I10)")       "myint =                  ",HUGE(netcut)  !nz = 10
    !write(*,"(a,I10)")       "myid =                  ",HUGE(myid)  !nz = 10
    !m = 4; n = 4; nblocks = 2
    !allocate(irn(nz),jcn(nz))
    !irn = (/4,4,4,2,2,2,1,1,3,3/)
    !jcn = (/4,3,1,4,3,1,2,1,2,1/)
    nz = indata(1)!indata%nz
    m = indata(2)!indata%m
    !allocate(irn1(nz))
    !do i = 1, nz
       !write(*,"(i10,i10)") irn1(i),jcn1(i)
    !end do
    n = m
    IF (indata(5).GE.2) then !%nsbbdblocks
      nblocks = indata(5)!%nsbbdblocks
    ELSE
      nblocks = 2
    ENDIF
    !allocate(irn(nz),jcn(nz))
    !do i = 1, nz
      !irn(i) = irn1(i)!indata%irn(i)
      !jcn(i) = jcn1(i)!indata%jcn(i)
    !end do
    !write(*,"(a)") "the original matrix is "
    !do i = 1, nz
    !   write(*,"(i4,i4)") irn(i),jcn(i)
    !end do
    !write(*,"(a,I10)")       "myint =                  ",irn(nz)  !nz = 10
    !write(*,"(a,I10)")       "myint =                  ",jcn(nz)  !nz = 10
    allocate(row_order(m),rowptr(nblocks+1), &
      column_order(n), colptr(nblocks+1))
    !write(*,"(a)") " "
    !write(*,"(a)") "generating the ordering ===== "
    control%COARSEN_SCHEME=2
    !control%COARSEST_SIZE=50
    !control%GRID_RDC_FAC=0.6
    write(*,"(i10,i10,i10,i10)") nz,m,n,nblocks
    call mc66(m,n,nz,irn1,jcn1,nblocks,control,seed, &
      row_order,info,rowptr,column_order,colptr,&
      netcut,rowdiff,kblocks)
    write(*,"(a)") "generating the ordering =====1 "
    if (info /= 0) then
      call mc66_print_message(info)
      if (info < 0) stop "mc66 failed"
    end if
    write(*,"(a,I10)")       "netcut =                  ",netcut
    write(*,"(a,f10.2,'%')") "row dimension imbalance = ", rowdiff
    ! dump block sizes
    !allocate(neqsb(nblocks))
    !do i = 1, nblocks
      !neqsb(i)=rowptr(i+1)-rowptr(i)
      !write(*,"('block ',i4,' of dimension  ',i10,' X ',i10)") &
      !i,rowptr(i+1)-rowptr(i),colptr(i+1)-colptr(i)
    !end do
    ! reorder the original matrix
    allocate(column_position(n))!, row_position(m))
    do i = 1, n
      column_position(column_order(i)) = i
    end do
    !goto 87
    !do i = 1, m
    !  row_position(row_order(i)) = i
    !end do
    !irn = row_position(irn)
    !jcn = column_position(jcn)
    !write(*,"(a)") " "
    !write(*,"(a)") "the reordered matrix is "
    !do i = 1, nz
    !  write(*,"(i5,i5)") irn(i),row_order(irn(i))!,jcn(i),row_order(irn(i)),column_order(jcn(i))
    !end do
    !write(*,"(a)") "the new reordered matrix is "


    !OPEN (UNIT=50,FILE='hsl_mp48ss.data')
    !READ (50,*) data%NEQ,data%NBLOCK,data%NE
    data%NEQ=m !indata%x
    data%NBLOCK=nblocks !indata%y
    data%NE=nz !indata%z
    ! Allocate arrays for matrix data
    ALLOCATE(data%NEQSB(1:data%NBLOCK),STAT=ST )
    ALLOCATE(data%EQPTR(1:data%NEQ+1),STAT=ST )
    ALLOCATE(data%EQVAR(1:data%NE),STAT=ST )
    ALLOCATE(data%VALUES(1:data%NE),STAT=ST )
    ALLOCATE(data%B(1:data%NEQ),STAT=ST )
    ! Read matrix data on host.
    !READ (50,*) data%NEQSB(1:data%NBLOCK)
    do i = 1, nblocks
      data%NEQSB(i)=rowptr(i+1)-rowptr(i)!neqsb(i)
      write(*,"('block ',i4,' of dimension  ',i10,' X ',i10)") &
      i,rowptr(i+1)-rowptr(i),colptr(i+1)-colptr(i)
    end do
    !READ (50,*) data%EQPTR(1:data%NEQ+1)
    h=1
    do j = 1, m
      data%EQPTR(j)=h
      do i = 1,neleperrow(row_order(j))
        data%EQVAR(h)=column_position(jcn1(ai1(row_order(j))+i-1))
        data%VALUES(h)=values1(ai1(row_order(j))+i-1)
        h=h+1
      end do
    end do
    !write(*,"(a)") "f"
    data%EQPTR(m+1)=nz+1
    !write(*,"(a)") "f"
    !READ (50,*) data%EQVAR(1:data%NE)
    !READ (50,*) data%VALUES(1:data%NE)
    ! Also read right hand side
    !READ (50,*) data%B(1:data%NEQ)
    do i = 1, m
      data%B(i)=b1(row_order(i))
    end do
    !write(*,"(a)") "f"
  END IF
  !do i = 1, nz
    !write(*,"(i5,f10.2)") data%EQVAR(i),data%VALUES(i)
  !end do
  !write(*,"(a)") "h"
  !write(*,"(a)") "the final matrix is "
  !do i = 1, nz
  !  write(*,"(i4,F10.2)") data%EQVAR(i),data%VALUES(i)
  !end do
  !write(*,"(a)") "the final B "
  !do i = 1, m
  !  write(*,"(i4,F10.2,F10.2,i4)") row_order(i),data%B(i),b1(i),data%EQPTR(i)
  !end do
  !write(*,"(i4)") data%EQPTR(m+1)
  CALL MPI_BARRIER(data%COMM,ERCODE)
  call MPI_COMM_RANK( fcomm1 , myid, ERCODE ) !MPI_COMM_WORLD
  call MPI_COMM_SIZE( fcomm1, numprocs, ERCODE ) !MPI_COMM_WORLD
  print *, "Process ", myid, " of ", numprocs, " is alive1"
  ! Call MP48AD/AD (combine analyse/factorize/solve)
  data%JOB = 25
  CALL MP48AD(data)
  IF (data%RANK.EQ.0) THEN
    IF (data%ERROR.LT.0) THEN
      WRITE (6,*) ' Unexpected error return'
    ELSE
      !WRITE (6,'(//A/,6ES11.3)') &
      !' The solution is: ',data%X(1:data%NEQ)
      do o=1,m
        x(column_order(o))=data%X(o)
      end do
    END IF
    deallocate(row_order,rowptr,column_order, colptr)
  END IF
  call MPI_COMM_RANK( fcomm1, myid, ERCODE ) !MPI_COMM_WORLD
  call MPI_COMM_SIZE( fcomm1, numprocs, ERCODE ) !MPI_COMM_WORLD
  print *, "Process ", myid, " of ", numprocs, " is alive2"
  data%JOB = 6
  CALL MP48AD(data)
  !deallocate(row_order)!,column_order,column_position)
!  CALL MPI_FINALIZE(ERCODE)
  !STOP
  !deallocate(irn,jcn,row_order,rowptr, &
  !  column_order, colptr,neqsb,column_position, &
  !  data%NEQSB,data%EQPTR,data%EQVAR,data%VALUES,data%B)
  !87 continue
END SUBROUTINE SPEC48_SINGLE

SUBROUTINE SPEC48_NOMC66(indata,jcn1,b1,values1,x,neleperrow,fcomm,rowptrin,colptrin)
  USE HSL_MP48_DOUBLE
  use constants
  IMPLICIT NONE

  !INTEGER, PARAMETER :: myreal = SELECTED_REAL_KIND(8)!kind(real*8)!KIND( 0.0D0 )
  !INTEGER, PARAMETER :: myint = selected_int_kind(8)!kind(interger*4)!KIND( 1)
  !INTEGER, PARAMETER :: myint64 = selected_int_kind(16)!kind(interger*4)!KIND( 1)
!  type  datawrp
!  integer (8) :: nz
!  integer (8) :: m
!  integer (8) :: mpisize
!  integer (8) :: nblock
!  integer (8) :: nsbbdblocks
!  end type datawrp
!  type(datawrp) indata
  integer myid, numprocs
  integer (8) :: i,j,o,h
  integer (8) :: nblocks,maxsbcols,ncols

  TYPE (MP48_DATA) data
  INTEGER ERCODE,ST
  LOGICAL FLAG
  integer (4) fcomm(*),fcomm1
  integer(4) neleperrow(*),jcn1(*)
  integer (8) indata(*),rowptrin(*),colptrin(*)
  real (kind=DPC) b1(*),values1(*),x(*)
  integer (8) :: nz
  integer (8) :: m
  fcomm1=fcomm(1)
  !integer (kind=myint64), pointer :: neqsb(:)

  ! rowptr: rowptr[1:nblocks] is the starting row index for
  !         each diagonal block. rowptr[nblocks+1]=m
  ! colptr: colptr[1:nblocks] is the starting column index
  !         for each diagonal block. colptr[nblocks+1]
  !         is the starting column index of the border.

  data%COMM = fcomm1 !MPI_COMM_WORLD
  !print *, "fcomm ", myid
  ! Initialize package
  data%JOB = 1
  CALL MP48AD(data)
  ! Reset control parameters (if required)
  ! Read all values on host
  data%ICNTL(7) = 3
  !data%ICNTL(11) = 1
  IF (data%RANK.EQ.0) THEN
    nz = indata(1)!%nz
    m = indata(2)!%m
    !write(*,"(i4,i4)") indata%nz,indata%m
    !allocate(irn1(nz))
    !do i = 1, nz
    !   write(*,"(i4,i4)") irn1(i),jcn1(i)
    !end do
    !n = m
    nblocks = indata(4)!%nblock
    !allocate(neqsb(nblocks))
    !do i = 1, nblocks
      !neqsb(i)=rowptrin(i+1)-rowptrin(i)
      !write(*,"('block ',i4,' of dimension  ',i10,' X ',i10)") &
      !i,neqsb(i),colptrin(i+1)-colptrin(i)
    !end do
    !OPEN (UNIT=50,FILE='hsl_mp48ss.data')
    !READ (50,*) data%NEQ,data%NBLOCK,data%NE
    data%NEQ=indata(2)!%m !indata%x
    data%NBLOCK=nblocks !indata%y
    data%NE=indata(1)!%nz !indata%z
    ! Allocate arrays for matrix data
    ALLOCATE(data%NEQSB(1:data%NBLOCK),STAT=ST )
    ALLOCATE(data%EQPTR(1:data%NEQ+1),STAT=ST )
    ALLOCATE(data%EQVAR(1:data%NE),STAT=ST )
    ALLOCATE(data%VALUES(1:data%NE),STAT=ST )
    ALLOCATE(data%B(1:data%NEQ),STAT=ST )
    ! Read matrix data on host.
    !READ (50,*) data%NEQSB(1:data%NBLOCK)
    maxsbcols=0
    do i = 1, nblocks
      data%NEQSB(i)=rowptrin(i+1)-rowptrin(i)!neqsb(i)
      write(*,"('nomc block ',i4,' of dimension  ',i10,' X ',i10)") &
      i,data%NEQSB(i),colptrin(i+1)-colptrin(i)
      ncols=colptrin(i+1)-colptrin(i)
      if(maxsbcols.LT.ncols)maxsbcols=ncols
    end do
    !print *, "row ", m,"nz ",nz,"maxcolsb ",m-colptrin(nblocks+1)
    maxsbcols=maxsbcols+m-colptrin(nblocks+1)
    !READ (50,*) data%EQPTR(1:data%NEQ+1)
    print *, "row ", m,"nz ",nz,"maxcolsb ",maxsbcols
    h=1
    data%MAXSBCOLS=maxsbcols
    do j = 1, m
      data%EQPTR(j)=h
      !if(neleperrow(j).LT.1) print *, "j",j,"row123 ", neleperrow(j), " of ", h
      do i = 1,neleperrow(j)
        !if(j.EQ.7071819)print *, "m ", j, "eq ", jcn1(h)
        data%EQVAR(h)=jcn1(h)
        !if(j.EQ.108494342)print *, "h ", h, " of ", jcn1(h)
        data%VALUES(h)=values1(h)
        h=h+1
      end do
    end do
    !write(*,"(a)") "f"
    data%EQPTR(m+1)=h
    !write(*,"(a)") "f"
    !READ (50,*) data%EQVAR(1:data%NE)
    !READ (50,*) data%VALUES(1:data%NE)
    ! Also read right hand side
    !READ (50,*) data%B(1:data%NEQ)
    do i = 1, m
      data%B(i)=b1(i)
    end do
    !write(*,"(a)") "f"
  END IF
  !do i = 1, nz
    !write(*,"(i5,f10.2)") data%EQVAR(i),data%VALUES(i)
  !end do
  !write(*,"(a)") "h"
  !write(*,"(a)") "the final matrix is "
  !do i = 1, nz
  !  write(*,"(i4,F10.2)") data%EQVAR(i),data%VALUES(i)
  !end do
  !write(*,"(a)") "the final B "
  !do i = 1, m
  !  write(*,"(i4,F10.2,F10.2,i4)") row_order(i),data%B(i),b1(i),data%EQPTR(i)
  !end do
  !write(*,"(i4)") data%EQPTR(m+1)
  CALL MPI_BARRIER(data%COMM,ERCODE)
  call MPI_COMM_RANK( fcomm1 , myid, ERCODE ) !MPI_COMM_WORLD
  call MPI_COMM_SIZE( fcomm1, numprocs, ERCODE ) !MPI_COMM_WORLD
  print *, "Process ", myid, " of ", numprocs, " is alive1"
  ! Call MP48AD/AD (combine analyse/factorize/solve)
  data%JOB = 25
  CALL MP48AD(data)
  print *, "Processa ", myid, " of ", numprocs, " is alive1"
  IF (data%RANK.EQ.0) THEN
    IF (data%ERROR.LT.0) THEN
      WRITE (6,*) ' Unexpected error return'
    ELSE
      !WRITE (6,'(//A/,6ES11.3)') &
      !' The solution is: ',data%X(1:data%NEQ)
      do o=1,m
        x(o)=data%X(o)
      end do
    END IF
  END IF
  !call MPI_COMM_RANK( fcomm1, myid, ERCODE ) !MPI_COMM_WORLD
  !call MPI_COMM_SIZE( fcomm1, numprocs, ERCODE ) !MPI_COMM_WORLD
  print *, "Process ", myid, " of ", numprocs, " is alive2"
  data%JOB = 6
  CALL MP48AD(data)
END SUBROUTINE SPEC48_NOMC66

SUBROUTINE SPEC51_RANK(INSIZE,CNTL6,IRN,JCN,VA)
  use constants
  IMPLICIT NONE

  !INTEGER, PARAMETER :: myreal = SELECTED_REAL_KIND(8)!kind(real*8)!KIND( 0.0D0 )
  !INTEGER, PARAMETER :: myint = selected_int_kind(16)!kind(interger*4)!KIND( 1)
  integer M,N,NE
  integer (4) INSIZE(*)
  integer (4) JCN(*),IRN(*)
  real (kind=DPC) VA(*),CNTL6(*)
  integer I,LA, MAXN,RANK1,SGNDET,T,NEFAC
  real (kind=DPC) LOGDET
  real (kind=DPC), pointer :: CNTL(:),RINFO(:),W(:)!A(:),
  integer, pointer :: COLS(:),ICNTL(:),INFO(:),IW(:),KEEP(:),ROWS(:)!JCN1(:),IRN1(:),
  M=INSIZE(1)
  N=INSIZE(2)
  NE=INSIZE(3)
  NEFAC=INSIZE(5)
  !write(*,"(A,i5,A,i5,A,i5,A,i5,A,f10.6)") 'M ',M,'NE ',NE,'DPC',DPC,'LA',LA,'A',VA(2)
  !write(*,"(A,i5,A,i5)") 'IRN ',IRN(2),'JCN ',JCN(2)
  !write(*,"(A,i5,A,i5)") 'IRN ',IRN(3),'JCN ',JCN(3)
  !write(*,"(A,i5,A,i5)") 'IRN ',IRN(4),'JCN ',JCN(4)
  if(NEFAC.EQ.0) then
    LA=2*NE
  else
    LA=ceiling((NEFAC/100.0)*NE)
  endif
  MAXN=N
  IF (N.LT.M) THEN
    MAXN=M
  END IF
  allocate(CNTL(10),RINFO(10),W(5*MAXN))!A(LA),
    !WRITE(6,'(/,A)')' Determinant is zero'

!  real (kind=myreal) A(LA),CNTL(10),LOGDET,RINFO(10),W(MAXN4)

  allocate(COLS(N),ICNTL(20),INFO(20),IW(6*M+3*N),ROWS(M))!JCN1(LA),IRN1(LA),
      !do I=1,NE
        !IRN1(I)=IRN(I)
        !JCN1(I)=JCN(I)
        !A(I)=VA(I)
      !end do

!     Factorize matrix
    IF (FSORD.EQ.1) THEN
      CALL MA48ID(CNTL,ICNTL)
    ELSE
      CALL MA48I(CNTL,ICNTL)
    ENDIF
    !WRITE(6,'(/,A,i5,F9.3)')' Determinant is zero2',FSORD,CNTL6(1)
    IF(CNTL6(1).EQ.0)THEN
      IF (FSORD.EQ.1) THEN
        CNTL(4)=1e-4
      else
        CNTL(4)=0.3
      ENDIF
    ELSE
      CNTL(4)=CNTL6(1)!1e-4!0.0000000001
    ENDIF
    !WRITE(6,'(/,A)')' Determinant is zero'
    T=M+5*N+4*N/ICNTL(6)+7
    allocate(KEEP(T))
    !WRITE(6,'(/,A)')' Determinant is zero1'
    ICNTL(7)=0
    IF (FSORD.EQ.1) THEN
      CALL MA48AD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
      !write(*,"(A,i5,A,i5,A,i5)") ' INFO =',INFO(1),'INFO10',INFO(10),'ct7',ICNTL(7)
      CALL MA48BD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,RINFO)
      !write(*,"(A,i5)") ' INFO1 =',INFO(1)

!     Compute the determinant
      CALL MA51CD(M,N,LA,VA,IRN,KEEP,SGNDET,LOGDET,W)
    ELSE
      CALL MA48A(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
      CALL MA48B(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,RINFO)
      CALL MA51C(M,N,LA,VA,IRN,KEEP,SGNDET,LOGDET,W)
    ENDIF

!     Print result
    IF(SGNDET.GT.0)THEN
      WRITE(6,'(/,A,F9.3)')&
          ' Determinant is positive; log(determinant) =',LOGDET
          INSIZE(4)=M
    ELSE IF(SGNDET.LT.0)THEN
        WRITE(6,'(/,A,F9.3)')&
          ' Determinant is negative; log(-determinant) =',LOGDET
          INSIZE(4)=M
      ELSE
        WRITE(6,'(/,A)')' Determinant is zero'
!     Determine the nonsingular submatrix of the factorization
        IF (FSORD.EQ.1) THEN
          CALL MA51AD(M,N,LA,IRN,KEEP,RANK1,ROWS,COLS,W)
        ELSE
          CALL MA51A(M,N,LA,IRN,KEEP,RANK1,ROWS,COLS,W)
        ENDIF
!     Print out rank and identifying vectors.
        INSIZE(4)=RANK1
        !write(*,"(A,i10,i10,i10,i10)") ' Rank of matrix Asub =',RANK1,M,N,COLS(3)
        !write(*,"(A,i10,i10,i10,i10)") ' Rank of matrix Asub =',RANK1,M,N,COLS(4)
        !do I=RANK+1,N
          !write(*,"(A,i5,A,i5)") ' ROWS =',ROWS(I),' COLS =',COLS(I)
          !write(*,"(A,i5,A,i5)") ' COLS =',COLS(I)
        !end do
        do I=1,M
          IRN(I)=ROWS(I)
          !write(*,"(A,i5)") ' IRN =',IRN(I)
        end do
        do I=1,N
          JCN(I)=COLS(I)
          !write(*,"(A,i5,i5)") ' JCN =',I,JCN(I)
        end do
    END IF
  deallocate(CNTL,RINFO,W)!A,
  deallocate(COLS,ICNTL,INFO,IW,KEEP,ROWS)!IRN1,JCN1,
END SUBROUTINE SPEC51_RANK

SUBROUTINE SPEC51M_RANK(INSIZE,CNTL6,IRN,JCN,VA,IRNA,JCNA,KEEP,W,IW)
  use constants
  IMPLICIT NONE

  !INTEGER, PARAMETER :: myreal = SELECTED_REAL_KIND(8)!kind(real*8)!KIND( 0.0D0 )
  !INTEGER, PARAMETER :: myint = selected_int_kind(16)!kind(interger*4)!KIND( 1)
  integer M,N,NE
  integer (4) INSIZE(*)
  integer (4) JCN(*),IRN(*),JCNA(*),IRNA(*),KEEP(*),IW(*)
  real (kind=DPC) VA(*),CNTL6(*),W(*)
  integer I,LA, MAXN,RANK1,SGNDET,T,NEFAC
  real (kind=DPC) LOGDET
  real (kind=DPC), pointer :: CNTL(:),RINFO(:)!,W(:),A(:),
  integer, pointer :: ICNTL(:),INFO(:)!,IW(:),KEEP(:),JCN1(:),IRN1(:),COLS(:),ROWS(:)
  M=INSIZE(1)
  N=INSIZE(2)
  NE=INSIZE(3)
  NEFAC=INSIZE(5)
  !write(*,"(A,i5,A,i5,A,i5,A,i5,A,f10.6)") 'M ',M,'NE ',NE,'DPC',DPC,'LA',LA,'A',VA(2)
  !write(*,"(A,i5,A,i5)") 'IRN ',IRN(2),'JCN ',JCN(2)
  !write(*,"(A,i5,A,i5)") 'IRN ',IRN(3),'JCN ',JCN(3)
  !write(*,"(A,i5,A,i5)") 'IRN ',IRN(4),'JCN ',JCN(4)
  if(NEFAC.EQ.0) then
    LA=2*NE
  else
    LA=ceiling((NEFAC/100.0)*NE)
  endif
  if(LA.NE.INSIZE(6)) then
    !write(*,"(A,i15,i15)") 'Note LA!!!!',LA,INSIZE(6)
    LA=INSIZE(6)
  end if
  MAXN=N
  IF (N.LT.M) THEN
    MAXN=M
  END IF
  allocate(CNTL(10),RINFO(10))!A(LA),,W(5*MAXN)
    !WRITE(6,'(/,A)')' Determinant is zero'

!  real (kind=myreal) A(LA),CNTL(10),LOGDET,RINFO(10),W(MAXN4)

  allocate(ICNTL(20),INFO(20))!JCN1(LA),IRN1(LA),,ROWS(M)COLS(N),,IW(6*M+3*N)
      !do I=1,NE
        !IRN1(I)=IRN(I)
        !JCN1(I)=JCN(I)
        !A(I)=VA(I)
      !end do

!     Factorize matrix
    IF (FSORD.EQ.1) THEN
      CALL MA48ID(CNTL,ICNTL)
    ELSE
      CALL MA48I(CNTL,ICNTL)
    ENDIF
    !WRITE(6,'(/,A,i5,F9.3)')' Determinant is zero2',FSORD,CNTL6(1)
    IF(CNTL6(1).EQ.0)THEN
      IF (FSORD.EQ.1) THEN
        CNTL(4)=1e-4
      else
        CNTL(4)=0.3
      ENDIF
    ELSE
      CNTL(4)=CNTL6(1)!1e-4!0.0000000001
    ENDIF
    !WRITE(6,'(/,A)')' Determinant is zero'
    !T=M+5*N+4*N/ICNTL(6)+7
    !allocate(KEEP(T))
    !WRITE(6,'(/,A)')' Determinant is zero1'
    ICNTL(7)=0
    IF (FSORD.EQ.1) THEN
      CALL MA48AD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
      !write(*,"(A,i5,A,i5,A,i5)") ' INFO =',INFO(1),'INFO10',INFO(10),'ct7',ICNTL(7)
      CALL MA48BD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,RINFO)
      !write(*,"(A,i5)") ' INFO1 =',INFO(1)

!     Compute the determinant
      CALL MA51CD(M,N,LA,VA,IRN,KEEP,SGNDET,LOGDET,W)
    ELSE
      CALL MA48A(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
      CALL MA48B(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,RINFO)
      CALL MA51C(M,N,LA,VA,IRN,KEEP,SGNDET,LOGDET,W)
    ENDIF

!     Print result
    IF(SGNDET.GT.0)THEN
      WRITE(6,'(/,A,F9.3)')&
          ' Determinant is positive; log(determinant) =',LOGDET
          INSIZE(4)=M
    ELSE IF(SGNDET.LT.0)THEN
        WRITE(6,'(/,A,F9.3)')&
          ' Determinant is negative; log(-determinant) =',LOGDET
          INSIZE(4)=M
      ELSE
        WRITE(6,'(/,A)')' Determinant is zero'
!     Determine the nonsingular submatrix of the factorization
        IF (FSORD.EQ.1) THEN
          CALL MA51AD(M,N,LA,IRN,KEEP,RANK1,IRNA,JCNA,W)!ROWS,COLS
        ELSE
          CALL MA51A(M,N,LA,IRN,KEEP,RANK1,IRNA,JCNA,W)!,ROWS,COLS
        ENDIF
!     Print out rank and identifying vectors.
        INSIZE(4)=RANK1
        !write(*,"(A,i10,i10,i10,i10)") ' Rank of matrix Asub =',RANK1,M,N,COLS(3)
        !write(*,"(A,i10,i10,i10,i10)") ' Rank of matrix Asub =',RANK1,M,N,COLS(4)
        !do I=RANK+1,N
          !write(*,"(A,i5,A,i5)") ' ROWS =',ROWS(I),' COLS =',COLS(I)
          !write(*,"(A,i5,A,i5)") ' COLS =',COLS(I)
        !end do
        !do I=1,M
          !IRN(I)=ROWS(I)
          !write(*,"(A,i5)") ' IRN =',IRN(I)
        !end do
        !do I=1,N
          !JCN(I)=COLS(I)
          !write(*,"(A,i5,i5)") ' JCN =',I,JCN(I)
        !end do
    END IF
  deallocate(CNTL,RINFO)!A,,W
  deallocate(ICNTL,INFO)!IRN1,JCN1,,ROWSCOLS,,KEEP,IW
END SUBROUTINE SPEC51M_RANK

!SUBROUTINE SPEC48_MSOL2(INSIZE,IRN,JCN,VA,B,X,IRNC,JCNC,VAC,IRNV,JCNV,VAV)
!  use constants
!  IMPLICIT NONE
!
!  INTEGER, PARAMETER :: myreal = SELECTED_REAL_KIND(8)!kind(real*8)!KIND( 0.0D0 )
!  INTEGER, PARAMETER :: myint = selected_int_kind(16)!kind(interger*4)!KIND( 1)
!  integer NEFAC,JOB
!  integer (kind=myint) INSIZE(*),JCN(*),IRN(*),JCNC(*),IRNC(*),JCNV(*),IRNV(*)
!  integer M,N,NE,T,NC,NV,MC,MV,NEC,NEV
!  real (kind=DPC) VA(*),VAC(*),VAV(*),B(*),X(*)
!  integer I,J,L,LA, MAXN,IFAIL,LP,II,SGNDET
!  LOGICAL TRANS,checksol
!  !DOUBLE PRECISION, pointer :: A(:),CNTL(:),RINFO(:),W(:),ERROR(:),RHS(:),SOL(:)
!  DOUBLE PRECISION, pointer :: A(:),CNTL(:),RINFO(:),W(:),ERROR(:),RHS(:),SOL(:),R(:),C(:),LOGDET
!  integer, pointer :: ICNTL(:),INFO(:),IW(:),KEEP(:),JCN1(:),IRN1(:)
!  !integer(kind=myint), pointer ::
!  M=INSIZE(1)
!  N=INSIZE(2)
!  NE=INSIZE(3)
!  MC=INSIZE(4)
!  NC=INSIZE(5)
!  NEC=INSIZE(6)
!  MV=INSIZE(7)
!  NV=INSIZE(8)
!  NEV=INSIZE(9)
!  NEFAC=INSIZE(10)
!  write(*,"(A,i10)") 'NC',NEC
!  write(*,"(A,i5,i5)") 'sizeof',sizeof(JOB),sizeof(I)
!!  write(*,"(i5,i5)") IRN(1),JCN(1)
!!  write(*,"(i5,i5)") IRN(2),JCN(2)
!!  write(*,"(i5,i5)") IRN(3),JCN(3)
!!  write(*,"(i5,i5)") IRN(4),JCN(4)
!  if(NEFAC.EQ.0) then
!    LA=3*NE
!  else
!    LA=ceiling((NEFAC/100.0)*NE)
!  endif
!  MAXN=N
!  IF (N.LT.M) THEN
!    MAXN=M
!  END IF
!  allocate(A(LA),CNTL(10),RINFO(10),W(5*MAXN),ERROR(3),RHS(M),SOL(M),R(MAXN),C(MAXN))
!
!!  real (kind=myreal) A(LA),CNTL(10),LOGDET,RINFO(10),W(MAXN4)
!!  write(*,"(i5,i5)") IRN(4),JCN(4)
!  allocate(JCN1(LA),IRN1(LA),ICNTL(20),INFO(20),IW(6*M+3*N))
!!  write(*,"(i5,i5)") IRN(4),JCN(4)
!  do I=1,NE
!    IRN1(I)=IRN(I)
!    JCN1(I)=JCN(I)
!    A(I)=VA(I)
!    !write(*,"(A,i5,i5,F20.2)") 'Ain',IRN1(I),JCN1(i),VA(I)
!  end do
!  write(*,"(A,i10,i10)") 'Matsize ',M,N
!! Scale matrix
!  LP = 6
!  CALL MC29AD(M,N,NE,A,IRN1,JCN1,R,C,W,LP,IFAIL)
!  IF (IFAIL.LT.0) THEN
!    WRITE (6,'(A)') ' Error in scaling routine'
!    STOP
!  END IF
!  DO I = 1,M
!    R(I) = EXP(R(I))
!  End Do
!  DO I = 1,N
!    C(I) = EXP(C(I))
!  End Do
!  DO II = 1,NE
!    I = IRN1(II)
!    J = JCN1(II)
!    !write(*,"(A,i5,i5,F10.2,F10.2,F10.2)") 'Ainbs',I,J,R(I),C(J),A(II)
!    A(II) = A(II)*R(I)*C(J)
!    !write(*,"(A,i5,i5,F20.2,F20.2)") 'Ain',I,J,R(I)*C(J),A(II)
!  End Do
!!     Factorize matrix
!  CALL MA48ID(CNTL,ICNTL)
!  !ICNTL(4)=0
!  !ICNTL(6)=1
!  T=M+5*N+4*N/ICNTL(6)+7
!  !CNTL(4)=0.00000000000000001
!  !CNTL(2)=0
!  allocate(KEEP(T))
!
!  JOB=1
!  CALL MA48AD(M,N,NE,1,LA,A,IRN1,JCN1,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
!  WRITE (6,FMT='(A,I3/A)') 'INFO(3) =',INFO(3)/NE
!  JOB=1
!  !do I=1,NE
!  !  write(*,"(A,i5,i5,i5,F20.2)") 'AAD',I,IRN1(I),JCN1(I),A(I)
!  !end do
!  CALL MA48BD(M,N,NE,1,LA,A,IRN1,JCN1,KEEP,CNTL,ICNTL,W,IW,INFO,&
!                RINFO)
!  WRITE (6,FMT='(A,I3/A)') 'INFO(4) =',INFO(4)/NE
!  IF (INFO(1).NE.0) THEN
!    WRITE (6,FMT='(A,I3/A)') 'STOP from MA48B/BD with INFO(1) =',&
!    INFO(1),'Solution not possible'
!    STOP
!  END IF
!
!  CALL MA51CD(M,N,LA,A,IRN1,KEEP,SGNDET,LOGDET,W)
!    IF(SGNDET.GT.0)THEN
!      WRITE(6,'(/,A,F9.3)')&
!          ' Determinant is positive; log(determinant) =',LOGDET
!    ELSE IF(SGNDET.LT.0)THEN
!        WRITE(6,'(/,A,F9.3)')&
!          ' Determinant is negative; log(-determinant) =',LOGDET
!      ELSE
!        WRITE(6,'(/,A)')' Determinant is zero'
!    END IF
!
!  do I=1,M
!    RHS(I)=B(I)*R(I)
!    !write(*,"(i5,F10.2)") I,B(I)
!  end do
!  JOB=1
!  TRANS = .FALSE.
!  CALL MA48CD(M,N,TRANS,JOB,LA,A,IRN1,KEEP,CNTL,ICNTL,&
!              RHS,SOL,ERROR,W,IW,INFO)
!  WRITE (6,FMT='(A,I3/A)') 'INFO(1) =',INFO(1)
!  do I=1,M
!    X(I)=SOL(I)*C(I)
!    !write(*,"(A,i5,F20.2)") 'UI',I,X(I)
!  end do
!  L=1
!  JOB=4
!  do J=1,NC
!    do I=1,M
!      RHS(I)=0
!    end do
!    !I=J
!    !write(*,"(A,i10)") 'UI',I
!    !I=NC
!    !write(*,"(A,i10)") 'NC',I
!    checksol=.FALSE.
!    do I=1,NEC
!      if(JCNC(I).EQ.J) then
!        RHS(IRNC(I))=VAC(I)*R(IRNC(I))
!        !if(I.EQ.NEC) then
!          !write(*,"(A,i5,i5,F20.2)") 'RHS',IRNC(I),JCNC(I),VAC(I)
!        !end if
!        !if(RHS(J).NE.0) then
!          checksol=.TRUE.
!        !end if
!      end if
!    end do
!    TRANS = .FALSE.
!    if(checksol) then
!    CALL MA48CD(M,N,TRANS,JOB,LA,A,IRN1,KEEP,CNTL,ICNTL,&
!              RHS,SOL,ERROR,W,IW,INFO)
!    do I=1,M
!      if(SOL(I).NE.0) then
!        VAV(L)=SOL(I)*C(I)
!        IRNV(L)=I
!        JCNV(L)=J
!        !if(L.GE.9450) then
!          !if(J.EQ.7) then
!            !write(*,"(A,i5,i5,F20.10)") 'V',J,I,VAV(L)
!          !end if
!        !end if
!        L=L+1
!      end if
!    end do
!    end if
!  end do
!  INSIZE(9)=L
!  deallocate(A,CNTL,RINFO,W,ERROR,RHS,SOL,R,C)
!  deallocate(IRN1,JCN1,ICNTL,INFO,IW,KEEP)
!END SUBROUTINE SPEC48_MSOL2

SUBROUTINE SPEC48_MSOL(INSIZE,IRN,JCN,VA,B,X,IRNC,JCNC,VAC,IRNB,JCNB,VALUESB,VECBIVI,BIVINZROW0,BIVINZCOL0)!,IRNV,JCNV,VAV
  use constants
  IMPLICIT NONE

  !INTEGER, PARAMETER :: myreal = SELECTED_REAL_KIND(8)!kind(real*8)!KIND( 0.0D0 )
  !INTEGER, PARAMETER :: myint = selected_int_kind(16)!kind(interger*4)!KIND( 1)
  !INTEGER, PARAMETER :: myint32 = selected_int_kind(16)!kind(interger*4)!KIND( 1)
  integer NEFAC,JOB
  !integer (kind=myint) JCN(*),IRN(*),JCNC(*),IRNC(*),JCNV(*),IRNV(*),INSIZE(*)
  integer(4) JCN(*),IRN(*),INSIZE(*),BIVINZROW0(*),BIVINZCOL0(*)!,JCNV(*),IRNV(*)
  integer(4)IRNB(*),JCNB(*)
  integer(4) JCNC(*),IRNC(*)
  !integer (kind=myint32) INSIZE(*)
  integer M,N,NE,T,NC,MC,NEC,RANK,J1!,NV,MV,NEV
  real(kind=DPC) VA(*),VAC(*),B(*),X(*),VALUESB(*),VECBIVI(*)
  !DOUBLE PRECISION LOGDET,SGNDET!,VAV(*)real (8)
  integer LA, MAXN,NBIVI,MBIVI
  !integer (kind=myint) I,J,L,L1
  integer(4) I,J,L,L1,L2,L3,L4,L5,M0,M1,M2,M3,M4,M5,MB,NB,NEB,J2,J3
  LOGICAL TRANS,checksol,isopen
  !DOUBLE PRECISION, pointer :: A(:),CNTL(:),RINFO(:),W(:),ERROR(:),RHS(:),SOL(:)
  real(kind=DPC), pointer :: CNTL(:),RINFO(:),W(:),ERROR1(:),SOL(:)!,VOUT(:)!,RHS(:)A(:),
  integer, pointer :: ICNTL(:),INFO(:),IW(:),KEEP(:),JCNB1(:)!,IRNOUT(:),JCNOUT(:)!,JCN1(:),IRN1(:)
  !integer(kind=myint), pointer ::
  character(len=1024) :: filename
  M=INSIZE(1)
  N=INSIZE(2)
  NE=INSIZE(3)
  MC=INSIZE(4)
  NC=INSIZE(5)
  NEC=INSIZE(6)
  MB=INSIZE(7)
  NB=INSIZE(8)
  NEB=INSIZE(9)
  NEFAC=INSIZE(10)
  RANK=INSIZE(11)
  J1=INSIZE(12)
  MBIVI=INSIZE(14)
  NBIVI=INSIZE(15)
  JCN(1:NE)=JCN(1:NE)+1
  allocate(JCNB1(NEB))
  JCNB1(1:NEB)=JCNB(1:NEB)+1
  do J=1,M
    L=J+NE
    IRN((IRN(L)+1):IRN(L+1))=J
  end do
  !write(*,"(A,i5,A,i5,i5,i5)") 'sizeof',sizeof(JOB),'RANK ',RANK,sizeof(I),NC
  !write(*,"(A,i5,i5)") 'MB NB',MB,NB
!  write(*,"(i5,i5)") IRN(2),JCN(2)
!  write(*,"(i5,i5)") IRN(3),JCN(3)
!  write(*,"(i5,i5)") IRN(4),JCN(4)
  if(NEFAC.EQ.0) then
    LA=2*NE
  else
    LA=ceiling((NEFAC/100.0)*NE)
  endif
  MAXN=N
  IF (N.LT.M) THEN
    MAXN=M
  END IF
  allocate(CNTL(10),RINFO(10),ERROR1(3),SOL(M))!,RHS(M),W(5*MAXN)A(LA),

!  real (kind=myreal) A(LA),CNTL(10),LOGDET,RINFO(10),W(MAXN4)
!  write(*,"(i5,i5)") IRN(4),JCN(4)
  allocate(ICNTL(20),INFO(20),IW(6*M+3*N))!JCN1(LA),IRN1(LA),
!  write(*,"(i5,i5)") IRN(4),JCN(4)
  !do I=1,NE
    !IRN1(I)=IRN(I)
    !JCN1(I)=JCN(I)
    !A(I)=VA(I)
    !write(*,"(A,i5,i5,F20.2)") 'Ain',IRN(I),JCN(i),VA(I)
  !end do
  !write(*,"(A,i10,i10)") 'Matsize ',M,N
!     Factorize matrix
  IF (FSORD.EQ.1) THEN
    CALL MA48ID(CNTL,ICNTL)
  ELSE
    CALL MA48I(CNTL,ICNTL)
  ENDIF
  !ICNTL(6)=1
  T=M+5*N+4*N/ICNTL(6)+7
  INSIZE(13)=T
!      IF (FSORD.EQ.1) THEN
!        CNTL(4)=1e-4
!      else
!        CNTL(4)=0.3
!      ENDIF
  !CNTL(4)=1e-10
  !CNTL(2)=0
  allocate(KEEP(T))
  !KEEP=0
  JOB=1
  IF (FSORD.EQ.1) THEN
    CALL MA48AD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  ELSE
    CALL MA48A(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  ENDIF
  IF (INFO(1).LT.0) THEN
  WRITE (6,'(A,I3)') 'Error STOP from MA48A/AD with INFO(1) =',INFO(1)
  STOP
  END IF
  deallocate(IW)
  allocate(W(M),IW(2*M+2*N))
  IF (FSORD.EQ.1) THEN
    CALL MA48BD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,RINFO)
  else
    CALL MA48B(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,RINFO)
  ENDIF
  IF (INFO(1).NE.0) THEN
    WRITE (6,FMT='(A,I3/A)') 'STOP from MA48B/BD with INFO(1) =',&
    INFO(1),'Solution not possible'
    write(*,"(A,i5)") 'RANK',INFO(5)
    STOP
  END IF
!     Compute the determinant
    !CALL MA51CD(M,N,LA,VA,IRN,KEEP,SGNDET,LOGDET,W)

!     Print result
    !IF(SGNDET.EQ.0)THEN
      !WRITE(6,'(/,A)')' Determinant is zero'
    !ELSE
      !write(*,"(A,F20.2)")' Log(determinant) =',LOGDET
    !END IF
  !do I=1,M
    !RHS(I)=B(I)
    !write(*,"(A,i5,F10.2)") 'B',I,B(I)
  !end do
  !deallocate(JCN1)
  JOB=1
  deallocate(W,IW)
  IF (JOB.EQ.1) THEN
    allocate(W(2*MAXN))
  else
    allocate(W(4*MAXN))
  END IF
  allocate(IW(MAXN))!A(LA),IRN1(LA),
  TRANS = .FALSE.
  !ICNTL(5)=0
  if(sum(abs(B(1:N))).GT.0) then
    IF (FSORD.EQ.1) THEN
      CALL MA48CD(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
              B,X,ERROR1,W,IW,INFO)
    ELSE
      CALL MA48C(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
              B,X,ERROR1,W,IW,INFO)
    ENDIF
  end if
  !WRITE (6,FMT='(A,I3/A)') 'INFO(1) =',INFO(1)
  !X(1:M)=SOL(1:M)
  !do I=1,M
    !X(I)=SOL(I)
    !write(*,"(A,i5,F20.2)") 'XI',I,SOL(I)
  !end do
  L=0
  TRANS = .FALSE.
  do J=1,MC-1!NC
    !do I=1,M
      !B(I)=0!RHS(I)=0
    !end do
    !write(*,"(A,i5,i5,i5)") 'JC',J,L,IRNC(J)
      J3=IRNC(J+1)-IRNC(J)
    if(J3.GT.0) then
      !checksol=.TRUE.
      B(1:M)=0
      do I=1,J3!NEC
        !write(*,"(A,i5,A,i5)") 'I',I,'JNC',JCNC(I+L)
        L2=I+L
        B(JCNC(L2)+1)=VAC(L2)
        !if(JCNC(I).EQ.J) then
          !B(IRNC(I))=VAC(I)!RHS(IRNC(I))=VAC(I)
          !if(J.EQ.7) then
          !write(*,"(A,i5,i5,F20.2)") 'RHS',J,JCNC(I+L),VAC(I+L)
        !end if
        !if(RHS(J).NE.0) then
          !checksol=.TRUE.
        !end if
      !end if
      end do
      L=L+J3
    !end if
    !if(checksol) then
      IF (FSORD.EQ.1) THEN
        CALL MA48CD(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
                  B,SOL,ERROR1,W,IW,INFO)
      else
        CALL MA48C(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
                  B,SOL,ERROR1,W,IW,INFO)
      ENDIF
      !WRITE (6,FMT='(A,I3/A)') 'INFO(1) =',INFO(1)
      do I=1,MB-1
        if(IRNB(I).NE.IRNB(I+1)) then !GE.1E-10.OR.SOL(I).LE.-1E-10) then
              L3=BIVINZCOL0(J)+BIVINZROW0(I)!(JCNB(j2)+1)!+NBIVI*(BIVINZROW0(JCNB(j2)+1)-1)
              DO j2=IRNB(I)+1,IRNB(I+1)
                !if(SOL(JCNB1(j2)).NE.0) then
                !$OMP ATOMIC
                VECBIVI(L3)=VECBIVI(L3)-SOL(JCNB1(j2))*VALUESB(j2)
                !if(L3.EQ.1) then
                !write(*,"(A,i10,i10,i10,F20.2,F10.2,F10.2)") 'L3',L3,J,j2,VECBIVI(L3),SOL(I),VALUESB(j2)
                !end if
                !end if
              end do
        end if
      end do
        if(IRNB(MB).NE.NEB) then !GE.1E-10.OR.SOL(I).LE.-1E-10) then
              L3=BIVINZCOL0(J)+BIVINZROW0(MB)!(JCNB(j2)+1)!+NBIVI*(BIVINZROW0(JCNB(j2)+1)-1)
              do j2=IRNB(MB)+1,NEB
                !if(VALUESB(j2).NE.0) then
                !L3=J+NB*JCNB(j2)
                !$OMP ATOMIC
                VECBIVI(L3)=VECBIVI(L3)-SOL(JCNB1(j2))*VALUESB(j2)
                !if(L3.EQ.1) then
                !write(*,"(A,i10,F20.2)") 'L3',L3,VECBIVI(L3)
                !end if
                !end if
              end do
        end if
    end if
    !write(*,"(A,i10)") 'L',JCNOUT(J)
  end do

 !do J=1,MC!NC
  J=MC
    !do I=1,M
      !B(I)=0!RHS(I)=0
    !end do
    !write(*,"(A,i5,i5,i5)") 'JC',J,L,IRNC(J)
      J3=NEC-IRNC(J)
    !write(*,"(A,i15,A,i15,A,i15)") 'J',J,'J3',J3,'IRNC',IRNC(J)
    if(J3.GT.0) then
      !checksol=.TRUE.
      B(1:M)=0
      do I=1,J3!NEC
        !write(*,"(A,i5,A,i5)") 'I',I,'JNC',JCNC(I+L)
        L2=I+L
        B(JCNC(L2)+1)=VAC(L2)
        !if(JCNC(I).EQ.J) then
          !B(IRNC(I))=VAC(I)!RHS(IRNC(I))=VAC(I)
          !if(J.EQ.7) then
          !write(*,"(A,i5,i5,F20.2)") 'RHS',J,JCNC(I+L),VAC(I+L)
        !end if
        !if(RHS(J).NE.0) then
          !checksol=.TRUE.
        !end if
      !end if
      end do
      L=L+J3
    !end if
    !if(checksol) then
      IF (FSORD.EQ.1) THEN
        CALL MA48CD(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
                  B,SOL,ERROR1,W,IW,INFO)
      ELSE
        CALL MA48C(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
                  B,SOL,ERROR1,W,IW,INFO)
      ENDIF
      do I=1,MB-1
        if(IRNB(I).NE.IRNB(I+1)) then !GE.1E-10.OR.SOL(I).LE.-1E-10) then
              L3=BIVINZCOL0(J)+BIVINZROW0(I)!(JCNB(j2)+1)!+NBIVI*(BIVINZROW0(JCNB(j2)+1)-1)
              DO j2=IRNB(I)+1,IRNB(I+1)
                !if(SOL(JCNB1(j2)).NE.0) then
                !$OMP ATOMIC
                VECBIVI(L3)=VECBIVI(L3)-SOL(JCNB1(j2))*VALUESB(j2)
                !if(L3.EQ.1) then
                !write(*,"(A,i10,i10,i10,F20.2,F10.2,F10.2)") 'L3',L3,J,j2,VECBIVI(L3),SOL(I),VALUESB(j2)
                !end if
                !end if
              end do
        end if
      end do
        if(IRNB(MB).NE.NEB) then !GE.1E-10.OR.SOL(I).LE.-1E-10) then
              L3=BIVINZCOL0(J)+BIVINZROW0(MB)!(JCNB(j2)+1)!+NBIVI*(BIVINZROW0(JCNB(j2)+1)-1)
              do j2=IRNB(MB)+1,NEB
                !if(VALUESB(j2).NE.0) then
                !L3=J+NB*JCNB(j2)
                !$OMP ATOMIC
                VECBIVI(L3)=VECBIVI(L3)-SOL(JCNB1(j2))*VALUESB(j2)
                !if(L3.EQ.1) then
                !write(*,"(A,i10,F20.2)") 'L3',L3,VECBIVI(L3)
                !end if
                !end if
              end do
        end if
    end if
    !write(*,"(A,i10)") 'L',JCNOUT(J)
  !end do
  J=INSIZE(16)
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A4,I4.4,I4.4,A4)") "_vav",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) VA(1:LA)
  CLOSE(J)
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A5,I4.4,I4.4,A4)") "_irnv",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) IRN(1:LA)
  CLOSE(J)
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A5,I4.4,I4.4,A4)") "_keep",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) KEEP(1:T)
  CLOSE(J)
  !INSIZE(9)=L
  !write(*,"(A,F10.2)") 'VIBI',VECBIVI(1)
  !deallocate(VOUT,IRNOUT)
  deallocate(CNTL,RINFO,W,ERROR1,SOL,JCNB1)!,RHSA,
  deallocate(ICNTL,INFO,IW,KEEP)!IRN1,,JCNOUT
END SUBROUTINE SPEC48_MSOL

SUBROUTINE SPEC48M_MSOL(INSIZE,IRN,JCN,VA,B,X,IRNC,JCNC,VAC,IRNB,JCNB,VALUESB,VECBIVI,BIVINZROW0,BIVINZCOL0)!,IRNV,JCNV,VAV
  use constants
  IMPLICIT NONE

  !INTEGER, PARAMETER :: myreal = SELECTED_REAL_KIND(8)!kind(real*8)!KIND( 0.0D0 )
  !INTEGER, PARAMETER :: myint = selected_int_kind(16)!kind(interger*4)!KIND( 1)
  !INTEGER, PARAMETER :: myint32 = selected_int_kind(16)!kind(interger*4)!KIND( 1)
  integer NEFAC,JOB
  !integer (kind=myint) JCN(*),IRN(*),JCNC(*),IRNC(*),JCNV(*),IRNV(*),INSIZE(*)
  integer(4) JCN(*),IRN(*),INSIZE(*),BIVINZROW0(*),BIVINZCOL0(*)!,JCNV(*),IRNV(*)
  integer(4)IRNB(*),JCNB(*)
  integer(4) JCNC(*),IRNC(*)
  !integer (kind=myint32) INSIZE(*)
  integer M,N,NE,T,NC,MC,NEC,RANK,J1!,NV,MV,NEV
  real(kind=DPC) VA(*),VAC(*),B(*),X(*),VALUESB(*),VECBIVI(*)
  !DOUBLE PRECISION LOGDET,SGNDET!,VAV(*)real (8)
  integer LA, MAXN,NBIVI,MBIVI
  !integer (kind=myint) I,J,L,L1
  integer(4) I,J,L,L1,L2,L3,L4,L5,M0,M1,M2,M3,M4,M5,MB,NB,NEB,J2,J3
  LOGICAL TRANS,checksol,isopen
  !DOUBLE PRECISION, pointer :: A(:),CNTL(:),RINFO(:),W(:),ERROR(:),RHS(:),SOL(:)
  real(kind=DPC), pointer :: CNTL(:),RINFO(:),W(:),ERROR1(:),SOL(:)!,VOUT(:)!,RHS(:)A(:),
  integer, pointer :: ICNTL(:),INFO(:),IW(:),KEEP(:),JCNB1(:)!,IRNOUT(:),JCNOUT(:)!,JCN1(:),IRN1(:)
  !integer(kind=myint), pointer ::
  character(len=1024) :: filename
  M=INSIZE(1)
  N=INSIZE(2)
  NE=INSIZE(3)
  MC=INSIZE(4)
  NC=INSIZE(5)
  NEC=INSIZE(6)
  MB=INSIZE(7)
  NB=INSIZE(8)
  NEB=INSIZE(9)
  NEFAC=INSIZE(10)
  RANK=INSIZE(11)
  J1=INSIZE(12)
  MBIVI=INSIZE(14)
  NBIVI=INSIZE(15)
  JCN(1:NE)=JCN(1:NE)+1
  allocate(JCNB1(NEB))
  JCNB1(1:NEB)=JCNB(1:NEB)+1
  do J=1,M
    L=J+NE
    IRN((IRN(L)+1):IRN(L+1))=J
  end do
  !write(*,"(A,i5,A,i5,i5,i5)") 'sizeof',sizeof(JOB),'RANK ',RANK,sizeof(I),NC
  !write(*,"(A,i5,i5)") 'MB NB',MB,NB
!  write(*,"(i5,i5)") IRN(2),JCN(2)
!  write(*,"(i5,i5)") IRN(3),JCN(3)
!  write(*,"(i5,i5)") IRN(4),JCN(4)
  if(NEFAC.EQ.0) then
    LA=2*NE
  else
    LA=ceiling((NEFAC/100.0)*NE)
  endif
  if(LA.NE.INSIZE(17)) then
    !write(*,"(A,i15,i15)") 'Note LA!!!!',LA,INSIZE(17)
    LA=INSIZE(17)
  end if
  MAXN=N
  IF (N.LT.M) THEN
    MAXN=M
  END IF
  allocate(CNTL(10),RINFO(10),ERROR1(3),SOL(M))!,RHS(M),W(5*MAXN)A(LA),

!  real (kind=myreal) A(LA),CNTL(10),LOGDET,RINFO(10),W(MAXN4)
!  write(*,"(i5,i5)") IRN(4),JCN(4)
  allocate(ICNTL(20),INFO(20),IW(6*M+3*N))!JCN1(LA),IRN1(LA),
!  write(*,"(i5,i5)") IRN(4),JCN(4)
  !do I=1,NE
    !IRN1(I)=IRN(I)
    !JCN1(I)=JCN(I)
    !A(I)=VA(I)
    !write(*,"(A,i5,i5,F20.2)") 'Ain',IRN(I),JCN(i),VA(I)
  !end do
  !write(*,"(A,i10,i10)") 'Matsize ',M,N
!     Factorize matrix
  IF (FSORD.EQ.1) THEN
    CALL MA48ID(CNTL,ICNTL)
  ELSE
    CALL MA48I(CNTL,ICNTL)
  ENDIF
  !ICNTL(6)=1
  T=M+5*N+4*N/ICNTL(6)+7
  INSIZE(13)=T
!      IF (FSORD.EQ.1) THEN
!        CNTL(4)=1e-4
!      else
!        CNTL(4)=0.3
!      ENDIF
  !CNTL(4)=1e-10
  !CNTL(2)=0
  allocate(KEEP(T))
  !KEEP=0
  JOB=1
  IF (FSORD.EQ.1) THEN
    CALL MA48AD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  ELSE
    CALL MA48A(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  ENDIF
  IF (INFO(1).LT.0) THEN
  WRITE (6,'(A,I3)') 'Error STOP from MA48A/AD with INFO(1) =',INFO(1)
  STOP
  END IF
  deallocate(IW)
  allocate(W(M),IW(2*M+2*N))
  IF (FSORD.EQ.1) THEN
    CALL MA48BD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,RINFO)
  else
    CALL MA48B(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,RINFO)
  ENDIF
  IF (INFO(1).NE.0) THEN
    WRITE (6,FMT='(A,I3/A)') 'STOP from MA48B/BD with INFO(1) =',&
    INFO(1),'Solution not possible'
    write(*,"(A,i5)") 'RANK',INFO(5)
    STOP
  END IF
!     Compute the determinant
    !CALL MA51CD(M,N,LA,VA,IRN,KEEP,SGNDET,LOGDET,W)

!     Print result
    !IF(SGNDET.EQ.0)THEN
      !WRITE(6,'(/,A)')' Determinant is zero'
    !ELSE
      !write(*,"(A,F20.2)")' Log(determinant) =',LOGDET
    !END IF
  !do I=1,M
    !RHS(I)=B(I)
    !write(*,"(A,i5,F10.2)") 'B',I,B(I)
  !end do
  !deallocate(JCN1)
  JOB=1
  deallocate(W,IW)
  IF (JOB.EQ.1) THEN
    allocate(W(2*MAXN))
  else
    allocate(W(4*MAXN))
  END IF
  allocate(IW(MAXN))!A(LA),IRN1(LA),
  TRANS = .FALSE.
  !ICNTL(5)=0
  if(sum(abs(B(1:N))).GT.0) then
    IF (FSORD.EQ.1) THEN
      CALL MA48CD(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
              B,X,ERROR1,W,IW,INFO)
    ELSE
      CALL MA48C(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
              B,X,ERROR1,W,IW,INFO)
    ENDIF
  end if
  !WRITE (6,FMT='(A,I3/A)') 'INFO(1) =',INFO(1)
  !X(1:M)=SOL(1:M)
  !do I=1,M
    !X(I)=SOL(I)
    !write(*,"(A,i5,F20.2)") 'XI',I,SOL(I)
  !end do
  L=0
  TRANS = .FALSE.
  do J=1,MC-1!NC
    !do I=1,M
      !B(I)=0!RHS(I)=0
    !end do
    !write(*,"(A,i5,i5,i5)") 'JC',J,L,IRNC(J)
      J3=IRNC(J+1)-IRNC(J)
    if(J3.GT.0) then
      !checksol=.TRUE.
      B(1:M)=0
      do I=1,J3!NEC
        !write(*,"(A,i5,A,i5)") 'I',I,'JNC',JCNC(I+L)
        L2=I+L
        B(JCNC(L2)+1)=VAC(L2)
        !if(JCNC(I).EQ.J) then
          !B(IRNC(I))=VAC(I)!RHS(IRNC(I))=VAC(I)
          !if(J.EQ.7) then
          !write(*,"(A,i5,i5,F20.2)") 'RHS',J,JCNC(I+L),VAC(I+L)
        !end if
        !if(RHS(J).NE.0) then
          !checksol=.TRUE.
        !end if
      !end if
      end do
      L=L+J3
    !end if
    !if(checksol) then
      IF (FSORD.EQ.1) THEN
        CALL MA48CD(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
                  B,SOL,ERROR1,W,IW,INFO)
      else
        CALL MA48C(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
                  B,SOL,ERROR1,W,IW,INFO)
      ENDIF
      !WRITE (6,FMT='(A,I3/A)') 'INFO(1) =',INFO(1)
      do I=1,MB-1
        if(IRNB(I).NE.IRNB(I+1)) then !GE.1E-10.OR.SOL(I).LE.-1E-10) then
              L3=BIVINZCOL0(J)+BIVINZROW0(I)!(JCNB(j2)+1)!+NBIVI*(BIVINZROW0(JCNB(j2)+1)-1)
              DO j2=IRNB(I)+1,IRNB(I+1)
                !if(SOL(JCNB1(j2)).NE.0) then
                !$OMP ATOMIC
                VECBIVI(L3)=VECBIVI(L3)-SOL(JCNB1(j2))*VALUESB(j2)
                !if(L3.EQ.1) then
                !write(*,"(A,i10,i10,i10,F20.2,F10.2,F10.2)") 'L3',L3,J,j2,VECBIVI(L3),SOL(I),VALUESB(j2)
                !end if
                !end if
              end do
        end if
      end do
        if(IRNB(MB).NE.NEB) then !GE.1E-10.OR.SOL(I).LE.-1E-10) then
              L3=BIVINZCOL0(J)+BIVINZROW0(MB)!(JCNB(j2)+1)!+NBIVI*(BIVINZROW0(JCNB(j2)+1)-1)
              do j2=IRNB(MB)+1,NEB
                !if(VALUESB(j2).NE.0) then
                !L3=J+NB*JCNB(j2)
                !$OMP ATOMIC
                VECBIVI(L3)=VECBIVI(L3)-SOL(JCNB1(j2))*VALUESB(j2)
                !if(L3.EQ.1) then
                !write(*,"(A,i10,F20.2)") 'L3',L3,VECBIVI(L3)
                !end if
                !end if
              end do
        end if
    end if
    !write(*,"(A,i10)") 'L',JCNOUT(J)
  end do

 !do J=1,MC!NC
  J=MC
    !do I=1,M
      !B(I)=0!RHS(I)=0
    !end do
    !write(*,"(A,i5,i5,i5)") 'JC',J,L,IRNC(J)
      J3=NEC-IRNC(J)
    !write(*,"(A,i15,A,i15,A,i15)") 'J',J,'J3',J3,'IRNC',IRNC(J)
    if(J3.GT.0) then
      !checksol=.TRUE.
      B(1:M)=0
      do I=1,J3!NEC
        !write(*,"(A,i5,A,i5)") 'I',I,'JNC',JCNC(I+L)
        L2=I+L
        B(JCNC(L2)+1)=VAC(L2)
        !if(JCNC(I).EQ.J) then
          !B(IRNC(I))=VAC(I)!RHS(IRNC(I))=VAC(I)
          !if(J.EQ.7) then
          !write(*,"(A,i5,i5,F20.2)") 'RHS',J,JCNC(I+L),VAC(I+L)
        !end if
        !if(RHS(J).NE.0) then
          !checksol=.TRUE.
        !end if
      !end if
      end do
      L=L+J3
    !end if
    !if(checksol) then
      IF (FSORD.EQ.1) THEN
        CALL MA48CD(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
                  B,SOL,ERROR1,W,IW,INFO)
      ELSE
        CALL MA48C(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
                  B,SOL,ERROR1,W,IW,INFO)
      ENDIF
      do I=1,MB-1
        if(IRNB(I).NE.IRNB(I+1)) then !GE.1E-10.OR.SOL(I).LE.-1E-10) then
              L3=BIVINZCOL0(J)+BIVINZROW0(I)!(JCNB(j2)+1)!+NBIVI*(BIVINZROW0(JCNB(j2)+1)-1)
              DO j2=IRNB(I)+1,IRNB(I+1)
                !if(SOL(JCNB1(j2)).NE.0) then
                !$OMP ATOMIC
                VECBIVI(L3)=VECBIVI(L3)-SOL(JCNB1(j2))*VALUESB(j2)
                !if(L3.EQ.1) then
                !write(*,"(A,i10,i10,i10,F20.2,F10.2,F10.2)") 'L3',L3,J,j2,VECBIVI(L3),SOL(I),VALUESB(j2)
                !end if
                !end if
              end do
        end if
      end do
        if(IRNB(MB).NE.NEB) then !GE.1E-10.OR.SOL(I).LE.-1E-10) then
              L3=BIVINZCOL0(J)+BIVINZROW0(MB)!(JCNB(j2)+1)!+NBIVI*(BIVINZROW0(JCNB(j2)+1)-1)
              do j2=IRNB(MB)+1,NEB
                !if(VALUESB(j2).NE.0) then
                !L3=J+NB*JCNB(j2)
                !$OMP ATOMIC
                VECBIVI(L3)=VECBIVI(L3)-SOL(JCNB1(j2))*VALUESB(j2)
                !if(L3.EQ.1) then
                !write(*,"(A,i10,F20.2)") 'L3',L3,VECBIVI(L3)
                !end if
                !end if
              end do
        end if
    end if
    !write(*,"(A,i10)") 'L',JCNOUT(J)
  !end do
  J=INSIZE(16)
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A4,I4.4,I4.4,A4)") "_vav",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) VA(1:LA)
  CLOSE(J)
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A5,I4.4,I4.4,A4)") "_irnv",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) IRN(1:LA)
  CLOSE(J)
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A5,I4.4,I4.4,A4)") "_keep",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) KEEP(1:T)
  CLOSE(J)
  !INSIZE(9)=L
  !write(*,"(A,F10.2)") 'VIBI',VECBIVI(1)
  !deallocate(VOUT,IRNOUT)
  deallocate(CNTL,RINFO,W,ERROR1,SOL,JCNB1)!,RHSA,
  deallocate(ICNTL,INFO,IW,KEEP)!IRN1,,JCNOUT
END SUBROUTINE SPEC48M_MSOL

SUBROUTINE SPEC48_ESOL(INSIZE,IRN,VA,KEEP,B,SOL)
  use constants
  IMPLICIT NONE

  integer NEFAC,JOB
  integer(4) INSIZE(*),IRN(*),KEEP(*)
  integer M,N,NE
  real (kind=DPC) VA(*),B(*),SOL(*)
  integer LA, MAXN,I
  real (kind=DPC), pointer :: CNTL(:),RINFO(:),W(:),ERROR1(:)
  integer, pointer :: ICNTL(:),INFO(:),IW(:)
  LOGICAL TRANS
  M=INSIZE(1)
  N=INSIZE(2)
  NE=INSIZE(3)
  NEFAC=INSIZE(10)
  !write(*,"(A,i10,i10)") 'Matsize ',M,N
  if(NEFAC.EQ.0) then
    LA=2*NE
  else
    LA=ceiling((NEFAC/100.0)*NE)
  endif
  MAXN=N
  IF (N.LT.M) THEN
    MAXN=M
  END IF
  allocate(CNTL(10),RINFO(10),ERROR1(3))
  allocate(ICNTL(20),INFO(20))

  IF (FSORD.EQ.1) THEN
    CALL MA48ID(CNTL,ICNTL)
  else
    CALL MA48I(CNTL,ICNTL)
  ENDIF
  JOB=1
  IF (JOB.EQ.1) THEN
    allocate(W(2*MAXN))
  else
    allocate(W(4*MAXN))
  END IF
  allocate(IW(MAXN))
  TRANS = .FALSE.
  IF (FSORD.EQ.1) THEN
    CALL MA48CD(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
              B,SOL,ERROR1,W,IW,INFO)
  else
    CALL MA48C(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
              B,SOL,ERROR1,W,IW,INFO)
  ENDIF
  !WRITE (6,FMT='(A,I3/A)') 'INFO(1) =',INFO(1)
  !do I=1,M
    !write(*,"(i5,F20.2)") I,SOL(I)
  !end do
  deallocate(CNTL,RINFO,W,ERROR1)
  deallocate(ICNTL,INFO,IW)
END SUBROUTINE SPEC48_ESOL

SUBROUTINE SPEC48M_ESOL(INSIZE,IRN,VA,KEEP,B,SOL)
  use constants
  IMPLICIT NONE

  integer NEFAC,JOB
  integer(4) INSIZE(*),IRN(*),KEEP(*)
  integer M,N,NE
  real (kind=DPC) VA(*),B(*),SOL(*)
  integer LA, MAXN,I
  real (kind=DPC), pointer :: CNTL(:),RINFO(:),W(:),ERROR1(:)
  integer, pointer :: ICNTL(:),INFO(:),IW(:)
  LOGICAL TRANS
  M=INSIZE(1)
  N=INSIZE(2)
  NE=INSIZE(3)
  NEFAC=INSIZE(10)
  !write(*,"(A,i10,i10)") 'Matsize ',M,N
  if(NEFAC.EQ.0) then
    LA=2*NE
  else
    LA=ceiling((NEFAC/100.0)*NE)
  endif
  if(LA.NE.INSIZE(17)) then
    !write(*,"(A,i15,i15)") 'Note LA!!!!',LA,INSIZE(17)
    LA=INSIZE(17)
  end if
  MAXN=N
  IF (N.LT.M) THEN
    MAXN=M
  END IF
  allocate(CNTL(10),RINFO(10),ERROR1(3))
  allocate(ICNTL(20),INFO(20))

  IF (FSORD.EQ.1) THEN
    CALL MA48ID(CNTL,ICNTL)
  else
    CALL MA48I(CNTL,ICNTL)
  ENDIF
  JOB=1
  IF (JOB.EQ.1) THEN
    allocate(W(2*MAXN))
  else
    allocate(W(4*MAXN))
  END IF
  allocate(IW(MAXN))
  TRANS = .FALSE.
  IF (FSORD.EQ.1) THEN
    CALL MA48CD(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
              B,SOL,ERROR1,W,IW,INFO)
  else
    CALL MA48C(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
              B,SOL,ERROR1,W,IW,INFO)
  ENDIF
  !WRITE (6,FMT='(A,I3/A)') 'INFO(1) =',INFO(1)
  !do I=1,M
    !write(*,"(i5,F20.2)") I,SOL(I)
  !end do
  deallocate(CNTL,RINFO,W,ERROR1)
  deallocate(ICNTL,INFO,IW)
END SUBROUTINE SPEC48M_ESOL

SUBROUTINE SPEC48_RPESOL(INSIZE,IRN,VA,KEEP,B,SOL,CNTL,RINFO,ERROR1,ICNTL,INFO,W,IW)
  use constants
  IMPLICIT NONE

  integer NEFAC,JOB
  integer(4) INSIZE(*),IRN(*),KEEP(*),ICNTL(*),INFO(*),IW(*)
  integer M,N,NE
  real(kind=DPC) VA(*),B(*),SOL(*)
  real(kind=DPC) CNTL(*),RINFO(*),W(*),ERROR1(*)
  integer LA!, MAXN,I
  !DOUBLE PRECISION, pointer :: CNTL(:),RINFO(:),W(:),ERROR1(:)
  !integer, pointer :: ICNTL(:),INFO(:),IW(:)
  LOGICAL TRANS
  !M=INSIZE(1)
  N=INSIZE(2)
  NE=INSIZE(3)
  NEFAC=INSIZE(10)
  !write(*,"(A,i10,i10)") 'Matsize ',M,N
  !if(NEFAC.EQ.0) then
    !LA=2*NE
  !else
    LA=ceiling((NEFAC/100.0)*NE)
  !endif
  !MAXN=N
  !IF (N.LT.M) THEN
    !MAXN=M
  !END IF
  !allocate(CNTL(10),RINFO(10),ERROR1(3))
  !allocate(ICNTL(20),INFO(20))

  IF (FSORD.EQ.1) THEN
    CALL MA48ID(CNTL,ICNTL)
  else
    CALL MA48I(CNTL,ICNTL)
  ENDIF
  JOB=1
  !IF (JOB.EQ.1) THEN
    !allocate(W(2*MAXN))
  !else
    !allocate(W(4*MAXN))
  !END IF
  !allocate(IW(MAXN))
  TRANS = .FALSE.
  IF (FSORD.EQ.1) THEN
    CALL MA48CD(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
              B,SOL,ERROR1,W,IW,INFO)
  else
    CALL MA48C(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
              B,SOL,ERROR1,W,IW,INFO)
  ENDIF
  !WRITE (6,FMT='(A,I3/A)') 'INFO(1) =',INFO(1)
  !do I=1,M
    !write(*,"(i5,F20.2)") I,SOL(I)
  !end do
  !deallocate(CNTL,RINFO,W,ERROR1)
  !deallocate(ICNTL,INFO,IW)
END SUBROUTINE SPEC48_RPESOL

SUBROUTINE SPEC48M_RPESOL(INSIZE,IRN,VA,KEEP,B,SOL,CNTL,RINFO,ERROR1,ICNTL,INFO,W,IW)
  use constants
  IMPLICIT NONE

  integer NEFAC,JOB
  integer(4) INSIZE(*),IRN(*),KEEP(*),ICNTL(*),INFO(*),IW(*)
  integer M,N,NE
  real(kind=DPC) VA(*),B(*),SOL(*)
  real(kind=DPC) CNTL(*),RINFO(*),W(*),ERROR1(*)
  integer LA!, MAXN,I
  !DOUBLE PRECISION, pointer :: CNTL(:),RINFO(:),W(:),ERROR1(:)
  !integer, pointer :: ICNTL(:),INFO(:),IW(:)
  LOGICAL TRANS
  !M=INSIZE(1)
  N=INSIZE(2)
  NE=INSIZE(3)
  NEFAC=INSIZE(10)
  !write(*,"(A,i10,i10)") 'Matsize ',M,N
  !if(NEFAC.EQ.0) then
    !LA=2*NE
  !else
    LA=ceiling((NEFAC/100.0)*NE)
  if(LA.NE.INSIZE(17)) then
    !write(*,"(A,i15,i15)") 'Note LA!!!!',LA,INSIZE(17)
    LA=INSIZE(17)
  end if
  !endif
  !MAXN=N
  !IF (N.LT.M) THEN
    !MAXN=M
  !END IF
  !allocate(CNTL(10),RINFO(10),ERROR1(3))
  !allocate(ICNTL(20),INFO(20))

  IF (FSORD.EQ.1) THEN
    CALL MA48ID(CNTL,ICNTL)
  else
    CALL MA48I(CNTL,ICNTL)
  ENDIF
  JOB=1
  !IF (JOB.EQ.1) THEN
    !allocate(W(2*MAXN))
  !else
    !allocate(W(4*MAXN))
  !END IF
  !allocate(IW(MAXN))
  TRANS = .FALSE.
  IF (FSORD.EQ.1) THEN
    CALL MA48CD(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
              B,SOL,ERROR1,W,IW,INFO)
  else
    CALL MA48C(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
              B,SOL,ERROR1,W,IW,INFO)
  ENDIF
  !WRITE (6,FMT='(A,I3/A)') 'INFO(1) =',INFO(1)
  !do I=1,M
    !write(*,"(i5,F20.2)") I,SOL(I)
  !end do
  !deallocate(CNTL,RINFO,W,ERROR1)
  !deallocate(ICNTL,INFO,IW)
END SUBROUTINE SPEC48M_RPESOL


SUBROUTINE SPEC48_SSOL(INSIZE,IRN,JCN,VA,B,X)
  use constants
  IMPLICIT NONE

  INTEGER, PARAMETER :: myreal = SELECTED_REAL_KIND(8)!kind(real*8)!KIND( 0.0D0 )
  INTEGER, PARAMETER :: myint = selected_int_kind(16)!kind(interger*4)!KIND( 1)
  INTEGER, PARAMETER :: myint32 = selected_int_kind(16)!kind(interger*4)!KIND( 1)
  integer NEFAC,JOB
  integer (kind=myint) JCN(*),IRN(*)
  integer (kind=myint32) INSIZE(*)
  integer M,N,NE,T
  real (kind=DPC) VA(*),B(*),X(*)
  integer I,J,L,LA, MAXN,SGNDET!,IFAIL,LP,II
  LOGICAL TRANS,checksol
  real (kind=DPC), pointer :: A(:),CNTL(:),RINFO(:),W(:),ERROR(:)!,LOGDET!,R(:),C(:)!,RHS(:),SOL(:)
  integer, pointer :: ICNTL(:),INFO(:),IW(:),KEEP(:),JCN1(:),IRN1(:)
  M=INSIZE(1)
  N=INSIZE(2)
  NE=INSIZE(3)
  NEFAC=INSIZE(4)
  if(NEFAC.EQ.0) then
    LA=2*NE
  else
    LA=ceiling((NEFAC/100.0)*NE)
  endif
  MAXN=N
  IF (N.LT.M) THEN
    MAXN=M
  END IF
  allocate(A(LA),CNTL(10),RINFO(10),ERROR(3))!,RHS(M),SOL(M),R(MAXN),C(MAXN),W(5*MAXN)
  allocate(JCN1(LA),IRN1(LA),ICNTL(20),INFO(20),IW(6*M+3*N))

  do I=1,NE
    IRN1(I)=IRN(I)
    JCN1(I)=JCN(I)
    A(I)=VA(I)
    !write(*,"(A,i5,i5,F20.2)") 'Ain',IRN1(I),JCN1(i),VA(I)
  end do
  write(*,"(A,i10,i10,i15)") 'Matsize ',M,N,NE
  !LP = 6
  !CALL MC29AD(M,N,NE,A,IRN1,JCN1,R,C,W,LP,IFAIL)
  !IF (IFAIL.LT.0) THEN
  !  WRITE (6,'(A)') ' Error in scaling routine'
  !  STOP
  !END IF
  !DO I = 1,M
  !  R(I) = EXP(R(I))
  !End Do
  !DO I = 1,N
  !  C(I) = EXP(C(I))
  !End Do
  !DO II = 1,NE
  !  I = IRN1(II)
  !  J = JCN1(II)
    !write(*,"(A,i5,i5,F10.2,F10.2,F10.2)") 'Ainbs',I,J,R(I),C(J),A(II)
  !  A(II) = A(II)*R(I)*C(J)
    !write(*,"(A,i5,i5,F20.2,F20.2)") 'Ain',I,J,R(I)*C(J),A(II)
  !End Do
!     Factorize matrix
  IF (FSORD.EQ.1) THEN
    CALL MA48ID(CNTL,ICNTL)
  else
    CALL MA48I(CNTL,ICNTL)
  endif
  !ICNTL(4)=0
  !ICNTL(6)=1
  T=M+5*N+4*N/ICNTL(6)+7
  !CNTL(4)=0.00000000000000001
  !CNTL(2)=0
  allocate(KEEP(T))

  JOB=1
  IF (FSORD.EQ.1) THEN
    CALL MA48AD(M,N,NE,1,LA,A,IRN1,JCN1,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  else
    CALL MA48A(M,N,NE,1,LA,A,IRN1,JCN1,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  endif
  IF (INFO(1).LT.0) THEN
  WRITE (6,'(A,I3)') 'Error STOP from MA48A/AD with INFO(1) =',INFO(1)
  WRITE (6,FMT='(A,I3/A)') 'INFO(3) =',INFO(3)/NE
  STOP
  END IF
  JOB=1
  !do I=1,NE
  !  write(*,"(A,i5,i5,i5,F20.2)") 'AAD',I,IRN1(I),JCN1(I),A(I)
  !end do
  deallocate(IW)
  allocate(W(M),IW(2*M+2*N))
  IF (FSORD.EQ.1) THEN
    CALL MA48BD(M,N,NE,1,LA,A,IRN1,JCN1,KEEP,CNTL,ICNTL,W,IW,INFO,&
                RINFO)
  else
    CALL MA48B(M,N,NE,1,LA,A,IRN1,JCN1,KEEP,CNTL,ICNTL,W,IW,INFO,&
                RINFO)
  endif
  WRITE (6,FMT='(A,I3/A)') 'INFO(4) =',INFO(4)/NE
  IF (INFO(1).NE.0) THEN
    WRITE (6,FMT='(A,I3/A)') 'STOP from MA48B/BD with INFO(1) =',&
    INFO(1),'Solution not possible'
    STOP
  END IF

  !do I=1,M
    !RHS(I)=B(I)
    !write(*,"(i5,F10.2)") I,B(I)
  !end do
  !do I=1,M
  !  B(I)=B(I)*R(I)
    !write(*,"(i5,F10.2)") I,B(I)
  !end do
  deallocate(JCN1)
  JOB=1
  deallocate(W,IW)
  IF (JOB.EQ.1) THEN
    allocate(W(2*MAXN))
  else
    allocate(W(4*MAXN))
  END IF
  allocate(IW(MAXN))!A(LA),IRN1(LA),
  TRANS = .FALSE.
  IF (FSORD.EQ.1) THEN
    CALL MA48CD(M,N,TRANS,JOB,LA,A,IRN1,KEEP,CNTL,ICNTL,&
              B,X,ERROR,W,IW,INFO)
  else
    CALL MA48C(M,N,TRANS,JOB,LA,A,IRN1,KEEP,CNTL,ICNTL,&
              B,X,ERROR,W,IW,INFO)
  endif
  WRITE (6,FMT='(A,I3/A)') 'INFO(1) =',INFO(1)
  !do I=1,M
  !  X(I)=X(I)*C(I)
    !write(*,"(A,i5,F20.2)") 'XD',I,X(I)
  !end do
  !do I=1,M
    !X(I)=SOL(I)
    !write(*,"(A,i5,F10.6)") 'XD',I,X(I)
  !end do
  deallocate(A,CNTL,RINFO,W,ERROR)!,RHS,SOL
  deallocate(IRN1,ICNTL,INFO,IW,KEEP)
END SUBROUTINE SPEC48_SSOL

SUBROUTINE SPEC48_SSOL2LA(INSIZE,IRN,JCN,VA,B,X)
  use constants
  IMPLICIT NONE

  integer NEFAC,JOB
  integer(4) JCN(*),IRN(*)
  integer(4) INSIZE(*)
  integer M,N,NE,T
  real(kind=DPC) VA(*),B(*),X(*)
  integer I,J,L,LA, MAXN,SGNDET!,IFAIL,LP,II
  LOGICAL TRANS,checksol
  real(kind=DPC), pointer :: CNTL(:),RINFO(:),W(:),ERROR1(:),LOGDET!A(:),,R(:),C(:)!,RHS(:),SOL(:)
  integer, pointer :: ICNTL(:),INFO(:),IW(:),KEEP(:)!,JCN1(:),IRN1(:)
  M=INSIZE(1)
  N=INSIZE(2)
  NE=INSIZE(3)
  NEFAC=INSIZE(4)
  if(NEFAC.EQ.0) then
    LA=2*NE
  else
    LA=ceiling((NEFAC/100.0)*NE)
  endif
  MAXN=N
  IF (N.LT.M) THEN
    MAXN=M
  END IF
  allocate(CNTL(10),RINFO(10),ERROR1(3))!,RHS(M),SOL(M),R(MAXN),C(MAXN),W(5*MAXN)A(LA),
  allocate(ICNTL(20),INFO(20),IW(6*M+3*N))!JCN1(LA),IRN1(LA),
  !write(*,"(A,i5,i5,F20.2)") 'Ain 1',IRN(1),JCN(1),VA(1)
  !do I=1,NE
    !IRN1(I)=IRN(I)
    !JCN1(I)=JCN(I)
    !A(I)=VA(I)
    !write(*,"(A,i5,i5,F20.2)") 'Ain',IRN(I),JCN(i),VA(I)
  !end do
  !write(*,"(A,i10,i10,i15)") 'Matsize ',M,N,NE
  !LP = 6
  !CALL MC29AD(M,N,NE,A,IRN1,JCN1,R,C,W,LP,IFAIL)
  !IF (IFAIL.LT.0) THEN
  !  WRITE (6,'(A)') ' Error in scaling routine'
  !  STOP
  !END IF
  !DO I = 1,M
  !  R(I) = EXP(R(I))
  !End Do
  !DO I = 1,N
  !  C(I) = EXP(C(I))
  !End Do
  !DO II = 1,NE
  !  I = IRN1(II)
  !  J = JCN1(II)
    !write(*,"(A,i5,i5,F10.2,F10.2,F10.2)") 'Ainbs',I,J,R(I),C(J),A(II)
  !  A(II) = A(II)*R(I)*C(J)
    !write(*,"(A,i5,i5,F20.2,F20.2)") 'Ain',I,J,R(I)*C(J),A(II)
  !End Do
!     Factorize matrix
  IF (FSORD.EQ.1) THEN
    CALL MA48ID(CNTL,ICNTL)
  else
    CALL MA48I(CNTL,ICNTL)
  endif
  !ICNTL(4)=0
  !ICNTL(6)=1
  T=M+5*N+4*N/ICNTL(6)+7
  !CNTL(4)=0.00000000000000001
  !CNTL(2)=0
  allocate(KEEP(T))

  JOB=1
  !ICNTL(8)=1
  !CNTL(2)=0.0
  IF (FSORD.EQ.1) THEN
    CALL MA48AD(M,N,NE,JOB,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  else
    CALL MA48A(M,N,NE,JOB,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  endif
  WRITE (6,FMT='(A,I3/A)') 'INFO(3) =',INFO(3)/NE
  !JOB=1
  !do I=1,NE
  !  write(*,"(A,i5,i5,i5,F20.2)") 'AAD',I,IRN1(I),JCN1(I),A(I)
  !end do
  deallocate(IW)
  allocate(W(M),IW(2*M+2*N))
  JOB=1
  !CNTL(3)=1e-5
  IF (FSORD.EQ.1) THEN
    CALL MA48BD(M,N,NE,JOB,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,&
                RINFO)
  else
    CALL MA48B(M,N,NE,JOB,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,&
                RINFO)
  endif
  WRITE (6,FMT='(A,I3/A)') 'INFO(4) =',INFO(4)/NE
  IF (INFO(1).NE.0) THEN
    WRITE (6,FMT='(A,I3/A)') 'STOP from MA48B/BD with INFO(1) =',&
    INFO(1),'Solution not possible'
    WRITE (6,FMT='(A,i10)') 'INFO(5) =',INFO(5)
    STOP
  END IF

  !do I=1,M
    !RHS(I)=B(I)
    !write(*,"(i5,F10.2)") I,B(I)
  !end do
  !do I=1,M
  !  B(I)=B(I)*R(I)
    !write(*,"(i5,F10.2)") I,B(I)
  !end do
  !deallocate(JCN1)
  JOB=1
  deallocate(W,IW)
  IF (JOB.EQ.1) THEN
    allocate(W(2*MAXN))
  else
    allocate(W(4*MAXN))
  END IF
  allocate(IW(MAXN))!A(LA),IRN1(LA),
  TRANS = .FALSE.
  IF (FSORD.EQ.1) THEN
    CALL MA48CD(M,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
              B,X,ERROR1,W,IW,INFO)
  else
    CALL MA48C(M,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
              B,X,ERROR1,W,IW,INFO)
  endif
  WRITE (6,FMT='(A,I3/A)') 'INFO(1) =',INFO(1)
  !do I=1,M
  !  X(I)=X(I)*C(I)
    !write(*,"(A,i5,F10.6)") 'XD',I,X(I)
  !end do
  !do I=1,M
    !X(I)=SOL(I)
    !write(*,"(A,i5,F10.6)") 'XD',I,X(I)
  !end do
  deallocate(CNTL,RINFO,W,ERROR1)!A,,RHS,SOL
  deallocate(ICNTL,INFO,IW,KEEP)!IRN1,
END SUBROUTINE SPEC48_SSOL2LA

SUBROUTINE SPEC48M_SSOL2LA(INSIZE,IRN,JCN,VA,B,X)
  use constants
  IMPLICIT NONE

  integer NEFAC,JOB
  integer(4) JCN(*),IRN(*)
  integer(4) INSIZE(*)
  integer M,N,NE,T
  real(kind=DPC) VA(*),B(*),X(*)
  integer I,J,L,LA, MAXN,SGNDET!,IFAIL,LP,II
  LOGICAL TRANS,checksol
  real(kind=DPC), pointer :: CNTL(:),RINFO(:),W(:),ERROR1(:),LOGDET!A(:),,R(:),C(:)!,RHS(:),SOL(:)
  integer, pointer :: ICNTL(:),INFO(:),IW(:),KEEP(:)!,JCN1(:),IRN1(:)
  M=INSIZE(1)
  N=INSIZE(2)
  NE=INSIZE(3)
  NEFAC=INSIZE(4)
  if(NEFAC.EQ.0) then
    LA=2*NE
  else
    LA=ceiling((NEFAC/100.0)*NE)
  endif
  if(LA.NE.INSIZE(6)) then
    !write(*,"(A,i15,i15)") 'Note LA!!!!',LA,INSIZE(6)
    LA=INSIZE(6)
  end if
  MAXN=N
  IF (N.LT.M) THEN
    MAXN=M
  END IF
  allocate(CNTL(10),RINFO(10),ERROR1(3))!,RHS(M),SOL(M),R(MAXN),C(MAXN),W(5*MAXN)A(LA),
  allocate(ICNTL(20),INFO(20),IW(6*M+3*N))!JCN1(LA),IRN1(LA),
  !write(*,"(A,i5,i5,F20.2)") 'Ain 1',IRN(1),JCN(1),VA(1)
  !do I=1,NE
    !IRN1(I)=IRN(I)
    !JCN1(I)=JCN(I)
    !A(I)=VA(I)
    !write(*,"(A,i5,i5,F20.2)") 'Ain',IRN(I),JCN(i),VA(I)
  !end do
  !write(*,"(A,i10,i10,i15)") 'Matsize ',M,N,NE
  !LP = 6
  !CALL MC29AD(M,N,NE,A,IRN1,JCN1,R,C,W,LP,IFAIL)
  !IF (IFAIL.LT.0) THEN
  !  WRITE (6,'(A)') ' Error in scaling routine'
  !  STOP
  !END IF
  !DO I = 1,M
  !  R(I) = EXP(R(I))
  !End Do
  !DO I = 1,N
  !  C(I) = EXP(C(I))
  !End Do
  !DO II = 1,NE
  !  I = IRN1(II)
  !  J = JCN1(II)
    !write(*,"(A,i5,i5,F10.2,F10.2,F10.2)") 'Ainbs',I,J,R(I),C(J),A(II)
  !  A(II) = A(II)*R(I)*C(J)
    !write(*,"(A,i5,i5,F20.2,F20.2)") 'Ain',I,J,R(I)*C(J),A(II)
  !End Do
!     Factorize matrix
  IF (FSORD.EQ.1) THEN
    CALL MA48ID(CNTL,ICNTL)
  else
    CALL MA48I(CNTL,ICNTL)
  endif
  !ICNTL(4)=0
  !ICNTL(6)=1
  T=M+5*N+4*N/ICNTL(6)+7
  !CNTL(4)=0.00000000000000001
  !CNTL(2)=0
  allocate(KEEP(T))

  JOB=1
  !ICNTL(8)=1
  !CNTL(2)=0.0
  IF (FSORD.EQ.1) THEN
    CALL MA48AD(M,N,NE,JOB,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  else
    CALL MA48A(M,N,NE,JOB,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  endif
  WRITE (6,FMT='(A,I3/A)') 'INFO(3) =',INFO(3)/NE
  !JOB=1
  !do I=1,NE
  !  write(*,"(A,i5,i5,i5,F20.2)") 'AAD',I,IRN1(I),JCN1(I),A(I)
  !end do
  deallocate(IW)
  allocate(W(M),IW(2*M+2*N))
  JOB=1
  !CNTL(3)=1e-5
  IF (FSORD.EQ.1) THEN
    CALL MA48BD(M,N,NE,JOB,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,&
                RINFO)
  else
    CALL MA48B(M,N,NE,JOB,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,&
                RINFO)
  endif
  WRITE (6,FMT='(A,I3/A)') 'INFO(4) =',INFO(4)/NE
  IF (INFO(1).NE.0) THEN
    WRITE (6,FMT='(A,I3/A)') 'STOP from MA48B/BD with INFO(1) =',&
    INFO(1),'Solution not possible'
    WRITE (6,FMT='(A,i10)') 'INFO(5) =',INFO(5)
    STOP
  END IF

  !do I=1,M
    !RHS(I)=B(I)
    !write(*,"(i5,F10.2)") I,B(I)
  !end do
  !do I=1,M
  !  B(I)=B(I)*R(I)
    !write(*,"(i5,F10.2)") I,B(I)
  !end do
  !deallocate(JCN1)
  JOB=1
  deallocate(W,IW)
  IF (JOB.EQ.1) THEN
    allocate(W(2*MAXN))
  else
    allocate(W(4*MAXN))
  END IF
  allocate(IW(MAXN))!A(LA),IRN1(LA),
  TRANS = .FALSE.
  IF (FSORD.EQ.1) THEN
    CALL MA48CD(M,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
              B,X,ERROR1,W,IW,INFO)
  else
    CALL MA48C(M,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
              B,X,ERROR1,W,IW,INFO)
  endif
  WRITE (6,FMT='(A,I3/A)') 'INFO(1) =',INFO(1)
  !do I=1,M
  !  X(I)=X(I)*C(I)
    !write(*,"(A,i5,F10.6)") 'XD',I,X(I)
  !end do
  !do I=1,M
    !X(I)=SOL(I)
    !write(*,"(A,i5,F10.6)") 'XD',I,X(I)
  !end do
  deallocate(CNTL,RINFO,W,ERROR1)!A,,RHS,SOL
  deallocate(ICNTL,INFO,IW,KEEP)!IRN1,
END SUBROUTINE SPEC48M_SSOL2LA

!SUBROUTINE SPEC48_SSOL0(INSIZE,IRN,JCN,VA,B,X)
!  use constants
!  IMPLICIT NONE
!
!  INTEGER, PARAMETER :: myreal = SELECTED_REAL_KIND(8)!kind(real*8)!KIND( 0.0D0 )
!  INTEGER, PARAMETER :: myint = selected_int_kind(16)!kind(interger*4)!KIND( 1)
!  integer NEFAC,JOB
!  integer (kind=myint) INSIZE(*),JCN(*),IRN(*)
!  integer M,N,NE,T
!  real (kind=DPC) VA(*),B(*),X(*)
!  integer I,J,L,LA, MAXN,IFAIL,LP,II,SGNDET
!  LOGICAL TRANS,checksol
!  DOUBLE PRECISION, pointer :: A(:),CNTL(:),RINFO(:),W(:),ERROR(:),RHS(:),SOL(:),R(:),C(:),LOGDET
!  integer, pointer :: ICNTL(:),INFO(:),IW(:),KEEP(:),JCN1(:),IRN1(:)
!  M=int(INSIZE(1))
!  N=int(INSIZE(2))
!  NE=int(INSIZE(3))
!  !NEFAC=int(INSIZE(4))
!  if(NEFAC.EQ.0) then
!    LA=3*NE
!  else
!    LA=ceiling((NEFAC/100.0)*NE)
!  endif
!  MAXN=N
!  IF (N.LT.M) THEN
!    MAXN=M
!  END IF
!  allocate(A(LA),CNTL(10),RINFO(10),W(5*MAXN),ERROR(3),RHS(M),SOL(M),R(MAXN),C(MAXN))
!  allocate(JCN1(LA),IRN1(LA),ICNTL(20),INFO(20),IW(6*M+3*N))
!
!  do I=1,NE
!    IRN1(I)=int(IRN(I))
!    JCN1(I)=int(JCN(I))
!    A(I)=dble(VA(I))
!    !write(*,"(A,i5,i5,F20.2)") 'Ain',IRN1(I),JCN1(i),VA(I)
!  end do
!  write(*,"(A,i5,i5)") 'Matsize ',M,N
!! Scale matrix
!  LP = 6
!  CALL MC29AD(M,N,NE,A,IRN1,JCN1,R,C,W,LP,IFAIL)
!  IF (IFAIL.LT.0) THEN
!    WRITE (6,'(A)') ' Error in scaling routine'
!    STOP
!  END IF
!  DO I = 1,M
!    R(I) = EXP(R(I))
!  End Do
!  DO I = 1,N
!    C(I) = EXP(C(I))
!  End Do
!  DO II = 1,NE
!    I = IRN1(II)
!    J = JCN1(II)
!    !write(*,"(A,i5,i5,F10.2,F10.2,F10.2)") 'Ainbs',I,J,R(I),C(J),A(II)
!    A(II) = A(II)*R(I)*C(J)
!    !write(*,"(A,i5,i5,F20.2,F20.2)") 'Ain',I,J,R(I)*C(J),A(II)
!  End Do
!!     Factorize matrix
!  CALL MA48ID(CNTL,ICNTL)
!  !ICNTL(4)=0
!  !ICNTL(6)=1
!  T=M+5*N+4*N/ICNTL(6)+7
!  !CNTL(4)=0.00000000000000001
!  !CNTL(2)=0
!  allocate(KEEP(T))
!
!  JOB=1
!  CALL MA48AD(M,N,NE,1,LA,A,IRN1,JCN1,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
!  WRITE (6,FMT='(A,I3/A)') 'INFO(3) =',INFO(3)/NE
!  JOB=1
!  !do I=1,NE
!  !  write(*,"(A,i5,i5,i5,F20.2)") 'AAD',I,IRN1(I),JCN1(I),A(I)
!  !end do
!  CALL MA48BD(M,N,NE,1,LA,A,IRN1,JCN1,KEEP,CNTL,ICNTL,W,IW,INFO,&
!                RINFO)
!  WRITE (6,FMT='(A,I3/A)') 'INFO(4) =',INFO(4)/NE
!  IF (INFO(1).NE.0) THEN
!    WRITE (6,FMT='(A,I3/A)') 'STOP from MA48B/BD with INFO(1) =',&
!    INFO(1),'Solution not possible'
!    STOP
!  END IF
!
!  do I=1,M
!    RHS(I)=B(I)*R(I)
!    !write(*,"(i5,F10.2)") I,B(I)
!  end do
!  JOB=1
!  TRANS = .FALSE.
!  CALL MA48CD(M,N,TRANS,JOB,LA,A,IRN1,KEEP,CNTL,ICNTL,&
!              RHS,SOL,ERROR,W,IW,INFO)
!  WRITE (6,FMT='(A,I3/A)') 'INFO(1) =',INFO(1)
!  do I=1,M
!    X(I)=SOL(I)*C(I)
!    !write(*,"(A,i5,F20.2)") 'XD',I,X(I)
!  end do
!  deallocate(A,CNTL,RINFO,W,ERROR,RHS,SOL,R,C)
!  deallocate(IRN1,JCN1,ICNTL,INFO,IW,KEEP)
!END SUBROUTINE SPEC48_SSOL0
!
!SUBROUTINE SPEC48_SSOL2(INSIZE,IRN,JCN,VA,B,X)
!  use constants
!  IMPLICIT NONE
!
!  INTEGER, PARAMETER :: myreal = SELECTED_REAL_KIND(8)!kind(real*8)!KIND( 0.0D0 )
!  INTEGER, PARAMETER :: myint = selected_int_kind(16)!kind(interger*4)!KIND( 1)
!  integer M,N,NE,T,NEFAC,JOB,J
!  integer (kind=myint) INSIZE(*),JCN(*),IRN(*)
!  real (kind=DPC) VA(*),B(*),X(*)
!  integer I,LA, MAXN,LP,II,IFAIL
!  LOGICAL TRANS
!  DOUBLE PRECISION, pointer :: A(:),CNTL(:),RINFO(:),W(:),ERROR(:),RHS(:),SOL(:),R(:),C(:)
!  integer, pointer :: JCN1(:),IRN1(:),ICNTL(:),INFO(:),IW(:),KEEP(:)
!  M=int(INSIZE(1))
!  N=int(INSIZE(2))
!  NE=int(INSIZE(3))
!  !NEFAC=int(INSIZE(4))
!  !write(*,"(i5,i5)") IRN(1),JCN(1)
!  !write(*,"(i5,i5)") IRN(2),JCN(2)
!  !write(*,"(i5,i5)") IRN(3),JCN(3)
!  !write(*,"(i5,i5)") IRN(4),JCN(4)
!  if(NEFAC.EQ.0) then
!    LA=3*NE
!  else
!    LA=ceiling((NEFAC/100.0)*NE)
!  endif
!  MAXN=N
!  IF (N.LT.M) THEN
!    MAXN=M
!  END IF
!  allocate(A(LA),CNTL(10),RINFO(10),W(5*MAXN),ERROR(3),RHS(M),SOL(M),R(MAXN),C(MAXN))
!!  real (kind=myreal) A(LA),CNTL(10),LOGDET,RINFO(10),W(MAXN4)
!  allocate(JCN1(LA),IRN1(LA),ICNTL(20),INFO(20),IW(6*M+3*N))
!  do I=1,NE
!    IRN1(I)=IRN(I)
!    JCN1(I)=JCN(I)
!    A(I)=VA(I)
!    !write(*,"(i5,i5,i5,F20.2)") I,IRN1(I),JCN1(I),A(I)
!  end do
!! Scale matrix
!  LP = 6
!  CALL MC29AD(M,N,NE,A,IRN1,JCN1,R,C,W,LP,IFAIL)
!  IF (IFAIL.LT.0) THEN
!    WRITE (6,'(A)') ' Error in scaling routine'
!    STOP
!  END IF
!  DO I = 1,M
!    R(I) = 1!EXP(R(I))
!  End Do
!  DO I = 1,N
!    C(I) = 1!EXP(C(I))
!  End Do
!  DO II = 1,NE
!    I = IRN1(II)
!    J = JCN1(II)
!    !write(*,"(A,i5,i5,F10.2,F10.2,F10.2)") 'Ainbs',I,J,R(I),C(J),A(II)
!    A(II) = A(II)*R(I)*C(J)
!    write(*,"(A,i5,i5,F20.2,F20.2)") 'Ain',I,J,R(I)*C(J),A(II)
!  End Do
!
!  do I=1,M
!    RHS(I)=B(I)*R(I)
!    write(*,"(i5,F20.2)") I,RHS(I)
!  end do
!!     Factorize matrix
!  CALL MA48ID(CNTL,ICNTL)
!  ICNTL(6)=1
!  !CNTL(2)=0.5
!  !CNTL(1)=0.0
!  T=M+5*N+4*N/ICNTL(6)+7
!  allocate(KEEP(T))
!  JOB=1
!  CALL MA48AD(M,N,NE,JOB,LA,A,IRN1,JCN1,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
!  CALL MA48BD(M,N,NE,JOB,LA,A,IRN1,JCN1,KEEP,CNTL,ICNTL,W,IW,INFO,RINFO)
!  WRITE (6,FMT='(A,I3/A)') 'INFO(1) =',INFO(1)
!  write(*,"(A,i5)") 'INFO(4) =',INFO(4)/NE
!  IF (INFO(1).NE.0) THEN
!    WRITE (6,FMT='(A,I3/A)') 'STOP from MA48B/BD with INFO(1) =',&
!    INFO(1),'Solution not possible'
!    STOP
!  END IF
!  TRANS = .FALSE.
!  !JOB=1
!  CALL MA48CD(M,N,TRANS,JOB,LA,A,IRN,KEEP,CNTL,ICNTL,&
!              RHS,SOL,ERROR,W,IW,INFO)
!  do I=1,N
!    X(I)=SOL(I)*C(I)
!    write(*,"(A,i5,F20.2)") 'SOL ',I,SOL(I)
!  end do
!  !deallocate(A,CNTL,RINFO,W,ERROR)
!  deallocate(CNTL,RINFO,ERROR,A,R,C)
!  deallocate(IRN1,JCN1,ICNTL,INFO,IW,KEEP)
!  deallocate(W)
!END SUBROUTINE SPEC48_SSOL2
!
!SUBROUTINE SPEC48_SSOL1(INSIZE,IRN,JCN,VA,B,X)
!  use constants
!  IMPLICIT NONE
!
!  INTEGER, PARAMETER :: myreal = SELECTED_REAL_KIND(8)!kind(real*8)!KIND( 0.0D0 )
!  INTEGER, PARAMETER :: myint = selected_int_kind(16)!kind(interger*4)!KIND( 1)
!  integer M,N,NE,T,NEFAC,JOB
!  integer (kind=myint) INSIZE(*),JCN(*),IRN(*)
!  real (kind=DPC) VA(*),B(*),X(*)
!  integer I,LA, MAXN
!  LOGICAL TRANS
!  DOUBLE PRECISION, pointer :: A(:),CNTL(:),RINFO(:),W(:),ERROR(:),RHS(:),SOL(:)
!  integer, pointer :: JCN1(:),IRN1(:),ICNTL(:),INFO(:),IW(:),KEEP(:)
!  M=INSIZE(1)
!  N=INSIZE(2)
!  NE=INSIZE(3)
!  NEFAC=INSIZE(4)
!  !write(*,"(i5,i5)") IRN(1),JCN(1)
!  !write(*,"(i5,i5)") IRN(2),JCN(2)
!  !write(*,"(i5,i5)") IRN(3),JCN(3)
!  !write(*,"(i5,i5)") IRN(4),JCN(4)
!  if(NEFAC.EQ.0) then
!    LA=3*NE
!  else
!    LA=ceiling((NEFAC/100.0)*NE)
!  endif
!  MAXN=N
!  IF (N.LT.M) THEN
!    MAXN=M
!  END IF
!  allocate(A(LA),CNTL(10),RINFO(10),W(5*MAXN),ERROR(3),RHS(M),SOL(M))
!
!!  real (kind=myreal) A(LA),CNTL(10),LOGDET,RINFO(10),W(MAXN4)
!  allocate(JCN1(LA),IRN1(LA),ICNTL(20),INFO(20),IW(6*M+3*N))
!  do I=1,NE
!    IRN1(I)=IRN(I)
!    JCN1(I)=JCN(I)
!    A(I)=VA(I)
!    write(*,"(i5,i5,i5,F20.2)") I,IRN1(I),JCN1(I),A(I)
!  end do
!  !B(8)=10
!  do I=1,M
!    RHS(I)=B(I)
!    write(*,"(i5,F20.2)") I,B(I)
!  end do
!!     Factorize matrix
!  CALL MA48ID(CNTL,ICNTL)
!  ICNTL(6)=N
!  !CNTL(2)=0.5
!  !CNTL(1)=0.0
!  T=M+5*N+4*N/ICNTL(6)+7
!  allocate(KEEP(T))
!  JOB=1
!  CALL MA48AD(M,N,NE,JOB,LA,A,IRN1,JCN1,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
!  CALL MA48BD(M,N,NE,JOB,LA,A,IRN1,JCN1,KEEP,CNTL,ICNTL,W,IW,INFO,RINFO)
!  WRITE (6,FMT='(A,I3/A)') 'INFO(1) =',INFO(1)
!  write(*,"(A,i5)") 'INFO(4) =',INFO(4)/NE
!  IF (INFO(1).NE.0) THEN
!    WRITE (6,FMT='(A,I3/A)') 'STOP from MA48B/BD with INFO(1) =',&
!    INFO(1),'Solution not possible'
!    STOP
!  END IF
!  TRANS = .FALSE.
!  !JOB=1
!  CALL MA48CD(N,N,TRANS,JOB,LA,A,IRN,KEEP,CNTL,ICNTL,&
!              RHS,SOL,ERROR,W,IW,INFO)
!  do I=1,N
!    X(I)=SOL(I)
!    write(*,"(i5,F20.2)") I,X(I)
!  end do
!  !deallocate(A,CNTL,RINFO,W,ERROR)
!  deallocate(CNTL,RINFO,ERROR,A)
!  deallocate(IRN1,JCN1,ICNTL,INFO,IW,KEEP)
!  deallocate(W)
!END SUBROUTINE SPEC48_SSOL1

SUBROUTINE MY_SPAR_ADD(VECBIVI1,BIVIINDX1,NZ1,VECBIVI0,BIVIINDX0,NZ0,NZ2)
  use constants
  IMPLICIT NONE
  real(kind=DPC) VECBIVI1(*),VECBIVI0(*)
  INTEGER NZ1(*),NZ0(*),NZ2(*)
  INTEGER BIVIINDX0(0:NZ0(1)),BIVIINDX1(0:NZ1(1))
  INTEGER j,j1,j1p,i,i1,j2,j3,j2p,ip,b1,e1,b2,e2,b3,e3,buffer,m,l
  buffer=100000!to avoid crash in fotran assignment
  j=NZ1(1)
  j1=NZ0(1)
  i=NZ2(1)
  do while(i.GT.0)!i=NZ2(1),1,-1
    i1=i
    j2=j
    j3=j1
    do while(BIVIINDX1(j).EQ.BIVIINDX0(j1).AND.BIVIINDX1(j).GT.-1)
          !VECBIVI1(i)=VECBIVI1(j)+VECBIVI0(j1)
          !BIVIINDX1(i)=BIVIINDX1(j)
          j1=j1-1
          j=j-1
          !i=i-1
          !GO TO 20
    end do
    i=i1-(j2-j)
    if(i.NE.i1) then
    ip=i+1
    j1p=j+1
    j2p=j1+1
    if(ip.GT.j2) then
      if(i1-ip.LT.buffer) then
        VECBIVI1(ip:i1)=VECBIVI1(j1p:j2)+VECBIVI0(j2p:j3)
        BIVIINDX1(ip:i1)=BIVIINDX1(j1p:j2)
      else
        l=(i1-ip)/buffer
        do m=0,l-1,1
          b1=ip+buffer*m
          e1=b1+buffer-1
          b2=j1p+buffer*m
          e2=b2+buffer-1
          b3=j2p+buffer*m
          e3=b3+buffer-1
          VECBIVI1(b1:e1)=VECBIVI1(b2:e2)+VECBIVI0(b3:e3)
          BIVIINDX1(b1:e1)=BIVIINDX1(b2:e2)
        end do
          VECBIVI1(e1+1:i1)=VECBIVI1(e2+1:j2)+VECBIVI0(e3+1:j3)
          BIVIINDX1(e1+1:i1)=BIVIINDX1(e2+1:j2)
      end if
    else
!      do m=j2,ip,-1
!        VECBIVI1(i1-j2+m)=VECBIVI1(m)+VECBIVI0(j3-j2+m)
!      end do
      if(j2-ip.LT.buffer) then
        VECBIVI1(i1-j2+ip:i1)=VECBIVI1(ip:j2)+VECBIVI0(j3-j2+ip:j3)
        BIVIINDX1(i1-j2+ip:i1)=BIVIINDX1(ip:j2)
      else
        l=(j2-ip)/buffer
        e1=i1-j2+ip+buffer*l
        e2=ip+buffer*l
        e3=j3-j2+ip+buffer*l
        VECBIVI1(e1:i1)=VECBIVI1(e2:j2)+VECBIVI0(e3:j3)
        BIVIINDX1(e1:i1)=BIVIINDX1(e2:j2)
        do m=l-1,0,-1
          b1=e1-buffer
          b2=e2-buffer
          b3=e3-buffer
          VECBIVI1(b1:e1-1)=VECBIVI1(b2:e2-1)+VECBIVI0(b3:e3-1)
          BIVIINDX1(b1:e1-1)=BIVIINDX1(b2:e2-1)
          e1=b1
          e2=b2
          e3=b3
        end do
      end if
      if(i-j1p.LT.buffer) then
        VECBIVI1(ip:ip+i-j1p)=VECBIVI1(j1p:i)+VECBIVI0(j2p:j3-j2+i)
        BIVIINDX1(ip:i1-j2+i)=BIVIINDX1(j1p:i)
      else
!      do m=i,j1p,-1
!        BIVIINDX1(i1-j2+m)=BIVIINDX1(m)
!      end do
!      do m=i,j1p,-1
!        VECBIVI1(i1-j2+m)=VECBIVI1(m)+VECBIVI0(j3-j2+m)
!      end do
        l=(i-j1p)/buffer
        e1=ip+buffer*l
        e2=j1p+buffer*l
        e3=j2p+buffer*l
        VECBIVI1(e1:ip+i-j1p)=VECBIVI1(e2:i)+VECBIVI0(e3:j3-j2+i)
        BIVIINDX1(e1:i1-j2+i)=BIVIINDX1(e2:i)
        do m=l-1,0,-1
          b1=e1-buffer
          b2=e2-buffer
          b3=e3-buffer
          !WRITE (*,"(i10,i10,i10,i10)"),b1,e1,b2,e2
          VECBIVI1(b1:e1-1)=VECBIVI1(b2:e2-1)+VECBIVI0(b3:e3-1)
          BIVIINDX1(b1:e1-1)=BIVIINDX1(b2:e2-1)
          e1=b1
          e2=b2
          e3=b3
        end do
      end if
      end if
    end if
    i1=i
    j2=j
    do while(BIVIINDX1(j).GT.BIVIINDX0(j1))
          !VECBIVI1(i)=VECBIVI1(j)
          !BIVIINDX1(i)=BIVIINDX1(j)
          j=j-1
          !i=i-1
          !GO TO 20
    end do
    i=i1-(j2-j)
    if(i.NE.i1) then
    ip=i+1
    j1p=j+1
    if(ip.GT.j2) then
      if(i1-ip.LT.buffer) then
        VECBIVI1(ip:i1)=VECBIVI1(j1p:j2)
        BIVIINDX1(ip:i1)=BIVIINDX1(j1p:j2)
      else
        l=(i1-ip)/buffer
        do m=0,l-1,1
          b1=ip+buffer*m
          e1=b1+buffer-1
          b2=j1p+buffer*m
          e2=b2+buffer-1
          VECBIVI1(b1:e1)=VECBIVI1(b2:e2)
          BIVIINDX1(b1:e1)=BIVIINDX1(b2:e2)
        end do
          VECBIVI1(e1+1:i1)=VECBIVI1(e2+1:j2)
          BIVIINDX1(e1+1:i1)=BIVIINDX1(e2+1:j2)
      end if
    else
      if(j2-ip.LT.buffer) then
        VECBIVI1(i1-j2+ip:i1)=VECBIVI1(ip:j2)
        BIVIINDX1(i1-j2+ip:i1)=BIVIINDX1(ip:j2)
      else
        l=(j2-ip)/buffer
        e1=i1-j2+ip+buffer*l
        e2=ip+buffer*l
        VECBIVI1(e1:i1)=VECBIVI1(e2:j2)
        BIVIINDX1(e1:i1)=BIVIINDX1(e2:j2)
        do m=l-1,0,-1
          b1=e1-buffer
          b2=e2-buffer
          VECBIVI1(b1:e1-1)=VECBIVI1(b2:e2-1)
          BIVIINDX1(b1:e1-1)=BIVIINDX1(b2:e2-1)
          e1=b1
          e2=b2
        end do
      end if
      if(i-j1p.LT.buffer) then
        VECBIVI1(ip:ip+i-j1p)=VECBIVI1(j1p:i)
        BIVIINDX1(ip:ip+i-j1p)=BIVIINDX1(j1p:i)
      else
        l=(i-j1p)/buffer
        e1=ip+buffer*l
        e2=j1p+buffer*l
        VECBIVI1(e1:ip+i-j1p)=VECBIVI1(e2:i)
        BIVIINDX1(e1:ip+i-j1p)=BIVIINDX1(e2:i)
        do m=l-1,0,-1
          b1=e1-buffer
          b2=e2-buffer
          !WRITE (*,"(i10,i10,i10,i10)"),b1,e1,b2,e2
          VECBIVI1(b1:e1-1)=VECBIVI1(b2:e2-1)
          BIVIINDX1(b1:e1-1)=BIVIINDX1(b2:e2-1)
          e1=b1
          e2=b2
        end do
      end if

!    if(ip.GT.j2) then
!      if(i1-ip.LT.buffer) then
!        VECBIVI1(ip:i1)=VECBIVI1(j1p:j2)
!        BIVIINDX1(ip:i1)=BIVIINDX1(j1p:j2)
!      else
!        l=(i1-ip)/buffer
!        !WRITE (*,"(i10,i10,i10,i10,i10)"),l,ip,i1,j1p,j2
!        do m=0,l-1,1
!          b1=ip+buffer*m
!          e1=b1+buffer-1
!          b2=j1p+buffer*m
!          e2=b2+buffer-1
!          !WRITE (*,"(i10,i10,i10,i10,i10)"),m,b1,e1,b2,e2
!          VECBIVI1(b1:e1)=VECBIVI1(b2:e2)
!          BIVIINDX1(b1:e1)=BIVIINDX1(b2:e2)
!        end do
!          VECBIVI1(e1+1:i1)=VECBIVI1(e2+1:j2)
!          BIVIINDX1(e1+1:i1)=BIVIINDX1(e2+1:j2)
!      end if
!    else
!      if(j2-ip.LT.buffer) then
!        VECBIVI1(i1-j2+ip:i1)=VECBIVI1(ip:j2)
!        BIVIINDX1(i1-j2+ip:i1)=BIVIINDX1(ip:j2)
!      else
!        l=(j2-ip)/buffer
!        !WRITE (*,"(i10,i10,i10,i10,i10)"),l,ip,i1,j1p,j2
!        do m=0,l-1,1
!          b1=i1-j2+ip+buffer*m
!          e1=b1+buffer-1
!          b2=ip+buffer*m
!          e2=b2+buffer-1
!          !WRITE (*,"(i10,i10,i10,i10,i10)"),m,b1,e1,b2,e2
!          VECBIVI1(b1:e1)=VECBIVI1(b2:e2)
!          BIVIINDX1(b1:e1)=BIVIINDX1(b2:e2)
!        end do
!          VECBIVI1(e1+1:i1)=VECBIVI1(e2+1:j2)
!          BIVIINDX1(e1+1:i1)=BIVIINDX1(e2+1:j2)
!      end if
!      if(i-j1p.LT.buffer) then
!        VECBIVI1(ip:ip+i-j1p)=VECBIVI1(j1p:i)
!        BIVIINDX1(ip:ip+i-j1p)=BIVIINDX1(j1p:i)
!      else
!        l=(i-j1p)/buffer
!        !WRITE (*,"(i10,i10,i10,i10,i10)"),l,ip,i1,j1p,j2
!        do m=0,l-1,1
!          b1=ip+buffer*m
!          e1=b1+buffer-1
!          b2=j1p+buffer*m
!          e2=b2+buffer-1
!          !WRITE (*,"(i10,i10,i10,i10,i10)"),m,b1,e1,b2,e2
!          VECBIVI1(b1:e1)=VECBIVI1(b2:e2)
!          BIVIINDX1(b1:e1)=BIVIINDX1(b2:e2)
!        end do
!          VECBIVI1(e1+1:ip+i-j1p)=VECBIVI1(e2+1:i)
!          BIVIINDX1(e1+1:ip+i-j1p)=BIVIINDX1(e2+1:i)
!      end if
    end if
    end if
    i1=i
    j2=j1
    do while (BIVIINDX1(j).LT.BIVIINDX0(j1))
    !VECBIVI1(i)=VECBIVI0(j1)
    !BIVIINDX1(i)=BIVIINDX0(j1)
    j1=j1-1
    !i=i-1
    end do
    i=i1-(j2-j1)
    if(i.NE.i1) then
    ip=i+1
    j1p=j1+1
      if(i1-ip.LT.buffer) then
        VECBIVI1(ip:i1)=VECBIVI0(j1p:j2)
        BIVIINDX1(ip:i1)=BIVIINDX0(j1p:j2)
      else
        l=(i1-ip)/buffer
        !WRITE (*,"(i10,i10,i10,i10,i10)"),l,ip,i1,j1p,j2
        do m=0,l-1,1
          b1=ip+buffer*m
          e1=b1+buffer-1
          b2=j1p+buffer*m
          e2=b2+buffer-1
          !WRITE (*,"(i10,i10,i10,i10,i10)"),m,b1,e1,b2,e2
          VECBIVI1(b1:e1)=VECBIVI0(b2:e2)
          BIVIINDX1(b1:e1)=BIVIINDX0(b2:e2)
        end do
          VECBIVI1(e1+1:i1)=VECBIVI0(e2+1:j2)
          BIVIINDX1(e1+1:i1)=BIVIINDX0(e2+1:j2)
      end if
    end if
    !20 CONTINUE
  end do
END SUBROUTINE MY_SPAR_ADD

SUBROUTINE MY_SPAR_ADD3L(VECBIVI1,BIVIINDX1,NZ1,VECBIVI0,BIVIINDX0,NZ0,NZ2)
  use constants
  IMPLICIT NONE
  real(kind=DPC) VECBIVI1(*),VECBIVI0(*)
  INTEGER(8) NZ1(*),NZ2(*)
  INTEGER(8) NZ0(*)
  INTEGER(8) BIVIINDX0(0:NZ0(1)),BIVIINDX1(0:NZ2(1))
  INTEGER(8) j,i,l
  j=NZ0(1)
  i=NZ1(1)
  do l=NZ2(1),1,-1
    if(BIVIINDX1(i).LT.BIVIINDX0(j)) then
      VECBIVI1(l)=VECBIVI0(j)
      BIVIINDX1(l)=BIVIINDX0(j)
      j=j-1
    else if(BIVIINDX1(i).EQ.BIVIINDX0(j)) then
      VECBIVI1(l)=VECBIVI1(i)+VECBIVI0(j)
      BIVIINDX1(l)=BIVIINDX0(j)
      j=j-1
      i=i-1
    else
      VECBIVI1(l)=VECBIVI1(i)
      BIVIINDX1(l)=BIVIINDX1(i)
      i=i-1
    end if
    !WRITE (*,"(A10,i10,i10,i10,i20,i20)")'IJJJJ',i,j,l,BIVIINDX1(i),BIVIINDX0(j)
  end do
  !WRITE (*,"(A10,i10,i10)")'IJJJJ',i,j
END SUBROUTINE MY_SPAR_ADD3L

SUBROUTINE MY_SPAR_ADDL(VECBIVI1,BIVIINDX1,NZ1,VECBIVI0,BIVIINDX0,NZ0,NZ2)
  use constants
  IMPLICIT NONE
  real(kind=DPC) VECBIVI1(*),VECBIVI0(*)
  INTEGER NZ1(*),NZ0(*),NZ2(*)
  INTEGER(8) BIVIINDX0(0:NZ0(1)),BIVIINDX1(0:NZ1(1))
  INTEGER j,j1,j1p,i,i1,j2,j3,j2p,ip,b1,e1,b2,e2,b3,e3,buffer,m,l
  buffer=100000!to avoid crash in fotran assignment
  j=NZ1(1)
  j1=NZ0(1)
  i=NZ2(1)
  do while(i.GT.0)!i=NZ2(1),1,-1
    i1=i
    j2=j
    j3=j1
    !WRITE (*,"(i10,i10,i10,i10,i10,i20,i20)"),i,i1,j2,j3,j1,BIVIINDX1(j),BIVIINDX0(j1)
    do while(BIVIINDX1(j).EQ.BIVIINDX0(j1).AND.BIVIINDX1(j).GT.-1)
          !VECBIVI1(i)=VECBIVI1(j)+VECBIVI0(j1)
          !BIVIINDX1(i)=BIVIINDX1(j)
          j1=j1-1
          j=j-1
          !i=i-1
          !GO TO 20
    end do
    !WRITE (*,"(i10,i10,i10,i10,i10,i20,i20)"),i,i1,j2,j3,j1,BIVIINDX1(j),BIVIINDX0(j1)
    i=i1-(j2-j)
    if(i.NE.i1) then
    ip=i+1
    j1p=j+1
    j2p=j1+1
    if(ip.GT.j2) then
      if(i1-ip.LT.buffer) then
        VECBIVI1(ip:i1)=VECBIVI1(j1p:j2)+VECBIVI0(j2p:j3)
        BIVIINDX1(ip:i1)=BIVIINDX1(j1p:j2)
      else
        l=(i1-ip)/buffer
        do m=0,l-1,1
          b1=ip+buffer*m
          e1=b1+buffer-1
          b2=j1p+buffer*m
          e2=b2+buffer-1
          b3=j2p+buffer*m
          e3=b3+buffer-1
          VECBIVI1(b1:e1)=VECBIVI1(b2:e2)+VECBIVI0(b3:e3)
          BIVIINDX1(b1:e1)=BIVIINDX1(b2:e2)
        end do
          VECBIVI1(e1+1:i1)=VECBIVI1(e2+1:j2)+VECBIVI0(e3+1:j3)
          BIVIINDX1(e1+1:i1)=BIVIINDX1(e2+1:j2)
      end if
    else
!      do m=j2,ip,-1
!        VECBIVI1(i1-j2+m)=VECBIVI1(m)+VECBIVI0(j3-j2+m)
!      end do
      if(j2-ip.LT.buffer) then
        VECBIVI1(i1-j2+ip:i1)=VECBIVI1(ip:j2)+VECBIVI0(j3-j2+ip:j3)
        BIVIINDX1(i1-j2+ip:i1)=BIVIINDX1(ip:j2)
      else
        l=(j2-ip)/buffer
        e1=i1-j2+ip+buffer*l
        e2=ip+buffer*l
        e3=j3-j2+ip+buffer*l
        VECBIVI1(e1:i1)=VECBIVI1(e2:j2)+VECBIVI0(e3:j3)
        BIVIINDX1(e1:i1)=BIVIINDX1(e2:j2)
        do m=l-1,0,-1
          b1=e1-buffer
          b2=e2-buffer
          b3=e3-buffer
          VECBIVI1(b1:e1-1)=VECBIVI1(b2:e2-1)+VECBIVI0(b3:e3-1)
          BIVIINDX1(b1:e1-1)=BIVIINDX1(b2:e2-1)
          e1=b1
          e2=b2
          e3=b3
        end do
      end if
      if(i-j1p.LT.buffer) then
        VECBIVI1(ip:ip+i-j1p)=VECBIVI1(j1p:i)+VECBIVI0(j2p:j3-j2+i)
        BIVIINDX1(ip:i1-j2+i)=BIVIINDX1(j1p:i)
      else
!      do m=i,j1p,-1
!        BIVIINDX1(i1-j2+m)=BIVIINDX1(m)
!      end do
!      do m=i,j1p,-1
!        VECBIVI1(i1-j2+m)=VECBIVI1(m)+VECBIVI0(j3-j2+m)
!      end do
        l=(i-j1p)/buffer
        e1=ip+buffer*l
        e2=j1p+buffer*l
        e3=j2p+buffer*l
        VECBIVI1(e1:ip+i-j1p)=VECBIVI1(e2:i)+VECBIVI0(e3:j3-j2+i)
        BIVIINDX1(e1:i1-j2+i)=BIVIINDX1(e2:i)
        do m=l-1,0,-1
          b1=e1-buffer
          b2=e2-buffer
          b3=e3-buffer
          !WRITE (*,"(i10,i10,i10,i10)"),b1,e1,b2,e2
          VECBIVI1(b1:e1-1)=VECBIVI1(b2:e2-1)+VECBIVI0(b3:e3-1)
          BIVIINDX1(b1:e1-1)=BIVIINDX1(b2:e2-1)
          e1=b1
          e2=b2
          e3=b3
        end do
      end if
      end if
    end if
    i1=i
    j2=j
    do while(BIVIINDX1(j).GT.BIVIINDX0(j1))
          !VECBIVI1(i)=VECBIVI1(j)
          !BIVIINDX1(i)=BIVIINDX1(j)
          j=j-1
          !i=i-1
          !GO TO 20
    end do
    !WRITE (*,"(i10,i10,i10,i10,i10,i10,i20,i20)"),i,i1,j2,j3,j1,j,BIVIINDX1(j),BIVIINDX0(j1)
    i=i1-(j2-j)
    if(i.NE.i1) then
    ip=i+1
    j1p=j+1
    if(ip.GT.j2) then
      if(i1-ip.LT.buffer) then
        VECBIVI1(ip:i1)=VECBIVI1(j1p:j2)
        BIVIINDX1(ip:i1)=BIVIINDX1(j1p:j2)
      else
        l=(i1-ip)/buffer
        do m=0,l-1,1
          b1=ip+buffer*m
          e1=b1+buffer-1
          b2=j1p+buffer*m
          e2=b2+buffer-1
          VECBIVI1(b1:e1)=VECBIVI1(b2:e2)
          BIVIINDX1(b1:e1)=BIVIINDX1(b2:e2)
        end do
          VECBIVI1(e1+1:i1)=VECBIVI1(e2+1:j2)
          BIVIINDX1(e1+1:i1)=BIVIINDX1(e2+1:j2)
      end if
    else
      if(j2-ip.LT.buffer) then
        VECBIVI1(i1-j2+ip:i1)=VECBIVI1(ip:j2)
        BIVIINDX1(i1-j2+ip:i1)=BIVIINDX1(ip:j2)
      else
        l=(j2-ip)/buffer
        e1=i1-j2+ip+buffer*l
        e2=ip+buffer*l
        VECBIVI1(e1:i1)=VECBIVI1(e2:j2)
        BIVIINDX1(e1:i1)=BIVIINDX1(e2:j2)
        do m=l-1,0,-1
          b1=e1-buffer
          b2=e2-buffer
          VECBIVI1(b1:e1-1)=VECBIVI1(b2:e2-1)
          BIVIINDX1(b1:e1-1)=BIVIINDX1(b2:e2-1)
          e1=b1
          e2=b2
        end do
      end if
      if(i-j1p.LT.buffer) then
        VECBIVI1(ip:ip+i-j1p)=VECBIVI1(j1p:i)
        BIVIINDX1(ip:ip+i-j1p)=BIVIINDX1(j1p:i)
      else
        l=(i-j1p)/buffer
        e1=ip+buffer*l
        e2=j1p+buffer*l
        VECBIVI1(e1:ip+i-j1p)=VECBIVI1(e2:i)
        BIVIINDX1(e1:ip+i-j1p)=BIVIINDX1(e2:i)
        do m=l-1,0,-1
          b1=e1-buffer
          b2=e2-buffer
          !WRITE (*,"(i10,i10,i10,i10)"),b1,e1,b2,e2
          VECBIVI1(b1:e1-1)=VECBIVI1(b2:e2-1)
          BIVIINDX1(b1:e1-1)=BIVIINDX1(b2:e2-1)
          e1=b1
          e2=b2
        end do
      end if

!    if(ip.GT.j2) then
!      if(i1-ip.LT.buffer) then
!        VECBIVI1(ip:i1)=VECBIVI1(j1p:j2)
!        BIVIINDX1(ip:i1)=BIVIINDX1(j1p:j2)
!      else
!        l=(i1-ip)/buffer
!        !WRITE (*,"(i10,i10,i10,i10,i10)"),l,ip,i1,j1p,j2
!        do m=0,l-1,1
!          b1=ip+buffer*m
!          e1=b1+buffer-1
!          b2=j1p+buffer*m
!          e2=b2+buffer-1
!          !WRITE (*,"(i10,i10,i10,i10,i10)"),m,b1,e1,b2,e2
!          VECBIVI1(b1:e1)=VECBIVI1(b2:e2)
!          BIVIINDX1(b1:e1)=BIVIINDX1(b2:e2)
!        end do
!          VECBIVI1(e1+1:i1)=VECBIVI1(e2+1:j2)
!          BIVIINDX1(e1+1:i1)=BIVIINDX1(e2+1:j2)
!      end if
!    else
!      if(j2-ip.LT.buffer) then
!        VECBIVI1(i1-j2+ip:i1)=VECBIVI1(ip:j2)
!        BIVIINDX1(i1-j2+ip:i1)=BIVIINDX1(ip:j2)
!      else
!        l=(j2-ip)/buffer
!        !WRITE (*,"(i10,i10,i10,i10,i10)"),l,ip,i1,j1p,j2
!        do m=0,l-1,1
!          b1=i1-j2+ip+buffer*m
!          e1=b1+buffer-1
!          b2=ip+buffer*m
!          e2=b2+buffer-1
!          !WRITE (*,"(i10,i10,i10,i10,i10)"),m,b1,e1,b2,e2
!          VECBIVI1(b1:e1)=VECBIVI1(b2:e2)
!          BIVIINDX1(b1:e1)=BIVIINDX1(b2:e2)
!        end do
!          VECBIVI1(e1+1:i1)=VECBIVI1(e2+1:j2)
!          BIVIINDX1(e1+1:i1)=BIVIINDX1(e2+1:j2)
!      end if
!      if(i-j1p.LT.buffer) then
!        VECBIVI1(ip:ip+i-j1p)=VECBIVI1(j1p:i)
!        BIVIINDX1(ip:ip+i-j1p)=BIVIINDX1(j1p:i)
!      else
!        l=(i-j1p)/buffer
!        !WRITE (*,"(i10,i10,i10,i10,i10)"),l,ip,i1,j1p,j2
!        do m=0,l-1,1
!          b1=ip+buffer*m
!          e1=b1+buffer-1
!          b2=j1p+buffer*m
!          e2=b2+buffer-1
!          !WRITE (*,"(i10,i10,i10,i10,i10)"),m,b1,e1,b2,e2
!          VECBIVI1(b1:e1)=VECBIVI1(b2:e2)
!          BIVIINDX1(b1:e1)=BIVIINDX1(b2:e2)
!        end do
!          VECBIVI1(e1+1:ip+i-j1p)=VECBIVI1(e2+1:i)
!          BIVIINDX1(e1+1:ip+i-j1p)=BIVIINDX1(e2+1:i)
!      end if
    end if
    end if
    i1=i
    j2=j1
    do while (BIVIINDX1(j).LT.BIVIINDX0(j1))
    !VECBIVI1(i)=VECBIVI0(j1)
    !BIVIINDX1(i)=BIVIINDX0(j1)
    j1=j1-1
    !i=i-1
    end do
    !WRITE (*,"(i10,i10,i10,i10,i10,i10,i20,i20)"),i,i1,j2,j3,j1,j,BIVIINDX1(j),BIVIINDX0(j1)
    i=i1-(j2-j1)
    if(i.NE.i1) then
    ip=i+1
    j1p=j1+1
      if(i1-ip.LT.buffer) then
        VECBIVI1(ip:i1)=VECBIVI0(j1p:j2)
        BIVIINDX1(ip:i1)=BIVIINDX0(j1p:j2)
      else
        l=(i1-ip)/buffer
        !WRITE (*,"(i10,i10,i10,i10,i10)"),l,ip,i1,j1p,j2
        do m=0,l-1,1
          b1=ip+buffer*m
          e1=b1+buffer-1
          b2=j1p+buffer*m
          e2=b2+buffer-1
          !WRITE (*,"(i10,i10,i10,i10,i10)"),m,b1,e1,b2,e2
          VECBIVI1(b1:e1)=VECBIVI0(b2:e2)
          BIVIINDX1(b1:e1)=BIVIINDX0(b2:e2)
        end do
          VECBIVI1(e1+1:i1)=VECBIVI0(e2+1:j2)
          BIVIINDX1(e1+1:i1)=BIVIINDX0(e2+1:j2)
      end if
    end if
    !20 CONTINUE
  end do
END SUBROUTINE MY_SPAR_ADDL

SUBROUTINE MY_SPAR_ADD1(VECBIVI1,BIVIINDX1,IRN,JCN,NZ1,VECBIVI0,BIVIINDX0,NZ0,NZ2,NCOL1)
  use constants
  IMPLICIT NONE
  real(kind=DPC) VECBIVI1(*),VECBIVI0(*)
  INTEGER IRN(*),JCN(*),NZ1(*),NZ0(*),NZ2(*),NCOL1(*)
  INTEGER BIVIINDX0(0:NZ0(1)),BIVIINDX1(0:NZ1(1))
  INTEGER j,j1,j1p,j2p,i,ncol,i1,j2,j3,ip,b1,e1,b2,b3,e3,e2,buffer,m,l
  buffer=100000!to avoid crash in fotran assignment
  j=NZ1(1)
  j1=NZ0(1)
  ncol=NCOL1(1)
  i=NZ2(1)
  do while(i.GT.0)!i=NZ2(1),1,-1
    i1=i
    j2=j
    j3=j1
    !WRITE (*,"(i10,i10,i10,i10,i10)"),i,i1,j2,j3,j1
    do while(BIVIINDX1(j).EQ.BIVIINDX0(j1).AND.BIVIINDX1(j).GT.-1)
          !VECBIVI1(i)=VECBIVI1(j)+VECBIVI0(j1)
          !IRN(i)=BIVIINDX1(j)/ncol+1
          !JCN(i)=MOD(BIVIINDX1(j),ncol)+1
          j1=j1-1
          j=j-1
          !i=i-1
          !GO TO 20
    end do
    !WRITE (*,"(i10,i10,i10,i10,i10,i10,i10)"),i,i1,j2,j3,j,BIVIINDX1(j),BIVIINDX0(j1)
    i=i1-(j2-j)
    if(i.NE.i1) then
    ip=i+1
    j1p=j+1
    j2p=j1+1
    if(ip.GT.j2) then
      if(i1-ip.LT.buffer) then
        VECBIVI1(ip:i1)=VECBIVI1(j1p:j2)+VECBIVI0(j2p:j3)
        IRN(ip:i1)=BIVIINDX1(j1p:j2)/ncol+1
        JCN(ip:i1)=MOD(BIVIINDX1(j1p:j2),ncol)+1
      else
        l=(i1-ip)/buffer
        do m=0,l-1,1
          b1=ip+buffer*m
          e1=b1+buffer-1
          b2=j1p+buffer*m
          e2=b2+buffer-1
          b3=j2p+buffer*m
          e3=b3+buffer-1
          VECBIVI1(b1:e1)=VECBIVI1(b2:e2)+VECBIVI0(b3:e3)
          IRN(b1:e1)=BIVIINDX1(b2:e2)/ncol+1
          JCN(b1:e1)=MOD(BIVIINDX1(b2:e2),ncol)+1
        end do
          VECBIVI1(e1+1:i1)=VECBIVI1(e2+1:j2)+VECBIVI0(e3+1:j3)
          IRN(e1+1:i1)=BIVIINDX1(e2+1:j2)/ncol+1
          JCN(e1+1:i1)=MOD(BIVIINDX1(e2+1:j2),ncol)+1
      end if
    else
!      do m=j2,ip,-1
!        VECBIVI1(i1-j2+m)=VECBIVI1(m)+VECBIVI0(j3-j2+m)
!      end do
      if(j2-ip.LT.buffer) then
        VECBIVI1(i1-j2+ip:i1)=VECBIVI1(ip:j2)+VECBIVI0(j3-j2+ip:j3)
        IRN(i1-j2+ip:i1)=BIVIINDX1(ip:j2)/ncol+1
        JCN(i1-j2+ip:i1)=MOD(BIVIINDX1(ip:j2),ncol)+1
      else
        l=(j2-ip)/buffer
        e1=i1-j2+ip+buffer*l
        e2=ip+buffer*l
        e3=j3-j2+ip+buffer*l
        VECBIVI1(e1:i1)=VECBIVI1(e2:j2)+VECBIVI0(e3:j3)
        IRN(e1:i1)=BIVIINDX1(e2:j2)/ncol+1
        JCN(e1:i1)=MOD(BIVIINDX1(e2:j2),ncol)+1
        do m=l-1,0,-1
          b1=e1-buffer
          b2=e2-buffer
          b3=e3-buffer
          VECBIVI1(b1:e1-1)=VECBIVI1(b2:e2-1)+VECBIVI0(b3:e3-1)
          IRN(b1:e1-1)=BIVIINDX1(b2:e2-1)/ncol+1
          JCN(b1:e1-1)=MOD(BIVIINDX1(b2:e2-1),ncol)+1
          e1=b1
          e2=b2
          e3=b3
        end do
      end if
      if(i-j1p.LT.buffer) then
        VECBIVI1(ip:ip+i-j1p)=VECBIVI1(j1p:i)+VECBIVI0(j2p:j3-j2+i)
        IRN(ip:ip+i-j1p)=BIVIINDX1(j1p:i)/ncol+1
        JCN(ip:ip+i-j1p)=MOD(BIVIINDX1(j1p:i),ncol)+1
      else
!      do m=i,j1p,-1
!        BIVIINDX1(i1-j2+m)=BIVIINDX1(m)
!      end do
!      do m=i,j1p,-1
!        VECBIVI1(i1-j2+m)=VECBIVI1(m)+VECBIVI0(j3-j2+m)
!      end do
        l=(i-j1p)/buffer
        e1=ip+buffer*l
        e2=j1p+buffer*l
        e3=j2p+buffer*l
        VECBIVI1(e1:ip+i-j1p)=VECBIVI1(e2:i)+VECBIVI0(e3:j3-j2+i)
        IRN(e1:ip+i-j1p)=BIVIINDX1(e2:i)/ncol+1
        JCN(e1:ip+i-j1p)=MOD(BIVIINDX1(e2:i),ncol)+1
        do m=l-1,0,-1
          b1=e1-buffer
          b2=e2-buffer
          b3=e3-buffer
          !WRITE (*,"(i10,i10,i10,i10)"),b1,e1,b2,e2
          VECBIVI1(b1:e1-1)=VECBIVI1(b2:e2-1)+VECBIVI0(b3:e3-1)
          IRN(b1:e1-1)=BIVIINDX1(b2:e2-1)/ncol+1
          JCN(b1:e1-1)=MOD(BIVIINDX1(b2:e2-1),ncol)+1
          e1=b1
          e2=b2
          e3=b3
        end do
      end if
      end if
    end if
    i1=i
    j2=j
    !WRITE (*,"(i10,i10,i10)"),i,j,j1
    do while(BIVIINDX1(j).GT.BIVIINDX0(j1))
          !VECBIVI1(i)=VECBIVI1(j)
          !IRN(i)=BIVIINDX1(j)/ncol+1
          !JCN(i)=MOD(BIVIINDX1(j),ncol)+1
          j=j-1
          !i=i-1
          !GO TO 20
    end do
    i=i1-(j2-j)
    if(i.NE.i1) then
    ip=i+1
    j1p=j+1
    !WRITE (*,"(i10,i10,i10,i10)"),i,j,j1p,i1

    if(ip.GT.j2) then
      if(i1-ip.LT.buffer) then
        VECBIVI1(ip:i1)=VECBIVI1(j1p:j2)
      else
        l=(i1-ip)/buffer
        do m=0,l-1,1
          b1=ip+buffer*m
          e1=b1+buffer-1
          b2=j1p+buffer*m
          e2=b2+buffer-1
          VECBIVI1(b1:e1)=VECBIVI1(b2:e2)
        end do
          VECBIVI1(e1+1:i1)=VECBIVI1(e2+1:j2)
      end if
    else
      if(j2-ip.LT.buffer) then
        VECBIVI1(i1-j2+ip:i1)=VECBIVI1(ip:j2)
      else
        l=(j2-ip)/buffer
        e1=i1-j2+ip+buffer*l
        e2=ip+buffer*l
        VECBIVI1(e1:i1)=VECBIVI1(e2:j2)
        do m=l-1,0,-1
          b1=e1-buffer
          b2=e2-buffer
          VECBIVI1(b1:e1-1)=VECBIVI1(b2:e2-1)
          e1=b1
          e2=b2
        end do
      end if
      if(i-j1p.LT.buffer) then
        VECBIVI1(ip:ip+i-j1p)=VECBIVI1(j1p:i)
      else
        l=(i-j1p)/buffer
        e1=ip+buffer*l
        e2=j1p+buffer*l
        VECBIVI1(e1:ip+i-j1p)=VECBIVI1(e2:i)
        do m=l-1,0,-1
          b1=e1-buffer
          b2=e2-buffer
          VECBIVI1(b1:e1-1)=VECBIVI1(b2:e2-1)
          e1=b1
          e2=b2
        end do
      end if
    end if

      if(i1-ip.LT.buffer) then
        IRN(ip:i1)=BIVIINDX1(j1p:j2)/ncol+1
        JCN(ip:i1)=MOD(BIVIINDX1(j1p:j2),ncol)+1
      else
        l=(i1-ip)/buffer
        !WRITE (*,"(i10,i10,i10,i10,i10)"),l,ip,i1,j1p,j2
        do m=0,l-1,1
          b1=ip+buffer*m
          e1=b1+buffer-1
          b2=j1p+buffer*m
          e2=b2+buffer-1
          !WRITE (*,"(i10,i10,i10,i10,i10)"),m,b1,e1,b2,e2
          IRN(b1:e1)=BIVIINDX1(b2:e2)/ncol+1
          JCN(b1:e1)=MOD(BIVIINDX1(b2:e2),ncol)+1
        end do
          IRN(e1+1:i1)=BIVIINDX1(e2+1:j2)/ncol+1
          JCN(e1+1:i1)=MOD(BIVIINDX1(e2+1:j2),ncol)+1
      end if
    end if
    !IRN(i+1:i1)=BIVIINDX1(j+1:j2)/ncol+1
    !JCN(i+1:i1)=MOD(BIVIINDX1(j+1:j2),ncol)+1
    i1=i
    j2=j1
    do while(BIVIINDX1(j).LT.BIVIINDX0(j1))
    !VECBIVI1(i)=VECBIVI0(j1)
    !IRN(i)=BIVIINDX0(j1)/ncol+1
    !JCN(i)=MOD(BIVIINDX0(j1),ncol)+1
    j1=j1-1
    !i=i-1
    end do
    i=i1-(j2-j1)
    ip=i+1
    j1p=j1+1
    if(i.NE.i1) then
      if(i1-ip.LT.buffer) then
        VECBIVI1(ip:i1)=VECBIVI0(j1p:j2)
        IRN(ip:i1)=BIVIINDX0(j1p:j2)/ncol+1
        JCN(ip:i1)=MOD(BIVIINDX0(j1p:j2),ncol)+1
      else
        l=(i1-ip)/buffer
        !WRITE (*,"(i10,i10,i10,i10,i10)"),l,ip,i1,j1p,j2
        do m=0,l-1,1
          b1=ip+buffer*m
          e1=b1+buffer-1
          b2=j1p+buffer*m
          e2=b2+buffer-1
          !WRITE (*,"(i10,i10,i10,i10,i10)"),m,b1,e1,b2,e2
          VECBIVI1(b1:e1)=VECBIVI0(b2:e2)
          IRN(b1:e1)=BIVIINDX0(b2:e2)/ncol+1
          JCN(b1:e1)=MOD(BIVIINDX0(b2:e2),ncol)+1
        end do
          VECBIVI1(e1+1:i1)=VECBIVI0(e2+1:j2)
          IRN(e1+1:i1)=BIVIINDX0(e2+1:j2)/ncol+1
          JCN(e1+1:i1)=MOD(BIVIINDX0(e2+1:j2),ncol)+1
      end if
    end if
    !20 CONTINUE
  end do
END SUBROUTINE MY_SPAR_ADD1

SUBROUTINE MY_SPAR_ADD1L(VECBIVI1,BIVIINDX1,IRN,JCN,NZ1,VECBIVI0,BIVIINDX0,NZ0,NZ2,NCOL1)
  use constants
  IMPLICIT NONE
  real(kind=DPC) VECBIVI1(*),VECBIVI0(*)
  INTEGER IRN(*),JCN(*),NZ1(*),NZ0(*),NZ2(*),NCOL1(*)
  INTEGER(8) BIVIINDX0(0:NZ0(1)),BIVIINDX1(0:NZ1(1))
  INTEGER j,j1,j1p,j2p,i,ncol,i1,j2,j3,ip,b1,e1,b2,b3,e3,e2,buffer,m,l
  buffer=100000!to avoid crash in fotran assignment
  j=NZ1(1)
  j1=NZ0(1)
  ncol=NCOL1(1)
  i=NZ2(1)
  do while(i.GT.0)!i=NZ2(1),1,-1
    i1=i
    j2=j
    j3=j1
    !WRITE (*,"(i10,i10,i10,i10,i10,i10,i10)"),i,i1,j2,j3,j1,BIVIINDX1(j),BIVIINDX0(j1)
    do while(BIVIINDX1(j).EQ.BIVIINDX0(j1).AND.BIVIINDX1(j).GT.-1)
          !VECBIVI1(i)=VECBIVI1(j)+VECBIVI0(j1)
          !IRN(i)=BIVIINDX1(j)/ncol+1
          !JCN(i)=MOD(BIVIINDX1(j),ncol)+1
          j1=j1-1
          j=j-1
          !i=i-1
          !GO TO 20
    end do
    !WRITE (*,"(A10,i10,i10,i10,i10,i10,i10,i10)"),'A',i,i1,j2,j3,j,BIVIINDX1(j),BIVIINDX0(j1)
    i=i1-(j2-j)
    if(i.NE.i1) then
    ip=i+1
    j1p=j+1
    j2p=j1+1
    if(ip.GT.j2) then
      if(i1-ip.LT.buffer) then
        !WRITE (*,"(i10,i10,i10,i10,i10,i10)"),ip,i1,j1p,j2,j2p,j3
        VECBIVI1(ip:i1)=VECBIVI1(j1p:j2)+VECBIVI0(j2p:j3)
        IRN(ip:i1)=BIVIINDX1(j1p:j2)/ncol+1
        JCN(ip:i1)=MOD(BIVIINDX1(j1p:j2),ncol)+1
      else
        l=(i1-ip)/buffer
        do m=0,l-1,1
          b1=ip+buffer*m
          e1=b1+buffer-1
          b2=j1p+buffer*m
          e2=b2+buffer-1
          b3=j2p+buffer*m
          e3=b3+buffer-1
          VECBIVI1(b1:e1)=VECBIVI1(b2:e2)+VECBIVI0(b3:e3)
          IRN(b1:e1)=BIVIINDX1(b2:e2)/ncol+1
          JCN(b1:e1)=MOD(BIVIINDX1(b2:e2),ncol)+1
        end do
          VECBIVI1(e1+1:i1)=VECBIVI1(e2+1:j2)+VECBIVI0(e3+1:j3)
          IRN(e1+1:i1)=BIVIINDX1(e2+1:j2)/ncol+1
          JCN(e1+1:i1)=MOD(BIVIINDX1(e2+1:j2),ncol)+1
      end if
    else
!      do m=j2,ip,-1
!        VECBIVI1(i1-j2+m)=VECBIVI1(m)+VECBIVI0(j3-j2+m)
!      end do
      if(j2-ip.LT.buffer) then
        VECBIVI1(i1-j2+ip:i1)=VECBIVI1(ip:j2)+VECBIVI0(j3-j2+ip:j3)
        IRN(i1-j2+ip:i1)=BIVIINDX1(ip:j2)/ncol+1
        JCN(i1-j2+ip:i1)=MOD(BIVIINDX1(ip:j2),ncol)+1
      else
        l=(j2-ip)/buffer
        e1=i1-j2+ip+buffer*l
        e2=ip+buffer*l
        e3=j3-j2+ip+buffer*l
        VECBIVI1(e1:i1)=VECBIVI1(e2:j2)+VECBIVI0(e3:j3)
        IRN(e1:i1)=BIVIINDX1(e2:j2)/ncol+1
        JCN(e1:i1)=MOD(BIVIINDX1(e2:j2),ncol)+1
        do m=l-1,0,-1
          b1=e1-buffer
          b2=e2-buffer
          b3=e3-buffer
          VECBIVI1(b1:e1-1)=VECBIVI1(b2:e2-1)+VECBIVI0(b3:e3-1)
          IRN(b1:e1-1)=BIVIINDX1(b2:e2-1)/ncol+1
          JCN(b1:e1-1)=MOD(BIVIINDX1(b2:e2-1),ncol)+1
          e1=b1
          e2=b2
          e3=b3
        end do
      end if
      if(i-j1p.LT.buffer) then
        VECBIVI1(ip:ip+i-j1p)=VECBIVI1(j1p:i)+VECBIVI0(j2p:j3-j2+i)
        IRN(ip:ip+i-j1p)=BIVIINDX1(j1p:i)/ncol+1
        JCN(ip:ip+i-j1p)=MOD(BIVIINDX1(j1p:i),ncol)+1
      else
!      do m=i,j1p,-1
!        BIVIINDX1(i1-j2+m)=BIVIINDX1(m)
!      end do
!      do m=i,j1p,-1
!        VECBIVI1(i1-j2+m)=VECBIVI1(m)+VECBIVI0(j3-j2+m)
!      end do
        l=(i-j1p)/buffer
        e1=ip+buffer*l
        e2=j1p+buffer*l
        e3=j2p+buffer*l
        VECBIVI1(e1:ip+i-j1p)=VECBIVI1(e2:i)+VECBIVI0(e3:j3-j2+i)
        IRN(e1:ip+i-j1p)=BIVIINDX1(e2:i)/ncol+1
        JCN(e1:ip+i-j1p)=MOD(BIVIINDX1(e2:i),ncol)+1
        do m=l-1,0,-1
          b1=e1-buffer
          b2=e2-buffer
          b3=e3-buffer
          !WRITE (*,"(i10,i10,i10,i10)"),b1,e1,b2,e2
          VECBIVI1(b1:e1-1)=VECBIVI1(b2:e2-1)+VECBIVI0(b3:e3-1)
          IRN(b1:e1-1)=BIVIINDX1(b2:e2-1)/ncol+1
          JCN(b1:e1-1)=MOD(BIVIINDX1(b2:e2-1),ncol)+1
          e1=b1
          e2=b2
          e3=b3
        end do
      end if
      end if
    end if
    i1=i
    j2=j
    !WRITE (*,"(i10,i10,i10)"),i,j,j1
    do while(BIVIINDX1(j).GT.BIVIINDX0(j1))
          !VECBIVI1(i)=VECBIVI1(j)
          !IRN(i)=BIVIINDX1(j)/ncol+1
          !JCN(i)=MOD(BIVIINDX1(j),ncol)+1
          j=j-1
          !i=i-1
          !GO TO 20
    end do
    i=i1-(j2-j)
    if(i.NE.i1) then
    ip=i+1
    j1p=j+1
    !WRITE (*,"(i10,i10,i10,i10)"),i,j,j1p,i1

    if(ip.GT.j2) then
      if(i1-ip.LT.buffer) then
        VECBIVI1(ip:i1)=VECBIVI1(j1p:j2)
      else
        l=(i1-ip)/buffer
        do m=0,l-1,1
          b1=ip+buffer*m
          e1=b1+buffer-1
          b2=j1p+buffer*m
          e2=b2+buffer-1
          VECBIVI1(b1:e1)=VECBIVI1(b2:e2)
        end do
          VECBIVI1(e1+1:i1)=VECBIVI1(e2+1:j2)
      end if
    else
      if(j2-ip.LT.buffer) then
        VECBIVI1(i1-j2+ip:i1)=VECBIVI1(ip:j2)
      else
        l=(j2-ip)/buffer
        e1=i1-j2+ip+buffer*l
        e2=ip+buffer*l
        VECBIVI1(e1:i1)=VECBIVI1(e2:j2)
        do m=l-1,0,-1
          b1=e1-buffer
          b2=e2-buffer
          VECBIVI1(b1:e1-1)=VECBIVI1(b2:e2-1)
          e1=b1
          e2=b2
        end do
      end if
      if(i-j1p.LT.buffer) then
        VECBIVI1(ip:ip+i-j1p)=VECBIVI1(j1p:i)
      else
        l=(i-j1p)/buffer
        e1=ip+buffer*l
        e2=j1p+buffer*l
        VECBIVI1(e1:ip+i-j1p)=VECBIVI1(e2:i)
        do m=l-1,0,-1
          b1=e1-buffer
          b2=e2-buffer
          VECBIVI1(b1:e1-1)=VECBIVI1(b2:e2-1)
          e1=b1
          e2=b2
        end do
      end if
    end if

      if(i1-ip.LT.buffer) then
        IRN(ip:i1)=BIVIINDX1(j1p:j2)/ncol+1
        JCN(ip:i1)=MOD(BIVIINDX1(j1p:j2),ncol)+1
      else
        l=(i1-ip)/buffer
        !WRITE (*,"(i10,i10,i10,i10,i10)"),l,ip,i1,j1p,j2
        do m=0,l-1,1
          b1=ip+buffer*m
          e1=b1+buffer-1
          b2=j1p+buffer*m
          e2=b2+buffer-1
          !WRITE (*,"(i10,i10,i10,i10,i10)"),m,b1,e1,b2,e2
          IRN(b1:e1)=BIVIINDX1(b2:e2)/ncol+1
          JCN(b1:e1)=MOD(BIVIINDX1(b2:e2),ncol)+1
        end do
          IRN(e1+1:i1)=BIVIINDX1(e2+1:j2)/ncol+1
          JCN(e1+1:i1)=MOD(BIVIINDX1(e2+1:j2),ncol)+1
      end if
    end if
    !IRN(i+1:i1)=BIVIINDX1(j+1:j2)/ncol+1
    !JCN(i+1:i1)=MOD(BIVIINDX1(j+1:j2),ncol)+1
    i1=i
    j2=j1
    do while(BIVIINDX1(j).LT.BIVIINDX0(j1))
    !VECBIVI1(i)=VECBIVI0(j1)
    !IRN(i)=BIVIINDX0(j1)/ncol+1
    !JCN(i)=MOD(BIVIINDX0(j1),ncol)+1
    j1=j1-1
    !i=i-1
    end do
    i=i1-(j2-j1)
    ip=i+1
    j1p=j1+1
    if(i.NE.i1) then
      if(i1-ip.LT.buffer) then
        VECBIVI1(ip:i1)=VECBIVI0(j1p:j2)
        IRN(ip:i1)=BIVIINDX0(j1p:j2)/ncol+1
        JCN(ip:i1)=MOD(BIVIINDX0(j1p:j2),ncol)+1
      else
        l=(i1-ip)/buffer
        !WRITE (*,"(i10,i10,i10,i10,i10)"),l,ip,i1,j1p,j2
        do m=0,l-1,1
          b1=ip+buffer*m
          e1=b1+buffer-1
          b2=j1p+buffer*m
          e2=b2+buffer-1
          !WRITE (*,"(i10,i10,i10,i10,i10)"),m,b1,e1,b2,e2
          VECBIVI1(b1:e1)=VECBIVI0(b2:e2)
          IRN(b1:e1)=BIVIINDX0(b2:e2)/ncol+1
          JCN(b1:e1)=MOD(BIVIINDX0(b2:e2),ncol)+1
        end do
          VECBIVI1(e1+1:i1)=VECBIVI0(e2+1:j2)
          IRN(e1+1:i1)=BIVIINDX0(e2+1:j2)/ncol+1
          JCN(e1+1:i1)=MOD(BIVIINDX0(e2+1:j2),ncol)+1
      end if
    end if
    !20 CONTINUE
  end do
END SUBROUTINE MY_SPAR_ADD1L

SUBROUTINE MY_SPAR_ADD4L(VECBIVI1,BIVIINDX1,IRN,JCN,NZ1,VECBIVI0,BIVIINDX0,NZ0,NZ2,NCOL1)
  use constants
  IMPLICIT NONE
  real(kind=DPC) VECBIVI1(*),VECBIVI0(*)
  INTEGER IRN(*),JCN(*),NCOL1(*)
  INTEGER(8) NZ0(*),NZ1(*),NZ2(*)
  INTEGER(8) BIVIINDX0(0:NZ0(1)),BIVIINDX1(0:NZ2(1))
  INTEGER(8) j,i,l,ncol
  ncol=NCOL1(1)
  j=NZ0(1)
  i=NZ1(1)
  do l=NZ2(1),1,-1
    if(BIVIINDX1(i).LT.BIVIINDX0(j)) then
      VECBIVI1(l)=VECBIVI0(j)
          IRN(l)=BIVIINDX0(j)/ncol+1
          JCN(l)=MOD(BIVIINDX0(j),ncol)+1
      j=j-1
    else if(BIVIINDX1(i).EQ.BIVIINDX0(j)) then
      VECBIVI1(l)=VECBIVI1(i)+VECBIVI0(j)
          IRN(l)=BIVIINDX0(j)/ncol+1
          JCN(l)=MOD(BIVIINDX0(j),ncol)+1
      j=j-1
      i=i-1
    else
      VECBIVI1(l)=VECBIVI1(i)
          IRN(l)=BIVIINDX1(i)/ncol+1
          JCN(l)=MOD(BIVIINDX1(i),ncol)+1
      i=i-1
    end if
    !WRITE (*,"(A10,i10,i10,i10,i20,i20)")'IJJJJ',i,j,l,BIVIINDX1(i),BIVIINDX0(j)
  end do
  !WRITE (*,"(A10,i10,i10)")'IJJJJ',i,j
END SUBROUTINE MY_SPAR_ADD4L


SUBROUTINE MY_SPAR_ADD1BKUP(VECBIVI1,BIVIINDX1,IRN,JCN,NZ1,VECBIVI0,BIVIINDX0,NZ0,NZ2,NCOL1)
  use constants
  IMPLICIT NONE
  real(kind=DPC) VECBIVI1(*),VECBIVI0(*)
  INTEGER IRN(*),JCN(*),BIVIINDX0(*),BIVIINDX1(*),NZ1(*),NZ0(*),NZ2(*),NCOL1(*)
  INTEGER j,j1,j1p,i,ncol,i1,j2,ip,b1,e1,b2,e2,buffer,m,l
  buffer=100000!to avoid crash in fotran assignment
  j=NZ1(1)
  j1=NZ0(1)
  ncol=NCOL1(1)
  i=NZ2(1)
  do while(i.GT.0)!i=NZ2(1),1,-1
    i1=i
    j2=j
    do while(BIVIINDX1(j).EQ.BIVIINDX0(j1))
          VECBIVI1(i)=VECBIVI1(j)+VECBIVI0(j1)
          !IRN(i)=BIVIINDX1(j)/ncol+1
          !JCN(i)=MOD(BIVIINDX1(j),ncol)+1
          j1=j1-1
          j=j-1
          i=i-1
          !GO TO 20
    end do
          IRN(i+1:i1)=BIVIINDX1(j+1:j2)/ncol+1
          JCN(i+1:i1)=MOD(BIVIINDX1(j+1:j2),ncol)+1
    i1=i
    j2=j
    do while(BIVIINDX1(j).GT.BIVIINDX0(j1))
          !VECBIVI1(i)=VECBIVI1(j)
          !IRN(i)=BIVIINDX1(j)/ncol+1
          !JCN(i)=MOD(BIVIINDX1(j),ncol)+1
          j=j-1
          !i=i-1
          !GO TO 20
    end do
    i=i1-(j2-j)
    if(i.NE.i1) then
    ip=i+1
    j1p=j+1
    if(ip.GT.j2) then
      if(i1-ip.LT.buffer) then
        VECBIVI1(ip:i1)=VECBIVI1(j1p:j2)
      else
        l=(i1-ip)/buffer
        !WRITE (*,"(i10,i10,i10,i10,i10)"),l,ip,i1,j1p,j2
        do m=0,l-1,1
          b1=ip+buffer*m
          e1=b1+buffer-1
          b2=j1p+buffer*m
          e2=b2+buffer-1
          !WRITE (*,"(i10,i10,i10,i10,i10)"),m,b1,e1,b2,e2
          VECBIVI1(b1:e1)=VECBIVI1(b2:e2)
        end do
          VECBIVI1(e1+1:i1)=VECBIVI1(e2+1:j2)
      end if
    else
      if(j2-ip.LT.buffer) then
        VECBIVI1(i1-j2+ip:i1)=VECBIVI1(ip:j2)
      else
        l=(j2-ip)/buffer
        !WRITE (*,"(i10,i10,i10,i10,i10)"),l,ip,i1,j1p,j2
        do m=0,l-1,1
          b1=i1-j2+ip+buffer*m
          e1=b1+buffer-1
          b2=ip+buffer*m
          e2=b2+buffer-1
          !WRITE (*,"(i10,i10,i10,i10,i10)"),m,b1,e1,b2,e2
          VECBIVI1(b1:e1)=VECBIVI1(b2:e2)
        end do
          VECBIVI1(e1+1:i1)=VECBIVI1(e2+1:j2)
      end if
      if(i-j1p.LT.buffer) then
        VECBIVI1(ip:ip+i-j1p)=VECBIVI1(j1p:i)
      else
        l=(i-j1p)/buffer
        !WRITE (*,"(i10,i10,i10,i10,i10)"),l,ip,i1,j1p,j2
        do m=0,l-1,1
          b1=ip+buffer*m
          e1=b1+buffer-1
          b2=j1p+buffer*m
          e2=b2+buffer-1
          !WRITE (*,"(i10,i10,i10,i10,i10)"),m,b1,e1,b2,e2
          VECBIVI1(b1:e1)=VECBIVI1(b2:e2)
        end do
          VECBIVI1(e1+1:ip+i-j1p)=VECBIVI1(e2+1:i)
      end if
    end if

      if(i1-ip.LT.buffer) then
        IRN(ip:i1)=BIVIINDX1(j1p:j2)/ncol+1
        JCN(ip:i1)=MOD(BIVIINDX1(j1p:j2),ncol)+1
      else
        l=(i1-ip)/buffer
        !WRITE (*,"(i10,i10,i10,i10,i10)"),l,ip,i1,j1p,j2
        do m=0,l-1,1
          b1=ip+buffer*m
          e1=b1+buffer-1
          b2=j1p+buffer*m
          e2=b2+buffer-1
          !WRITE (*,"(i10,i10,i10,i10,i10)"),m,b1,e1,b2,e2
          IRN(b1:e1)=BIVIINDX1(b2:e2)/ncol+1
          JCN(b1:e1)=MOD(BIVIINDX1(b2:e2),ncol)+1
        end do
          IRN(e1+1:i1)=BIVIINDX1(e2+1:j2)/ncol+1
          JCN(e1+1:i1)=MOD(BIVIINDX1(e2+1:j2),ncol)+1
      end if
    end if
    i1=i
    j2=j1
    do while(BIVIINDX1(j).LT.BIVIINDX0(j1))
    !VECBIVI1(i)=VECBIVI0(j1)
    !IRN(i)=BIVIINDX0(j1)/ncol+1
    !JCN(i)=MOD(BIVIINDX0(j1),ncol)+1
    j1=j1-1
    !i=i-1
    end do
    i=i1-(j2-j1)
    ip=i+1
    j1p=j1+1
    if(i.NE.i1) then
      if(i1-ip.LT.buffer) then
        VECBIVI1(ip:i1)=VECBIVI0(j1p:j2)
        IRN(ip:i1)=BIVIINDX0(j1p:j2)/ncol+1
        JCN(ip:i1)=MOD(BIVIINDX0(j1p:j2),ncol)+1
      else
        l=(i1-ip)/buffer
        !WRITE (*,"(i10,i10,i10,i10,i10)"),l,ip,i1,j1p,j2
        do m=0,l-1,1
          b1=ip+buffer*m
          e1=b1+buffer-1
          b2=j1p+buffer*m
          e2=b2+buffer-1
          !WRITE (*,"(i10,i10,i10,i10,i10)"),m,b1,e1,b2,e2
          VECBIVI1(b1:e1)=VECBIVI0(b2:e2)
          IRN(b1:e1)=BIVIINDX0(b2:e2)/ncol+1
          JCN(b1:e1)=MOD(BIVIINDX0(b2:e2),ncol)+1
        end do
          VECBIVI1(e1+1:i1)=VECBIVI0(e2+1:j2)
          IRN(e1+1:i1)=BIVIINDX0(e2+1:j2)/ncol+1
          JCN(e1+1:i1)=MOD(BIVIINDX0(e2+1:j2),ncol)+1
      end if
    end if
    !20 CONTINUE
  end do
END SUBROUTINE MY_SPAR_ADD1BKUP

SUBROUTINE MY_SPAR_ADD2(VECBIVI1,BIVIINDX1,IRN,JCN,NZ1,VECBIVI0,BIVIINDX0,NZ0,NZ2,NCOL1,VECBIVI2,IRN2,JCN2,J2,CNTL3)
  use constants
  IMPLICIT NONE
  real(kind=DPC) VECBIVI1(*),VECBIVI0(*),VECBIVI2(*),CNTL3(*),cntl
  INTEGER IRN(*),JCN(*),IRN2(*),JCN2(*),NZ1(*),NZ0(*),NZ2(*),NCOL1(*),j,j1,i,ncol,l,J2(*)
  INTEGER BIVIINDX0(0:NZ0(1)),BIVIINDX1(0:NZ1(1))
  j=NZ1(1)
  j1=NZ0(1)
  ncol=NCOL1(1)
  l=NZ2(1)
  cntl=cntl3(1)
  do i=NZ2(1),1,-1
    IF(BIVIINDX1(j).EQ.BIVIINDX0(j1)) then
          VECBIVI1(i)=VECBIVI1(j)+VECBIVI0(j1)
          IRN(i)=BIVIINDX1(j)/ncol+1
          JCN(i)=MOD(BIVIINDX1(j),ncol)+1
          j1=j1-1
          j=j-1
          GO TO 20
    END IF
    IF(BIVIINDX1(j).GT.BIVIINDX0(j1)) then
          VECBIVI1(i)=VECBIVI1(j)
          IRN(i)=BIVIINDX1(j)/ncol+1
          JCN(i)=MOD(BIVIINDX1(j),ncol)+1
          j=j-1
          GO TO 20
    END IF
    VECBIVI1(i)=VECBIVI0(j1)
    IRN(i)=BIVIINDX0(j1)/ncol+1
    JCN(i)=MOD(BIVIINDX0(j1),ncol)+1
    j1=j1-1
    20 CONTINUE
    IF(VECBIVI1(i).GT.cntl.OR.VECBIVI1(i).LT.-cntl) then
          VECBIVI2(l)=VECBIVI1(i)
          IRN2(l)=IRN(i)
          JCN2(l)=JCN(i)
          l=l-1
    END IF
  end do
  J2(1)=l-1
END SUBROUTINE MY_SPAR_ADD2

SUBROUTINE MY_SPAR_COMP(BIVIINDX1,NZ1,BIVIINDX0,NZ0,NZ2)
  INTEGER BIVIINDX0(*),BIVIINDX1(*),NZ1(*),NZ0(*),NZ2(*),j,j1,i
  if(BIVIINDX1(NZ1(1)).LT.BIVIINDX0(NZ0(1))) then
      j=1
      j1=0
      do i=1,NZ1(1),1
        do while (BIVIINDX1(i).GT.BIVIINDX0(j))
          j=j+1
        end do
        if(BIVIINDX1(i).LT.BIVIINDX0(j)) then
          j1=j1+1
        end if
      end do
      NZ2(1)=NZ0(1)+j1
  else
      j=1
      j1=0
      do i=1,NZ0(1),1
        do while (BIVIINDX0(i).GT.BIVIINDX1(j))
          j=j+1
          !WRITE (*,"(i10,i10,i10,i10)"),j,NZ1(1),i,NZ0(1)
        end do
        if(BIVIINDX0(i).LT.BIVIINDX1(j)) then
          j1=j1+1
        end if
      end do
      NZ2(1)=NZ1(1)+j1
  end if
END SUBROUTINE MY_SPAR_COMP

SUBROUTINE MY_SPAR_COMPL(BIVIINDX1,NZ1,BIVIINDX0,NZ0,NZ2)
  INTEGER(8) BIVIINDX0(*),BIVIINDX1(*),NZ0(*),NZ1(*),NZ2(*)
  INTEGER(8) j,j1,i
  !WRITE (*,"(A,i10)"),'in F ',BIVIINDX0(2)
  if(BIVIINDX1(NZ1(1)).LT.BIVIINDX0(NZ0(1))) then
      j=1
      j1=0
      do i=1,NZ1(1),1
        do while (BIVIINDX1(i).GT.BIVIINDX0(j))
          j=j+1
        end do
        if(BIVIINDX1(i).LT.BIVIINDX0(j)) then
          j1=j1+1
        end if
      end do
      NZ2(1)=NZ0(1)+j1
  else
      j=1
      j1=0
      do i=1,NZ0(1),1
        do while (BIVIINDX0(i).GT.BIVIINDX1(j))
          j=j+1
          !WRITE (*,"(i10,i10,i10,i10)"),j,NZ1(1),i,NZ0(1)
        end do
        if(BIVIINDX0(i).LT.BIVIINDX1(j)) then
          j1=j1+1
        end if
      end do
      NZ2(1)=NZ1(1)+j1
  end if
END SUBROUTINE MY_SPAR_COMPL

SUBROUTINE MY_VEC_COMZ(VECBIVI,BIVIINDX,col,row,colsize,NZ0,NZ1)
  use constants
  IMPLICIT NONE
  INTEGER, PARAMETER :: myint = selected_int_kind(4)
  real(kind=DPC) VECBIVI(*)
  INTEGER BIVIINDX(*),NZ1(*),col(*),row(*),colsize(*),i,j,l
  integer (kind=myint) NZ0(*)
  j=1
  do i=1,NZ0(1)
    if(VECBIVI(i).NE.0) then
      VECBIVI(j)=VECBIVI(i)
      l=i-1
      BIVIINDX(j)=col(MOD(l,colsize(1))+1)+row(l/colsize(1)+1)
      j=j+1
    end if
  end do
  NZ1(1)=j-1
END SUBROUTINE MY_VEC_COMZ

SUBROUTINE PREP48_ALU(INSIZE,IRN,JCN,VA)
  use constants
  IMPLICIT NONE
  integer NEFAC,JOB
  integer(4) JCN(*),IRN(*),INSIZE(*)
  integer M,N,NE,T,RANK,J1,J
  real(kind=DPC) VA(*)
  logical isopen
  !DOUBLE PRECISION LOGDET,SGNDET
  integer LA, MAXN
  real(kind=DPC), pointer :: CNTL(:),RINFO(:),W(:)
  integer, pointer :: ICNTL(:),INFO(:),IW(:),KEEP(:)
  character(len=1024) :: filename
  M=INSIZE(1)
  N=INSIZE(2)
  NE=INSIZE(3)
  NEFAC=INSIZE(10)
  RANK=INSIZE(11)
  J1=INSIZE(12)
  !JCN(1:NE)=JCN(1:NE)+1
  if(NEFAC.EQ.0) then
    LA=2*NE
  else
    LA=ceiling((NEFAC/100.0)*NE)
  endif
  MAXN=N
  IF (N.LT.M) THEN
    MAXN=M
  END IF
  allocate(CNTL(10),RINFO(10))
  allocate(ICNTL(20),INFO(20),IW(6*M+3*N))
  IF (FSORD.EQ.1) THEN
    CALL MA48ID(CNTL,ICNTL)
  else
    CALL MA48I(CNTL,ICNTL)
  endif
  T=M+5*N+4*N/ICNTL(6)+7
  INSIZE(13)=T
  allocate(KEEP(T))
  KEEP=0
  JOB=1
  IF (FSORD.EQ.1) THEN
    CALL MA48AD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  else
    CALL MA48A(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  endif
  IF (INFO(1).LT.0) THEN
  WRITE (6,'(A,I3)') 'Error STOP from MA48A/AD with INFO(1) =',INFO(1)
  STOP
  END IF
  deallocate(IW)
  allocate(W(M),IW(2*M+2*N))
  IF (FSORD.EQ.1) THEN
    CALL MA48BD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,&
                RINFO)
  else
    CALL MA48B(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,&
                RINFO)
  endif
  IF (INFO(1).NE.0) THEN
    WRITE (6,FMT='(A,I3/A)') 'STOP from MA48B/BD with INFO(1) =',&
    INFO(1),'Solution not possible'
    write(*,"(A,i5)") 'RANK',INFO(5)
    STOP
  END IF
  !write(*,"(A,i5)") 'RANK',5555
  J=INSIZE(16)
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A4,I4.4,I4.4,A4)") "_vav",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) VA(1:LA)
  CLOSE(J)
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A5,I4.4,I4.4,A4)") "_irnv",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) IRN(1:LA)
  CLOSE(J)
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A5,I4.4,I4.4,A4)") "_keep",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) KEEP(1:T)
  CLOSE(J)
  deallocate(CNTL,RINFO,W)
  deallocate(ICNTL,INFO,IW,KEEP)
END SUBROUTINE PREP48_ALU

SUBROUTINE PREP48_ALU1(INSIZE,IRN,JCN,VA,W,IW,KEEP)
  use constants
  IMPLICIT NONE
  integer NEFAC,JOB
  integer(4) JCN(*),IRN(*),INSIZE(*),IW(*),KEEP(*)
  integer M,N,NE,T,RANK,J1,J
  real(kind=DPC) VA(*),W(*)
  real(kind=DPC) LA1
  logical isopen
  !DOUBLE PRECISION LOGDET,SGNDET
  integer LA, MAXN
  real(kind=DPC), pointer :: CNTL(:),RINFO(:)
  integer, pointer :: ICNTL(:),INFO(:)
  character(len=1024) :: filename
  M=INSIZE(1)
  N=INSIZE(2)
  NE=INSIZE(3)
  NEFAC=INSIZE(10)
  RANK=INSIZE(11)
  J1=INSIZE(12)
  !JCN(1:NE)=JCN(1:NE)+1
  if(NEFAC.EQ.0) then
    LA=2*NE
  else
    LA1=(NEFAC/100.0)*NE
    LA=ceiling(LA1)
  endif
  if(LA.NE.INSIZE(17)) then
    !write(*,"(A,i15,i15)") 'Note LA!!!!',LA,INSIZE(17)
    LA=INSIZE(17)
  end if
  MAXN=N
  IF (N.LT.M) THEN
    MAXN=M
  END IF
  allocate(CNTL(10),RINFO(10))
  allocate(ICNTL(20),INFO(20))!,IW(6*M+3*N))
  IF (FSORD.EQ.1) THEN
    CALL MA48ID(CNTL,ICNTL)
  else
    CALL MA48I(CNTL,ICNTL)
  endif
  T=M+5*N+4*N/ICNTL(6)+7
  INSIZE(13)=T
  !allocate(KEEP(T))
  !KEEP=0
  !write(*,"(A,i15,f15.5)") 'LA',LA,LA1
  JOB=1
  IF (FSORD.EQ.1) THEN
    CALL MA48AD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  else
    CALL MA48A(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  endif
  IF (INFO(1).LT.0) THEN
  WRITE (6,'(A,I3)') 'Error STOP from MA48A/AD with INFO(1) =',INFO(1)
  STOP
  END IF
  !deallocate(IW)
  !allocate(W(M),IW(2*M+2*N))
  IF (FSORD.EQ.1) THEN
    CALL MA48BD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,&
                RINFO)
  else
    CALL MA48B(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,&
                RINFO)
  endif
  IF (INFO(1).NE.0) THEN
    WRITE (6,FMT='(A,I3/A)') 'STOP from MA48B/BD with INFO(1) =',&
    INFO(1),'Solution not possible'
    write(*,"(A,i5)") 'RANK',INFO(5)
    STOP
  END IF
  !write(*,"(A,i5)") 'RANK',5555
  J=INSIZE(16)
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A4,I4.4,I4.4,A4)") "_vav",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) VA(1:LA)
  CLOSE(J)
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A5,I4.4,I4.4,A4)") "_irnv",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) IRN(1:LA)
  CLOSE(J)
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A5,I4.4,I4.4,A4)") "_keep",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) KEEP(1:T)
  CLOSE(J)
  deallocate(CNTL,RINFO)!,W)
  deallocate(ICNTL,INFO)!,IW,KEEP)
END SUBROUTINE PREP48_ALU1

SUBROUTINE PREP48_MSOL(INSIZE,IRN,JCN,VA,IRNC,JCNC,VAC,IRNB,JCNB,VALUESB,VECBIVI,BIVINZROW0,BIVINZCOL0)!,IRNV,JCNV,VAV
  use constants
  IMPLICIT NONE
  integer NEFAC,JOB
  integer(8) BIVINZROW0(*),L3
  integer(4) JCN(*),IRN(*),INSIZE(*),BIVINZCOL0(*)!,JCNV(*),IRNV(*)
  integer(4)IRNB(*),JCNB(*)
  integer(4) JCNC(*),IRNC(*)
  integer M,N,NE,T,NC,MC,NEC,RANK,J1!,NV,MV,NEV
  real(kind=DPC) VA(*),VAC(*),VALUESB(*),VECBIVI(*)
  !DOUBLE PRECISION LOGDET,SGNDET!,VAV(*)real (8)
  integer LA, MAXN,NBIVI,MBIVI
  integer(4) I,J,L,L1,L2,L4,L5,M0,M1,M2,M3,M4,M5,MB,NB,NEB,J2,J3
  LOGICAL TRANS,checksol,isopen
  real(kind=DPC), pointer :: CNTL(:),RINFO(:),W(:),ERROR1(:),SOL(:),B(:)!,VOUT(:)!,RHS(:)A(:),
  integer, pointer :: ICNTL(:),INFO(:),IW(:),KEEP(:),JCNB1(:)!,IRNOUT(:),JCNOUT(:)!,JCN1(:),IRN1(:)
  character(len=1024) :: filename
  M=INSIZE(1)
  N=INSIZE(2)
  NE=INSIZE(3)
  MC=INSIZE(4)
  NC=INSIZE(5)
  NEC=INSIZE(6)
  MB=INSIZE(7)
  NB=INSIZE(8)
  NEB=INSIZE(9)
  NEFAC=INSIZE(10)
  RANK=INSIZE(11)
  J1=INSIZE(12)
  MBIVI=INSIZE(14)
  NBIVI=INSIZE(15)
  JCN(1:NE)=JCN(1:NE)+1
  allocate(JCNB1(NEB))
  JCNB1(1:NEB)=JCNB(1:NEB)+1
  do J=1,M
    L=J+NE
    IRN((IRN(L)+1):IRN(L+1))=J
  end do
  if(NEFAC.EQ.0) then
    LA=2*NE
  else
    LA=ceiling((NEFAC/100.0)*NE)
  endif
  MAXN=N
  IF (N.LT.M) THEN
    MAXN=M
  END IF
  allocate(CNTL(10),RINFO(10),ERROR1(3),SOL(M),B(M))!,RHS(M),W(5*MAXN)A(LA),
  allocate(ICNTL(20),INFO(20),IW(6*M+3*N))!JCN1(LA),IRN1(LA),
  IF (FSORD.EQ.1) THEN
    CALL MA48ID(CNTL,ICNTL)
  else
    CALL MA48I(CNTL,ICNTL)
  endif
  T=M+5*N+4*N/ICNTL(6)+7
  INSIZE(13)=T
  allocate(KEEP(T))
  KEEP=0
  JOB=1
  IF (FSORD.EQ.1) THEN
    CALL MA48AD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  else
    CALL MA48A(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  endif
  IF (INFO(1).LT.0) THEN
  WRITE (6,'(A,I3)') 'Error STOP from MA48A/AD with INFO(1) =',INFO(1)
  STOP
  END IF
  deallocate(IW)
  allocate(W(M),IW(2*M+2*N))
  IF (FSORD.EQ.1) THEN
    CALL MA48BD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,&
                RINFO)
  else
    CALL MA48B(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,&
                RINFO)
  endif
  IF (INFO(1).NE.0) THEN
    WRITE (6,FMT='(A,I3/A)') 'STOP from MA48B/BD with INFO(1) =',&
    INFO(1),'Solution not possible'
    write(*,"(A,i5)") 'RANK',INFO(5)
    STOP
  END IF
  JOB=1
  deallocate(W,IW)
  IF (JOB.EQ.1) THEN
    allocate(W(2*MAXN))
  else
    allocate(W(4*MAXN))
  END IF
  allocate(IW(MAXN))!A(LA),IRN1(LA),
  L=0
  TRANS = .FALSE.
  do J=1,MC-1!NC
      J3=IRNC(J+1)-IRNC(J)
    if(J3.GT.0) then
      B(1:M)=0
      do I=1,J3!NEC
        L2=I+L
        B(JCNC(L2)+1)=VAC(L2)
      end do
      L=L+J3
      IF (FSORD.EQ.1) THEN
        CALL MA48CD(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
                  B,SOL,ERROR1,W,IW,INFO)
      else
        CALL MA48C(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
                  B,SOL,ERROR1,W,IW,INFO)
      endif
      do I=1,MB-1
        if(IRNB(I).NE.IRNB(I+1)) then !GE.1E-10.OR.SOL(I).LE.-1E-10) then
              L3=BIVINZCOL0(J)+BIVINZROW0(I)!(JCNB(j2)+1)!+NBIVI*(BIVINZROW0(JCNB(j2)+1)-1)
              DO j2=IRNB(I)+1,IRNB(I+1)
                VECBIVI(L3)=VECBIVI(L3)-SOL(JCNB1(j2))*VALUESB(j2)
              end do
        end if
      end do
        if(IRNB(MB).NE.NEB) then !GE.1E-10.OR.SOL(I).LE.-1E-10) then
              L3=BIVINZCOL0(J)+BIVINZROW0(MB)!(JCNB(j2)+1)!+NBIVI*(BIVINZROW0(JCNB(j2)+1)-1)
              do j2=IRNB(MB)+1,NEB
                VECBIVI(L3)=VECBIVI(L3)-SOL(JCNB1(j2))*VALUESB(j2)
              end do
        end if
    end if
  end do
  J=MC
      J3=NEC-IRNC(J)
    if(J3.GT.0) then
      B(1:M)=0
      do I=1,J3!NEC
        L2=I+L
        B(JCNC(L2)+1)=VAC(L2)
      end do
      L=L+J3
      IF (FSORD.EQ.1) THEN
        CALL MA48CD(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
                  B,SOL,ERROR1,W,IW,INFO)
      else
        CALL MA48C(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
                  B,SOL,ERROR1,W,IW,INFO)
      endif
      do I=1,MB-1
        if(IRNB(I).NE.IRNB(I+1)) then !GE.1E-10.OR.SOL(I).LE.-1E-10) then
              L3=BIVINZCOL0(J)+BIVINZROW0(I)!(JCNB(j2)+1)!+NBIVI*(BIVINZROW0(JCNB(j2)+1)-1)
              DO j2=IRNB(I)+1,IRNB(I+1)
                VECBIVI(L3)=VECBIVI(L3)-SOL(JCNB1(j2))*VALUESB(j2)
              end do
        end if
      end do
        if(IRNB(MB).NE.NEB) then !GE.1E-10.OR.SOL(I).LE.-1E-10) then
              L3=BIVINZCOL0(J)+BIVINZROW0(MB)!(JCNB(j2)+1)!+NBIVI*(BIVINZROW0(JCNB(j2)+1)-1)
              do j2=IRNB(MB)+1,NEB
                VECBIVI(L3)=VECBIVI(L3)-SOL(JCNB1(j2))*VALUESB(j2)
              end do
        end if
    end if
  J=INSIZE(16)
    if(J.GT.98) then 
      J=7
    end if
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A4,I4.4,I4.4,A4)") "_vav",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) VA(1:LA)
  CLOSE(J)
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A5,I4.4,I4.4,A4)") "_irnv",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) IRN(1:LA)
  CLOSE(J)
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A5,I4.4,I4.4,A4)") "_keep",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) KEEP(1:T)
  CLOSE(J)
  deallocate(CNTL,RINFO,W,ERROR1,SOL,JCNB1,B)!,RHSA,
  deallocate(ICNTL,INFO,IW,KEEP)!IRN1,,JCNOUT
END SUBROUTINE PREP48_MSOL

SUBROUTINE PREP48M_MSOL(INSIZE,IRN,JCN,VA,IRNC,JCNC,VAC,IRNB,JCNB,VALUESB,VECBIVI,BIVINZROW0,BIVINZCOL0,JCNB1,SOL,B,W,IW,KEEP)!,IRNV,JCNV,VAV
  use constants
  IMPLICIT NONE
  integer NEFAC,JOB
  integer(8) BIVINZROW0(*),L3
  integer(4) JCN(*),IRN(*),INSIZE(*),BIVINZCOL0(*)!,JCNV(*),IRNV(*)
  integer(4)IRNB(*),JCNB(*),KEEP(*),IW(*)
  integer(4) JCNC(*),IRNC(*),JCNB1(*)
  integer M,N,NE,T,NC,MC,NEC,RANK,J1!,NV,MV,NEV
  real(kind=DPC) VA(*),VAC(*),VALUESB(*),VECBIVI(*),SOL(*),B(*),W(*)
  !DOUBLE PRECISION LOGDET,SGNDET!,VAV(*)real (8)
  integer LA, MAXN,NBIVI,MBIVI
  integer(4) I,J,L,L1,L2,L4,L5,M0,M1,M2,M3,M4,M5,MB,NB,NEB,J2,J3
  LOGICAL TRANS,checksol,isopen
  real(kind=DPC), pointer :: CNTL(:),RINFO(:),ERROR1(:)!,VOUT(:)!,RHS(:)A(:),,SOL(:),B(:),W(:)
  integer, pointer :: ICNTL(:),INFO(:)!,IRNOUT(:),JCNOUT(:)!,JCN1(:),IRN1(:),JCNB1(:),IW(:),KEEP(:)
  character(len=1024) :: filename
  M=INSIZE(1)
  N=INSIZE(2)
  NE=INSIZE(3)
  MC=INSIZE(4)
  NC=INSIZE(5)
  NEC=INSIZE(6)
  MB=INSIZE(7)
  NB=INSIZE(8)
  NEB=INSIZE(9)
  NEFAC=INSIZE(10)
  RANK=INSIZE(11)
  J1=INSIZE(12)
  MBIVI=INSIZE(14)
  NBIVI=INSIZE(15)
  JCN(1:NE)=JCN(1:NE)+1
  !allocate(JCNB1(NEB))
  JCNB1(1:NEB)=JCNB(1:NEB)+1
  do J=1,M
    L=J+NE
    IRN((IRN(L)+1):IRN(L+1))=J
  end do
  if(NEFAC.EQ.0) then
    LA=2*NE
  else
    LA=ceiling((NEFAC/100.0)*NE)
  endif
  if(LA.NE.INSIZE(17)) then
    !write(*,"(A,i15,i15)") 'Note LA!!!!',LA,INSIZE(17)
    LA=INSIZE(17)
  end if
  MAXN=N
  IF (N.LT.M) THEN
    MAXN=M
  END IF
  allocate(CNTL(10),RINFO(10),ERROR1(3))!,RHS(M),W(5*MAXN)A(LA),,SOL(M),B(M)
  allocate(ICNTL(20),INFO(20))!JCN1(LA),IRN1(LA),,IW(6*M+3*N)
  IF (FSORD.EQ.1) THEN
    CALL MA48ID(CNTL,ICNTL)
  else
    CALL MA48I(CNTL,ICNTL)
  endif
  T=M+5*N+4*N/ICNTL(6)+7
  INSIZE(13)=T
  !allocate(KEEP(T))
  !KEEP=0
  JOB=1
  IF (FSORD.EQ.1) THEN
    CALL MA48AD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  else
    CALL MA48A(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  endif
  IF (INFO(1).LT.0) THEN
  WRITE (6,'(A,I3)') 'Error STOP from MA48A/AD with INFO(1) =',INFO(1)
  STOP
  END IF
  !deallocate(IW)
  !allocate(W(M),IW(2*M+2*N))
  IF (FSORD.EQ.1) THEN
    CALL MA48BD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,&
                RINFO)
  else
    CALL MA48B(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,&
                RINFO)
  endif
  IF (INFO(1).NE.0) THEN
    WRITE (6,FMT='(A,I3/A)') 'STOP from MA48B/BD with INFO(1) =',&
    INFO(1),'Solution not possible'
    write(*,"(A,i5)") 'RANK',INFO(5)
    STOP
  END IF
  JOB=1
  !deallocate(W,IW)
  !IF (JOB.EQ.1) THEN
  !  allocate(W(2*MAXN))
  !else
  !  allocate(W(4*MAXN))
  !END IF
  !allocate(IW(MAXN))!A(LA),IRN1(LA),
  L=0
  TRANS = .FALSE.
  do J=1,MC-1!NC
      J3=IRNC(J+1)-IRNC(J)
    if(J3.GT.0) then
      B(1:M)=0
      do I=1,J3!NEC
        L2=I+L
        B(JCNC(L2)+1)=VAC(L2)
      end do
      L=L+J3
      IF (FSORD.EQ.1) THEN
        CALL MA48CD(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
                  B,SOL,ERROR1,W,IW,INFO)
      else
        CALL MA48C(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
                  B,SOL,ERROR1,W,IW,INFO)
      endif
      do I=1,MB-1
        if(IRNB(I).NE.IRNB(I+1)) then !GE.1E-10.OR.SOL(I).LE.-1E-10) then
              L3=BIVINZCOL0(J)+BIVINZROW0(I)!(JCNB(j2)+1)!+NBIVI*(BIVINZROW0(JCNB(j2)+1)-1)
              DO j2=IRNB(I)+1,IRNB(I+1)
                VECBIVI(L3)=VECBIVI(L3)-SOL(JCNB1(j2))*VALUESB(j2)
              end do
        end if
      end do
        if(IRNB(MB).NE.NEB) then !GE.1E-10.OR.SOL(I).LE.-1E-10) then
              L3=BIVINZCOL0(J)+BIVINZROW0(MB)!(JCNB(j2)+1)!+NBIVI*(BIVINZROW0(JCNB(j2)+1)-1)
              do j2=IRNB(MB)+1,NEB
                VECBIVI(L3)=VECBIVI(L3)-SOL(JCNB1(j2))*VALUESB(j2)
              end do
        end if
    end if
  end do
  J=MC
      J3=NEC-IRNC(J)
    if(J3.GT.0) then
      B(1:M)=0
      do I=1,J3!NEC
        L2=I+L
        B(JCNC(L2)+1)=VAC(L2)
      end do
      L=L+J3
      IF (FSORD.EQ.1) THEN
        CALL MA48CD(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
                  B,SOL,ERROR1,W,IW,INFO)
      else
        CALL MA48C(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
                  B,SOL,ERROR1,W,IW,INFO)
      endif
      do I=1,MB-1
        if(IRNB(I).NE.IRNB(I+1)) then !GE.1E-10.OR.SOL(I).LE.-1E-10) then
              L3=BIVINZCOL0(J)+BIVINZROW0(I)!(JCNB(j2)+1)!+NBIVI*(BIVINZROW0(JCNB(j2)+1)-1)
              DO j2=IRNB(I)+1,IRNB(I+1)
                VECBIVI(L3)=VECBIVI(L3)-SOL(JCNB1(j2))*VALUESB(j2)
              end do
        end if
      end do
        if(IRNB(MB).NE.NEB) then !GE.1E-10.OR.SOL(I).LE.-1E-10) then
              L3=BIVINZCOL0(J)+BIVINZROW0(MB)!(JCNB(j2)+1)!+NBIVI*(BIVINZROW0(JCNB(j2)+1)-1)
              do j2=IRNB(MB)+1,NEB
                VECBIVI(L3)=VECBIVI(L3)-SOL(JCNB1(j2))*VALUESB(j2)
              end do
        end if
    end if
  J=INSIZE(16)
    if(J.GT.98) then 
      J=7
    end if
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A4,I4.4,I4.4,A4)") "_vav",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) VA(1:LA)
  CLOSE(J)
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A5,I4.4,I4.4,A4)") "_irnv",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) IRN(1:LA)
  CLOSE(J)
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A5,I4.4,I4.4,A4)") "_keep",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) KEEP(1:T)
  CLOSE(J)
  deallocate(CNTL,RINFO,ERROR1)!,RHSA,,SOL,JCNB1,B,W
  deallocate(ICNTL,INFO)!IRN1,,JCNOUT,IW,KEEP
END SUBROUTINE PREP48M_MSOL

SUBROUTINE SPAR_MULMIN(SOL,NROW,NZ,IRN,JCN,VA,RES)
  use constants
  IMPLICIT NONE
  integer I,J
  integer(4) JCN(*),IRN(*),NROW(*),NZ(*)
  real(kind=DPC) VA(*),RES(*),SOL(0:NROW(1)-1)
      do I=1,NROW(1)-1
              DO J=IRN(I)+1,IRN(I+1)
                !$OMP ATOMIC
                RES(I)=RES(I)-VA(j)*SOL(JCN(J))
              end do
      end do
              DO J=IRN(NROW(1))+1,NZ(1)
                !$OMP ATOMIC
                RES(I)=RES(I)-VA(j)*SOL(JCN(J))
              end do
END SUBROUTINE SPAR_MULMIN

SUBROUTINE SPAR_MULADD(SOL,NROW,NZ,IRN,JCN,VA,RES)
  use constants
  IMPLICIT NONE
  integer I,J
  integer(4) JCN(*),IRN(*),NROW(*),NZ(*)
  real(kind=DPC) VA(*),RES(*),SOL(0:NROW(1)-1)
      do I=1,NROW(1)-1
              DO J=IRN(I)+1,IRN(I+1)
                RES(I)=RES(I)+VA(j)*SOL(JCN(J))
              end do
      end do
              DO J=IRN(NROW(1))+1,NZ(1)
                RES(NROW(1))=RES(NROW(1))+VA(j)*SOL(JCN(J))
              end do
END SUBROUTINE SPAR_MULADD

SUBROUTINE SPAR_MULNOADD(SOL,NROW,NZ,IRN,JCN,VA,RES)
  use constants
  IMPLICIT NONE
  integer I,J
  integer(4) JCN(*),IRN(*),NROW(*),NZ(*)
  real(kind=DPC) VA(*),RES(*),SOL(0:NROW(1)-1)
  RES(1:NROW(1))=0
      do I=1,NROW(1)-1
              DO J=IRN(I)+1,IRN(I+1)
                RES(I)=RES(I)+VA(j)*SOL(JCN(J))
              end do
      end do
              DO J=IRN(NROW(1))+1,NZ(1)
                RES(NROW(1))=RES(NROW(1))+VA(j)*SOL(JCN(J))
              end do
END SUBROUTINE SPAR_MULNOADD

SUBROUTINE SPAR_VBIVIADD(SOL,BVCOL,BVROW,BVSIZE,NROW,NCOL,NZ,IRN,JCN,VA,RES)
  use constants
  IMPLICIT NONE
  integer I,J
  integer(4) BVCOL(*),JCN(*),IRN(*),NROW(*),NZ(*),NCOL(*)
  integer(8) BVROW(*),BVSIZE(*),L
  real(kind=DPC) VA(*),RES(0:(BVSIZE(1)-1)),SOL(0:NCOL(1))
      do I=1,NROW(1)-1
        L=BVCOL(1)+BVROW(I)
              DO J=IRN(I)+1,IRN(I+1)
                !$OMP ATOMIC
                RES(L)=RES(L)-VA(J)*SOL(JCN(J))
              end do
      end do
        L=BVCOL(1)+BVROW(NROW(1))
              DO J=IRN(NROW(1))+1,NZ(1)
                !$OMP ATOMIC
                RES(L)=RES(L)-VA(J)*SOL(JCN(J))
              end do
END SUBROUTINE SPAR_VBIVIADD

SUBROUTINE PATIO_MAT(insizeda,IRN,JCN,VBIVI,IRN1A,JCN1A)
  use constants
  IMPLICIT NONE
  integer I,J,colcut,rowcut
  integer(4) JCN(*),IRN(*),IRN1A(*),insizeda(*),JCN1A(*)
  real(kind=DPC) VBIVI(*)
  LOGICAL iterstop
        do I=1,insizeda(3)
          iterstop=.TRUE.
          colcut=0
          rowcut=0
          if(insizeda(1).NE.insizeda(4)) then
          do J=insizeda(1),insizeda(4)+1,-1
            if(irn1a(J).LT.irn(I)) then
              rowcut=rowcut+1
            else if(irn1a(j).EQ.irn(i)) then
              JCN(i)=1
              IRN(i)=1
              VBIVI(i)=0
              iterstop=.false.
              exit
            endif
          end do
          endif
          if(iterstop.EQV..true.) then
            if(insizeda(2).NE.insizeda(4)) then
            do J=insizeda(2),insizeda(4)+1,-1
            if(jcn1a(J).LT.jcn(I)) then
              colcut=colcut+1
            else if(JCN1A(j).EQ.JCN(i)) then
              JCN(i)=1
              IRN(i)=1
              VBIVI(i)=0
              iterstop=.false.
              exit
            endif
            end do
            endif
          endif
          if(iterstop.EQV..true.) then
            IRN(i)=IRN(i)-rowcut
            JCN(i)=JCN(i)-colcut
          endif
        end do
END SUBROUTINE PATIO_MAT
