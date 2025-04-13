! COPYRIGHT (c) 2001 Council for the Central Laboratory
!                    of the Research Councils
!
! Version 2.1.0
! See ChangeLog for version history
!
MODULE HSL_mc66_double
  use HSL_fa14_double
  use HSL_mc65_double
  USE HSL_ZD11_double
  implicit none
  private
  public mc66_control, mc66, mc66_print_message, fa14_seed

  INTEGER, PARAMETER :: myreal = KIND( 1.0D+0 )
  INTEGER, PARAMETER :: myint = KIND( 1)

  integer (kind = myint), parameter, public :: &
       WARN_MATRIX_TOO_SMALL = 1

  ! IRN is not within the  range of {\tt [1,M]}. In
  ! matrix_construct using coordinate formatted input data.
  ! related entries excluded
  integer (kind = myint), parameter :: WARN_RANGE_IRN = 2

  ! JCN is not within the  range of {\tt [1,N]}. In
  ! matrix_construct using coordinate formatted input data.
  ! related entries excluded
  integer (kind = myint), parameter :: WARN_RANGE_JCN = 3

  ! if duplicate entries were
  !  found when cleaning the matrix. These entries are ignored.
  integer (kind = myint), parameter :: WARN_DUP_ENTRY = 4

  ! if empty column(s) was found
  integer (kind = myint), parameter :: WARN_EMPTY_COL = 5

  integer (kind = myint), parameter :: ERR_MEMORY_ALLOC = -1
  integer (kind = myint), parameter :: ERR_MEMORY_DEALLOC = -2
  integer (kind = myint), parameter :: ERR_NEG_NBLOCKS = -3
  integer (kind = myint), parameter :: ERR_M_NONPOSITIVE = -4 ! m<= 0
  integer (kind = myint), parameter :: ERR_N_NONPOSITIVE = -5 ! n <= 0
  integer (kind = myint), parameter :: ERR_NZ_NEGATIVE = -6 ! nz < 0

  type mc66_control
     ! lp: error message channel. negative if suppressed
     ! wp: warning message channel. negative if suppressed
     ! mp: diagonostic  message channel. negative if suppressed
     integer (kind = myint) :: lp=6,wp=6,mp=6

     !
     ! print level: < 0 FOR QUIET WORKING, = 0 FOR MODERATE
     !              > 0 FOR DETAILED PRINTING.
     integer (kind = myint) :: print_level = -1

     ! kl_AGGRESSIVE: WHETHER AGGRESSIVE (COSTLY) OR MODERATE kl
     ! REFINEMENT.
     !          <= 0 IF MODERATE, > 0 IF AGGRESIVE. dEFAULT -1
     ! IMBALANCE: AMOUNT OF IMBALANCE ALLOWED PER BISECTION.
     ! 0.01 <=> 1% DEFAULT 0.01
     real (kind = myreal) :: kl_AGGRESSIVE = -1
     real (kind = myreal) :: MAX_IMBALANCE = 0.01

     ! number of levels in the multilevel
     integer (kind = myint) :: mglevel = huge(0)

     ! grid_rdc_fac: minimum grid reduction factor which must
     ! be achieved during coarsening
     real (kind = myreal) :: grid_rdc_fac = 0.75_myreal

     ! smallest problem size to stop coarsening
     integer (kind = myint) :: coarsest_size = 100

     ! coarsening scheme: 1 for heavy edge coarsening
     ! 2 for a more aggressive coarsening scheme (faster, but slighly
     ! lower quality)
     integer (kind = myint) :: coarsen_scheme = 1

     ! number of random KL initial partition to be carried out on the
     ! coarsest matrix
     integer (kind = myint) :: NUM_COARSEST_KL = 4
  end type mc66_control
  type block
     type (ZD11_type), pointer :: mat
     ! super-column weights
     integer (kind = myint), dimension (:), pointer :: col_wgt
     ! how many blocks this particular submatrix will be divided into
     integer (kind = myint) :: npart
  end type block

  ! used in breadth first search: a list of blocks
  type block_list
     ! total number of blocks to achieve
     integer (kind = myint) :: max_blocks
     ! number of matrix in the list so far
     integer (kind = myint) :: num_blocks
     ! list of matrices
     type (block), pointer :: blocks(:)
     ! pointer to tell which row each of the matrix in the list
     ! start and end
     integer (kind = myint), pointer :: row_sta(:)
  end type block_list

   type multigrid
     ! m: size of this level (number of rows)
     integer (kind = myint) :: m
     ! this level of matrix
     type (ZD11_type), pointer :: graph
     ! the graph
     type (ZD11_type), pointer :: matrix
     ! where each row of this level of matrix will go
     integer (kind = myint), pointer, dimension (:) :: where
     ! row weight: how many rows does this row of the coarse grid
     ! matrix stands for
     integer (kind = myint), pointer, dimension (:) :: row_wgt
     ! the weight of each column,
     ! which is the number of original columns this column stands for
     ! (columns are merged if they have exactly the same pattern,
     ! deleted if they have only one element)
     integer (kind = myint), pointer, dimension (:) :: col_wgt
     ! the level
     integer (kind = myint) :: level
     ! pointer to the fine and coarse grid (matrix)
     type (multigrid), pointer :: coarse, fine
     ! the cut given by this partitioning 'where'
     integer (kind = myint), pointer :: mincut
     ! the prolongation operator
     type (ZD11_type), pointer :: p
  end type multigrid

  type extendable_row
     integer (kind = myint) :: num_column
     integer (kind = myint), pointer :: column(:)
     integer (kind = myint), pointer :: val(:)
  end type extendable_row

  interface expand
     module procedure iexpand1
  end interface

  type kl_type
     type (ZD11_type), pointer :: KL_matrix, KL_graph
     ! two foot_prints of the two subdomain and their depth
     integer (kind = myint), dimension (:,:), pointer :: foot_print
     ! the column weight (>1 for super-columns)
     integer (kind = myint), dimension (:), pointer :: colwgt
  end type kl_type

  ! node for double linked list
  type list_node_type
     integer (kind = myint) :: id
     type (list_node_type), pointer :: next,prev
  end type List_node_type

  ! this is the first node in the bucket, any new node will be added
  ! before this node and become the first node
  ! thus this is a last-n-first-out scheme
  type list_node_pointer
     type (list_node_type), pointer :: current_node
  end type list_node_pointer

  type queue_type
     integer (kind = myint) :: low_bound,up_bound
     ! maximum number of nodes allowed (as dimension of mynode)
     ! maximum gain so far in the buckets
     ! (give the highest non-empty bucket)
     ! total number of nodes in the buckets
     integer (kind = myint) :: maxnodes,maxgain,num_nodes
     ! the starting node of each bucket (point to NULL if empty)
     type (list_node_pointer), pointer, dimension(:) :: buckets
     ! the location of each node in the buckets
     type (list_node_type), pointer, dimension (:) :: mynode
  end type queue_type

  interface mc66
     module procedure monet
  end interface
  interface mc66_print_message
     module procedure monet_print_message
  end interface

contains

  subroutine monet_print_message(info,lp,wp,context)
    ! printing error message
    ! info: is an integer scaler of INTENT (IN).
    !       It is the information flag
    !       whose corresponding error message is to be printed.
    ! lp: is an OPTIONAL integer scaler of INTENT (IN).
    !         It is the unit number the user wish to print the
    !         error message to. If this number
    !         is negative, printing is supressed. If not supplied,
    !         unit 6 is the default unit to print the error message.
    ! wp: is an OPTIONAL integer scaler of INTENT (IN).
    !         It is the unit number the user wish to print the
    !         warning message to. If this number
    !         is negative, printing is supressed. If not supplied,
    !         unit 6 is the default unit to print the error message.
    ! context: is an OPTIONAL assumed size CHARACTER array of
    !          INTENT (IN). It describes the context under
    !          which the error occured.
    integer (kind = myint), intent (in) :: info
    integer (kind = myint), intent (in), optional :: lp,wp
    character (len = *), optional, intent (in) :: context
    integer (kind = myint) :: unit,length
!    integer (kind = myint) :: warning_tag,rem,div

    unit = 6

    if (info > 0) then
       if (present(wp)) then
          unit = wp
       end if
       if (unit <= 0) return
       write(unit,advance = "yes", fmt = "('  ')")
       write(unit,advance = "yes", fmt = "(' WARNING: ')")
    else if (info < 0) then
       if (present(lp)) then
          unit = lp
       end if
       if (unit <= 0) return
       write(unit,advance = "yes", fmt = "('  ')")
       write(unit,advance = "yes", fmt = "(' ERROR: ')")
    end if

    if (present(context)) then
       length = len_trim(context)
       write(unit,advance = "yes", fmt = "(a,' : ')") context(1:length)
    end if

    ! decompose warning message
!!$    if (info > 0) then
!!$       warning_tag = 1
!!$       div = info
!!$       do while (div > 0)
!!$          rem = mod(div,2)
!!$          div = div/2
!!$          if (rem > 0) call monet_print_single_message(unit,&
!!$               warning_tag)
!!$          warning_tag = warning_tag*2
!!$       end do
!!$    else
!!$       call monet_print_single_message(unit,info)
!!$    end if
       call monet_print_single_message(unit,info)

  end subroutine monet_print_message


  subroutine monet_print_single_message(unit,info)
    integer (kind = myint), intent (in) :: info,unit
    select case (info)
    case (0)
       write(unit,"(a)") "successful completion"
    case (WARN_MATRIX_TOO_SMALL)
       write(unit,"(A)") "matrix size is relatively small"
       write(unit,"(A)") "compared with number of blocks"
    case(WARN_RANGE_IRN)
       write(unit,"(A)") "IRN is not within the &
            &range of [1,M] in MATRIX_CONSTRUCT"
       write(unit,"(A)") "using coordinate formatted input data. &
            &Such entries are excluded"
    case(WARN_RANGE_JCN)
       write(unit,"(A)") "JCN is not within the &
            &range of [1,N] in matrix_construct"
       write(unit,"(A)") " using coordinate formatted &
            &input data. Such entries are excluded"
    case (WARN_DUP_ENTRY)
       write(unit,"(A)") "duplicate entries were found&
            & and ignored"
    case (WARN_EMPTY_COL)
       write(unit,"(A)") "empty column(s) was found"
       write(unit,"(A)") "the matrix is singular, the empty column(s) "
       write(unit,"(A)") "was placed at the last diagonal block"
    case (ERR_MEMORY_ALLOC)
       write(unit,"(A)") "memory allocation failed"
    case (ERR_MEMORY_DEALLOC)
       write(unit,"(A)") "memory deallocation failed"
    case (ERR_NEG_NBLOCKS)
       write(unit,"(A)") "number of blocks nonpositive or greater than&
            & min(M,N)"
    case(ERR_M_NONPOSITIVE)
       write(unit,"(A)") "M <= 0"
    case(ERR_N_NONPOSITIVE)
       write(unit,"(A)") "N <= 0"
    case(ERR_NZ_NEGATIVE)
       write(unit,"(A)") "NZ < 0"
    case default
       write(unit,"(a,i10,a)") "you have supplied an info flag of ",&
            info," this is not a recognized info flag"
    end select

  end subroutine monet_print_single_message


  subroutine block_list_init(blklist,maxblks,info)
    ! initialise a list of matrix pointers

    ! list of matrices
    type (block_list), intent (out) :: blklist

    ! maxiumum number of blocks in the matrix list
    integer (kind = myint), intent (in) :: maxblks

    ! i: loop index
    ! ierr: error flag for memory allocation
    integer (kind = myint) :: i,ierr

    ! info: info tag
    integer (kind = myint), intent (out) :: info

    info = 0
    ! assign enough space for the matrix list
    allocate(blklist%blocks(maxblks),stat = ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_ALLOC
       return
    end if

    blklist%max_blocks = maxblks

    ! nullify all the matrix pointer associated with each block
    do i = 1, maxblks
       nullify(blklist%blocks(i)%mat,blklist%blocks(i)%col_wgt)
    end do

    blklist%num_blocks = 1
    ! column starting and stopping index, the border
    ! is counted as the last block
    allocate(blklist%row_sta(maxblks+1),stat = ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_ALLOC
       return
    end if
    blklist%row_sta(1) = 1

  end subroutine block_list_init


  subroutine block_list_deallocate(blklist,info)
    ! list of matrices
    type (block_list), intent (inout) :: blklist

    ! loop index etc
    integer (kind = myint) :: i, ierr

    ! info: integer scaler of intent (out). Information tag
    integer (kind = myint), intent (out) :: info

    info = 0
    ! assign space for matrix type, and let matrix pointers
    ! point to them
    do i = 1, blklist%max_blocks
       if (associated(blklist%blocks(i)%mat)) then
          call mc65_matrix_destruct(blklist%blocks(i)%mat,ierr)
          if (ierr /= 0) then
             info = ERR_MEMORY_DEALLOC
             return
          end if
          deallocate(blklist%blocks(i)%col_wgt,&
               blklist%blocks(i)%mat,stat = ierr)
          if (ierr /= 0) then
             info = ERR_MEMORY_DEALLOC
             return
          end if
       end if
    end do
    blklist%num_blocks = 0
    deallocate(blklist%blocks,blklist%row_sta,stat = ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_DEALLOC
       return
    end if
  end subroutine block_list_deallocate


  function num_orig_row(grid)
    ! find the size of finest-level matrix

    integer (kind = myint) :: num_orig_row

    ! the coarse grid
    type (multigrid), target :: grid
    ! grid iterator
    type (multigrid), pointer :: grid_temp

    num_orig_row = grid%m
    grid_temp => grid%fine
    do while (associated(grid_temp))
       num_orig_row = grid_temp%m
       grid_temp => grid_temp%fine
    end do
  end function num_orig_row


  subroutine random_permute(order,seed)
    ! generate a random permutation of order
    integer (kind = myint), intent (inout) :: order(:)
    ! random number seed
    type (fa14_seed), intent (inout) :: seed
    integer (kind = myint):: n,left,ii,tmp
    real (kind = myreal) :: x

    n = size(order)
    left = n
    do left = n, 2, -1
       call fa14_random_real(seed,.true.,x)
       ii = min(int(x*left)+1,left)
       tmp = order(ii)
       order(ii) = order(left)
       order(left) = tmp
    end do
  end subroutine random_permute


  subroutine coarsen(grid,controller,info)
! ------------------------------------
! coarsen the grid and set up the
! coarse grid equation, the prolongator
! and restrictor
! ------------------------------------

    ! grid: the grid to be coarsened.
    type (multigrid), intent (inout), target :: grid

    ! controller: the controller containing control parameters
    type (mc66_control), intent (in) :: controller

    ! info: the info tag
    integer (kind = myint), intent (out) :: info

    type (multigrid), pointer :: cgrid

    integer (kind = myint) :: ierr,estream,wstream

    estream = controller%lp
    wstream = controller%wp

    info = 0
    allocate(cgrid,stat = ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_ALLOC
       return
    end if

    cgrid%fine => grid

    grid%coarse => cgrid

    if (controller%coarsen_scheme == 2) then
       ! find the prolongator using more aggressive coarsening
       call prolng_heavy_edges(grid,ierr)
       if (ierr /= 0) then
          info = ierr
          call monet_print_message(info,estream,wstream,&
               "prolng_heavy_edges")
          if (info < 0) return
       end if
    else
       ! find the prolongator using heavy edge coarsening
       call prolng_heavy_edge(grid,ierr)
       if (ierr /= 0) then
          info = ierr
          call monet_print_message(info,estream,wstream,&
               "prolng_heavy_edge")
          if (info < 0) return
       end if
    end if

    cgrid%level = grid%level + 1

  end subroutine coarsen


  subroutine prolng_heavy_edge(grid,info)
! ------------------------------------
!    calculate the prolongator:
!    match the vertices of the heaviest
!    edges
! ------------------------------------

    ! input fine grid
    type (multigrid), intent (inout), target :: grid

    ! info: information
    integer (kind = myint), intent (out) :: info

    ! coarse grid based on the fine grid
    type (multigrid), pointer :: cgrid

    ! the fine grid row connectivity graph
    type (ZD11_type), pointer :: graph
    ! the coarse grid prolngator
    type (ZD11_type), pointer :: p

    ! the number of fine and coarse grid vertices
    integer (kind = myint) :: nvtx, cnvtx

    ! working variables
    integer (kind = myint) :: v,u,ierr,j,i

    integer (kind = myint), pointer, dimension (:) :: the_row
    real (kind = myreal), pointer, dimension (:) :: the_row_val

    ! whether a vertex is matched already
    integer (kind = myint), parameter :: unmatched=-1

    ! matching status of each vertex
    integer (kind = myint), allocatable, dimension (:) :: match

    ! maximum weight and index of edges connected to
    ! the current vertex
    real (kind = myreal) :: maxwgt
    integer (kind = myint) :: maxind

    ! order of which vertex is visited for matching
    integer (kind = myint), allocatable :: order(:),nrow(:)
    integer (kind = myint) :: nz
    integer (kind = myint), dimension (:), pointer :: ia

    info = 0

    ! Allocate local arrays
    allocate(order(grid%m),nrow(grid%m),stat = ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_ALLOC
       return
    end if


    ! allocate the prolongation matrix pointers
    cgrid => grid%coarse
    graph => grid%graph
    allocate(cgrid%p,stat = ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_ALLOC
       return
    end if

    ! allocate the graph and matrix pointer and the mincut pointer
    ! so that everything is defined
    allocate(cgrid%graph, cgrid%matrix, cgrid%mincut,stat = ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_ALLOC
       return
    end if


    p => cgrid%p
    nvtx = graph%m

    ! prolongator start here ================================

    ! initialise the matching status
    allocate(match(nvtx),stat = ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_ALLOC
       return
    end if

    match(:) = unmatched

    ! randomly permute the vertex order
    do i = 1, nvtx
       order(i) = i
    end do
!    call random_permute(order)

    ! loop over each vertex and match along the heaviest edge
    cnvtx = 0
    nz = 0
    nrow = 0
    do i = 1, nvtx
       v = order(i)
       ! already matched, next vertex please
       if (match(v) /= unmatched) cycle
!       the_row => mc65_matrix_getrow(graph,v)
!       the_row_val => mc65_matrix_getrowval(graph,v)
       call mc65_matrix_getrow(graph,v,the_row)
       call mc65_matrix_getrowval(graph,v,the_row_val)
       maxwgt=-huge(0.0_myreal)
       maxind=v ! in the case no match is found then match itself
       do j = 1, size(the_row)
          u = the_row(j)
          ! heavy edge matching
          if (match(u)==unmatched.and.maxwgt < &
               abs(the_row_val(j))) then
             maxwgt=abs(the_row_val(j))
             maxind=u
          end if
       end do
       ! the neighbor with heaviest weight
       match(v) = maxind
       match(maxind) = v
       cnvtx = cnvtx + 1
       ! construct the prolongation matrix:
       ! find vertex v and maxind is linked
       ! with the coarse grid vertex cnvtx
       nz = nz + 1
       nrow(v) = nrow(v) + 1
       if (maxind /= v) then
          nz = nz + 1
          nrow(maxind) = nrow(maxind) + 1
       end if
    end do

    call mc65_matrix_construct(p,nvtx,nz,ierr,cnvtx,type = "pattern")
    if (ierr /= 0) then
       info = ERR_MEMORY_ALLOC
       return
    end if

    ia => p%ptr
    ia(1) = 1
    do i = 1, nvtx
       ia(i+1) = ia(i)+nrow(i)
    end do
    ! the row counter: nothing in matrix P filled yet
    nrow = 0

    match = unmatched
    ! loop over each vertex and match along the heaviest edge
    cnvtx = 0
    do i = 1, nvtx
       v = order(i)
       ! already matched, next vertex please
       if (match(v) /= unmatched) cycle
!      the_row => mc65_matrix_getrow(graph,v)
!      the_row_val => mc65_matrix_getrowval(graph,v)
       call mc65_matrix_getrow(graph,v,the_row)
       call mc65_matrix_getrowval(graph,v,the_row_val)
       maxwgt=-huge(0.0_myreal)
       maxind=v ! in the case no match is found then match itself
       do j = 1, size(the_row)
          u = the_row(j)
          ! heavy edge matching
          if (match(u)==unmatched.and.maxwgt < &
               abs(the_row_val(j))) then
             maxwgt=abs(the_row_val(j))
             maxind=u
          end if
       end do
       ! the neighbor with heaviest weight
       match(v) = maxind
       match(maxind) = v
       cnvtx = cnvtx + 1
       ! construct the prolongation matrix:
       ! find vertex v and maxind is linked
       ! with the coarse grid vertex cnvtx
       p%col(p%ptr(v)+nrow(v)) = cnvtx
       nrow(v) = nrow(v)+1
       if (maxind /= v) then
          p%col(p%ptr(maxind)+nrow(maxind)) = cnvtx
          nrow(maxind) = nrow(maxind)+1
       end if

    end do

    ! deallocate the match pointer for vertices
    deallocate(match,STAT = IERR)
    if (ierr /= 0) then
       info = ERR_MEMORY_DEALLOC
       return
    end if

    ! prolongator end ================================

    ! size of coarse grid
    cgrid%m = cnvtx

  end subroutine prolng_heavy_edge


  subroutine prolng_heavy_edges(grid,info)
    ! ------------------------------------
    !    calculate the prolongator:
    !    match the vertices of the heaviest
    !    edges. However unlike heavy edge,
    !    in this scheme a vertex can merge with an
    !    already merged (matched) vertex
    ! ------------------------------------

     ! input fine grid
    type (multigrid), intent (inout), target :: grid

    ! info: information tag
    integer (kind = myint), intent (out) :: info

    ! coarse grid based on the fine grid
    type (multigrid), pointer :: cgrid

    ! the fine grid row connectivity graph
    type (ZD11_type), pointer :: graph
    ! the coarse grid prolngator
    type (ZD11_type), pointer :: p

    ! the number of fine and coarse grid vertices
    integer (kind = myint) :: nvtx, cnvtx

    ! working variables
    integer (kind = myint) :: v,u,ierr,j,i

    integer (kind = myint), pointer, dimension (:) :: the_row
    real (kind = myreal), pointer, dimension (:) :: the_row_val

    ! whether a vertex is matched already
    integer (kind = myint), parameter :: unmatched=-1

    ! matching status of each vertex
    integer (kind = myint), allocatable, dimension (:) :: match

    ! maximum weight and index of edges connected
    ! to the current vertex
    real (kind = myreal) :: maxwgt
    integer (kind = myint) :: maxind

    ! order of which vertex is visited for matching
    integer (kind = myint), allocatable :: order(:),nrow(:)
    integer (kind = myint) :: nz
    integer (kind = myint), dimension (:), pointer :: ia

    ! crow_wgt: coarse vertices weight
    ! row_wgt: fine vertices weight
    ! cluster weight: weight of a potential new cluster
    integer (kind = myint), dimension (:), pointer :: row_wgt
    integer (kind = myint), allocatable, dimension (:) :: crow_wgt
    real (kind = myreal) :: cluster_wgt

    info = 0

    ! Allocate local arrays
    allocate(order(grid%m),nrow(grid%m),crow_wgt(grid%m),stat = ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_ALLOC
       return
    end if

    ! allocate the prolongation matrix pointers
    cgrid => grid%coarse

    graph => grid%graph

    ! allocate the graph and matrix pointer and the mincut pointer
    ! so that everything is defined

    allocate(cgrid%p,cgrid%graph, cgrid%matrix, &
         cgrid%mincut,stat = ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_ALLOC
    end if

    p => cgrid%p
    nvtx = graph%m

    ! prolongator start here ================================

    ! initialise the matching status
    allocate(match(nvtx),stat = ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_ALLOC
    end if

    match(:) = unmatched

    ! randomly permute the vertex order
    do i = 1, nvtx
       order(i) = i
    end do
!    call random_permute(order)
    ! loop over each vertex and match along the heaviest edge
    cnvtx = 0
    nz = 0
    nrow(:) = 0
    row_wgt => grid%row_wgt
    do i = 1, nvtx
       v = order(i)
       ! already matched, next vertex please
       if (match(v) /= unmatched) cycle
!       the_row => mc65_matrix_getrow(graph,v)
!       the_row_val => mc65_matrix_getrowval(graph,v)
       call mc65_matrix_getrow(graph,v,the_row)
       call mc65_matrix_getrowval(graph,v,the_row_val)
       maxwgt=-huge(0.0_myreal)
       maxind=v ! in the case no match is found then match itself
       ! heavy edge matching, but penalise forming heavy
       ! clusters
       do j = 1, size(the_row)
          u = the_row(j)
          ! heavy edge matching
          if (match(u) /= unmatched) then
             cluster_wgt = row_wgt(v) + crow_wgt(match(u))
          else
             cluster_wgt = row_wgt(v) + row_wgt(u)
          end if
 !         if (maxwgt < abs(the_row_val(j))/(cluster_wgt).and.&
         if (maxwgt < abs(the_row_val(j)).and.&
             cluster_wgt < 0.15*num_orig_row(grid)) then
             maxwgt=abs(the_row_val(j))
             maxind=u
          end if
       end do
       ! testing to control the excessive coarsening
       !       if (cnvtx < i*0.5) maxind = v
       ! a new coarse node is added only
       ! if the neighbor is not already matched
       if (match(maxind) == unmatched) then
          ! match the two end node to the new coarse node
          cnvtx = cnvtx + 1
          match(v) = cnvtx
          match(maxind) = cnvtx
          ! a new entry linking v with cnvtx
          nz = nz + 1
          nrow(v) = nrow(v) + 1
          ! another new entry if v is not matched to itself
          if (maxind /= v) then
             nz = nz + 1
             nrow(maxind) = nrow(maxind) + 1
             crow_wgt(cnvtx) = row_wgt(maxind) + row_wgt(v)
          else
             crow_wgt(cnvtx) = row_wgt(v)
          end if
       else
          ! match v to the old coarse node that its
          ! neighbor is matched to
          match(v)= match(maxind)
          ! a new entry linking v with the old coarse node
          nz = nz + 1
          nrow(v) = nrow(v) + 1
          crow_wgt(match(maxind)) = crow_wgt(match(maxind)) + &
               row_wgt(v)
       end if
       ! construct the prolongation matrix:
       ! find vertex v and maxind is linked
       ! with the coarse grid vertex cnvtx
    end do

    call mc65_matrix_construct(p,nvtx,nz,ierr,cnvtx,type = "pattern")
    if (ierr /= 0) then
       info = ERR_MEMORY_ALLOC
    end if

    ia => p%ptr
    ia(1) = 1
    do i = 1, nvtx
       ia(i+1) = ia(i)+nrow(i)
    end do
    ! the row counter: nothing in matrix P filled yet
    nrow = 0

    match = unmatched
    ! loop over each vertex and match along the heaviest edge
    cnvtx = 0
    do i = 1, nvtx
       v = order(i)
       ! already matched, next vertex please
       if (match(v) /= unmatched) cycle
!       the_row => mc65_matrix_getrow(graph,v)
!       the_row_val => mc65_matrix_getrowval(graph,v)
       call mc65_matrix_getrow(graph,v,the_row)
       call mc65_matrix_getrowval(graph,v,the_row_val)
       maxwgt=-huge(0.0_myreal)
       maxind=v ! in the case no match is found then match itself
       do j = 1, size(the_row)
          u = the_row(j)
          ! heavy edge matching
          if (match(u) /= unmatched) then
             cluster_wgt = row_wgt(v) + crow_wgt(match(u))
          else
             cluster_wgt = row_wgt(v) + row_wgt(u)
          end if
!          if (maxwgt < abs(the_row_val(j))/(cluster_wgt).and.&
          if (maxwgt < abs(the_row_val(j)).and.&
             cluster_wgt < 0.15*num_orig_row(grid)) then
             maxwgt=abs(the_row_val(j))
             maxind=u
          end if
       end do
       ! testing to control the excessive coarsening
!       if (cnvtx < i*0.5) maxind = v
       ! a new coarse node is added only
       ! if the neighbor is not already matched
       if (match(maxind) == unmatched) then
          ! match the two end node to the new coarse node
          cnvtx = cnvtx + 1
          match(v) = cnvtx
          match(maxind) = cnvtx
          ! construct the prolongation matrix:
          ! fine vertex v and maxind is linked
          ! with the coarse grid vertex cnvtx
          p%col(p%ptr(v)+nrow(v)) = cnvtx
          nrow(v) = nrow(v)+1
          if (maxind /= v) then
             p%col(p%ptr(maxind)+nrow(maxind)) = cnvtx
             nrow(maxind) = nrow(maxind)+1
             crow_wgt(cnvtx) = row_wgt(maxind) + row_wgt(v)
          else
             crow_wgt(cnvtx) = row_wgt(v)
          end if
       else
          ! match v to the old coarse node that its neighbor
          ! is matched to
          match(v)= match(maxind)
          ! a new entry linking v with the old coarse node
          p%col(p%ptr(v)+nrow(v)) = match(maxind)
          nrow(v) = nrow(v)+1
          crow_wgt(match(maxind)) = crow_wgt(match(maxind)) + &
               row_wgt(v)
       end if
    end do

    ! deallocate the match pointer for vertices
    deallocate(match,stat = ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_DEALLOC
    end if

! prolongator end ================================

! size of coarse grid
    cgrid%m = cnvtx

  end subroutine prolng_heavy_edges


  subroutine extrow_construct(row,n,info)
    type ( extendable_row), intent (out) :: row
    integer (kind = myint), intent (in) :: n
    integer (kind = myint), intent (out) :: info
    integer (kind = myint) :: ierr
    info = 0
    allocate(row%column(n),row%val(n),stat = ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_ALLOC
       return
    end if
    row%num_column = 0 ! no content yet
  end subroutine extrow_construct


  subroutine extrow_destruct(row,info)
    type ( extendable_row), intent (inout) :: row
    integer (kind = myint), intent (out) :: info
    integer (kind = myint) :: ierr
    info = 0
    deallocate(row%column,row%val,stat = ierr)
    if (ierr /= 0) info = ERR_MEMORY_DEALLOC
  end subroutine extrow_destruct


  subroutine extrow(the_row,column_index, column_value, info)
    ! extend an csr_matrix_row by upto 10% in memory
    ! danger --
    ! only used with a stand alone csr_matrix_row, do not use
    ! with a row in the csr_matrix since that
    ! row is pointed to a small part of a large
    ! contiguous space and cannot be "expanded"

    ! the_row: the row to be extended
    type (extendable_row), intent (inout) :: the_row
    ! column_index: integer content of the row
    integer (kind = myint), intent (in) :: column_index
    ! column_value: real (double) content of the row
    integer (kind = myint), intent (in) :: column_value
    ! info: info flag. info /= 0 if memory allocation failed
    integer (kind = myint), intent (out) :: info
    integer (kind = myint) :: sp

    info = 0
    the_row%num_column = the_row%num_column + 1
    if (the_row%num_column > size(the_row%column)) then
       sp = the_row%num_column + max(10,the_row%num_column/10)
       call expand(the_row%column,sp,info)
       if (info /= 0) then
          return
       end if
       call expand(the_row%val,sp,info)
       if (info /= 0) then
          return
       end if
    end if
    the_row%column(the_row%num_column) = column_index
    the_row%val(the_row%num_column) = column_value

  end subroutine extrow


  ! =========== for integer (kind = myint) arrays ================
  subroutine iexpand1(p1,upbound_new1,info)
    integer (kind = myint), pointer:: p2(:),p1(:)
    integer (kind = myint) upbound_new1
    integer (kind = myint) ub(1),lb(1),upb(1)
    integer (kind = myint) info
    integer (kind = myint) i

    upb(1)=upbound_new1

    lb=lbound(p1)
    ub=ubound(p1)

     allocate(p2(lb(1):upb(1)),stat=info)
    if (info == 0) then
       p2=0
       do i=lb(1),min(ub(1),upb(1))
          p2(i)=p1(i)
       end do
       deallocate(p1,stat = info)
       if (info == 0) then
          p1 => p2
       else
          info = ERR_MEMORY_DEALLOC
       end if
    else
       info = ERR_MEMORY_ALLOC
    end if

  end subroutine iexpand1


  ! given a row order and the row grouping,
  ! output the column ordering
  ! so that using the row and column ordering,
  ! the matrix will become
  ! bordered block diagonal.
  subroutine get_column_order(matrix,num_blocks,max_nblks,&
       row_order,row_sta,column_order,column_sta,info)
     ! matrix: the matrix to be ordered
    type (ZD11_type), target, intent (in) :: matrix

    ! max_nblks: number of blocks to parttion the matrix into
    integer (kind = myint), intent (in) :: max_nblks

    ! num_blocks: number of blocks actually achieved
    integer (kind = myint), intent (in) :: num_blocks

    ! row ordering: row_order(i) is the original row index
    ! of the i-row of the reordered matrix
    integer (kind = myint), intent (in), dimension (:) :: row_order

    ! row_sta: the starting row index for each block.
    ! row_sta(i) is the starting row index of the i-th block.
    integer (kind = myint), intent (in), dimension (:)  :: row_sta

    ! column order: column_order(i) is the original column index
    ! of the i-column of the reordered matrix
    integer (kind = myint), intent (out), dimension (:) :: &
         column_order

    ! column_sta: the starting column index for each block.
    ! column_sta(i) is the starting column index of the i-th block.
    integer (kind = myint), intent (out), dimension (:) :: column_sta

    ! info: < 0 if allocation failed
    integer (kind = myint), intent (out) :: info

    ! m : column size
    integer (kind = myint) :: m

    ! i,j,k: loop index
    ! col: column index
    ! ierr: error tag
    integer (kind = myint) :: i,j,k,col,ierr

    ! counter: counter(j) = 0 initially. If
    !      counter(j) < 0 --- j is a shared column (netcut)
    !      counter(j) = 0 --- j is an empty column
    !      counter(j) = i --- column j only appear in block i
    integer (kind = myint), allocatable :: counter(:)

    ! iwork: temperatory working array which contain
    ! the number of columns in each block,
    ! with column_sta(max_nblks+1) the number of columns in the
    ! interface block
    integer (kind = myint), allocatable :: iwork(:)
    integer (kind = myint) :: irow

    info = 0
    m = matrix%n

    allocate(counter(m),stat = ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_ALLOC
       return
    end if

    counter(:) = 0
    ! loop over each block and increment a counter
    do i = 1, num_blocks
       do j = row_sta(i), row_sta(i+1)-1
          irow = row_order(j)
          do k = matrix%ptr(irow), matrix%ptr(irow+1) - 1
             col = matrix%col(k)
             if (counter(col) == 0) then
                counter(col) = i
             else if (counter(col) > 0.and. counter(col) /= i) then
                ! if this column has been see in another
                ! block before,
                ! then put in interface block
                counter(col) = -counter(col)
             end if
          end do
       end do
    end do

    allocate(iwork(max_nblks+1), stat = ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_ALLOC
       return
    end if

    ! iwork now will contain the number of columns in each block,
    ! with iwork(max_nblks+1) the number of columns in the
    ! interface block
    iwork(:) = 0
    do i = 1, m
       col = counter(i)
       if (col>0) then
          iwork(col) = iwork(col) + 1
       else if (col == 0) then ! put empty column at the last block
!          iwork(max_nblks+1) = iwork(max_nblks+1) + 1
          info = WARN_EMPTY_COL
          iwork(max_nblks) = iwork(max_nblks) + 1
       end if ! do nothing for the interface block
    end do

    ! set column_sta to be the starting column index for each block.
    ! column_sta(i) is the starting column index of the i-th block.
    column_sta(1) = 1
    do i = 1, max_nblks
       column_sta(i+1) = column_sta(i) + iwork(i)
    end do

    ! now collect column indices that belongs to each block
    do i = 1, m
       col = counter(i)
       if (col > 0) then
          column_order(column_sta(col)) = i
          column_sta(col) = column_sta(col) + 1
       else if (col == 0) then ! put empty column at the last block
          column_order(column_sta(max_nblks)) = i
          column_sta(max_nblks) = column_sta(max_nblks) + 1
       else
          column_order(column_sta(max_nblks+1)) = i
          column_sta(max_nblks+1) = column_sta(max_nblks+1) + 1
       end if
    end do

    ! set again column_sta to be the starting column
    ! index for each block.
    ! column_sta(i) is the starting column index of the i-th block.
    column_sta(1) = 1
    do i = 1, max_nblks
       column_sta(i+1) = column_sta(i) + iwork(i)
    end do

    deallocate(counter,iwork,stat = ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_ALLOC
       return
    end if

  end subroutine get_column_order


! misc routines used by KL in calculating cuts
! and managing boundary list etc

! initialise the in and out degree used in KL
  subroutine KL_init(kl_obj,matrix,graph,where,col_wgt,&
       out_deg,in_deg,mincut,ierr)
    ! kl_obj: the matrix and connectivity graph
    ! KL is working on
    type (kl_type), intent (out) :: kl_obj

    ! the matrix
    type (ZD11_type), target, intent (in) :: matrix
    ! the connectivity of rows
    type (ZD11_type), target, intent (in) :: graph


    ! where each nodes are at the moment (either 0 or 1)
    integer (kind = myint), intent (in) :: where(:)

    ! the column weight
    integer (kind = myint), dimension (:), target :: col_wgt

    ! the number of connections of a node to its
    ! subdomain and the opposite
    ! subdomain
    integer (kind = myint), pointer :: in_deg(:),out_deg(:)

    ! the current net cut
    integer (kind = myint), intent (out) :: mincut

    ! ierr: error flag. not zero only if allocation failed
    integer (kind = myint), intent (out) :: ierr

    ! size of the problem (number of rows)
    integer (kind = myint) :: nvtx

    ! loop index
    integer (kind = myint) :: i

    ! matrix row pointer
    integer (kind = myint), dimension (:), pointer :: the_row

    ierr = 0

    ! initialise the foot_print
    kl_obj%KL_matrix => matrix
    KL_obj%KL_graph => graph

    ! initialise column weight
    KL_obj%colwgt => col_wgt

    allocate(KL_obj%foot_print(matrix%n,0:1), stat = ierr)
    if (ierr /= 0) return

    KL_obj%foot_print(:,:) = 0
    do i = 1, matrix%m
!       the_row => mc65_matrix_getrow(matrix,i)
       call mc65_matrix_getrow(matrix,i,the_row)
       if (where(i) == 0) then
          KL_obj%foot_print(the_row(:),0) = &
               KL_obj%foot_print(the_row,0) + 1
       else if (where(i) == 1) then
          KL_obj%foot_print(the_row(:),1) = &
               KL_obj%foot_print(the_row,1) + 1
       end if
    end do

    nvtx = kl_obj%KL_matrix%m
    allocate(in_deg(nvtx),out_deg(nvtx), stat = ierr)
    if (ierr /= 0) return

    do i = 1, nvtx
       call degree_calc(kl_obj,i,where(i),out_deg(i),in_deg(i))
    end do
    mincut = 0
    do i = 1, size(KL_obj%foot_print,1)
       if (KL_obj%foot_print(i,0) > 0 .and. &
            KL_obj%foot_print(i,1) > 0 ) then
          mincut = mincut + KL_obj%colwgt(i)
       end if
    end do

  end subroutine KL_init


  subroutine KL_deallocate(kl_obj,out_deg, in_deg, ierr)
    ! kl_obj: the matrix and connectivity graph
    ! KL is working on
    type (kl_type), intent (inout) :: kl_obj

    ! the number of connections of a node to its
    ! subdomain and the opposite
    ! subdomain
    integer (kind = myint), pointer :: in_deg(:),out_deg(:)

    ! ierr: error flag, nonzero only if deallocation failed.
    integer (kind = myint), intent (out) :: ierr

    ierr = 0

    deallocate(out_deg, in_deg, stat = ierr)
    if (ierr /= 0) return

    deallocate(KL_obj%foot_print, stat = ierr)
    if (ierr /= 0) return

    nullify(kl_obj%KL_matrix,KL_obj%KL_graph)

  end subroutine KL_deallocate


  subroutine boundary_init(nbd,bd_mask)
    integer (kind = myint), intent (out) :: nbd
    integer (kind = myint), intent (out), dimension (:) :: bd_mask
    bd_mask(:) = -1
    nbd = 0
  end subroutine boundary_init


  subroutine boundary_delete(node,nbd,bd_list,bd_mask)
    ! delete node from the boundary list
    integer (kind = myint), intent (in) :: node
    ! number of boundary nodes
    integer (kind = myint), intent (inout) :: nbd
    ! the list and mask (also act as pointer)
    integer (kind = myint), intent (inout), dimension (:) :: &
         bd_mask,bd_list

!!$    if (bd_mask(node) == -1) then
!!$       stop "the node to be deleted from the bdry is not in there!"
!!$    end if

    ! move the last node in the boundary list to
    ! the position of this deleted node
    bd_list(bd_mask(node)) = bd_list(nbd)
    bd_mask(bd_list(nbd)) = bd_mask(node)
    nbd = nbd - 1
    bd_mask(node) = -1
  end subroutine boundary_delete


  subroutine boundary_add(node,nbd,bd_list,bd_mask)
    ! add node to the boundary list
    integer (kind = myint), intent (in) :: node
    ! number of boundary nodes
    integer (kind = myint), intent (inout) :: nbd
    ! the list and mask (also act as pointer)
    integer (kind = myint), intent (inout), dimension (:) :: &
         bd_mask,bd_list

!!$    if (bd_mask(node) /= -1) then
!!$       stop "the node to be added to the bdry is already in there!"
!!$    end if

    ! add the node
    nbd = nbd + 1
    bd_list(nbd) = node
    bd_mask(node) = nbd
  end subroutine boundary_add


! return a list of neighbors in the row connectivity table
!  function get_neighbor(kl_obj,i) result(neighbor)
  subroutine get_neighbor(kl_obj,i,neighbor)
    ! kl_obj: the matrix and connectivity graph
    ! KL is working on
    type (kl_type), intent (in), target :: kl_obj
    integer (kind = myint), intent (in) :: i
    integer (kind = myint), dimension (:), pointer :: neighbor

!   neighbor => mc65_matrix_getrow(KL_obj%KL_graph,i)
    call mc65_matrix_getrow(KL_obj%KL_graph,i,neighbor)
  end subroutine get_neighbor


! return a row in the matrix
!  function get_row(kl_obj,i) result(neighbor)
   subroutine get_row(kl_obj,i,neighbor)
    ! kl_obj: the matrix and connectivity graph
    ! KL is working on
    type (kl_type), intent (in) :: kl_obj
    integer (kind = myint) :: i
    integer (kind = myint), dimension (:), pointer :: neighbor

!    neighbor => mc65_matrix_getrow(kl_obj%KL_matrix,i)
    call mc65_matrix_getrow(KL_obj%KL_matrix,i,neighbor)
  end subroutine get_row


  subroutine footprint_change(kl_obj,node,from)
    ! kl_obj: the matrix and connectivity graph
    ! KL is working on
    type (kl_type), intent (inout) :: kl_obj
    ! calculate the cut increment due to the move
    ! of node from subdomain "from"
    ! to the other side
    ! also update the "foot_print" of the two subdomain

    ! the node to move from "from"
    integer (kind = myint), intent (in) :: node
    integer (kind = myint), intent (in) :: from
    ! the increment of cut
!    integer (kind = myint), intent (out) :: cut_inc

    !local variables:
    ! neighbors of node
    integer (kind = myint), dimension (:), pointer :: neighbor


    ! the domain to move to
    integer (kind = myint) :: to

    to = mod(from + 1, 2)
!    neighbor => get_row(kl_obj,node)
    call get_row(kl_obj,node,neighbor)

    KL_obj%foot_print(neighbor,from) = &
         KL_obj%foot_print(neighbor,from) - 1
    KL_obj%foot_print(neighbor,to) = &
         KL_obj%foot_print(neighbor,to) + 1

  end subroutine footprint_change


  subroutine degree_calc(kl_obj,node,domain,my_out_deg,my_in_deg)
    ! kl_obj: the matrix and connectivity graph
    ! KL is working on
    type (kl_type), intent (inout) :: kl_obj
    ! calculate the degree of node "node" in subdomain "domain"

    ! the node in "domain"
    integer (kind = myint), intent (in) :: node
    integer (kind = myint), intent (in) :: domain
    ! new out and in degree
    integer (kind = myint), intent (out) :: my_out_deg,my_in_deg

    !local variables:
    ! neighbors of node
    integer (kind = myint), dimension (:), pointer :: neighbor
    ! loop index
    integer (kind = myint) :: i
    ! the domain on the other side
    integer (kind = myint) :: other

     my_out_deg = 0; my_in_deg = 0
    other = mod(domain + 1, 2)
!    neighbor => get_row(kl_obj,node)
    call get_row(kl_obj,node,neighbor)
    do i = 1, size(neighbor)
       if (KL_obj%foot_print(neighbor(i),other) /= 0) then
          my_out_deg = my_out_deg + KL_obj%colwgt(neighbor(i))
       end if
       ! if this column index is not the only one in the subdomain!
       if (KL_obj%foot_print(neighbor(i),domain) > 1) then
          my_in_deg = my_in_deg +  KL_obj%colwgt(neighbor(i))
       end if
    end do

  end subroutine degree_calc


  subroutine queue_init(queue,maxnodes,maxgain,mingain,ierr)
    ! initialise the queue
    type (queue_type), intent (out) :: queue
    integer (kind = myint), intent (in) :: maxnodes, maxgain,mingain
    integer (kind = myint), intent (out) :: ierr
    integer (kind = myint) :: i

    ierr = 0
! COMMENTED OUT AS THIS IS NOT GOING TO HAPPEN
! WITH IN MC66
!    if (mingain > maxgain) then
!       ierr = -2
!       return
!    end if

! allocate buckets
    allocate(queue%buckets(mingain:maxgain), stat = ierr)
    if (ierr /= 0) then
       return
    end if
    queue%low_bound = mingain
    queue%up_bound = maxgain

! no maxgain yet (all buckets empty)
    queue%maxgain = queue%low_bound - 1

! no nodes
    queue%num_nodes = 0

! nullify the buckets
    do i = queue%low_bound, queue%up_bound
       nullify(queue%buckets(i)%current_node)
    end do

! allocate node pointers which say where in the buckets
! it locates
    allocate(queue%mynode(maxnodes), stat = ierr)
    if (ierr /= 0) then
       return
    end if
    queue%maxnodes = maxnodes

  end subroutine queue_init


  subroutine queue_deallocate(queue,ierr)
    ! deallocate a queue
    type (queue_type), intent (inout) :: queue
    ! if not zero: deallocation failed
    integer (kind = myint), intent (out) :: ierr

    ierr = 0

! deallocate buckets array
    deallocate(queue%buckets, stat = ierr)
    if (ierr /= 0) then
       return
    end if

! deallocate node pointers which say where in the buckets
! it locates
    deallocate(queue%mynode, stat = ierr)
    if (ierr /= 0) then
       return
    end if

  end subroutine queue_deallocate


  subroutine queue_add(queue,id,gain)
    ! add into the queue a node with gain "gain"
    type (queue_type), intent (inout) :: queue
    integer (kind = myint), intent (in) :: id,gain

    type (list_node_type), pointer :: new_node

!!$    if (gain > queue%up_bound.or. gain < queue%low_bound) then
!!$       write(*,*) "GAIN = ",GAIN, " BOUND = ",&
!!$          queue%low_bound,queue%up_bound
!!$       stop "in queue_add, gain out of range"
!!$    end if
!!$
!!$    if (id < 0 .or. id > queue%maxnodes) then
!!$       stop "in queue_add, id out of range"
!!$    end if

    ! node add one
    queue%num_nodes = queue%num_nodes + 1

    ! update the current maxgain
    if (gain > queue%maxgain) queue%maxgain = gain

    new_node => queue%mynode(id)

    ! fill its content
    new_node%id = id

    ! add this node at the front of the queue
    new_node%next => queue%buckets(gain)%current_node
    nullify(new_node%prev)

    ! creat the link provided that new_node is not
    ! the only one in the list
    if (associated(new_node%next)) then
       new_node%next%prev => new_node
    end if

    ! the bucket now start from this node
    queue%buckets(gain)%current_node => new_node

  end subroutine queue_add


  subroutine queue_delete(queue,id,gain)
    ! delete the node id from bucket gain
    type (queue_type), intent (inout) :: queue
    integer (kind = myint), intent (in) :: id,gain


    type (list_node_type), pointer :: new_node

!!$    if (gain > queue%up_bound.or. gain < queue%low_bound) then
!!$       stop "in queue_delete, gain out of range"
!!$    end if
!!$
!!$    if (id < 0 .or. id > queue%maxnodes) then
!!$       stop "in queue_delete, id out of range"
!!$    end if

! node less one
    queue%num_nodes = queue%num_nodes - 1

    new_node => queue%mynode(id)

! if the node to be deleted is not the last in the bucket
    if (associated(new_node%next)) then
       new_node%next%prev => new_node%prev
    end if

! if the node to be deleted is not the first in the bucket
    if (associated(new_node%prev)) then
       new_node%prev%next => new_node%next
! if it is the first then bucket will start from the next
    else
       queue%buckets(gain)%current_node => new_node%next
    end if

! update the max gain
    if (gain == queue%maxgain) then
       ! find the highest none empty buckets
       if (queue%num_nodes == 0) then
          queue%maxgain = queue%low_bound - 1
          return
       end if
       do while (.not.&
            associated(queue%buckets(queue%maxgain)%current_node))
          queue%maxgain = queue%maxgain - 1
       end do
    end if
  end subroutine queue_delete


   subroutine queue_getmax(queue,id,gain)
    ! give the node with the largest gain in the queue
    ! last in first out!
    type(queue_type), intent(in) :: queue
    integer (kind = myint), intent (out) :: id,gain

    if (queue%num_nodes == 0) then
       gain = 0
       id = -1
       return
    end if
    gain = queue%maxgain
    id = queue%buckets(gain)%current_node%id
  end subroutine queue_getmax


  subroutine kl(plevel,ostream,estream,wstream,&
       matrix,graph,row_wgt,col_wgt,&
       where,target_nnodes,&
       aggressive,unbalance,mincut,nnodes,info)

    ! PLEVEL: print level
    integer (kind = myint), intent (in) :: plevel
    ! ostream: output stream
    ! estream: error stream
    ! wstream: warning stream
    integer (kind = myint), intent (in) :: ostream,estream,wstream

    ! matrix: the matrix
    type (ZD11_type), intent (in), target :: matrix
    ! graph: the graph of row connectivity
    type (ZD11_type), intent (in), target :: graph

    ! the weight of each row, which is the number of rows the coarse
    ! row stands for. Also the weight of each column,
    ! which is the number of original columns this column stands for
    ! (columns are merged if they have exactly the same pattern,
    ! deleted if they have only one element)
    integer (kind = myint), dimension (:), intent (in) :: row_wgt
    integer (kind = myint), dimension (:), intent (in), target :: col_wgt

    ! IN: where each nodes are at the moment (either 0 or 1)
    ! OUT: where each node will be after KL
    integer (kind = myint), intent (inout) :: where(matrix%m)

    ! target_nnodes: the targeting size of the sum of row (node)
    ! weights in each subdomain
    ! sum(target_nnodes) = sum(row_wgt)
    integer (kind = myint), intent (in) :: target_nnodes(0:1)

    ! whether aggressive (>0) or not (<0)
    real (kind = myreal), intent (in) :: aggressive

    ! maximum allowed imbalance in a bisection in percentage term/
    ! this is defined as
    ! abs(nnodes(1)-target_nnodes(1))/sum(target_nnodes)
    !  = abs(nnodes(0)-target_nnodes(0))/sum(target_nnodes)
    !     <= unbalance
    real (kind = myreal), intent (in) :: unbalance

    ! the minimum cut after KL
    integer (kind = myint), intent (out) :: mincut

     ! sum of node weights in each of the two subdomain.
    ! sum(nnodes) = sum(target_nnodes)
    integer (kind = myint), intent (out) :: nnodes(0:1)

    ! info : information tag
    integer (kind = myint), intent (out) :: info

    ! total number of nodes on this level (= matrix%m)
    integer (kind = myint) :: nvtx

    ! maximum possible gain of moving a node
    integer (kind = myint) :: maxgain

    ! the number of connections of a node to its subdomain
    ! and the opposite subdomain
    integer (kind = myint), pointer :: in_deg(:),out_deg(:)


    ! the queue which holds boundary nodes along the
    ! current partitioning one for each subdomain
    type (queue_type) :: queue(0:1)

    ! the array which says whether a node is locked
    integer (kind = myint) :: locked(matrix%m)

    ! the current net cut, and the best so far in the round
    integer (kind = myint) :: mincut_cur, mincut_cur_best

    ! the current domain adn the destined domain
    integer (kind = myint) :: to,from

    ! gain and node id with the best gain
    integer (kind = myint) :: gain,best, gain1

    ! the array which records the previous moves, used
    ! to undo them when needed
    integer (kind = myint) :: swaps(matrix%m)

    ! list of boundary nodes and a mask which says
    ! whether a given node
    ! is in the boundary (and if so where in the list it is)
    integer (kind = myint) :: bd_list(matrix%m),bd_mask(matrix%m)

    ! the number of nodes in the boundary
    integer (kind = myint) :: nbd

    ! a neigbbor and its in out degree
    integer (kind = myint) :: friend,my_in_deg,my_out_deg

    ! best inbalance and current imbalance
    integer (kind = myint) :: balance, balance_cur, balance_cur_best

    ! the last swap that has give improvement
    integer (kind = myint) :: last_goodswp

    ! neighbor indices
    integer (kind = myint), pointer, dimension (:) :: neighbor

    ! loop index
    integer (kind = myint) :: i,j,iswaps,ikl

    ! cut change after a swap
    integer (kind = myint) :: cut_inc

    ! maximum of extra swapping allowed that does
    ! not give improvement
    integer (kind = myint) :: max_extra_swaps

    ! maximum allowed KL loop
    integer (kind = myint) :: max_kl

    ! average row weights (node weight)
    !    = sum(target_nnodes)/2
    real (kind = myreal) :: avg_nnodes

    ! kl_obj: the matrix and connectivity graph
    ! KL is working on
    type (kl_type) :: kl_obj

    ! ierr: error flags
    integer (kind = myint) :: ierr

    info = 0

    nvtx = matrix%m
    avg_nnodes = 0.5*sum(target_nnodes)

    call KL_init(kl_obj,matrix,graph,where,col_wgt,&
         out_deg,in_deg,mincut,ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_ALLOC
       call monet_print_message(info,estream,wstream,&
               "in KL_init")
       return
    end if

    ! very timid
    !    max_kl = 1
    !    max_extra_swaps = 0

    if (aggressive <= 0) then
       ! moderate
       max_kl = 10
       max_extra_swaps = max(15,min(nvtx/100,100))
    else
       ! aggressive
       max_kl = 1000
       max_extra_swaps = huge(0)
    end if
    if (plevel > 1) write(ostream,*) 'unbalance = ',&
         unbalance,' max_kl = ',max_kl

    if (plevel > 1) write(ostream,*) "in KL"

    ! work out the number of nodes in each of the partitioning
    nnodes = 0
    do i = 1, nvtx
       nnodes(where(i)) = nnodes(where(i)) + row_wgt(i)
    end do
    if (plevel > 1) then
       write(ostream,*)
       write(ostream,*) " ==============="
       write(ostream,*) "before KL mincut = ", &
         mincut," problem size",matrix%m," nnodes = ",nnodes
    end if


    ! imbalance
    balance = abs(nnodes(0) - target_nnodes(0))

    ! boundary node list initialisation
    call boundary_init(nbd,bd_mask)
    do i = 1, nvtx
       ! find the boundary nodes and their gain,
       ! a boundary node is one with
       ! external link, or an isolated node
       if (out_deg(i) > 0 .or. out_deg(i) + in_deg(i) == 0) then
          call boundary_add(i,nbd,bd_list,bd_mask)
       end if
    end do

    ! maximum gain is the sum of all super-column weights
    maxgain = sum(col_wgt)

    ! start of the KL sweeps
    KLLOOP: do ikl = 1, max_kl

       ! initialisation ======================
       ! initialise the queue which records the boundary
       ! nodes and their gains
       call queue_init(queue(0),nvtx,MAXGAIN = maxgain, &
            MINGAIN = -maxgain, ierr = ierr)
       if (ierr /= 0) then
          info = ERR_MEMORY_ALLOC
          call monet_print_message(info,estream,wstream,&
               "in queue_init")
          return
       end if
       call queue_init(queue(1),nvtx,MAXGAIN = maxgain, &
            MINGAIN = -maxgain, ierr = ierr)
       if (ierr /= 0) then
          info = ERR_MEMORY_ALLOC
          call monet_print_message(info,estream,wstream,&
               "in queue_init")
          return
       end if
       if (nbd == 0) exit KLLOOP

       do i = 1, nbd
          j = bd_list(i)
          ! find the boundary nodes and their gain
          gain = out_deg(j) - in_deg(j)
          call queue_add(queue(where(j)),j,gain)
       end do

       ! initialise the locking array
       locked = -1

       ! current mincut and current best mincut initialised,
       ! also initialise the best move so far
       mincut_cur = mincut
       mincut_cur_best = mincut
       balance_cur_best = balance
       last_goodswp = 0
       ! ================ end initialisation

!       write(ostream,*) 'before SWAPPING'

       ! do a number of swaps
       SWAPPING: do iswaps = 1, nvtx + 1
          ! decide which domain to choose the node from:
          ! get from the domain with more nodes if quite unbalanced,
          ! else get it from the domain with best gain
!        if (abs(nnodes(0) - target_nnodes(0))/real(nvtx) &
!            >= unbalance) then
         if (abs(nnodes(0) - target_nnodes(0))/avg_nnodes &
              >= unbalance) then
             if (nnodes(0) < target_nnodes(0)) then
                from = 1; to = 0
             else
                from = 0; to = 1
             end if
          else
             call queue_getmax(queue(0),best,gain)
             call queue_getmax(queue(1),best,gain1)
             if (gain1 > gain) then
                from = 1; to = 0
             else
                from = 0; to = 1
             end if
          end if

          ! get a node in domain "from" with the largest
          ! gain and delete this node from the queue
          call queue_getmax(queue(from),best,gain)
          if (best == -1) then
             exit SWAPPING
          end if
          if (gain < 0 .and. iswaps - last_goodswp > &
               max_extra_swaps) then
             exit SWAPPING
          endif
          ! do not swap if this will result in one
          ! domain having zero weight!
          if (nnodes(from) - row_wgt(best) == 0) exit SWAPPING

          call queue_delete(queue(from),best,gain)

          ! move it to the other side
          where(best) = to
          nnodes(to) = nnodes(to) + row_wgt(best)
          nnodes(from) = nnodes(from) - row_wgt(best)

          ! record this move and lock this node
          swaps(iswaps) = best
          locked(best) = iswaps

          ! update the current mincut by adding the cut
          ! increment resulted from this move
          call footprint_change(kl_obj,best,from)
          cut_inc = in_deg(best) - out_deg(best)
          mincut_cur = mincut_cur + cut_inc
          if (plevel > 2) then
             write(ostream,'(a,i6,a,i6,a,i6,a,i6)') &
                  "move ",best," to ",to," gain = ",gain,&
                  " cutinc = ",cut_inc," mincut = ",mincut_cur
             write(ostream,'(a,2i6,a,2i6)') 'domain nodes from ',&
               nnodes(to) - row_wgt(best),&
               nnodes(from) + row_wgt(best), &
               ' to ',nnodes(to),nnodes(from)
          end if

          ! if this is a good move (negative cut_inc,
          ! or improved balance)
          ! over the best so far then record it
          ! so that we may go back to it when necessary
          balance_cur = abs(nnodes(0) - target_nnodes(0))
          !          if (mincut_cur < mincut_cur_best .or. &
          if ((mincut_cur < mincut_cur_best.and.&
               balance_cur/avg_nnodes <= unbalance).or. &
               (mincut_cur == mincut_cur_best .and. &
               balance_cur < balance_cur_best)) then
             last_goodswp = iswaps
             mincut_cur_best = mincut_cur
             balance_cur_best = balance_cur
             if (plevel > 1) write(ostream,*) &
                  "current best = ",mincut_cur_best,&
                  " balance = ",balance_cur_best
          end if

          ! update the in and out degree itself:
          ! if not isolated and external degree
          ! changed to zero, delete from boundary
          call degree_calc(kl_obj,best,to,my_out_deg,my_in_deg)
          out_deg(best) = my_out_deg; in_deg(best) = my_in_deg

          if (my_out_deg == 0 .and. my_in_deg > 0) then
             call boundary_delete(best,nbd,bd_list,bd_mask)
          end if

          ! update the in and out degree neighbors
 !         neighbor => get_neighbor(kl_obj,best)
          call get_neighbor(kl_obj,best,neighbor)
          do i = 1, size(neighbor)
             friend = neighbor(i)
             call degree_calc(kl_obj,friend,where(friend),&
                  my_out_deg,my_in_deg)

             ! if this was a boundary node
             if (bd_mask(friend) /= -1) then
                ! and not a boundary node any more then
                ! remove from queue if were in the queue
                ! also delete from boundary
                if (my_out_deg <= 0) then
                   if (locked(friend) == -1) then
                      call queue_delete(queue(where(friend)),friend,&
                           out_deg(friend) - in_deg(friend))
                   end if
                   call boundary_delete(friend,nbd,bd_list,bd_mask)
                   ! else update the position in the
                   ! queue since gain changed, provided it is in
                   ! the queue (not locked) and the gain has changed
                else if (locked(friend) == -1.and. &
                     out_deg(friend) - in_deg(friend) /= &
                     my_out_deg - my_in_deg) then
                   call queue_delete(queue(where(friend)),friend,&
                        out_deg(friend) - in_deg(friend))
                   call queue_add(queue(where(friend)),friend,&
                        my_out_deg - my_in_deg)
                end if
             else
                ! if was not a boundary node and not locked,
                ! but now become a boundary node,
                ! add into queue
                if (my_out_deg > 0) then
                   if (locked(friend) == -1) then
                      call queue_add(queue(where(friend)),friend,&
                           my_out_deg - my_in_deg)
                   end if
                   ! add to the boundary
                   call boundary_add(friend,nbd,bd_list,bd_mask)
                end if
             end if
             out_deg(friend) = my_out_deg; in_deg(friend) = my_in_deg
          end do

       end do SWAPPING

       ! undo bad moved

       UNROLL: do i = iswaps - 1, last_goodswp + 1, -1
          best = swaps(i)
          from = where(best)
          to = mod(from + 1,2)
          where(best) = to
          nnodes(to) = nnodes(to) + row_wgt(best)
          nnodes(from) = nnodes(from) - row_wgt(best)

          ! update the cuts
          call footprint_change(kl_obj,best,from)
          cut_inc = in_deg(best) - out_deg(best)

          mincut_cur = mincut_cur + cut_inc

          ! update the in and out degree of itself
          call degree_calc(kl_obj,best,to,my_out_deg,my_in_deg)
          out_deg(best) = my_out_deg; in_deg(best) = my_in_deg

          ! if was not boundary but now become one
          if (bd_mask(best) == -1 .and. my_out_deg > 0) then
             call boundary_add(best,nbd,bd_list,bd_mask)
             ! or if was in the boundary and now not one
! COMMENTED OUT AS THIS IS NOT GOING TO HAPPEN
! WITH IN MC66. IN FACT SHOULD NOT HAPPEN
! AT ALL SINCE IF BEST IS NOW AN INTERNAL NODE,
! THEN IT SHOULD NOT BE MOVE DURING SWAPPING
! IN THE FIRST PLACE: ONLY BOUNDARY NODES
! ARE SWAPPED DURING SWAPPING PHASE
!          else if (bd_mask(best) /= -1 .and. my_out_deg == 0 &
!               .and. my_in_deg > 0) then
!             call boundary_delete(best,nbd,bd_list,bd_mask)
          end if

          ! update the in and out degree of the neighbors
!          neighbor => get_neighbor(kl_obj,best)
          call get_neighbor(kl_obj,best,neighbor)
          if (plevel>10) write(ostream,*) &
               "number of neighbors = ",size(neighbor)
          do j = 1, size(neighbor)
             friend = neighbor(j)
             call degree_calc(kl_obj,friend,where(friend),&
                  my_out_deg,my_in_deg)
             ! if this was a boundary node
             if (bd_mask(friend) /= -1) then
                ! and not a boundary node any more then
                if (my_out_deg <= 0) then
                   ! delete from boundary
                   call boundary_delete(friend,nbd,bd_list,bd_mask)
                end if
             else
                ! if was not a boundary node but now become
                ! a boundary node,
                if (my_out_deg > 0) then
                   ! add to the boundary
                   call boundary_add(friend,nbd,bd_list,bd_mask)
                end if
             end if
             out_deg(friend) = my_out_deg; in_deg(friend) = my_in_deg
          end do

       end do UNROLL

!       write(ostream,*) 'after UNROLL'

       balance_cur = abs(nnodes(0) - target_nnodes(0))
!!$       if (mincut_cur /= mincut_cur_best .or. &
!!$             balance_cur /= balance_cur_best) then
!!$          write(ostream,*) mincut_cur,mincut_cur_best
!!$          write(ostream,*) balance_cur, balance_cur_best
!!$          write(ostream,*) " nnodes = ",nnodes
!!$          stop " inconsistency in the cut after rolling back!"
!!$       end if

       ! deallocate queue
       call queue_deallocate(queue(0),ierr)
       if (ierr /= 0) then
          info = ERR_MEMORY_DEALLOC
          call monet_print_message(info,estream,wstream,&
               "in queue_deallocate")
          return
       end if
       call queue_deallocate(queue(1),ierr)
       if (ierr /= 0) then
          info = ERR_MEMORY_DEALLOC
          call monet_print_message(info,estream,wstream,&
               "in queue_deallocate")
          return
       end if

       ! KL finished if neither cut improved nor balance improved,
       if (plevel > 1) write(ostream,&
            "('best_cut = ',i10,' swps = ',i10,' nunrol =',i10)") &
            mincut,iswaps,iswaps - last_goodswp
       if (mincut <= mincut_cur_best .and. &
            balance_cur >= balance) then
          exit KLLOOP
       else
          mincut = mincut_cur_best
          balance = balance_cur_best
       end if

    end do KLLOOP
!       write(ostream,*) 'where = ',where

    if (plevel > 1) write(ostream,'(a,i10,i10,a,i10,a,i10)') &
         ' KL finish with domain nnodes = ',nnodes, &
         ' cut = ',mincut," ikl = ",ikl

    call KL_deallocate(kl_obj,out_deg,in_deg,ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_DEALLOC
       call monet_print_message(info,estream,wstream,&
               "in KL_deallocate")
       return
    end if

  end subroutine kl


  subroutine level_print(ostream,title1,level,title2,res)
    ! print the title in a indented fasion depend on which level.
    ! will give output like this
    !   ===== title on level   1 =====
    !      ===== title on level   2 =====
    !         ===== title on level   3 =====
    !            ===== title on level   4 =====
    !               ===== title on level   5 =====
    !                  ===== title on level   6 =====
    !                      .......
    !
    integer (kind = myint) :: ostream
    character (len = *), intent(in) :: title1
    integer (kind = myint), intent(in) :: level
    real (kind = myreal), optional, intent (in) :: res
    character (len = *), optional, intent(in) :: title2
    character (len=80) fmt
    integer (kind = myint) :: char_len1,char_len2

    char_len1=len_trim(title1)

    if (present(res).and.present(title2)) then
       char_len2=len_trim(title2)
       write(fmt, &
            "('(',i4,'X,''===== '',a',i4,',i4,a',i4,',g14.3,'' ====='')')") &
            level*3, char_len1,char_len2
       write(ostream,fmt) title1,level,title2,res
    else if (present(res)) then

       write(fmt, &
            "('(',i4,'X,''===== '',a',i4,',i4,3h is,g14.3,'' ====='')')") &
            level*3, char_len1
       write(ostream,fmt) title1,level,res
    else
       write(fmt,"('(',i4,'X,''===== '',a',i4,',i4,'' ====='')')") &
            level*3, char_len1
       write(ostream,fmt) title1,level
    end if
  end subroutine level_print


  subroutine matrix_to_graph(estream,wstream,matrix,graph,info,col_wgt)
    ! subroutine which gives a graph  A*A_trans
    ! faster!
    ! estream, wstream: error and warning stream
    integer (kind = myint) :: estream, wstream
    type (ZD11_type), intent (in), target :: matrix
    type (ZD11_type), intent (out), target :: graph

     ! the super-column weight.
    integer (kind = myint), dimension (:), intent (in), optional :: &
         col_wgt
    ! info: information tag
    integer (kind = myint), intent (out) :: info

    type (ZD11_type), target :: graph_temp
    integer (kind = myint) :: ierr

    info = 0

    call mc65_matrix_transpose(matrix,graph_temp,ierr)
    if (ierr /= 0) then
       if (ierr == -1) info = ERR_MEMORY_ALLOC
       if (ierr == -2) info = ERR_MEMORY_DEALLOC
       call monet_print_message(info,estream,wstream,&
            "in matrix_to_graph mc65_matrix_transpose")
       return
    end if

    if (present(col_wgt)) then
       call mc65_matrix_multiply_graph(matrix,graph_temp,graph,ierr,&
            col_wgt = real(col_wgt,kind = myreal))
    else
       call mc65_matrix_multiply_graph(matrix,graph_temp,graph,ierr)
    end if
    if (ierr /= 0) then
       if (ierr == -1) info = ERR_MEMORY_ALLOC
       if (ierr == -2) info = ERR_MEMORY_DEALLOC
       call monet_print_message(info,estream,wstream,&
            "in matrix_to_graph mc65_matrix_multiply_graph")
       return
    end if


    call mc65_matrix_destruct(graph_temp,ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_DEALLOC
       call monet_print_message(info,estream,wstream,&
            "in matrix_to_graph mc65_matrix_destruct")
       return
    end if

  end subroutine matrix_to_graph


  recursive subroutine multilevel(grid,controller,mglevel,seed,&
       target_nnodes,info)

    ! this level of matrix (grid)
    type (multigrid), intent (inout), target :: grid


    ! controller: the controller containing control parameters
    type (mc66_control), intent (in) :: controller


    ! number of multilevel allowed in this multilevel run,
    ! this may be smaller than in the controller%mglevel
    ! so as to allow termination of the multilevel
    ! in case coarsening is not making enough reduction in
    ! problem size
    integer (kind = myint), intent (in) :: mglevel
    ! random number seed
    type (fa14_seed), intent (inout) :: seed

    ! target_nnodes: the targeting size of the sum of row (node)
    ! weights in each subdomain
    ! sum(target_nnodes) = sum(row_wgt)
    integer (kind = myint), intent (in) :: target_nnodes(0:1)

    ! info: integer scalar of intent (inout). Info tag
    integer (kind = myint), intent (out) :: info

    ! cgrid: the coarse level grid
    type (multigrid), pointer :: cgrid

    ! p: the coarse grid prolongator
    type (ZD11_type), pointer :: p

    ! r: the fine grid restrictor
    type (ZD11_type), target :: r

    ! fwhere, cwhere: the partition on fine and coarse grid
    integer (kind = myint), dimension (:), pointer :: fwhere, cwhere

    ! cgraph: the coarse grid graph of row connectivity
    ! cmatrix: the coarse grid matrix
    type (ZD11_type), pointer :: cgraph, cmatrix

    ! matrix: the fine grid matrix
    type (ZD11_type), pointer ::  matrix

    ! row_wgt, crow_wgt: fine and coarse matrix row weight
    ! col_wgt: fine grid column weight
    integer (kind = myint), dimension (:), pointer :: row_wgt, &
         crow_wgt
    integer (kind = myint), dimension (:), pointer :: col_wgt

    ! cnvtx: number of vertex (rows) in the coarse matrix
    integer (kind = myint) :: cnvtx

    ! i,j: the loop index
    integer (kind = myint) :: i,j

    ! order: a random permutation
    integer (kind = myint), allocatable, dimension (:) :: order

    ! where: the best partition
    ! mincut: the best netcut so far
    integer (kind = myint), allocatable, dimension (:) :: where
    integer (kind = myint) :: mincut

    ! row_sum: sum of row weights of rows assigned to the partition 0
    ! for initial partition
    ! imin: the row with the smallest row weight
    ! min_rowwgt: minimum row weights
    ! rowwgt: a row weight
    integer (kind = myint) :: row_sum,imin,min_rowwgt,rowwgt

    ! agg: the aggressiveness
    real (kind = myreal) :: agg


    ! nnodes: number of nodes in subdomains 0 and 1
    ! nnodes_copy: back up of nnodes
    integer (kind = myint) :: nnodes(0:1),nnodes_copy(0:1)

    ! ierr: error tag
    integer (kind = myint) :: ierr

    ! ostream: channel for diagonalstic print
    ! estream: error stream
    ! wstream: warning stream
    ! plevel: print level for diagnostic print
    integer (kind = myint) :: estream,wstream,ostream, plevel

    info = 0

    ! Allocate local arrays
    allocate(order(grid%m),where(grid%m),stat = ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_ALLOC
       return
    end if

    ostream = controller%mp; wstream = controller%wp
    estream = controller%lp
    plevel = controller%print_level

    if (plevel >= 1) &
         call level_print(ostream,"size of grid on level ",&
         grid%level,&
         " is ",real(grid%m,myreal))
    if (plevel > 0) &
         write(ostream,"('matrix size=',3i12)") grid%m,grid%matrix%n,&
         sum(grid%col_wgt(1:grid%matrix%n))

    ! if this is the last level or matrix size too small,
    ! partition and return
    if (grid%level >= mglevel .or. &
         grid%m <= controller%coarsest_size) then
       nullify(grid%coarse)
       ! partitioning randomly subject to load balance,
       ! but only do this if this is not the finest level,
       ! on finest level the partitioning is
       ! input from outside
       if (grid%level /= 1) then
          mincut = huge(0)
          ! on the coarest grid, do NUM_COARSEST_KL random partition
          ! and pick the best result.
          do j = 1, controller%NUM_COARSEST_KL
             do i = 1, grid%m
                order(i) = i
             end do
             call random_permute(order,seed)
             ! assign rows to partition 0 until target_nnodes(0
             min_rowwgt = huge(0)
             row_sum = 0
             grid%where(order(1:grid%m)) = 1
             ! add nodes into domain 0 until full
             do i = 1, grid%m
                rowwgt = grid%row_wgt(order(i))
                if (rowwgt < min_rowwgt) then
                   min_rowwgt = rowwgt
                   imin = i
                end if
                if (row_sum + rowwgt <= &
                     target_nnodes(0)) then
                   row_sum = row_sum + rowwgt
                   grid%where(order(i)) = 0
                else
                   exit
                end if
             end do
             agg = 1.0 ! aggressive coarsening on the coarest level
             call KL(PLEVEL = plevel,OSTREAM = ostream,&
                  ESTREAM = estream, WSTREAM = wstream,&
                  MATRIX = grid%matrix,GRAPH = grid%graph, &
                  ROW_WGT = grid%row_wgt, COL_WGT = grid%col_wgt,&
                  WHERE = grid%where, target_nnodes = target_nnodes, &
                  AGGRESSIVE = agg, &
                  UNBALANCE =controller%max_imbalance, &
                  MINCUT = grid%mincut, NNODES = nnodes,INFO = ierr)
             if (ierr /= 0) then
                info = ierr
                if (info < 0) return
             end if

             if (mincut > grid%mincut) then
                nnodes_copy = nnodes
                mincut = grid%mincut
                where = grid%where
             end if
             if (plevel > 1) write(ostream,&
                  "(a,i3,a,i10,a,2i12,a,i12)") &
                  "after ini. part. on ",grid%level,&
                  " mincut ",grid%mincut," nnodes = ",nnodes,&
                  " max_row_wgt = ",maxval(grid%row_wgt)
          end do
          nnodes = nnodes_copy
          grid%where = where
          grid%mincut = mincut
          if (plevel >= 1) call level_print(ostream,&
               "on coarsest level ",grid%level,&
               " mincut =  ",real(mincut,myreal))
       else
          ! one level algorithm
          row_sum = 0
          do i = 1, grid%m
             order(i) = i
          end do
!          call random_permute(order)
          grid%where(order(1:grid%m)) = 1
          ! assign rows to partition 0 until target_nnodes(0)
          ! assign rows to partition 0 until target_nnodes(0
          min_rowwgt = huge(0)
          row_sum = 0
          grid%where(order(1:grid%m)) = 1
          ! add nodes into domain 0 until full
          do i = 1, grid%m
             rowwgt = grid%row_wgt(order(i))
             if (rowwgt < min_rowwgt) then
                min_rowwgt = rowwgt
                imin = i
             end if
             if (row_sum + rowwgt <= &
                  target_nnodes(0)) then
                row_sum = row_sum + rowwgt
                grid%where(order(i)) = 0
             else
                exit
             end if
          end do
          ! make sure domain 0 has at least one row
          if (row_sum == 0) then
             grid%where(order(imin)) = 0
          end if
          call KL(PLEVEL = plevel,OSTREAM = ostream,&
               ESTREAM = estream, WSTREAM = wstream,&
               MATRIX = grid%matrix,GRAPH = grid%graph, &
               ROW_WGT = grid%row_wgt, COL_WGT = grid%col_wgt,&
               WHERE = grid%where, target_nnodes = target_nnodes, &
               AGGRESSIVE = controller%kl_aggressive, &
               UNBALANCE = controller%max_imbalance, &
               MINCUT = grid%mincut, NNODES = nnodes,INFO = ierr)
          if (ierr /= 0) then
             info = ierr
             if (info < 0) return
          end if

       end if
       if (plevel > 1) write(ostream,"(6i12)") grid%m,grid%matrix%n,&
            sum(grid%col_wgt(1:grid%matrix%n)),nnodes,grid%mincut


       return
    end if

    ! coarsening
    if (plevel > 1) call level_print(ostream,"before coarsening",&
         grid%level)
    call coarsen(grid,controller,ierr)
    if (ierr /= 0) then
       info = ierr
       if (info < 0) return
    end if
    if (plevel > 1) call level_print(ostream,"after coarsening",&
         grid%level)
    cgrid => grid%coarse

    ! see if the grid reduction is not achieved,
    ! or the size of grid is <= 2, if so, set the allowed
    ! maximum level to current level and partition this level
    ! deallocate the coarse grid quantities that has been
    ! allocated so far
    if (cgrid%m/real(grid%m,myreal) >  controller%grid_rdc_fac.or.&
         cgrid%m <= 2) then
       if (plevel >= 1) then
          write(ostream,"(a,i10)") &
            " next coarse grid size = ",cgrid%m
          write(ostream,"(a,g12.4,a)") &
               " next coarse grid size/current grid size >  &
               &controller%grid_rdc_fac (",&
               controller%grid_rdc_fac,")"
          write(ostream,&
               "(' current level is taken as the coarest level')")
       end if
       call multilevel(grid,controller,grid%level,seed,target_nnodes,ierr)
       if (ierr /= 0) then
          info = ierr
          if (info < 0) return
       end if
       call multilevel_deallocate_last(cgrid,ierr)
       if (ierr /= 0) then
          info = ierr
          call monet_print_message(info,estream,wstream,&
               "inside multilevel_deallocate_last")
          if (info < 0) return
       end if
       return
    end if

    ! restriction ================

    ! form the coarse grid graph and matrix
    ! cmatrix = P^T*matrix = R*matrix
    p => cgrid%p
    matrix => grid%matrix
    cmatrix => cgrid%matrix
    cnvtx = cgrid%m

    ! r = P^T
    call mc65_matrix_transpose(p,r,ierr)
    if (ierr /= 0) then
       if (ierr == -1) info = ERR_MEMORY_ALLOC
       if (ierr == -2) info = ERR_MEMORY_DEALLOC
       call monet_print_message(info,estream,wstream,&
            "multilevel mc65_matrix_transpose")
       return
    end if

    ! cmatrix=r*matrix
    call mc65_matrix_multiply(r,matrix,cmatrix,ierr)
    if (ierr /= 0) then
       if (ierr == -1) info = ERR_MEMORY_ALLOC
       if (ierr == -2) info = ERR_MEMORY_DEALLOC
       call monet_print_message(info,estream,wstream,&
            "multilevel mc65_matrix_multiply")
       return
    end if


    ! alloate coarse grid quantities
    allocate(cgrid%where(cnvtx),cgrid%row_wgt(cnvtx),&
         cgrid%col_wgt(matrix%n), stat = ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_ALLOC
       call monet_print_message(info,estream,wstream,"multilevel")
       return
    end if

    ! condense the cmatrix, get rid of columns of only one entry
    ! and also merge columns of the same pattern
    col_wgt => cgrid%col_wgt
    ! initial column weight inherated from fine grid
    col_wgt(1:matrix%n) = grid%col_wgt(1:matrix%n)
    call mc65_matrix_condense(cmatrix,ierr,col_wgt_int=col_wgt,&
         iremove=1)
    if (ierr /= 0) then
       if (ierr == -1) info = ERR_MEMORY_ALLOC
       if (ierr == -2) info = ERR_MEMORY_DEALLOC
       call monet_print_message(info,estream,wstream,&
            "multilevel mc65_matrix_condense")
       return
    end if


    ! row weight cw = R*w
    row_wgt => grid%row_wgt
    crow_wgt => cgrid%row_wgt
    call mc65_matrix_multiply_vector(r,row_wgt,crow_wgt,ierr)

    ! deallocate the restrction matrix
    call mc65_matrix_destruct(r,ierr)

    ! extract the coarse grid by
    !      cgraph = cmatrix^T*cmatrix
    ! for the moment, inctead of the Gerlikin fomulation:
    !      cgraph = R*graph*P
    cgraph => cgrid%graph

    call matrix_to_graph(estream,wstream,cmatrix,cgraph,ierr,col_wgt)
    if (ierr /= 0) then
       info = ierr
       if (info < 0) return
    end if

    call multilevel(cgrid,controller,mglevel,seed,target_nnodes,ierr)
    if (ierr /= 0) then
       info = ierr
       if (info < 0) return
    end if

    ! prolongation ================

    ! injection of the partition on coarse grid to the fine grid:
    !    grid%where = P*cgrid%where
    ! here P is a special matrix with only one non-zero entry per row
    fwhere => grid%where
    cwhere => cgrid%where
    call mc65_matrix_multiply_vector(p,cwhere,fwhere,ierr)

    ! post smoothing
    call KL(PLEVEL = plevel,OSTREAM = ostream,&
         ESTREAM = estream, WSTREAM = wstream,&
         MATRIX = grid%matrix,GRAPH = grid%graph, &
         ROW_WGT = grid%row_wgt, COL_WGT = grid%col_wgt,&
         WHERE = grid%where, target_nnodes = target_nnodes, &
         AGGRESSIVE = controller%kl_aggressive,&
         UNBALANCE = controller%max_imbalance, &
         MINCUT = grid%mincut, NNODES = nnodes,INFO = ierr)
    if (ierr /= 0) then
       info = ierr
       if (info < 0) return
    end if


    if (plevel > 1) write(ostream,"(6i12)") grid%m,grid%matrix%n,&
         sum(grid%col_wgt(1:grid%matrix%n)),nnodes,grid%mincut


    if (plevel >= 1) then
       call level_print(ostream," after post smoothing ",&
            grid%level," mincut ",&
            real(grid%mincut,kind = myreal))
       call level_print(ostream," after post smoothing ",&
            grid%level," imbalance = ",&
            abs(nnodes(0)-nnodes(1))*real(100.0,myreal)/&
            (nnodes(0)+nnodes(1)))
    end if

    ! deallocate the previous level
    call multilevel_deallocate(cgrid,ierr)
    if (ierr /= 0) then
       info = ierr
       call monet_print_message(info,estream,wstream,&
            "multilevel_deallocate")
       if (info < 0) return
    end if

  end subroutine multilevel

  subroutine multilevel_deallocate(grid,info)
    type (multigrid), pointer :: grid
    INTEGER (KIND = myint), intent (out) :: INFO
    integer (kind = myint) :: ierr

    info = 0
    call mc65_matrix_destruct(grid%graph,ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_DEALLOC
       return
    end if
    call mc65_matrix_destruct(grid%matrix,ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_DEALLOC
       return
    end if
    call mc65_matrix_destruct(grid%p,ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_DEALLOC
       return
    end if
    deallocate(grid%graph, grid%p, &
         grid%matrix, grid%where, grid%mincut, grid%row_wgt, &
         grid%col_wgt, stat = ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_DEALLOC
       return
    end if
    deallocate(grid, stat = ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_DEALLOC
       return
    end if

  end subroutine multilevel_deallocate


  subroutine multilevel_deallocate_last(grid,info)
    INTEGER (KIND = myint), intent (out) :: INFO
    type (multigrid), pointer :: grid
    integer (kind = myint) :: ierr

    info = 0
    call mc65_matrix_destruct(grid%p,ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_DEALLOC
       return
    end if

    deallocate(grid%graph, grid%p, &
         grid%matrix, grid%mincut, stat = ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_DEALLOC
       return
    end if

    deallocate(grid,stat = ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_DEALLOC
       return
    end if
  end subroutine multilevel_deallocate_last


  recursive subroutine recursive_order(blklist,controller,seed,&
       num_bisect,order_global,excluded_col,info)

    ! blklist: the current list of blocks
    type (block_list), intent (inout) :: blklist

    ! controller: the controller containing control parameters
    type (mc66_control), intent (in) :: controller

    ! random number seed
    type (fa14_seed), intent (inout) :: seed

    ! this level of matrix (grid)
    type (multigrid), target :: grid

    integer (kind = myint), intent (inout) :: num_bisect
    integer (kind = myint), intent (inout) :: order_global(:)
    type (extendable_row), intent (inout) :: excluded_col
    ! info: information tag
    integer (kind = myint), intent (out) :: info

    ! the matrix to be partitioned
    type (ZD11_type), pointer :: matrix
    ! the graph correspnds to the row connectivity of matrix
    type (ZD11_type), target :: graph
    ! which row is on which processor (0 or 1)
    integer (kind = myint), pointer, dimension(:) :: row_proc

    ! super-column weight for the two new submatrices
    integer (kind = myint), dimension (:), pointer :: col_wgt

    type (ZD11_type),pointer :: matrix_up,matrix_down
    integer (kind = myint) :: num_repeated_column,ierr
    integer (kind = myint), pointer :: order(:),order_temp(:)
    integer (kind = myint) :: i

    ! the net cut for this bisection, and the
    ! total mincut
    integer (kind = myint), target :: mincut
    integer (kind = myint) :: mincut_all = 0

    ! imax: the matrix with the largest size/number of sub. domains
    !  in the current list
    integer (kind = myint) :: imax,maxsize,maxm


    ! target_nnodes: target sum of weights (number of nodes)
    ! in each of the two submatrices from the partition
    integer (kind = myint) :: target_nnodes(0:1)

    ! npart: number of parts this current matrix is supposed
    ! to be partitioned into
    ! avg: avarage size of the submatrix
    ! rem: reminder
    ! m: size of current matrix
    ! n: column dim. of current matrix
    ! npart0: number of parts the top half of the matrix will contain
    ! npart1:number of parts the bottom half of the matrix
    !        will contain
    ! col_size: column size of the original (finest) matrix
    integer (kind = myint) :: npart,avg,rem,m,n,npart0,npart1
    integer (kind = myint) :: col_size = -1


    ! ostream: channel for diagonalstic print
    ! estream: error stream
    ! wstream: warning stream
    ! plevel: print level for diagnostic print
    integer (kind = myint) :: estream,wstream,ostream, plevel

    info = 0
    ostream = controller%mp; wstream = controller%wp
    estream = controller%lp
    plevel = controller%print_level

    ! check through the matrix list to get the one with the
    ! largest number of submatrices to be partitioned into,
    ! and bisect that one
    maxsize = -huge(0_myint); maxm = -huge(0_myint)
    do i = 1, blklist%num_blocks
       if (blklist%blocks(i)%npart > maxsize .or.&
            (blklist%blocks(i)%npart == maxsize.and.&
            blklist%blocks(i)%mat%m > maxm)) then
          maxsize = blklist%blocks(i)%npart
          maxm = blklist%blocks(i)%mat%m
          imax = i
       end if
    end do

    ! bisect this matrix
    matrix => blklist%blocks(imax)%mat
    m = matrix%m
    if (col_size < 0) col_size = m
    if (blklist%num_blocks >= blklist%max_blocks.or.m < 2) then
       if (blklist%num_blocks /= blklist%max_blocks) then
          info = WARN_MATRIX_TOO_SMALL
          call monet_print_message(info,estream,wstream,&
               "in recursive ordering:")
       end if
       return
    end if

    ! work out the size of each of the twp submatrices
    ! and also the number of submatrices each will in turn
    ! be partitioned into
    npart = blklist%blocks(imax)%npart
    avg = m/npart; rem = m-avg*npart
    npart0 = (npart+1)/2; npart1=npart-npart0
    target_nnodes(0) = avg*npart0 + min(rem,npart0)
    target_nnodes(1) = m - target_nnodes(0)

    ! if one of the target nnodes is zero,
    ! we have too many blocks compared with matrix dimension.
    if (min(target_nnodes(0), target_nnodes(1)) <= 0) then
       info = WARN_MATRIX_TOO_SMALL
       call monet_print_message(info,estream,wstream,&
            "in recursive ordering:")
       return
    end if

    if (plevel > 2) then
       write(ostream,*) "target = ",target_nnodes," npart = ",npart
    end if

    ! First condense it
    col_wgt => blklist%blocks(imax)%col_wgt
    call mc65_matrix_condense(matrix,ierr,col_wgt_int=col_wgt,&
         iremove=1)
    if (ierr /= 0) then
       if (ierr == -1) info = ERR_MEMORY_ALLOC
       if (ierr == -2) info = ERR_MEMORY_DEALLOC
       call MONET_print_message(info,estream,wstream,&
            "recursive order mc65_matrix_condense")
       return
    end if
    m = matrix%m
    n = matrix%n


    allocate(order(m), stat = ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_ALLOC
       call monet_print_message(info,estream,wstream,&
            "recursive ordering")
       return
    end if

    num_bisect = num_bisect + 1
    blklist%num_blocks = blklist%num_blocks + 1

    ! work out the graph which describe the connectivity of rows

    call matrix_to_graph(estream,wstream,matrix,graph,ierr)
    if (ierr /= 0) then
       info = ierr
       if (info < 0) return
       return
    end if

    ! smoothing
    nullify(grid%fine)
    grid%m = graph%m
    grid%level = 1
    grid%graph => graph
    grid%matrix => matrix
    allocate(row_proc(graph%m),stat = ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_ALLOC
       call monet_print_message(info,estream,wstream,&
            "recursive ordering")
       return
    end if

    ! will be assigned inside multilebel
    grid%where => row_proc
    ! initialise row weight to 1
    allocate(grid%row_wgt(graph%m),grid%col_wgt(n),stat = ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_ALLOC
       call monet_print_message(info,estream,wstream,&
            "recursive ordering")
       return
    end if
    grid%row_wgt = 1
    grid%col_wgt(1:n) = col_wgt(1:n)
    grid%mincut => mincut
    if (plevel >= 0) write(ostream,"(a,i10,a,i10)") &
         "partition matrix of size ",m," X ",n
    call multilevel(grid,controller,controller%mglevel,seed,&
         target_nnodes,ierr)
    if (ierr /= 0) then
       info = ierr
       if (info < 0) return
    end if

    deallocate(grid%row_wgt,grid%col_wgt,stat = ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_DEALLOC
       call monet_print_message(info,estream,wstream,&
            "recursive ordering")
       return
    end if

    mincut_all = mincut_all + mincut
    if (plevel > 3) write(ostream,*) "after KL mincut = ",mincut
    if (plevel > 3) write(ostream,*) "after KL mincut_all = ",&
         mincut_all

    ! deallocate the graph concerning the connectivity of rows
    call mc65_matrix_destruct(graph,ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_DEALLOC
       call monet_print_message(info,estream,wstream,&
            "recursive ordering")
       return
    end if


    ! calculate the net cut and get rid of the repeated ones
    ! generate two matrices
    allocate(matrix_up,matrix_down,stat = ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_ALLOC
       call monet_print_message(info,estream,wstream,&
            "recursive ordering")
       return
    end if
    call repeated_column(controller,matrix,row_proc,&
         matrix_up,matrix_down,num_repeated_column, &
         order,excluded_col,col_wgt,ierr)
    if (ierr /= 0) then
       info = ierr
       if (info < 0) return
    end if
    deallocate(row_proc,stat = ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_DEALLOC
       call monet_print_message(info,estream,wstream,&
            "recursive ordering")
       return
    end if

    ! merge the local order into the global ordering
    allocate(order_temp(size(order)), stat = ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_ALLOC
       call monet_print_message(info,estream,wstream,&
            "recursive ordering")
       return
    end if

    do i = 1,size(order)
       order_temp(i) = order_global(order(i)+blklist%row_sta(imax)-1)
    end do
    order_global(blklist%row_sta(imax):blklist%row_sta(imax+1) - 1) &
         = order_temp
    deallocate(order,order_temp,stat = ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_DEALLOC
       call monet_print_message(info,estream,wstream,&
            "recursive ordering")
       return
    end if


    ! update the matrix list by adding in the two new matrix and take
    ! away the original parent matrix
    ! first move the list from imax+1 forward a space,
    ! then deallocate the old matrix
    ! (except if it is the root matrix,
    ! which will be deallocated on order.f90)
    ! to give room to two new matrix
    do i = blklist%num_blocks,imax+2,-1
       blklist%blocks(i)%mat => blklist%blocks(i-1)%mat
       blklist%blocks(i)%col_wgt => blklist%blocks(i-1)%col_wgt
       blklist%blocks(i)%npart = blklist%blocks(i-1)%npart
    end do

    ! destruct the already partitioned matrix
    call mc65_matrix_destruct(blklist%blocks(imax)%mat,ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_DEALLOC
       call monet_print_message(info,estream,wstream,&
            "recursive ordering")
       return
    end if

    deallocate(blklist%blocks(imax)%mat, stat = ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_DEALLOC
       call monet_print_message(info,estream,wstream,&
            "recursive ordering")
       return
    end if


    ! add two new matrix to the block list.
    blklist%blocks(imax)%mat => matrix_up
    blklist%blocks(imax)%npart = npart0
    blklist%blocks(imax+1)%mat => matrix_down
    blklist%blocks(imax+1)%npart = npart1

    ! inherate the column weights from that of the parent.
    allocate(col_wgt(matrix_down%n), stat = ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_ALLOC
       call monet_print_message(info,estream,wstream,&
            "recursive ordering")
       return
    end if

    col_wgt(1:matrix_down%n) = &
         blklist%blocks(imax)%col_wgt(1:matrix_down%n)
    blklist%blocks(imax+1)%col_wgt => col_wgt

    ! update the size pointer
    do i = blklist%num_blocks+1,imax+2,-1
       blklist%row_sta(i) = blklist%row_sta(i-1)
    end do
    blklist%row_sta(imax+1) = blklist%row_sta(imax) + matrix_up%m

    if (plevel >= 0) then
       write(ostream,"('blocks = ',i6,' netcut = ',i10,&
         &' netcut% = ',f10.2,'%')")&
         blklist%num_blocks,num_repeated_column,&
         100*real(num_repeated_column,kind = myreal)/col_size
!       write(ostream,"('block size ',100i10)") &
!            (blklist%row_sta(i+1)- &
!            blklist%row_sta(i), i=1, blklist%num_blocks)
       write(ostream,"(a)") " "
    end if

    ! resursivelly recalaulate the ordering
    call recursive_order(blklist,controller,seed,num_bisect,order_global,&
         excluded_col,ierr)
    if (ierr /= 0) then
       info = ierr
       if (info < 0) return
    end if

  end subroutine recursive_order


  subroutine repeated_column(controller,matrix,row_proc,matrix_up,&
       matrix_down,num_repeated_column,order,excluded_col,&
       col_wgt,info)
    ! calculating repeated entries (that shares the same
    ! column) in the two subdomain
    ! controller: the controller containing control parameters
    type (mc66_control), intent (in) :: controller
    type (ZD11_type), intent (in), target :: matrix
    integer (kind = myint), intent (out) :: num_repeated_column
    type (ZD11_type), intent (out), target :: matrix_up,matrix_down
    integer (kind = myint), intent (out) :: order(:)
    type (extendable_row), intent (inout) :: excluded_col
    ! super-column weights
    integer (kind = myint), intent (in) :: col_wgt(:)
    ! info: info tag
    integer (kind = myint), intent (out) :: info

    integer (kind = myint), pointer, dimension(:) :: row_proc
    integer (kind = myint) :: i,j,k,ierr,n0,n1
    integer (kind = myint), dimension (:), pointer :: the_row,the_row2
    integer (kind = myint), pointer, dimension (:) :: &
         foot_print0, foot_print1
    integer  (kind = myint):: num_up_row, num_down_row, &
         num_column,num_row
    integer  (kind = myint) :: nz_up, nz_down
    integer (kind = myint), pointer, dimension (:) :: ia_up,ia_down
    ! estream: error stream
    ! wstream: warning stream
    integer (kind = myint) :: estream,wstream

    info = 0
    estream = controller%lp; wstream = controller%wp

    ! the new order will let all rows on proc 0 be first,
    ! then the rest be second
    num_row = 0
    do i = 1, matrix%m
       if (row_proc(i) == 0) then
          num_row = num_row + 1
          order(num_row) = i
       end if
    end do
    do i = 1, matrix%m
       if (row_proc(i) == 1) then
          num_row = num_row + 1
          order(num_row) = i
       end if
    end do

    allocate(foot_print0(matrix%n),foot_print1(matrix%n),stat = ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_ALLOC
       call monet_print_message(info,estream,&
            wstream,"repeated column")
       return
    end if


    foot_print0 = 0; foot_print1 = 0
    do i = 1, matrix%m
!       the_row => mc65_matrix_getrow(matrix,i)
       call mc65_matrix_getrow(matrix,i,the_row)
       if (row_proc(i) == 0) then
          foot_print0(the_row) = &
               foot_print0(the_row) + 1
       elseif (row_proc(i) == 1) then
          foot_print1(the_row) = &
               foot_print1(the_row) + 1
       end if
    end do

    ! find the number of repeated entry,
    ! and also record the column names
    num_repeated_column  = 0
    do i = 1, matrix%n
       if (foot_print0(i) /= 0 .and. foot_print1(i) /= 0) then
          num_repeated_column = num_repeated_column + col_wgt(i)
          call extrow(excluded_col,i,col_wgt(i),ierr)
          if (ierr /= 0) then
             info = ierr
             call monet_print_message(info,estream,&
                  wstream,"extrow")
             if (info < 0) return
          end if
       else
          foot_print0(i) = 0
       end if
    end do

    ! number of rows that are in the first processor
    ! and on the second processor
    n1 = sum(row_proc)
    n0 = matrix%m - n1

    ! FIND THE NUMBER OF ZEROS

    ! find the number of nonzeros etc for the two matrix
    nz_up = 0; nz_down = 0
    num_up_row = 0
    num_down_row = 0
    do i = 1, matrix%m
!       the_row => mc65_matrix_getrow(matrix,i)
       call mc65_matrix_getrow(matrix,i,the_row)
       num_column = 0
       do j = 1, size(the_row)
          k = the_row(j)
          ! only count columns that does not
          ! appear in both partitioning
          if (foot_print0(k) == 0) then
             num_column = num_column + 1
          end if
       end do
       if (row_proc(i) == 0) then
          num_up_row = num_up_row + 1
          nz_up = nz_up + num_column
       else
          num_down_row = num_down_row + 1
          nz_down = nz_down + num_column
       end if
    end do

    call mc65_matrix_construct(matrix_up,n0,nz_up,ierr,&
         matrix%n,"pattern")
    if (ierr /= 0) then
       info = ERR_MEMORY_ALLOC
       call monet_print_message(info,estream,&
            wstream,"in repeated column mc65_matrix_construct")
       return
    end if
    call mc65_matrix_construct(matrix_down,n1,nz_down,ierr,&
         matrix%n,"pattern")
    if (ierr /= 0) then
       info = ERR_MEMORY_ALLOC
       call monet_print_message(info,estream,&
            wstream,"in repeated column mc65_matrix_construct")
       return
    end if
    ia_up => matrix_up%ptr
    ia_down => matrix_down%ptr

    ! UPDATE THE ROW POINTER
    num_up_row = 0
    num_down_row = 0
    ia_up(1) = 1; ia_down(1) = 1
    do i = 1, matrix%m
!       the_row => mc65_matrix_getrow(matrix,i)
       call mc65_matrix_getrow(matrix,i,the_row)
       num_column = 0
       do j = 1, size(the_row)
          k = the_row(j)
          ! only count columns that does not appear
          ! in both partitioning
          if (foot_print0(k) == 0) then
             num_column = num_column + 1
          end if
       end do
       if (row_proc(i) == 0) then
          num_up_row = num_up_row + 1
          ia_up(num_up_row+1) = ia_up(num_up_row) + num_column
       else
          num_down_row = num_down_row + 1
          ia_down(num_down_row+1) = ia_down(num_down_row) + num_column
       end if
    end do


    ! get rid of the repeated entry
    num_up_row = 0
    num_down_row = 0
    do i = 1, matrix%m
!      the_row => mc65_matrix_getrow(matrix,i)
       call mc65_matrix_getrow(matrix,i,the_row)
       ! which new matrix this row belongs to: up or down?
       if (row_proc(i) == 0) then
          num_up_row = num_up_row + 1
!          the_row2 => MC65_matrix_getrow(matrix_up,num_up_row)
          call mc65_matrix_getrow(matrix_up,num_up_row,the_row2)
       else
          num_down_row = num_down_row + 1
!          the_row2 => MC65_matrix_getrow(matrix_down,num_down_row)
          call mc65_matrix_getrow(matrix_down,num_down_row,the_row2)
       end if
       ! number of column entries in this row of
       ! the up or down matrix
       num_column = 0
       do j = 1, size(the_row)
          k = the_row(j)
          ! only count columns that does not appear
          ! in both partitioning
          if (foot_print0(k) == 0) then
             num_column = num_column + 1
             the_row2(num_column) = k
          end if
       end do
    end do
    deallocate(foot_print0,foot_print1,STAT = IERR)
    if (ierr /= 0) then
       info = ERR_MEMORY_DEALLOC
       call monet_print_message(info,estream,&
            wstream,"in repeated column deallocate")
       return
    end if

  end subroutine repeated_column


  subroutine monet(m,n,nz,irn,jcn,nblocks,control,seed, &
       row_order,info,rowptr,column_order,colptr,netcut,&
       rowdiff,kblocks)

    ! ======================== INPUT =========================
    ! m: INTEGER scalar of INTENT (IN). On entry it must
    !    be set by the user to hold the row dimension of the matrix
    ! n: INTEGER scalar of INTENT (IN). On entry it must be
    !    set by the user to hold the column dimension of the matrix
    ! nz: INTEGER scalar of INTENT (IN). On entry it must be
    !    set by the user to hold the number of nonzeros
    integer (kind=myint), intent (in) ::  m,n,nz

    ! irn: is an INTEGER rank one array  of size NZ with INTENT (IN).
    !      It must be set by the user to hold the row
    !      indices of the matrix.
    ! jcn: is an INTEGER rank one array  of size NZ with INTENT (IN).
    !      It must be set by the user to hold the column indices of
    !      the matrix.
    integer (kind=myint), intent (in) :: irn(nz)
    integer (kind=myint), intent (in) :: jcn(nz)

    ! nblocks:  INTEGER scalar of INTENT (IN). On entry it must
    !           be set by the user to hold the number of blocks
    !           in the Bordered Blocked Diagonal form that the
    !           sparse matrix is to be reordered into
    integer (kind = myint), intent (in) :: nblocks

    ! control: which controls the multilevel algorithm
    type (mc66_control), intent (inout) :: control

    ! the random number seed
    type (fa14_seed), intent (inout) :: seed

    ! ========================== OUTPUT =======================
    ! row ordering:  is an INTEGER rank one array  of size N
    !                with INTENT (OUT). row_order(i) is
    !                the original row index of the i-row of
    !                the reordered matrix.
    integer (kind = myint), intent (out) :: row_order(m)

    ! info: is an INTEGER scalar of INTENT (OUT). Info tag
    integer (kind = myint), intent (out) :: info



    ! ====================== optional output =======================
    ! rowptr: is an OPTIONAL INTEGER array of rank one of
    ! size nblocks+1 with INTENT (OUT).
    ! On exit, rowptr(I) gives the starting row index for
    ! block I in the reordered bordered blocked diagonal form.
    integer (kind = myint), optional, intent (out), &
         dimension (nblocks+1) :: rowptr

    ! column ordering:  is an INTEGER rank one array
    !    of size N with INTENT (OUT). column_order(i) is the
    !    original column index of the i-column of
    !    the reordered matrix.
    integer (kind = myint), intent (out), optional, &
         dimension (n) :: column_order

    ! colptr: is an OPTIONAL INTEGER array of rank one of
    ! size nblocks+1 with INTENT (OUT).
    ! On exit, colptr(I) gives the starting column index for
    ! block I in the reordered bordered blocked diagonal form.
    integer (kind = myint), intent (out), optional, &
         dimension (nblocks+1) :: colptr

    ! netcut: is an OPTIONAL INTEGER  scalar of INTENT (OUT).
    ! On exit, it holds the netcut  (the number of columns in the
    ! border of the SBBD form).
    integer (kind = myint), intent (out), optional :: netcut

    ! rowdiff: is an OPTIONAL REAL (double precision REAL
    ! for HSL_MC66_DOUBLE) scalar of INTENT (OUT). On exit,
    ! it holds the amount of rowdiff in row dimension of the
    ! blocks, defined as the difference between the
    ! maximum row dimension among all blocks and the
    ! average row dimension (= M/NBLOCKS), devided by the
    ! average row dimension, and expressed in percentage term.
    ! Thus rowdiff = 10 means that there is a 10% imbalance.
    real (kind = myreal), intent (OUT), optional :: rowdiff

    ! kblocks: is an OPTIONAL INTEGER  scalar of INTENT (OUT).
    ! on exit, it holds the actual number of blocks
    ! with nonzero row dimension
    ! in the SBBD form. kblocks = nblocks unless
    ! info = WARN_MATRIX_TOO_SMALL, which happens
    ! when the user set nblocks to be
    ! very close to or greater than M. In this case
    ! kblocks < nblocks.
    integer (kind = myint), intent (OUT), optional :: kblocks

    ! ==================== local variables =======================
    integer (kind = myint) :: i,ierr,num_bisect,ncut
    real (kind = myreal) :: imb
    type (ZD11_type), pointer :: matrix
    type (ZD11_type), target :: matrix_tmp
    type (extendable_row)  :: excluded_col

    ! a list of blocks
    type (block_list) :: blklist


    ! ostream: channel for diagonalstic print
    ! estream: error stream
    ! wstream: warning stream
    ! plevel: print level for diagnostic print
    integer (kind = myint) :: estream,wstream,ostream,plevel

    ! max_bk_size maximum block size after partitioning
    integer (kind = myint) :: max_bk_size

    info = 0

    ! check control parameters
    if (control%max_imbalance < 0) control%max_imbalance = 0
    if (control%max_imbalance > 1) control%max_imbalance = 1

    if (control%mglevel <= 0) control%mglevel = 1

    if (control%coarsen_scheme > 2) control%coarsen_scheme = 2
    if (control%coarsen_scheme < 1) control%coarsen_scheme = 1

    if (control%coarsest_size < 2) control%coarsest_size = 2

    if (control%num_coarsest_kl < 1) control%num_coarsest_kl = 1



    if (control%print_level > 1) control%print_level = 1
    if (control%mp < 0) control%print_level  = -huge(0)

    ! initialise the matrix list
    ostream = control%mp; wstream = control%wp
    estream = control%lp; plevel = control%print_level

    if (m <= 0) then
       info = ERR_M_NONPOSITIVE
       call monet_print_message(info,estream,wstream,"mc66")
       return
    end if
    if (n <= 0) then
       info = ERR_N_NONPOSITIVE
       call monet_print_message(info,estream,wstream,"mc66")
       return
    end if
    if (nz < 0) then
       info = ERR_NZ_NEGATIVE
       call monet_print_message(info,estream,wstream,"mc66")
       return
    end if

    if (nblocks <= 0 .or. nblocks > min(m,n)) then
       info = ERR_NEG_NBLOCKS
       call monet_print_message(info,estream,wstream,"mc66")
       return
    end if

    call block_list_init(blklist,nblocks,ierr)
    if (ierr /= 0) then
       info = ierr
       call monet_print_message(info,estream,wstream,&
            "in block_list_init")
       if (ierr < 0) return
    end if

    ! convert to internal matrix format to the first matrix
    ! on block_list

    ! this storage is deallocated in blklist
    allocate(matrix,stat = ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_ALLOC
       call monet_print_message(info,estream,wstream,"in mc66")
       if (info < 0) return
    end if
    call mc65_matrix_construct(matrix,m,n,nz,irn,jcn,ierr,&
         checking = 1)
    if (ierr /= 0) then
       if (ierr < 0) then
          if (ierr == -1) then
             info = ERR_MEMORY_ALLOC
          else
             info = ERR_MEMORY_DEALLOC
          end if
          call monet_print_message(info,estream,wstream,&
               "in mc66a mc65_matrix_construct")
          return
       else
          if (ierr == 2) then
             info = WARN_RANGE_IRN
          else if (ierr == 3) then
             info = WARN_RANGE_JCN
          else
             info = WARN_DUP_ENTRY
          end if
          call monet_print_message(info,estream,wstream,&
               "in mc66a mc65_matrix_construct")
       end if
    end if


    ! the first block point to the original matrix
    blklist%blocks(1)%mat => matrix

    !to be partitioned into this many submatrics
    blklist%blocks(1)%npart = nblocks

    ! allocate and set the column weights for this first block
    allocate(blklist%blocks(1)%col_wgt(matrix%n), stat = ierr)
    if (ierr /= 0) then
       info = ERR_MEMORY_ALLOC
       call monet_print_message(info,estream,wstream,"in mc66")
       return
    end if

    ! uniform super-column weights on entry
    blklist%blocks(1)%col_wgt = 1

    ! initialise the starting position of the second block
    ! needed in case nblocks = 1 and no partitioning is done.
    blklist%row_sta(2) = 1+matrix%m

    ! initially no excluded columns
    !
    call extrow_construct(excluded_col,0,ierr)
    if (ierr /= 0) then
       info = ierr
       call monet_print_message(info,estream,wstream,&
            "in extrow_construct")
       if (info < 0) return
    end if

    !
    ! recursive row_order the matrix
    !
    num_bisect  = 0
    do i = 1, matrix%m
       row_order(i) = i
    end do

    call recursive_order(blklist,control,seed,num_bisect,row_order,&
         excluded_col,ierr)
    if (ierr /= 0) then
       info = ierr
       if (info < 0) return
    end if

    ! copy the row and column starting index
    if (present(rowptr)) then
       ! in case we have less blocks then we hope to achieve
       ! the rest of row pointers is the same as the last in
       ! existence
       rowptr(1:blklist%num_blocks+1) = &
            blklist%row_sta(1:blklist%num_blocks+1)
       do i = blklist%num_blocks+2,nblocks+1
          rowptr(i) = rowptr(blklist%num_blocks+1)
       end do
    end if

    if (present(KBLOCKS)) then
       KBLOCKS = blklist%num_blocks
    end if

    ncut = nint(0.1+sum(excluded_col%val(1:excluded_col%num_column)))
    max_bk_size = -huge(0)
    do i = 1, blklist%num_blocks
       max_bk_size = max(max_bk_size,blklist%row_sta(i+1)&
            -blklist%row_sta(i))
    end do
    imb = 100*(max_bk_size-real(m,myreal)/nblocks)&
         /(real(m,myreal)/nblocks)
    if (plevel > -1) then
       write(ostream,"('netcut = ',i10,' netcut% = ',f12.2)") &
            ncut,100*ncut/real(n,myreal)
        write(ostream,"(a,f10.2,'%')") "rowdiff = ", imb
    end if

    ! if need column order we need to regenerate the matrix
    ! since the original matrix has been condensed
    if (present(column_order).and.present(colptr)) then
       call mc65_matrix_construct(matrix_tmp,m,n,nz,irn,jcn,ierr,&
         checking = 1)
       if (ierr /= 0) then
          if (ierr < 0) then
             if (ierr == -1) then
                info = ERR_MEMORY_ALLOC
             else
                info = ERR_MEMORY_DEALLOC
             end if
             call monet_print_message(info,estream,wstream,&
                  "in mc66b mc65_matrix_construct")
             return
          else
             if (ierr == 2) then
                info = WARN_RANGE_IRN
             else if (ierr == 3) then
                info = WARN_RANGE_JCN
             else
                info = WARN_DUP_ENTRY
             end if
             call monet_print_message(info,estream,wstream,&
                  "in mc66")
          end if
       end if
       call get_column_order(matrix_tmp,blklist%num_blocks,nblocks,&
            row_order,blklist%row_sta,column_order,colptr,ierr)
       if (ierr /= 0) then
          info = ierr
          call monet_print_message(info,estream,wstream,&
               "in MC66: get_column_order")
          if (info < 0) return
       end if

       call mc65_matrix_destruct(matrix_tmp,ierr)
       if (ierr /= 0) then
          info = ERR_MEMORY_DEALLOC
          call monet_print_message(info,estream,wstream,&
               "in mc66 matrix_destruct")
          return
       end if

    end if

    call extrow_destruct(excluded_col,ierr)
    if (ierr /= 0) then
       info = ierr
       call monet_print_message(info,estream,wstream,&
            "extrow_destruct")
       if (info < 0) return
    end if

    call block_list_deallocate(blklist,ierr)
    if (ierr /= 0) then
       info = ierr
       call monet_print_message(info,estream,wstream,&
            "block_list_deallocate")
       if (ierr < 0) return
    end if

    if (present(netcut)) then
       netcut = ncut
    end if
    if (present(rowdiff)) then
       rowdiff = imb
    end if

  end subroutine monet


end module HSL_mc66_double



