module parallel_mod
   !
   !  This module contains information about the MPI partitioning and 
   !  its communicators.
   !  
   ! There are 3 spaces in MagIC
   !    - Grid Space (φ, θ, r)
   !      - φ:    contiguous, local
   !      - θ:    contiguous, distributed
   !      - r:    contiguous, distributed
   !    - LM-Space (l, m, r) radia Loop
   !      - l:    contiguous, local
   !      - m: discontiguous, distributed
   !      - r:    contiguous, distributed
   !    - ML-Space (r, m, l) MLLoop
   !      - m: discontiguous, distributed
   !      - l: "sort of" contiguous, distributed
   !      - r:    contiguous, local
   !  
   ! The available MPI ranks are split into a 2D cartesian grid. 
   ! Notice that the cartesian for (l,m,r) *hast* to be the same 
   ! as (φ,θ,r). This is because the Legendre Transform (which 
   ! converts from (φ,θ) -> (l,m)) is called inside of the radial
   ! loop for each r independently. It would be very inefficient 
   ! (and messy) to reorganize the r's in the middle of the radial loop.
   ! 
   ! The communicators are:
   ! 
   ! - comm_gs: cartesian grid for (φ, θ, r)
   !   - comm_r: only the r direction of comm_gs
   !   - comm_theta: only the θ direction of comm_gs
   !   
   ! - comm_lm: copy of comm_gs
   !   - comm_r: same as above
   !   - comm_m: copy of comm_theta
   ! 
   ! - comm_mlo: cartesian grid for (m, l, r) used in the MLloop
   !   - comm_mo: only the m direction of comm_mlo_space
   !   - comm_lo: only the l direction of comm_mlo_space
   ! 
   ! The duplications kinda make the code more readable, however.
   ! 
   !   Author: Rafael Lago (MPCDF) May 2018
   ! 

   use MPI
   use logic, only: l_save_out, lVerbose

   implicit none

   !   MPI_COMM_WORLD
   integer :: rank, n_ranks
   
   !   Intra-node Information
   integer :: comm_intra
   integer :: intra_rank, n_ranks_intra
   integer, allocatable :: rank2intra(:)
   integer, allocatable :: intra2rank(:)
   
   !   Grid Space
   integer ::    comm_gs
   integer ::    comm_r,    comm_theta
   integer :: n_ranks_r, n_ranks_theta
   integer ::   coord_r,   coord_theta
   integer, allocatable :: gs2rank(:,:), &
                           rank2theta(:), &
                           rank2r(:)
   
   !   LM-Space (radial Loop)
   !   Those will be just copies of variables above for theta
   integer ::    comm_lm
   integer ::    comm_m
   integer :: n_ranks_m
   integer ::   coord_m
   integer, allocatable :: lm2rank(:,:), &   ! Maybe rename "lm2rank" to "mr2rank"? Because it takes lm2rank(coord_m, coord_r) as arguments...
                           rank2m(:)
   
   !   ML-Space (ML Loop)
   integer ::   comm_mlo
   integer ::    comm_lo,    comm_mo
   integer :: n_ranks_lo, n_ranks_mo, n_ranks_mlo
   integer ::   coord_lo,   coord_mo, coord_mlo
   integer, allocatable :: mlo2rank(:,:), &
                           rank2mo(:), &
                           rank2lo(:)
   
   !   Others (might be deprecated)
   integer :: nThreads
   integer :: nLMBs_per_rank
   integer :: rank_with_l1m0
   integer :: chunksize
   integer :: ierr
   logical :: l_verb_paral=.false.
   
   !   Maps coordinates from one cartesian grid to another.
   !   Acronyms convention:
   !   gsp:  grid space as (φ,θ,r). 
   !         Partition in (θ,r), all φ are local
   !   lmr:  spherical-harmonic domain, with tuples ordered as (l,m,r); used in the radial loop
   !         Partition in (m,r), all l are local
   !   mlo:  spherical-harmonic domain, with tuples ordered as (m,l,r); used in the "lmloop"
   !         Partition in (m,l), all r are local; however, the (m,l) tuples are treated as 1D 
   !         instead of 2D.
   !   rnk: rank according to MPI_COMM_WORLD
   type, public :: cart_map
      integer, allocatable :: gsp2rnk(:,:)
      integer, allocatable :: rnk2gsp(:,:)  ! rnk2gsp(:,1) => θ, rnk2gsp(:,2) => r
      integer, allocatable :: rnk2lmr(:,:)  ! rnk2lmr(:,1) => m, rnk2lmr(:,2) => r
      integer, allocatable :: rnk2mlo(:)   
      integer, allocatable :: lmr2mlo(:,:)
      integer, allocatable :: lmr2rnk(:,:)
      integer, allocatable :: mlo2lmr(:,:)  ! mlo2lmr(:,1) => m, mlo2lmr(:,2) => r
   end type cart_map
   
   type(cart_map) :: cart

contains
   
   !------------------------------------------------------------------------------
   subroutine initialize_mpi_world
      !
      !-- Get number (name) of processor
      !
      integer :: world_group, intra_group, i, j
      integer, allocatable :: tmp(:)
      call mpi_comm_rank(MPI_COMM_WORLD,rank,    ierr)
      call mpi_comm_size(MPI_COMM_WORLD,n_ranks, ierr)

      nThreads  = 1
      chunksize = 16

      if (rank .ne. 0) l_save_out = .false.
      if (rank .ne. 0) lVerbose   = .false.
      
      !
      !-- Figures out which ranks are local within their nodes
      !
      call mpi_comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, comm_intra, ierr)
      call mpi_comm_rank(comm_intra, intra_rank, ierr)
      call mpi_comm_size(comm_intra, n_ranks_intra, ierr)
      allocate(intra2rank(0:n_ranks_intra-1))
      allocate(rank2intra(0:n_ranks-1))
      
      call mpi_comm_group(MPI_COMM_WORLD, world_group, ierr)
      call mpi_comm_group(comm_intra, intra_group, ierr)
      
      allocate(tmp(0:n_ranks-1))
      do i=0,n_ranks-1
         tmp(i)=i
      end do
      
      call mpi_group_translate_ranks(world_group, n_ranks, tmp, intra_group, &
                  rank2intra, ierr)
      call mpi_group_translate_ranks(intra_group, n_ranks_intra, tmp(0:n_ranks_intra-1),&
                  world_group, intra2rank, ierr)
                  
      deallocate(tmp)
      
   end subroutine initialize_mpi_world
   
   !------------------------------------------------------------------------------
   subroutine finalize_cart_map
      
      deallocate(cart%gsp2rnk)
      deallocate(cart%rnk2gsp)
      deallocate(cart%rnk2lmr)
      deallocate(cart%rnk2mlo)
      deallocate(cart%lmr2mlo)
      deallocate(cart%mlo2lmr)
   
   end subroutine finalize_cart_map
   
   !------------------------------------------------------------------------------
   subroutine initialize_mpi_decomposition
      !
      !   Author: Rafael Lago (MPCDF) May 2018
      ! 
      integer :: dims(2), coords(2), i, irank
      logical :: periods(2)
      
      if (rank == 0) write(*,*) ' '
      call optimize_decomposition_simple
      if (rank == 0) then
         write(*,*) '! MPI Domain Decomposition Info: '
         write(*,'(A,I0,A,I0)') ' ! Grid Space (θ,r): ', n_ranks_theta, "x", n_ranks_r
         write(*,'(A,I0,A,I0)') ' !   LM Space (m,r): ', n_ranks_m,  "x", n_ranks_r
         write(*,'(A,I0,A,I0)') ' !   ML Space (l,m): ', n_ranks_lo, "x", n_ranks_mo
         write(*,'(A,I0)')      ' !      Total Ranks: ', n_ranks
      end if
      call check_decomposition
      
      !-- Grid Space
      !
      ! All arrays will be allocated as (θ,r), meaning that neighbouring ranks 
      ! should have contiguous θ. For that end, I need to flip the dimensions 
      ! when creating communicators using MPI because it uses row major instead 
      ! of Fortran's column major.
      dims    = (/n_ranks_r, n_ranks_theta/)
      periods = (/.true., .false./)
            
      call MPI_Cart_Create(MPI_COMM_WORLD, 2, dims, periods, .true., comm_gs, ierr)
      call check_MPI_error(ierr)
      
      call MPI_Comm_Rank(comm_gs, irank, ierr) 
      call check_MPI_error(ierr)
      
      call MPI_Cart_Coords(comm_gs, irank, 2, coords, ierr)
      call check_MPI_error(ierr)
      
      call MPI_Comm_Split(comm_gs, coords(2), irank, comm_r, ierr)
      call check_MPI_error(ierr)
      call MPI_Comm_Rank(comm_r, coord_r, ierr) 
      call check_MPI_error(ierr)
      
      call MPI_Comm_Split(comm_gs, coords(1), irank, comm_theta, ierr)
      call MPI_Comm_Rank(comm_theta, coord_theta, ierr) 
      call check_MPI_error(ierr)
      
      ! old
      allocate(rank2theta(0:n_ranks-1))
      allocate(rank2r(0:n_ranks-1))
      allocate(gs2rank(0:n_ranks_theta-1,0:n_ranks_r-1))
      ! new
      allocate(cart%gsp2rnk(0:n_ranks_theta-1,0:n_ranks_r-1))
      allocate(cart%rnk2gsp(0:n_ranks-1, 2))
      
      do i=0,n_ranks-1
         call mpi_cart_coords(comm_gs, i, 2, coords, ierr)
         
         cart%gsp2rnk(coords(2),coords(1)) = i
         cart%rnk2gsp(i,1) = coords(2)  ! θ/m
         cart%rnk2gsp(i,2) = coords(1)  ! r
         
         rank2theta(i) = coords(2)
         rank2r(i)     = coords(1)
         gs2rank(coords(2),coords(1)) = i
      end do
      
      !-- LM Space (radial Loop)
      !   
      !   This is just a copy
      !   
      !   PS: is it a problem that I'm using a periodic communicator for 
      !       lm-space as well?
      comm_lm   = comm_gs
      comm_m    = comm_theta
      n_ranks_m = n_ranks_theta
      coord_m   = coord_theta
      
      ! old
      allocate(rank2m(0:n_ranks-1))
      allocate(lm2rank(0:n_ranks_m-1,0:n_ranks_r-1))
      ! new
      allocate(cart%rnk2lmr(0:n_ranks-1, 2))
      allocate(cart%lmr2rnk(0:n_ranks_m-1,0:n_ranks_r-1))
      cart%rnk2lmr = cart%rnk2gsp
      cart%lmr2rnk = cart%gsp2rnk
      
      lm2rank   = gs2rank
      rank2m    = rank2theta
      
      !-- ML Space (ML Loop)
      !   
      !   
      n_ranks_mlo = n_ranks
      n_ranks_mo  = n_ranks_m
      n_ranks_lo  = n_ranks_r
      coord_mo    = coord_m
      coord_lo    = coord_r
      coord_mlo   = rank
      
      comm_mlo = MPI_COMM_WORLD
      
      allocate(rank2mo(0:n_ranks-1))
      allocate(rank2lo(0:n_ranks-1))
      allocate(mlo2rank(0:n_ranks_mo-1,0:n_ranks_lo-1))
      
      rank2lo  = rank2r
      rank2mo  = rank2m
      mlo2rank = gs2rank
      
      ! new
      allocate(cart%rnk2mlo(0:n_ranks-1))
      allocate(cart%mlo2lmr(0:n_ranks_mlo-1, 2))
      allocate(cart%lmr2mlo(0:n_ranks_m-1,0:n_ranks_r-1))
      ! I know, it is unnecessary. But it helps with the reading of the code imho
      do i=0,n_ranks-1
         cart%rnk2mlo(i) = i
      end do
      cart%mlo2lmr = cart%rnk2gsp
      cart%lmr2mlo = cart%lmr2rnk
      
      if (l_verb_paral) call print_mpi_distribution
      
   end subroutine initialize_mpi_decomposition
   
   !------------------------------------------------------------------------------
   subroutine check_MPI_error(code)

      integer, intent(in) :: code
      character(len=MPI_MAX_ERROR_STRING) :: error_str
      integer :: ierr, strlen

      if (code /= MPI_SUCCESS) then
          call MPI_Error_string(code, error_str, strlen, ierr)
          write(*, '(A, A)') 'MPI error: ', trim(error_str)
          call MPI_Abort(MPI_COMM_WORLD, code, ierr)
      endif

   end subroutine check_MPI_error
   
   !------------------------------------------------------------------------------
   subroutine finalize_mpi_decomposition

      call MPI_Comm_Free(comm_gs,    ierr) 
      call MPI_Comm_Free(comm_r,     ierr) 
      call MPI_Comm_Free(comm_theta, ierr) 
      
!       call MPI_Comm_Free(comm_mlo, ierr)
!       call MPI_Comm_Free(comm_mo,  ierr)
!       call MPI_Comm_Free(comm_lo,  ierr)
      
      deallocate(rank2theta, rank2r, gs2rank)
      deallocate(rank2m, lm2rank)
      deallocate(rank2mo, rank2lo, mlo2rank)
      
      call finalize_cart_map
      
   end subroutine finalize_mpi_decomposition

   !------------------------------------------------------------------------------
   subroutine optimize_decomposition_simple
      !   
      !   This is a *very* simple function to split all n_ranks into a 2D 
      !   grid. In the future we might want to make this function much more 
      !   sophisticated.
      !
      !   Author: Rafael Lago (MPCDF) May 2018
      ! 
      real     :: log2
      integer  :: nPow2
      
      !-- If the grid dimensions were not given
      if (n_ranks_theta == 0 .and. n_ranks_r == 0) then
         if (rank==0) write(*,*) '! Automatic splitting of MPI ranks for Grid Space...'
         
         log2 = log(real(n_ranks)) / log(2.0) ! Gets log2(n_ranks)
         n_ranks_theta = 2**FLOOR(log2/2)
         n_ranks_r = n_ranks/n_ranks_theta
      else if (n_ranks_theta == 0) then
         n_ranks_theta = n_ranks/n_ranks_r
      else if (n_ranks_r == 0) then
         n_ranks_r = n_ranks/n_ranks_theta
      end if
      
      n_ranks_m = n_ranks_theta
         
      n_ranks_mo = n_ranks_r
      n_ranks_lo = n_ranks_theta
      
   end subroutine optimize_decomposition_simple
   
   !------------------------------------------------------------------------------
   subroutine check_decomposition

      if (n_ranks_r * n_ranks_theta .NE. n_ranks) then
        print *, '! Invalid parallelization partition! Make sure that n_ranks_r'//&
                 '* n_ranks_theta equals the number of available processes!'
        stop
      end if
      if (n_ranks_r * n_ranks_m .NE. n_ranks) then
        print *, '! Invalid parallelization partition! Make sure that n_ranks_r'//&
                 '* n_ranks_m equals the number of available processes!'
        stop
      end if
      if (n_ranks_lo * n_ranks_mo .NE. n_ranks) then
        print *, '! Invalid parallelization partition! Make sure that n_ranks_lo'//&
                 '* n_ranks_mo equals the number of available processes!'
        stop
      end if

   end subroutine check_decomposition
   
   !------------------------------------------------------------------------------
   subroutine print_mpi_distribution
   
      integer :: world_group, tmp_group, i, j
      integer, allocatable :: tmp_ranks(:), world_ranks(:)
      
      if (rank==0) then
         write(*,'(A)') '! MPI ranks distribution:'
         write(*,'(A)') '! ----------------------'
      end if
      
      call mpi_comm_group(MPI_COMM_WORLD, world_group, ierr)
      allocate(tmp_ranks(0:n_ranks-1))
      allocate(world_ranks(0:n_ranks-1))
      do i=0,n_ranks-1
         tmp_ranks(i)=i
      end do
      
      ! Prints intranode ranks
      do i=0,n_ranks-1
         if ((rank==i) .and. (intra_rank == 0)) then
            write(*,'(A,I0,A)') '! Intranode [@',rank,']'
            call print_array(intra2rank)
!             write(*,'(A,I0)',ADVANCE="NO") '! [',intra2rank(0)
!             do j=1,n_ranks_intra-1
!                write(*,'(A,I0)',ADVANCE="NO") ",",intra2rank(j)
!             end do
!             write(*,'(A,I0)') ']'
         end if
         call mpi_barrier(mpi_comm_world,ierr)
      end do
      
      
      ! Translate ranks from r communicator
      call mpi_comm_group(comm_r, tmp_group, ierr)
      call mpi_group_translate_ranks(tmp_group, n_ranks_r, &
              tmp_ranks(0:n_ranks_r-1), world_group, &
              world_ranks(0:n_ranks_r-1), ierr)
      
      ! Prints radial ranks
      do i=0,n_ranks-1
         if ((rank==i) .and. (coord_r == 0)) then
            write(*,'(A,I0,A)') '! Radial [@',rank,']'
            call print_array(world_ranks(0:n_ranks_r-1))
         end if
         call mpi_barrier(mpi_comm_world,ierr)
      end do
      
      
      ! Translate ranks from θ communicator
      call mpi_comm_group(comm_theta, tmp_group, ierr)
      call mpi_group_translate_ranks(tmp_group, n_ranks_theta, &
              tmp_ranks(0:n_ranks_theta-1), world_group, &
              world_ranks(0:n_ranks_theta-1), ierr)
      
      ! Prints theta ranks
      do i=0,n_ranks-1
         if ((rank==i) .and. (coord_theta == 0)) then
            write(*,'(A,I0,A)') '! Theta [@',rank,']'
            call print_array(world_ranks(0:n_ranks_theta-1))
         end if
         call mpi_barrier(mpi_comm_world,ierr)
      end do
      
      
      deallocate(tmp_ranks)
      deallocate(world_ranks)
   
   end subroutine print_mpi_distribution
   
   !------------------------------------------------------------------------------
   subroutine print_array(arr)
      integer, intent(in) :: arr(:)
      integer :: i
   
      write(*,'(A,I0)',ADVANCE="NO") '! [',arr(1)
      do i=2,size(arr)
         write(*,'(A,I0)',ADVANCE="NO") ",",arr(i)
      end do
      write(*,'(A,I0)') ']'
   
   end subroutine

!-------------------------------------------------------------------------------
end module parallel_mod
