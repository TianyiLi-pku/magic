! #include "perflib_preproc.cpp"
! !-- Holds the main abstract transposition object, as well as the default
! !   transposition object. Others can be implemented later on.
! !
! !   Author: Rafael Lago (MPCDF) October 2019
! !   
! module transposition
! 
! 
!    !-- Abstract transposition object
!    !
!    type, abstract, public :: transp_t
!       complex(cp), target :: container_from(:,:,:)
!       complex(cp), target :: container_to(:,:,:)
!       integer, allocatable :: rq(:)
!       integer :: nRq
!       integer :: nfields
!       contains
!         procedure :: wait  => transp_wait
!         procedure :: start => transp_start_abs
!    end type transp_t
!    
!    !-- ML to Radial transposition object using 
!    !   non-blocking send-receive
!    !
!    type, public, extends(transp_t) :: ml2r_sendrecv_t
!       integer, pointer :: sends(:), recvs(:)
!       integer, pointer :: s_type(:), r_type(:)
!       contains
!         final :: ml2r_sendrecv_finalize
!         procedure :: start => ml2r_sendrecv_start
!    end type ml2r_sendrecv_t
!    
!    abstract  interface
!       subroutine transp_start_abs(self)
!          import :: transp_t
!          class(transp_t), intent(inout) :: self
!       end subroutine transp_start_abs
!    end interface
! 
!    contains
!    
!    !-- Wait function
!    !   
!    !   This is very generic since it just calls MPI's waitall. 
!    !   It can, however, be overriden if something more advanced is 
!    !   needed
!    !   
!    subroutine transp_wait(self)
!       class(transp_t), intent(inout) :: self
!       integer :: ierr
!       call mpi_waitall(this%nEq,self%rq,MPI_STATUSES_IGNORE,ierr)
!    end subroutine transp_wait
!    
!    !-- Constructor
!    !
!    function new_ml2r_sendrecv(container_from,container_to, s_type, r_type) result(self)
!       complex(cp), target, intent(in) :: container_from(n_mlo_loc, n_r_max, *)
!       complex(cp), target, intent(in) :: container_to(n_lm_loc, nRstart:nRstop, *)
!       integer, target, intent(in) :: s_type(:)
!       integer, target, intent(in) :: r_type(:)
!       type(ml2r_sendrecv_t) :: self
!       integer :: nsends, nrecvs
!       
!       self%nfields = size(container_from,3)
!       self%container_from => container_from
!       self%container_to => container_to
!       self%s_type => s_type
!       self%r_type => r_type
!       
!       nsends = size(ml2r_dests)
!       nrecvs = size(ml2r_sources)
!       
!       allocate(self%rq(nsends+nrecvs))
!       self%sends => self%rq(1:nsends)
!       self%recvs => self%rq(nsends+1:nrecvs)
!       
!       self%rq = MPI_REQUEST_NULL
!    end function new_ml2r_sendrecv
!    
!    !-- Destructor
!    !
!    subroutine ml2r_sendrecv_finalize(self)
!       type(ml2r_sendrecv_t), intent(inout) :: self
!       integer :: i, ierr
!       do i=1,self%nRq
!          if (self%rq(i) /= MPI_REQUEST_NULL) call mpi_cancel(self%rq(i), ierr)
!       end do
!       
!       nullify(self%container_from)
!       nullify(self%container_to)
!       nullify(self%sends)
!       nullify(self%recvs)
!       deallocate(self%rq)
!    end subroutine ml2r_sendrecv_finalize
!    
!    !----------------------------------------------------------------------------
!    subroutine ml2r_sendrecv_start(self)
!       !   
!       !   This is supposed to be a general-purpose transposition. All the 
!       !   customization should happen during the creation of the types!
!       !
!       !   Author: Rafael Lago (MPCDF) January 2018
!       !   
!       !   
!       class(ml2r_sendrecv), intent(inout) :: self
!       
!       integer :: i, j, k, icoord_mlo, icoord_r, il_r, ierr
!       
!       !-- Starts the sends
!       do j=1,size(ml2r_dests)
!          icoord_mlo = ml2r_dests(j)
!          icoord_r = cart%mlo2lmr(icoord_mlo,2)
!          il_r = dist_r(icoord_r,1)
!          call mpi_isend(self%container_from(1,il_r,1), self%nfields, self%s_type(j), icoord_mlo, 1, comm_mlo, self%sends(j), ierr)
!       end do
!       
!       !-- Starts the receives
!       do j=1,size(ml2r_sources)
!          icoord_mlo = ml2r_sources(j)
!          call mpi_irecv(self%container_to, self%nfields, self%r_type(j), icoord_mlo, 1, comm_mlo, self%recvs(j), ierr)
!       end do
!       
!       !-- Copies data which is already local
!       do i=1,size(ml2r_loc_dspl,1)
!          k = ml2r_loc_dspl(i,1)
!          j = ml2r_loc_dspl(i,2)
!          self%container_to(j,nRstart:nRstop,1:self%nfields) = self%container_from(k,nRstart:nRstop,1:self%nfields)
!       end do
!    end subroutine ml2r_sendrecv_start
!    
! 
!    
!    
! end module transposition
