#include "perflib_preproc.cpp"
module communications

#ifdef WITH_MPI
   use mpimod
#endif
   use precision_mod
   use mem_alloc, only: memWrite, bytes_allocated
   use parallel_mod, only: coord_r, n_ranks_r, ierr, & 
                           comm_r, rank, comm_gs,   &
                           comm_m, n_ranks_theta
   use LMLoop_data, only: llm, ulm
   use geometry 
   use blocking, only: lo_map, lmStartB, lmStopB
   use logic, only: l_mag, l_conv, l_heat, l_chemical_conv, &
       &            l_mag_kin, l_TP_form, l_double_curl
   use useful, only: abortRun
   use LMmapping, only: map_dist_st, map_glbl_st, map_mlo
 
 
use blocking
 
 
   implicit none
 
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !    private uncomment this
   public 
 
   interface get_global_sum
      module procedure get_global_sum_cmplx_2d, get_global_sum_cmplx_1d, &
                       get_global_sum_real_2d
   end interface
   
   interface get_global_sum_dist
      module procedure get_global_sum_cmplx_2d_dist, get_global_sum_cmplx_1d_dist, &
                       get_global_sum_real_2d_dist
   end interface
 
   type, public :: lm2r_type
      integer, allocatable :: final_wait_array(:), s_request(:), r_request(:)
      complex(cp), pointer :: temp_Rloc(:,:,:), arr_Rloc(:,:,:)
      integer :: count
   end type lm2r_type

   type, public :: r2lm_type
      integer, allocatable :: final_wait_array(:), s_request(:), r_request(:)
      complex(cp), pointer :: temp_Rloc(:,:,:)
      integer :: count
   end type r2lm_type
 
   type, public :: gather_type
      integer, allocatable :: gather_mpi_type(:)
      integer :: dim2
   end type gather_type
   
 
   ! MPI datatypes for the redistribution of the d?dt arrays
   integer, save, allocatable :: s_transfer_type(:),s_transfer_type_nr_end(:)
   integer, save, allocatable :: r_transfer_type(:),r_transfer_type_nr_end(:)
   integer, save, allocatable :: s_transfer_type_cont(:,:)
   integer, save, allocatable :: s_transfer_type_nr_end_cont(:,:)
   integer, save, allocatable :: r_transfer_type_cont(:,:)
   integer, save, allocatable :: r_transfer_type_nr_end_cont(:,:)
   integer :: r_lm_gather_type, r_lm_gather_type_lm_end
   integer, allocatable :: s_request(:),r_request(:),final_wait_array(:)
   integer, allocatable :: array_of_statuses(:,:)
 
   public :: gather_from_lo_to_rank0,scatter_from_rank0_to_lo, &
   &         gather_all_from_lo_to_rank0
   public :: get_global_sum, finalize_communications,                      &
   &         r2lo_redist_wait,r2lo_redist_start,initialize_communications, &
   &         create_lm2r_type!,lo2r_redist,lm2r_redist
   public :: lo2r_redist_start,lo2r_redist_wait, create_r2lm_type,         &
   &         destroy_r2lm_type,destroy_lm2r_type
   public :: get_global_sum_dist, r2lo_redist_start_dist,                  &
   &         lo2r_redist_start_dist, lo2r_redist_wait_dist

   
 
   ! declaration of the types for the redistribution
   type(lm2r_type), public :: lo2r_flow, lo2r_s
   type(lm2r_type), public :: lo2r_field, lo2r_xi

   type(r2lm_type), public :: r2lo_flow, r2lo_s, r2lo_xi, r2lo_b
 
   type(gather_type), public :: gt_OC,gt_IC,gt_cheb
 
   complex(cp), allocatable :: temp_gather_lo(:)
   complex(cp), allocatable :: temp_r2lo(:,:)
   
   !-- ML to Radial transposition object
   !
   !   Author: Rafael Lago (MPCDF) January 2018
   !   
   type, public :: ml2r_transp
      integer, pointer :: sends(:), recvs(:)
      integer, pointer :: rq(:)
      contains
        final :: deallocate_ml2r_transp
        procedure :: wait  => transpose_ml2r_wait
        procedure :: start => transpose_ml2r_start
   end type ml2r_transp
   
   integer, allocatable :: ml2r_s_type(:), ml2r_r_type(:)
   integer, allocatable :: ml2r_dests(:),  ml2r_sources(:)
   integer, allocatable :: ml2r_loc_dspl(:,:)
   
   type(ml2r_transp), public :: ml2r_flow, ml2r_s
   type(ml2r_transp), public :: ml2r_field, ml2r_xi
   
   !-- Radial to ML transposition object
   !
   !   Author: Rafael Lago (MPCDF) October 2019
   !   
   type, public :: r2lm_transp
      integer, pointer :: sends(:), recvs(:)
      integer, pointer :: rq(:)
      contains
        final :: deallocate_r2lm_transp
        procedure :: wait  => transpose_r2lm_wait
        procedure :: start => transpose_r2lm_start
   end type r2lm_transp
   
   type(r2lm_transp), public :: r2lm_dflowdt, r2lm_dsdt, &
        &             r2lm_dxidt, r2lm_dbdt
   
   !-- Slice/gather interface
   !   
   !-- TODO: These are not writen for performance right now. Dunno if it will 
   !   ever be needed.
   !   
   interface slice_Flm
      module procedure slice_Flm_cmplx, slice_Flm_real
   end interface slice_Flm
   
   interface slice_FlmP
      module procedure slice_FlmP_cmplx, slice_FlmP_real
   end interface slice_FlmP
   
   !-- ?????????
   !
   public :: myAllGather, slice_f, slice_Flm, slice_FlmP, gather_f, &
             gather_FlmP, gather_Flm, transpose_m_theta, transpose_theta_m, &
             transform_new2old, transform_old2new, printMatrix, printTriplets,&
             printMatrixInt, printArray
             
contains

   !-- Constructor
   !   
   !   Initialize the MPI types for the transposition from ML to Radial
   !   
   !   Author: Rafael Lago (MPCDF) January 2018
   !   
   function allocate_ml2r_transp() result(self)
      type(ml2r_transp) :: self
      integer :: nsends, nrecvs
      nsends = size(ml2r_dests)
      nrecvs = size(ml2r_sources)
      
      allocate(self%rq(nsends+nrecvs))
      self%sends => self%rq(1:nsends)
      self%recvs => self%rq(nsends+1:nrecvs)
      
      self%rq = MPI_REQUEST_NULL
   end function allocate_ml2r_transp
   
   subroutine deallocate_ml2r_transp(self)
      type(ml2r_transp), intent(inout) :: self
      nullify(self%sends)
      nullify(self%recvs)
      nullify(self%rq)
   end subroutine deallocate_ml2r_transp
   
   !-- Constructor
   !   
   !   Initialize the MPI types for the transposition from Radial to LM
   !   
   !   Author: Rafael Lago (MPCDF) January 2018
   !   
   function allocate_r2lm_transp() result(self)
      type(r2lm_transp) :: self
      integer :: nsends, nrecvs
      nsends = size(ml2r_dests)
      nrecvs = size(ml2r_sources)
      
      allocate(self%rq(nsends+nrecvs))
      self%sends => self%rq(1:nsends)
      self%recvs => self%rq(nsends+1:nrecvs)
      
      self%rq = MPI_REQUEST_NULL
   end function allocate_r2lm_transp
   
   subroutine deallocate_r2lm_transp(self)
      !   
      !   Author: Rafael Lago (MPCDF) October 2019
      !    
      type(r2lm_transp), intent(inout) :: self
      nullify(self%sends)
      nullify(self%recvs)
      nullify(self%rq)
   end subroutine deallocate_r2lm_transp

   !----------------------------------------------------------------------------
   subroutine initialize_ml2r_tranposition
      !   
      !   Initialize the MPI types for the transposition from ML to Radial
      !   The types for the requests must be created somewhere else!
      !   This is supposed to be a general-purpose transposition. It might be 
      !   possible to optimize it for further specialized data structures.
      !   The MPI datatype here is built as indexed->resized->vector->resized
      !   
      !   Author: Rafael Lago (MPCDF) January 2019
      !    
      integer :: i, j, l, m, icoord_m, icoord_mlo, icoord_r, in_r, lm, ierr
      integer :: col_type, ext_type, mtx_type
      integer :: send_displacements(n_mlo_array, 0:n_ranks_mlo-1)
      integer :: recv_displacements(n_lm_loc, 0:n_ranks_mlo-1)
      integer :: send_counter_i(0:n_ranks_mlo-1), recv_counter_i(0:n_ranks_mlo-1), inblocks
      integer :: blocklenghts(max(n_mlo_array, n_lm_loc))
      
      integer (kind=mpi_address_kind) :: lb, extend, bytesCMPLX
      integer :: nsends, nrecvs
      
      call mpi_type_get_extent(MPI_DOUBLE_COMPLEX, lb, bytesCMPLX, ierr) 
      
      send_displacements = -1
      recv_displacements = -1
      send_counter_i     = 0
      recv_counter_i     = 0
      blocklenghts       = 1 ! This is actually fixed to 1!
      lb = 0
      
      !-- Loops over each (m,l) tuple to determine which rank in lmr will need it
      !   There is no good way of doing this; I could loop over the *local* tuples,
      !   but then I would get displacements in a funny order. The "easiest" 
      !   solution I see is to simply loop over *all* (m,l) tuples and skip those
      !   which do not belong here.
      do lm=1,lm_max
         l = map_glbl_st%lm2l(lm)
         m = map_glbl_st%lm2m(lm)
         i = map_mlo%ml2i(m,l)
         if (i<1) cycle
         
         icoord_m = m_tsid(m)
         
         do icoord_r=0, n_ranks_r-1
            icoord_mlo = cart%lmr2mlo(icoord_m, icoord_r)
            inblocks = send_counter_i(icoord_mlo) + 1
            send_counter_i(icoord_mlo) = inblocks
            send_displacements(inblocks,icoord_mlo) = i-1
         end do
      end do
      
      !-- Counts how many ranks will receive from me, and then allocate the (global) 
      !   array of types. 
      !   We subtract 1, if this rank has to "send data to itself", because we will
      !   not create a datatype nor an MPI request for that. 
      !   I cannot fathom a realitic situation in which a rank would *not* have to
      !   send something to itself, but who knows, maybe I'm missing something
      nsends = count(send_counter_i>0)
      if (send_counter_i(coord_mlo)>0) nsends = nsends - 1
      allocate(ml2r_s_type(nsends))
      allocate(ml2r_dests(nsends))
      
      !-- Build the Send MPI types:
      !   First we create indexed type; it is a vector which looks more or less like
      !      send_buf({i1,i2,i3,...}, j, k)
      !   Then we resize it, so that it has the length of the first dimension. Next
      !   we build an vector type out of it obtaining
      !      send_buf({i1,i2,i3,...}, l_r:u_r, k)
      !   Finally, we resize that, so that it comprises the whole data for a single 
      !   variable. Since we send the whole containers at once, we will be sending
      !   and array of those last types
      j = 0
      do icoord_mlo=0,n_ranks_mlo-1
         
         inblocks = send_counter_i(icoord_mlo)
         if (inblocks>0 .and. icoord_mlo /= coord_mlo) then
            j=j+1
            ml2r_dests(j) = icoord_mlo
            icoord_r = cart%mlo2lmr(icoord_mlo,2)
            in_r = dist_r(icoord_r,0)
            
            call MPI_Type_indexed(inblocks,blocklenghts(1:inblocks), send_displacements(1:inblocks,icoord_mlo), &
               MPI_DOUBLE_COMPLEX, col_type, ierr)
               
            extend = int(n_mlo_loc*bytesCMPLX,kind=mpi_address_kind)
            call MPI_Type_create_resized(col_type, lb, extend, ext_type, ierr)
            
            call MPI_Type_vector(in_r, 1, 1, ext_type, mtx_type, ierr)
            
            extend = int(n_mlo_loc*n_r_max*bytesCMPLX,kind=mpi_address_kind)
            call MPI_Type_create_resized(mtx_type, lb, extend, ml2r_s_type(j), ierr)
            
            call MPI_Type_commit(ml2r_s_type(j),ierr)
            call MPI_Type_free(col_type,ierr)
            call MPI_Type_free(mtx_type,ierr)
            call MPI_Type_free(ext_type,ierr)
         end if
      end do
      
      !-- Loops over each local (l,m) tuple to figure out which rank will send it to me!
      do i=1,n_lm_loc
         m  = map_dist_st%lm2m(i)
         l  = map_dist_st%lm2l(i)
         icoord_mlo = mlo_tsid(m,l)
         
         inblocks = recv_counter_i(icoord_mlo) + 1
         recv_counter_i(icoord_mlo) = inblocks
         recv_displacements(inblocks,icoord_mlo) = i-1
      end do
      
      !-- Counts how many ranks will send to me (similar to nsends)
      nrecvs = count(recv_counter_i>0)
      if (send_counter_i(coord_mlo)>0) nrecvs = nrecvs - 1
      allocate(ml2r_r_type(nrecvs))
      allocate(ml2r_sources(nrecvs))
      
      !-- Build the Recv MPI types (similar to Send MPI Types)
      !
      j = 0
      do icoord_mlo=0,n_ranks_mlo-1

         inblocks = recv_counter_i(icoord_mlo)
         if (inblocks>0 .and. icoord_mlo /= coord_mlo) then
            j = j + 1
            ml2r_sources(j) = icoord_mlo
            call MPI_Type_indexed(inblocks,blocklenghts(1:inblocks),&
               recv_displacements(1:inblocks,icoord_mlo), MPI_DOUBLE_COMPLEX, col_type, ierr)
            
            extend = int(n_lm_loc*bytesCMPLX,kind=mpi_address_kind)
            call MPI_Type_create_resized(col_type, lb, extend, ext_type, ierr)
            call MPI_Type_vector(n_r_loc, 1, 1, ext_type, mtx_type, ierr)
            
            extend = int(n_lm_loc*n_r_loc*bytesCMPLX,kind=mpi_address_kind)
            call MPI_Type_create_resized(mtx_type, lb, extend, ml2r_r_type(j), ierr)
            
            call MPI_Type_commit(ml2r_r_type(j),ierr)
            call MPI_Type_free(col_type,ierr)
            call MPI_Type_free(mtx_type,ierr)
            call MPI_Type_free(ext_type,ierr)
         end if
      end do
      
      !-- Finally, we need to keep the info pertaining to the copying of 
      !   data which is already local
      inblocks = send_counter_i(coord_mlo)
      allocate(ml2r_loc_dspl(inblocks,2))
      ml2r_loc_dspl(:,1) = send_displacements(1:inblocks,coord_mlo) + 1
      ml2r_loc_dspl(:,2) = recv_displacements(1:inblocks,coord_mlo) + 1
      
   end subroutine initialize_ml2r_tranposition
   
   !----------------------------------------------------------------------------
   subroutine finalize_ml2r_tranposition
      !
      !   Author: Rafael Lago (MPCDF) January 2018
      !   
      !       call mpi_type_free() #TODO
      deallocate(ml2r_s_type)
      deallocate(ml2r_dests)
      deallocate(ml2r_loc_dspl)
   end subroutine finalize_ml2r_tranposition
   
   !----------------------------------------------------------------------------
   subroutine finalize_communications_new
      !
      !   Author: Rafael Lago (MPCDF) January 2018
      !   
      call finalize_ml2r_tranposition
      call deallocate_ml2r_transp(ml2r_s)
      call deallocate_ml2r_transp(ml2r_flow)
      call deallocate_ml2r_transp(ml2r_field)
      call deallocate_ml2r_transp(ml2r_xi)
   end subroutine finalize_communications_new
  
   !----------------------------------------------------------------------------
   subroutine initialize_communications_new
   !-- Will replace the function below
   !
   !   Author: Rafael Lago (MPCDF) January 2018
   !   
      call initialize_ml2r_tranposition
      ml2r_s     = allocate_ml2r_transp()
      ml2r_flow  = allocate_ml2r_transp()
      ml2r_field = allocate_ml2r_transp()
      ml2r_xi    = allocate_ml2r_transp()
      
      r2lm_dflowdt = allocate_r2lm_transp()
      r2lm_dsdt    = allocate_r2lm_transp()
      r2lm_dxidt   = allocate_r2lm_transp()
      r2lm_dbdt    = allocate_r2lm_transp()
   end subroutine initialize_communications_new
  
   subroutine initialize_communications
      integer :: proc,my_lm_per_rank
      integer(lip) :: local_bytes_used
#ifdef WITH_MPI
      integer(kind=MPI_ADDRESS_KIND) :: zerolb, extent, sizeof_double_complex
!       integer(kind=MPI_ADDRESS_KIND) :: lb_marker, myextent, true_lb, true_extent
      integer :: base_col_type,temptype
      integer :: blocklengths(8),blocklengths_on_last(8),displs(8),displs_on_last(8)
      integer :: i, n_r_on_last_rank, n_r_per_rank

      ! first setup the datatype. It is not equal for all ranks. The n_ranks_r-1 coord_r can
      ! have a smaller datatype.
      ! Due to the different number of radial and lm points for the ranks, 
      ! we need essentially three different datatypes
      ! transfer_type: Standard for the transfers between ranks 0-(n_ranks_r-2)
      ! transfer_type_nr_end: for transfers involving coord_r (n_ranks_r-1) as receiver
      ! transfer_type_lm_end: for transfers involving coord_r (n_ranks_r-1) as sender
      ! +----+----+----+----+
      ! |    |    |    |    |
      ! |    |    |    |    |
      ! +----+----+----+----+
      ! |    |    |    |    |
      ! |    |    |    |    |
      ! +----+----+----+----+
      ! lm_per_rank is set here
      ! ATTENTION: for the last coord_r, the numbers are different and are
      !            stored in n_r_on_last_rank and lm_on_last_rank

      local_bytes_used = bytes_allocated
      allocate(s_transfer_type(n_ranks_r))
      allocate(s_transfer_type_nr_end(n_ranks_r))
      allocate(r_transfer_type(n_ranks_r))
      allocate(r_transfer_type_nr_end(n_ranks_r))
      allocate(s_transfer_type_cont(n_ranks_r,8))
      allocate(s_transfer_type_nr_end_cont(n_ranks_r,8))
      allocate(r_transfer_type_cont(n_ranks_r,8))
      allocate(r_transfer_type_nr_end_cont(n_ranks_r,8))
      bytes_allocated = bytes_allocated + 32*n_ranks_r*SIZEOF_INTEGER
      
      ! -- TODO:
      !    Those were for the old code, which had n_r_loc per rank fixed, 
      !    except for the last rank. Now things are much more flexible.
      !    Eliminate these two variables as soon as the parallelization of the 
      !    ML Loop is completed!
      ! 
      n_r_per_rank = dist_r(0,0)
      n_r_on_last_rank  = dist_r(n_ranks_r-1,0)

      do proc=0,n_ranks_r-1
         my_lm_per_rank=lmStopB(proc+1)-lmStartB(proc+1)+1
         !write(*,"(2(A,I4))") "lm_per_rank on coord_r ", proc," is ",my_lm_per_rank
         call MPI_Type_vector(n_r_per_rank,my_lm_per_rank,&
              &lm_max,MPI_DEF_COMPLEX,s_transfer_type(proc+1),ierr)
         call MPI_Type_commit(s_transfer_type(proc+1),ierr)
         if (ierr /= MPI_SUCCESS) print *,"~~~~~>",rank,__LINE__,__FILE__

         ! The same for the last coord_r for nR
         call MPI_Type_vector(n_r_on_last_rank,my_lm_per_rank,&
              &lm_max,MPI_DEF_COMPLEX,s_transfer_type_nr_end(proc+1),ierr)
         call MPI_Type_commit(s_transfer_type_nr_end(proc+1),ierr)
         if (ierr /= MPI_SUCCESS) print *,"~~~~~>",rank,__LINE__,__FILE__

         ! we do not need special receive datatypes, as the buffers are 
         ! contiguous in memory but for ease of reading, we define the 
         ! receive datatypes explicitly
         call MPI_Type_contiguous(my_lm_per_rank*n_r_per_rank,&
              & MPI_DEF_COMPLEX,r_transfer_type(proc+1),ierr)
         call MPI_Type_commit(r_transfer_type(proc+1),ierr)
         if (ierr /= MPI_SUCCESS) print *,"~~~~~>",rank,__LINE__,__FILE__
         call MPI_Type_contiguous(my_lm_per_rank*n_r_on_last_rank,&
              &MPI_DEF_COMPLEX,r_transfer_type_nr_end(proc+1),ierr)
         call MPI_Type_commit(r_transfer_type_nr_end(proc+1),ierr)
         if (ierr /= MPI_SUCCESS) print *,"~~~~~>",rank,__LINE__,__FILE__


         ! define the transfer types for the containers
         ! same schema as for the other types
         ! some temporary datatypes, not needed for communication
         ! but only for constructing the final datatypes
         call MPI_Type_get_extent(MPI_DEF_COMPLEX,zerolb,sizeof_double_complex,ierr)
         call MPI_Type_contiguous(my_lm_per_rank,MPI_DEF_COMPLEX,temptype,ierr)
         zerolb=0
         extent=lm_max*sizeof_double_complex
         call MPI_Type_create_resized(temptype,zerolb,extent,base_col_type,ierr)
         !call MPI_type_get_extent(base_col_type,lb_marker,myextent,ierr)
         !write(*,"(2(A,I10))") "base_col_type: lb = ",lb_marker,", extent = ",myextent
         blocklengths = [ n_r_per_rank, n_r_per_rank, n_r_per_rank, n_r_per_rank, &
                       &  n_r_per_rank, n_r_per_rank, n_r_per_rank, n_r_per_rank ]
         displs       = [ 0,           n_r_per_rank, 2*n_r_per_rank,    &
                       &  3*n_r_per_rank, 4*n_r_per_rank, 5*n_r_per_rank,&
                       &  6*n_r_per_rank, 7*n_r_per_rank ] 
         blocklengths_on_last = [ n_r_on_last_rank,n_r_on_last_rank,n_r_on_last_rank,&
                               &  n_r_on_last_rank,n_r_on_last_rank,n_r_on_last_rank,&
                               &  n_r_on_last_rank, n_r_on_last_rank ]
         displs_on_last     = [ 0,          n_r_on_last_rank, 2*n_r_on_last_rank,   &
                             &  3*n_r_on_last_rank, 4*n_r_on_last_rank,             &
                             &  5*n_r_on_last_rank, 6*n_r_on_last_rank,             &
                             &  7*n_r_on_last_rank ]
         do i=1,8
            call MPI_Type_vector(i,n_r_per_rank*my_lm_per_rank,n_r_max*my_lm_per_rank,&
                 & MPI_DEF_COMPLEX,r_transfer_type_cont(proc+1,i),ierr)
            call MPI_Type_commit(r_transfer_type_cont(proc+1,i),ierr)
            if (ierr /= MPI_SUCCESS) print *,"~~~~~>",rank,__LINE__,__FILE__
            call MPI_Type_vector(i,n_r_on_last_rank*my_lm_per_rank, &
                 & n_r_max*my_lm_per_rank,MPI_DEF_COMPLEX,      &
                 & r_transfer_type_nr_end_cont(proc+1,i),ierr)
            call MPI_Type_commit(r_transfer_type_nr_end_cont(proc+1,i),ierr)
            if (ierr /= MPI_SUCCESS) print *,"~~~~~>",rank,__LINE__,__FILE__

            call MPI_Type_indexed(i,blocklengths(1:i),&
                 & displs(1:i),base_col_type,s_transfer_type_cont(proc+1,i),ierr)
            call MPI_Type_commit(s_transfer_type_cont(proc+1,i),ierr)
            if (ierr /= MPI_SUCCESS) print *,"~~~~~>",rank,__LINE__,__FILE__
            call MPI_Type_indexed(i,blocklengths_on_last(1:i),&
                 & displs_on_last(1:i),base_col_type,         &
                 & s_transfer_type_nr_end_cont(proc+1,i),ierr)
            call MPI_Type_commit(s_transfer_type_nr_end_cont(proc+1,i),ierr)
            if (ierr /= MPI_SUCCESS) print *,"~~~~~>",rank,__LINE__,__FILE__

#if 0
            if (i == 8) then
               call MPI_type_get_extent(r_transfer_type_cont(proc+1,i), &
                                        lb_marker,myextent,ierr)
               call MPI_type_get_true_extent(r_transfer_type_cont(proc+1,i), &
                                             true_lb,true_extent,ierr)
               write(*,"(2(A,I3),3(A,I10))")                                    &
                    & "r_transfer_type_cont(",proc+1,",",i,"): lb = ",lb_marker,&
                    & ", extent = ",myextent,", true extent = ",true_extent
               call MPI_type_get_extent(s_transfer_type_cont(proc+1,i), &
                                        lb_marker,myextent,ierr)
               call MPI_type_get_true_extent(s_transfer_type_cont(proc+1,i), &
                                             true_lb,true_extent,ierr)
               write(*,"(2(A,I3),3(A,I10))")                                    &
                    & "s_transfer_type_cont(",proc+1,",",i,"): lb = ",lb_marker,&
                    & ", extent = ",myextent,", true extent = ",true_extent
               
               call MPI_type_get_extent(r_transfer_type_nr_end_cont(proc+1,i), &
                                        lb_marker,myextent,ierr)
               call MPI_type_get_true_extent(r_transfer_type_nr_end_cont(proc+1,i),&
                                             true_lb,true_extent,ierr)
               write(*,"(2(A,I3),3(A,I10))")                                           &
                    & "r_transfer_type_nr_end_cont(",proc+1,",",i,"): lb = ",lb_marker,&
                    & ", extent = ",myextent,", true extent = ",true_extent
               call MPI_type_get_extent(s_transfer_type_nr_end_cont(proc+1,i),&
                                        lb_marker,myextent,ierr)
               call MPI_type_get_true_extent(s_transfer_type_nr_end_cont(proc+1,i),&
                                             true_lb,true_extent,ierr)
               write(*,"(2(A,I3),3(A,I10))")                                           &
                    & "s_transfer_type_nr_end_cont(",proc+1,",",i,"): lb = ",lb_marker,&
                    & ", extent = ",myextent,", true extent = ",true_extent
            end if

#endif
         end do
      end do
#else
      local_bytes_used=bytes_allocated
#endif


      call create_gather_type(gt_OC,n_r_max)
      call create_gather_type(gt_IC,n_r_ic_max)

#ifdef WITH_MPI
      allocate(s_request(n_ranks_r-1),r_request(n_ranks_r-1))
      bytes_allocated = bytes_allocated + 2*(n_ranks_r-1)*SIZEOF_INTEGER
      allocate(array_of_statuses(MPI_STATUS_SIZE,2*(n_ranks_r-1)))
      bytes_allocated = bytes_allocated + &
                        2*(n_ranks_r-1)*MPI_STATUS_SIZE*SIZEOF_INTEGER
      allocate(final_wait_array(2*(n_ranks_r-1)))
      bytes_allocated = bytes_allocated + 2*(n_ranks_r-1)*SIZEOF_INTEGER
#endif

      if ( l_heat ) then
         call create_lm2r_type(lo2r_s,2)
         if ( l_TP_form ) then
            call create_r2lm_type(r2lo_s,3) ! since we also need u\grad P
         else
            call create_r2lm_type(r2lo_s,2)
         end if
      end if
      if ( l_chemical_conv ) then
         call create_lm2r_type(lo2r_xi,2)
         call create_r2lm_type(r2lo_xi,2)
      end if
      if ( l_conv .or. l_mag_kin) then
         call create_lm2r_type(lo2r_flow,7)
         if ( l_double_curl ) then
            call create_r2lm_type(r2lo_flow,4)
         else
            call create_r2lm_type(r2lo_flow,3)
         end if
      end if

      if ( l_mag ) then
         call create_lm2r_type(lo2r_field,5)
         call create_r2lm_type(r2lo_b,3)
      end if


      ! allocate a temporary array for the gather operations.
      allocate(temp_r2lo(lm_max,l_r:u_r))
      bytes_allocated = bytes_allocated + &
                        lm_max*(u_r-l_r+1)*SIZEOF_DEF_COMPLEX
      if ( coord_r == 0 ) then
         allocate(temp_gather_lo(1:lm_max))
         bytes_allocated = bytes_allocated + lm_max*SIZEOF_DEF_COMPLEX
      else
         allocate(temp_gather_lo(1))
      end if

      local_bytes_used = bytes_allocated - local_bytes_used
      call memWrite('communications.f90', local_bytes_used)

   end subroutine initialize_communications
!-------------------------------------------------------------------------------
   subroutine finalize_communications

      call destroy_gather_type(gt_OC)
      call destroy_gather_type(gt_IC)

#ifdef WITH_MPI
      deallocate( s_request, r_request, array_of_statuses, final_wait_array )
      deallocate( s_transfer_type, s_transfer_type_nr_end )
      deallocate( r_transfer_type, r_transfer_type_nr_end )
      deallocate( s_transfer_type_cont, s_transfer_type_nr_end_cont )
      deallocate( r_transfer_type_cont, r_transfer_type_nr_end_cont )
#endif

      if ( l_heat ) then
         call destroy_lm2r_type(lo2r_s)
         if ( l_TP_form ) then
            call destroy_r2lm_type(r2lo_s) ! since we also need u\grad P
         else
            call destroy_r2lm_type(r2lo_s)
         end if
      end if
      if ( l_chemical_conv ) then
         call destroy_lm2r_type(lo2r_xi)
         call destroy_r2lm_type(r2lo_xi)
      end if
      if ( l_conv .or. l_mag_kin) then
         call destroy_lm2r_type(lo2r_flow)
         call destroy_r2lm_type(r2lo_flow)
      end if

      if ( l_mag ) then
         call destroy_lm2r_type(lo2r_field)
         call destroy_r2lm_type(r2lo_b)
      end if

      deallocate( temp_r2lo, temp_gather_lo )

   end subroutine finalize_communications
!-------------------------------------------------------------------------------
   complex(cp) function get_global_sum_cmplx_2d(dwdt_local) result(global_sum)

      complex(cp), intent(in) :: dwdt_local(:,:)
      
#ifdef WITH_MPI
      integer :: ierr
      complex(cp) :: local_sum
      
      local_sum = sum( dwdt_local )
      call MPI_Allreduce(local_sum,global_sum,1,MPI_DEF_COMPLEX, &
                         MPI_SUM,comm_r,ierr)
#else
      global_sum= sum(dwdt_local)
#endif

   end function get_global_sum_cmplx_2d
!-------------------------------------------------------------------------------
   complex(cp) function get_global_sum_cmplx_2d_dist(dwdt_local) result(global_sum)

      complex(cp), intent(in) :: dwdt_local(:,:)
      
#ifdef WITH_MPI
      integer :: ierr
      complex(cp) :: local_sum
      
      local_sum = sum( dwdt_local )
      call MPI_Allreduce(local_sum,global_sum,1,MPI_DEF_COMPLEX, &
                         MPI_SUM,comm_gs,ierr)
#else
      global_sum= sum(dwdt_local)
#endif

   end function get_global_sum_cmplx_2d_dist
!-------------------------------------------------------------------------------
   real(cp) function get_global_sum_real_2d(dwdt_local) result(global_sum)

      real(cp), intent(in) :: dwdt_local(:,:)
      
#ifdef WITH_MPI
      integer :: ierr
      real(cp) :: local_sum
      
      local_sum = sum( dwdt_local )
      call MPI_Allreduce(local_sum,global_sum,1,MPI_DEF_REAL,MPI_SUM,comm_r,ierr)
#else
      global_sum= sum(dwdt_local)
#endif

   end function get_global_sum_real_2d
!-------------------------------------------------------------------------------
   real(cp) function get_global_sum_real_2d_dist(dwdt_local) result(global_sum)

      real(cp), intent(in) :: dwdt_local(:,:)
      
#ifdef WITH_MPI
      integer :: ierr
      real(cp) :: local_sum
      
      local_sum = sum( dwdt_local )
      call MPI_Allreduce(local_sum,global_sum,1,MPI_DEF_REAL,MPI_SUM,&
                         comm_gs,ierr)
#else
      global_sum= sum(dwdt_local)
#endif

   end function get_global_sum_real_2d_dist
!-------------------------------------------------------------------------------
   complex(cp) function get_global_sum_cmplx_1d(arr_local) result(global_sum)
      !
      ! Kahan summation algorithm
      !
      ! .. code-block:: c
      !
      !    function KahanSum(input)
      !    var sum = 0.0
      !    var c = 0.0             //A running compensation for lost low-order bits.
      !    for i = 1 to input.length do
      !       y = input[i] - c    //So far, so good: c is zero.
      !       t = sum + y         //Alas, sum is big, y small, 
      !                           //so low-order digits of y are lost.
      !       c = (t - sum) - y   //(t - sum) recovers the high-order part of y; 
      !                           //subtracting y recovers -(low part of y)
      !       sum = t             //Algebraically, c should always be zero. 
      !                           //Beware eagerly optimising compilers!
      !       //Next time around, the lost low part will be added to y in a fresh attempt.
      !    return sum
      !
      ! ..
      !
      
      complex(cp), intent(in) :: arr_local(:)
      
#ifdef WITH_MPI
      integer :: lb,ub,ierr,i
      complex(cp) :: local_sum,c,y,t

      lb = lbound(arr_local,1)
      ub = ubound(arr_local,1)
      local_sum = 0.0_cp
      c = 0.0_cp          !A running compensation for lost low-order bits.
      do i=lb,ub
         y = arr_local(i) - c ! So far, so good: c is zero.
         t = local_sum + y    ! Alas, sum is big, y small, 
                              ! so low-order digits of y are lost.
         c = (t - local_sum) - y ! (t - sum) recovers the high-order part of y; 
                                 ! subtracting y recovers -(low part of y)
         local_sum = t           ! Algebraically, c should always be zero. 
                                 ! Beware eagerly optimising compilers
         ! Next time around, the lost low part will be added to y in a fresh attempt.
      end do

      !local_sum = sum( arr_local )
      call MPI_Allreduce(local_sum,global_sum,1,MPI_DEF_COMPLEX, &
                          MPI_SUM,comm_r,ierr)
#else
      global_sum = sum( arr_local )
#endif

   end function get_global_sum_cmplx_1d
!-------------------------------------------------------------------------------
   complex(cp) function get_global_sum_cmplx_1d_dist(arr_local) result(global_sum)
      !
      ! Kahan summation algorithm
      ! 
      ! Read the comment of the previous function
      ! 
      
      complex(cp), intent(in) :: arr_local(:)
      
#ifdef WITH_MPI
      integer :: lb,ub,ierr,i
      complex(cp) :: local_sum,c,y,t

      lb = lbound(arr_local,1)
      ub = ubound(arr_local,1)
      local_sum = 0.0_cp
      c = 0.0_cp          !A running compensation for lost low-order bits.
      do i=lb,ub
         y = arr_local(i) - c ! So far, so good: c is zero.
         t = local_sum + y    ! Alas, sum is big, y small, 
                              ! so low-order digits of y are lost.
         c = (t - local_sum) - y ! (t - sum) recovers the high-order part of y; 
                                 ! subtracting y recovers -(low part of y)
         local_sum = t           ! Algebraically, c should always be zero. 
                                 ! Beware eagerly optimising compilers
         ! Next time around, the lost low part will be added to y in a fresh attempt.
      end do

      !local_sum = sum( arr_local )
      call MPI_Allreduce(local_sum,global_sum,1,MPI_DEF_COMPLEX, &
                          MPI_SUM,comm_gs,ierr)
#else
      global_sum = sum( arr_local )
#endif

   end function get_global_sum_cmplx_1d_dist
!-------------------------------------------------------------------------------
   subroutine gather_all_from_lo_to_rank0(self,arr_lo,arr_full)

      type(gather_type) :: self
      complex(cp) :: arr_lo(llm:ulm,self%dim2)
      complex(cp) :: arr_full(1:lm_max,self%dim2)
      
      integer :: l,m,nR
#ifdef WITH_MPI
      integer :: ierr,irank
      !complex(cp) :: temp_lo((1:lm_max,self%dim2)
      complex(cp), allocatable :: temp_lo(:,:)
      integer :: gather_tag,status(MPI_STATUS_SIZE)

      if ( coord_r == 0 ) allocate(temp_lo(1:lm_max,self%dim2))
      if (n_ranks_r == 1) then
         ! copy the data on coord_r 0
         do nR=1,self%dim2
            temp_lo(llm:ulm,nR)=arr_lo(:,nR)
         end do
      else
         !call MPI_Barrier(comm_r,ierr)
         gather_tag=1990
         if ( coord_r == 0 ) then
            do irank=1,n_ranks_r-1
               call MPI_Recv(temp_lo(lmStartB(irank+1),1),1,       &
                    & self%gather_mpi_type(irank),irank,gather_tag,&
                    & comm_r,status,ierr)
            end do

            ! copy the data on coord_r 0
            do nR=1,self%dim2
               temp_lo(llm:ulm,nR)=arr_lo(:,nR)
            end do
            !write(*,"(A,I3,A,2ES22.14)") "recving temp_lo(",1+irank*lm_per_rank,") &
            !     & = ",sum(temp_lo(1+irank*lm_per_rank:1+irank*lm_per_rank+        &
            !     & lm_on_last_rank,:))
         else
            ! Now send the data to coord_r 0
            !write(*,"(A,I5,A,I2)") "Sending ",(ulm-llm+1)*self%dim2," &
            !   &    dc from coord_r ",coord_r
            !write(*,"(A,2ES22.14)") "sending arr_lo = ", sum(arr_lo)
            call MPI_Send(arr_lo,self%dim2*(ulm-llm+1),MPI_DEF_COMPLEX, &
                          0,gather_tag,comm_r,ierr)
         end if
         !call MPI_Barrier(comm_r,ierr)
      end if

      if ( coord_r == 0 ) then    
         ! reorder
         if ( .not. l_axi ) then
            do nR=1,self%dim2
               do l=0,l_max
                  do m=0,l,minc
                     arr_full(map_glbl_st%lm2(l,m),nR) = temp_lo(lo_map%lm2(l,m),nR)
                  end do
               end do
            end do
         else
            do nR=1,self%dim2
               do l=0,l_max
                  arr_full(map_glbl_st%lm2(l,0),nR) = temp_lo(lo_map%lm2(l,0),nR)
               end do
            end do
         end if
         deallocate(temp_lo)
      end if
#else
      if ( .not. l_axi ) then
         do nR=1,self%dim2
            do l=0,l_max
               do m=0,l,minc
                  arr_full(map_glbl_st%lm2(l,m),nR) = arr_lo(lo_map%lm2(l,m),nR)
               end do
            end do
         end do
      else
         do nR=1,self%dim2
            do l=0,l_max
               arr_full(map_glbl_st%lm2(l,0),nR) = arr_lo(lo_map%lm2(l,0),nR)
            end do
         end do
      end if
#endif

   end subroutine gather_all_from_lo_to_rank0
!-------------------------------------------------------------------------------
   subroutine create_gather_type(self,dim2)
      !
      ! Define the datatypes for gather_all_from_lo_to_rank0
      ! the sending array has dimension (llm:ulm,1:dim2)
      ! receiving array has dimension (1:lm_max,1:dim2)
      !

      type(gather_type) :: self
      integer :: dim2

      integer :: proc

#ifdef WITH_MPI
      allocate(self%gather_mpi_type(0:n_ranks_r-1))
      ! 1. Datatype for the data on one coord_r 
      do proc=0,n_ranks_r-1
         call MPI_type_vector(dim2,lmStopB(proc+1)-lmStartB(proc+1)+1,&
              &               lm_max,MPI_DEF_COMPLEX,              &
              &               self%gather_mpi_type(proc),ierr)
         call MPI_Type_commit(self%gather_mpi_type(proc),ierr)
         if (ierr /= MPI_SUCCESS) print *,"~~~~~>",rank,__LINE__,__FILE__
      end do
#endif
      ! 2. Datatype for the data on the last coord_r
      !call MPI_Type_vector(dim2,lmStopB(n_ranks_r)-lmStartB(n_ranks_r)+1,&
      !     &lm_max,MPI_DEF_COMPLEX,&
      !     & self%gather_mpi_type_end,ierr)
      !call MPI_Type_commit(self%gather_mpi_type_end,ierr)
      self%dim2=dim2

   end subroutine create_gather_type
!-------------------------------------------------------------------------------
   subroutine destroy_gather_type(self)

      type(gather_type) :: self

      integer :: proc

#ifdef WITH_MPI
      do proc=0,n_ranks_r-1
         call MPI_Type_free(self%gather_mpi_type(proc),ierr)
      end do

      deallocate(self%gather_mpi_type)
#endif

   end subroutine destroy_gather_type
!-------------------------------------------------------------------------------
   subroutine gather_from_lo_to_rank0(arr_lo,arr_full)

      complex(cp) :: arr_lo(llm:ulm)
      complex(cp) :: arr_full(1:lm_max)

      integer :: l,m
#ifdef WITH_MPI
      integer :: sendcounts(0:n_ranks_r-1),displs(0:n_ranks_r-1)
      integer :: irank
      !complex(cp) :: temp_lo(1:lm_max)

      do irank=0,n_ranks_r-1
         sendcounts(irank) = lmStopB(irank+1)-lmStartB(irank+1)+1
         displs(irank) = lmStartB(irank+1)-1 !irank*lm_per_rank
      end do
      !sendcounts(n_ranks_r-1) = lm_on_last_rank
      
      call MPI_GatherV(arr_lo,sendcounts(coord_r),MPI_DEF_COMPLEX,&
           &           temp_gather_lo,sendcounts,displs,          &
           &           MPI_DEF_COMPLEX,0,comm_r,ierr)

      if ( coord_r == 0 ) then
         ! reorder
         if ( .not. l_axi ) then
            do l=0,l_max
               do m=0,l,minc
                  arr_full(map_glbl_st%lm2(l,m)) = temp_gather_lo(lo_map%lm2(l,m))
               end do
            end do
         else
            do l=0,l_max
               arr_full(map_glbl_st%lm2(l,0)) = temp_gather_lo(lo_map%lm2(l,0))
            end do
         end if
      end if
#else
      if ( .not. l_axi ) then
         do l=0,l_max
            do m=0,l,minc
               arr_full(map_glbl_st%lm2(l,m)) = arr_lo(lo_map%lm2(l,m))
            end do
         end do
      else
         do l=0,l_max
            arr_full(map_glbl_st%lm2(l,0)) = arr_lo(lo_map%lm2(l,0))
         end do
      end if
#endif
    
   end subroutine gather_from_lo_to_rank0
!-------------------------------------------------------------------------------
   subroutine scatter_from_rank0_to_lo(arr_full,arr_lo)

      complex(cp) :: arr_full(1:lm_max)
      complex(cp) :: arr_lo(llm:ulm)

      integer :: l,m
#ifdef WITH_MPI
      integer :: sendcounts(0:n_ranks_r-1),displs(0:n_ranks_r-1)
      integer :: irank
      !complex(cp) :: temp_lo(1:lm_max)

      do irank=0,n_ranks_r-1
         sendcounts(irank) = lmStopB(irank+1)-lmStartB(irank+1)+1
         displs(irank) = lmStartB(irank+1)-1
      end do

      if ( coord_r == 0 ) then
         ! reorder
         if ( .not. l_axi ) then
            do l=0,l_max
               do m=0,l,minc
                  temp_gather_lo(lo_map%lm2(l,m)) = arr_full(map_glbl_st%lm2(l,m))
               end do
            end do
         else
            do l=0,l_max
               temp_gather_lo(lo_map%lm2(l,0)) = arr_full(map_glbl_st%lm2(l,0))
            end do
         end if
      end if

      call MPI_ScatterV(temp_gather_lo,sendcounts,displs,MPI_DEF_COMPLEX,&
           &            arr_lo,sendcounts(coord_r),MPI_DEF_COMPLEX,0,       &
           &            comm_r,ierr)
#else
      if ( .not. l_axi ) then
         do l=0,l_max
            do m=0,l,minc
               arr_lo(lo_map%lm2(l,m)) = arr_full(map_glbl_st%lm2(l,m))
            end do
         end do
      else
         do l=0,l_max
            arr_lo(lo_map%lm2(l,0)) = arr_full(map_glbl_st%lm2(l,0))
         end do
      end if
#endif

   end subroutine scatter_from_rank0_to_lo
!-------------------------------------------------------------------------------
   subroutine create_lm2r_type(self,count)

      type(lm2r_type) :: self
      integer, optional, intent(in) :: count

      if (.not. present(count)) then
         self%count=1
      else
         self%count = count
      end if
#ifdef WITH_MPI
      allocate(self%s_request(n_ranks_r-1))
      allocate(self%r_request(n_ranks_r-1))
      allocate(self%final_wait_array(2*(n_ranks_r-1)))
      bytes_allocated = bytes_allocated+4*(n_ranks_r-1)*SIZEOF_INTEGER
#endif
      allocate(self%temp_Rloc(1:lm_max,l_r:u_r,1:self%count))
      bytes_allocated = bytes_allocated+&
                        lm_max*(u_r-l_r+1)*self%count*SIZEOF_DEF_COMPLEX

   end subroutine create_lm2r_type
!-------------------------------------------------------------------------------
   subroutine destroy_lm2r_type(self)

      type(lm2r_type) :: self

#ifdef WITH_MPI
      deallocate(self%final_wait_array)
      deallocate(self%s_request)
      deallocate(self%r_request)
#endif
      deallocate(self%temp_Rloc)

   end subroutine destroy_lm2r_type
!-------------------------------------------------------------------------------
   subroutine create_r2lm_type(self,count)

      type(r2lm_type) :: self
      integer, optional, intent(in) :: count

      if (.not. present(count)) then
         self%count=1
      else
         self%count = count
      end if
#ifdef WITH_MPI
      allocate(self%s_request(n_ranks_r-1))
      allocate(self%r_request(n_ranks_r-1))
      allocate(self%final_wait_array(2*(n_ranks_r-1)))
      bytes_allocated = bytes_allocated+4*(n_ranks_r-1)*SIZEOF_INTEGER
#endif
      allocate(self%temp_Rloc(1:lm_max,l_r:u_r,1:self%count))
      bytes_allocated = bytes_allocated+&
                        lm_max*(u_r-l_r+1)*self%count*SIZEOF_DEF_COMPLEX

   end subroutine create_r2lm_type
!-------------------------------------------------------------------------------
   subroutine destroy_r2lm_type(self)

      type(r2lm_type) :: self

#ifdef WITH_MPI
      deallocate(self%s_request)
      deallocate(self%r_request)
      deallocate(self%final_wait_array)
#endif
      deallocate(self%temp_Rloc)

   end subroutine destroy_r2lm_type
!-------------------------------------------------------------------------------
  ! --------------------- NONBLOCKING ---------------------
  ! Here comes the nonblocking variant
   subroutine lm2r_redist_start(self,arr_LMloc,arr_Rloc)

      type(lm2r_type) :: self
      complex(cp), intent(in)  :: arr_LMloc(llm:ulm,n_r_max,*)
      complex(cp), intent(out) :: arr_Rloc(lm_max,l_r:u_r,*)

      integer :: i
#ifdef WITH_MPI
      ! Local variables
      integer :: send_pe,recv_pe,irank, n_r_per_rank
      integer :: transfer_tag=1111

      !PERFON('lm2r_st')
      
      ! -- TODO:
      !    Those were for the old code, which had n_r_loc per rank fixed, 
      !    except for the last rank. Now things are much more flexible.
      !    Eliminate these two variables as soon as the parallelization of the 
      !    ML Loop is completed!
      ! 
      n_r_per_rank = dist_r(0,0)
!       n_r_on_last_rank = dist_r(n_ranks_r-1,0)
      
      if ( coord_r < n_ranks_r-1 ) then
         ! all the ranks from [0,n_ranks_r-2]
         do irank=0,n_ranks_r-1
            !if (coord_r == irank) then
            ! just copy
            !   arr_LMLoc(llm:ulm,l_r:u_r)=arr_Rloc(llm:ulm,l_r:u_r)
            !else
            ! send_pe: send to this coord_r
            ! recv_pe: receive from this coord_r
            send_pe = modulo(coord_r+irank,n_ranks_r)
            recv_pe = modulo(coord_r-irank+n_ranks_r,n_ranks_r)
            !print*,"send to ",send_pe,",     recv from ",recv_pe
            if (coord_r == send_pe) then
               !PERFON('loc_copy')
               ! just copy
               do i=1,self%count
                  arr_Rloc(llm:ulm,l_r:u_r,i)= &
                       arr_LMloc(llm:ulm,l_r:u_r,i)
               end do
               !PERFOFF
            else
               !PERFON('irecv')
               !if (recv_pe == n_ranks_r-1) then
               !   call MPI_Irecv(arr_Rloc(lmStartB(recv_pe+1),l_r,1),   &
               !        &         1,s_transfer_type_cont(n_ranks_r,self%count),&
               !        &         recv_pe,transfer_tag,comm_r,       &
               !        &         self%r_request(irank),ierr)
               !else
                  call MPI_Irecv(arr_Rloc(lmStartB(recv_pe+1),l_r,1),     &
                       &         1,s_transfer_type_cont(recv_pe+1,self%count),&
                       &         recv_pe,transfer_tag,comm_r,         &
                       &         self%r_request(irank),ierr)
               !end if
               !PERFOFF
               !PERFON('isend')
               if (send_pe == n_ranks_r-1) then
                  call MPI_Isend(arr_LMloc(llm,1+n_r_per_rank*send_pe,1),          &
                       &         1,r_transfer_type_nr_end_cont(coord_r+1,self%count),&
                       &         send_pe,transfer_tag,comm_r,             &
                       &         self%s_request(irank),ierr)
               else
                  call MPI_Isend(arr_LMloc(llm,1+n_r_per_rank*send_pe,1),   &
                       &         1,r_transfer_type_cont(coord_r+1,self%count),&
                       &         send_pe,transfer_tag,comm_r,      &
                       &         self%s_request(irank),ierr)
               end if
               !PERFOFF
            end if
         end do

         i=1
         do irank=1,n_ranks_r-1
            self%final_wait_array(i)=self%s_request(irank)
            self%final_wait_array(i+1)=self%r_request(irank)
            i = i + 2
         end do
         !print*,"Waiting for completion of nonblocking communication 1"
         !call mpi_waitall(2*n_ranks_r,final_wait_array,array_of_statuses,ierr)
         !print*,"Nonblocking communication 1 is done."
      else
         ! coord_r  ==  n_ranks_r-1
         ! all receives are with the s_transfer_type_nr_end
         ! all sends are done with r_transfer_type_lm_end
         do irank=0,n_ranks_r-1
            ! send_pe: send to this coord_r
            ! recv_pe: receive from this coord_r
            send_pe = modulo(coord_r+irank,n_ranks_r)
            recv_pe = modulo(coord_r-irank+n_ranks_r,n_ranks_r)
            !print*,"send to ",send_pe,",     recv from ",recv_pe
            if (coord_r == send_pe) then
               !PERFON('loc_copy')
               ! just copy
               do i=1,self%count
                  arr_Rloc(llm:ulm,l_r:u_r,i)= &
                        arr_LMloc(llm:ulm,l_r:u_r,i)
               end do
               !PERFOFF
            else
               !PERFON('irecv')
               call MPI_Irecv(arr_Rloc(lmStartB(recv_pe+1),l_r,1),            &
                    &         1,s_transfer_type_nr_end_cont(recv_pe+1,self%count),&
                    &         recv_pe,transfer_tag,comm_r,                &
                    &         self%r_request(irank),ierr)
               !PERFOFF
               !PERFON('isend')
               call MPI_Isend(arr_LMloc(llm,1+n_r_per_rank*send_pe,1),    &
                    &         1,r_transfer_type_cont(coord_r+1,self%count), &
                    &         send_pe,transfer_tag,comm_r,       &
                    &         self%s_request(irank),ierr)
               !PERFOFF
            end if
         end do
         i=1
         do irank=1,n_ranks_r-1
            self%final_wait_array(i)=self%s_request(irank)
            self%final_wait_array(i+1)=self%r_request(irank)
            i = i + 2
         end do
         !print*,"Waiting for completion of nonblocking communication 1"
         !call mpi_waitall(2*(n_ranks_r-1),final_wait_array,array_of_statuses,ierr)
         !print*,"Nonblocking communication 1 is done."
      end if
      !write(*,"(A,I3)") "lm2r_redist_start on n_ranks_r=",n_ranks_r
      !PERFOFF


#else
      do i=1,self%count
         arr_Rloc(llm:ulm,l_r:u_r,i)= arr_LMloc(llm:ulm,l_r:u_r,i)
      end do
#endif

   end subroutine lm2r_redist_start
!-------------------------------------------------------------------------------
   subroutine lm2r_redist_wait(self)

      type(lm2r_type) :: self
#ifdef WITH_MPI
      integer :: ierr
      integer :: array_of_statuses(MPI_STATUS_SIZE,2*n_ranks_r)

      !PERFON('lm2r_wt')
      !write(*,"(A,I3)") "n_ranks_r = ",n_ranks_r
      !write(*,"(2(A,I3))") "Waiting for ",2*(n_ranks_r-1)," requests,", &
      !   &             size(self%final_wait_array)
      call MPI_Waitall(2*(n_ranks_r-1),self%final_wait_array,array_of_statuses,ierr)
      !PERFOFF
#endif

   end subroutine lm2r_redist_wait
!-------------------------------------------------------------------------------
   subroutine lo2r_redist_start(self,arr_lo,arr_Rloc)

      type(lm2r_type) :: self
      complex(cp), intent(in) :: arr_lo(llm:ulm,1:n_r_max,*)
      complex(cp), target, intent(out) :: arr_Rloc(1:lm_max,l_r:u_r,*)
  
  
      PERFON('lo2r_st')
      !call lm2r_redist(arr_lo,temp_lo)
      self%arr_Rloc(1:,l_r:,1:) => arr_Rloc(1:lm_max,l_r:u_r,1:self%count)
      !self%arr_Rloc(1:,l_r:) => arr_Rloc(1:,l_r:)
      call lm2r_redist_start(self,arr_lo,self%temp_Rloc)
      PERFOFF

   end subroutine lo2r_redist_start
!-------------------------------------------------------------------------------
   subroutine lo2r_redist_start_dist(self,arr_lo,arr_dist)

      type(lm2r_type) :: self
      complex(cp), intent(in) :: arr_lo(llm:ulm,1:n_r_max,*)
      complex(cp), target, intent(out) :: arr_dist(1:n_lm_loc,l_r:u_r,*)
  
  
      PERFON('lo2r_st')
      !call lm2r_redist(arr_lo,temp_lo)
      self%arr_Rloc(1:,l_r:,1:) => arr_dist(1:n_lm_loc,l_r:u_r,1:self%count)
      call lm2r_redist_start(self,arr_lo,self%temp_Rloc)
      PERFOFF

   end subroutine lo2r_redist_start_dist
!-------------------------------------------------------------------------------
   subroutine lo2r_redist_wait(self)

      type(lm2r_type) :: self
  
      ! Local variables
      integer :: nR,l,m,i
  
      !PERFON("lo2r_wt")
      call lm2r_redist_wait(self)
      ! now in self%temp_Rloc we do have the lo_ordered r-local part
      ! now reorder to the original ordering
      if ( .not. l_axi ) then
         do i=1,self%count
            do nR=l_r,u_r
               do l=0,l_max
                  do m=0,l,minc
                     self%arr_Rloc(map_glbl_st%lm2(l,m),nR,i) = &
                            self%temp_Rloc(lo_map%lm2(l,m),nR,i)
                  end do
               end do
            end do
         end do
      else
         do i=1,self%count
            do nR=l_r,u_r
               do l=0,l_max
                  self%arr_Rloc(map_glbl_st%lm2(l,0),nR,i) = &
                         self%temp_Rloc(lo_map%lm2(l,0),nR,i)
               end do
            end do
         end do
      end if
      !PERFOFF

   end subroutine lo2r_redist_wait
!-------------------------------------------------------------------------------
   subroutine lo2r_redist_wait_dist(self)

      type(lm2r_type) :: self
  
      ! Local variables
      integer :: nR,l,m,i,lm
  
      !PERFON("lo2r_wt")
      call lm2r_redist_wait(self)
      ! now in self%temp_Rloc we do have the lo_ordered r-local part
      ! now reorder to the original ordering
      if ( .not. l_axi ) then
         do i=1,self%count
            do nR=l_r,u_r
               do lm=1,n_lm_loc
                  l = map_dist_st%lm2l(lm)
                  m = map_dist_st%lm2m(lm)
                  self%arr_Rloc(lm,nR,i) = self%temp_Rloc(lo_map%lm2(l,m),nR,i)
               end do
            end do
         end do
      else
         do i=1,self%count
            do nR=l_r,u_r
               if (map_dist_st%lm2(0,0) > 0) then
                  do l=0,l_max
                     self%arr_Rloc(map_dist_st%lm2(l,0),nR,i) = &
                           self%temp_Rloc(lo_map%lm2(l,0),nR,i)
                  end do
               end if
            end do
         end do
      end if
      !PERFOFF

   end subroutine lo2r_redist_wait_dist
!-------------------------------------------------------------------------------
   subroutine r2lm_redist_start(self,arr_rloc,arr_LMloc)

      type(r2lm_type) :: self
      complex(cp), intent(in) :: arr_Rloc(lm_max,l_r:u_r,*)
      complex(cp), intent(out) :: arr_LMloc(llm:ulm,n_r_max,*)
  
      integer :: i
#ifdef WITH_MPI
      ! Local variables
      integer :: send_pe, recv_pe, irank
      integer :: transfer_tag=1111
  
      !write(*,"(A)") "----------- start r2lm_redist -------------"
      !PERFON('r2lm_dst')
      if (coord_r < n_ranks_r-1) then
         ! all the ranks from [0,n_ranks_r-2]
         do irank=0,n_ranks_r-1
            !if (coord_r == irank) then
            ! just copy
            !   arr_LMLoc(llm:ulm,l_r:u_r)=arr_Rloc(llm:ulm,l_r:u_r)
            !else
            ! send_pe: send to this coord_r
            ! recv_pe: receive from this coord_r
            send_pe = modulo(coord_r+irank,n_ranks_r)
            recv_pe = modulo(coord_r-irank+n_ranks_r,n_ranks_r)
            if (coord_r == send_pe) then
               do i=1,self%count
                  arr_LMLoc(llm:ulm,l_r:u_r,i)= &
                            arr_Rloc(llm:ulm,l_r:u_r,i)
               end do
            else
               call MPI_Isend(arr_Rloc(lmStartB(send_pe+1),l_r,1),     &
                    &         1,s_transfer_type_cont(send_pe+1,self%count),&
                    &         send_pe,transfer_tag,comm_r,         &
                    &         self%s_request(irank),ierr)
               if (recv_pe == n_ranks_r-1) then
                  !-- TODO
                  !   This is using dist_r(0,0) instead of nR_per_rank. 
                  !   dunno for how long this will work
                  call MPI_Irecv(arr_LMloc(llm,1+dist_r(0,0)*recv_pe,1),          &
                       &         1,r_transfer_type_nr_end_cont(coord_r+1,self%count),&
                       &         recv_pe,transfer_tag,comm_r,             &
                       &         self%r_request(irank),ierr)
               else
                  call MPI_Irecv(arr_LMloc(llm,1+dist_r(0,0)*recv_pe,1),     &
                       &         1,r_transfer_type_cont(coord_r+1,self%count),  &
                       &         recv_pe,transfer_tag,comm_r,        &
                       &         self%r_request(irank),ierr)
               end if
            end if
         end do
  
         i=1
         do irank=1,n_ranks_r-1
            self%final_wait_array(i)=self%s_request(irank)
            self%final_wait_array(i+1)=self%r_request(irank)
            i = i + 2
         end do
      else
         ! coord_r  ==  n_ranks_r-1
         ! all receives are with the r_transfer_type_lm_end
         ! all sends are done with s_transfer_type_nr_end
         do irank=0,n_ranks_r-1
            ! send_pe: send to this coord_r
            ! recv_pe: receive from this coord_r
            send_pe = modulo(coord_r+irank,n_ranks_r)
            recv_pe = modulo(coord_r-irank+n_ranks_r,n_ranks_r)
            !print*,"send to ",send_pe,",     recv from ",recv_pe
            if (coord_r == send_pe) then
               ! just copy
               do i=1,self%count
                  arr_LMLoc(llm:ulm,l_r:u_r,i)= &
                  &        arr_Rloc(llm:ulm,l_r:u_r,i)
               end do
            else
               !-- TODO
               !   This is using dist_r(0,0) instead of nR_per_rank. 
               !   dunno for how long this will work
               call MPI_Irecv(arr_LMloc(llm,1+dist_r(0,0)*recv_pe,1),    &
                    &         1,r_transfer_type_cont(coord_r+1,self%count), &
                    &         recv_pe,transfer_tag,comm_r,       &
                    &         self%r_request(irank),ierr)
               call MPI_Isend(arr_Rloc(lmStartB(send_pe+1),l_r,1),            &
                    &         1,s_transfer_type_nr_end_cont(send_pe+1,self%count),&
                    &         send_pe,transfer_tag,comm_r,                &
                    &         self%s_request(irank),ierr)
  
            end if
         end do
         i=1
         do irank=1,n_ranks_r-1
            self%final_wait_array(i)=self%s_request(irank)
            self%final_wait_array(i+1)=self%r_request(irank)
            i = i + 2
         end do
      end if
  
      !PERFOFF
      !write(*,"(A)") "----------- end   r2lm_redist -------------"
#else
      do i=1,self%count
         arr_LMLoc(llm:ulm,l_r:u_r,i)=arr_Rloc(llm:ulm,l_r:u_r,i)
      end do
#endif

   end subroutine r2lm_redist_start
!-------------------------------------------------------------------------------
   subroutine r2lm_redist_wait(self)

      type(r2lm_type) :: self
#ifdef WITH_MPI
      integer :: ierr
      integer :: array_of_statuses(MPI_STATUS_SIZE,2*n_ranks_r)

      !PERFON('lm2r_wt')
      !write(*,"(A,I3)") "n_ranks_r = ",n_ranks_r
      !write(*,"(2(A,I3))") "Waiting for ",2*(n_ranks_r-1)," requests,", &
      !   &             size(self%final_wait_array)
      call MPI_Waitall(2*(n_ranks_r-1),self%final_wait_array,array_of_statuses,ierr)
      !PERFOFF
#endif

   end subroutine r2lm_redist_wait
!-------------------------------------------------------------------------------
   subroutine r2lo_redist_start(self,arr_Rloc,arr_lo)

      type(r2lm_type) :: self
      complex(cp), intent(in) :: arr_Rloc(1:lm_max,l_r:u_r,*)
      complex(cp), intent(out) :: arr_lo(llm:ulm,1:n_r_max,*)
  
      ! Local variables
      integer :: nR,l,m,i

      self%temp_Rloc(1:,l_r:,1:) = arr_Rloc(1:lm_max,l_r:u_r,1:self%count)
  
      ! Just copy the array with permutation
      !PERFON('r2lo_dst')
      if ( .not. l_axi ) then
         do i=1,self%count
            do nR=l_r,u_r
               do l=0,l_max
                  do m=0,l,minc
                     self%temp_Rloc(lo_map%lm2(l,m),nR,i) = & 
                                      arr_Rloc(map_glbl_st%lm2(l,m),nR,i)
                  end do
               end do
            end do
         end do
      else
         do i=1,self%count
            do nR=l_r,u_r
               do l=0,l_max
                  self%temp_Rloc(lo_map%lm2(l,0),nR,i) = & 
                                   arr_Rloc(map_glbl_st%lm2(l,0),nR,i)
               end do
            end do
         end do
      end if
  
      call r2lm_redist_start(self,self%temp_Rloc,arr_lo)
      !PERFOFF

   end subroutine r2lo_redist_start
!-------------------------------------------------------------------------------
   subroutine r2lo_redist_start_dist(self,arr_dist,arr_lo)
      !
      ! This is just a temporary cheat for performing the transposition.
      !>@TODO Do a true paralellization in θ of this function
      ! 

      type(r2lm_type) :: self
      complex(cp), intent(in) :: arr_dist(1:n_lm_loc,l_r:u_r,*)
      complex(cp), intent(out) :: arr_lo(llm:ulm,1:n_r_max,*)
  
      ! Local variables
      integer :: nR,l,m,i
      complex(cp) :: tmp_glb(1:lm_max)
      

!       self%temp_Rloc(1:,l_r:,1:) = arr_Rloc(1:lm_max,l_r:u_r,1:self%count)
  
      ! Just copy the array with permutation
      !PERFON('r2lo_dst')
      if ( .not. l_axi ) then
         do i=1,self%count
            do nR=l_r,u_r
               call gather_Flm(arr_dist(1:n_lm_loc,nR,i), tmp_glb(1:lm_max))
               self%temp_Rloc(1:lm_max,nR,i) = tmp_glb(1:lm_max)
               do l=0,l_max
                  do m=0,l,minc
                     self%temp_Rloc(lo_map%lm2(l,m),nR,i) = tmp_glb(map_glbl_st%lm2(l,m))
                  end do
               end do
            end do
         end do
!          stop
      else
         do i=1,self%count
            do nR=l_r,u_r
               call gather_Flm(arr_dist(1:n_lm_loc,nR,i), tmp_glb(1:lm_max))
               self%temp_Rloc(1:lm_max,nR,i) = tmp_glb(1:lm_max)
               do l=0,l_max
                  self%temp_Rloc(lo_map%lm2(l,0),nR,i) = tmp_glb(map_glbl_st%lm2(l,0))
               end do
            end do
         end do
      end if
  
      call r2lm_redist_start(self,self%temp_Rloc,arr_lo)
      !PERFOFF

   end subroutine r2lo_redist_start_dist
!-------------------------------------------------------------------------------
   subroutine r2lo_redist_wait(self)

      type(r2lm_type) :: self
  
      !PERFON("r2lo_wt")
      call r2lm_redist_wait(self)
      !PERFOFF

   end subroutine r2lo_redist_wait
!-------------------------------------------------------------------------------
   subroutine lm2lo_redist(arr_LMloc,arr_lo)

      complex(cp), intent(in) :: arr_LMloc(llm:ulm,1:n_r_max)
      complex(cp), intent(out) :: arr_lo(llm:ulm,1:n_r_max)
  
      ! Local variables
      integer :: nR,l,m
  
      if (n_ranks_r > 1) then
         call abortRun('lm2lo not yet parallelized')
      end if
  
      if ( .not. l_axi ) then
         do nR=1,n_r_max
            do l=0,l_max
               do m=0,l,minc
                  arr_lo(lo_map%lm2(l,m),nR) = arr_LMloc(map_glbl_st%lm2(l,m),nR)
               end do
            end do
         end do
      else
         do nR=1,n_r_max
            do l=0,l_max
               arr_lo(lo_map%lm2(l,0),nR) = arr_LMloc(map_glbl_st%lm2(l,0),nR)
            end do
         end do
      end if
    
   end subroutine lm2lo_redist
!-------------------------------------------------------------------------------
   subroutine lo2lm_redist(arr_lo,arr_LMloc)

      complex(cp), intent(in) :: arr_lo(llm:ulm,1:n_r_max)
      complex(cp), intent(out) :: arr_LMloc(llm:ulm,1:n_r_max)
  
      ! Local variables
      integer :: nR,l,m
  
      if (n_ranks_r > 1) then
         call abortRun('lo2lm not yet parallelized')
      end if
  
      if ( .not. l_axi ) then
         do nR=1,n_r_max
            do l=0,l_max
               do m=0,l,minc
                  arr_LMloc(map_glbl_st%lm2(l,m),nR) = arr_lo(lo_map%lm2(l,m),nR)
               end do
            end do
         end do
      else
         do nR=1,n_r_max
            do l=0,l_max
               arr_LMloc(map_glbl_st%lm2(l,0),nR) = arr_lo(lo_map%lm2(l,0),nR)
            end do
         end do
      end if
      
   end subroutine lo2lm_redist
!-------------------------------------------------------------------------------
#ifdef WITH_MPI
#define STANDARD 1001
#define EXTENDED 1002
#define DATATYPE 1003

#define ALLGATHER STANDARD
#if (ALLGATHER==STANDARD)
!-------------------------------------------------------------------------------
   subroutine myAllGather(arr,dim1,dim2)

      use blocking
      use parallel_mod
  
      implicit none
  
      integer,     intent(in) :: dim1,dim2
      complex(cp), intent(inout) :: arr(dim1,dim2)
  
      integer :: lmStart_on_rank,lmStop_on_rank
      integer :: sendcount,recvcounts(0:n_ranks_r-1),displs(0:n_ranks_r-1)
      integer :: irank,nR
      !double precision :: local_sum, global_sum, recvd_sum
  
      !write(*,"(A,ES15.8)") "before: arr = ",sum(real(conjg(arr)*arr))
  
      lmStart_on_rank = lmStartB(1+coord_r*nLMBs_per_rank)
      lmStop_on_rank  = lmStopB(min((coord_r+1)*nLMBs_per_rank,nLMBs)) 
      sendcount  = lmStop_on_rank-lmStart_on_rank+1
      do irank=0,n_ranks_r-1
         recvcounts(irank) = lmStopB ( min((irank+1)*nLMBs_per_rank,nLMBs) ) &
                           - lmStartB( 1+irank*nLMBs_per_rank ) + 1
      end do
      do irank=0,n_ranks_r-1
         !displs(irank) = (nR-1)*dim1 + sum(recvcounts(0:irank-1))
         displs(irank) = sum(recvcounts(0:irank-1))
      end do
      !write(*,"(4X,2I4,A,I6)") coord_r,nR," displs = ",displs(coord_r)
      !write(*,"(5(A,I4))") "LMBlocks ",1+coord_r*nLMBs_per_rank,"->", &
      !     &                (coord_r+1)*nLMBs_per_rank,&
      !     &", lm runs from ",lmStart_on_rank," to ",lmStop_on_rank,", &
      !     &  recvcounts = ",recvcounts(coord_r)
      do nR=1,dim2
         !local_sum = sum( real( conjg(arr(lmStart_on_rank:lmStop_on_rank,nR))*&
         !     &           arr(lmStart_on_rank:lmStop_on_rank,nR) ) )
         call MPI_AllGatherV(MPI_IN_PLACE,sendcount,MPI_DEF_COMPLEX,     &
              &              arr(1,nR),recvcounts,displs,MPI_DEF_COMPLEX,&
              &              comm_r,ierr)
         !recvd_sum = sum( real( &
         !     & conjg(arr( lmStartB(1+&
         !     &(modulo(coord_r+1,2)*nLMBs_per_rank)):&
         !     &lmStopB(1+(modulo(coord_r+1,2)+1)*nLMBs_per_rank-1),nR ))*&
         !     &       arr( lmStartB(1+(modulo(coord_r+1,2)*nLMBs_per_rank)):&
         !     &lmStopB(1+(modulo(coord_r+1,2)+1)*nLMBs_per_rank-1),nR ) ) )
         !global_sum = sum( real( conjg(arr(:,nR))*arr(:,nR) ) )
         !write(*,"(4X,A,I4,3(A,ES20.13))") "nR = ",nR,": l_sum = ",&
         !&      local_sum,", r_sum = ",recvd_sum,", g_sum = ", global_sum
      end do
      !write(*,"(A)") "---------------------------"

   end subroutine myAllGather
!-------------------------------------------------------------------------------
#elif (ALLGATHER==EXTENDED)
   subroutine myAllGather(arr,edim1,dim2)

      use blocking
      use parallel_mod

      implicit none

      integer,     intent(in) :: edim1,dim2
      complex(cp), intent(inout) :: arr(edim1, dim2)

      integer :: sendcount,recvcount
      integer :: irank,nR

      recvcount = edim1/n_ranks_r
      do nR=1,dim2
         call MPI_AllGather(MPI_IN_PLACE,sendcount,MPI_DEF_COMPLEX,&
              &             arr(1,nR),recvcount,MPI_DEF_COMPLEX,   &
              &             comm_r,ierr)
      end do

   end subroutine myAllGather
!-------------------------------------------------------------------------------
#elif (ALLGATHER==DATATYPE)
   ! now with a datatype to have one big transfer, not many small ones
   subroutine myAllGather(arr,dim1,dim2)

      use blocking
      use parallel_mod

      implicit none

      integer,     intent(in) :: dim1,dim2
      complex(cp), intent(inout) :: arr(dim1,dim2)

      integer :: lmStart_on_rank,lmStop_on_rank
      integer :: sendcount,recvcounts(0:n_ranks_r-1),displs(0:n_ranks_r-1)
      integer :: irank,nR
      integer :: sendtype, new_sendtype
      integer(kind=MPI_ADDRESS_KIND) :: lb,extent,extent_dcmplx

      PERFON('mk_dt')
      ! definition of the datatype (will later be pulled out of here)
      ! we assume dim1=lm_max and dim2=n_r_max
      call mpi_type_get_extent(MPI_DEF_COMPLEX,lb,extent_dcmplx,ierr)
      call mpi_type_vector(dim2,1,dim1,MPI_DEF_COMPLEX,sendtype,ierr)
      lb=0
      extent=extent_dcmplx
      call mpi_type_create_resized(sendtype,lb,extent,new_sendtype,ierr)
      call mpi_type_commit(new_sendtype,ierr)


      lmStart_on_rank = lmStartB(1+coord_r*nLMBs_per_rank)
      lmStop_on_rank  = lmStopB(1+(coord_r+1)*nLMBs_per_rank-1) 
      sendcount  = lmStop_on_rank-lmStart_on_rank+1
      do irank=0,n_ranks_r-1
         recvcounts(irank) = lmStopB( (irank+1)*nLMBs_per_rank ) - &
                            lmStartB( 1+irank*nLMBs_per_rank ) + 1
      end do
      do irank=0,n_ranks_r-1
         !displs(irank) = (nR-1)*dim1 + sum(recvcounts(0:irank-1))
         displs(irank) = sum(recvcounts(0:irank-1))
      end do
      PERFOFF
      PERFON('comm')
      call MPI_AllGatherV(MPI_IN_PLACE,sendcount,new_sendtype,&
           &              arr,recvcounts,displs,new_sendtype, &
           &              comm_r,ierr)
      PERFOFF

   end subroutine myAllGather
!------------------------------------------------------------------------------
#endif
#endif

   !----------------------------------------------------------------------------
   subroutine slice_f(f_global, f_local)
      !
      !   Author: Rafael Lago (MPCDF) August 2017
      !
      real(cp),  intent(in)  :: f_global(n_phi_max, n_theta_max)
      real(cp),  intent(out) :: f_local(n_phi_max, n_theta_loc)
      
      f_local = f_global(:,l_theta:u_theta)
   end subroutine slice_f
   
   !----------------------------------------------------------------------------
   subroutine slice_Flm_cmplx(Flm_global, Flm_local)
      !
      !   Author: Rafael Lago (MPCDF) August 2017
      !
      complex(cp),  intent(in)  :: Flm_global(lm_max)
      complex(cp),  intent(out) :: Flm_local(n_lm_loc)
      
      integer :: i, l_lm, u_lm, l_lm_g, u_lm_g, m
      
      do i = 1, n_m_loc
        m = dist_m(coord_m,i)
        l_lm  = map_dist_st%lm2(m, m)
        u_lm  = map_dist_st%lm2(l_max, m)
        l_lm_g = map_glbl_st%lm2(m, m)
        u_lm_g = map_glbl_st%lm2(l_max, m)
        Flm_local(l_lm:u_lm) = Flm_global(l_lm_g:u_lm_g)
      end do
   end subroutine slice_Flm_cmplx
   
   !----------------------------------------------------------------------------
   subroutine slice_Flm_real(Flm_global, Flm_local)
      !
      !   Author: Rafael Lago (MPCDF) August 2017
      !
      real(cp),  intent(in)  :: Flm_global(lm_max)
      real(cp),  intent(out) :: Flm_local(n_lm_loc)
      
      integer :: i, l_lm, u_lm, l_lm_g, u_lm_g, m
      
      do i = 1, n_m_loc
        m = dist_m(coord_m,i)
        l_lm  = map_dist_st%lm2(m, m)
        u_lm  = map_dist_st%lm2(l_max, m)
        l_lm_g = map_glbl_st%lm2(m, m)
        u_lm_g = map_glbl_st%lm2(l_max, m)
        Flm_local(l_lm:u_lm) = Flm_global(l_lm_g:u_lm_g)
      end do
   end subroutine slice_Flm_real

   !----------------------------------------------------------------------------
   subroutine slice_FlmP_cmplx(FlmP_global, Flm_local)
      !
      !   Author: Rafael Lago (MPCDF) August 2017
      !
      complex(cp),  intent(in)  :: FlmP_global(lmP_max)
      complex(cp),  intent(out) :: Flm_local(n_lmP_loc)
      
      integer :: i, l_lm, u_lm, l_lm_g, u_lm_g, m
      
      do i = 1, n_m_loc
        m = dist_m(coord_m, i)
        l_lm  = map_dist_st%lmP2(m, m)
        u_lm  = map_dist_st%lmP2(l_max, m)
        l_lm_g = map_glbl_st%lmP2(m,   m)
        u_lm_g = map_glbl_st%lmP2(l_max+1, m)
        Flm_local(l_lm:u_lm) = FlmP_global(l_lm_g:u_lm_g)
      end do
   end subroutine slice_FlmP_cmplx
   
   !----------------------------------------------------------------------------
   subroutine slice_FlmP_real(FlmP_global, Flm_local)
      !
      !   Author: Rafael Lago (MPCDF) August 2017
      !
      real(cp),  intent(in)  :: FlmP_global(lmP_max)
      real(cp),  intent(out) :: Flm_local(n_lmP_loc)
      
      integer :: i, l_lm, u_lm, l_lm_g, u_lm_g, m
      
      do i = 1, n_m_loc
        m = dist_m(coord_m, i)
        l_lm  = map_dist_st%lm2(m, m)
        u_lm  = map_dist_st%lm2(l_max, m)
        l_lm_g = map_glbl_st%lmP2(m,   m)
        u_lm_g = map_glbl_st%lmP2(l_max+1, m)
        Flm_local(l_lm:u_lm) = FlmP_global(l_lm_g:u_lm_g)
      end do
   end subroutine slice_FlmP_real

   !----------------------------------------------------------------------------
   subroutine gather_FlmP(Flm_local, FlmP_global)
      !
      !   Author: Rafael Lago (MPCDF) August 2017
      !
      complex(cp),  intent(in)  :: Flm_local(n_lmP_loc)
      complex(cp),  intent(out) :: FlmP_global(lmP_max)
      
      complex(cp) ::  buffer(lmP_max)
      integer :: irank, j, m, l_lm_loc, u_lm_loc, l_lm_glb, u_lm_glb
      integer :: pos, ilen, in_m, Rq(n_ranks_m), ierr
      
      !-- buffer will receive all messages, but they are ordered by ranks,
      !   not by m.
      pos = 1
      do irank=0,n_ranks_m-1
         in_m = dist_m(irank,0)
         ilen = in_m*(l_max+2) - sum(dist_m(irank,1:in_m))
         if (coord_m == irank) buffer(pos:pos+ilen-1) = Flm_local(1:n_lmP_loc)
         CALL MPI_IBCAST(buffer(pos:pos+ilen-1), ilen, MPI_DOUBLE_COMPLEX, &
                         irank, comm_m, Rq(irank+1), ierr)
         pos = pos + ilen
      end do
      
      CALL MPI_WAITALL(n_ranks_theta, Rq, MPI_STATUSES_IGNORE, ierr)
      
      !-- Reorders the buffer
      l_lm_loc = 1
      do irank=0,n_ranks_m-1
         do j = 1, dist_m(irank,0)
            m = dist_m(irank,j)
            u_lm_loc = l_lm_loc + l_max+1 - m   ! (l_max+2 - m) points for lmP
            l_lm_glb = map_glbl_st%lmP2(m  ,m)
            u_lm_glb = map_glbl_st%lmP2(l_max+1,m)
            FlmP_global(l_lm_glb:u_lm_glb) = buffer(l_lm_loc:u_lm_loc)
            l_lm_loc = u_lm_loc + 1
         end do
      end do
      
   end subroutine gather_FlmP
   
   !----------------------------------------------------------------------------
   subroutine gather_Flm(Flm_local, Flm_global)
      !
      !   IMPORTANT: this function will only work if the all l points are stored
      !   locally and grouped together in "chunks". There is no requirements 
      !   for the m's. 
      !   
      !   For instance, if m1 and m2 are two consecutive m points in a given
      !   rank, then Flm must be stored such that
      !   [(m1:l_max,m1) (m2:l_max,m2)]
      !   are consecutive.
      !   Notice that there is not assumptions about m1 or m2. Also, the l's do
      !   not need to be stored in ascending or descending order. And example 
      !   for l_max=6 and dist_m=(/5,3/) is:
      !   [(5,3) (6,3) (3,3) (4,3) (6,5) (5,5)]
      !   This is perfectly valid grouping.
      !      
      !   Author: Rafael Lago (MPCDF) August 2017
      !
      complex(cp),  intent(in)  :: Flm_local(n_lm_loc)
      complex(cp),  intent(out) :: Flm_global(lm_max)
      
      complex(cp) ::  buffer(lm_max)
      integer :: irank, j, m, l_lm_loc, u_lm_loc, l_lm_glb, u_lm_glb
      integer :: in_m, pos, ilen, Rq(n_ranks_m), ierr
      
      !-- The buffer will receive all messages, but they are ordered by ranks,
      !   not by m.
      
      pos = 1
      do irank=0,n_ranks_m-1
         in_m = dist_m(irank,0)
         ilen = in_m*(l_max+1) - sum(dist_m(irank,1:in_m))
         if (coord_m == irank) buffer(pos:pos+ilen-1) = Flm_local(1:n_lm_loc)
         CALL MPI_IBCAST(buffer(pos:pos+ilen-1), ilen, MPI_DOUBLE_COMPLEX, &
                         irank, comm_m, Rq(irank+1), ierr)
         pos = pos + ilen
      end do
      
      CALL MPI_WAITALL(n_ranks_m, Rq, MPI_STATUSES_IGNORE, ierr)
      
      !-- Reorders the buffer
      l_lm_loc = 1
      do irank=0,n_ranks_m-1
         do j = 1, dist_m(irank, 0)
            m = dist_m(irank, j)
            u_lm_loc = l_lm_loc + l_max - m   ! (l_max+1 - m) points for lm
            l_lm_glb = map_glbl_st%lm2(m , m)
            u_lm_glb = map_glbl_st%lm2(l_max , m)
            Flm_global(l_lm_glb:u_lm_glb) = buffer(l_lm_loc:u_lm_loc)
            l_lm_loc = u_lm_loc + 1
         end do
      end do
      
   end subroutine gather_Flm
   
   !----------------------------------------------------------------------------
   subroutine gather_f(f_local, f_global)
      !
      !   Author: Rafael Lago (MPCDF) August 2017
      !
      real(cp),  intent(inout) :: f_local( n_phi_max, l_theta:u_theta)
      real(cp),  intent(out)   :: f_global(n_phi_max, n_theta_max)
      
      integer :: i, ierr
      integer :: Rq(n_ranks_theta) 
      
      !-- Copies local content to f_global
      f_global = 0.0
      f_global(:,l_theta:u_theta) = f_local(:,l_theta:u_theta)
      
      do i=0,n_ranks_theta-1
         CALL MPI_IBCAST(f_global(:,dist_theta(i,1):dist_theta(i,2)),        &
                         n_phi_max*dist_theta(i,0), MPI_DOUBLE_PRECISION, i, &
                         comm_theta, Rq(i+1), ierr)
      end do
      
      CALL MPI_WAITALL(n_ranks_theta, Rq, MPI_STATUSES_IGNORE, ierr)
      
   end subroutine gather_f
   
   !-- General-purpose transposition. All the customization should happen 
   !   during the creation of the types!
   !
   !   Author: Rafael Lago (MPCDF) January 2018
   !   
   subroutine transpose_ml2r_start(self, container_ml, container_rm, n_fields)
      class(ml2r_transp), intent(inout) :: self
      integer,            intent(in)    :: n_fields
      complex(cp),        intent(inout) :: container_rm(n_lm_loc, l_r:u_r, n_fields)
      complex(cp),        intent(in)    :: container_ml(n_mlo_loc, n_r_max, n_fields)
      
      integer :: i, j, k, icoord_mlo, icoord_r, il_r, ierr
      
      !-- Starts the sends
      do j=1,size(ml2r_dests)
         icoord_mlo = ml2r_dests(j)
         icoord_r = cart%mlo2lmr(icoord_mlo,2)
         il_r = dist_r(icoord_r,1)
         call mpi_isend(container_ml(1,il_r,1), n_fields, ml2r_s_type(j), icoord_mlo, 1, comm_mlo, self%sends(j), ierr)
      end do
      
      !-- Starts the receives
      do j=1,size(ml2r_sources)
         icoord_mlo = ml2r_sources(j)
         call mpi_irecv(container_rm, n_fields, ml2r_r_type(j), icoord_mlo, 1, comm_mlo, self%recvs(j), ierr)
      end do
      
      !-- Copies data which is already local
      do i=1,size(ml2r_loc_dspl,1)
         k = ml2r_loc_dspl(i,1)
         j = ml2r_loc_dspl(i,2)
         container_rm(j,l_r:u_r,1:n_fields) = container_ml(k,l_r:u_r,1:n_fields)
      end do
      
   end subroutine transpose_ml2r_start
   
   !
   !   Author: Rafael Lago (MPCDF) January 2018
   !   
   subroutine transpose_ml2r_wait(self)
      class(ml2r_transp), intent(inout) :: self
      call mpi_waitall(size(self%rq),self%rq,MPI_STATUSES_IGNORE,ierr)
   end subroutine transpose_ml2r_wait
   
   
   !-- General-purpose transposition. All the customization should happen 
   !   during the creation of the types!
   !
   !   Author: Rafael Lago (MPCDF) October 2019
   !   
   subroutine transpose_r2lm_start(self, container_rm, container_ml, n_fields)
      class(r2lm_transp), intent(inout) :: self
      integer,            intent(in)    :: n_fields
      complex(cp),        intent(inout) :: container_ml(n_mlo_loc, n_r_max, n_fields)
      complex(cp),        intent(in)    :: container_rm(n_lm_loc, l_r:u_r, n_fields)
      
      integer :: i, j, k, icoord_mlo, icoord_r, il_r, ierr
      
      !-- Starts the receives
      do j=1,size(ml2r_sources)
         icoord_mlo = ml2r_sources(j)
         call mpi_isend(container_rm, n_fields, ml2r_r_type(j), icoord_mlo, 1, comm_mlo, self%recvs(j), ierr)
      end do
      
      !-- Starts the sends
      do j=1,size(ml2r_dests)
         icoord_mlo = ml2r_dests(j)
         icoord_r = cart%mlo2lmr(icoord_mlo,2)
         il_r = dist_r(icoord_r,1)
         call mpi_irecv(container_ml(1,il_r,1), n_fields, ml2r_s_type(j), icoord_mlo, 1, comm_mlo, self%sends(j), ierr)
      end do
      
      !-- Copies data which is already local
      do i=1,size(ml2r_loc_dspl,1)
         k = ml2r_loc_dspl(i,1)
         j = ml2r_loc_dspl(i,2)
         container_ml(k,l_r:u_r,1:n_fields) = container_rm(j,l_r:u_r,1:n_fields)
      end do
   end subroutine transpose_r2lm_start
   
   !
   !   Author: Rafael Lago (MPCDF) October 2019
   !   
   subroutine transpose_r2lm_wait(self)
      class(r2lm_transp), intent(inout) :: self
      call mpi_waitall(size(self%rq),self%rq,MPI_STATUSES_IGNORE,ierr)
   end subroutine transpose_r2lm_wait
   
   !-- Transposition from (m_loc,θ_glb) to (θ_loc,m_glb).
   !   
   !   
   !   Author: Rafael Lago (MPCDF) August 2017
   !
   !-- TODO this with mpi_type to stride the data
   !
   subroutine transpose_m_theta(f_m_theta, f_theta_m)
      complex(cp), intent(inout) :: f_m_theta(n_m_max, n_theta_loc)
      complex(cp), intent(inout) :: f_theta_m(n_theta_max, n_m_loc)
      
      complex(cp) :: sendbuf(n_m_max * n_theta_loc)
      complex(cp) :: recvbuf(n_m_loc, n_theta_max)
      
      integer :: sendcount(0:n_ranks_m-1)
      integer :: recvcount(0:n_ranks_m-1)
      integer :: senddispl(0:n_ranks_m-1)
      integer :: recvdispl(0:n_ranks_m-1)
      integer :: irank, j, itheta, m, pos
      
      pos = 1
      do irank=0,n_ranks_m-1
         !-- Copy each m which belongs to the irank-th rank into the send buffer
         !   column-wise. That will simplify a lot things later
         !
         !-- TODO check performance of this; implementing this with mpi_type
         !   striding the data will probably be faster
         senddispl(irank) = pos-1
         do itheta=1,n_theta_loc
            do j=1,dist_m(irank,0)
               m = dist_m(irank,j)/minc
               sendbuf(pos) = f_m_theta(m+1,itheta)
               pos = pos + 1
            end do
         end do
         
         sendcount(irank) = pos - senddispl(irank) - 1
         recvdispl(irank) = irank*n_m_loc*dist_theta(irank,0)
         recvcount(irank) =   n_m_loc*dist_theta(irank,0)
      end do
      
      call MPI_ALLTOALLV(sendbuf, sendcount, senddispl, MPI_DOUBLE_COMPLEX, &
                         recvbuf, recvcount, recvdispl, MPI_DOUBLE_COMPLEX, &
                         comm_m, irank)
      f_theta_m = transpose(recvbuf)
      
   end subroutine transpose_m_theta
   
   !-- Transposition from (θ_loc,m_glb) to (m_loc,θ_glb)
   !   
   !   Author: Rafael Lago (MPCDF) August 2017
   !
   !-- TODO this with mpi_type to stride the data
   !
   subroutine transpose_theta_m(f_theta_m, f_m_theta)
      complex(cp), intent(inout) :: f_theta_m(n_theta_max, n_m_loc)
      complex(cp), intent(inout) :: f_m_theta(n_m_max, n_theta_loc)
      
      complex(cp) :: sendbuf(n_m_loc * n_theta_max)
      complex(cp) :: recvbuf(n_theta_loc,  n_m_max)
      
      integer :: sendcount(0:n_ranks_theta-1)
      integer :: recvcount(0:n_ranks_theta-1)
      integer :: senddispl(0:n_ranks_theta-1)
      integer :: recvdispl(0:n_ranks_theta-1)
      integer :: irank, j, pos, n_t, l_t, u_t
      integer :: m_arr(n_ranks_theta*n_m_array) 
      
      recvcount = 0
      pos = 1
      do irank=0,n_ranks_theta-1
         !-- Copy each theta chunk so that the send buffer is contiguous
         !-- TODO check performance of this; implementing this with mpi_type
         !   striding the data will probably be faster
         senddispl(irank) = pos-1
         n_t = dist_theta(irank,0)
         l_t = dist_theta(irank,1)
         u_t = dist_theta(irank,2)
         do j=1, n_m_loc
            sendbuf(pos:pos + n_t - 1) = f_theta_m(l_t:u_t,j)
            pos = pos + n_t
         end do
         
         sendcount(irank) = pos - senddispl(irank) - 1
         recvdispl(irank) = sum(recvcount)
         recvcount(irank) = dist_m(irank,0) * n_t
      end do
      
      call MPI_ALLTOALLV(sendbuf, sendcount, senddispl, MPI_DOUBLE_COMPLEX, &
                         recvbuf, recvcount, recvdispl, MPI_DOUBLE_COMPLEX, &
                         comm_theta, irank)
      
      !-- Now we reorder the receiver buffer. If the m distribution looks like:
      !   rank 0: 0, 4,  8, 12, 16
      !   rank 1: 1, 5,  9, 13
      !   rank 2: 2, 6, 10, 14
      !   rank 3: 3, 7, 11, 15
      !   then the columns of recvbuf are ordered as 0,4,8,12,16,1,5,9,13(...)
      !   and so forth. m_arr will contain this ordering (+1):
      m_arr = reshape(transpose(dist_m(:,1:)), &
                      (/n_ranks_m*n_m_array/))/minc + 1
      j = 1
      do pos = 1, n_ranks_theta*n_m_array
         if (m_arr(pos) < 1) cycle
         f_m_theta(m_arr(pos),:) = recvbuf(:,j)
         j = j + 1
      end do
   end subroutine transpose_theta_m
   

!  /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
!  /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
!
!  From here on, are only temporary functions. They are not optimized. They 
!  should be deleted once the transition is over
!  
!  LM Loop transposes and Gathers and Etc
!
!-------------------------------------------------------------------------------
   subroutine transform_new2old(Fmlo_new, Fmlo_old)
      complex(cp), intent(in) :: Fmlo_new(n_mlo_loc, n_r_max)
      complex(cp), intent(inout) :: Fmlo_old(llm:ulm,   n_r_max)
      
      complex(cp) :: recvbuff(n_r_max)
      integer :: irank, ierr, lm, l, m
      
      do lm=1,lm_max
!          m = map_glbl_st%lm2m(lm)
!          l = map_glbl_st%lm2l(lm)
         m = lo_map%lm2m(lm)
         l = lo_map%lm2l(lm)
         irank = map_mlo%ml2coord(m,l)
         if (irank==coord_mlo) recvbuff = Fmlo_new(map_mlo%ml2i(m,l),:)
         call mpi_bcast(recvbuff, n_r_max, MPI_DOUBLE_COMPLEX, irank, comm_mlo, ierr)

         if (lm>=llm .and. lm<=ulm) Fmlo_old(lm,:) = recvbuff
      end do

   end subroutine transform_new2old
   
!-------------------------------------------------------------------------------
!  A funny thing about the function below is that it is not actually necessary
!  for MagIC itself, but only for the transition. And it is quite the
!  complicated one if I want to have reasonable performance.
!  I'll just use regular-extremely-expensive P2P communications. Do not even
!  try to use this in a large scale simulations. 
   subroutine transform_old2new(Fmlo_old, Fmlo_new)
!-------------------------------------------------------------------------------
      complex(cp), intent(in) :: Fmlo_old(llm:ulm,   n_r_max)
      complex(cp), intent(inout) :: Fmlo_new(n_mlo_loc, n_r_max)
      
      complex(cp) :: recvbuff(n_r_max)
      integer :: old2coord(l_max, l_max)
      integer :: irank, ierr, lm, l, m
      integer :: nLMB_start, nLMB_end
      
      old2coord = -1
      do irank=0,n_ranks_r-1
         nLMB_start = 1+irank*nLMBs_per_rank
         nLMB_end   = min((irank+1)*nLMBs_per_rank,nLMBs)
         do lm=lmStartB(nLMB_start),lmStopB(nLMB_end)
            m = lo_map%lm2m(lm)
            l = lo_map%lm2l(lm)
            
            if (irank==coord_r) recvbuff = Fmlo_old(lm,:)
            call mpi_bcast(recvbuff, n_r_max, MPI_DOUBLE_COMPLEX, irank, comm_r, ierr)
            if (map_mlo%ml2coord(m,l)==coord_mlo) Fmlo_new(map_mlo%ml2i(m,l),:) = recvbuff
         end do
      end do
      
   end subroutine transform_old2new
!--------------------------------------------------------------------------------
!@> Delete me after conversion!
   subroutine test_field(newfield, oldfield, name)
      character(len=*), intent(in) :: name
      complex(cp), intent(in) :: newfield(n_mlo_loc,n_r_max)
      complex(cp), intent(in) :: oldfield(llm:ulm,n_r_max)
      complex(cp) :: test_old(llm:ulm,n_r_max)
      complex(cp) :: test_new(n_mlo_loc,n_r_max)
      real(cp)    :: test_norm, error_threshold
      
      error_threshold = 1e-6 !EPSILON(1.0_cp)
      
      call transform_new2old(newfield, test_old)
      test_norm = ABS(SUM(oldfield - test_old))
      IF (test_norm>error_threshold) print *, "||",name,"n2o|| : ", test_norm

      call transform_old2new(oldfield, test_new)
      test_norm = ABS(SUM(newfield - test_new))
      IF (test_norm>error_threshold) print *, "||",name,"o2n|| : ", test_norm
   end subroutine
   
!-------------------------------------------------------------------------------
   subroutine printArray(inMat, o_prefix, o_fmtString)
      complex(cp), intent(in) :: inMat(:)
      character(*), optional, intent(in) :: o_fmtString
      character(*), optional, intent(in) :: o_prefix
      
      character(:), allocatable :: fmtString, prefix
      integer :: nrow, irow
      
      if (present(o_fmtString)) then
         fmtString = o_fmtString
      else
         fmtString = "(F)"
      end if
      
      if (present(o_prefix)) then
         prefix = o_prefix
      else
         prefix = ""
      end if
      
      nrow = size(inMat,1)
      
      write(*,"(A)", advance="NO") prefix//"["
      do irow=1,nrow-1
         write(*,fmtString, advance="NO") real(inMat(irow))
         write(*,"(A)", advance="NO") "+"
         write(*,fmtString, advance="NO") aimag(inMat(irow))
         write(*,"(A)", advance="NO") "i "
      end do
      write(*,fmtString, advance="NO") real(inMat(irow))
      write(*,"(A)", advance="NO") "+"
      write(*,fmtString, advance="NO") aimag(inMat(irow))
      write(*,"(A)", advance="NO") "i]"//NEW_LINE("A")
      flush(6)
   end subroutine printArray
!-------------------------------------------------------------------------------
   subroutine printArrayReal(inMat, o_fmtString)
      real(cp), intent(in) :: inMat(:)
      character(*), optional, intent(in) :: o_fmtString
      
      character(:), allocatable :: fmtString
      integer :: nrow, irow
      
      if (present(o_fmtString)) then
         fmtString = o_fmtString
      else
         fmtString = "(F)"
      end if
      
      nrow = size(inMat)
      
      write(*,"(A)", advance="NO") "["
      do irow=1,nrow-1
         write(*,fmtString, advance="NO") real(inMat(irow))
         write(*,"(A)", advance="NO") " "
      end do
      write(*,fmtString, advance="NO") real(inMat(irow))
      write(*,"(A)", advance="NO") "]"//NEW_LINE("A")
      flush(6)
   end subroutine printArrayReal
!-------------------------------------------------------------------------------
   subroutine printMatrix(inMat, o_fmtString)
      complex(cp), intent(in) :: inMat(:,:)
      character(*), optional, intent(in) :: o_fmtString
      
      character(:), allocatable :: fmtString
      integer :: nrow, ncol, irow, icol
      
      if (present(o_fmtString)) then
         fmtString = o_fmtString
      else
         fmtString = "(F00.0)"
      end if
      
      nrow = size(inMat,1)
      ncol = size(inMat,2)
      
      do irow=1,nrow
         write(*,"(A)", advance="NO") "["
         do icol=1,ncol-1
            write(*,fmtString, advance="NO") real(inMat(irow,icol))
            write(*,"(A)", advance="NO") "+"
            write(*,fmtString, advance="NO") aimag(inMat(irow,icol))
            write(*,"(A)", advance="NO") "i "
         end do
         write(*,fmtString, advance="NO") real(inMat(irow,ncol))
         write(*,"(A)", advance="NO") "+"
         write(*,fmtString, advance="NO") aimag(inMat(irow,ncol))
         write(*,"(A)", advance="NO") "i]"//NEW_LINE("A")
      end do
      flush(6)
   end subroutine printMatrix
!-------------------------------------------------------------------------------
   subroutine printMatrixInt(inMat, o_fmtString)
      integer, intent(in) :: inMat(:,:)
      character(*), optional, intent(in) :: o_fmtString
      
      character(:), allocatable :: fmtString
      integer :: nrow, ncol, irow, icol
      
      if (present(o_fmtString)) then
         fmtString = o_fmtString
      else
         fmtString = "(I0)"
      end if
      
      nrow = size(inMat,1)
      ncol = size(inMat,2)
      
      do irow=1,nrow
         write(*,"(A)", advance="NO") "["
         do icol=1,ncol-1
            write(*,fmtString, advance="NO") inMat(irow,icol)
            write(*,"(A)", advance="NO") " "
         end do
         write(*,fmtString, advance="NO") inMat(irow,ncol)
         write(*,"(A)", advance="NO") "]"//NEW_LINE("A")
      end do
      flush(6)
   end subroutine printMatrixInt
!-------------------------------------------------------------------------------
   subroutine printMatrixReal(inMat, o_fmtString)
      real(cp), intent(in) :: inMat(:,:)
      character(*), optional, intent(in) :: o_fmtString
      
      character(:), allocatable :: fmtString
      integer :: nrow, ncol, irow, icol
      
      if (present(o_fmtString)) then
         fmtString = o_fmtString
      else
         fmtString = "(F00.0)"
      end if
      
      nrow = size(inMat,1)
      ncol = size(inMat,2)
      
      do irow=1,nrow
         write(*,"(A)", advance="NO") "["
         do icol=1,ncol-1
            write(*,fmtString, advance="NO") inMat(irow,icol)
            write(*,"(A)", advance="NO") " "
         end do
         write(*,fmtString, advance="NO") inMat(irow,ncol)
         write(*,"(A)", advance="NO") "]"//NEW_LINE("A")
      end do
      flush(6)
   end subroutine printMatrixReal
!-------------------------------------------------------------------------------
   subroutine printTriplets(inMat)
      complex(cp), intent(in) :: inMat(:,:)
      
      character(:), allocatable :: fmtString
      integer :: nrow, ncol, irow, icol, l_idx, m_idx, r_idx
      
      fmtString = "(I2)"
      
      nrow = size(inMat,1)
      ncol = size(inMat,2)
      
      do irow=1,nrow
         write(*,"(A)", advance="NO") "[("
         do icol=1,ncol-1
            if (real(inMat(irow,icol))<=0.0) then
               m_idx = -int(real(inMat(irow,icol)))
               l_idx = -1
            else
               m_idx = map_glbl_st%lm2m(int(real(inMat(irow,icol))))
               l_idx = map_glbl_st%lm2l(int(real(inMat(irow,icol))))
            end if
            r_idx = int(aimag(inMat(irow,icol)))
            write(*,fmtString, advance="NO") m_idx
            write(*,"(A)", advance="NO") ","
            write(*,fmtString, advance="NO") l_idx
            write(*,"(A)", advance="NO") ","
            write(*,fmtString, advance="NO") r_idx
            write(*,"(A)", advance="NO") ") | ("
         end do
         if (real(inMat(irow,icol))<=0.0) then
            m_idx = -int(real(inMat(irow,icol)))
            l_idx = -1
         else
            m_idx = map_glbl_st%lm2m(int(real(inMat(irow,ncol))))
            l_idx = map_glbl_st%lm2l(int(real(inMat(irow,ncol))))
         end if
         r_idx = int(aimag(inMat(irow,ncol)))
         write(*,fmtString, advance="NO") m_idx
         write(*,"(A)", advance="NO") ","
         write(*,fmtString, advance="NO") l_idx
         write(*,"(A)", advance="NO") ","
         write(*,fmtString, advance="NO") r_idx
         write(*,"(A)", advance="NO") ")]"//NEW_LINE("A")
      end do
      flush(6)
   end subroutine printTriplets

end module communications
