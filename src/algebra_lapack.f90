module algebra

   use precision_mod, only: cp

   implicit none

   private

   public :: prepare_mat, solve_mat
   
   ! There might be a large difference between the solutions when using 
   ! the multiple rhs in the lapack call (for some reason). It is hard 
   ! to debug the code sometimes because of these differences. Set this flag 
   ! for changing which method to use
   ! 0: loops over each RHS individually (slowest, more precise)
   ! 1: solves once for the real and once for the imaginary part (original version)
   ! anything else: solves both real and imaginary at once in one 
   !     block (optimized, different precision)
   ! 
   integer, parameter :: multiple_rhs_method = 2

   interface solve_mat
      module procedure solve_mat_real_rhs
      module procedure solve_mat_complex_rhs
      module procedure solve_mat_complex_rhs_multi
   end interface solve_mat

contains

   subroutine solve_mat_complex_rhs(a,len_a,n,pivot,rhs)
      !
      !  This routine does the backward substitution into a lu-decomposed real 
      !  matrix a (to solve a * x = bc1) were bc1 is the right hand side  
      !  vector. On return x is stored in bc1.                            
      !                                                                     

      !-- Input variables:
      integer,  intent(in) :: n          ! dimension of problem
      integer,  intent(in) :: len_a      ! first dim of a
      integer,  intent(in) :: pivot(n)   ! pivot pointer of legth n
      real(cp), intent(in) :: a(len_a,n) ! real n X n matrix

      !-- Output variables
      complex(cp), intent(inout) :: rhs(n) ! on input RHS of problem

      !-- Local variables:
      real(cp) :: tmpr(n), tmpi(n)
      integer :: info, i

      do i=1,n
         tmpr(i) = real(rhs(i))
         tmpi(i) = aimag(rhs(i))
      end do

#if (DEFAULT_PRECISION==sngl)
      call sgetrs('N',n,1,a,len_a,pivot,tmpr,n,info)
      call sgetrs('N',n,1,a,len_a,pivot,tmpi,n,info)
#elif (DEFAULT_PRECISION==dble)
      call dgetrs('N',n,1,a,len_a,pivot,tmpr,n,info)
      call dgetrs('N',n,1,a,len_a,pivot,tmpi,n,info)
#endif

      do i=1,n
         rhs(i)=cmplx(tmpr(i),tmpi(i),kind=cp)
      end do

   end subroutine solve_mat_complex_rhs
!-----------------------------------------------------------------------------
   subroutine solve_mat_complex_rhs_multi(a,len_a,n,pivot,rhs,nRHSs)
      !
      !  This routine does the backward substitution into a lu-decomposed real
      !  matrix a (to solve a * x = bc ) simultaneously for nRHSs complex 
      !  vectors bc. On return the results are stored in the bc.                  
      !

      !-- Input variables:
      integer,  intent(in) :: n           ! dimension of problem
      integer,  intent(in) :: len_a       ! leading dimension of a
      integer,  intent(in) :: pivot(n)       ! pivot pointer of length n
      real(cp), intent(in) :: a(len_a,n)  ! real n X n matrix
      integer,  intent(in) :: nRHSs       ! number of right-hand sides

      complex(cp), intent(inout) :: rhs(:,:) ! on input RHS of problem

      !-- Local variables:
      real(cp), allocatable :: tmpr(:,:), tmpi(:,:)
      real(cp) :: norm_b, norm_diff
      integer :: info, i, j
      
      ! There is a hack here to get deterministic solves - read the 
      ! description of the multiple_rhs_method flag at the top of this 
      ! module - Lago
      ! 
      !   Precise variant:
      ! ------------------------------------------------------------------
      if (multiple_rhs_method==0) then
         allocate(tmpr(n,1), tmpi(n,1))
         do j=1,nRHSs
            tmpr(:,1) = real(rhs(:,j))
            tmpi(:,1) = aimag(rhs(:,j))

#if (DEFAULT_PRECISION==sngl)
            call sgetrs('N',n,1,a(1:n,1:n),n,pivot(1:n),tmpr(1:n,1),n,info)
            call sgetrs('N',n,1,a(1:n,1:n),n,pivot(1:n),tmpi(1:n,1),n,info)
#elif (DEFAULT_PRECISION==dble)
            call dgetrs('N',n,1,a(1:n,1:n),n,pivot(1:n),tmpr(1:n,1),n,info)
            call dgetrs('N',n,1,a(1:n,1:n),n,pivot(1:n),tmpi(1:n,1),n,info)
#endif

            rhs(:,j)=cmplx(tmpr(:,1),tmpi(:,1),kind=cp)
         end do
      
      !   Original Variant:
      ! ------------------------------------------------------------------
      else if (multiple_rhs_method==1) then
         allocate(tmpr(n,nRHSs), tmpi(n,nRHSs))
         do j=1,nRHSs
            do i=1,n
               tmpr(i,j) = real(rhs(i,j))
               tmpi(i,j) = aimag(rhs(i,j))
            end do
         end do

#if (DEFAULT_PRECISION==sngl)
         call sgetrs('N',n,nRHSs,a(1:n,1:n),n,pivot(1:n),tmpr(1:n,:),n,info)
         call sgetrs('N',n,nRHSs,a(1:n,1:n),n,pivot(1:n),tmpi(1:n,:),n,info)
#elif (DEFAULT_PRECISION==dble)
         call dgetrs('N',n,nRHSs,a(1:n,1:n),n,pivot(1:n),tmpr(1:n,:),n,info)
         call dgetrs('N',n,nRHSs,a(1:n,1:n),n,pivot(1:n),tmpi(1:n,:),n,info)
#endif

         do j=1,nRHSs
            do i=1,n
               rhs(i,j)=cmplx(tmpr(i,j),tmpi(i,j),kind=cp)
            end do
         end do
         deallocate(tmpr,tmpi)
      
      !   Optimized Variant:
      ! ------------------------------------------------------------------
      else
         allocate(tmpr(n,2*nRHSs))
         tmpr(1:n,1:nRHSs) = real(rhs(1:n,1:nRHSs))
         tmpr(1:n,nRHSs+1:2*nRHSs) = aimag(rhs(1:n,1:nRHSs))
         
#if (DEFAULT_PRECISION==sngl)
         call sgetrs('N',n,2*nRHSs,a(1:n,1:n),n,pivot(1:n),tmpr(1:n,1:2*nRHSs),n,info)
#elif (DEFAULT_PRECISION==dble)
         call dgetrs('N',n,2*nRHSs,a(1:n,1:n),n,pivot(1:n),tmpr(1:n,1:2*nRHSs),n,info)
#endif

         rhs = cmplx(tmpr(1:n,1:nRHSs),tmpr(1:n,nRHSs+1:2*nRHSs))
         deallocate(tmpr)
      end if       
         

   end subroutine solve_mat_complex_rhs_multi
!-----------------------------------------------------------------------------
   subroutine solve_mat_real_rhs(a,len_a,n,pivot,rhs)
      !
      !     like the linpack routine
      !     backward substitution of vector b into lu-decomposed matrix a
      !     to solve  a * x = b for a single real vector b
      !
      !     sub prepare_mat must be called once first to initialize a and pivot
      !

      !-- Input variables:
      integer,  intent(in) :: n         ! dim of problem
      integer,  intent(in) :: len_a     ! first dim of a
      integer,  intent(in) :: pivot(n)  ! pivot information
      real(cp), intent(in) :: a(len_a,n)

      !-- Output: solution stored in rhs(n)
      real(cp), intent(inout) :: rhs(n)
      integer :: info

#if (DEFAULT_PRECISION==sngl)
      call sgetrs('N',n,1,a,len_a,pivot,rhs,n,info)
#elif (DEFAULT_PRECISION==dble)
      call dgetrs('N',n,1,a,len_a,pivot,rhs,n,info)
#endif

   end subroutine solve_mat_real_rhs
!-----------------------------------------------------------------------------
   subroutine prepare_mat(a,len_a,n,pivot,info)
      !
      !     like the linpack routine
      !
      !     lu decomposes the real matrix a(n,n) via gaussian elimination
      !

      !-- Input variables:
      integer,  intent(in) :: len_a,n
      real(cp), intent(inout) :: a(len_a,n)

      !-- Output variables:
      integer, intent(out) :: pivot(n)   ! pivoting information
      integer, intent(out) :: info

#if (DEFAULT_PRECISION==sngl)
      call sgetrf(n,n,a(1:n,1:n),n,pivot(1:n),info)
#elif (DEFAULT_PRECISION==dble)
      call dgetrf(n,n,a(1:n,1:n),n,pivot(1:n),info)
#endif

   end subroutine prepare_mat
!-----------------------------------------------------------------------------
end module algebra
