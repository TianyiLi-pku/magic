#include "perflib_preproc.cpp"
#include "mkl_dfti.f90"

module fft

   use precision_mod
   use constants, only: one
   use geometry, only: nrp, n_phi_max, n_theta_max, m_max, n_theta_loc
   use blocking, only: nfs
   use mkl_dfti
 
   implicit none
   
   private
 
   type(DFTI_DESCRIPTOR), pointer :: c2r_handle, r2c_handle
   type(DFTI_DESCRIPTOR), pointer :: phi2m_dhandle, m2phi_dhandle
   integer :: status
 
   public :: fft_thetab, init_fft, fft_to_real, finalize_fft
   public :: finalize_fft_phi, fft_phi_loc
   public :: initialize_fft_phi

contains

   subroutine init_fft(number_of_points)
      !
      !   MKL equivalent of init_fft
      !

      !-- Input variable
      integer, intent(in) :: number_of_points ! number of points

      ! Fourier transformation complex->REAL with MKL DFTI interface
      ! init FFT
      status = DftiCreateDescriptor( c2r_handle, DFTI_DOUBLE, DFTI_REAL, &
                                     1, number_of_points )
      status = DftiSetValue( c2r_handle, DFTI_NUMBER_OF_TRANSFORMS, nfs )
      status = DftiSetValue( c2r_handle, DFTI_INPUT_DISTANCE, nrp )
      status = DftiSetValue( c2r_handle, DFTI_OUTPUT_DISTANCE, nrp )
      !status = DftiSetValue( c2r_handle, DFTI_CONJUGATE_EVEN_STORAGE, &
      !                       DFTI_COMPLEX_COMPLEX )
      status = DftiSetValue( c2r_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE )
      ! This is supposed to fix the issue for ifort < 14 but it increases the
      ! walltime
      !status = DftiSetValue( c2r_handle, DFTI_NUMBER_OF_USER_THREADS, nThreads)
      status = DftiCommitDescriptor( c2r_handle )
  
      ! Fourier transformation REAL->complex with MKL DFTI interface
      ! init FFT
      status = DftiCreateDescriptor( r2c_handle, DFTI_DOUBLE, DFTI_REAL, &
                                     1, number_of_points )
      status = DftiSetValue( r2c_handle, DFTI_NUMBER_OF_TRANSFORMS, nfs )
      status = DftiSetValue( r2c_handle, DFTI_INPUT_DISTANCE, nrp )
      status = DftiSetValue( r2c_handle, DFTI_OUTPUT_DISTANCE, nrp )
      !status = DftiSetValue( r2c_handle, DFTI_CONJUGATE_EVEN_STORAGE, &

      status = DftiSetValue( r2c_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE )
      status = DftiSetValue( r2c_handle, DFTI_FORWARD_SCALE, &
                             one/real(number_of_points,cp) )
      !status = DftiSetValue( r2c_handle, DFTI_NUMBER_OF_USER_THREADS, nThreads)
      status = DftiCommitDescriptor( r2c_handle )
      
   end subroutine init_fft
   
   !------------------------------------------------------------------------------
   subroutine finalize_fft

      status = DftiFreeDescriptor( c2r_handle )
      status = DftiFreeDescriptor( r2c_handle )

   end subroutine finalize_fft

   !----------------------------------------------------------------------------
   subroutine initialize_fft_phi
      !
      !   The phi2m_dhandle is the handle for the real transform.
      !   It will store only Nφ/2+1 modes, because the others are the conjugate 
      !   of these. This is for the distributed case as well (and thus, will 
      !   repeat n_theta_loc times)
      !   
      !   Author: Rafael Lago, MPCDF, July 2017
      !   
      integer :: status
      
      status = DftiCreateDescriptor( phi2m_dhandle, DFTI_DOUBLE, DFTI_REAL, 1, n_phi_max )
      status = DftiSetValue( phi2m_dhandle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
      status = DftiSetValue( phi2m_dhandle, DFTI_NUMBER_OF_TRANSFORMS, n_theta_loc)
      status = DftiSetValue( phi2m_dhandle, DFTI_INPUT_DISTANCE, n_phi_max )
      status = DftiSetValue( phi2m_dhandle, DFTI_OUTPUT_DISTANCE, n_phi_max/2+1 )
      status = DftiSetValue( phi2m_dhandle, DFTI_PLACEMENT, DFTI_NOT_INPLACE )
      status = DftiSetValue( phi2m_dhandle, DFTI_FORWARD_SCALE, 1.0_cp/real(n_phi_max,cp) )
      status = DftiCommitDescriptor( phi2m_dhandle )
      
      status = DftiCreateDescriptor( m2phi_dhandle, DFTI_DOUBLE, DFTI_REAL, 1, n_phi_max )
      status = DftiSetValue( m2phi_dhandle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
      status = DftiSetValue( m2phi_dhandle, DFTI_NUMBER_OF_TRANSFORMS, n_theta_loc)
      status = DftiSetValue( m2phi_dhandle, DFTI_INPUT_DISTANCE, n_phi_max/2+1 )
      status = DftiSetValue( m2phi_dhandle, DFTI_OUTPUT_DISTANCE, n_phi_max )
      status = DftiSetValue( m2phi_dhandle, DFTI_PLACEMENT, DFTI_NOT_INPLACE )
      status = DftiSetValue( m2phi_dhandle, DFTI_FORWARD_SCALE, 1.0_cp/real(n_phi_max,cp) )
      status = DftiCommitDescriptor( m2phi_dhandle )
    
  end subroutine initialize_fft_phi

   !----------------------------------------------------------------------------
   subroutine finalize_fft_phi
      !
      !   Author: Rafael Lago, MPCDF, July 2017
      !   
   
      integer :: status
      
      status = DftiFreeDescriptor(phi2m_dhandle)
      status = DftiFreeDescriptor(m2phi_dhandle)
   end subroutine
  
   !----------------------------------------------------------------------------
   subroutine fft_phi_loc(f, g, dir)
      !   FFT meant to be used for distributed vectors. It does not perform 
      !   the transposition. That can be achieved with a call to
      !   transpose_m_theta or transpose_theta_m.
      !   
      !   Forward:  converts f(φ,θ_loc) to g(m,θ_loc) 
      !   Backward: converts g(m,θ_loc) to f(φ,θ_loc)
      !   
      !   OBS: g has (n_phi_max/2+1)*n_theta_loc points, not m_max*n_theta_loc points.
      !       Only the first m_max rows are required by MagIC.
      !  
      !   Author: Rafael Lago, MPCDF, November 2017
      ! 
      real(cp),    intent(inout)  :: f(n_phi_max*n_theta_loc)
      complex(cp), intent(inout)  :: g((n_phi_max/2+1)*n_theta_loc)
      integer,     intent(in)     :: dir
      integer :: status
      
      PERFON('fftphi')      
      if (dir == 1) then
         status = DftiComputeForward( phi2m_dhandle, f(:), g(:) )
      else if (dir == -1) then
         status = DftiComputeBackward( m2phi_dhandle, g(:), f(:) )
      else
         print *, "Unknown direction in fft_phi_loc: ", dir
      end if
      PERFOFF
    
  end subroutine fft_phi_loc

  
   !----------------------------------------------------------------------------
   subroutine fft_thetab(f,dir)
      !  
      !   Author: ???
      ! 
      real(cp), intent(inout) :: f(nrp,nfs)
      integer,  intent(in) :: dir            ! back or forth transform

      if (dir == -1) then
         ! run FFT
         status = DftiComputeForward( r2c_handle, f(:,1) )
      else if (dir == 1) then
         ! we want a backward transform, and the real array f 
         ! is to be interpreted as complex array

         ! run FFT
         status = DftiComputeBackward( c2r_handle, f(:,1) )
      end if

   end subroutine fft_thetab
   !----------------------------------------------------------------------------
   subroutine fft_to_real(f,ld_f,nrep)
      !  
      !   Author: ???
      ! 
      integer,  intent(in) :: ld_f,nrep
      real(cp), intent(inout) :: f(ld_f,nrep)

      type(DFTI_DESCRIPTOR), pointer :: local_c2r_handle
      real(cp), allocatable :: work_array(:, :)
      allocate(work_array(ld_f+2, nrep))

      ! Fourier transformation complex->REAL with MKL DFTI interface
      ! init FFT
      status = DftiCreateDescriptor( local_c2r_handle, DFTI_DOUBLE, &
                                     DFTI_REAL, 1, n_phi_max )
      status = DftiSetValue( local_c2r_handle, DFTI_NUMBER_OF_TRANSFORMS, nrep )
      status = DftiSetValue( local_c2r_handle, DFTI_OUTPUT_DISTANCE, ld_f )
      status = DftiSetValue( local_c2r_handle, DFTI_INPUT_DISTANCE, ld_f+2 )
      status = DftiSetValue( local_c2r_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE )
      status = DftiCommitDescriptor( local_c2r_handle )

      ! run FFT
      work_array(1:n_phi_max, :) = f(1:n_phi_max, :)
      work_array(n_phi_max+1:n_phi_max+2, :) = 0.0_cp
      status = DftiComputeBackward( local_c2r_handle, work_array(:, 1), f(:,1) )
      deallocate(work_array)

      status = DftiFreeDescriptor( local_c2r_handle )

   end subroutine fft_to_real
end module fft
