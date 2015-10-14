module shtns

   use precision_mod, only: cp
   use truncation, only: ncp, m_max, l_max, n_theta_max, n_phi_max, &
                         minc, lm_max, n_m_max
   use blocking, only: nfs
   use horizontal_data, only: dLh, gauss, theta_ord, D_m, O_sin_theta_E2
   use radial_functions, only: r
   use parallel_mod

   implicit none

   include "shtns.f"

   private

   public :: init_shtns, scal_to_spat, scal_to_grad_spat, pol_to_grad_spat, & 
             torpol_to_spat, pol_to_curlr_spat, torpol_to_curl_spat,        &
             torpol_to_dphspat

contains

   subroutine init_shtns()

      integer :: it, ip
      integer :: nlm
      integer :: norm
      integer :: nThreads

      if ( rank == 0 ) then
         call shtns_verbose(1)
      end if

      call shtns_use_threads(0)

      norm = SHT_ORTHONORMAL + SHT_NO_CS_PHASE

      call shtns_set_size(l_max, m_max/minc, minc, norm)
      call shtns_calc_nlm(nlm, l_max, m_max/minc, minc)
      call shtns_precompute(SHT_GAUSS, SHT_PHI_CONTIGUOUS, 1.e-12_cp, &
                            n_theta_max, n_phi_max)
      if (lm_max /= nlm) then
         print*, "error: nlm /= lm_max", nlm, lm_max
#ifdef WITH_MPI
         call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
#else
         stop "error: nlm /= lm_max"
#endif
      end if

   end subroutine
!------------------------------------------------------------------------------
   subroutine scal_to_spat(Slm, fieldc)
      ! transform a spherical harmonic field into grid space
      complex(cp), intent(in) :: Slm(:)
      real(cp), intent(out) :: fieldc(:, :)

      call shtns_SH_to_spat(Slm, fieldc)

   end subroutine scal_to_spat
!------------------------------------------------------------------------------
   subroutine scal_to_grad_spat(Slm, gradtc, gradpc)
      ! transform a scalar spherical harmonic field into it's gradient
      ! on the grid
      complex(cp), intent(in) :: Slm(:)
      real(cp), intent(out) :: gradtc(:, :), gradpc(:, :)

      call shtns_sph_to_spat(Slm, gradtc, gradpc)

   end subroutine scal_to_grad_spat
!------------------------------------------------------------------------------
   subroutine pol_to_grad_spat(Slm, gradtc, gradpc)

      complex(cp), intent(in) :: Slm(:)
      real(cp), intent(out) :: gradtc(:, :), gradpc(:, :)

      ! local
      complex(cp), allocatable :: Qlm(:)
      integer :: ip, it, lm

      allocate(Qlm(lm_max))

      do lm = 1, lm_max
         Qlm(lm) = dLh(lm) * Slm(lm)
      end do

      call shtns_sph_to_spat(Qlm, gradtc, gradpc)
      
      deallocate(Qlm)

   end subroutine pol_to_grad_spat
!------------------------------------------------------------------------------
   subroutine torpol_to_spat(Wlm, dWlm, Zlm, vrc, vtc, vpc)
      complex(cp), intent(in) :: Wlm(:), dWlm(:), Zlm(:)
      real(cp), intent(out) :: vrc(:, :), vtc(:, :), vpc(:, :)

      ! local
      complex(cp), allocatable :: Qlm(:)
      integer :: lm

      allocate(Qlm(lm_max))

      do lm = 1, lm_max
         Qlm(lm) = dLh(lm) * Wlm(lm)
      end do

      call shtns_qst_to_spat(Qlm, dWlm, Zlm, vrc, vtc, vpc)

      deallocate(Qlm)
   end subroutine torpol_to_spat
!------------------------------------------------------------------------------
   subroutine torpol_to_dphspat(dWlm, Zlm, dvtdp, dvpdp)
      !
      ! Computes horizontal phi derivative of a toroidal/poloidal field
      !
      complex(cp), intent(in) :: dWlm(:), Zlm(:)
      real(cp), intent(out) :: dvtdp(:, :), dvpdp(:, :)

      ! local
      complex(cp), allocatable :: Slm(:), Tlm(:)
      complex(cp) :: i = cmplx(0.0, 1.0)
      integer :: lm, it, ip
      real(cp) :: m

      allocate(Slm(lm_max), Tlm(lm_max))
      do lm = 1, lm_max
         m = D_m(lm)
         Slm(lm) = i*m*dWlm(lm)
         Tlm(lm) = i*m*Zlm(lm)
      end do

      call shtns_sphtor_to_spat(Slm, Tlm, dvtdp, dvpdp)

      do it=1, n_theta_max
         do ip=1, n_phi_max
            dvtdp(ip, it) = dvtdp(ip, it) * O_sin_theta_E2(it)
            dvpdp(ip, it) = dvpdp(ip, it) * O_sin_theta_E2(it)
         end do
      end do

      deallocate(Slm, Tlm)
   end subroutine torpol_to_dphspat
!------------------------------------------------------------------------------
   subroutine pol_to_curlr_spat(Qlm, cvrc)
      complex(cp), intent(in) :: Qlm(:)
      real(cp), intent(out) :: cvrc(:, :)

      ! local
      complex(cp), allocatable :: dQlm(:)
      integer :: lm

      allocate(dQlm(lm_max))

      do lm = 1, lm_max
         dQlm(lm) = dLh(lm) * Qlm(lm)
      end do

      call shtns_SH_to_spat(dQlm, cvrc)

      deallocate(dQlm)
   end subroutine pol_to_curlr_spat
!------------------------------------------------------------------------------
   subroutine torpol_to_curl_spat(Blm, dBlm, ddBlm, Jlm, dJlm, nR, &
                                 cvrc, cvtc, cvpc)
      complex(cp), intent(in) :: Blm(:), dBlm(:), ddBlm(:)
      complex(cp), intent(in) :: Jlm(:), dJlm(:)
      integer, intent(in) :: nR
      real(cp), intent(out) :: cvrc(:, :), cvtc(:, :), cvpc(:, :)

      ! local
      complex(cp), allocatable :: Qlm(:), Tlm(:)
      integer :: it, ip, lm

      allocate(Qlm(lm_max), Tlm(lm_max))

      do lm = 1, lm_max
         Qlm(lm) = dLh(lm) * Jlm(lm)
         Tlm(lm) = 1/r(nR)**2 * dLh(lm) * Blm(lm) - ddBlm(lm)
      end do

      call shtns_qst_to_spat(Qlm, dJlm, Tlm, cvrc, cvtc, cvpc)

      deallocate(Qlm, Tlm)
   end subroutine torpol_to_curl_spat
!------------------------------------------------------------------------------
end module shtns
