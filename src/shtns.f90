#include "perflib_preproc.cpp"
module shtns

   use precision_mod, only: cp
   use constants, only: ci
   use communications
   use geometry
   use LMmapping
   use horizontal_data, only: dLh, O_sin_theta_E2, dLh_loc, D_m_loc
   use radial_functions, only: or2
   use parallel_mod
   use fft, only: fft_phi_loc
   use omp_lib

   implicit none

   include "shtns.f"

   private

   public :: initialize_shtns, finalize_shtns,                                    &
   &         scal_to_spat, scal_to_grad_spat, pol_to_grad_spat,             &
   &         torpol_to_spat, pol_to_curlr_spat, torpol_to_curl_spat,        &
   &         torpol_to_dphspat, spat_to_SH, spat_to_sphertor,               &
   &         spat_to_SH_dist, sh_to_spat_dist,                              &
   &         sph_to_spat_dist, torpol_to_spat_dist,                         &
   &         pol_to_curlr_spat_dist, pol_to_grad_spat_dist,                 &
   &         torpol_to_dphspat_dist, torpol_to_curl_spat_dist,              &
   &         spat_to_sphtor_dist

contains

   !----------------------------------------------------------------------------
   subroutine initialize_shtns
      !
      !   Lago: This is important! minc here is sort of the opposite of mres in 
      !   SHTns. In MagIC, the number of m modes stored is n_m_max=m_max/minc+1.
      !   In SHTns, it is simply m_max. Conversely, there are m_max m modes in 
      !   MagIC, though not all of them are effectivelly stored/computed. In 
      !   SHTns, the number of m modes is m_max*mres+1. This makes things very
      !   confusing specially because lm2 and lmP2 arrays (and their relatives)
      !   store m_max m points, and sets those which are not multiple of minc 
      !   to 0.
      !   
      !   In other words, this:
      !   call shtns_lmidx(lm_idx, l_idx, m/minc)
      !   returns the same as 
      !   lm_idx = lm2(l_idx, m)
      !   provided that m is a multiple of minc.
      !

      integer :: norm

      call shtns_verbose(0)
      if ( rank == 0 ) call shtns_verbose(1)

      call shtns_use_threads(0)

      norm = SHT_ORTHONORMAL + SHT_NO_CS_PHASE

      call shtns_set_size(l_max, m_max/minc, minc, norm)
      call shtns_precompute(SHT_GAUSS, SHT_PHI_CONTIGUOUS, &
                            1.e-10_cp, n_theta_max, n_phi_max)
      call shtns_save_cfg(0)

      if ( rank == 0 ) call shtns_verbose(0)

      call shtns_set_size(l_max+1, m_max/minc, minc, norm)
      call shtns_precompute(SHT_QUICK_INIT, SHT_PHI_CONTIGUOUS, &
           &                1.e-10_cp, n_theta_max, n_phi_max)
      call shtns_save_cfg(1)

      call shtns_load_cfg(0)
      
   end subroutine initialize_shtns
   
   !----------------------------------------------------------------------------
   subroutine finalize_shtns
      continue
   end subroutine finalize_shtns
   
   !----------------------------------------------------------------------------
   subroutine scal_to_spat(Slm, fieldc)
      ! transform a spherical harmonic field into grid space
      complex(cp), intent(in) :: Slm(lm_max)
      real(cp), intent(out) :: fieldc(n_phi_max, n_theta_max)

      call shtns_SH_to_spat(Slm, fieldc)

   end subroutine scal_to_spat

   !----------------------------------------------------------------------------
   subroutine scal_to_grad_spat(Slm, gradtc, gradpc)
      !   
      !   Transform a scalar spherical harmonic field into it's gradient
      !   on the grid
      !   
      complex(cp), intent(in) :: Slm(lm_max)
      real(cp), intent(out) :: gradtc(n_phi_max, n_theta_max)
      real(cp), intent(out) :: gradpc(n_phi_max, n_theta_max)

      call shtns_sph_to_spat(Slm, gradtc, gradpc)

   end subroutine scal_to_grad_spat

   !----------------------------------------------------------------------------
   subroutine pol_to_grad_spat(Slm, gradtc, gradpc)

      complex(cp), intent(in) :: Slm(lm_max)
      real(cp), intent(out) :: gradtc(n_phi_max, n_theta_max)
      real(cp), intent(out) :: gradpc(n_phi_max, n_theta_max)

      !-- Local variables
      complex(cp) :: Qlm(lm_max)
      integer :: lm

      do lm = 1, lm_max
         Qlm(lm) = dLh(lm) * Slm(lm)
      end do

      call shtns_sph_to_spat(Qlm, gradtc, gradpc)

   end subroutine pol_to_grad_spat

   !----------------------------------------------------------------------------
   subroutine torpol_to_spat(Wlm, dWlm, Zlm, vrc, vtc, vpc)
      complex(cp), intent(in) :: Wlm(lm_max), dWlm(lm_max), Zlm(lm_max)
      real(cp), intent(out) :: vrc(n_phi_max, n_theta_max)
      real(cp), intent(out) :: vtc(n_phi_max, n_theta_max)
      real(cp), intent(out) :: vpc(n_phi_max, n_theta_max)
      
      !-- Local variables
      complex(cp) :: Qlm(lm_max)
      integer :: lm

      do lm = 1, lm_max
         Qlm(lm) = dLh(lm) * Wlm(lm)
      end do

      call shtns_qst_to_spat(Qlm, dWlm, Zlm, vrc, vtc, vpc)

   end subroutine torpol_to_spat

   !----------------------------------------------------------------------------
   subroutine torpol_to_dphspat(dWlm, Zlm, dvtdp, dvpdp)
      !
      !   Computes horizontal phi derivative of a toroidal/poloidal field
      !
      complex(cp), intent(in) :: dWlm(lm_max), Zlm(lm_max)
      real(cp), intent(out) :: dvtdp(n_phi_max, n_theta_max)
      real(cp), intent(out) :: dvpdp(n_phi_max, n_theta_max)

      !-- Local variables
      complex(cp) :: Slm(lm_max), Tlm(lm_max)
      integer :: lm, it, ip
      real(cp) :: m_rp

      do lm = 1, lm_max
         m_rp = D_m_loc(lm)
         Slm(lm) = ci*m_rp*dWlm(lm)
         Tlm(lm) = ci*m_rp*Zlm(lm)
      end do

      call shtns_sphtor_to_spat(Slm, Tlm, dvtdp, dvpdp)

      do it=1, n_theta_max
         do ip=1, n_phi_max
            dvtdp(ip, it) = dvtdp(ip, it) * O_sin_theta_E2(it)
            dvpdp(ip, it) = dvpdp(ip, it) * O_sin_theta_E2(it)
         end do
      end do

   end subroutine torpol_to_dphspat
   
   !----------------------------------------------------------------------------
   subroutine pol_to_curlr_spat(Qlm, cvrc)
      complex(cp), intent(in) :: Qlm(lm_max)
      real(cp), intent(out) :: cvrc(n_phi_max, n_theta_max)

      !-- Local variables
      complex(cp) :: dQlm(lm_max)
      integer :: lm


      do lm = 1, lm_max
         dQlm(lm) = dLh(lm) * Qlm(lm)
      end do

      call shtns_SH_to_spat(dQlm, cvrc)

   end subroutine pol_to_curlr_spat

   !----------------------------------------------------------------------------
   subroutine torpol_to_curl_spat(Blm, ddBlm, Jlm, dJlm, nR, &
              &                  cvrc, cvtc, cvpc)
      complex(cp), intent(in) :: Blm(lm_max), ddBlm(lm_max)
      complex(cp), intent(in) :: Jlm(lm_max), dJlm(lm_max)
      integer, intent(in) :: nR
      real(cp), intent(out) :: cvrc(n_phi_max, n_theta_max)
      real(cp), intent(out) :: cvtc(n_phi_max, n_theta_max)
      real(cp), intent(out) :: cvpc(n_phi_max, n_theta_max)

      !-- Local variables
      complex(cp) :: Qlm(lm_max), Tlm(lm_max)
      integer :: lm

      do lm = 1, lm_max
         Qlm(lm) = dLh(lm) * Jlm(lm)
         Tlm(lm) = or2(nR) * dLh(lm) * Blm(lm) - ddBlm(lm)
      end do

      call shtns_qst_to_spat(Qlm, dJlm, Tlm, cvrc, cvtc, cvpc)

   end subroutine torpol_to_curl_spat

   !----------------------------------------------------------------------------
   subroutine spat_to_SH(f, fLM)

      real(cp), intent(in) :: f(n_phi_max, n_theta_max)
      complex(cp), intent(out) :: fLM(lmP_max)

      call shtns_load_cfg(1)
      call shtns_spat_to_sh(f, fLM)
      call shtns_load_cfg(0)

   end subroutine spat_to_SH

   !----------------------------------------------------------------------------
   subroutine spat_to_sphertor(f,g,fLM,gLM)

      real(cp), intent(in) :: f(n_phi_max,n_theta_max)
      real(cp), intent(in) :: g(n_phi_max,n_theta_max)
      complex(cp), intent(out) :: fLM(lmP_max)
      complex(cp), intent(out) :: gLM(lmP_max)

      call shtns_load_cfg(1)
      call shtns_spat_to_sphtor(f,g,fLM,gLM)
      call shtns_load_cfg(0)

   end subroutine spat_to_sphertor


   !----------------------------------------------------------------------------
   ! 
   ! 
   !                       Distributed functions
   !   
   !
   !----------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   subroutine spat_to_SH_dist(Vr_loc, Qn_lmP_loc)
      !   
      !   Transform the spherical harmonic coefficients Qlm into its spatial 
      !   representation Vr.
      !  
      !   Author: Rafael Lago, MPCDF, July 2017
      !
      real(cp),     intent(inout) :: Vr_loc(n_phi_max,n_theta_loc)
      complex(cp),  intent(out)   :: Qn_lmP_loc(n_lmP_loc)
      
      !-- Local variables
      complex(cp) ::  QlP_loc(n_theta_max,n_m_loc)
      complex(cp) ::  transpose_loc(n_m_max,n_theta_loc)
      complex(cp) ::  Fr_loc(n_phi_max/2+1,n_theta_loc)
      
      integer :: m, l_lm, u_lm, i
      
      call shtns_load_cfg(1) ! l_max + 1
      
      !-- TODO: The FFT must be performed for an array with the dimensions of 
      !   Fr_loc which may end up paded with zeroes.
      !   Is there any way to tell MKL to perform a "truncated" FFT?
      !-- TODO: Terrible performance here!
      call fft_phi_loc(Vr_loc, Fr_loc, 1)
      transpose_loc(1:n_m_max,1:n_theta_loc) = Fr_loc(1:n_m_max,1:n_theta_loc)
      
      call transpose_m_theta(transpose_loc, QlP_loc)
      
      !-- Now do the Legendre transform using the new function in a loop
      do i = 1, n_m_loc
        m = dist_m(coord_m, i)
        l_lm = map_dist_st%lmP2(m, m)
        u_lm = map_dist_st%lmP2(l_max, m)
        call shtns_spat_to_sh_ml(m/minc,QlP_loc(:,i),Qn_lmP_loc(l_lm:u_lm),l_max+1)
      end do
      
      call shtns_load_cfg(0) ! l_max

   end subroutine spat_to_SH_dist
   
   !----------------------------------------------------------------------------
   subroutine spat_to_sphtor_dist(Vr_loc, Vt_loc, Sn_lmP_loc, Tn_lmP_loc)
      !   
      !   Transform the spherical harmonic coefficients Qlm into its spatial 
      !   representation Vr.
      !  
      !   Author: Rafael Lago, MPCDF, July 2017
      !
      real(cp),     intent(inout) :: Vr_loc(n_phi_max,n_theta_loc)
      real(cp),     intent(inout) :: Vt_loc(n_phi_max,n_theta_loc)
      complex(cp),  intent(out)   :: Sn_lmP_loc(n_lmP_loc)
      complex(cp),  intent(out)   :: Tn_lmP_loc(n_lmP_loc)
      
      !-- Local variables
      complex(cp) ::  SlP_loc(n_theta_max,n_m_loc)
      complex(cp) ::  TlP_loc(n_theta_max,n_m_loc)
      complex(cp) ::  transpose_loc(n_m_max,n_theta_loc)
      complex(cp) ::  F_loc(n_phi_max/2+1,n_theta_loc)
      
      integer :: m, l_lm, u_lm, i
      
      call shtns_load_cfg(1) ! l_max + 1
      
      !-- TODO: The FFT must be performed for an array with the dimensions of 
      !   F_loc which may end up paded with zeroes.
      !   Is there any way to tell MKL to perform a "truncated" FFT?
      !-- TODO: Terrible performance here!
      call fft_phi_loc(Vr_loc, F_loc, 1)
      transpose_loc(1:n_m_max,1:n_theta_loc) = F_loc(1:n_m_max,1:n_theta_loc)
      call transpose_m_theta(transpose_loc, SlP_loc)
      
      call fft_phi_loc(Vt_loc, F_loc, 1)
      transpose_loc(1:n_m_max,1:n_theta_loc) = F_loc(1:n_m_max,1:n_theta_loc)
      call transpose_m_theta(transpose_loc, TlP_loc)
      
      !-- Now do the Legendre transform using the new function in a loop
      do i = 1, n_m_loc
        m = dist_m(coord_m, i)
        l_lm  = map_dist_st%lmP2(m, m)
        u_lm  = map_dist_st%lmP2(l_max, m)
        call shtns_spat_to_sphtor_ml(m/minc,SlP_loc(:,i), TlP_loc(:,i), &
                           Sn_lmP_loc(l_lm:u_lm),Tn_lmP_loc(l_lm:u_lm), l_max+1)
      end do
      
      call shtns_load_cfg(0) ! l_max

   end subroutine spat_to_sphtor_dist

   !----------------------------------------------------------------------------
   subroutine sh_to_spat_dist(Qn_lm_loc, Vr_loc)
      !   
      !   Transform the spherical harmonic coefficients Qlm into its spatial 
      !   representation Vr
      !  
      !   Author: Rafael Lago, MPCDF, July 2017
      !
      complex(cp),  intent(inout)   :: Qn_lm_loc(n_lm_loc)
      real(cp),     intent(out)     :: Vr_loc(n_phi_max, n_theta_loc)
      
      !-- Local variables
      complex(cp) :: Ql_loc(n_theta_max,n_m_loc)
      complex(cp) :: transposed_loc(n_m_max,n_theta_loc)
      complex(cp) :: F_loc(n_phi_max/2+1,n_theta_loc)
      
      
      
      complex(cp) :: Qn_lm_loc_tmp(n_lm_loc)
      complex(cp) :: Ql_loc_tmp(n_theta_max,n_m_loc)
      
      integer :: i, l_lm, u_lm, m
      
      Qn_lm_loc_tmp = Qn_lm_loc
      Ql_loc_tmp    = Ql_loc
      
      PERFON('dbg1')
      !$OMP PARALLEL DO PRIVATE(i,m,l_lm,u_lm) &
      !$OMP &   SHARED(Qn_lm_loc,Ql_loc,n_m_loc,coord_m, dist_m, map_dist_st, l_max, minc) &
      !$OMP &   DEFAULT(None) &
      !$OMP &   SCHEDULE(Dynamic,4)
      do i = 1, n_m_loc
        m = dist_m(coord_m, i)
        l_lm = map_dist_st%lm2(m, m)
        u_lm = map_dist_st%lm2(l_max, m)
        call shtns_sh_to_spat_ml(m/minc, Qn_lm_loc(l_lm:u_lm), Ql_loc(:,i), &
                                 l_max)
      end do
      !$OMP END PARALLEL DO
      PERFOFF
      
      PERFON('dbg2')
      do i = 1, n_m_loc
        m = dist_m(coord_m, i)
        l_lm = map_dist_st%lm2(m, m)
        u_lm = map_dist_st%lm2(l_max, m)
        call shtns_sh_to_spat_ml(m/minc, Qn_lm_loc_tmp(l_lm:u_lm), Ql_loc_tmp(:,i),l_max)
      end do
      PERFOFF
      
      call transpose_theta_m(Ql_loc, transposed_loc)
      
      !-- TODO: The FFT must be performed for an array with the dimensions of 
      !   F_loc which may end up paded with zeroes.
      !   Is there any way to tell MKL to perform a "truncated" FFT?
      !-- TODO: Terrible performance here!
      F_loc = 0.0
      F_loc(1:n_m_max,1:n_theta_loc) = transposed_loc
      
      call fft_phi_loc(Vr_loc, F_loc, -1)
      
   end subroutine sh_to_spat_dist

   !----------------------------------------------------------------------------
   subroutine sph_to_spat_dist(Sn_lm_loc, gVt_loc, gVp_loc)
      !
      !   Transform spheroidal spherical harmonic coefficients Slm to the 
      !   spatial theta and phi components (Vt,Vp), effectively computing the 
      !   gradient of S. This is the distributed counterpart of [scal_to_grad_spat]
      !
      !   Author: Rafael Lago, MPCDF, July 2017
      !
      complex(cp),  intent(inout)   :: Sn_lm_loc(n_lm_loc)
      real(cp),     intent(out)     :: gVt_loc(n_phi_max, n_theta_loc)
      real(cp),     intent(out)     :: gVp_loc(n_phi_max, n_theta_loc)
      
      !-- Local variables
      complex(cp) :: Sl_t(n_theta_max,n_m_loc), Sl_p(n_theta_max,n_m_loc)
      complex(cp) :: tranposed_loc(n_m_max,n_theta_loc)
      complex(cp) :: F_loc(n_phi_max/2+1,n_theta_loc)
      
      integer :: i, l_lm, u_lm, m
      
      
      do i = 1, n_m_loc
        m = dist_m(coord_m, i)
        l_lm  = map_dist_st%lm2(m, m)
        u_lm  = map_dist_st%lm2(l_max, m)
        call shtns_sph_to_spat_ml(m/minc, Sn_lm_loc(l_lm:u_lm), &
                                         Sl_t(:,i), Sl_p(:,i), l_max)
      end do
      
      !-- TODO: write the rest of this function such that both the poloidal and 
      !   the toroidal part can be transposed and Fourier-transformed 
      !   simultanouesly
      call transpose_theta_m(Sl_t, tranposed_loc)
      F_loc = 0.0
      F_loc(1:n_m_max,1:n_theta_loc) = tranposed_loc
      call fft_phi_loc(gVt_loc, F_loc, -1)
      
      
      call transpose_theta_m(Sl_p, tranposed_loc)
      F_loc = 0.0
      F_loc(1:n_m_max,1:n_theta_loc) = tranposed_loc
      call fft_phi_loc(gVp_loc, F_loc, -1)
            
   end subroutine
   
      !----------------------------------------------------------------------------
   subroutine qst_to_spat_dist(Qn_lm_loc, Sn_lm_loc, Tn_lm_loc, &
                                  Vr_loc, Vt_loc, Vp_loc)
      !  
      !   3D vector transform from spherical coordinates to radial-spheroidal-
      !   toroidal spectral components. 
      !
      !   Author: Rafael Lago, MPCDF, July 2017
      !
      complex(cp),  intent(inout) :: Qn_lm_loc(n_lm_loc)
      complex(cp),  intent(inout) :: Sn_lm_loc(n_lm_loc)
      complex(cp),  intent(inout) :: Tn_lm_loc(n_lm_loc)
      real(cp),     intent(out)   :: Vr_loc(n_phi_max, n_theta_loc)
      real(cp),     intent(out)   :: Vt_loc(n_phi_max, n_theta_loc)
      real(cp),     intent(out)   :: Vp_loc(n_phi_max, n_theta_loc)
      
      !-- Local variables
      complex(cp) :: Ql_loc(n_theta_max,n_m_loc)
      complex(cp) :: Sl_loc(n_theta_max,n_m_loc)
      complex(cp) :: Tl_loc(n_theta_max,n_m_loc)
      complex(cp) :: transposed_loc(n_m_max,n_theta_loc)
      complex(cp) :: F_loc(n_phi_max/2+1,n_theta_loc)
      
      integer :: i, l_lm, u_lm, m
      
      do i = 1, n_m_loc
        m = dist_m(coord_m, i)
        l_lm  = map_dist_st%lm2(m, m)
        u_lm  = map_dist_st%lm2(l_max, m)
        call shtns_qst_to_spat_ml(m/minc, Qn_lm_loc(l_lm:u_lm), Sn_lm_loc(l_lm:u_lm), &
                                  Tn_lm_loc(l_lm:u_lm), Ql_loc(:,i), &
                                  Sl_loc(:,i), Tl_loc(:,i), l_max)
      end do
      
      !-- TODO: write the rest of this function such that both the poloidal and 
      !   the toroidal part can be transposed and Fourier-transformed 
      !   simultanouesly
      call transpose_theta_m(Ql_loc, transposed_loc)
      F_loc = 0.0_cp
      F_loc(1:n_m_max,1:n_theta_loc) = transposed_loc
      call fft_phi_loc(Vr_loc, F_loc, -1)

      call transpose_theta_m(Sl_loc, transposed_loc)
      F_loc = 0.0_cp
      F_loc(1:n_m_max,1:n_theta_loc) = transposed_loc
      call fft_phi_loc(Vt_loc, F_loc, -1)

      call transpose_theta_m(Tl_loc, transposed_loc)
      F_loc = 0.0_cp
      F_loc(1:n_m_max,1:n_theta_loc) = transposed_loc
      call fft_phi_loc(Vp_loc, F_loc, -1)
      
   end subroutine qst_to_spat_dist
   
   !----------------------------------------------------------------------------
   subroutine pol_to_curlr_spat_dist(Pn_lm_loc, cVr_loc)
      !   
      !   Computes the spherical harmonic coefficients Qlm throught the 
      !   poloidal scalar Plm, and then transforms it into the curl of its spatial 
      !   representation cVr
      !  
      !   Author: Rafael Lago, MPCDF, August 2017
      !   
      complex(cp),  intent(inout)   :: Pn_lm_loc(n_lm_loc)
      real(cp),     intent(out)     :: cVr_loc(n_phi_max, n_theta_loc)
      
      !-- Local variables
      complex(cp) :: Qn_lm_loc(n_lm_loc)
      complex(cp) :: Ql_loc(n_theta_max,n_m_loc)
      complex(cp) :: transpose_loc(n_m_max,n_theta_loc)
      complex(cp) :: F_loc(n_phi_max/2+1,n_theta_loc)
      
      integer :: i, l_lm, u_lm, m
      
      Qn_lm_loc(1:n_lm_loc) = dLh_loc(1:n_lm_loc) * Pn_lm_loc(1:n_lm_loc)
      
      call sh_to_spat_dist(Qn_lm_loc,cVr_loc)
      
   end subroutine pol_to_curlr_spat_dist

   !----------------------------------------------------------------------------
   subroutine pol_to_grad_spat_dist(Pn_lm_loc, gVt_loc, gVp_loc)
      !   
      !   Computes Qlm (spheroidal spherical harmonic coefficients) from Plm
      !   (poloidal scalar) and then transforms them into the spatial theta and 
      !   phi components (Vt,Vp), effectively computing the gradient of S
      ! 
      !   Author: Rafael Lago, MPCDF, November 2017
      !
      complex(cp),  intent(inout)   :: Pn_lm_loc(lm_max)
      real(cp),     intent(out)     :: gVt_loc(n_phi_max, n_theta_max)
      real(cp),     intent(out)     :: gVp_loc(n_phi_max, n_theta_max)
      
      !-- Local variables
      complex(cp) :: Qn_lm_loc(n_lm_loc)
      
      Qn_lm_loc(1:n_lm_loc) = dLh_loc(1:n_lm_loc) * Pn_lm_loc(1:n_lm_loc)
      
      call sph_to_spat_dist(Qn_lm_loc, gVt_loc, gVp_loc)
            
   end subroutine pol_to_grad_spat_dist
   
   !----------------------------------------------------------------------------
   subroutine torpol_to_spat_dist(Pn_lm_loc, Sn_lm_loc, Tn_lm_loc, &
                                  Vr_loc, Vt_loc, Vp_loc)
      !   
      !   Computes the spherical harmonic coefficients Qlm throught the 
      !   poloidal scalar Plm, then transforms the 3D vector from spherical 
      !   coordinates to radial-spheroidal-toroidal spectral components. 
      !
      !   Author: Rafael Lago, MPCDF, December 2017
      !   
      complex(cp),  intent(inout) :: Pn_lm_loc(n_lm_loc)
      complex(cp),  intent(inout) :: Sn_lm_loc(n_lm_loc)
      complex(cp),  intent(inout) :: Tn_lm_loc(n_lm_loc)
      real(cp),     intent(out)   :: Vr_loc(n_phi_max, n_theta_loc)
      real(cp),     intent(out)   :: Vt_loc(n_phi_max, n_theta_loc)
      real(cp),     intent(out)   :: Vp_loc(n_phi_max, n_theta_loc)
      
      !-- Local variables
      complex(cp) :: Qn_lm_loc(n_lm_loc)
      
      Qn_lm_loc(1:n_lm_loc) = dLH_loc(1:n_lm_loc)*Pn_lm_loc(1:n_lm_loc)
      
      call qst_to_spat_dist(Qn_lm_loc,Sn_lm_loc,Tn_lm_loc,Vr_loc,Vt_loc,Vp_loc)
      
   end subroutine torpol_to_spat_dist
   
   !----------------------------------------------------------------------------
   subroutine torpol_to_dphspat_dist(dWn_lm_loc, Zn_lm_loc, dVt_loc, dVp_loc)
      !   
      !   Transform spheroidal-toroidal spherical harmonic coefficients &
      !   to the derivative of the spatial theta and phi components (Vt,Vp).
      ! 
      !   Author: Rafael Lago, MPCDF, November 2017
      !
      complex(cp),  intent(inout)   :: dWn_lm_loc(n_lm_loc)
      complex(cp),  intent(inout)   :: Zn_lm_loc(n_lm_loc)
      real(cp),     intent(out)     :: dVt_loc(n_phi_max, n_theta_max)
      real(cp),     intent(out)     :: dVp_loc(n_phi_max, n_theta_max)
      
      !-- Local variables
      complex(cp) :: Sn_lm_loc(n_lm_loc)
      complex(cp) :: Tn_lm_loc(n_lm_loc)
      complex(cp) :: Sl_loc(n_theta_max,n_m_loc)
      complex(cp) :: Tl_loc(n_theta_max,n_m_loc)
      complex(cp) :: transpose_loc(n_m_max,n_theta_loc)
      complex(cp) :: F_loc(n_phi_max/2+1,n_theta_loc)
      
      integer  :: i, j, l_lm, u_lm, m
      complex(cp) :: m_cp
      
      !-- TODO: do we need those temporary arrays?
      !   This is D_m_loc is just m...
      do i = 1, n_lm_loc
         m_cp = ci*D_m_loc(i)
         Sn_lm_loc(i) = m_cp*dWn_lm_loc(i)
         Tn_lm_loc(i) = m_cp* Zn_lm_loc(i)
      end do
      
      do i = 1, n_m_loc
        m = dist_m(coord_m, i)
        l_lm = map_dist_st%lm2(m, m)
        u_lm = map_dist_st%lm2(l_max, m)
        call shtns_sphtor_to_spat_ml(m/minc, Sn_lm_loc(l_lm:u_lm), Tn_lm_loc(l_lm:u_lm), &
                                     Sl_loc(:,i), Tl_loc(:,i), l_max)
      end do
      
      !-- TODO: write the rest of this function such that both the poloidal and 
      !   the toroidal part can be transposed and Fourier-transformed simultanouesly
      call transpose_theta_m(Sl_loc, transpose_loc)
      F_loc = 0.0_cp
      F_loc(1:n_m_max,1:n_theta_loc) = transpose_loc
      call fft_phi_loc(dVt_loc, F_loc, -1)
      do j=1, n_theta_loc
         do i=1, n_phi_max
            dVt_loc(i, j) = dVt_loc(i, j) * O_sin_theta_E2(l_theta+j-1)
         end do
      end do

      call transpose_theta_m(Tl_loc, transpose_loc)
      F_loc = 0.0_cp
      F_loc(1:n_m_max,1:n_theta_loc) = transpose_loc
      call fft_phi_loc(dVp_loc, F_loc, -1)
      
      do j=1, n_theta_loc
         do i=1, n_phi_max
            dVp_loc(i, j) = dVp_loc(i, j) * O_sin_theta_E2(l_theta+j-1)
         end do
      end do
            
   end subroutine torpol_to_dphspat_dist
   
   !----------------------------------------------------------------------------
   subroutine torpol_to_curl_spat_dist(Bn_lm_loc, ddBn_lm_loc, Jn_lm_loc, Sn_lm_loc, &
                                       nR, cVr_loc, cVt_loc, cVp_loc)
      !   3D vector transform from spherical coordinates to the curl of the 
      !   radial-spheroidal-toroidal spectral components. 
      ! 
      !   Author: Rafael Lago, MPCDF, November 2017
      !   
      !-- TODO: Review the notation here
      
      complex(cp),  intent(inout) :: Bn_lm_loc(n_lm_loc)
      complex(cp),  intent(inout) :: ddBn_lm_loc(n_lm_loc)
      complex(cp),  intent(inout) :: Jn_lm_loc(n_lm_loc)
      complex(cp),  intent(inout) :: Sn_lm_loc(n_lm_loc)
      integer,      intent(in)    :: nR
      real(cp),     intent(out)   :: cVr_loc(n_phi_max, n_theta_max)
      real(cp),     intent(out)   :: cVt_loc(n_phi_max, n_theta_max)
      real(cp),     intent(out)   :: cVp_loc(n_phi_max, n_theta_max)
      
      !-- Local variables
      complex(cp) :: Qn_lm_loc(n_lm_loc)
      complex(cp) :: Tn_lm_loc(n_lm_loc)
      
      Qn_lm_loc(1:n_lm_loc) = dLH_loc(1:n_lm_loc)*Jn_lm_loc(1:n_lm_loc)
      Tn_lm_loc(1:n_lm_loc) = or2(nR) * dLH_loc(1:n_lm_loc) * Bn_lm_loc(1:n_lm_loc) &
                        - ddBn_lm_loc(1:n_lm_loc) 
      
      call qst_to_spat_dist(Qn_lm_loc,Sn_lm_loc,Tn_lm_loc,cVr_loc,cVt_loc,cVp_loc)
      
   end subroutine torpol_to_curl_spat_dist

end module shtns
