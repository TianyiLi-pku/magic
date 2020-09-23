module nonlinear_bcs

   use iso_fortran_env, only: output_unit
   use precision_mod
   use LMmapping, only: map_dist_st
   use truncation, only: n_phi_max, l_axi, l_max, n_r_cmb, n_r_icb, &
       &                 nThetaStart, nThetaStop, n_lmP_loc, n_lm_loc
   use radial_functions, only: r_cmb, r_icb, rho0
   use physical_parameters, only: sigma_ratio, conductance_ma, prmag, oek
   use horizontal_data, only: dTheta1S_loc, dTheta1A_loc, dPhi_loc, O_sin_theta, &
       &                      dLh_loc, sn2, cosTheta
   use constants, only: two
   use sht, only: scal_to_SH
   use useful, only: abortRun
   use sht, only: scal_to_SH

   implicit none

   private

   public :: get_br_v_bcs, get_b_nl_bcs, v_rigid_boundary

contains

   subroutine get_br_v_bcs(br,vt,vp,omega,O_r_E_2,O_rho,br_vt_lm,br_vp_lm)
      !
      !  Purpose of this subroutine is to calculate the nonlinear term
      !  of the magnetic boundary condition for a conducting mantle or
      !  inner core in space (r,lm).
      !  Calculation is performed for the theta block:
      !
      !  .. code-block:: fortran
      !
      !                  n_theta_min<=n_theta<=n_theta_min+n_theta_block-1
      !
      !  On input ``br``, ``vt`` and ``vp`` are given on all phi points and
      !  thetas in the specific block.
      !  On output the contribution of these grid points to all
      !  degree and orders is stored in ``br_vt_lm`` and ``br_vp_lm``.
      !  Output is ``[r/sin(theta)*Br*U]=[(0,br_vt_lm,br_vp_lm)]``
      !

      !-- input:
      real(cp), intent(in) :: br(nThetaStart:nThetaStop,n_phi_max) ! :math:`r^2 B_r`
      real(cp), intent(in) :: vt(nThetaStart:nThetaStop,n_phi_max) ! :math:`r \sin\theta u_\theta`
      real(cp), intent(in) :: vp(nThetaStart:nThetaStop,n_phi_max) ! :math:`r \sin\theta u_\phi`
      real(cp), intent(in) :: omega          ! rotation rate of mantle or IC
      real(cp), intent(in) :: O_r_E_2        ! :math:`1/r^2`
      real(cp), intent(in) :: O_rho          ! :math:`1/\tilde{\rho}` (anelastic)

      !-- Output variables:
      complex(cp), intent(inout) :: br_vt_lm(n_lmP_loc) ! br*vt/(sin(theta)**2*r**2)
      complex(cp), intent(inout) :: br_vp_lm(n_lmP_loc) ! br*(vp/(sin(theta)**2*r**2)-omega_ma)

      !-- Local variables:
      integer :: n_theta      ! number of theta position
      integer :: n_phi        ! number of longitude
      real(cp) :: br_vt(nThetaStart:nThetaStop,n_phi_max)
      real(cp) :: br_vp(nThetaStart:nThetaStop,n_phi_max)
      real(cp) :: fac          ! 1/( r**2 sin(theta)**2 )

      !$omp parallel do default(shared) private(n_phi,fac,n_theta)
      do n_phi=1,n_phi_max
         do n_theta=nThetaStart,nThetaStop
            fac=O_sin_theta(n_theta)*O_sin_theta(n_theta)*O_r_E_2*O_rho
            br_vt(n_theta,n_phi)= fac*br(n_theta,n_phi)*vt(n_theta,n_phi)
            br_vp(n_theta,n_phi)= br(n_theta,n_phi) * &
            &                        ( fac*vp(n_theta,n_phi) - omega )
         end do
      end do
      !$omp end parallel do

      call scal_to_SH(br_vt, br_vt_lm, l_max)
      call scal_to_SH(br_vp, br_vp_lm, l_max)

   end subroutine get_br_v_bcs
!----------------------------------------------------------------------------
   subroutine get_b_nl_bcs(bc,br_vt_lm,br_vp_lm,b_nl_bc,aj_nl_bc)
      !
      !  Purpose of this subroutine is to calculate the nonlinear term
      !  of the magnetic boundary condition for a conducting mantle in
      !  physical space (theta,phi), assuming that the conductance
      !  of the mantle is much smaller than that of the core.
      !  Calculation is performed for the theta block:
      !
      !  .. code-block:: fortran
      !
      !                  n_theta_min<=n_theta<=n_theta_min+n_theta_block-1
      !

      !-- Input variables:
      character(len=3), intent(in) :: bc                 ! Distinguishes 'CMB' and 'ICB'
      complex(cp),      intent(in) :: br_vt_lm(n_lmP_loc) ! :math:`B_r u_\theta/(r^2\sin^2\theta)`
      complex(cp),      intent(in) :: br_vp_lm(n_lmP_loc) ! :math:`B_r u_\phi/(r^2\sin^2\theta)`

      !-- Output variables:
      complex(cp), intent(out) :: b_nl_bc(n_lm_loc)  ! nonlinear bc for b
      complex(cp), intent(out) :: aj_nl_bc(n_lm_loc) ! nonlinear bc for aj

      !-- Local variables:
      integer :: l,m       ! degree and order
      integer :: lm        ! position of degree and order
      integer :: lmP       ! same as lm but for l running to l_max+1
      integer :: lmPS,lmPA ! lmP for l-1 and l+1
      real(cp) :: fac

      if ( bc == 'CMB' ) then

         fac=conductance_ma*prmag

         !$omp parallel do default(shared) private(lm,l,m,lmP,lmPS,lmPA)
         do lm=1,n_lm_loc
            l   =map_dist_st%lm2l(lm)
            if ( l == 0 ) cycle
            m   =map_dist_st%lm2m(lm)
            lmP =map_dist_st%lm2lmP(lm)
            lmPS=map_dist_st%lmP2lmPS(lmP)
            lmPA=map_dist_st%lmP2lmPA(lmP)
            if ( l > m ) then
               b_nl_bc(lm)= fac/dLh_loc(lm) * (  dTheta1S_loc(lm)*br_vt_lm(lmPS)  &
               &                               - dTheta1A_loc(lm)*br_vt_lm(lmPA)  &
               &                               +     dPhi_loc(lm)*br_vp_lm(lmP)   )
               aj_nl_bc(lm)=-fac/dLh_loc(lm) * ( dTheta1S_loc(lm)*br_vp_lm(lmPS)  &
               &                               - dTheta1A_loc(lm)*br_vp_lm(lmPA)  &
               &                               -     dPhi_loc(lm)*br_vt_lm(lmP)    )
            else if ( l == m ) then
               b_nl_bc(lm)= fac/dLh_loc(lm) * ( - dTheta1A_loc(lm)*br_vt_lm(lmPA)  &
               &                                    + dPhi_loc(lm)*br_vp_lm(lmP)   )
               aj_nl_bc(lm)=-fac/dLh_loc(lm) * ( - dTheta1A_loc(lm)*br_vp_lm(lmPA) &
               &                                     - dPhi_loc(lm)*br_vt_lm(lmP)   )
            end if
         end do
         !$omp end parallel do

      else if ( bc == 'ICB' ) then

         fac=sigma_ratio*prmag
         !$omp parallel do default(shared) private(lm,l,m,lmP,lmPS,lmPA)
         do lm=1,n_lm_loc
            l   =map_dist_st%lm2l(lm)
            if ( l == 0 ) cycle
            m   =map_dist_st%lm2m(lm)
            lmP =map_dist_st%lm2lmP(lm)
            lmPS=map_dist_st%lmP2lmPS(lmP)
            lmPA=map_dist_st%lmP2lmPA(lmP)
            if ( l > m ) then
               aj_nl_bc(lm)=-fac/dLh_loc(lm) * ( dTheta1S_loc(lm)*br_vp_lm(lmPS)   &
               &                               - dTheta1A_loc(lm)*br_vp_lm(lmPA)   &
               &                               -     dPhi_loc(lm)*br_vt_lm(lmP)    )
            else if ( l == m ) then
               aj_nl_bc(lm)=-fac/dLh_loc(lm) * (- dTheta1A_loc(lm)*br_vp_lm(lmPA) &
               &                                    - dPhi_loc(lm)*br_vt_lm(lmP)    )
            end if
         end do
         !$omp end parallel do

      else
         call abortRun('Wrong input of bc into get_b_nl_bcs')
      end if

   end subroutine get_b_nl_bcs
!-------------------------------------------------------------------------
   subroutine v_rigid_boundary(nR,omega,lDeriv,vrr,vtr,vpr,cvrr,dvrdtr, &
              &                dvrdpr,dvtdpr,dvpdpr)
      !
      !  Purpose of this subroutine is to set the velocities and their
      !  derivatives at a fixed boundary.
      !  While ``vt`` is zero, since we only allow for rotation about the
      !  :math:`z`-axis, ``vp= r \sin(theta) v_phi = r**2 sin(theta)**2 omega``
      !  and ``cvr= r**2 * radial component of (\curl v) = r**2  2 cos(theta) omega``
      !

      !-- Input of variables:
      integer,  intent(in) :: nR            ! no of radial grid point
      logical,  intent(in) :: lDeriv        ! derivatives required ?

      !-- Input of boundary rotation rate
      real(cp), intent(in) :: omega         ! boundary rotation rate

      !-- output:
      real(cp), intent(out) :: vrr(nThetaStart:nThetaStop,n_phi_max)
      real(cp), intent(out) :: vpr(nThetaStart:nThetaStop,n_phi_max)
      real(cp), intent(out) :: vtr(nThetaStart:nThetaStop,n_phi_max)
      real(cp), intent(out) :: cvrr(nThetaStart:nThetaStop,n_phi_max)
      real(cp), intent(out) :: dvrdtr(nThetaStart:nThetaStop,n_phi_max)
      real(cp), intent(out) :: dvrdpr(nThetaStart:nThetaStop,n_phi_max)
      real(cp), intent(out) :: dvtdpr(nThetaStart:nThetaStop,n_phi_max)
      real(cp), intent(out) :: dvpdpr(nThetaStart:nThetaStop,n_phi_max)

      !-- Local variables:
      real(cp) :: r2
      integer :: nTheta,nThetaNHS
      integer :: nPhi

      if ( nR == n_r_cmb ) then
         r2=r_cmb*r_cmb
      else if ( nR == n_r_icb ) then
         r2=r_icb*r_icb
      else
         write(output_unit,*)
         write(output_unit,*) '! v_rigid boundary called for grid'
         write(output_unit,*) '! point which is not a boundary !  '
         return
      end if

      !$omp parallel do default(shared) private(nTheta,nThetaNHS,nPhi)
      do nPhi=1,n_phi_max
         do nTheta=nThetaStart,nThetaStop
            nThetaNHS =(nTheta+1)/2 ! northern hemisphere,sn2 has size n_theta_max/2
            vrr(nTheta,nPhi)=0.0_cp
            vtr(nTheta,nPhi)=0.0_cp
            vpr(nTheta,nPhi)=r2*rho0(nR)*sn2(nThetaNHS)*omega
            if ( lDeriv ) then
               cvrr(nTheta,nPhi)  =r2*rho0(nR)*two*cosTheta(nTheta)*omega
               dvrdtr(nTheta,nPhi)=0.0_cp
               dvrdpr(nTheta,nPhi)=0.0_cp
               dvtdpr(nTheta,nPhi)=0.0_cp
               dvpdpr(nTheta,nPhi)=0.0_cp
            end if
         end do
      end do
      !$omp end parallel do

   end subroutine v_rigid_boundary
!-------------------------------------------------------------------------
end module nonlinear_bcs
