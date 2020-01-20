#include "perflib_preproc.cpp"
module rIterThetaBlocking_shtns_mod
#ifdef WITHOMP
   use omp_lib
#endif
   use precision_mod
   use rIterThetaBlocking_mod, only: rIterThetaBlocking_t
   use geometry
   use logic, only: l_mag, l_conv, l_mag_kin, l_heat, l_ht, l_anel,  &
       &            l_mag_LF, l_conv_nl, l_mag_nl, l_b_nl_cmb,       &
       &            l_b_nl_icb, l_rot_ic, l_cond_ic, l_rot_ma,       &
       &            l_cond_ma, l_dtB, l_store_frame, l_movie_oc,     &
       &            l_TO, l_chemical_conv, l_TP_form, l_probe,       &
       &            l_precession, l_double_curl
   use radial_functions, only: or2, orho1
   use constants, only: zero
   use leg_helper_mod, only: leg_helper_t
   use nonlinear_lm_mod, only:nonlinear_lm_t
   use grid_space_arrays_mod, only: grid_space_arrays_t
   use TO_arrays_mod, only: TO_arrays_t
   use dtB_arrays_mod, only: dtB_arrays_t
   use torsional_oscillations, only: getTO, getTOnext, getTOfinish
#ifdef WITH_MPI
   use graphOut_mod, only: graphOut_mpi
#else
   use graphOut_mod, only: graphOut
#endif
   use dtB_mod, only: get_dtBLM, get_dH_dtBLM
   use out_movie, only: store_movie_frame
   use outRot, only: get_lorentz_torque
   use courant_mod, only: courant
   use nonlinear_bcs, only: get_br_v_bcs, v_rigid_boundary
   use nl_special_calc
   use shtns
   use horizontal_data
   use fields, only: s_Rdist,ds_Rdist, z_Rdist,dz_Rdist, p_Rdist,dp_Rdist, &
       &             b_Rdist,db_Rdist,ddb_Rdist, aj_Rdist,dj_Rdist,       &
       &             w_Rdist,dw_Rdist,ddw_Rdist, xi_Rdist
   use physical_parameters, only: ktops, kbots, n_r_LCR
   use probe_mod
   use parallel_mod
   use blocking, only: nfs

   implicit none

   private

   type, public, extends(rIterThetaBlocking_t) :: rIterThetaBlocking_shtns_t
      integer :: nThreads
      type(grid_space_arrays_t) :: gsa
      type(TO_arrays_t) :: TO_arrays
      type(dtB_arrays_t) :: dtB_arrays
      type(nonlinear_lm_t) :: nl_lm
      real(cp) :: lorentz_torque_ic,lorentz_torque_ma
      
   contains
      procedure :: initialize => initialize_rIterThetaBlocking_shtns
      procedure :: finalize => finalize_rIterThetaBlocking_shtns
      procedure :: do_iteration => do_iteration_ThetaBlocking_shtns
      procedure :: getType => getThisType
      
      ! Distributed Update - Lago
      procedure :: transform_to_grid_space_shtns
      procedure :: transform_to_lm_space_shtns
      
   end type rIterThetaBlocking_shtns_t

contains

   function getThisType(this)

      class(rIterThetaBlocking_shtns_t) :: this
      character(len=100) :: getThisType
      getThisType="rIterThetaBlocking_shtns_t"

   end function getThisType
!------------------------------------------------------------------------------
   subroutine initialize_rIterThetaBlocking_shtns(this)

      class(rIterThetaBlocking_shtns_t) :: this

      if (n_ranks_theta < 2) print *, "n_ranks_theta is too small!!!!!!!!!!!"
      if (n_ranks_theta < 2) stop
      
      call this%allocate_common_arrays()
      call this%gsa%initialize(n_phi_max, l_theta, u_theta)
      if ( l_TO ) call this%TO_arrays%initialize()
      call this%dtB_arrays%initialize()
      
      call this%nl_lm%initialize(n_lmP_loc)
      

   end subroutine initialize_rIterThetaBlocking_shtns
!------------------------------------------------------------------------------
   subroutine finalize_rIterThetaBlocking_shtns(this)

      class(rIterThetaBlocking_shtns_t) :: this

      call this%deallocate_common_arrays()
!       call this%gsa_glb%finalize()
      if ( l_TO ) call this%TO_arrays%finalize()
      call this%dtB_arrays%finalize()
!       call this%nl_lm_glb%finalize()
      
      ! Distributed Update - Lago
      call this%gsa%finalize()
      call this%nl_lm%finalize()

   end subroutine finalize_rIterThetaBlocking_shtns
!------------------------------------------------------------------------------
   subroutine do_iteration_ThetaBlocking_shtns(this,nR,nBc,time,dt,dtLast, &
        &                 dsdt,dwdt,dzdt,dpdt,dxidt,dbdt,djdt,             &
        &                 dVxVhLM,dVxBhLM,dVSrLM,dVPrLM,dVXirLM,           &
        &                 br_vt_lm_cmb,br_vp_lm_cmb,                       &
        &                 br_vt_lm_icb,br_vp_lm_icb,                       &
        &                 lorentz_torque_ic, lorentz_torque_ma,            &
        &                 HelLMr,Hel2LMr,HelnaLMr,Helna2LMr,viscLMr,       &
        &                 uhLMr,duhLMr,gradsLMr,fconvLMr,fkinLMr,fviscLMr, &
        &                 fpoynLMr,fresLMr,EperpLMr,EparLMr,EperpaxiLMr,   &
        &                 EparaxiLMr)

      class(rIterThetaBlocking_shtns_t) :: this
      integer,  intent(in) :: nR,nBc
      real(cp), intent(in) :: time,dt,dtLast

      complex(cp), intent(out) :: dwdt(:),dzdt(:),dpdt(:),dsdt(:),dVSrLM(:)
      complex(cp), intent(out) :: dxidt(:),dVPrLM(:),dVXirLM(:)
      complex(cp), intent(out) :: dbdt(:),djdt(:),dVxVhLM(:),dVxBhLM(:)
      !---- Output of nonlinear products for nonlinear
      !     magnetic boundary conditions (needed in s_updateB.f):
      complex(cp), intent(out) :: br_vt_lm_cmb(:) ! product br*vt at CMB
      complex(cp), intent(out) :: br_vp_lm_cmb(:) ! product br*vp at CMB
      complex(cp), intent(out) :: br_vt_lm_icb(:) ! product br*vt at ICB
      complex(cp), intent(out) :: br_vp_lm_icb(:) ! product br*vp at ICB
      real(cp),    intent(out) :: lorentz_torque_ma, lorentz_torque_ic
      real(cp),    intent(out) :: HelLMr(:),Hel2LMr(:),HelnaLMr(:),Helna2LMr(:)
      real(cp),    intent(out) :: viscLMr(:)
      real(cp),    intent(out) :: uhLMr(:), duhLMr(:) ,gradsLMr(:)
      real(cp),    intent(out) :: fconvLMr(:), fkinLMr(:), fviscLMr(:)
      real(cp),    intent(out) :: fpoynLMr(:), fresLMr(:)
      real(cp),    intent(out) :: EperpLMr(:), EparLMr(:), EperpaxiLMr(:), EparaxiLMr(:)

      integer :: lm
      logical :: lGraphHeader=.false.
      logical :: DEBUG_OUTPUT=.false.
      real(cp) :: c, lorentz_torques_ic
      
      this%nR=nR
      this%nBc=nBc
      this%isRadialBoundaryPoint=(nR == n_r_cmb).or.(nR == n_r_icb)

      if ( this%l_cour ) then
         this%dtrkc=1.e10_cp
         this%dthkc=1.e10_cp
      end if
      if ( this%lTOCalc ) then
         !------ Zero lm coeffs for first theta block:
         call this%TO_arrays%set_zero()
      end if

      call this%leg_helper%legPrepG(this%nR,this%nBc,this%lDeriv,this%lRmsCalc, &
           &                        this%lPressCalc,this%l_frame,this%lTOnext,  &
           &                        this%lTOnext2,this%lTOcalc)

      if (DEBUG_OUTPUT) then
         write(*,"(I3,A,I1,2(A,L1))") this%nR,": nBc = ", &
              & this%nBc,", lDeriv = ",this%lDeriv,", l_mag = ",l_mag
      end if

      this%lorentz_torque_ma = 0.0_cp
      this%lorentz_torque_ic = 0.0_cp
      lorentz_torques_ic = 0.0_cp
      c = 0.0_cp

      br_vt_lm_cmb=zero
      br_vp_lm_cmb=zero
      br_vt_lm_icb=zero
      br_vp_lm_icb=zero
      
      HelLMr     =0.0_cp
      Hel2LMr    =0.0_cp
      HelnaLMr   =0.0_cp
      Helna2LMr  =0.0_cp
      viscLMr    =0.0_cp
      uhLMr      =0.0_cp
      duhLMr     =0.0_cp
      gradsLMr   =0.0_cp
      fconvLMr   =0.0_cp
      fkinLMr    =0.0_cp
      fviscLMr   =0.0_cp
      fpoynLMr   =0.0_cp
      fresLMr    =0.0_cp
      EperpLMr   =0.0_cp
      EparLMr    =0.0_cp
      EperpaxiLMr=0.0_cp
      EparaxiLMr =0.0_cp

      
      call this%nl_lm%set_zero()

      call this%transform_to_grid_space_shtns(time)

      !--------- Calculation of nonlinear products in grid space:
      if ( (.not.this%isRadialBoundaryPoint) .or. this%lMagNlBc .or. &
            this%lRmsCalc ) then

         call this%gsa%get_nl_shtns(time, this%nR, this%nBc, this%lRmsCalc)
         
         call this%transform_to_lm_space_shtns
      
      else if ( l_mag ) then
         do lm=1,n_lmP_loc
            this%nl_lm%VxBtLM(lm)=0.0_cp
            this%nl_lm%VxBpLM(lm)=0.0_cp
         end do
      end if

      !---- Calculation of nonlinear products needed for conducting mantle or
      !     conducting inner core if free stress BCs are applied:
      !     input are brc,vtc,vpc in (theta,phi) space (plus omegaMA and ..)
      !     ouput are the products br_vt_lm_icb, br_vt_lm_cmb, br_vp_lm_icb,
      !     and br_vp_lm_cmb in lm-space, respectively the contribution
      !     to these products from the points theta(nThetaStart)-theta(nThetaStop)
      !     These products are used in get_b_nl_bcs.
      if ( this%nR == n_r_cmb .and. l_b_nl_cmb ) then
         call get_br_v_bcs(this%gsa%brc,this%gsa%vtc,          &
              &            this%gsa%vpc,this%leg_helper%omegaMA,    &
              &            or2(this%nR),orho1(this%nR),                  &
              &            br_vt_lm_cmb,br_vp_lm_cmb)
      else if ( this%nR == n_r_icb .and. l_b_nl_icb ) then
         call get_br_v_bcs(this%gsa%brc,this%gsa%vtc,          &
              &            this%gsa%vpc,this%leg_helper%omegaIC,    &
              &            or2(this%nR),orho1(this%nR),                  &
              &            br_vt_lm_icb,br_vp_lm_icb)
      end if
      
      !--------- Calculate Lorentz torque on inner core:
      !          each call adds the contribution of the theta-block to
      !          lorentz_torque_ic
      if ( this%nR == n_r_icb .and. l_mag_LF .and. l_rot_ic .and. l_cond_ic  ) then
         call get_lorentz_torque(lorentz_torques_ic,                &
              &                  this%gsa%brc,                 &
              &                  this%gsa%bpc,this%nR)
      end if

      !--------- Calculate Lorentz torque on mantle:
      !          note: this calculates a torque of a wrong sign.
      !          sign is reversed at the end of the theta blocking.
      if ( this%nR == n_r_cmb .and. l_mag_LF .and. l_rot_ma .and. l_cond_ma ) then
         call get_lorentz_torque(this%lorentz_torque_ma,   &
              &                  this%gsa%brc,        &
              &                  this%gsa%bpc,this%nR)
      end if
      
      !--------- Calculate courant condition parameters:
      if ( this%l_cour ) then
         !PRINT*,"Calling courant with this%nR=",this%nR
         call courant(this%nR,this%dtrkc,this%dthkc,this%gsa%vrc, &
              &       this%gsa%vtc,this%gsa%vpc,             &
              &       this%gsa%brc,this%gsa%btc,             &
              &       this%gsa%bpc)
      end if

      !--------- Since the fields are given at gridpoints here, this is a good
      !          point for graphical output:
      !< parallelization postponed >
!       if ( this%l_graph ) then
! #ifdef WITH_MPI
!!             call graphOut_mpi(time,this%nR,this%gsa%vrc,           &
!                  &            this%gsa%vtc,this%gsa%vpc,           &
!                  &            this%gsa%brc,this%gsa%btc,           &
!                  &            this%gsa%bpc,this%gsa%sc,            &
!                  &            this%gsa%pc,this%gsa%xic,            &
!                  &            1 ,this%sizeThetaB, lGraphHeader)
!! #else
!             call graphOut(time,this%nR,this%gsa%vrc,           &
!                  &        this%gsa%vtc,this%gsa%vpc,           &
!                  &        this%gsa%brc,this%gsa%btc,           &
!                  &        this%gsa%bpc,this%gsa%sc,            &
!                  &        this%gsa%pc,this%gsa%xic,            &
!                  &        1 ,this%sizeThetaB,lGraphHeader)
! #endif
!       end if
      

      if ( this%l_probe_out ) then
         print *, "PANIC l_probe_out"
         call probe_out(time,this%nR,this%gsa%vpc, 1,this%sizeThetaB)
      end if
      !< / parallelization postponed >
      
      !--------- Helicity output:
      if ( this%lHelCalc ) then
         print *, "PANIC lHelCalc"
         call get_helicity(this%gsa%vrc,this%gsa%vtc,          &
              &        this%gsa%vpc,this%gsa%cvrc,             &
              &        this%gsa%dvrdtc,                        &
              &        this%gsa%dvrdpc,                        &
              &        this%gsa%dvtdrc,                        &
              &        this%gsa%dvpdrc,HelLMr,Hel2LMr,         &
              &        HelnaLMr,Helna2LMr,this%nR,1 )
      end if

      !--------- Viscous heating:
      if ( this%lPowerCalc ) then
         print *, "PANIC lPowerCalc"
         call get_visc_heat(this%gsa%vrc,this%gsa%vtc,this%gsa%vpc,     &
              &        this%gsa%cvrc,this%gsa%dvrdrc,this%gsa%dvrdtc,   &
              &        this%gsa%dvrdpc,this%gsa%dvtdrc,this%gsa%dvtdpc, &
              &        this%gsa%dvpdrc,this%gsa%dvpdpc,viscLMr,         &
              &        this%nR,1)
      end if
  
      !--------- horizontal velocity :
      if ( this%lViscBcCalc ) then
         print *, "PANIC lViscBcCalc"
         call get_nlBLayers(this%gsa%vtc,    &
              &             this%gsa%vpc,    &
              &             this%gsa%dvtdrc, &
              &             this%gsa%dvpdrc, &
              &             this%gsa%drSc,   &
              &             this%gsa%dsdtc,  &
              &             this%gsa%dsdpc,  &
              &             uhLMr,duhLMr,gradsLMr,nR,1 )
      end if


      if ( this%lFluxProfCalc ) then
         print *, "PANIC lFluxProfCalc"
          call get_fluxes(this%gsa%vrc,this%gsa%vtc,             &
                 &        this%gsa%vpc,this%gsa%dvrdrc,          &
                 &        this%gsa%dvtdrc,                       &
                 &        this%gsa%dvpdrc,                       &
                 &        this%gsa%dvrdtc,                       &
                 &        this%gsa%dvrdpc,this%gsa%sc,           &
                 &        this%gsa%pc,this%gsa%brc,              &
                 &        this%gsa%btc,this%gsa%bpc,             &
                 &        this%gsa%cbtc,this%gsa%cbpc,           &
                 &        fconvLMr,fkinLMr,fviscLMr,fpoynLMr,    &
                 &        fresLMr,nR,1 )
      end if

      if ( this%lPerpParCalc ) then
         print *, "PANIC lPerpParCalc"
          call get_perpPar(this%gsa%vrc,this%gsa%vtc,       &
                 &         this%gsa%vpc,EperpLMr,EparLMr,   &
                 &         EperpaxiLMr,EparaxiLMr,nR,1 )
      end if


      !--------- Movie output:
      if ( this%l_frame .and. l_movie_oc .and. l_store_frame ) then
         print *, "PANIC l_frame"
         call store_movie_frame(this%nR,this%gsa%vrc,                &
              &                 this%gsa%vtc,this%gsa%vpc,           &
              &                 this%gsa%brc,this%gsa%btc,           &
              &                 this%gsa%bpc,this%gsa%sc,            &
              &                 this%gsa%drSc,                       &
              &                 this%gsa%dvrdpc,                     &
              &                 this%gsa%dvpdrc,                     &
              &                 this%gsa%dvtdrc,                     &
              &                 this%gsa%dvrdtc,                     &
              &                 this%gsa%cvrc,                       &
              &                 this%gsa%cbrc,                       &
              &                 this%gsa%cbtc,1 ,                    &
              &                 this%sizeThetaB,this%leg_helper%bCMB)
      end if


      !--------- Stuff for special output:
      !--------- Calculation of magnetic field production and advection terms
      !          for graphic output:
      if ( l_dtB ) then
         print *, "PANIC dtBLM"
         call get_dtBLM(this%nR,this%gsa%vrc,this%gsa%vtc,                    &
              &         this%gsa%vpc,this%gsa%brc,                            &
              &         this%gsa%btc,this%gsa%bpc,                            &
              &         1 ,this%sizeThetaB,this%dtB_arrays%BtVrLM,            &
              &         this%dtB_arrays%BpVrLM,this%dtB_arrays%BrVtLM,        &
              &         this%dtB_arrays%BrVpLM,this%dtB_arrays%BtVpLM,        &
              &         this%dtB_arrays%BpVtLM,this%dtB_arrays%BrVZLM,        &
              &         this%dtB_arrays%BtVZLM,this%dtB_arrays%BtVpCotLM,     &
              &         this%dtB_arrays%BpVtCotLM,this%dtB_arrays%BtVZcotLM,  &
              &         this%dtB_arrays%BtVpSn2LM,this%dtB_arrays%BpVtSn2LM,  &
              &         this%dtB_arrays%BtVZsn2LM)
      end if


      !--------- Torsional oscillation terms:
      if ( ( this%lTONext .or. this%lTONext2 ) .and. l_mag ) then
         print *, "PANIC lTONext"
         call getTOnext(this%leg_helper%zAS,this%gsa%brc,   &
              &         this%gsa%btc,this%gsa%bpc,&
              &         this%lTONext,this%lTONext2,dt,dtLast,this%nR, &
              &         1 ,this%sizeThetaB,this%BsLast,      &
              &         this%BpLast,this%BzLast)
      end if

      if ( this%lTOCalc ) then
         print *, "PANIC lTOCalc"
         call getTO(this%gsa%vrc,this%gsa%vtc,    &
              &     this%gsa%vpc,this%gsa%cvrc,   &
              &     this%gsa%dvpdrc,this%gsa%brc, &
              &     this%gsa%btc,this%gsa%bpc,    &
              &     this%gsa%cbrc,this%gsa%cbtc,  &
              &     this%BsLast,this%BpLast,this%BzLast,              &
              &     this%TO_arrays%dzRstrLM,this%TO_arrays%dzAstrLM,  &
              &     this%TO_arrays%dzCorLM,this%TO_arrays%dzLFLM,     &
              &     dtLast,this%nR,1,this%sizeThetaB)
      end if
      

      lorentz_torque_ic = lorentz_torques_ic
      this%lorentz_torque_ic = lorentz_torques_ic
      lorentz_torque_ma = this%lorentz_torque_ma

      if (DEBUG_OUTPUT) then
         print *, "PANIC DEBUG_OUTPUT"
         call this%nl_lm%output()
      end if

      !-- Partial calculation of time derivatives (horizontal parts):
      !   input flm...  is in (l,m) space at radial grid points this%nR !
      !   Only dVxBh needed for boundaries !
      !   get_td finally calculates the d*dt terms needed for the
      !   time step performed in s_LMLoop.f . This should be distributed
      !   over the different models that s_LMLoop.f parallelizes over.
      !write(*,"(A,I4,2ES20.13)") "before_td: ", &
      !     &  this%nR,sum(real(conjg(VxBtLM)*VxBtLM)),sum(real(conjg(VxBpLM)*VxBpLM))
      PERFON('getTd')
      call this%nl_lm%get_td(this%nR, this%nBc, this%lRmsCalc, &
           &                 this%lPressCalc, dVSrLM, dVPrLM, dVXirLM,   &
           &                 dVxVhLM, dVxBhLM, dwdt, dzdt, dpdt, dsdt,   &
           &                 dxidt, dbdt, djdt, this%leg_helper)
      PERFOFF
           
      !write(*,"(A,I4,ES20.13)") "after_td:  ", &
      !     & this%nR,sum(real(conjg(dVxBhLM(:,this%nR_Mag))*dVxBhLM(:,this%nR_Mag)))
      !-- Finish calculation of TO variables:
      if ( this%lTOcalc ) then
         print *, "PANIC lTOcalc"
         call getTOfinish(this%nR, dtLast, this%leg_helper%zAS,             &
              &           this%leg_helper%dzAS, this%leg_helper%ddzAS,      &
              &           this%TO_arrays%dzRstrLM, this%TO_arrays%dzAstrLM, &
              &           this%TO_arrays%dzCorLM, this%TO_arrays%dzLFLM)
      end if
      
      !--- Form partial horizontal derivaties of magnetic production and
      !    advection terms:
      if ( l_dtB ) then
         print *, "PANIC l_dtB"
         call get_dH_dtBLM(this%nR,this%dtB_arrays%BtVrLM,this%dtB_arrays%BpVrLM,&
              &            this%dtB_arrays%BrVtLM,this%dtB_arrays%BrVpLM,        &
              &            this%dtB_arrays%BtVpLM,this%dtB_arrays%BpVtLM,        &
              &            this%dtB_arrays%BrVZLM,this%dtB_arrays%BtVZLM,        &
              &            this%dtB_arrays%BtVpCotLM,this%dtB_arrays%BpVtCotLM,  &
              &            this%dtB_arrays%BtVpSn2LM,this%dtB_arrays%BpVtSn2LM)
      end if
      
    end subroutine do_iteration_ThetaBlocking_shtns
!-------------------------------------------------------------------------------
!
! Distributed Update - Lago
! 
!-------------------------------------------------------------------------------
   subroutine transform_to_grid_space_shtns(this,time)
!@>author Rafael Lago, MPCDF, December 2017
!-------------------------------------------------------------------------------

      class(rIterThetaBlocking_shtns_t) :: this
      real(cp), intent(in) :: time
      integer :: nR

      PERFON('sh2g')
      nR = this%nR

      if ( l_conv .or. l_mag_kin ) then
         if ( l_heat ) then
            call sh_to_spat_dist(s_Rdist(:, nR), this%gsa%sc)
            if ( this%lViscBcCalc ) then
               call sph_to_spat_dist(s_Rdist(:, nR), this%gsa%dsdtc, this%gsa%dsdpc)
               if (this%nR == n_r_cmb .and. ktops==1) then
                  this%gsa%dsdtc=0.0_cp
                  this%gsa%dsdpc=0.0_cp
               end if
               if (this%nR == n_r_icb .and. kbots==1) then
                  this%gsa%dsdtc=0.0_cp
                  this%gsa%dsdpc=0.0_cp
               end if
            end if
         end if

         if ( this%lRmsCalc ) then
            call sph_to_spat_dist(p_Rdist(:, nR), this%gsa%dpdtc, this%gsa%dpdpc)
         end if

         if ( this%lPressCalc ) then ! Pressure
            call sh_to_spat_dist(p_Rdist(:, nR), this%gsa%pc)
         end if

         if ( l_chemical_conv ) then ! Chemical composition
            call sh_to_spat_dist(xi_Rdist(:, nR), this%gsa%xic)
         end if

         if ( l_HT .or. this%lViscBcCalc ) then
            call sh_to_spat_dist(ds_Rdist(:, nR), this%gsa%drsc)
         endif
         if ( this%nBc == 0 ) then
            call torpol_to_spat_dist(w_Rdist(:, nR), dw_Rdist(:, nR),  z_Rdist(:, nR), &
                 &              this%gsa%vrc, this%gsa%vtc, this%gsa%vpc)
            if ( this%lDeriv ) then
               call torpol_to_spat_dist(dw_Rdist(:, nR), ddw_Rdist(:, nR), dz_Rdist(:, nR), &
                    &              this%gsa%dvrdrc, this%gsa%dvtdrc, this%gsa%dvpdrc)

               call pol_to_curlr_spat_dist(z_Rdist(:, nR), this%gsa%cvrc)

               call pol_to_grad_spat_dist(w_Rdist(:, nR), this%gsa%dvrdtc, this%gsa%dvrdpc)
               call torpol_to_dphspat_dist(dw_Rdist(:, nR),  z_Rdist(:, nR), &
                    &                 this%gsa%dvtdpc, this%gsa%dvpdpc)

            end if
         else if ( this%nBc == 1 ) then ! Stress free
             ! TODO don't compute vrc as it is set to 0 afterward
            call torpol_to_spat_dist(w_Rdist(:, nR), dw_Rdist(:, nR),  z_Rdist(:, nR), &
                 &              this%gsa%vrc, this%gsa%vtc, this%gsa%vpc)
            this%gsa%vrc = 0.0_cp
            if ( this%lDeriv ) then
               this%gsa%dvrdtc = 0.0_cp
               this%gsa%dvrdpc = 0.0_cp
               call torpol_to_spat_dist(dw_Rdist(:, nR), ddw_Rdist(:, nR), dz_Rdist(:, nR), &
                    &              this%gsa%dvrdrc, this%gsa%dvtdrc, this%gsa%dvpdrc)
               call pol_to_curlr_spat_dist(z_Rdist(:, nR), this%gsa%cvrc)
               call torpol_to_dphspat_dist(dw_Rdist(:, nR),  z_Rdist(:, nR), &
                    &                 this%gsa%dvtdpc, this%gsa%dvpdpc)
            end if
         else if ( this%nBc == 2 ) then
            if ( this%nR == n_r_cmb ) then
               call v_rigid_boundary(this%nR,this%leg_helper%omegaMA,this%lDeriv, &
                    &                this%gsa%vrc,this%gsa%vtc,this%gsa%vpc,this%gsa%cvrc,this%gsa%dvrdtc, &
                    &                this%gsa%dvrdpc,this%gsa%dvtdpc,this%gsa%dvpdpc,time)
            else if ( this%nR == n_r_icb ) then
               call v_rigid_boundary(this%nR,this%leg_helper%omegaIC,this%lDeriv, &
                    &                this%gsa%vrc,this%gsa%vtc,this%gsa%vpc,this%gsa%cvrc,this%gsa%dvrdtc, &
                    &                this%gsa%dvrdpc,this%gsa%dvtdpc,this%gsa%dvpdpc,time)
            end if
            if ( this%lDeriv ) then
               call torpol_to_spat_dist(dw_Rdist(:, nR), ddw_Rdist(:, nR), dz_Rdist(:, nR), &
                    &              this%gsa%dvrdrc, this%gsa%dvtdrc, this%gsa%dvpdrc)
            end if
         end if
      end if

      if ( l_mag .or. l_mag_LF ) then
         call torpol_to_spat_dist(b_Rdist(:, nR), db_Rdist(:, nR),  aj_Rdist(:, nR),    &
              &              this%gsa%brc, this%gsa%btc, this%gsa%bpc)

         if ( this%lDeriv ) then
            call torpol_to_curl_spat_dist(b_Rdist(:, nR), ddb_Rdist(:, nR),        &
                 &                   aj_Rdist(:, nR), dj_Rdist(:, nR), nR,    &
                 &                   this%gsa%cbrc, this%gsa%cbtc, this%gsa%cbpc)
         end if
      end if
      PERFOFF

   end subroutine transform_to_grid_space_shtns
!-------------------------------------------------------------------------------
   subroutine transform_to_lm_space_shtns(this)
!@>author Rafael Lago, MPCDF, December 2017
!-------------------------------------------------------------------------------
      class(rIterThetaBlocking_shtns_t) :: this
      
      integer :: nPhi, nTheta
      
      PERFON('g2sh')
      if ( (.not.this%isRadialBoundaryPoint .or. this%lRmsCalc) &
            .and. ( l_conv_nl .or. l_mag_LF ) ) then
         if ( l_conv_nl .and. l_mag_LF ) then
            if ( this%nR>n_r_LCR ) then
               do nTheta=l_theta,u_theta
                  do nPhi=1, n_phi_max
                     this%gsa%Advr(nPhi, nTheta)=this%gsa%Advr(nPhi, nTheta) + this%gsa%LFr(nPhi, nTheta)
                     this%gsa%Advt(nPhi, nTheta)=this%gsa%Advt(nPhi, nTheta) + this%gsa%LFt(nPhi, nTheta)
                     this%gsa%Advp(nPhi, nTheta)=this%gsa%Advp(nPhi, nTheta) + this%gsa%LFp(nPhi, nTheta)
                  end do
               end do
            end if
         else if ( l_mag_LF ) then
            if ( this%nR > n_r_LCR ) then
               do nTheta=l_theta, u_theta
                  do nPhi=1, n_phi_max
                     this%gsa%Advr(nPhi, nTheta) = this%gsa%LFr(nPhi, nTheta)
                     this%gsa%Advt(nPhi, nTheta) = this%gsa%LFt(nPhi, nTheta)
                     this%gsa%Advp(nPhi, nTheta) = this%gsa%LFp(nPhi, nTheta)
                  end do
               end do
            else
               do nTheta=l_theta, u_theta
                  do nPhi=1, n_phi_max
                     this%gsa%Advr(nPhi,nTheta)=0.0_cp
                     this%gsa%Advt(nPhi,nTheta)=0.0_cp
                     this%gsa%Advp(nPhi,nTheta)=0.0_cp
                  end do
               end do
            end if
         end if

         if ( l_precession ) then
            do nTheta=l_theta,u_theta
               do nPhi=1, n_phi_max
                  this%gsa%Advr(nPhi, nTheta)=this%gsa%Advr(nPhi, nTheta) + this%gsa%PCr(nPhi, nTheta)
                  this%gsa%Advt(nPhi, nTheta)=this%gsa%Advt(nPhi, nTheta) + this%gsa%PCt(nPhi, nTheta)
                  this%gsa%Advp(nPhi, nTheta)=this%gsa%Advp(nPhi, nTheta) + this%gsa%PCp(nPhi, nTheta)
               end do
            end do

         end if

         call spat_to_SH_dist(this%gsa%Advr, this%nl_lm%AdvrLM)
         call spat_to_SH_dist(this%gsa%Advt, this%nl_lm%AdvtLM)
         call spat_to_SH_dist(this%gsa%Advp, this%nl_lm%AdvpLM)

         if ( this%lRmsCalc .and. l_mag_LF .and. this%nR>n_r_LCR ) then
            ! LF treated extra:
            call spat_to_SH_dist(this%gsa%LFr, this%nl_lm%LFrLM)
            call spat_to_SH_dist(this%gsa%LFt, this%nl_lm%LFtLM)
            call spat_to_SH_dist(this%gsa%LFp, this%nl_lm%LFpLM)
         end if
      end if
      if ( (.not.this%isRadialBoundaryPoint) .and. l_heat ) then
         call spat_to_SH_dist(this%gsa%VSr, this%nl_lm%VSrLM)
         call spat_to_SH_dist(this%gsa%VSt, this%nl_lm%VStLM)
         call spat_to_SH_dist(this%gsa%VSp, this%nl_lm%VSpLM)

         if (l_anel) then ! anelastic stuff
            if ( l_mag_nl .and. this%nR>n_r_LCR ) then
               call spat_to_SH_dist(this%gsa%ViscHeat, this%nl_lm%ViscHeatLM)
               call spat_to_SH_dist(this%gsa%OhmLoss, this%nl_lm%OhmLossLM)
            else
               call spat_to_SH_dist(this%gsa%ViscHeat, this%nl_lm%ViscHeatLM)
            end if
         end if
      end if
      if ( (.not.this%isRadialBoundaryPoint) .and. l_TP_form ) then
         call spat_to_SH_dist(this%gsa%VPr, this%nl_lm%VPrLM)
      end if
      if ( (.not.this%isRadialBoundaryPoint) .and. l_chemical_conv ) then
         call spat_to_SH_dist(this%gsa%VXir, this%nl_lm%VXirLM)
         call spat_to_SH_dist(this%gsa%VXit, this%nl_lm%VXitLM)
         call spat_to_SH_dist(this%gsa%VXip, this%nl_lm%VXipLM)
      end if
      if ( l_mag_nl ) then
         if ( .not.this%isRadialBoundaryPoint .and. this%nR>n_r_LCR ) then
            call spat_to_SH_dist(this%gsa%VxBr, this%nl_lm%VxBrLM)
            call spat_to_SH_dist(this%gsa%VxBt, this%nl_lm%VxBtLM)
            call spat_to_SH_dist(this%gsa%VxBp, this%nl_lm%VxBpLM)
         else
            call spat_to_SH_dist(this%gsa%VxBt, this%nl_lm%VxBtLM)
            call spat_to_SH_dist(this%gsa%VxBp, this%nl_lm%VxBpLM)
         end if
      end if

      if ( this%lRmsCalc ) then
         call spat_to_sphtor_dist(this%gsa%dpdtc, this%gsa%dpdpc,    &
                                 this%nl_lm%PFt2LM, this%nl_lm%PFp2LM)
         call spat_to_sphtor_dist(this%gsa%CFt2, this%gsa%CFp2,      &
                                 this%nl_lm%CFt2LM, this%nl_lm%CFp2LM)
         if ( l_conv_nl ) then
            call spat_to_sphtor_dist(this%gsa%Advt2, this%gsa%Advp2, &
                               this%nl_lm%Advt2LM, this%nl_lm%Advp2LM)
         end if
         if ( l_mag_nl .and. this%nR>n_r_LCR ) then
            call spat_to_sphtor_dist(this%gsa%LFt2, this%gsa%LFp2,   &
                                 this%nl_lm%LFt2LM, this%nl_lm%LFp2LM)
         end if
      end if
      PERFOFF
      
   end subroutine transform_to_lm_space_shtns
end module rIterThetaBlocking_shtns_mod
