#include "perflib_preproc.cpp"
module radialLoop

   use precision_mod
   use mem_alloc, only: memWrite, bytes_allocated
   use geometry
   use physical_parameters, only: ktopv, kbotv
   use blocking, only: nThetaBs, sizeThetaB
   use logic, only: l_dtB, l_mag, l_mag_LF, lVerbose, l_rot_ma, l_rot_ic,    &
       &            l_cond_ic, l_mag_kin, l_cond_ma, l_mag_nl,               &
       &            l_single_matrix, l_double_curl, l_chemical_conv,         &
       &            l_b_nl_icb, l_b_nl_cmb, l_TP_form
   use constants, only: zero
   use parallel_mod, only: coord_r, n_ranks_r, rank
   use rIteration_mod, only: rIteration_t
   use rIterThetaBlocking_mod, only: rIterThetaBlocking_t
   use rIterThetaBlocking_shtns_mod, only: rIterThetaBlocking_shtns_t

#ifdef WITH_MPI
   use graphOut_mod, only: graphOut_mpi_header
#else
   use graphOut_mod, only: graphOut_header
#endif

   implicit none

   private

   public :: initialize_radialLoop, finalize_radialLoop, radialLoopG

   class(rIteration_t), pointer :: this_rIteration

contains

   subroutine initialize_radialLoop

      character(len=100) :: this_type
      integer(lip) :: local_bytes_used

      local_bytes_used = bytes_allocated
      
      allocate( rIterThetaBlocking_shtns_t :: this_rIteration )

      this_type = this_rIteration%getType()
      if (rank == 0) write(*,"(2A)") "Using rIteration type: ",trim(this_type)
      call this_rIteration%initialize()
      select type (this_rIteration)
         class is (rIterThetaBlocking_t)
           call this_rIteration%set_ThetaBlocking(nThetaBs,sizeThetaB)
         class default
           print*,"this_rIteration has no matching type in radialLoop.f90"
      end select

      local_bytes_used = bytes_allocated-local_bytes_used

      call memWrite('radialLoop.f90', local_bytes_used)

   end subroutine initialize_radialLoop
!----------------------------------------------------------------------------
   subroutine finalize_radialLoop

      call this_rIteration%finalize()
      deallocate(this_rIteration)

   end subroutine finalize_radialLoop
!----------------------------------------------------------------------------
   subroutine radialLoopG(l_graph,l_cour,l_frame,time,dt,dtLast,         &
              &          lTOCalc,lTONext,lTONext2,lHelCalc,lPowerCalc,   &
              &          lRmsCalc,lPressCalc,lViscBcCalc,lFluxProfCalc,  &
              &          lPerpParCalc,l_probe_out,dsdt,dwdt,dzdt,dpdt,   &
              &          dxidt,dbdt,djdt,dVxVhLM,dVxBhLM,dVSrLM,dVPrLM,  &
              &          dVXirLM,lorentz_torque_ic,lorentz_torque_ma,    &
              &          br_vt_lm_cmb,br_vp_lm_cmb,                      &
              &          br_vt_lm_icb,br_vp_lm_icb,                      &
              &          HelLMr,Hel2LMr,HelnaLMr,Helna2LMr,viscLMr,uhLMr,&
              &          duhLMr,gradsLMr,fconvLMr,fkinLMr,fviscLMr,      &
              &          fpoynLMr,fresLMr,EperpLMr,EparLMr,              &
              &          EperpaxiLMr,EparaxiLMr,dtrkc,dthkc)
      !
      !  This subroutine performs the actual time-stepping.
      !

      !--- Input of variables:
      logical,      intent(in) :: l_graph,l_cour,l_frame
      logical,      intent(in) :: lTOcalc,lTONext,lTONext2,lHelCalc
      logical,      intent(in) :: lPowerCalc
      logical,      intent(in) :: lViscBcCalc,lFluxProfCalc,lPerpParCalc
      logical,      intent(in) :: lRmsCalc
      logical,      intent(in) :: l_probe_out
      logical,      intent(in) :: lPressCalc
      real(cp),     intent(in) :: time,dt,dtLast

      !---- Output of explicit time step:
      !---- dVSrLM and dVxBhLM are output of contributions to explicit time step that
      !     need a further treatment (radial derivatives required):
      complex(cp), intent(out) :: dwdt(n_lm_loc,nRstart:nRstop)
      complex(cp), intent(out) :: dzdt(n_lm_loc,nRstart:nRstop)
      complex(cp), intent(out) :: dpdt(n_lm_loc,nRstart:nRstop)
      complex(cp), intent(out) :: dsdt(n_lm_loc,nRstart:nRstop)
      complex(cp), intent(out) :: dxidt(n_lmChe_loc,nRstartChe:nRstopChe)
      complex(cp), intent(out) :: dVSrLM(n_lm_loc,nRstart:nRstop)
      complex(cp), intent(out) :: dVPrLM(n_lmTP_loc,nRstartTP:nRstopTP)
      complex(cp), intent(out) :: dVXirLM(n_lmChe_loc,nRstartChe:nRstopChe)
      complex(cp), intent(out) :: dbdt(n_lmMag_loc,nRstartMag:nRstopMag)
      complex(cp), intent(out) :: djdt(n_lmMag_loc,nRstartMag:nRstopMag)
      complex(cp), intent(out) :: dVxVhLM(n_lmDC_loc,nRstartDC:nRstopDC)
      complex(cp), intent(out) :: dVxBhLM(n_lmMag_loc,nRstartMag:nRstopMag)
      real(cp),    intent(out) :: lorentz_torque_ma,lorentz_torque_ic

      !---- Output for axisymmetric helicity:
      real(cp),    intent(out) :: HelLMr(l_max+1,nRstart:nRstop)
      real(cp),    intent(out) :: Hel2LMr(l_max+1,nRstart:nRstop)
      real(cp),    intent(out) :: HelnaLMr(l_max+1,nRstart:nRstop)
      real(cp),    intent(out) :: Helna2LMr(l_max+1,nRstart:nRstop)
      real(cp),    intent(out) :: uhLMr(l_max+1,nRstart:nRstop)
      real(cp),    intent(out) :: duhLMr(l_max+1,nRstart:nRstop)
      real(cp),    intent(out) :: viscLMr(l_max+1,nRstart:nRstop)
      real(cp),    intent(out) :: gradsLMr(l_max+1,nRstart:nRstop)
      real(cp),    intent(out) :: fkinLMr(l_max+1,nRstart:nRstop)
      real(cp),    intent(out) :: fconvLMr(l_max+1,nRstart:nRstop)
      real(cp),    intent(out) :: fviscLMr(l_max+1,nRstart:nRstop)
      real(cp),    intent(out) :: fresLMr(l_maxMag+1,nRstartMag:nRstopMag)
      real(cp),    intent(out) :: fpoynLMr(l_maxMag+1,nRstartMag:nRstopMag)
      real(cp),    intent(out) :: EperpLMr(l_max+1,nRstart:nRstop)
      real(cp),    intent(out) :: EparLMr(l_max+1,nRstart:nRstop)
      real(cp),    intent(out) :: EperpaxiLMr(l_max+1,nRstart:nRstop)
      real(cp),    intent(out) :: EparaxiLMr(l_max+1,nRstart:nRstop)

      !---- Output of nonlinear products for nonlinear
      !     magnetic boundary conditions (needed in s_updateB.f):
      complex(cp), intent(out) :: br_vt_lm_cmb(n_lmP_loc) ! product br*vt at CMB
      complex(cp), intent(out) :: br_vp_lm_cmb(n_lmP_loc) ! product br*vp at CMB
      complex(cp), intent(out) :: br_vt_lm_icb(n_lmP_loc) ! product br*vt at ICB
      complex(cp), intent(out) :: br_vp_lm_icb(n_lmP_loc) ! product br*vp at ICB

      !---- Output for Courant criteria:
      real(cp),intent(out) :: dtrkc(nRstart:nRstop),dthkc(nRstart:nRstop)


      !--- Local variables:
      integer :: nBc, lm, nR
      integer :: nR_Mag, nRDC, nRTP, nRChe
      integer :: nThetaStart!,nThetaStop
      logical :: lDeriv,lOutBc,lMagNlBc
      logical :: lGraphHeader    ! Write header into graph file
      logical :: isRadialBoundaryPoint

      lGraphHeader=l_graph
      if ( lGraphHeader ) then
#ifdef WITH_MPI
         call graphOut_mpi_header(time,nR,nThetaStart,sizeThetaB)
#else
         call graphOut_header(time)
#endif
      end if
      
      if ( l_cour ) then
         if ( coord_r == 0 ) then
            dtrkc(n_r_cmb)=1.e10_cp
            dthkc(n_r_cmb)=1.e10_cp
         elseif (coord_r == n_ranks_r-1) then
            dtrkc(n_r_icb)=1.e10_cp
            dthkc(n_r_icb)=1.e10_cp
         end if
      end if

      !------ Set nonlinear terms that are possibly needed at the boundaries.
      !       They may be overwritten by get_td later.
      do lm=1,n_lm_loc
         if ( coord_r == 0 ) then
            dVSrLM(lm,n_r_cmb) =zero
            if ( l_chemical_conv ) dVXirLM(lm,n_r_cmb)=zero
            if ( l_mag ) dVxBhLM(lm,n_r_cmb)=zero
            if ( l_double_curl ) dVxVhLM(lm,n_r_cmb)=zero
         elseif (coord_r == n_ranks_r-1) then
            dVSrLM(lm,n_r_icb) =zero
            if ( l_chemical_conv ) dVXirLM(lm,n_r_icb)=zero
            if ( l_mag ) dVxBhLM(lm,n_r_icb)=zero
            if ( l_double_curl ) dVxVhLM(lm,n_r_icb)=zero
         end if
      end do

      !------ Having to calculate non-linear boundary terms?
      lMagNlBc=.false.
      if ( ( l_mag_nl .or. l_mag_kin ) .and.                          &
           &       ( ktopv == 1 .or. l_cond_ma .or.                   &
           &          ( ktopv == 2 .and. l_rot_ma ) ) .or.            &
           &       ( kbotv == 1 .or. l_cond_ic .or.                   &
           &          ( kbotv == 2 .and. l_rot_ic ) ) )               &
           &     lMagNlBc=.true.

      !------ When boundary output, Courant criterion, or non-magnetic 
      !       boundary conditions are required I have to calculate 
      !       the fields at the boundaries. This is done in one thread and 
      !       is triggered by lOutBc=.true.
      lOutBc=.false.
      if ( lTOCalc .or. lHelCalc .or. l_frame .or.         &
           & l_cour .or. l_dtB .or. lMagNlBc .or. l_graph  &
           & .or. lPerpParCalc .or. lViscBcCalc .or.       &
           & lFluxProfCalc .or. lRmsCalc .or. lPowerCalc   &
           & .or. l_single_matrix ) lOutBc=.true.

      !nRstart=n_r_cmb
      !nRstop =n_r_icb-1

      !--- Start the big do loop over the radial threads:
      
      nR_Mag = 1
      nRDC   = 1
      nRTP   = 1
      nRChe  = 1
      
      !nThreadsRmax=1
      do nR=nRstart,nRstop
         !IF( nTh > nThreadsRmax ) nThreadsRmax=nTh
         if ( lVerbose ) then
            write(*,'(/," ! Starting radial level ",i4)') nR
         end if

         !nR = nRC
         nBc = 0
         lDeriv = .true.
         isRadialBoundaryPoint=(nR == n_r_cmb).or.(nR == n_r_icb)
         
         if ((l_mag .or. l_mag_LF)) nR_Mag = nR
         if (l_double_curl        ) nRDC   = nR
         if (l_TP_form            ) nRTP   = nR
         if (l_chemical_conv      ) nRChe  = nR

         if ( nR == n_r_cmb ) then 
            if ( lOutBc ) then
               !nR  = n_r_cmb
               nBc = ktopv
               lDeriv= lTOCalc .or. lHelCalc .or. l_frame .or. lPerpParCalc   &
            &          .or. lViscBcCalc .or. lFluxProfCalc .or. lRmsCalc .or. &
            &          lPowerCalc
            else
               cycle   ! Nothing needs to be done by thread one !
            end if
         elseif ( nR == n_r_icb ) then
            if ( lOutBc ) then
               !nR = n_r_icb
               nBc = kbotv
               lDeriv= lTOCalc .or. lHelCalc .or. l_frame  .or. lPerpParCalc  &
            &          .or. lViscBcCalc .or. lFluxProfCalc .or. lRmsCalc .or. &
            &          lPowerCalc
            else
               cycle
            end if
         end if

         call this_rIteration%set_steering_variables(l_cour,lTOCalc,lTOnext, &
              & lTOnext2,lDeriv,lRmsCalc,lHelCalc,lPowerCalc,l_frame,        &
              & lMagNlBc,l_graph,lViscBcCalc,lFluxProfCalc,lPerpParCalc,     &
              & lPressCalc, l_probe_out)

         call this_rIteration%do_iteration(nR,nBc,time,dt,dtLast,              &
              & dsdt(:,nR),dwdt(:,nR),dzdt(:,nR),dpdt(:,nR),dxidt(:,nR),            &
              & dbdt(:,nR_Mag),djdt(:,nR_Mag),dVxVhLM(:,nRDC),dVxBhLM(:,nR_Mag),                 &
              & dVSrLM(:,nR),dVPrLM(:,nR),dVXirLM(:,nR),br_vt_lm_cmb,        &
              & br_vp_lm_cmb,br_vt_lm_icb,br_vp_lm_icb,lorentz_torque_ic,      &
              & lorentz_torque_ma,HelLMr(:,nR),Hel2LMr(:,nR),HelnaLMr(:,nR),   &
              & Helna2LMr(:,nR),viscLMr(:,nR),uhLMr(:,nR),duhLMr(:,nR),        &
              & gradsLMr(:,nR),fconvLMr(:,nR),fkinLMr(:,nR),fviscLMr(:,nR),    &
              & fpoynLMr(:,nR_Mag),fresLMr(:,nR_Mag),EperpLMr(:,nR),           &
              & EparLMr(:,nR),EperpaxiLMr(:,nR),EparaxiLMr(:,nR))
         
         dtrkc(nR)=this_rIteration%dtrkc      
         dthkc(nR)=this_rIteration%dthkc
      end do
         
      !----- Correct sign of mantel Lorentz torque (see above):
      lorentz_torque_ma=-lorentz_torque_ma
      
   end subroutine radialLoopG
!----------------------------------------------------------------------------
end module radialLoop
