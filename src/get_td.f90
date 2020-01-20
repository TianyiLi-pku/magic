module nonlinear_lm_mod

   use, intrinsic :: iso_c_binding
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use geometry, only: lm_max, l_max, lm_maxMag, lmP_max, n_lm_loc,         &
       &             n_lmMag_loc, n_lmP_loc, n_lmDC_loc, n_lmChe_loc, n_lmTP_loc
   use communications, only: slice_FlmP, gather_FlmP
   use logic, only : l_anel, l_conv_nl, l_corr, l_heat, l_anelastic_liquid, &
       &             l_mag_nl, l_mag_kin, l_mag_LF, l_conv, l_mag, l_RMS,   &
       &             l_chemical_conv, l_TP_form, l_single_matrix, l_double_curl
   use radial_functions, only: r, or2, or1, beta, rho0, rgrav, epscProf,    &
       &             or4, temp0, alpha0, ogrun, orho1
   use physical_parameters, only: CorFac, ra, epsc, ViscHeatFac,            &
       &             OhmLossFac, n_r_LCR, epscXi, BuoFac, ThExpNb
   use horizontal_data, only: dLh_loc,  dTheta1S_loc, dTheta1A_loc,         &
       &             dPhi_loc, dTheta2A_loc, dTheta3A_loc, dTheta4A_loc,    &
       &             dPhi0_loc, dTheta2S_loc, dTheta3S_loc, dTheta4S_loc,   &
       &             hdif_V_loc, hdif_B_loc
   use RMS, only: Adv2hInt, Pre2hInt, Buo2hInt, Cor2hInt, LF2hInt,  &
       &          Geo2hInt, Mag2hInt, ArcMag2hInt, CLF2hInt, PLF2hInt, &
       &          CIA2hInt, Arc2hInt
   use leg_helper_mod, only: leg_helper_t
   use constants, only: zero, two
   use fields, only: w_Rdist, dw_Rdist, ddw_Rdist, z_Rdist, dz_Rdist
   use RMS_helpers, only: hIntRms
   use LMmapping, only: map_dist_st, map_glbl_st

   implicit none
   
   type :: nonlinear_lm_t
      !----- Nonlinear terms in lm-space: 
      complex(cp), allocatable :: AdvrLM(:), AdvtLM(:), AdvpLM(:)
      complex(cp), allocatable :: LFrLM(:),  LFtLM(:),  LFpLM(:)
      complex(cp), allocatable :: VxBrLM(:), VxBtLM(:), VxBpLM(:)
      complex(cp), allocatable :: VSrLM(:),  VStLM(:),  VSpLM(:)
      complex(cp), allocatable :: VXirLM(:),  VXitLM(:),  VXipLM(:)
      complex(cp), allocatable :: VPrLM(:)
      complex(cp), allocatable :: ViscHeatLM(:), OhmLossLM(:)
      !----- RMS calculations
      complex(cp), allocatable :: Advt2LM(:), Advp2LM(:)
      complex(cp), allocatable :: LFt2LM(:), LFp2LM(:)
      complex(cp), allocatable :: CFt2LM(:), CFp2LM(:)
      complex(cp), allocatable :: PFt2LM(:), PFp2LM(:)
 
   contains
 
      procedure :: initialize
      procedure :: finalize
      procedure :: output
      procedure :: set_zero
      procedure :: get_td
      procedure :: slice_all  => slice_all_nonlinear_lm
      procedure :: gather_all  => gather_all_nonlinear_lm
 
   end type nonlinear_lm_t

contains

   subroutine initialize(this,lmP_length)

      class(nonlinear_lm_t) :: this
      integer, intent(in)   :: lmP_length

      allocate( this%AdvrLM(lmP_length) )   
      allocate( this%AdvtLM(lmP_length) )   
      allocate( this%AdvpLM(lmP_length) )   
      allocate( this%LFrLM(lmP_length) )    
      allocate( this%LFtLM(lmP_length) )    
      allocate( this%LFpLM(lmP_length) )    
      allocate( this%VxBrLM(lmP_length) )   
      allocate( this%VxBtLM(lmP_length) )   
      allocate( this%VxBpLM(lmP_length) )   
      allocate( this%VSrLM(lmP_length) )    
      allocate( this%VStLM(lmP_length) )    
      allocate( this%VSpLM(lmP_length) )    
      allocate( this%ViscHeatLM(lmP_length) )
      allocate( this%OhmLossLM(lmP_length) )
      bytes_allocated = bytes_allocated + 14*lmP_length*SIZEOF_DEF_COMPLEX

      if ( l_TP_form ) then
         allocate( this%VPrLM(lmP_length) )    
         bytes_allocated = bytes_allocated + lmP_length*SIZEOF_DEF_COMPLEX
      end if
      
      if ( l_chemical_conv ) then
         allocate( this%VXirLM(lmP_length) )    
         allocate( this%VXitLM(lmP_length) )    
         allocate( this%VXipLM(lmP_length) )    
         bytes_allocated = bytes_allocated + 3*lmP_length*SIZEOF_DEF_COMPLEX
      end if

      !-- RMS calculations
      if ( l_RMS ) then
         allocate( this%Advt2LM(lmP_length) )
         allocate( this%Advp2LM(lmP_length) )
         allocate( this%LFt2LM(lmP_length) )
         allocate( this%LFp2LM(lmP_length) )
         allocate( this%CFt2LM(lmP_length) )
         allocate( this%CFp2LM(lmP_length) )
         allocate( this%PFt2LM(lmP_length) )
         allocate( this%PFp2LM(lmP_length) )
         bytes_allocated = bytes_allocated + 8*lmP_length*SIZEOF_DEF_COMPLEX
      end if

   end subroutine initialize
!----------------------------------------------------------------------------
   subroutine finalize(this)

      class(nonlinear_lm_t) :: this

      deallocate( this%AdvrLM )   
      deallocate( this%AdvtLM )   
      deallocate( this%AdvpLM )   
      deallocate( this%LFrLM )    
      deallocate( this%LFtLM )    
      deallocate( this%LFpLM )    
      deallocate( this%VxBrLM )   
      deallocate( this%VxBtLM )   
      deallocate( this%VxBpLM )   
      deallocate( this%VSrLM )    
      deallocate( this%VStLM )    
      deallocate( this%VSpLM )    
      deallocate( this%ViscHeatLM )
      deallocate( this%OhmLossLM )

      if ( l_TP_form ) then
         deallocate( this%VPrLM )    
      end if

      if ( l_chemical_conv ) then
         deallocate( this%VXirLM )    
         deallocate( this%VXitLM )    
         deallocate( this%VXipLM )    
      end if

      !-- RMS calculations
      if ( l_RMS ) then
         deallocate( this%Advt2LM )
         deallocate( this%Advp2LM )
         deallocate( this%LFt2LM )
         deallocate( this%LFp2LM )
         deallocate( this%CFt2LM )
         deallocate( this%CFp2LM )
         deallocate( this%PFt2LM )
         deallocate( this%PFp2LM )
      end if

   end subroutine finalize
!----------------------------------------------------------------------------
   subroutine set_zero(this)

      class(nonlinear_lm_t) :: this
      
      this%AdvrLM    =zero
      this%AdvtLM    =zero
      this%AdvpLM    =zero
      this%LFrLM     =zero
      this%LFtLM     =zero
      this%LFpLM     =zero
      this%VxBrLM    =zero
      this%VxBtLM    =zero
      this%VxBpLM    =zero
      this%VSrLM     =zero
      this%VStLM     =zero
      this%VSpLM     =zero
      this%ViscHeatLM=zero
      this%OhmLossLM =zero

      if ( l_TP_form ) then
         this%VPrLM =zero
      end if

      if ( l_chemical_conv ) then
         this%VXirLM =zero
         this%VXitLM =zero
         this%VXipLM =zero
      end if

      if ( l_RMS ) then
         this%Advt2LM=zero
         this%Advp2LM=zero
         this%LFp2LM =zero
         this%LFt2LM =zero
         this%CFt2LM =zero
         this%CFp2LM =zero
         this%PFt2LM =zero
         this%PFp2LM =zero
      end if

   end subroutine set_zero
!----------------------------------------------------------------------------
   subroutine output(this)

      class(nonlinear_lm_t) :: this
      
      write(*,"(A,6ES20.12)") "AdvrLM,AdvtLM,AdvpLM = ",&
           & sum(this%AdvrLM), sum(this%AdvtLM),        &
           & sum(this%AdvpLM)

   end subroutine output

!----------------------------------------------------------------------------
subroutine get_td(this,nR,nBc,lRmsCalc,lPressCalc,dVSrLM,dVPrLM,dVXirLM, &
              &      dVxVhLM,dVxBhLM,dwdt,dzdt,dpdt,dsdt,dxidt,dbdt,djdt,   &
              &      leg_helper)
      !
      !  Purpose of this to calculate time derivatives
      !  dwdt,dzdt,dpdt,dsdt,dxidt,dbdt,djdt
      !  and auxiliary arrays dVSrLM, dVXirLM and dVxBhLM, dVxVhLM
      !  from non-linear terms in spectral form,
      !  contained in flmw1-3,flms1-3, flmb1-3 (input)
    
      !-- Input of variables:
      class(nonlinear_lm_t) :: this

      integer,            intent(in) :: nR
      integer,            intent(in) :: nBc ! signifies boundary conditions
      logical,            intent(in) :: lRmsCalc
      logical,            intent(in) :: lPressCalc
      type(leg_helper_t), intent(in) :: leg_helper
    
      !-- Output of variables:
      complex(cp), intent(out) :: dwdt(n_lm_loc),dzdt(n_lm_loc)
      complex(cp), intent(out) :: dpdt(n_lm_loc),dsdt(n_lm_loc)
      complex(cp), intent(out) :: dxidt(n_lmChe_loc)
      complex(cp), intent(out) :: dbdt(n_lmMag_loc),djdt(n_lmMag_loc)
      complex(cp), intent(out) :: dVxBhLM(n_lmMag_loc)
      complex(cp), intent(out) :: dVxVhLM(n_lmDC_loc)
      complex(cp), intent(out) :: dVSrLM(n_lm_loc)
      complex(cp), intent(out) :: dVXirLM(n_lmChe_loc)
      complex(cp), intent(out) :: dVPrLM(n_lmTP_loc)
    
      !-- Local variables:
      integer :: l,m,lm,lmS,lmA,lmP,lmPS,lmPA
      complex(cp) :: CorPol(n_lm_loc)
      complex(cp) :: AdvPol(n_lm_loc),AdvTor(n_lm_loc)
      complex(cp) :: LFPol(n_lm_loc),LFTor(n_lm_loc)
      complex(cp) :: Geo(n_lm_loc),CLF(n_lm_loc),PLF(n_lm_loc)
      complex(cp) :: ArcMag(n_lm_loc),Mag(n_lm_loc),CIA(n_lm_loc),Arc(n_lm_loc)
      complex(cp) :: Buo(n_lm_loc)
      complex(cp) :: AdvPol_loc,CorPol_loc,AdvTor_loc,CorTor_loc
      complex(cp) :: dsdt_loc, dxidt_loc
    
      integer, parameter :: DOUBLE_COMPLEX_PER_CACHELINE=4
      
      integer :: lm_maybe_skip_first
      
      lm_maybe_skip_first = 1
      if (map_dist_st%lm2(0,0) > 0) lm_maybe_skip_first = 2
    
      if (nBc == 0 .or. lRmsCalc ) then
    
         if ( l_conv ) then  ! Convection
            
            if (map_dist_st%lm2(0,0) > 0) then  ! if m=0 is in this rank
               lm  = map_dist_st%lm2(0,0)
               lmA = map_dist_st%lm2lmA(lm)
               lmP = map_dist_st%lm2lmP(lm)
               lmPA= map_dist_st%lmP2lmPA(lmP)
               
               if ( l_conv_nl ) then
                  AdvPol_loc=      or2(nR)*this%AdvrLM(lm)
                  AdvTor_loc=-dTheta1A_loc(lm)*this%AdvpLM(lmPA)
               else
                  AdvPol_loc=zero
                  AdvTor_loc=zero
               end if
               if ( l_corr ) then
                  CorPol_loc=two*CorFac*or1(nR) * dTheta2A_loc(lm)* z_Rdist(lmA,nR)
                  CorTor_loc= two*CorFac*or2(nR) * (                 &
                  &                dTheta3A_loc(lm)*dw_Rdist(lmA,nR) +    &
                  &        or1(nR)*dTheta4A_loc(lm)* w_Rdist(lmA,nR) )
               else
                  CorPol_loc=zero
                  CorTor_loc=zero
               end if

               if ( l_single_matrix ) then
                  dwdt(lm)=AdvPol_loc!+CorPol_loc
               else
                  dwdt(lm)=AdvPol_loc+CorPol_loc
               end if

               dzdt(lm)=AdvTor_loc+CorTor_loc

               if ( lRmsCalc ) then

                  Buo(lm) =BuoFac*rgrav(nR)*rho0(nR)*leg_helper%sR(lm)
                  if ( l_mag_LF .and. nR>n_r_LCR ) then
                     LFPol(lm) =      or2(nR)*this%LFrLM(lm)
                     LFTor(lm) =-dTheta1A_loc(lm)*this%LFpLM(lmPA)
                     AdvPol(lm)=AdvPol_loc-LFPol(lm)
                     AdvTor(lm)=AdvTor_loc-LFTor(lm)
                  else
                     AdvPol(lm)=AdvPol_loc
                     AdvTor(lm)=AdvTor_loc
                  end if
                  CorPol(lm)=CorPol_loc
               end if
            end if
            
            do lm=lm_maybe_skip_first,n_lm_loc
               l   =map_dist_st%lm2l(lm)
               m   =map_dist_st%lm2m(lm)
               lmS =map_dist_st%lm2lmS(lm)
               lmA =map_dist_st%lm2lmA(lm)
               lmP =map_dist_st%lm2lmP(lm)
               lmPS=map_dist_st%lmP2lmPS(lmP)
               lmPA=map_dist_st%lmP2lmPA(lmP)
               
               if ( l_double_curl ) then ! Pressure is not needed

                  if ( l_corr ) then
                     if ( l < l_max .and. l > m ) then
                        CorPol_loc =two*CorFac*or2(nR)*orho1(nR)*(               &
                        &                   dPhi0_loc(lm)*(                      &
                        &         -ddw_Rdist(lm,nR)+beta(nR)*dw_Rdist(lm,nR)     + &
                        &             ( beta(nR)*or1(nR)+or2(nR))*               &
                        &                         leg_helper%dLhw(lm) )    + &
                        &             dTheta3A_loc(lm)*( dz_Rdist(lmA,nR)-        &
                        &                            beta(nR)*z_Rdist(lmA,nR) ) + &
                        &             dTheta3S_loc(lm)*( dz_Rdist(lmS,nR)-        &
                        &                            beta(nR)*z_Rdist(lmS,nR) ) + &
                        &          or1(nR)* (                                    &
                        &             dTheta4A_loc(lm)* z_Rdist(lmA,nR)           &
                        &            -dTheta4S_loc(lm)* z_Rdist(lmS,nR) ) )
                     else if ( l == l_max ) then
                        CorPol_loc =two*CorFac*or2(nR)*orho1(nR)*(               &
                        &                   dPhi0_loc(lm)*(                      &
                        &         -ddw_Rdist(lm,nR)+beta(nR)*dw_Rdist(lm,nR)     + &
                        &             ( beta(nR)*or1(nR)+or2(nR))*               &
                        &                         leg_helper%dLhw(lm) ) )
                     else if ( l == m ) then
                        CorPol_loc =two*CorFac*or2(nR)*orho1(nR)*(               &
                        &                   dPhi0_loc(lm)*(                      &
                        &         -ddw_Rdist(lm,nR)+beta(nR)*dw_Rdist(lm,nR)     + &
                        &             ( beta(nR)*or1(nR)+or2(nR))*               &
                        &                         leg_helper%dLhw(lm) )    + &
                        &             dTheta3A_loc(lm)*( dz_Rdist(lmA,nR)-        &
                        &                            beta(nR)*z_Rdist(lmA,nR) ) + &
                        &          or1(nR)* (                                    &
                        &             dTheta4A_loc(lm)* z_Rdist(lmA,nR) ) )
                     end if
                  else
                     CorPol_loc=zero
                  end if

                  if ( l_conv_nl ) then

                     if ( l > m ) then
                        dVxVhLM(lm)=      orho1(nR)*r(nR)*r(nR)* ( &
                        &        dTheta1S_loc(lm)*this%AdvtLM(lmPS) -  &
                        &        dTheta1A_loc(lm)*this%AdvtLM(lmPA) +  &
                        &             dPhi_loc(lm)*this%AdvpLM(lmP)  )
                     else if ( l == m ) then
                        dVxVhLM(lm)=      orho1(nR)*r(nR)*r(nR)* ( &
                        &      - dTheta1A_loc(lm)*this%AdvtLM(lmPA) +  &
                        &        dPhi_loc(lm)*this%AdvpLM(lmP)  )
                     end if

                     AdvPol_loc=dLh_loc(lm)*or4(nR)*orho1(nR)*this%AdvrLM(lmP)

                  else

                     AdvPol_loc =zero
                     dVxVhLM(lm)=zero

                  endif

               else ! We don't use the double curl

                  if ( l_corr .and. nBc /= 2 ) then
                     if ( l < l_max .and. l > m ) then
                        CorPol_loc =two*CorFac*or1(nR) * (      &
                        &       dPhi0_loc(lm)*dw_Rdist(lm,nR) +  & ! phi-deriv of dw/dr
                        &    dTheta2A_loc(lm)*z_Rdist(lmA,nR) -  & ! sin(theta) dtheta z
                        &    dTheta2S_loc(lm)*z_Rdist(lmS,nR) )
                     else if ( l == l_max ) then
                        CorPol_loc= two*CorFac*or1(nR) * (      &
                        &            dPhi0_loc(lm)*dw_Rdist(lm,nR)  )
                     else if ( l == m ) then
                        CorPol_loc = two*CorFac*or1(nR) * (      &
                        &        dPhi0_loc(lm)*dw_Rdist(lm,nR)  + &
                        &     dTheta2A_loc(lm)*z_Rdist(lmA,nR) )
                     end if
                  else
                     CorPol_loc=zero
                  end if

                  if ( l_conv_nl ) then
                     AdvPol_loc=or2(nR)*this%AdvrLM(lmP)
                  else
                     AdvPol_loc=zero
                  endif

               end if ! Double curl or not for the poloidal equation

               dwdt(lm)=AdvPol_loc+CorPol_loc

               if ( lRmsCalc ) then ! RMS force balance

                  if ( l_TP_form .or. l_anelastic_liquid ) then
                     Buo(lm) =BuoFac*alpha0(nR)*rgrav(nR)*(              &
                     &        rho0(nR)*leg_helper%sR(lm)-ViscHeatFac*&
                     &        (ThExpNb*alpha0(nR)*temp0(nR)+ogrun(nR))*  &
                     &        leg_helper%preR(lm) )
                  else
                     Buo(lm) =BuoFac*rho0(nR)*rgrav(nR)*leg_helper%sR(lm)
                  end if

                  if ( l_double_curl ) then 
                     ! In that case we have to recompute the Coriolis force
                     ! since we also want the pressure gradient
                     if ( l_corr .and. nBc /= 2 ) then
                        if ( l < l_max .and. l > m ) then
                           CorPol_loc =two*CorFac*or1(nR) * (      &
                           &       dPhi0_loc(lm)*dw_Rdist(lm,nR) +  & ! phi-deriv of dw/dr
                           &    dTheta2A_loc(lm)*z_Rdist(lmA,nR) -  & ! sin(theta) dtheta z
                           &    dTheta2S_loc(lm)*z_Rdist(lmS,nR) )
                        else if ( l == l_max ) then
                           CorPol_loc= two*CorFac*or1(nR) * ( &
                           &            dPhi0_loc(lm)*dw_Rdist(lm,nR)  )
                        else if ( l == m ) then
                           CorPol_loc = two*CorFac*or1(nR) * (      &
                           &        dPhi0_loc(lm)*dw_Rdist(lm,nR)  + &
                           &     dTheta2A_loc(lm)*z_Rdist(lmA,nR) )
                        end if
                     else
                        CorPol_loc=zero
                     end if

                     ! We also need to recompute AdvPol_loc here
                     if ( l_conv_nl ) then
                        AdvPol_loc=or2(nR)*this%AdvrLM(lmP)
                     else
                        AdvPol_loc=zero
                     endif

                  end if

                  if ( l_mag_LF .and. nR>n_r_LCR ) then
                     LFPol(lm) =or2(nR)*this%LFrLM(lmP)
                     AdvPol(lm)=AdvPol_loc-LFPol(lm)
                  else
                     AdvPol(lm)=AdvPol_loc
                  end if
                  CorPol(lm)=CorPol_loc

               end if

               if ( l_corr ) then
                  if ( l < l_max .and. l > m ) then
                     CorTor_loc=          two*CorFac*or2(nR) * (      &
                     &                dPhi0_loc(lm)*z_Rdist(lm,nR)   + &
                     &            dTheta3A_loc(lm)*dw_Rdist(lmA,nR)  + &
                     &    or1(nR)*dTheta4A_loc(lm)* w_Rdist(lmA,nR)  + &
                     &            dTheta3S_loc(lm)*dw_Rdist(lmS,nR)  - &
                     &    or1(nR)*dTheta4S_loc(lm)* w_Rdist(lmS,nR)  )
                  else if ( l == l_max ) then
                     CorTor_loc=two*CorFac*or2(nR) * ( &
                     &            dPhi0_loc(lm)*z_Rdist(lm,nR)   )
                  else if ( l == m ) then
                     CorTor_loc=  two*CorFac*or2(nR) * (  &
                     &        dPhi0_loc(lm)*z_Rdist(lm,nR)   + &
                     &    dTheta3A_loc(lm)*dw_Rdist(lmA,nR)  + &
                     &    or1(nR)*dTheta4A_loc(lm)* w_Rdist(lmA,nR)  )
                  end if
               else
                  CorTor_loc=zero
               end if
    
               if ( l_conv_nl ) then
                  if ( l > m ) then
                     AdvTor_loc=   -dPhi_loc(lm)*this%AdvtLM(lmP)  + &
                     &          dTheta1S_loc(lm)*this%AdvpLM(lmPS) - &
                     &          dTheta1A_loc(lm)*this%AdvpLM(lmPA)
                  else if ( l == m ) then
                     AdvTor_loc=   -dPhi_loc(lm)*this%AdvtLM(lmP)  - &
                     &          dTheta1A_loc(lm)*this%AdvpLM(lmPA)
                  end if
               else
                  AdvTor_loc=zero
               end if
    
               dzdt(lm)=CorTor_loc+AdvTor_loc
               ! until here
    
               if ( lRmsCalc ) then
                  if ( l_mag_LF .and. nR>n_r_LCR ) then
                     !------ When RMS values are required, the Lorentz force is treated
                     !       separately:
       
                     if ( l > m ) then
                        !------- LFTor= 1/(E*Pm) * curl( curl(B) x B )_r
                        LFTor(lm) =   -dPhi_loc(lm)*this%LFtLM(lmP)  + &
                        &          dTheta1S_loc(lm)*this%LFpLM(lmPS) - &
                        &          dTheta1A_loc(lm)*this%LFpLM(lmPA)
                     else if ( l == m ) then
                        LFTor(lm) =   -dPhi_loc(lm)*this%LFtLM(lmP)  - &
                        &          dTheta1A_loc(lm)*this%LFpLM(lmPA)
                     end if
                     AdvTor(lm)=AdvTor_loc-LFTor(lm)
                  else
                     AdvTor(lm)=AdvTor_loc
                  end if
               end if
    
            end do

!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!
!           BEGIN OF POSTPONED
!           Postponed. This has to do with Diagnostics!
!           
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------    
            if ( lRmsCalc ) then
               print *, "lRmsCalc not yet parallelized! @", __LINE__, __FILE__
               stop
    
               if ( l_conv_nl ) then
                  call hIntRms(AdvPol,nR,1,lm_max,0,Adv2hInt(:,nR),map_glbl_st, .false.)
                  call hIntRms(this%Advt2LM,nR,1,lmP_max,1,Adv2hInt(:,nR),map_glbl_st, &
                       &       .true.)
                  call hIntRms(this%Advp2LM,nR,1,lmP_max,1,Adv2hInt(:,nR),map_glbl_st, &
                       &       .true.)
               end if

               if ( l_TP_form .or. l_anelastic_liquid ) then
                  call hIntRms(leg_helper%dpR,nR,1,lm_max,0, &
                       &       Pre2hInt(:,nR),map_glbl_st,.false.)
               else
                  call hIntRms(leg_helper%dpR-beta(nR)*leg_helper%preR,&
                       &       nR,1,lm_max,0,Pre2hInt(:,nR),map_glbl_st,.false.)
               end if
               call hIntRms(this%PFt2LM,nR,1,lmP_max,1,Pre2hInt(:,nR),map_glbl_st,.true.)
               call hIntRms(this%PFp2LM,nR,1,lmP_max,1,Pre2hInt(:,nR),map_glbl_st,.true.)

               ! rho* grad(p/rho) = grad(p) - beta*p
               if ( ra /= 0.0_cp ) &
                  call hIntRms(Buo,nR,1,lm_max,0,Buo2hInt(:,nR),map_glbl_st,.false.)
               if ( l_corr ) then
                  call hIntRms(CorPol,nR,1,lm_max,0,Cor2hInt(:,nR),map_glbl_st,.false.)
                  call hIntRms(this%CFt2LM,nR,1,lmP_max,1,Cor2hInt(:,nR), &
                       &       map_glbl_st,.true.)
                  calL hIntRms(this%CFp2LM,nR,1,lmP_max,1,Cor2hInt(:,nR), &
                       &       map_glbl_st,.true.)
               end if
               if ( l_mag_LF .and. nR>n_r_LCR ) then
                  call hIntRms(LFPol,nR,1,lm_max,0,LF2hInt(:,nR),map_glbl_st,.false.)
                  call hIntRms(this%LFt2LM,nR,1,lmP_max,1,LF2hInt(:,nR),map_glbl_st,.true.)
                  call hIntRms(this%LFp2LM,nR,1,lmP_max,1,LF2hInt(:,nR),map_glbl_st,.true.)
               end if

      !---------------------------------------------------------------------------------------------
      ! This is the original:
      !---------------------------------------------------------------------------------------------
               do lm=1,lm_max
                  Geo(lm)=CorPol(lm)-leg_helper%dpR(lm)+beta(nR)*leg_helper%preR(lm)
                  CLF(lm)=CorPol(lm)+LFPol(lm)
                  PLF(lm)=LFPol(lm)-leg_helper%dpR(lm)+beta(nR)*leg_helper%preR(lm)
                  Mag(lm)=Geo(lm)+LFPol(lm)
                  Arc(lm)=Geo(lm)+Buo(lm)
                  ArcMag(lm)=Mag(lm)+Buo(lm)
                  CIA(lm)=ArcMag(lm)+AdvPol(lm)
                  !CIA(lm)=CorPol(lm)+Buo(lm)+AdvPol(lm)
               end do
               call hIntRms(Geo,nR,1,lm_max,0,Geo2hInt(:,nR),map_glbl_st,.false.)
               call hIntRms(CLF,nR,1,lm_max,0,CLF2hInt(:,nR),map_glbl_st,.false.)
               call hIntRms(PLF,nR,1,lm_max,0,PLF2hInt(:,nR),map_glbl_st,.false.)
               call hIntRms(Mag,nR,1,lm_max,0,Mag2hInt(:,nR),map_glbl_st,.false.)
               call hIntRms(Arc,nR,1,lm_max,0,Arc2hInt(:,nR),map_glbl_st,.false.)
               call hIntRms(ArcMag,nR,1,lm_max,0,ArcMag2hInt(:,nR),map_glbl_st,.false.)
               call hIntRms(CIA,nR,1,lm_max,0,CIA2hInt(:,nR),map_glbl_st,.false.)
      !---------------------------------------------------------------------------------------------
      ! This is (more or less) how this specific piece will look like afterwards: 
      !---------------------------------------------------------------------------------------------
      !                do lm=1,n_lm_loc
      !                   l = map_dist_st%lm2l(lm)
      !                   m = map_dist_st%lm2m(lm)
      !                   lm_glb = lm2(l,m)
      !                   
      !                   Geo(lm)=CorPol(lm)-leg_helper%dpR(lm)+beta(nR)*leg_helper%preR(lm)
      !                   CLF(lm)=CorPol(lm)+LFPol(lm)
      !                   PLF(lm)=LFPol(lm)-leg_helper%dpR(lm)+beta(nR)*leg_helper%preR(lm)
      !                   Mag(lm)=Geo(lm)+LFPol(lm)
      !                   Arc(lm)=Geo(lm)+Buo(lm)
      !                   ArcMag(lm)=Mag(lm)+Buo(lm)
      !                   CIA(lm)=ArcMag(lm)+AdvPol(lm)
      !                   !CIA(lm)=CorPol(lm)+Buo(lm)+AdvPol(lm)
      !                end do
      !
      !                call hIntRms(Geo(1:n_lm_loc),nR,1,n_lm_loc,0,Geo2hInt(0:l_max,nR),map_dist_st,.false.)
      !                call hIntRms(CLF(1:n_lm_loc),nR,1,n_lm_loc,0,CLF2hInt(0:l_max,nR),map_dist_st,.false.)
      !                call hIntRms(PLF(1:n_lm_loc),nR,1,n_lm_loc,0,PLF2hInt(0:l_max,nR),map_dist_st,.false.)
      !                call hIntRms(Mag(1:n_lm_loc),nR,1,n_lm_loc,0,Mag2hInt(0:l_max,nR),map_dist_st,.false.)
      !                call hIntRms(Arc(1:n_lm_loc),nR,1,n_lm_loc,0,Arc2hInt(0:l_max,nR),map_dist_st,.false.)
      !                call hIntRms(ArcMag(1:n_lm_loc),nR,1,n_lm_loc,0,ArcMag2hInt(0:l_max,nR),map_dist_st,.false.)
      !                call hIntRms(CIA(1:n_lm_loc),nR,1,n_lm_loc,0,CIA2hInt(0:l_max,nR),map_dist_st,.false.)               
      !                
      !                call mpi_iallreduce(MPI_IN_PLACE, Geo2hInt(0:l_max,nR), l_max+1, MPI_DEF_REAL, MPI_SUM, comm_theta, Rq(1), ierr)
      !                call mpi_iallreduce(MPI_IN_PLACE, CLF2hInt(0:l_max,nR), l_max+1, MPI_DEF_REAL, MPI_SUM, comm_theta, Rq(2), ierr)
      !                call mpi_iallreduce(MPI_IN_PLACE, PLF2hInt(0:l_max,nR), l_max+1, MPI_DEF_REAL, MPI_SUM, comm_theta, Rq(3), ierr)
      !                call mpi_iallreduce(MPI_IN_PLACE, Mag2hInt(0:l_max,nR), l_max+1, MPI_DEF_REAL, MPI_SUM, comm_theta, Rq(4), ierr)
      !                call mpi_iallreduce(MPI_IN_PLACE, Arc2hInt(0:l_max,nR), l_max+1, MPI_DEF_REAL, MPI_SUM, comm_theta, Rq(5), ierr)
      !                call mpi_iallreduce(MPI_IN_PLACE, CIA2hInt(0:l_max,nR), l_max+1, MPI_DEF_REAL, MPI_SUM, comm_theta, Rq(6), ierr)
      !                call mpi_iallreduce(MPI_IN_PLACE, ArcMag2hInt(0:l_max,nR), l_max+1, MPI_DEF_REAL, MPI_SUM, comm_theta, Rq(7), ierr)
      !                
      !                call mpi_waitall(7, Rq(1:7), MPI_STATUSES_IGNORE, ierr )
      !                
      !                call gather_Flm(Geo   , Geo)
      !                call gather_Flm(CLF   , CLF)
      !                call gather_Flm(PLF   , PLF)
      !                call gather_Flm(Mag   , Mag)
      !                call gather_Flm(Arc   , Arc)
      !                call gather_Flm(ArcMag, ArcMag)
      !                call gather_Flm(CIA   , CIA)
      !---------------------------------------------------------------------------------------------
      ! Not very charming, but the gathers will drop (those are local variables)
      ! I am not sure yet, but maybe I can afford to do only a single mpi_iallreduce after ALL those 
      ! lm loops which appear next. I think it is just a sum, but I need to double check it.
      !---------------------------------------------------------------------------------------------

               do lm=1,lm_max
                  lmP = map_glbl_st%lm2lmP(lm)
                  Geo(lm)=-this%CFt2LM(lmP)-this%PFt2LM(lmP)
                  CLF(lm)=-this%CFt2LM(lmP)+this%LFt2LM(lmP)
                  PLF(lm)=this%LFt2LM(lmP)-this%PFt2LM(lmP)
                  Mag(lm)=Geo(lm)+this%LFt2LM(lmP)
                  Arc(lm)=Geo(lm)
                  ArcMag(lm)=Mag(lm)
                  CIA(lm)=ArcMag(lm)+this%Advt2LM(lmP)
                  !CIA(lm)=-this%CFt2LM(lmP)+this%Advt2LM(lmP)
               end do
               call hIntRms(Geo,nR,1,lm_max,0,Geo2hInt(:,nR),map_glbl_st,.true.)
               call hIntRms(CLF,nR,1,lm_max,0,CLF2hInt(:,nR),map_glbl_st,.true.)
               call hIntRms(PLF,nR,1,lm_max,0,PLF2hInt(:,nR),map_glbl_st,.true.)
               call hIntRms(Mag,nR,1,lm_max,0,Mag2hInt(:,nR),map_glbl_st,.true.)
               call hIntRms(Arc,nR,1,lm_max,0,Arc2hInt(:,nR),map_glbl_st,.true.)
               call hIntRms(ArcMag,nR,1,lm_max,0,ArcMag2hInt(:,nR),map_glbl_st,.true.)
               call hIntRms(CIA,nR,1,lm_max,0,CIA2hInt(:,nR),map_glbl_st,.true.)
    
               do lm=1,lm_max
                  lmP = map_glbl_st%lm2lmP(lm)
                  Geo(lm)=-this%CFp2LM(lmP)-this%PFp2LM(lmP)
                  CLF(lm)=-this%CFp2LM(lmP)+this%LFp2LM(lmP)
                  PLF(lm)=this%LFp2LM(lmP)-this%PFp2LM(lmP)
                  Mag(lm)=Geo(lm)+this%LFp2LM(lmP)
                  Arc(lm)=Geo(lm)
                  ArcMag(lm)=Mag(lm)
                  CIA(lm)=ArcMag(lm)+this%Advp2LM(lmP)
                  !CIA(lm)=-this%CFp2LM(lmP)+this%Advp2LM(lmP)
               end do
               call hIntRms(Geo,nR,1,lm_max,0,Geo2hInt(:,nR),map_glbl_st,.true.)
               call hIntRms(CLF,nR,1,lm_max,0,CLF2hInt(:,nR),map_glbl_st,.true.)
               call hIntRms(PLF,nR,1,lm_max,0,PLF2hInt(:,nR),map_glbl_st,.true.)
               call hIntRms(Mag,nR,1,lm_max,0,Mag2hInt(:,nR),map_glbl_st,.true.)
               call hIntRms(Arc,nR,1,lm_max,0,Arc2hInt(:,nR),map_glbl_st,.true.)
               call hIntRms(ArcMag,nR,1,lm_max,0,ArcMag2hInt(:,nR),map_glbl_st,.true.)
               call hIntRms(CIA,nR,1,lm_max,0,CIA2hInt(:,nR),map_glbl_st,.true.)

            end if
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!
!           END OF POSTPONED
!
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------

            ! In case double curl is calculated dpdt is useless
            if ( (.not. l_double_curl) .or. lPressCalc ) then 
               do lm=lm_maybe_skip_first,n_lm_loc
                  l   =map_dist_st%lm2l(lm)
                  m   =map_dist_st%lm2m(lm)
                  lmS =map_dist_st%lm2lmS(lm)
                  lmA =map_dist_st%lm2lmA(lm)
                  lmP =map_dist_st%lm2lmP(lm)
                  lmPS=map_dist_st%lmP2lmPS(lmP)
                  lmPA=map_dist_st%lmP2lmPA(lmP)
       
                  !------ Recycle CorPol and AdvPol:
                  if ( l_corr ) then
                     if ( l < l_max .and. l > m ) then
                        CorPol_loc=                    two*CorFac*or2(nR) *      &
                        &                    ( -dPhi0_loc(lm) * ( dw_Rdist(lm,nR) &
                        &                       +or1(nR)*leg_helper%dLhw(lm) &
                        &                                                  )     &
                        &                       +dTheta3A_loc(lm)*z_Rdist(lmA,nR) &
                        &                       +dTheta3S_loc(lm)*z_Rdist(lmS,nR) &
                        &                    )
       
                     else if ( l == l_max ) then
                        CorPol_loc=  two*CorFac*or2(nR) * ( -dPhi0_loc(lm) *  &
                                    ( dw_Rdist(lm,nR) + or1(nR)*leg_helper%dLhw(lm) ) )
       
                     else if ( l == m ) then
                        CorPol_loc=                    two*CorFac*or2(nR) *      &
                        &                    ( -dPhi0_loc(lm) * ( dw_Rdist(lm,nR) &
                        &                       +or1(nR)*leg_helper%dLhw(lm) &
                        &                                                   )    &
                        &                      +dTheta3A_loc(lm)*z_Rdist(lmA,nR)  &
                        &                    )
       
                     end if
                  else
                     CorPol_loc=zero
                  end if
                  if ( l_conv_nl ) then
                     if ( l > m ) then
                        AdvPol_loc= dTheta1S_loc(lm)*this%AdvtLM(lmPS) - &
                        &           dTheta1A_loc(lm)*this%AdvtLM(lmPA) + &
                        &               dPhi_loc(lm)*this%AdvpLM(lmP)
                     else if ( l == m ) then
                        AdvPol_loc=-dTheta1A_loc(lm)*this%AdvtLM(lmPA) + &
                        &               dPhi_loc(lm)*this%AdvpLM(lmP)
                     end if
                  else
                     AdvPol_loc=zero
                  end if
                  dpdt(lm)=AdvPol_loc+CorPol_loc
       
               end do ! lm loop
            end if
         else
            do lm=lm_maybe_skip_first,n_lm_loc
               dwdt(lm) =0.0_cp
               dzdt(lm) =0.0_cp
               dpdt(lm) =0.0_cp
            end do
         end if ! l_conv ?

      end if

      if ( nBc == 0 ) then
         
         if ( l_heat ) then
            
            if (map_dist_st%lm2(0,0) > 0) then  ! if m=0 is in this rank
               lm  = map_dist_st%lm2(0,0)

               dsdt_loc  =epsc*epscProf(nR) !+opr/epsS*divKtemp0(nR)
               dVSrLM(lm)=this%VSrLM(lm)
               if ( l_TP_form ) dVPrLM(lm)=this%VPrLM(lm)
               if ( l_anel ) then
                  if ( l_anelastic_liquid .or. l_TP_form ) then
                     if ( l_mag_nl ) then
                        dsdt_loc=dsdt_loc+                                        &
                        &    ViscHeatFac*hdif_V_loc(lm)*temp0(nR)*this%ViscHeatLM(lm)+  &
                        &     OhmLossFac*hdif_B_loc(lm)*temp0(nR)*this%OhmLossLM(lm)
                     else
                        dsdt_loc=dsdt_loc+ &
                        &    ViscHeatFac*hdif_V_loc(lm)*temp0(nR)*this%ViscHeatLM(lm)
                     end if
                  else
                     if ( l_mag_nl ) then
                        dsdt_loc=dsdt_loc+ViscHeatFac*hdif_V_loc(lm)*this%ViscHeatLM(lm)+ &
                        &                  OhmLossFac*hdif_B_loc(lm)*this%OhmLossLM(lm)
                     else
                        dsdt_loc=dsdt_loc+ViscHeatFac*hdif_V_loc(lm)*this%ViscHeatLM(lm)
                     end if
                  end if
               end if
               dsdt(lm)=dsdt_loc
            end if
    
            do lm=lm_maybe_skip_first,n_lm_loc
               l   =map_dist_st%lm2l(lm)
               m   =map_dist_st%lm2m(lm)
               lmP =map_dist_st%lm2lmP(lm)
               lmPS=map_dist_st%lmP2lmPS(lmP)
               lmPA=map_dist_st%lmP2lmPA(lmP)
               
               !------ This is horizontal heat advection:
    
               if ( l > m ) then
                  dsdt_loc= -dTheta1S_loc(lm)*this%VStLM(lmPS) &
                  &         +dTheta1A_loc(lm)*this%VStLM(lmPA) &
                  &         -dPhi_loc(lm)*this%VSpLM(lmP)
               else if ( l == m ) then
                  dsdt_loc=  dTheta1A_loc(lm)*this%VStLM(lmPA) &
                  &          -dPhi_loc(lm)*this%VSpLM(lmP)
               end if
               if ( l_anel ) then
                  if ( l_anelastic_liquid .or. l_TP_form ) then
                     dsdt_loc = dsdt_loc+ &
                     &          ViscHeatFac*hdif_V_loc(lm)*temp0(nR)*this%ViscHeatLM(lmP)
                     if ( l_mag_nl ) then
                        dsdt_loc = dsdt_loc+ &
                        &          OhmLossFac*hdif_B_loc(lm)*temp0(nR)*this%OhmLossLM(lmP)
                     end if
                  else
                     dsdt_loc = dsdt_loc+ &
                     &          ViscHeatFac*hdif_V_loc(lm)*this%ViscHeatLM(lmP)
                     if ( l_mag_nl ) then
                        dsdt_loc = dsdt_loc+ &
                        &          OhmLossFac*hdif_B_loc(lm)*this%OhmLossLM(lmP)
                     end if
                  end if
               end if
               !-----   simplified form for linear onset !
               !        not ds not saved in the current program form!
               !                 dsdt(lm)=
               !                    -dLh(lm)*w(lm,nR)*or2(nR)*dsR(1)
               dVSrLM(lm)=this%VSrLM(lmP)
               dsdt(lm) = dsdt_loc
               if ( l_TP_form ) dVPrLM(lm)=this%VPrLM(lmP)
            end do
         else
            
            do lm=lm_maybe_skip_first,n_lm_loc
               dsdt(lm)  =0.0_cp
               dVSrLM(lm)=0.0_cp
            end do
         end if

         if ( l_chemical_conv ) then
            if (map_dist_st%lm2(0,0) > 0) then  ! if m=0 is in this rank
               lm  = map_dist_st%lm2(0,0)
               dVXirLM(lm)=this%VXirLM(lm)
               dxidt(lm)  =epscXi
            end if
    
            do lm=lm_maybe_skip_first,n_lm_loc
               l   =map_dist_st%lm2l(lm)
               m   =map_dist_st%lm2m(lm)
               lmP =map_dist_st%lm2lmP(lm)
               lmPS=map_dist_st%lmP2lmPS(lmP)
               lmPA=map_dist_st%lmP2lmPA(lmP)
               !------ This is horizontal heat advection:
    
               if ( l > m ) then
                  dxidt_loc= -dTheta1S_loc(lm)*this%VXitLM(lmPS) &
                  &          +dTheta1A_loc(lm)*this%VXitLM(lmPA) &
                  &          -dPhi_loc(lm)*this%VXipLM(lmP)
               else if ( l == m ) then
                  dxidt_loc=  dTheta1A_loc(lm)*this%VXitLM(lmPA) &
                  &          -dPhi_loc(lm)*this%VXipLM(lmP)
               end if
               dVXirLM(lm)=this%VXirLM(lmP)
               dxidt(lm)  =dxidt_loc
            end do
         end if
    
         if ( l_mag_nl .or. l_mag_kin  ) then
            do lm=1,n_lm_loc
               l   =map_dist_st%lm2l(lm)
               m   =map_dist_st%lm2m(lm)
               lmP =map_dist_st%lm2lmP(lm)
               lmPS=map_dist_st%lmP2lmPS(lmP)
               lmPA=map_dist_st%lmP2lmPA(lmP)
               
               if ((l == 0) .and. (m == 0)) then
                  dVxBhLM(lm)= -r(nR)*r(nR)* dTheta1A_loc(lm)*this%VxBtLM(lmPA)
                  dbdt(lm)   = -dTheta1A_loc(lm)*this%VxBpLM(lmPA)
                  djdt(lm)   = zero
                  cycle ! <---------------------------

               !------- This is the radial part of the dynamo terms \curl(VxB)
               else if ( l > m ) then

                  dbdt(lm)=  dTheta1S_loc(lm)*this%VxBpLM(lmPS) &
                  &         -dTheta1A_loc(lm)*this%VxBpLM(lmPA) &
                  &         -dPhi_loc(lm)    *this%VxBtLM(lmP)

               else if ( l == m ) then
                  dbdt(lm)= -dTheta1A_loc(lm)*this%VxBpLM(lmPA) &
                  &         -dPhi_loc(lm)    *this%VxBtLM(lmP)
               end if
    
               !------- Radial component of
               !           \curl\curl(UxB) = \grad\div(UxB) - \laplace(VxB)
    
               !------- This is the radial part of \laplace (UxB)
               djdt(lm)=dLh_loc(lm)*or4(nR)*this%VxBrLM(lmP)
    
               !------- This is r^2 * horizontal divergence of (UxB)
               !        Radial derivative performed in get_dr_td
               if ( l > m ) then
                  dVxBhLM(lm)=            r(nR)*r(nR)* ( &
                  &    dTheta1S_loc(lm)*this%VxBtLM(lmPS) -  &
                  &    dTheta1A_loc(lm)*this%VxBtLM(lmPA) +  &
                  &    dPhi_loc(lm)*this%VxBpLM(lmP)  )
               else if ( l == m ) then
                  dVxBhLM(lm)=              r(nR)*r(nR)* ( &
                  &    - dTheta1A_loc(lm)*this%VxBtLM(lmPA) +  &
                  &    dPhi_loc(lm)*this%VxBpLM(lmP)  )
               end if
            end do
         else if ( l_mag ) then
            do lm=1,n_lm_loc
               dbdt(lm)   =zero
               djdt(lm)   =zero
               dVxBhLM(lm)=zero
            end do
         end if

      else   ! boundary !

         if ( l_mag_nl .or. l_mag_kin ) then
    
            !----- Stress free boundary, only nl mag. term for poloidal field needed.
            !      Because the radial derivative will be taken, this will contribute to
            !      the other radial grid points.
            do lm=1,n_lm_loc
               l   = map_dist_st%lm2l(lm)
               m   = map_dist_st%lm2m(lm)
               lmP = map_dist_st%lm2lmP(lm)
               lmPS= map_dist_st%lmP2lmPS(lmP)   ! l-1
               lmPA= map_dist_st%lmP2lmPA(lmP)   ! l+1
               
               if ((l == 0) .and. (m == 0)) then
                  dVxBhLM(lm)=zero
                  dVSrLM(lm) =zero
                  cycle ! <---------------------------
               else if ( l > m ) then
                  dVxBhLM(lm)=r(nR)*r(nR)* (               &
                  &      dTheta1S_loc(lm)*this%VxBtLM(lmPS) -  &
                  &      dTheta1A_loc(lm)*this%VxBtLM(lmPA) +  &
                  &          dPhi_loc(lm)*this%VxBpLM(lmP)  )
               else if ( l == m ) then ! (l-1) not allowed !
                  dVxBhLM(lm)=r(nR)*r(nR)* (               &
                  &    - dTheta1A_loc(lm)*this%VxBtLM(lmPA) +  &
                  &    dPhi_loc(lm)*this%VxBpLM(lmP)  )
               end if
               dVSrLM(lm)=zero
            end do
    
         else
            do lm=1,n_lm_loc
               if ( l_mag ) dVxBhLM(lm)=zero
               dVSrLM(lm) =zero
            end do
         end if
         if ( l_double_curl ) then
            do lm=1,n_lm_loc
               dVxVhLM(lm)=zero
            end do
         end if
         if ( l_chemical_conv ) then
            do lm=1,n_lm_loc
               dVXirLM(lm)=zero
            end do
         end if
         if ( l_TP_form ) then
            do lm=1,n_lm_loc
               dVPrLM(lm)=zero
            end do
         end if
         
      end if  ! boundary ? lvelo ?
      
   end subroutine get_td
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
   subroutine slice_all_nonlinear_lm(this, nl_lm_glb)
      !
      !@>author Rafael Lago, MPCDF, December 2017
      !
      class(nonlinear_lm_t)      :: this, nl_lm_glb
      
      call slice_FlmP(nl_lm_glb%AdvrLM, this%AdvrLM)   
      call slice_FlmP(nl_lm_glb%AdvtLM, this%AdvtLM)   
      call slice_FlmP(nl_lm_glb%AdvpLM, this%AdvpLM)   
      call slice_FlmP(nl_lm_glb%LFrLM, this%LFrLM)    
      call slice_FlmP(nl_lm_glb%LFtLM, this%LFtLM)    
      call slice_FlmP(nl_lm_glb%LFpLM, this%LFpLM)    
      call slice_FlmP(nl_lm_glb%VxBrLM, this%VxBrLM)   
      call slice_FlmP(nl_lm_glb%VxBtLM, this%VxBtLM)   
      call slice_FlmP(nl_lm_glb%VxBpLM, this%VxBpLM)   
      call slice_FlmP(nl_lm_glb%VSrLM, this%VSrLM)    
      call slice_FlmP(nl_lm_glb%VStLM, this%VStLM)    
      call slice_FlmP(nl_lm_glb%VSpLM, this%VSpLM)    
      call slice_FlmP(nl_lm_glb%ViscHeatLM, this%ViscHeatLM)
      call slice_FlmP(nl_lm_glb%OhmLossLM, this%OhmLossLM)

      if ( l_TP_form ) then
         call slice_FlmP(nl_lm_glb%VPrLM, this%VPrLM)    
      end if
      
      if ( l_chemical_conv ) then
         call slice_FlmP(nl_lm_glb%VXirLM, this%VXirLM)    
         call slice_FlmP(nl_lm_glb%VXitLM, this%VXitLM)    
         call slice_FlmP(nl_lm_glb%VXipLM, this%VXipLM)    
      end if

      !-- RMS calculations
      if ( l_RMS ) then
         call slice_FlmP(nl_lm_glb%Advt2LM, this%Advt2LM)
         call slice_FlmP(nl_lm_glb%Advp2LM, this%Advp2LM)
         call slice_FlmP(nl_lm_glb%LFt2LM, this%LFt2LM)
         call slice_FlmP(nl_lm_glb%LFp2LM, this%LFp2LM)
         call slice_FlmP(nl_lm_glb%CFt2LM, this%CFt2LM)
         call slice_FlmP(nl_lm_glb%CFp2LM, this%CFp2LM)
         call slice_FlmP(nl_lm_glb%PFt2LM, this%PFt2LM)
         call slice_FlmP(nl_lm_glb%PFp2LM, this%PFp2LM)
      end if
      
   end subroutine slice_all_nonlinear_lm
!------------------------------------------------------------------------------
   subroutine gather_all_nonlinear_lm(this, nl_lm_glb)
      !
      !@>author Rafael Lago, MPCDF, December 2017
      !
      class(nonlinear_lm_t)      :: this, nl_lm_glb
      
      call gather_FlmP(this%AdvrLM, nl_lm_glb%AdvrLM)   
      call gather_FlmP(this%AdvtLM, nl_lm_glb%AdvtLM)   
      call gather_FlmP(this%AdvpLM, nl_lm_glb%AdvpLM)   
      call gather_FlmP(this%LFrLM, nl_lm_glb%LFrLM)    
      call gather_FlmP(this%LFtLM, nl_lm_glb%LFtLM)    
      call gather_FlmP(this%LFpLM, nl_lm_glb%LFpLM)    
      call gather_FlmP(this%VxBrLM, nl_lm_glb%VxBrLM)   
      call gather_FlmP(this%VxBtLM, nl_lm_glb%VxBtLM)   
      call gather_FlmP(this%VxBpLM, nl_lm_glb%VxBpLM)   
      call gather_FlmP(this%VSrLM, nl_lm_glb%VSrLM)    
      call gather_FlmP(this%VStLM, nl_lm_glb%VStLM)    
      call gather_FlmP(this%VSpLM, nl_lm_glb%VSpLM)    
      call gather_FlmP(this%ViscHeatLM, nl_lm_glb%ViscHeatLM)
      call gather_FlmP(this%OhmLossLM, nl_lm_glb%OhmLossLM)

      if ( l_TP_form ) then
         call gather_FlmP(this%VPrLM, nl_lm_glb%VPrLM)    
      end if
      
      if ( l_chemical_conv ) then
         call gather_FlmP(this%VXirLM, nl_lm_glb%VXirLM)    
         call gather_FlmP(this%VXitLM, nl_lm_glb%VXitLM)    
         call gather_FlmP(this%VXipLM, nl_lm_glb%VXipLM)    
      end if

      !-- RMS calculations
      if ( l_RMS ) then
         call gather_FlmP(this%Advt2LM, nl_lm_glb%Advt2LM)
         call gather_FlmP(this%Advp2LM, nl_lm_glb%Advp2LM)
         call gather_FlmP(this%LFt2LM, nl_lm_glb%LFt2LM)
         call gather_FlmP(this%LFp2LM, nl_lm_glb%LFp2LM)
         call gather_FlmP(this%CFt2LM, nl_lm_glb%CFt2LM)
         call gather_FlmP(this%CFp2LM, nl_lm_glb%CFp2LM)
         call gather_FlmP(this%PFt2LM, nl_lm_glb%PFt2LM)
         call gather_FlmP(this%PFp2LM, nl_lm_glb%PFp2LM)
      end if
      
   end subroutine gather_all_nonlinear_lm
!------------------------------------------------------------------------------
end module nonlinear_lm_mod
