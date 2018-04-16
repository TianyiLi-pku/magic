!-------------------------------------------------------------------------------
module arrays_dist
!@>details Distributed implementation of grid_space_arrays_t and nonlinear_lm_t
!> for coord_theta>1.
!
!@>author Rafael Lago, MPCDF, December 2017
!-------------------------------------------------------------------------------
   use grid_space_arrays_mod, only: grid_space_arrays_t
   use mem_alloc, only: bytes_allocated
   use truncation
   use nonlinear_lm_mod, only: nonlinear_lm_t
   use logic, only: l_conv_nl, l_heat_nl, l_mag_nl, l_anel, l_mag_LF, &
       &            l_RMS, l_chemical_conv, l_TP_form, l_HT, l_heat, l_conv, &
       &            l_mag_kin, l_precession, l_mag, l_anelastic_liquid, &
       &            l_corr, l_double_curl, l_single_matrix
       
   use fields, only: s_Rloc,ds_Rloc, z_Rloc,dz_Rloc, p_Rloc,dp_Rloc, &
       &             b_Rloc,db_Rloc,ddb_Rloc, aj_Rloc,dj_Rloc,       &
       &             w_Rloc,dw_Rloc,ddw_Rloc, xi_Rloc,               &
       &             s_dist,ds_dist, z_dist,dz_dist, p_dist,dp_dist, &
       &             b_dist,db_dist,ddb_dist, aj_dist,dj_dist,       &
       &             w_dist,dw_dist,ddw_dist, xi_dist
   use precision_mod
   use radial_data, only: nRstart, nRstop
   
   use radial_functions, only: r, or2, or1, beta, rho0, rgrav, epscProf, &
       &                       or4, temp0, alpha0, ogrun, orho1
   use physical_parameters, only: CorFac, ra, epsc, ViscHeatFac, &
       &                          OhmLossFac, n_r_LCR, epscXi,   &
       &                          BuoFac, ThExpNb
   use blocking, only: lm2l, lm2m, lm2lmP, lmP2lmPS, lmP2lmPA, lm2lmA, &
       &               lm2lmS, st_map, lm2
   use horizontal_data, only: dLh, dTheta1S, dTheta1A, dPhi, dTheta2A, &
       &                      dTheta3A, dTheta4A, dPhi0, dTheta2S,     &
       &                      dTheta3S, dTheta4S, hdif_V, hdif_B,      &
       &                      dTheta1A_loc, dTheta2A_loc, dTheta3A_loc,&
       &                      dTheta4A_loc, dTheta1S_loc, dTheta2S_loc,&
       &                      dTheta3S_loc, dTheta4S_loc, dLH_loc,     &
       &                      dPhi_loc, dPhi0_loc, hdif_B_loc, &
       &                      hdif_V_loc, hdif_S_loc, hdif_Xi_loc
       
   use RMS, only: Adv2hInt, Pre2hInt, Buo2hInt, Cor2hInt, LF2hInt,  &
       &          Geo2hInt, Mag2hInt, ArcMag2hInt, CLF2hInt, PLF2hInt, &
       &          CIA2hInt, Arc2hInt
   use leg_helper_mod, only: leg_helper_t
   use constants, only: zero, two
   use RMS_helpers, only: hIntRms
   use distributed_theta, only: dist_map

   implicit none

   private
   
   type, public, extends(grid_space_arrays_t) :: grid_space_arrays_dist_t
   contains
      procedure :: initialize => initialize_grid_space             ! override
      procedure :: slice_all  => slice_all_grid_space
      procedure :: gather_all  => gather_all_grid_space
   end type grid_space_arrays_dist_t
   
   type, public, extends(nonlinear_lm_t) :: nonlinear_lm_dist_t
   contains
      procedure :: initialize => initialize_nonlinear_lm           ! override
      procedure :: get_td_dist
      procedure :: slice_all  => slice_all_nonlinear_lm
      procedure :: gather_all  => gather_all_nonlinear_lm
   end type nonlinear_lm_dist_t

contains

!-------------------------------------------------------------------------------
   subroutine gather_all_grid_space(this, gsa_glb)
!@>author Rafael Lago, MPCDF, December 2017
!-------------------------------------------------------------------------------
      class(grid_space_arrays_dist_t) :: this
      class(grid_space_arrays_t)      :: gsa_glb
      
      call gather_f(this%Advr, gsa_glb%Advr)
      call gather_f(this%Advt, gsa_glb%Advt)
      call gather_f(this%Advp, gsa_glb%Advp)
      call gather_f(this%LFr, gsa_glb%LFr)
      call gather_f(this%LFt, gsa_glb%LFt)
      call gather_f(this%LFp, gsa_glb%LFp)
      call gather_f(this%VxBr, gsa_glb%VxBr)
      call gather_f(this%VxBt, gsa_glb%VxBt)
      call gather_f(this%VxBp, gsa_glb%VxBp)
      call gather_f(this%VSr, gsa_glb%VSr)
      call gather_f(this%VSt, gsa_glb%VSt)
      call gather_f(this%VSp, gsa_glb%VSp)
      call gather_f(this%ViscHeat, gsa_glb%ViscHeat)
      call gather_f(this%OhmLoss, gsa_glb%OhmLoss)

      if ( l_TP_form ) then
         call gather_f(this%VPr, gsa_glb%VPr)
      end if
      
      if ( l_precession ) then
         call gather_f(this%PCr, gsa_glb%PCr)
         call gather_f(this%PCt, gsa_glb%PCt)
         call gather_f(this%PCp, gsa_glb%PCp)
      end if

      if ( l_chemical_conv ) then
         call gather_f(this%VXir, gsa_glb%VXir)
         call gather_f(this%VXit, gsa_glb%VXit)
         call gather_f(this%VXip, gsa_glb%VXip)
      end if

      !----- Fields calculated from these help arrays by legtf:
      call gather_f(this%vrc, gsa_glb%vrc)
      call gather_f(this%vtc, gsa_glb%vtc)
      call gather_f(this%vpc, gsa_glb%vpc)
      call gather_f(this%dvrdrc, gsa_glb%dvrdrc)
      call gather_f(this%dvtdrc, gsa_glb%dvtdrc)
      call gather_f(this%dvpdrc, gsa_glb%dvpdrc)
      call gather_f(this%cvrc, gsa_glb%cvrc)
      call gather_f(this%dvrdtc, gsa_glb%dvrdtc)
      call gather_f(this%dvrdpc, gsa_glb%dvrdpc)
      call gather_f(this%dvtdpc, gsa_glb%dvtdpc)
      call gather_f(this%dvpdpc, gsa_glb%dvpdpc)
      call gather_f(this%brc, gsa_glb%brc)
      call gather_f(this%btc, gsa_glb%btc)
      call gather_f(this%bpc, gsa_glb%bpc)
      call gather_f(this%cbrc, gsa_glb%cbrc)
      call gather_f(this%cbtc, gsa_glb%cbtc)
      call gather_f(this%cbpc, gsa_glb%cbpc)
      call gather_f(this%sc, gsa_glb%sc)
      call gather_f(this%drSc, gsa_glb%drSc)
      call gather_f(this%pc, gsa_glb%pc)
      call gather_f(this%dsdtc, gsa_glb%dsdtc)
      call gather_f(this%dsdpc, gsa_glb%dsdpc)

      if ( l_chemical_conv ) then
         call gather_f(this%xic, gsa_glb%xic)
      else
         gsa_glb%xic(1,1) = this%xic(1,1)
      end if

      !-- RMS Calculations
      if ( l_RMS ) then
         call gather_f(this%Advt2, gsa_glb%Advt2)
         call gather_f(this%Advp2, gsa_glb%Advp2)
         call gather_f(this%LFt2, gsa_glb%LFt2)
         call gather_f(this%LFp2, gsa_glb%LFp2)
         call gather_f(this%CFt2, gsa_glb%CFt2)
         call gather_f(this%CFp2, gsa_glb%CFp2)
         call gather_f(this%dpdtc, gsa_glb%dpdtc)
         call gather_f(this%dpdpc, gsa_glb%dpdpc)
      end if
      
   end subroutine gather_all_grid_space

!-------------------------------------------------------------------------------
   subroutine gather_all_nonlinear_lm(this, nl_lm_glb)
!@>author Rafael Lago, MPCDF, December 2017
!-------------------------------------------------------------------------------
      class(nonlinear_lm_dist_t) :: this
      class(nonlinear_lm_t)      :: nl_lm_glb
      
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

!-------------------------------------------------------------------------------
   subroutine slice_all_grid_space(this, gsa_glb)
!@>author Rafael Lago, MPCDF, December 2017
!-------------------------------------------------------------------------------
      class(grid_space_arrays_dist_t) :: this
      class(grid_space_arrays_t)      :: gsa_glb
      
      call slice_f(gsa_glb%Advr, this%Advr)
      call slice_f(gsa_glb%Advt, this%Advt)
      call slice_f(gsa_glb%Advp, this%Advp)
      call slice_f(gsa_glb%LFr, this%LFr)
      call slice_f(gsa_glb%LFt, this%LFt)
      call slice_f(gsa_glb%LFp, this%LFp)
      call slice_f(gsa_glb%VxBr, this%VxBr)
      call slice_f(gsa_glb%VxBt, this%VxBt)
      call slice_f(gsa_glb%VxBp, this%VxBp)
      call slice_f(gsa_glb%VSr, this%VSr)
      call slice_f(gsa_glb%VSt, this%VSt)
      call slice_f(gsa_glb%VSp, this%VSp)
      call slice_f(gsa_glb%ViscHeat, this%ViscHeat)
      call slice_f(gsa_glb%OhmLoss, this%OhmLoss)

      if ( l_TP_form ) then
         call slice_f(gsa_glb%VPr, this%VPr)
      end if
      
      if ( l_precession ) then
         call slice_f(gsa_glb%PCr, this%PCr)
         call slice_f(gsa_glb%PCt, this%PCt)
         call slice_f(gsa_glb%PCp, this%PCp)
      end if

      if ( l_chemical_conv ) then
         call slice_f(gsa_glb%VXir, this%VXir)
         call slice_f(gsa_glb%VXit, this%VXit)
         call slice_f(gsa_glb%VXip, this%VXip)
      end if

      !----- Fields calculated from these help arrays by legtf:
      call slice_f(gsa_glb%vrc, this%vrc )
      call slice_f(gsa_glb%vtc, this%vtc)
      call slice_f(gsa_glb%vpc, this%vpc)
      call slice_f(gsa_glb%dvrdrc, this%dvrdrc)
      call slice_f(gsa_glb%dvtdrc, this%dvtdrc)
      call slice_f(gsa_glb%dvpdrc, this%dvpdrc)
      call slice_f(gsa_glb%cvrc, this%cvrc)
      call slice_f(gsa_glb%dvrdtc, this%dvrdtc)
      call slice_f(gsa_glb%dvrdpc, this%dvrdpc)
      call slice_f(gsa_glb%dvtdpc, this%dvtdpc)
      call slice_f(gsa_glb%dvpdpc, this%dvpdpc)
      call slice_f(gsa_glb%brc, this%brc)
      call slice_f(gsa_glb%btc, this%btc)
      call slice_f(gsa_glb%bpc, this%bpc)
      call slice_f(gsa_glb%cbrc, this%cbrc)
      call slice_f(gsa_glb%cbtc, this%cbtc)
      call slice_f(gsa_glb%cbpc, this%cbpc)
      call slice_f(gsa_glb%sc, this%sc)
      call slice_f(gsa_glb%drSc, this%drSc)
      call slice_f(gsa_glb%pc, this%pc)
      call slice_f(gsa_glb%dsdtc, this%dsdtc)
      call slice_f(gsa_glb%dsdpc, this%dsdpc)

      if ( l_chemical_conv ) then
         call slice_f(gsa_glb%xic, this%xic)
      else
         this%xic(1,1) = gsa_glb%xic(1,1)
      end if

      !-- RMS Calculations
      if ( l_RMS ) then
         call slice_f(gsa_glb%Advt2, this%Advt2)
         call slice_f(gsa_glb%Advp2, this%Advp2)
         call slice_f(gsa_glb%LFt2, this%LFt2)
         call slice_f(gsa_glb%LFp2, this%LFp2)
         call slice_f(gsa_glb%CFt2, this%CFt2)
         call slice_f(gsa_glb%CFp2, this%CFp2)
         call slice_f(gsa_glb%dpdtc, this%dpdtc)
         call slice_f(gsa_glb%dpdpc, this%dpdpc)
      end if
      
   end subroutine slice_all_grid_space

!-------------------------------------------------------------------------------
   subroutine slice_all_nonlinear_lm(this, nl_lm_glb)
!@>author Rafael Lago, MPCDF, December 2017
!-------------------------------------------------------------------------------
      class(nonlinear_lm_dist_t) :: this
      class(nonlinear_lm_t)      :: nl_lm_glb
      
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

!-------------------------------------------------------------------------------   
   subroutine initialize_grid_space(this)
!@>author Rafael Lago, MPCDF, December 2017
!-------------------------------------------------------------------------------
      class(grid_space_arrays_dist_t) :: this
      
      allocate( this%Advr(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%Advt(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%Advp(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%LFr(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%LFt(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%LFp(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%VxBr(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%VxBt(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%VxBp(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%VSr(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%VSt(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%VSp(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%ViscHeat(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%OhmLoss(n_phi_max,n_theta_beg:n_theta_end) )
      bytes_allocated=bytes_allocated + 14*n_phi_max*n_theta_loc*SIZEOF_DEF_REAL

      if ( l_TP_form ) then
         allocate( this%VPr(n_phi_max,n_theta_beg:n_theta_end) )
         bytes_allocated=bytes_allocated + n_phi_max*n_theta_loc*SIZEOF_DEF_REAL
      end if
      
      if ( l_precession ) then
         allocate( this%PCr(n_phi_max,n_theta_beg:n_theta_end) )
         allocate( this%PCt(n_phi_max,n_theta_beg:n_theta_end) )
         allocate( this%PCp(n_phi_max,n_theta_beg:n_theta_end) )
         bytes_allocated=bytes_allocated + 3*n_phi_max*n_theta_loc*SIZEOF_DEF_REAL
      end if

      if ( l_chemical_conv ) then
         allocate( this%VXir(n_phi_max,n_theta_beg:n_theta_end) )
         allocate( this%VXit(n_phi_max,n_theta_beg:n_theta_end) )
         allocate( this%VXip(n_phi_max,n_theta_beg:n_theta_end) )
         bytes_allocated=bytes_allocated + 3*n_phi_max*n_theta_loc*SIZEOF_DEF_REAL
      end if

      !----- Fields calculated from these help arrays by legtf:
      allocate( this%vrc(n_phi_max,n_theta_beg:n_theta_end) ) 
      allocate( this%vtc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%vpc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%dvrdrc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%dvtdrc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%dvpdrc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%cvrc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%dvrdtc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%dvrdpc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%dvtdpc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%dvpdpc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%brc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%btc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%bpc(n_phi_max,n_theta_beg:n_theta_end) )
      this%btc=1.0e50_cp
      this%bpc=1.0e50_cp
      allocate( this%cbrc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%cbtc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%cbpc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%sc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%drSc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%pc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%dsdtc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%dsdpc(n_phi_max,n_theta_beg:n_theta_end) )
      bytes_allocated=bytes_allocated + 22*n_phi_max*n_theta_loc*SIZEOF_DEF_REAL

      if ( l_chemical_conv ) then
         allocate( this%xic(n_phi_max,n_theta_beg:n_theta_end) )
         bytes_allocated=bytes_allocated + n_phi_max*n_theta_loc*SIZEOF_DEF_REAL
      else
         allocate( this%xic(1,1) )
      end if

      !-- RMS Calculations
      if ( l_RMS ) then
         allocate ( this%Advt2(n_phi_max,n_theta_beg:n_theta_end) )
         allocate ( this%Advp2(n_phi_max,n_theta_beg:n_theta_end) )
         allocate ( this%LFt2(n_phi_max,n_theta_beg:n_theta_end) )
         allocate ( this%LFp2(n_phi_max,n_theta_beg:n_theta_end) )
         allocate ( this%CFt2(n_phi_max,n_theta_beg:n_theta_end) )
         allocate ( this%CFp2(n_phi_max,n_theta_beg:n_theta_end) )
         allocate ( this%dpdtc(n_phi_max,n_theta_beg:n_theta_end) )
         allocate ( this%dpdpc(n_phi_max,n_theta_beg:n_theta_end) )
         bytes_allocated=bytes_allocated + 8*n_phi_max*n_theta_loc*SIZEOF_DEF_REAL
      end if
      !write(*,"(A,I15,A)") "grid_space_arrays: allocated ",bytes_allocated,"B."

   end subroutine initialize_grid_space
   
!-------------------------------------------------------------------------------
   subroutine initialize_nonlinear_lm(this)
!@>author Rafael Lago, MPCDF, December 2017
!-------------------------------------------------------------------------------

      class(nonlinear_lm_dist_t) :: this

      allocate( this%AdvrLM(lmP_loc) )   
      allocate( this%AdvtLM(lmP_loc) )   
      allocate( this%AdvpLM(lmP_loc) )   
      allocate( this%LFrLM(lmP_loc) )    
      allocate( this%LFtLM(lmP_loc) )    
      allocate( this%LFpLM(lmP_loc) )    
      allocate( this%VxBrLM(lmP_loc) )   
      allocate( this%VxBtLM(lmP_loc) )   
      allocate( this%VxBpLM(lmP_loc) )   
      allocate( this%VSrLM(lmP_loc) )    
      allocate( this%VStLM(lmP_loc) )    
      allocate( this%VSpLM(lmP_loc) )    
      allocate( this%ViscHeatLM(lmP_loc) )
      allocate( this%OhmLossLM(lmP_loc) )
      bytes_allocated = bytes_allocated + 14*lmP_loc*SIZEOF_DEF_COMPLEX

      if ( l_TP_form ) then
         allocate( this%VPrLM(lmP_loc) )    
         bytes_allocated = bytes_allocated + lmP_loc*SIZEOF_DEF_COMPLEX
      end if
      
      if ( l_chemical_conv ) then
         allocate( this%VXirLM(lmP_loc) )    
         allocate( this%VXitLM(lmP_loc) )    
         allocate( this%VXipLM(lmP_loc) )    
         bytes_allocated = bytes_allocated + 3*lmP_loc*SIZEOF_DEF_COMPLEX
      end if

      !-- RMS calculations
      if ( l_RMS ) then
         allocate( this%Advt2LM(lmP_loc) )
         allocate( this%Advp2LM(lmP_loc) )
         allocate( this%LFt2LM(lmP_loc) )
         allocate( this%LFp2LM(lmP_loc) )
         allocate( this%CFt2LM(lmP_loc) )
         allocate( this%CFp2LM(lmP_loc) )
         allocate( this%PFt2LM(lmP_loc) )
         allocate( this%PFp2LM(lmP_loc) )
         bytes_allocated = bytes_allocated + 8*lmP_loc*SIZEOF_DEF_COMPLEX
      end if

   end subroutine initialize_nonlinear_lm
!----------------------------------------------------------------------------
   subroutine get_td_dist(this_dist,nR,nBc,lRmsCalc,lPressCalc,dVSrLM,dVPrLM,dVXirLM, &
              &      dVxVhLM,dVxBhLM,dwdt,dzdt,dpdt,dsdt,dxidt,dbdt,djdt,   &
              &      leg_helper,this)
      !
      !  Purpose of this to calculate time derivatives
      !  dwdt,dzdt,dpdt,dsdt,dxidt,dbdt,djdt
      !  and auxiliary arrays dVSrLM, dVXirLM and dVxBhLM, dVxVhLM
      !  from non-linear terms in spectral form,
      !  contained in flmw1-3,flms1-3, flmb1-3 (input)
      !
    
      !-- Input of variables:
      class(nonlinear_lm_dist_t) :: this_dist
      class(nonlinear_lm_t) :: this

      integer,            intent(in) :: nR
      integer,            intent(in) :: nBc ! signifies boundary conditions
      logical,            intent(in) :: lRmsCalc
      logical,            intent(in) :: lPressCalc
      type(leg_helper_t), intent(in) :: leg_helper
    
      !-- Output of variables:
      complex(cp), intent(out) :: dwdt(lm_max),dzdt(lm_max)
      complex(cp), intent(out) :: dpdt(lm_max),dsdt(lm_max)
      complex(cp), intent(out) :: dxidt(lm_max)
      complex(cp), intent(out) :: dbdt(lm_maxMag),djdt(lm_maxMag)
      complex(cp), intent(out) :: dVxBhLM(lm_maxMag)
      complex(cp), intent(out) :: dVxVhLM(lm_max)
      complex(cp), intent(out) :: dVSrLM(lm_max)
      complex(cp), intent(out) :: dVXirLM(lm_max)
      complex(cp), intent(out) :: dVPrLM(lm_max)
    
      !-- Local variables:
      integer :: l,m,lm,lmS,lmA,lmP,lmPS,lmPA
      complex(cp) :: CorPol(lm_max)
      complex(cp) :: AdvPol(lm_max),AdvTor(lm_max)
      complex(cp) :: LFPol(lm_max),LFTor(lm_max)
      complex(cp) :: Geo(lm_max),CLF(lm_max),PLF(lm_max)
      complex(cp) :: ArcMag(lm_max),Mag(lm_max),CIA(lm_max),Arc(lm_max)
      complex(cp) :: Buo(lm_max)
      complex(cp) :: AdvPol_loc,CorPol_loc,AdvTor_loc,CorTor_loc
      complex(cp) :: dsdt_loc, dxidt_loc
    
      integer, parameter :: DOUBLE_COMPLEX_PER_CACHELINE=4
      
      !--- Duplications - yay -.-"
      complex(cp) :: CorPol_dist(lm_loc)
      complex(cp) :: AdvPol_dist(lm_loc),AdvTor_dist(lm_loc)
      complex(cp) :: LFPol_dist(lm_loc),LFTor_dist(lm_loc)
      complex(cp) :: Geo_dist(lm_loc),CLF_dist(lm_loc),PLF_dist(lm_loc)
      complex(cp) :: ArcMag_dist(lm_loc),Mag_dist(lm_loc),CIA_dist(lm_loc),Arc_dist(lm_loc)
      complex(cp) :: Buo_dist(lm_loc)
      
      complex(cp) :: dwdt_dist(lm_loc),dzdt_dist(lm_loc)
      complex(cp) :: dpdt_dist(lm_loc),dsdt_dist(lm_loc)
      complex(cp) :: dxidt_dist(lm_loc)
      complex(cp) :: dbdt_dist(lm_loc),djdt_dist(lm_loc)
      complex(cp) :: dVxBhLM_dist(lm_loc)
      complex(cp) :: dVxVhLM_dist(lm_loc)
      complex(cp) :: dVSrLM_dist(lm_loc)
      complex(cp) :: dVXirLM_dist(lm_loc)
      complex(cp) :: dVPrLM_dist(lm_loc)    
      
      integer :: Rq(10)
      integer :: lm_maybe_skip_first
      integer :: lm_dist
!       integer :: lmS_dist 
!       integer :: lmA_dist 
!       integer :: lmP_dist 
!       integer :: lmPS_dist
!       integer :: lmPA_dist
      integer :: lm_glb 
!       integer :: lmS_glb
!       integer :: lmA_glb
!       integer :: lmP_glb
!       integer :: lmPS_glb
!       integer :: lmPA_glb
      
      dwdt_dist    = zero
      dzdt_dist    = zero
      AdvPol_dist  = zero
      AdvTor_dist  = zero
      CorPol_dist  = zero
      dVxVhLM_dist = zero
      Buo_dist     = zero
      LFPol_dist   = zero
      LFTor_dist   = zero
      
      lm_maybe_skip_first = 1
      if (dist_map%lm2(0,0) > 0) lm_maybe_skip_first = 2
    
      if (nBc == 0 .or. lRmsCalc ) then
    
         if ( l_conv ) then  ! Convection
            
            if (dist_map%lm2(0,0) > 0) then  ! if m=0 is in this rank
               lm  = dist_map%lm2(0,0)
               lmA = dist_map%lm2lmA(lm)
               lmP = dist_map%lm2lmP(lm)
               lmPA= dist_map%lmP2lmPA(lmP)
               lm_glb = lm2(0,0) ! needed for leg_helper only!
               
               if ( l_conv_nl ) then
                  AdvPol_loc=      or2(nR)*this_dist%AdvrLM(lm)
                  AdvTor_loc=-dTheta1A_loc(lm)*this_dist%AdvpLM(lmPA)
               else
                  AdvPol_loc=zero
                  AdvTor_loc=zero
               end if
               if ( l_corr ) then
                  CorPol_loc=two*CorFac*or1(nR) * dTheta2A_loc(lm)* z_dist(lmA,nR)
                  CorTor_loc= two*CorFac*or2(nR) * (                 &
                  &                dTheta3A_loc(lm)*dw_dist(lmA,nR) +    &
                  &        or1(nR)*dTheta4A_loc(lm)* w_dist(lmA,nR) )
               else
                  CorPol_loc=zero
                  CorTor_loc=zero
               end if

               if ( l_single_matrix ) then
                  dwdt_dist(lm)=AdvPol_loc!+CorPol_loc
               else
                  dwdt_dist(lm)=AdvPol_loc+CorPol_loc
               end if

               dzdt_dist(lm)=AdvTor_loc+CorTor_loc

               if ( lRmsCalc ) then

                  Buo_dist(lm) =BuoFac*rgrav(nR)*rho0(nR)*leg_helper%sR(lm_glb)
                  if ( l_mag_LF .and. nR>n_r_LCR ) then
                     LFPol_dist(lm) =      or2(nR)*this_dist%LFrLM(lm)
                     LFTor_dist(lm) =-dTheta1A_loc(lm)*this_dist%LFpLM(lmPA)
                     AdvPol_dist(lm)=AdvPol_loc-LFPol_dist(lm)
                     AdvTor_dist(lm)=AdvTor_loc-LFTor_dist(lm)
                  else
                     AdvPol_dist(lm)=AdvPol_loc
                     AdvTor_dist(lm)=AdvTor_loc
                  end if
                  CorPol_dist(lm)=CorPol_loc
               end if
            end if
            
            do lm=lm_maybe_skip_first,lm_loc
               l   =dist_map%lm2l(lm)
               m   =dist_map%lm2m(lm)
               lmS =dist_map%lm2lmS(lm)
               lmA =dist_map%lm2lmA(lm)
               lmP =dist_map%lm2lmP(lm)
               lmPS=dist_map%lmP2lmPS(lmP)
               lmPA=dist_map%lmP2lmPA(lmP)
               lm_glb = lm2(l,m) ! needed for leg_helper only!
               
               if ( l_double_curl ) then ! Pressure is not needed

                  if ( l_corr ) then
                     if ( l < l_max .and. l > m ) then
                        CorPol_loc =two*CorFac*or2(nR)*orho1(nR)*(               &
                        &                   dPhi0_loc(lm)*(                          &
                        &         -ddw_dist(lm,nR)+beta(nR)*dw_dist(lm,nR)     + &
                        &             ( beta(nR)*or1(nR)+or2(nR))*               &
                        &                         leg_helper%dLhw(lm_glb) )        + &
                        &             dTheta3A_loc(lm)*( dz_dist(lmA,nR)-            &
                        &                            beta(nR)*z_dist(lmA,nR) ) + &
                        &             dTheta3S_loc(lm)*( dz_dist(lmS,nR)-            &
                        &                            beta(nR)*z_dist(lmS,nR) ) + &
                        &          or1(nR)* (                                    &
                        &             dTheta4A_loc(lm)* z_dist(lmA,nR)               &
                        &            -dTheta4S_loc(lm)* z_dist(lmS,nR) ) )
                     else if ( l == l_max ) then
                        CorPol_loc =two*CorFac*or2(nR)*orho1(nR)*(               &
                        &                   dPhi0_loc(lm)*(                          &
                        &         -ddw_dist(lm,nR)+beta(nR)*dw_dist(lm,nR)     + &
                        &             ( beta(nR)*or1(nR)+or2(nR))*               &
                        &                         leg_helper%dLhw(lm_glb) ) )
                     else if ( l == m ) then
                        CorPol_loc =two*CorFac*or2(nR)*orho1(nR)*(               &
                        &                   dPhi0_loc(lm)*(                          &
                        &         -ddw_dist(lm,nR)+beta(nR)*dw_dist(lm,nR)     + &
                        &             ( beta(nR)*or1(nR)+or2(nR))*               &
                        &                         leg_helper%dLhw(lm_glb) )        + &
                        &             dTheta3A_loc(lm)*( dz_dist(lmA,nR)-            &
                        &                            beta(nR)*z_dist(lmA,nR) ) + &
                        &          or1(nR)* (                                    &
                        &             dTheta4A_loc(lm)* z_dist(lmA,nR) ) )
                     end if
                  else
                     CorPol_loc=zero
                  end if

                  if ( l_conv_nl ) then

                     if ( l > m ) then
                        dVxVhLM_dist(lm)=      orho1(nR)*r(nR)*r(nR)* ( &
                        &        dTheta1S_loc(lm)*this_dist%AdvtLM(lmPS) -  &
                        &        dTheta1A_loc(lm)*this_dist%AdvtLM(lmPA) +  &
                        &             dPhi_loc(lm)*this_dist%AdvpLM(lmP)  )
                     else if ( l == m ) then
                        dVxVhLM_dist(lm)=      orho1(nR)*r(nR)*r(nR)* ( &
                        &      - dTheta1A_loc(lm)*this_dist%AdvtLM(lmPA) +  &
                        &        dPhi_loc(lm)*this_dist%AdvpLM(lmP)  )
                     end if

                     AdvPol_loc=dLh_loc(lm)*or4(nR)*orho1(nR)*this_dist%AdvrLM(lmP)

                  else

                     AdvPol_loc =zero
                     dVxVhLM_dist(lm)=zero

                  endif

               else ! We don't use the double curl

                  if ( l_corr .and. nBc /= 2 ) then
                     if ( l < l_max .and. l > m ) then
                        CorPol_loc =two*CorFac*or1(nR) * (  &
                        &       dPhi0_loc(lm)*dw_dist(lm,nR) +  & ! phi-deriv of dw/dr
                        &    dTheta2A_loc(lm)*z_dist(lmA,nR) -  & ! sin(theta) dtheta z
                        &    dTheta2S_loc(lm)*z_dist(lmS,nR) )
                     else if ( l == l_max ) then
                        CorPol_loc= two*CorFac*or1(nR) * ( &
                        &            dPhi0_loc(lm)*dw_dist(lm,nR)  )
                     else if ( l == m ) then
                        CorPol_loc = two*CorFac*or1(nR) * (  &
                        &        dPhi0_loc(lm)*dw_dist(lm,nR)  + &
                        &     dTheta2A_loc(lm)*z_dist(lmA,nR) )
                     end if
                  else
                     CorPol_loc=zero
                  end if

                  if ( l_conv_nl ) then
                     AdvPol_loc=or2(nR)*this_dist%AdvrLM(lmP)
                  else
                     AdvPol_loc=zero
                  endif

               end if ! Double curl or not for the poloidal equation

               dwdt_dist(lm)=AdvPol_loc+CorPol_loc

               if ( lRmsCalc ) then ! RMS force balance

                  if ( l_TP_form .or. l_anelastic_liquid ) then
                     Buo_dist(lm) =BuoFac*alpha0(nR)*rgrav(nR)*(              &
                     &        rho0(nR)*leg_helper%sR(lm_glb)-ViscHeatFac*    &
                     &        (ThExpNb*alpha0(nR)*temp0(nR)+ogrun(nR))*  &
                     &        leg_helper%preR(lm_glb) )
                  else
                     Buo_dist(lm) =BuoFac*rho0(nR)*rgrav(nR)*leg_helper%sR(lm_glb)
                  end if

                  if ( l_double_curl ) then 
                     ! In that case we have to recompute the Coriolis force
                     ! since we also want the pressure gradient
                     if ( l_corr .and. nBc /= 2 ) then
                        if ( l < l_max .and. l > m ) then
                           CorPol_loc =two*CorFac*or1(nR) * (  &
                           &       dPhi0_loc(lm)*dw_dist(lm,nR) +  & ! phi-deriv of dw/dr
                           &    dTheta2A_loc(lm)*z_dist(lmA,nR) -  & ! sin(theta) dtheta z
                           &    dTheta2S_loc(lm)*z_dist(lmS,nR) )
                        else if ( l == l_max ) then
                           CorPol_loc= two*CorFac*or1(nR) * ( &
                           &            dPhi0_loc(lm)*dw_dist(lm,nR)  )
                        else if ( l == m ) then
                           CorPol_loc = two*CorFac*or1(nR) * (  &
                           &        dPhi0_loc(lm)*dw_dist(lm,nR)  + &
                           &     dTheta2A_loc(lm)*z_dist(lmA,nR) )
                        end if
                     else
                        CorPol_loc=zero
                     end if

                     ! We also need to recompute AdvPol_loc here
                     if ( l_conv_nl ) then
                        AdvPol_loc=or2(nR)*this_dist%AdvrLM(lmP)
                     else
                        AdvPol_loc=zero
                     endif

                  end if

                  if ( l_mag_LF .and. nR>n_r_LCR ) then
                     LFPol_dist(lm) =or2(nR)*this_dist%LFrLM(lmP)
                     AdvPol_dist(lm)=AdvPol_loc-LFPol_dist(lm)
                  else
                     AdvPol_dist(lm)=AdvPol_loc
                  end if
                  CorPol_dist(lm)=CorPol_loc

               end if

               if ( l_corr ) then
                  if ( l < l_max .and. l > m ) then
                     CorTor_loc=          two*CorFac*or2(nR) * (  &
                     &                dPhi0_loc(lm)*z_dist(lm,nR)   + &
                     &            dTheta3A_loc(lm)*dw_dist(lmA,nR)  + &
                     &    or1(nR)*dTheta4A_loc(lm)* w_dist(lmA,nR)  + &
                     &            dTheta3S_loc(lm)*dw_dist(lmS,nR)  - &
                     &    or1(nR)*dTheta4S_loc(lm)* w_dist(lmS,nR)  )
                  else if ( l == l_max ) then
                     CorTor_loc=two*CorFac*or2(nR) * ( &
                     &            dPhi0_loc(lm)*z_dist(lm,nR)   )
                  else if ( l == m ) then
                     CorTor_loc=  two*CorFac*or2(nR) * (  &
                     &        dPhi0_loc(lm)*z_dist(lm,nR)   + &
                     &    dTheta3A_loc(lm)*dw_dist(lmA,nR)  + &
                     &    or1(nR)*dTheta4A_loc(lm)* w_dist(lmA,nR)  )
                  end if
               else
                  CorTor_loc=zero
               end if
    
               if ( l_conv_nl ) then
                  if ( l > m ) then
                     AdvTor_loc=   -dPhi_loc(lm)*this_dist%AdvtLM(lmP)  + &
                     &          dTheta1S_loc(lm)*this_dist%AdvpLM(lmPS) - &
                     &          dTheta1A_loc(lm)*this_dist%AdvpLM(lmPA)
                  else if ( l == m ) then
                     AdvTor_loc=   -dPhi_loc(lm)*this_dist%AdvtLM(lmP)  - &
                     &          dTheta1A_loc(lm)*this_dist%AdvpLM(lmPA)
                  end if
               else
                  AdvTor_loc=zero
               end if
    
               dzdt_dist(lm)=CorTor_loc+AdvTor_loc
               ! until here
    
               if ( lRmsCalc ) then
                  if ( l_mag_LF .and. nR>n_r_LCR ) then
                     !------ When RMS values are required, the Lorentz force is treated
                     !       separately:
       
                     if ( l > m ) then
                        !------- LFTor= 1/(E*Pm) * curl( curl(B) x B )_r
                        LFTor_dist(lm) =   -dPhi_loc(lm)*this_dist%LFtLM(lmP)  + &
                        &          dTheta1S_loc(lm)*this_dist%LFpLM(lmPS) - &
                        &          dTheta1A_loc(lm)*this_dist%LFpLM(lmPA)
                     else if ( l == m ) then
                        LFTor_dist(lm) =   -dPhi_loc(lm)*this_dist%LFtLM(lmP)  - &
                        &          dTheta1A_loc(lm)*this_dist%LFpLM(lmPA)
                     end if
                     AdvTor_dist(lm)=AdvTor_loc-LFTor_dist(lm)
                  else
                     AdvTor_dist(lm)=AdvTor_loc
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
                  call hIntRms(AdvPol,nR,1,lm_max,0,Adv2hInt(:,nR),st_map, .false.)
                  call hIntRms(this%Advt2LM,nR,1,lmP_max,1,Adv2hInt(:,nR),st_map, &
                       &       .true.)
                  call hIntRms(this%Advp2LM,nR,1,lmP_max,1,Adv2hInt(:,nR),st_map, &
                       &       .true.)
               end if

               if ( l_TP_form .or. l_anelastic_liquid ) then
                  call hIntRms(leg_helper%dpR,nR,1,lm_max,0, &
                       &       Pre2hInt(:,nR),st_map,.false.)
               else
                  call hIntRms(leg_helper%dpR-beta(nR)*leg_helper%preR,&
                       &       nR,1,lm_max,0,Pre2hInt(:,nR),st_map,.false.)
               end if
               call hIntRms(this%PFt2LM,nR,1,lmP_max,1,Pre2hInt(:,nR),st_map,.true.)
               call hIntRms(this%PFp2LM,nR,1,lmP_max,1,Pre2hInt(:,nR),st_map,.true.)

               ! rho* grad(p/rho) = grad(p) - beta*p
               if ( ra /= 0.0_cp ) &
                  call hIntRms(Buo,nR,1,lm_max,0,Buo2hInt(:,nR),st_map,.false.)
               if ( l_corr ) then
                  call hIntRms(CorPol,nR,1,lm_max,0,Cor2hInt(:,nR),st_map,.false.)
                  call hIntRms(this%CFt2LM,nR,1,lmP_max,1,Cor2hInt(:,nR), &
                       &       st_map,.true.)
                  calL hIntRms(this%CFp2LM,nR,1,lmP_max,1,Cor2hInt(:,nR), &
                       &       st_map,.true.)
               end if
               if ( l_mag_LF .and. nR>n_r_LCR ) then
                  call hIntRms(LFPol,nR,1,lm_max,0,LF2hInt(:,nR),st_map,.false.)
                  call hIntRms(this%LFt2LM,nR,1,lmP_max,1,LF2hInt(:,nR),st_map,.true.)
                  call hIntRms(this%LFp2LM,nR,1,lmP_max,1,LF2hInt(:,nR),st_map,.true.)
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
               call hIntRms(Geo,nR,1,lm_max,0,Geo2hInt(:,nR),st_map,.false.)
               call hIntRms(CLF,nR,1,lm_max,0,CLF2hInt(:,nR),st_map,.false.)
               call hIntRms(PLF,nR,1,lm_max,0,PLF2hInt(:,nR),st_map,.false.)
               call hIntRms(Mag,nR,1,lm_max,0,Mag2hInt(:,nR),st_map,.false.)
               call hIntRms(Arc,nR,1,lm_max,0,Arc2hInt(:,nR),st_map,.false.)
               call hIntRms(ArcMag,nR,1,lm_max,0,ArcMag2hInt(:,nR),st_map,.false.)
               call hIntRms(CIA,nR,1,lm_max,0,CIA2hInt(:,nR),st_map,.false.)
      !---------------------------------------------------------------------------------------------
      ! This is (more or less) how this specific piece will look like afterwards: 
      !---------------------------------------------------------------------------------------------
      !                do lm=1,lm_loc
      !                   l = dist_map%lm2l(lm)
      !                   m = dist_map%lm2m(lm)
      !                   lm_glb = lm2(l,m)
      !                   
      !                   Geo_dist(lm)=CorPol_dist(lm)-leg_helper%dpR(lm_glb)+beta(nR)*leg_helper%preR(lm_glb)
      !                   CLF_dist(lm)=CorPol_dist(lm)+LFPol_dist(lm)
      !                   PLF_dist(lm)=LFPol_dist(lm)-leg_helper%dpR(lm_glb)+beta(nR)*leg_helper%preR(lm_glb)
      !                   Mag_dist(lm)=Geo_dist(lm)+LFPol_dist(lm)
      !                   Arc_dist(lm)=Geo_dist(lm)+Buo_dist(lm)
      !                   ArcMag_dist(lm)=Mag_dist(lm)+Buo_dist(lm)
      !                   CIA_dist(lm)=ArcMag_dist(lm)+AdvPol_dist(lm)
      !                   !CIA_dist(lm)=CorPol(lm_glb)+Buo(lm_glb)+AdvPol(lm_glb)
      !                end do
      !
      !                call hIntRms(Geo_dist(1:lm_loc),nR,1,lm_loc,0,Geo2hInt(0:l_max,nR),dist_map,.false.)
      !                call hIntRms(CLF_dist(1:lm_loc),nR,1,lm_loc,0,CLF2hInt(0:l_max,nR),dist_map,.false.)
      !                call hIntRms(PLF_dist(1:lm_loc),nR,1,lm_loc,0,PLF2hInt(0:l_max,nR),dist_map,.false.)
      !                call hIntRms(Mag_dist(1:lm_loc),nR,1,lm_loc,0,Mag2hInt(0:l_max,nR),dist_map,.false.)
      !                call hIntRms(Arc_dist(1:lm_loc),nR,1,lm_loc,0,Arc2hInt(0:l_max,nR),dist_map,.false.)
      !                call hIntRms(ArcMag_dist(1:lm_loc),nR,1,lm_loc,0,ArcMag2hInt(0:l_max,nR),dist_map,.false.)
      !                call hIntRms(CIA_dist(1:lm_loc),nR,1,lm_loc,0,CIA2hInt(0:l_max,nR),dist_map,.false.)               
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
      !                call gather_Flm(Geo_dist   , Geo)
      !                call gather_Flm(CLF_dist   , CLF)
      !                call gather_Flm(PLF_dist   , PLF)
      !                call gather_Flm(Mag_dist   , Mag)
      !                call gather_Flm(Arc_dist   , Arc)
      !                call gather_Flm(ArcMag_dist, ArcMag)
      !                call gather_Flm(CIA_dist   , CIA)
      !---------------------------------------------------------------------------------------------
      ! Not very charming, but the gathers will drop (those are local variables)
      ! I am not sure yet, but maybe I can afford to do only a single mpi_iallreduce after ALL those 
      ! lm loops which appear next. I think it is just a sum, but I need to double check it.
      !---------------------------------------------------------------------------------------------

               do lm=1,lm_max
                  lmP =lm2lmP(lm)
                  Geo(lm)=-this%CFt2LM(lmP)-this%PFt2LM(lmP)
                  CLF(lm)=-this%CFt2LM(lmP)+this%LFt2LM(lmP)
                  PLF(lm)=this%LFt2LM(lmP)-this%PFt2LM(lmP)
                  Mag(lm)=Geo(lm)+this%LFt2LM(lmP)
                  Arc(lm)=Geo(lm)
                  ArcMag(lm)=Mag(lm)
                  CIA(lm)=ArcMag(lm)+this%Advt2LM(lmP)
                  !CIA(lm)=-this%CFt2LM(lmP)+this%Advt2LM(lmP)
               end do
               call hIntRms(Geo,nR,1,lm_max,0,Geo2hInt(:,nR),st_map,.true.)
               call hIntRms(CLF,nR,1,lm_max,0,CLF2hInt(:,nR),st_map,.true.)
               call hIntRms(PLF,nR,1,lm_max,0,PLF2hInt(:,nR),st_map,.true.)
               call hIntRms(Mag,nR,1,lm_max,0,Mag2hInt(:,nR),st_map,.true.)
               call hIntRms(Arc,nR,1,lm_max,0,Arc2hInt(:,nR),st_map,.true.)
               call hIntRms(ArcMag,nR,1,lm_max,0,ArcMag2hInt(:,nR),st_map,.true.)
               call hIntRms(CIA,nR,1,lm_max,0,CIA2hInt(:,nR),st_map,.true.)
    
               do lm=1,lm_max
                  lmP =lm2lmP(lm)
                  Geo(lm)=-this%CFp2LM(lmP)-this%PFp2LM(lmP)
                  CLF(lm)=-this%CFp2LM(lmP)+this%LFp2LM(lmP)
                  PLF(lm)=this%LFp2LM(lmP)-this%PFp2LM(lmP)
                  Mag(lm)=Geo(lm)+this%LFp2LM(lmP)
                  Arc(lm)=Geo(lm)
                  ArcMag(lm)=Mag(lm)
                  CIA(lm)=ArcMag(lm)+this%Advp2LM(lmP)
                  !CIA(lm)=-this%CFp2LM(lmP)+this%Advp2LM(lmP)
               end do
               call hIntRms(Geo,nR,1,lm_max,0,Geo2hInt(:,nR),st_map,.true.)
               call hIntRms(CLF,nR,1,lm_max,0,CLF2hInt(:,nR),st_map,.true.)
               call hIntRms(PLF,nR,1,lm_max,0,PLF2hInt(:,nR),st_map,.true.)
               call hIntRms(Mag,nR,1,lm_max,0,Mag2hInt(:,nR),st_map,.true.)
               call hIntRms(Arc,nR,1,lm_max,0,Arc2hInt(:,nR),st_map,.true.)
               call hIntRms(ArcMag,nR,1,lm_max,0,ArcMag2hInt(:,nR),st_map,.true.)
               call hIntRms(CIA,nR,1,lm_max,0,CIA2hInt(:,nR),st_map,.true.)

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
               do lm=lm_maybe_skip_first,lm_loc
                  l   =dist_map%lm2l(lm)
                  m   =dist_map%lm2m(lm)
                  lmS =dist_map%lm2lmS(lm)
                  lmA =dist_map%lm2lmA(lm)
                  lmP =dist_map%lm2lmP(lm)
                  lmPS=dist_map%lmP2lmPS(lmP)
                  lmPA=dist_map%lmP2lmPA(lmP)
                  lm_glb = lm2(l,m) ! needed for leg_helper only!
       
                  !------ Recycle CorPol and AdvPol:
                  if ( l_corr ) then
                     !PERFON('td_cv2c')
                     if ( l < l_max .and. l > m ) then
                        CorPol_loc=                    two*CorFac*or2(nR) *  &
                        &                    ( -dPhi0_loc(lm) * ( dw_dist(lm,nR) &
                        &                       +or1(nR)*leg_helper%dLhw(lm_glb) &
                        &                                                  ) &
                        &                       +dTheta3A_loc(lm)*z_dist(lmA,nR) &
                        &                       +dTheta3S_loc(lm)*z_dist(lmS,nR) &
                        &                    )
       
                     else if ( l == l_max ) then
                        CorPol_loc=  two*CorFac*or2(nR) * ( -dPhi0_loc(lm) *  &
                                    ( dw_dist(lm,nR) + or1(nR)*leg_helper%dLhw(lm_glb) ) )
       
                     else if ( l == m ) then
                        CorPol_loc=                    two*CorFac*or2(nR) *  &
                        &                    ( -dPhi0_loc(lm) * ( dw_dist(lm,nR) &
                        &                       +or1(nR)*leg_helper%dLhw(lm_glb) &
                        &                                                   )&
                        &                      +dTheta3A_loc(lm)*z_dist(lmA,nR)  &
                        &                    )
       
                     end if
                     !PERFOFF
                  else
                     CorPol_loc=zero
                  end if
                  if ( l_conv_nl ) then
                     !PERFON('td_cv2nl')
                     if ( l > m ) then
                        AdvPol_loc= dTheta1S_loc(lm)*this_dist%AdvtLM(lmPS) - &
                        &           dTheta1A_loc(lm)*this_dist%AdvtLM(lmPA) + &
                        &               dPhi_loc(lm)*this_dist%AdvpLM(lmP)
                     else if ( l == m ) then
                        AdvPol_loc=-dTheta1A_loc(lm)*this_dist%AdvtLM(lmPA) + &
                        &               dPhi_loc(lm)*this_dist%AdvpLM(lmP)
                     end if
                     !PERFOFF
                  else
                     AdvPol_loc=zero
                  end if
                  dpdt_dist(lm)=AdvPol_loc+CorPol_loc
       
               end do ! lm loop
            end if
         else
            do lm=lm_maybe_skip_first,lm_loc
               dwdt_dist(lm) =0.0_cp
               dzdt_dist(lm) =0.0_cp
               dpdt_dist(lm) =0.0_cp
            end do
         end if ! l_conv ?

      end if
!---------------------------------------------------------------------------------------------------
!                    GATHER
!---------------------------------------------------------------------------------------------------
      !@>TODO: lots of things are set to zero (e.g. if some l_* flags are deactive)
      !        those are, of course, unnecessary to gather. Others will bug everything
      !        out if gathered e.g. dVxVhLM.
      call gather_Flm(dwdt_dist   (1:lm_loc),dwdt   (1:lm_max)) 
      call gather_Flm(dzdt_dist   (1:lm_loc),dzdt   (1:lm_max)) 
      call gather_Flm(dpdt_dist   (1:lm_loc),dpdt   (1:lm_max)) 
      call gather_Flm(AdvPol_dist (1:lm_loc),AdvPol (1:lm_max)) 
      call gather_Flm(AdvTor_dist (1:lm_loc),AdvTor (1:lm_max)) 
      call gather_Flm(CorPol_dist (1:lm_loc),CorPol (1:lm_max)) 
      if (l_double_curl) call gather_Flm(dVxVhLM_dist(1:lm_loc),dVxVhLM(1:lm_max))
      call gather_Flm(Buo_dist    (1:lm_loc),Buo    (1:lm_max)) 
      call gather_Flm(LFPol_dist  (1:lm_loc),LFPol  (1:lm_max)) 
      call gather_Flm(LFTor_dist  (1:lm_loc),LFTor  (1:lm_max)) 
      
!---------------------------------------------------------------------------------------------------
!                    GATHER
!---------------------------------------------------------------------------------------------------    
      if ( nBc == 0 ) then
         
         if ( l_heat ) then
            
            if (dist_map%lm2(0,0) > 0) then  ! if m=0 is in this rank
               lm  = dist_map%lm2(0,0)

               dsdt_loc  =epsc*epscProf(nR) !+opr/epsS*divKtemp0(nR)
               dVSrLM_dist(lm)=this_dist%VSrLM(lm)
               if ( l_TP_form ) dVPrLM_dist(lm)=this_dist%VPrLM(lm)
               if ( l_anel ) then
                  if ( l_anelastic_liquid .or. l_TP_form ) then
                     if ( l_mag_nl ) then
                        dsdt_loc=dsdt_loc+                                        &
                        &    ViscHeatFac*hdif_V_loc(lm)*temp0(nR)*this_dist%ViscHeatLM(lm)+  &
                        &     OhmLossFac*hdif_B_loc(lm)*temp0(nR)*this_dist%OhmLossLM(lm)
                     else
                        dsdt_loc=dsdt_loc+ &
                        &    ViscHeatFac*hdif_V_loc(lm)*temp0(nR)*this_dist%ViscHeatLM(lm)
                     end if
                  else
                     if ( l_mag_nl ) then
                        dsdt_loc=dsdt_loc+ViscHeatFac*hdif_V_loc(lm)*this_dist%ViscHeatLM(lm)+ &
                        &                  OhmLossFac*hdif_B_loc(lm)*this_dist%OhmLossLM(lm)
                     else
                        dsdt_loc=dsdt_loc+ViscHeatFac*hdif_V_loc(lm)*this_dist%ViscHeatLM(lm)
                     end if
                  end if
               end if
               dsdt_dist(lm)=dsdt_loc
            end if
    
            do lm=lm_maybe_skip_first,lm_loc
               l   =dist_map%lm2l(lm)
               m   =dist_map%lm2m(lm)
               lmP =dist_map%lm2lmP(lm)
               lmPS=dist_map%lmP2lmPS(lmP)
               lmPA=dist_map%lmP2lmPA(lmP)
               
               !------ This is horizontal heat advection:
               !PERFON('td_h1')
    
               if ( l > m ) then
                  dsdt_loc= -dTheta1S_loc(lm)*this_dist%VStLM(lmPS) &
                  &         +dTheta1A_loc(lm)*this_dist%VStLM(lmPA) &
                  &         -dPhi_loc(lm)*this_dist%VSpLM(lmP)
               else if ( l == m ) then
                  dsdt_loc=  dTheta1A_loc(lm)*this_dist%VStLM(lmPA) &
                  &          -dPhi_loc(lm)*this_dist%VSpLM(lmP)
               end if
               !PERFOFF
               !PERFON('td_h2')
               if ( l_anel ) then
                  if ( l_anelastic_liquid .or. l_TP_form ) then
                     dsdt_loc = dsdt_loc+ &
                     &          ViscHeatFac*hdif_V_loc(lm)*temp0(nR)*this_dist%ViscHeatLM(lmP)
                     if ( l_mag_nl ) then
                        dsdt_loc = dsdt_loc+ &
                        &          OhmLossFac*hdif_B_loc(lm)*temp0(nR)*this_dist%OhmLossLM(lmP)
                     end if
                  else
                     dsdt_loc = dsdt_loc+ &
                     &          ViscHeatFac*hdif_V_loc(lm)*this_dist%ViscHeatLM(lmP)
                     if ( l_mag_nl ) then
                        dsdt_loc = dsdt_loc+ &
                        &          OhmLossFac*hdif_B_loc(lm)*this_dist%OhmLossLM(lmP)
                     end if
                  end if
               end if
               !PERFOFF
               !-----   simplified form for linear onset !
               !        not ds not saved in the current program form!
               !                 dsdt_dist(lm)=
               !                    -dLh(lm)*w(lm,nR)*or2(nR)*dsR(1)
               dVSrLM_dist(lm)=this_dist%VSrLM(lmP)
               dsdt_dist(lm) = dsdt_loc
               if ( l_TP_form ) dVPrLM_dist(lm)=this_dist%VPrLM(lmP)
            end do
         else
            
            do lm=lm_maybe_skip_first,lm_loc
               dsdt_dist(lm)  =0.0_cp
               dVSrLM_dist(lm)=0.0_cp
            end do
         end if

         if ( l_chemical_conv ) then
            if (dist_map%lm2(0,0) > 0) then  ! if m=0 is in this rank
               lm  = dist_map%lm2(0,0)
               dVXirLM_dist(lm)=this_dist%VXirLM(lm)
               dxidt_dist(lm)  =epscXi
            end if
    
            do lm=lm_maybe_skip_first,lm_loc
               l   =dist_map%lm2l(lm)
               m   =dist_map%lm2m(lm)
               lmP =dist_map%lm2lmP(lm)
               lmPS=dist_map%lmP2lmPS(lmP)
               lmPA=dist_map%lmP2lmPA(lmP)
               !------ This is horizontal heat advection:
    
               if ( l > m ) then
                  dxidt_loc= -dTheta1S_loc(lm)*this_dist%VXitLM(lmPS) &
                  &          +dTheta1A_loc(lm)*this_dist%VXitLM(lmPA) &
                  &          -dPhi_loc(lm)*this_dist%VXipLM(lmP)
               else if ( l == m ) then
                  dxidt_loc=  dTheta1A_loc(lm)*this_dist%VXitLM(lmPA) &
                  &          -dPhi_loc(lm)*this_dist%VXipLM(lmP)
               end if
               dVXirLM_dist(lm)=this_dist%VXirLM(lmP)
               dxidt_dist(lm)  =dxidt_loc
            end do
         end if
    
         if ( l_mag_nl .or. l_mag_kin  ) then
            do lm=1,lm_loc
               l   =dist_map%lm2l(lm)
               m   =dist_map%lm2m(lm)
               lmP =dist_map%lm2lmP(lm)
               lmPS=dist_map%lmP2lmPS(lmP)
               lmPA=dist_map%lmP2lmPA(lmP)
               
               if ((l == 0) .and. (m == 0)) then
                  dVxBhLM_dist(lm)= -r(nR)*r(nR)* dTheta1A_loc(lm)*this_dist%VxBtLM(lmPA)
                  dbdt_dist(lm)   = -dTheta1A_loc(lm)*this_dist%VxBpLM(lmPA)
                  djdt_dist(lm)   = zero
                  cycle ! <---------------------------

               !------- This is the radial part of the dynamo terms \curl(VxB)
               else if ( l > m ) then

                  dbdt_dist(lm)=  dTheta1S_loc(lm)*this_dist%VxBpLM(lmPS) &
                  &         -dTheta1A_loc(lm)*this_dist%VxBpLM(lmPA) &
                  &         -dPhi_loc(lm)    *this_dist%VxBtLM(lmP)

               else if ( l == m ) then
                  dbdt_dist(lm)= -dTheta1A_loc(lm)*this_dist%VxBpLM(lmPA) &
                  &         -dPhi_loc(lm)    *this_dist%VxBtLM(lmP)
               end if
    
               !------- Radial component of
               !           \curl\curl(UxB) = \grad\div(UxB) - \laplace(VxB)
    
               !------- This is the radial part of \laplace (UxB)
               djdt_dist(lm)=dLh_loc(lm)*or4(nR)*this_dist%VxBrLM(lmP)
    
               !------- This is r^2 * horizontal divergence of (UxB)
               !        Radial derivative performed in get_dr_td
               if ( l > m ) then
                  dVxBhLM_dist(lm)=            r(nR)*r(nR)* ( &
                  &    dTheta1S_loc(lm)*this_dist%VxBtLM(lmPS) -  &
                  &    dTheta1A_loc(lm)*this_dist%VxBtLM(lmPA) +  &
                  &    dPhi_loc(lm)*this_dist%VxBpLM(lmP)  )
               else if ( l == m ) then
                  dVxBhLM_dist(lm)=              r(nR)*r(nR)* ( &
                  &    - dTheta1A_loc(lm)*this_dist%VxBtLM(lmPA) +  &
                  &    dPhi_loc(lm)*this_dist%VxBpLM(lmP)  )
               end if
            end do
         else if ( l_mag ) then
            do lm=1,lm_loc
               dbdt_dist(lm)   =zero
               djdt_dist(lm)   =zero
               dVxBhLM_dist(lm)=zero
            end do
         end if
!---------------------------------------------------------------------------------------------------
!                    GATHER
!---------------------------------------------------------------------------------------------------
         call gather_Flm(dsdt_dist   (1:lm_loc),dsdt   (1:lm_max)) 
         call gather_Flm(dVSrLM_dist (1:lm_loc),dVSrLM (1:lm_max)) 
         if ( l_heat .and. l_TP_form ) call gather_Flm(dVPrLM_dist (1:lm_loc),dVPrLM (1:lm_max)) 
         if ( l_chemical_conv ) call gather_Flm(dVXirLM_dist(1:lm_loc), dVXirLM (1:lm_max)) 
         if ( l_chemical_conv ) call gather_Flm(dxidt_dist  (1:lm_loc), dxidt   (1:lm_max)) 
         if ( l_mag_nl .or. l_mag_kin .or. l_mag ) call gather_Flm(dbdt_dist   (1:lm_loc), dbdt   (1:lm_max)) 
         if ( l_mag_nl .or. l_mag_kin .or. l_mag ) call gather_Flm(djdt_dist   (1:lm_loc), djdt   (1:lm_max)) 
         if ( l_mag_nl .or. l_mag_kin .or. l_mag ) call gather_Flm(dVxBhLM_dist(1:lm_loc), dVxBhLM(1:lm_max)) 
!---------------------------------------------------------------------------------------------------
!                    GATHER
!---------------------------------------------------------------------------------------------------   

!---------------------------------------------------------------------------------------------------
!                    DONE
!---------------------------------------------------------------------------------------------------
      else   ! boundary !

         if ( l_mag_nl .or. l_mag_kin ) then
    
            !----- Stress free boundary, only nl mag. term for poloidal field needed.
            !      Because the radial derivative will be taken, this will contribute to
            !      the other radial grid points.
            do lm=1,lm_loc
               l   = dist_map%lm2l(lm)
               m   = dist_map%lm2m(lm)
               lmP = dist_map%lm2lmP(lm)
               lmPS= dist_map%lmP2lmPS(lmP)   ! l-1
               lmPA= dist_map%lmP2lmPA(lmP)   ! l+1
               
               if ((l == 0) .and. (m == 0)) then
                  dVxBhLM_dist(lm)=zero
                  dVSrLM_dist(lm) =zero
                  cycle ! <---------------------------
               else if ( l > m ) then
                  dVxBhLM_dist(lm)=r(nR)*r(nR)* (               &
                  &      dTheta1S_loc(lm)*this_dist%VxBtLM(lmPS) -  &
                  &      dTheta1A_loc(lm)*this_dist%VxBtLM(lmPA) +  &
                  &          dPhi_loc(lm)*this_dist%VxBpLM(lmP)  )
               else if ( l == m ) then ! (l-1) not allowed !
                  dVxBhLM_dist(lm)=r(nR)*r(nR)* (               &
                  &    - dTheta1A_loc(lm)*this_dist%VxBtLM(lmPA) +  &
                  &    dPhi_loc(lm)*this_dist%VxBpLM(lmP)  )
               end if
               dVSrLM_dist(lm)=zero
            end do
    
         else
            do lm=1,lm_loc
               if ( l_mag ) dVxBhLM_dist(lm)=zero
               dVSrLM_dist(lm) =zero
            end do
         end if
         if ( l_double_curl ) then
            do lm=1,lm_loc
               dVxVhLM_dist(lm)=zero
            end do
         end if
         if ( l_chemical_conv ) then
            do lm=1,lm_loc
               dVXirLM_dist(lm)=zero
            end do
         end if
         if ( l_TP_form ) then
            do lm=1,lm_loc
               dVPrLM_dist(lm)=zero
            end do
         end if
         
         if ( l_mag_nl .or. l_mag_kin .or. l_mag ) call gather_Flm(dVxBhLM_dist(1:lm_loc), dVxBhLM(1:lm_max)) 
         call gather_Flm(dVSrLM_dist (1:lm_loc), dVSrLM(1:lm_max)) 
         if ( l_double_curl ) call gather_Flm(dVxVhLM_dist (1:lm_loc), dVxVhLM(1:lm_max))
         if ( l_chemical_conv ) call gather_Flm(dVXirLM_dist (1:lm_loc), dVXirLM(1:lm_max)) 
         if ( l_TP_form ) call gather_Flm(dVPrLM_dist (1:lm_loc), dVPrLM(1:lm_max)) 
      end if  ! boundary ? lvelo ?

   end subroutine get_td_dist
!-----------------------------------------------------------------------------
   
end module arrays_dist