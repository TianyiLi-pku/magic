#include "perflib_preproc.cpp"
module LMLoop_mod

#ifdef WITH_LIKWID
#include "likwid_f90.h"
#endif

   use fields
   use fieldsLast
   use omp_lib
   use precision_mod
   use parallel_mod
   use mem_alloc, only: memWrite, bytes_allocated
   use geometry, only: l_max, lm_max, n_r_max, n_r_maxMag, n_r_icb,    &
       &            n_r_cmb, n_mlo_loc
!    use blocking, only: lmStartB, lmStopB, lo_map
   use logic, only: l_mag, l_conv, l_anelastic_liquid, lVerbose, l_heat, &
       &            l_single_matrix, l_chemical_conv, l_TP_form,         &
       &            l_save_out
   use output_data, only: n_log_file, log_file
   use timing, only: wallTime,subTime,writeTime
   use LMLoop_data, only: llm, ulm, llmMag, ulmMag
   use debugging,  only: debug_write
   use communications, only: GET_GLOBAL_SUM, lo2r_redist_start, lo2r_xi, &
       &                    lo2r_s, lo2r_flow, lo2r_field,               &
       &                    lo2r_redist_start_dist, transform_new2old, transform_old2new, &
       &                    lo2r_redist_wait_dist, ml2r_s, get_global_sum_dist, test_field
   use updateS_mod
   use updateZ_mod
   use updateWP_mod
   use updateWPT_mod
   use updateWPS_mod
   use updateB_mod
   use updateXi_mod
   use radial_functions
   
   use radial_der          ! Added by LAGO, remove later
   use lmmapping
   use blocking  !!! uncomment above

   implicit none

   private

   public :: LMLoop, initialize_LMLoop, finalize_LMLoop

contains

   subroutine initialize_LMLoop

      integer(lip) :: local_bytes_used

      local_bytes_used = bytes_allocated

      if ( l_single_matrix ) then
         if ( l_TP_form ) then
            call initialize_updateWPT
         else
            call initialize_updateWPS
         end if
      else
         call initialize_updateS
         call initialize_updateWP
      end if

      if ( l_chemical_conv ) call initialize_updateXi

      call initialize_updateZ
      if ( l_mag ) call initialize_updateB

      local_bytes_used = bytes_allocated-local_bytes_used

      call memWrite('LMLoop.f90',local_bytes_used)

   end subroutine initialize_LMLoop
!----------------------------------------------------------------------------
   subroutine finalize_LMLoop

      if ( l_single_matrix ) then
         if ( l_TP_form ) then
            call finalize_updateWPT
         else
            call finalize_updateWPS
         end if
      else
         call finalize_updateS
         call finalize_updateWP
      end if

      if ( l_chemical_conv ) call finalize_updateXi

      call finalize_updateZ
      ! There is a strange bug here! This finalize_updateB will
      ! make MagIC freeze depending on the distribution of the MPI
      ! ranks! Uncomment as soon as this problem is solved...
!       if ( l_mag ) call finalize_updateB

   end subroutine finalize_LMLoop
!----------------------------------------------------------------------------
   subroutine LMLoop(w1,coex,time,dt,lMat,lRmsNext,lPressNext,dVxVhLM, &
              &      dVxBhLM,dVSrLM,dVPrLM,dVXirLM,dsdt,dwdt,          &
              &      dzdt,dpdt,dxidt,dbdt,djdt,lorentz_torque_ma,      &
              &      lorentz_torque_ic,b_nl_cmb,aj_nl_cmb,             &
              &      aj_nl_icb)
      !
      !  This subroutine performs the actual time-stepping.
      !
      !

      !-- Input of variables:
      real(cp),    intent(in) :: w1,coex
      real(cp),    intent(in) :: dt,time
      logical,     intent(in) :: lMat
      logical,     intent(in) :: lRmsNext
      logical,     intent(in) :: lPressNext

      !--- Input from radialLoop:
      !    These fields are provided in the R-distributed space!
      ! for djdt in update_b
      complex(cp), intent(inout) :: dVxBhLM(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(inout) :: dVxVhLM(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dVSrLM(llm:ulm,n_r_max)   ! for dsdt  in update_s
      complex(cp), intent(inout) :: dVPrLM(llm:ulm,n_r_max)   ! for dsdt  in update_s
      complex(cp), intent(inout) :: dVXirLM(llm:ulm,n_r_max)  ! for dxidt in update_xi
      !integer,     intent(in) :: n_time_step

      !--- Input from radialLoop and then redistributed:
      complex(cp), intent(inout) :: dsdt(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dxidt(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dwdt(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dzdt(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dpdt(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dbdt(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(inout) :: djdt(llmMag:ulmMag,n_r_maxMag)
      real(cp),    intent(in) :: lorentz_torque_ma,lorentz_torque_ic
      complex(cp), intent(in) :: b_nl_cmb(lm_max)   ! nonlinear bc for b at CMB
      complex(cp), intent(in) :: aj_nl_cmb(lm_max)  ! nonlinear bc for aj at CMB
      complex(cp), intent(in) :: aj_nl_icb(lm_max)  ! nonlinear bc for dr aj at ICB

      !--- Local counter
      integer :: nLMB
      integer :: l,nR,ierr
      integer :: tStart(4),tStop(4),tPassed(4)

      logical,parameter :: DEBUG_OUTPUT=.false.
      !--- Inner core rotation from last time step
      real(cp), save :: omega_icLast
      real(cp) :: z10(n_r_max)
      complex(cp) :: sum_dwdt
      
      ! Duplicatas ------------------------------------------
      real(cp)    :: d_omega_ma_dtLast_dist       ! Time derivative of OC rotation of previous step
      real(cp)    :: d_omega_ic_dtLast_dist       ! Time derivative of IC rotation of previous step
      complex(cp) :: dzdt_dist(n_mlo_loc,n_r_max)
      real(cp)    :: omega_ma_dist                ! Calculated OC rotation
      real(cp)    :: omega_ic_dist                ! Calculated IC rotation
      
      complex(cp) :: dVxVhLM_dist(n_mlo_loc,n_r_max)
      complex(cp) :: dwdt_dist(n_mlo_loc,n_r_max) 
      complex(cp) :: dpdt_dist(n_mlo_loc,n_r_max)
      
      !!! New layout  TMP
      !!!----------------------------------
      !!!! SAME !!!
      
      complex(cp) :: dVSrLM_dist(n_mlo_loc,n_r_max)
      complex(cp) :: dsdt_dist(n_mlo_loc,n_r_max)
      
      real(cp)    :: test_norm, error_threshold, test_normi
      integer :: nLMB_start, nLMB_end   
      integer :: m, lm, i, j, k

      call transform_old2new(s_LMloc, s_LMdist)
      call transform_old2new(ds_LMloc, ds_LMdist)
      call transform_old2new(w_LMloc, w_LMdist)
      call transform_old2new(dsdt, dsdt_dist)
      call transform_old2new(dVSrLM, dVSrLM_dist)
      call transform_old2new(dsdtLast_LMloc, dsdtLast_LMdist)
      
      !!! END New layout  TMP
      !!!----------------------------------
      

      PERFON('LMloop')
      !LIKWID_ON('LMloop')
      if ( lVerbose .and. l_save_out ) then
         open(newunit=n_log_file, file=log_file, status='unknown', &
         &    position='append')
      end if

      omega_icLast=omega_ic

      if ( lMat ) then ! update matrices:
      !---- The following logicals tell whether the respective inversion
      !     matrices have been updated. lMat=.true. when a general
      !     update is necessary. These logicals are THREADPRIVATE and
      !     stored in the module matrices in m_mat.F90:
         lZ10mat=.false.
         do l=0,l_max
            if ( l_single_matrix ) then
               if ( l_TP_form .or. l_anelastic_liquid ) then
                  lWPTmat(l)=.false.
               else
                  lWPSmat(l)=.false.
               end if
            else
               lWPmat(l)=.false.
               lSmat(l) =.false.
            end if
            lZmat(l) =.false.
            if ( l_mag ) lBmat(l) =.false.
            if ( l_chemical_conv ) lXimat(l)=.false.
         end do
         
         ! NeEEEEEEEEEWWWWWWWWWW
         if ( l_single_matrix ) then
            continue
         else
            lZ10mat_new = .false.
            lSmat_new(:) =.false.
            lZmat_new(:) =.false.
            lWPmat_new(:)=.false.
         end if
      end if
      

      !nThreadsLMmax = 1
      nLMB=1+coord_r
      !nTh=1
      if ( lVerbose ) then
         write(*,'(/," ! lm block no:",i3)') nLMB
         call wallTime(tStart)
      end if

      if ( l_heat ) then ! dp,workA usead as work arrays
         if ( DEBUG_OUTPUT ) then
            write(*,"(A,I2,8ES20.12)") "s_before: ",nLMB,   &
                 & GET_GLOBAL_SUM( s_LMloc(:,:) ),          &
                 & GET_GLOBAL_SUM( ds_LMloc(:,:) ),         &
                 & GET_GLOBAL_SUM( dsdt(:,:) ),             &
                 & GET_GLOBAL_SUM( dsdtLast_LMloc(:,:) )
         end if
         if ( .not. l_single_matrix ) then
! ! !             PERFON('up_S')
            if ( l_anelastic_liquid ) then
            
               print *, "Not Parallelized!", __LINE__, __FILE__
               stop
               call updateS_ala(s_LMloc, ds_LMloc, w_LMloc, dVSrLM,dsdt,  & 
                    &           dsdtLast_LMloc, w1, coex, dt, nLMB)
            else
               
               PERFON('par_old')
               PERFON('upS_old')
               call updateS( s_LMloc, ds_LMloc, w_LMloc, dVSrLM,dsdt, &
                    &        dsdtLast_LMloc, w1, coex, dt, nLMB )
               PERFOFF
               PERFON('trp_old')
               call lo2r_redist_start_dist(lo2r_s,s_LMloc_container,s_Rdist_test)
               call lo2r_redist_wait_dist(lo2r_s)
               PERFOFF
               PERFOFF
               
               PERFON('par_new')
               PERFON('upS_new')
               call updateS_new( s_LMdist, ds_LMdist, w_LMdist, dVSrLM_dist, dsdt_dist, &
                    &             dsdtLast_LMdist, w1, coex, dt )
               PERFOFF
               PERFON('trp_new')
               call ml2r_s%start(s_LMdist_container, s_Rdist_container, 2)
               call ml2r_s%wait()
               PERFOFF
               PERFOFF
               
               call test_field(s_LMdist   , s_LMloc , "s")
               call test_field(ds_LMdist  , ds_LMloc, "ds")
               call test_field(w_LMdist   , w_LMloc , "w")
               call test_field(dsdt_dist  , dsdt    , "dsdt")
               call test_field(dVSrLM_dist, dVSrLM  , "dVSrLM")
               call test_field(dsdtLast_LMdist, dsdtLast_LMloc, "dsdtLast")
            end if
! ! ! !             PERFOFF
            
!             Here one could start the redistribution of s_LMloc,ds_LMloc etc. with a 
!             nonblocking send
!             PERFON('rdstSst')

!             PERFON('trp_old')
!             call lo2r_redist_start_dist(lo2r_s,s_LMloc_container,s_Rdist_test)
!             call lo2r_redist_wait_dist(lo2r_s)
!             PERFOFF
            
!             PERFON('trp_new')
!             call ml2r_s%start(s_LMdist_container, s_Rdist_container, 2)
!             call ml2r_s%wait()
!             PERFOFF
!             test_norm  = SUM(ABS(REAL(s_Rdist_test) - REAL(s_Rdist_container)))
!             test_normi  = SUM(ABS(AIMAG(s_Rdist_test) - AIMAG(s_Rdist_container)))
!             IF (test_norm+test_normi>error_threshold) print *, "|| cont || : ", test_norm + test_normi
            
!             PERFOFF
         end if

         if ( DEBUG_OUTPUT ) then
            write(*,"(A,I2,4ES20.12)") "s_after : ",nLMB,  &
                 & get_global_sum_dist( s_LMdist(:,:) ),         &
                 & get_global_sum_dist( ds_LMdist(:,:) )
            write(*,"(A,I2,8ES22.14)") "s_after(bnd_r): ",nLMB, &
                 & get_global_sum_dist( s_LMdist(:,n_r_icb) ),        &
                 & get_global_sum_dist( s_LMdist(:,n_r_cmb) ),        &
                 & get_global_sum_dist( ds_LMdist(:,n_r_icb) ),       &
                 & get_global_sum_dist( ds_LMdist(:,n_r_cmb) )
         end if
      end if

      if ( l_chemical_conv ) then ! dp,workA usead as work arrays
         call updateXi(xi_LMloc,dxi_LMloc,dVXirLM,dxidt,dxidtLast_LMloc, &
              &        w1,coex,dt,nLMB)

         call lo2r_redist_start_dist(lo2r_xi,xi_LMloc_container,xi_Rdist_container)
      end if
      
      if ( l_conv ) then
         if ( DEBUG_OUTPUT ) then
            write(*,"(A,I2,6ES20.12)") "z_before: ",nLMB,   &
                 & GET_GLOBAL_SUM( z_LMloc(:,:) ),          &
                 & GET_GLOBAL_SUM( dz_LMloc(:,:) ),         &
                 & GET_GLOBAL_SUM( dzdtLast_lo(:,:) )
         end if
         
         
         call transform_old2new(z_LMloc, z_LMdist)
         call transform_old2new(dz_LMloc, dz_LMdist)
         call transform_old2new(dzdt, dzdt_dist)
         call transform_old2new(dzdtLast_lo, dzdtLast_lodist)

         omega_ma_dist = omega_ma
         omega_ic_dist = omega_ic
         d_omega_ma_dtLast_dist = d_omega_ma_dtLast
         d_omega_ic_dtLast_dist = d_omega_ic_dtLast
         
         ! dp, dVSrLM, workA used as work arrays
!          PERFON('upZ_old')
!          call updateZ( z_LMloc, dz_LMloc, dzdt, dzdtLast_lo, time, &
!               &        omega_ma,d_omega_ma_dtLast,                 &
!               &        omega_ic,d_omega_ic_dtLast,                 &
!               &        lorentz_torque_ma,lorentz_torque_maLast,    &
!               &        lorentz_torque_ic,lorentz_torque_icLast,    &
!               &        w1,coex,dt,lRmsNext )
!          PERFOFF
         
         PERFON('upZ_new')
         call updateZ_new( z_LMdist, dz_LMdist, dzdt_dist, dzdtLast_lodist, time, &
              &        omega_ma_dist,d_omega_ma_dtLast_dist,                 &
              &        omega_ic_dist,d_omega_ic_dtLast_dist,                 &
              &        lorentz_torque_ma,lorentz_torque_maLast,    &
              &        lorentz_torque_ic,lorentz_torque_icLast,    &
              &        w1,coex,dt,lRmsNext ) 
         PERFOFF
         
         call test_field(z_LMdist        , z_LMloc    , "z_")
         call test_field(dz_LMdist       , dz_LMloc   , "dz_")
         call test_field(dzdtLast_lodist, dzdtLast_lo, "dzdtLast_lo_")
         
         call transform_new2old(z_LMdist, z_LMloc)
         call transform_new2old(dz_LMdist, dz_LMloc)
         call transform_new2old(dzdtLast_lodist,dzdtLast_lo)
              

         !call MPI_Barrier(comm_r,ierr)

         if ( DEBUG_OUTPUT ) then
            write(*,"(A,I2,6ES20.12)") "z_after: ",nLMB,  &
                 & GET_GLOBAL_SUM( z_LMloc(:,:) ),        &
                 & GET_GLOBAL_SUM( dz_LMloc(:,:) ),       &
                 & GET_GLOBAL_SUM( dzdtLast_lo(:,:) )
         end if
         ! dVSrLM, workA used as work arrays
         if ( DEBUG_OUTPUT ) then
            sum_dwdt=GET_GLOBAL_SUM( dwdt(:,:) )
            write(*,"(A,I2,8ES22.14,4(I3,F19.16))") "wp_before: ",nLMB,&
                 & GET_GLOBAL_SUM( w_LMloc(:,:) ),                     &
                 & GET_GLOBAL_SUM( p_LMloc(:,:) ),                     &
                 & GET_GLOBAL_SUM( dwdtLast_LMloc(:,:) ),              &
                 & GET_GLOBAL_SUM( dpdtLast_LMloc(:,:) ),              &
                 & exponent(real(sum_dwdt)),fraction(real(sum_dwdt)),  &
                 & exponent(aimag(sum_dwdt)),fraction(aimag(sum_dwdt))
         end if

         if ( l_single_matrix ) then
            if ( coord_r == rank_with_l1m0 ) then
               do nR=1,n_r_max
                  z10(nR)=real(z_LMloc(lo_map%lm2(1,0),nR))
               end do
            end if
#ifdef WITH_MPI
            call MPI_Bcast(z10,n_r_max,MPI_DEF_REAL,rank_with_l1m0, &
                 &         comm_r,ierr)
#endif
            if ( l_TP_form ) then
               call updateWPT( w_LMloc, dw_LMloc, ddw_LMloc, z10, dwdt,     &
                 &             dwdtLast_LMloc, p_LMloc, dp_LMloc, dpdt,     &
                 &             dpdtLast_LMloc, s_LMloc, ds_LMloc, dVSrLM,   &
                 &             dVPrLM, dsdt, dsdtLast_LMloc, w1, coex, dt,  &
                 &             nLMB, lRmsNext )
            else
               call updateWPS( w_LMloc, dw_LMloc, ddw_LMloc, z10, dwdt,    &
                 &             dwdtLast_LMloc, p_LMloc, dp_LMloc, dpdt,    &
                 &             dpdtLast_LMloc, s_LMloc, ds_LMloc, dVSrLM,  &
                 &             dsdt, dsdtLast_LMloc, w1, coex, dt, nLMB,   &
                 &             lRmsNext )
            end if

            call lo2r_redist_start_dist(lo2r_s,s_LMloc_container,s_Rdist_container)
         else
            dVxVhLM_dist = 0.0
            call transform_old2new( w_LMloc       , w_LMdist        )
            call transform_old2new( dw_LMloc      , dw_LMdist       )
            call transform_old2new( ddw_LMloc     , ddw_LMdist      )
            call transform_old2new( dwdt          , dwdt_dist       )    !
            call transform_old2new( dwdtLast_LMloc, dwdtLast_LMdist )    !
            call transform_old2new( p_LMloc       , p_LMdist        )
            call transform_old2new( dp_LMloc      , dp_LMdist       )
            call transform_old2new( dpdt          , dpdt_dist       )    !
            call transform_old2new( dpdtLast_LMloc, dpdtLast_LMdist )    !
            call transform_old2new( s_LMloc       , s_LMdist        )
            if (l_double_curl)   call transform_old2new( dVxVhLM       , dVxVhLM_dist    )    !
            if (l_chemical_conv) call transform_old2new( xi_LMloc      , xi_LMdist       )
            
            PERFON('upWP_new')
            call updateWP_new( w_LMdist, dw_LMdist, ddw_LMdist, dVxVhLM_dist, dwdt_dist,     &
                 &         dwdtLast_LMdist, p_LMdist, dp_LMdist, dpdt_dist,         &
                 &         dpdtLast_LMdist, s_LMdist, xi_LMdist, w1, coex, dt, &
                 &         lRmsNext, lPressNext)
            PERFOFF
                 
                 
!             PERFON('upWP_old')
!             call updateWP( w_LMloc, dw_LMloc, ddw_LMloc, dVxVhLM, dwdt,     &
!                  &         dwdtLast_LMloc, p_LMloc, dp_LMloc, dpdt,         &
!                  &         dpdtLast_LMloc, s_LMloc, xi_LMloc, w1, coex, dt, &
!                  &         nLMB, lRmsNext, lPressNext)
!             PERFOFF

            call transform_new2old( w_LMdist        , w_LMloc       )
            call transform_new2old( dw_LMdist       , dw_LMloc      )
            call transform_new2old( ddw_LMdist      , ddw_LMloc     )
            call transform_new2old( dwdt_dist       , dwdt          )    !
            call transform_new2old( dwdtLast_LMdist , dwdtLast_LMloc)    !
            call transform_new2old( p_LMdist        , p_LMloc       )
            call transform_new2old( dp_LMdist       , dp_LMloc      )
            call transform_new2old( dpdtLast_LMdist , dpdtLast_LMloc)    !
            call transform_new2old( s_LMdist        , s_LMloc       )
            if (l_double_curl)   call transform_new2old( dVxVhLM_dist, dVxVhLM)    !
            if (l_chemical_conv) call transform_new2old( xi_LMdist   , xi_LMloc )
            
            call test_field(w_LMdist       , w_LMloc       , "w_")
            call test_field(dw_LMdist      , dw_LMloc      , "dw_")
            call test_field(ddw_LMdist     , ddw_LMloc     , "ddw_")
            call test_field(dwdt_dist      , dwdt          , "dwdt_")
            call test_field(dwdtLast_LMdist, dwdtLast_LMloc, "dwdtLast_")
            call test_field(p_LMdist       , p_LMloc       , "p_")
            call test_field(dp_LMdist      , dp_LMloc      , "dp_")
            call test_field(dpdtLast_LMdist, dpdtLast_LMloc, "dpdtLast_")
            call test_field(s_LMdist       , s_LMloc       , "s_")
            if (l_double_curl)   call test_field(dVxVhLM_dist   , dVxVhLM       , "dVxVhLM_")
            if (l_chemical_conv) call test_field(xi_LMdist      , xi_LMloc      , "xi_")
            

            if ( DEBUG_OUTPUT ) then
               write(*,"(A,I2,12ES22.14)") "wp_after: ",nLMB,  &
                    & GET_GLOBAL_SUM( w_LMloc(:,:) ),          &
                    & GET_GLOBAL_SUM( p_LMloc(:,:) ),          &
                    & GET_GLOBAL_SUM( dwdtLast_LMloc(:,:) ),   &
                    & GET_GLOBAL_SUM( dpdtLast_LMloc(:,:) ),   &
                    &GET_GLOBAL_SUM( dw_LMloc(:,:) )
               write(*,"(A,I2,4ES22.14)") "wp_after(bnd_r): ",nLMB, &
                    & GET_GLOBAL_SUM( w_LMloc(:,n_r_icb) ),         &
                    & GET_GLOBAL_SUM( w_LMloc(:,n_r_cmb) )
            end if
         end if
         call lo2r_redist_start_dist(lo2r_flow,flow_LMloc_container,flow_Rdist_container)
      end if
      if ( l_mag ) then ! dwdt,dpdt used as work arrays
         if ( DEBUG_OUTPUT ) then
            write(*,"(A,I2,14ES20.12)") "b_before: ",nLMB,   &
                 & GET_GLOBAL_SUM(  b_LMloc(:,:) ),          & 
                 & GET_GLOBAL_SUM( aj_LMloc(:,:) ),          &
                 & GET_GLOBAL_SUM( b_ic_LMloc(:,:) ),        &
                 & GET_GLOBAL_SUM( aj_ic_LMloc(:,:) ),       &
                 & GET_GLOBAL_SUM( dbdt_icLast_LMloc(:,:) ), &
                 & GET_GLOBAL_SUM( djdt_icLast_LMloc(:,:) ), &
                 & GET_GLOBAL_SUM( dVxBhLM(:,:) )
         end if
         !LIKWID_ON('up_B')
         PERFON('up_B')
         call updateB( b_LMloc,db_LMloc,ddb_LMloc,aj_LMloc,dj_LMloc,ddj_LMloc, &
              &        dVxBhLM, dbdt, dbdtLast_LMloc, djdt, djdtLast_LMloc,    &
              &        b_ic_LMloc, db_ic_LMloc, ddb_ic_LMloc, aj_ic_LMloc,     &
              &        dj_ic_LMloc, ddj_ic_LMloc, dbdt_icLast_LMloc,           &
              &        djdt_icLast_LMloc, b_nl_cmb, aj_nl_cmb, aj_nl_icb,      &
              &        omega_icLast, w1, coex, dt, time, nLMB, lRmsNext )
         PERFOFF
         !LIKWID_OFF('up_B')
         call lo2r_redist_start_dist(lo2r_field,field_LMloc_container,field_Rdist_container)

         if ( DEBUG_OUTPUT ) then
            write(*,"(A,I2,8ES20.12)") "b_after: ",nLMB, &
                 & GET_GLOBAL_SUM(  b_LMloc(:,:) ),      & 
                 & GET_GLOBAL_SUM( aj_LMloc(:,:) ),      &
                 & GET_GLOBAL_SUM( dbdtLast_LMloc(:,:) ),&
                 & GET_GLOBAL_SUM( djdtLast_LMloc(:,:) )
            write(*,"(A,I2,8ES20.12)") "b_ic_after: ",nLMB, &
                 & GET_GLOBAL_SUM( b_ic_LMloc(:,:) ),       &
                 & GET_GLOBAL_SUM( aj_ic_LMloc(:,:) ),      &
                 & GET_GLOBAL_SUM( dbdt_icLast_LMloc(:,:) ),&
                 & GET_GLOBAL_SUM( djdt_icLast_LMloc(:,:) )
         end if
      end if

      if ( lVerbose ) then
         call wallTime(tStop)
         call subTime(tStart,tStop,tPassed)
         call writeTime(n_log_file,'! Time for thread:',tPassed)
      end if


      lorentz_torque_maLast=lorentz_torque_ma
      lorentz_torque_icLast=lorentz_torque_ic

      if ( lVerbose .and. l_save_out ) close(n_log_file)

      !LIKWID_OFF('LMloop')
      PERFOFF
   end subroutine LMLoop
   
   
!--------------------------------------------------------------------------------
end module LMLoop_mod
