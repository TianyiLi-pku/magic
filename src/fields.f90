module fields
   !
   !  This module contains the potential fields and their radial
   !  derivatives
   !
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use geometry, only: lm_max, n_r_max, lm_maxMag, n_r_maxMag, &
       &                 n_r_ic_maxMag, n_mlo_loc, n_mloMag_loc, &
       &                 n_lm_loc, n_lmMag_loc, l_r, u_r 
   use logic, only: l_chemical_conv
   use LMLoop_data, only: llm, ulm, llmMag, ulmMag
   use parallel_mod, only: coord_r, n_ranks_theta
 
   implicit none

   private
 
   !-- Velocity potentials:
   complex(cp), public, allocatable, target :: flow_LMloc_container(:,:,:)
   complex(cp), public, pointer :: w_LMloc(:,:),dw_LMloc(:,:),ddw_LMloc(:,:)
   complex(cp), public, allocatable, target :: flow_Rdist_container(:,:,:)
   complex(cp), public, pointer :: w_Rdist(:,:), dw_Rdist(:,:), ddw_Rdist(:,:)
 
   complex(cp), public, pointer :: z_LMloc(:,:),dz_LMloc(:,:)
   complex(cp), public, pointer :: z_Rdist(:,:), dz_Rdist(:,:)
 
   !-- Entropy:
   complex(cp), public, allocatable, target :: s_LMloc_container(:,:,:)
   complex(cp), public, allocatable, target :: s_Rdist_container(:,:,:)
   complex(cp), public, allocatable, target :: s_Rdist_test(:,:,:)  !DELETEME
   complex(cp), public, pointer :: s_LMloc(:,:), ds_LMloc(:,:)
   complex(cp), public, pointer :: s_Rdist(:,:), ds_Rdist(:,:)
 
   !-- Chemical composition:
   complex(cp), public, allocatable, target :: xi_LMloc_container(:,:,:)
   complex(cp), public, allocatable, target :: xi_Rdist_container(:,:,:)
   complex(cp), public, pointer :: xi_LMloc(:,:), dxi_LMloc(:,:)
   complex(cp), public, pointer :: xi_Rdist(:,:), dxi_Rdist(:,:)

   !-- Pressure:
   complex(cp), public, pointer :: p_LMloc(:,:), dp_LMloc(:,:)
   complex(cp), public, pointer :: p_Rdist(:,:), dp_Rdist(:,:)
 
   !-- Magnetic field potentials:
   complex(cp), public, allocatable :: b(:,:)
   complex(cp), public, allocatable, target :: field_LMloc_container(:,:,:)
   complex(cp), public, allocatable, target :: field_Rdist_container(:,:,:)
   complex(cp), public, pointer :: b_LMloc(:,:), db_LMloc(:,:), ddb_LMloc(:,:)
   complex(cp), public, pointer :: b_Rdist(:,:), db_Rdist(:,:), ddb_Rdist(:,:)
   complex(cp), public, pointer :: aj_LMloc(:,:), dj_LMloc(:,:), ddj_LMloc(:,:)
   complex(cp), public, pointer :: aj_Rdist(:,:), dj_Rdist(:,:)

   !-- Magnetic field potentials in inner core:
   !   NOTE: n_r_loc-dimension may be smaller once CHEBFT is addopted
   !         for even chebs
   complex(cp), public, allocatable :: b_ic(:,:)  
   complex(cp), public, allocatable :: db_ic(:,:)
   complex(cp), public, allocatable :: ddb_ic(:,:)
   complex(cp), public, allocatable :: aj_ic(:,:) 
   complex(cp), public, allocatable :: dj_ic(:,:)
   complex(cp), public, allocatable :: ddj_ic(:,:)
   complex(cp), public, allocatable :: b_ic_LMloc(:,:)  
   complex(cp), public, allocatable :: db_ic_LMloc(:,:)
   complex(cp), public, allocatable :: ddb_ic_LMloc(:,:)
   complex(cp), public, allocatable :: aj_ic_LMloc(:,:) 
   complex(cp), public, allocatable :: dj_ic_LMloc(:,:)
   complex(cp), public, allocatable :: ddj_ic_LMloc(:,:)

   complex(cp), public, allocatable :: work_LMloc(:,:) ! Needed in update routines
   
   
   !![NEW LAYOUT]
   !-- Velocity potentials:
   complex(cp), public, allocatable, target :: flow_LMdist_container(:,:,:)
   complex(cp), public, pointer :: w_LMdist(:,:),dw_LMdist(:,:),ddw_LMdist(:,:)
 
   complex(cp), public, pointer :: z_LMdist(:,:),dz_LMdist(:,:)
 
   !-- Entropy:
   complex(cp), public, allocatable, target :: s_LMdist_container(:,:,:)
   complex(cp), public, pointer :: s_LMdist(:,:), ds_LMdist(:,:)
 
   !-- Chemical composition:
   complex(cp), public, allocatable, target :: xi_LMdist_container(:,:,:)
   complex(cp), public, pointer :: xi_LMdist(:,:), dxi_LMdist(:,:)

   !-- Pressure:
   complex(cp), public, pointer :: p_LMdist(:,:), dp_LMdist(:,:)
 
   !-- Magnetic field potentials:
   complex(cp), public, allocatable, target :: field_LMdist_container(:,:,:)
   complex(cp), public, pointer :: b_LMdist(:,:), db_LMdist(:,:), ddb_LMdist(:,:)
   complex(cp), public, pointer :: aj_LMdist(:,:), dj_LMdist(:,:), ddj_LMdist(:,:)
 
   !-- Magnetic field potentials in inner core:
   !   NOTE: n_r_loc-dimension may be smaller once CHEBFT is addopted
   !         for even chebs
   complex(cp), public, allocatable :: b_ic_LMdist(:,:)  
   complex(cp), public, allocatable :: db_ic_LMdist(:,:)
   complex(cp), public, allocatable :: ddb_ic_LMdist(:,:)
   complex(cp), public, allocatable :: aj_ic_LMdist(:,:) 
   complex(cp), public, allocatable :: dj_ic_LMdist(:,:)
   complex(cp), public, allocatable :: ddj_ic_LMdist(:,:)

   complex(cp), public, allocatable :: work_LMdist(:,:) ! Needed in update routines
   !![END NEW LAYOUT]
   
   
   !-- Rotation rates:
   real(cp), public :: omega_ic,omega_ma

   public :: initialize_fields, finalize_fields

contains

   subroutine initialize_fields

      !-- Velocity potentials:
      if ( coord_r == 0 ) then
         allocate( b(lm_maxMag,n_r_maxMag) )
         bytes_allocated = bytes_allocated +  &
                           lm_maxMag*n_r_maxMag*SIZEOF_DEF_COMPLEX
         allocate( b_ic(lm_maxMag,n_r_ic_maxMag) )  
         allocate( db_ic(lm_maxMag,n_r_ic_maxMag) )
         allocate( ddb_ic(lm_maxMag,n_r_ic_maxMag) )
         allocate( aj_ic(lm_maxMag,n_r_ic_maxMag) ) 
         allocate( dj_ic(lm_maxMag,n_r_ic_maxMag) )
         allocate( ddj_ic(lm_maxMag,n_r_ic_maxMag) )
         bytes_allocated = bytes_allocated + &
                           6*lm_maxMag*n_r_ic_maxMag*SIZEOF_DEF_COMPLEX

      else
         allocate( b(1,n_r_maxMag) )
         bytes_allocated = bytes_allocated + n_r_maxMag*SIZEOF_DEF_COMPLEX
         allocate( b_ic(1,n_r_ic_maxMag) )  
         allocate( db_ic(1,n_r_ic_maxMag) )
         allocate( ddb_ic(1,n_r_ic_maxMag) )
         allocate( aj_ic(1,n_r_ic_maxMag) ) 
         allocate( dj_ic(1,n_r_ic_maxMag) )
         allocate( ddj_ic(1,n_r_ic_maxMag) )
         bytes_allocated = bytes_allocated + 6*n_r_ic_maxMag*SIZEOF_DEF_COMPLEX
      end if
      allocate( flow_LMloc_container(llm:ulm,n_r_max,1:7) )
      w_LMloc(llm:ulm,1:n_r_max)   => flow_LMloc_container(:,:,1)
      dw_LMloc(llm:ulm,1:n_r_max)  => flow_LMloc_container(:,:,2)
      ddw_LMloc(llm:ulm,1:n_r_max) => flow_LMloc_container(:,:,3)
      z_LMloc(llm:ulm,1:n_r_max)   => flow_LMloc_container(:,:,4)
      dz_LMloc(llm:ulm,1:n_r_max)  => flow_LMloc_container(:,:,5)
      p_LMloc(llm:ulm,1:n_r_max)   => flow_LMloc_container(:,:,6)
      dp_LMloc(llm:ulm,1:n_r_max)  => flow_LMloc_container(:,:,7)

      allocate( flow_Rdist_container(n_lm_loc,l_r:u_r,1:7) )
      w_Rdist(1:n_lm_loc,l_r:u_r)   => flow_Rdist_container(:,:,1)
      dw_Rdist(1:n_lm_loc,l_r:u_r)  => flow_Rdist_container(:,:,2)
      ddw_Rdist(1:n_lm_loc,l_r:u_r) => flow_Rdist_container(:,:,3)
      z_Rdist(1:n_lm_loc,l_r:u_r)   => flow_Rdist_container(:,:,4)
      dz_Rdist(1:n_lm_loc,l_r:u_r)  => flow_Rdist_container(:,:,5)
      p_Rdist(1:n_lm_loc,l_r:u_r)   => flow_Rdist_container(:,:,6)
      dp_Rdist(1:n_lm_loc,l_r:u_r)  => flow_Rdist_container(:,:,7)

      !-- Entropy:
      allocate( s_LMloc_container(llm:ulm,n_r_max,1:2) )
      s_LMloc(llm:ulm,1:n_r_max)  => s_LMloc_container(:,:,1)
      ds_LMloc(llm:ulm,1:n_r_max) => s_LMloc_container(:,:,2)
      allocate( s_Rdist_container(n_lm_loc,l_r:u_r,1:2) )
      s_Rdist(1:n_lm_loc,l_r:u_r)   => s_Rdist_container(:,:,1)
      ds_Rdist(1:n_lm_loc,l_r:u_r)  => s_Rdist_container(:,:,2)

      allocate( s_Rdist_test(n_lm_loc,l_r:u_r,1:2) )  !DELETEME
      
      bytes_allocated = bytes_allocated + &
                        9*(ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
      bytes_allocated = bytes_allocated + &
                        9*n_lm_loc*(u_r-l_r+1)*SIZEOF_DEF_COMPLEX

      !-- Chemical composition:
      if ( l_chemical_conv ) then
         allocate( xi_LMloc_container(llm:ulm,n_r_max,1:2) )
         xi_LMloc(llm:ulm,1:n_r_max)  => xi_LMloc_container(:,:,1)
         dxi_LMloc(llm:ulm,1:n_r_max) => xi_LMloc_container(:,:,2)
         allocate( xi_Rdist_container(n_lm_loc,l_r:u_r,1:2) )
         xi_Rdist(1:n_lm_loc,l_r:u_r)   => xi_Rdist_container(:,:,1)
         dxi_Rdist(1:n_lm_loc,l_r:u_r)  => xi_Rdist_container(:,:,2)
         bytes_allocated = bytes_allocated + &
                           2*(ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
         bytes_allocated = bytes_allocated + &
                           2*n_lm_loc*(u_r-l_r+1)*SIZEOF_DEF_COMPLEX
      else
         allocate( xi_LMloc_container(1,1,2) ) ! For debugging
         xi_LMloc(1:1,1:1)  => xi_LMloc_container(:,:,1)
         dxi_LMloc(1:1,1:1) => xi_LMloc_container(:,:,2)
         allocate( xi_Rdist_container(1,1,2) )
         xi_Rdist(1:1,1:1)   => xi_Rdist_container(:,:,1)
         dxi_Rdist(1:1,1:1)  => xi_Rdist_container(:,:,2)
      end if

      !-- Magnetic field potentials:
      allocate( field_LMloc_container(llmMag:ulmMag,n_r_maxMag,1:6) )
      b_LMloc(llmMag:ulmMag,1:n_r_maxMag)   => field_LMloc_container(:,:,1)
      db_LMloc(llmMag:ulmMag,1:n_r_maxMag)  => field_LMloc_container(:,:,2)
      ddb_LMloc(llmMag:ulmMag,1:n_r_maxMag) => field_LMloc_container(:,:,3)
      aj_LMloc(llmMag:ulmMag,1:n_r_maxMag)  => field_LMloc_container(:,:,4)
      dj_LMloc(llmMag:ulmMag,1:n_r_maxMag)  => field_LMloc_container(:,:,5)
      ddj_LMloc(llmMag:ulmMag,1:n_r_maxMag) => field_LMloc_container(:,:,6)
      bytes_allocated = bytes_allocated + &
                        6*(ulmMag-llmMag+1)*n_r_maxMag*SIZEOF_DEF_COMPLEX

      allocate( field_Rdist_container(n_lmMag_loc,l_r:u_r,1:5) )
      b_Rdist(1:n_lmMag_loc,l_r:u_r)   => field_Rdist_container(:,:,1)
      db_Rdist(1:n_lmMag_loc,l_r:u_r)  => field_Rdist_container(:,:,2)
      ddb_Rdist(1:n_lmMag_loc,l_r:u_r) => field_Rdist_container(:,:,3)
      aj_Rdist(1:n_lmMag_loc,l_r:u_r)  => field_Rdist_container(:,:,4)
      dj_Rdist(1:n_lmMag_loc,l_r:u_r)  => field_Rdist_container(:,:,5)
      bytes_allocated = bytes_allocated + &
                        5*n_lmMag_loc*(u_r-l_r+1)*SIZEOF_DEF_COMPLEX

      !-- Magnetic field potentials in inner core:
      !   NOTE: n_r_loc-dimension may be smaller once CHEBFT is adopted
      !         for even chebs
      allocate( b_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )  
      allocate( db_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )
      allocate( ddb_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )
      allocate( aj_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) ) 
      allocate( dj_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )
      allocate( ddj_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )
      bytes_allocated = bytes_allocated + &
                        6*(ulmMag-llmMag+1)*n_r_ic_maxMag*SIZEOF_DEF_COMPLEX
      
      allocate( work_LMloc(llm:ulm,1:n_r_max) )
      bytes_allocated = bytes_allocated + (ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
      
      call initialize_fields_dist(bytes_allocated)

   end subroutine initialize_fields
!----------------------------------------------------------------------------
   subroutine initialize_fields_dist(bytes)
   
      integer(lip), intent(inout) :: bytes

      allocate( flow_LMdist_container(1:n_mlo_loc,n_r_max,1:7) )
      w_LMdist(1:n_mlo_loc,1:n_r_max)   => flow_LMdist_container(:,:,1)
      dw_LMdist(1:n_mlo_loc,1:n_r_max)  => flow_LMdist_container(:,:,2)
      ddw_LMdist(1:n_mlo_loc,1:n_r_max) => flow_LMdist_container(:,:,3)
      z_LMdist(1:n_mlo_loc,1:n_r_max)   => flow_LMdist_container(:,:,4)
      dz_LMdist(1:n_mlo_loc,1:n_r_max)  => flow_LMdist_container(:,:,5)
      p_LMdist(1:n_mlo_loc,1:n_r_max)   => flow_LMdist_container(:,:,6)
      dp_LMdist(1:n_mlo_loc,1:n_r_max)  => flow_LMdist_container(:,:,7)

      !-- Entropy:
      allocate( s_LMdist_container(1:n_mlo_loc,n_r_max,1:2) )
      s_LMdist(1:n_mlo_loc,1:n_r_max)  => s_LMdist_container(:,:,1)
      ds_LMdist(1:n_mlo_loc,1:n_r_max) => s_LMdist_container(:,:,2)

      bytes = bytes + 9*(n_mlo_loc)*n_r_max*SIZEOF_DEF_COMPLEX

      !-- Chemical composition:
      if ( l_chemical_conv ) then
         allocate( xi_LMdist_container(1:n_mlo_loc,n_r_max,1:2) )
         xi_LMdist(1:n_mlo_loc,1:n_r_max)  => xi_LMdist_container(:,:,1)
         dxi_LMdist(1:n_mlo_loc,1:n_r_max) => xi_LMdist_container(:,:,2)
         bytes = bytes + 2*(n_mlo_loc)*n_r_max*SIZEOF_DEF_COMPLEX
      else
         allocate( xi_LMdist_container(1,1,2) ) ! For debugging
         xi_LMdist(1:1,1:1)  => xi_LMdist_container(:,:,1)
         dxi_LMdist(1:1,1:1) => xi_LMdist_container(:,:,2)
      end if

      !-- Magnetic field potentials:
      allocate( field_LMdist_container(1:n_mloMag_loc,n_r_maxMag,1:6) )
      b_LMdist(1:n_mloMag_loc,1:n_r_maxMag)   => field_LMdist_container(:,:,1)
      db_LMdist(1:n_mloMag_loc,1:n_r_maxMag)  => field_LMdist_container(:,:,2)
      ddb_LMdist(1:n_mloMag_loc,1:n_r_maxMag) => field_LMdist_container(:,:,3)
      aj_LMdist(1:n_mloMag_loc,1:n_r_maxMag)  => field_LMdist_container(:,:,4)
      dj_LMdist(1:n_mloMag_loc,1:n_r_maxMag)  => field_LMdist_container(:,:,5)
      ddj_LMdist(1:n_mloMag_loc,1:n_r_maxMag) => field_LMdist_container(:,:,6)
      bytes = bytes + 6*(n_mloMag_loc)*n_r_maxMag*SIZEOF_DEF_COMPLEX

      !-- Magnetic field potentials in inner core:
      !   NOTE: n_r_loc-dimension may be smaller once CHEBFT is adopted
      !         for even chebs
      allocate( b_ic_LMdist(1:n_mloMag_loc,n_r_ic_maxMag) )  
      allocate( db_ic_LMdist(1:n_mloMag_loc,n_r_ic_maxMag) )
      allocate( ddb_ic_LMdist(1:n_mloMag_loc,n_r_ic_maxMag) )
      allocate( aj_ic_LMdist(1:n_mloMag_loc,n_r_ic_maxMag) ) 
      allocate( dj_ic_LMdist(1:n_mloMag_loc,n_r_ic_maxMag) )
      allocate( ddj_ic_LMdist(1:n_mloMag_loc,n_r_ic_maxMag) )
      bytes = bytes + 6*(n_mloMag_loc)*n_r_ic_maxMag*SIZEOF_DEF_COMPLEX
      
      allocate( work_LMdist(1:n_mlo_loc,1:n_r_max) )
      bytes = bytes + (n_mlo_loc)*n_r_max*SIZEOF_DEF_COMPLEX

   end subroutine initialize_fields_dist
!----------------------------------------------------------------------------
   subroutine finalize_fields

      nullify( xi_Rdist )
      nullify( dxi_Rdist)
      nullify( s_Rdist  )
      nullify( ds_Rdist )
      nullify( z_Rdist  )
      nullify( dz_Rdist )
      nullify( p_Rdist  )
      nullify( dp_Rdist )
      nullify( b_Rdist  )
      nullify( db_Rdist )
      nullify( ddb_Rdist)
      nullify( aj_Rdist )
      nullify( dj_Rdist )
      nullify( w_Rdist  )
      nullify( dw_Rdist )
      nullify( ddw_Rdist)
      
      !!!!!!!!!!!! I think I need to make sure that there are no hanging
      !!!!!!!!!!!! MPI requests involving these arrays as buffers!
      !!!!!!!!!!!! I'm getting weird errors
      deallocate( b, b_ic, db_ic, ddb_ic )
      deallocate( aj_ic, dj_ic, ddj_ic)
      deallocate( flow_LMloc_container )
      deallocate( flow_LMdist_container )
      deallocate( flow_Rdist_container )
      deallocate( s_LMloc_container )
      deallocate( s_LMdist_container )
      deallocate( s_Rdist_container )
!       deallocate( b_ic_LMloc )
!       deallocate( db_ic_LMloc )
!       deallocate( ddb_ic_LMloc )
!       deallocate( aj_ic_LMloc )
!       deallocate( dj_ic_LMloc, ddj_ic_LMloc )
!       deallocate( xi_LMloc_container, xi_Rdist_container )
!       deallocate( field_LMloc_container, field_Rdist_container )
!       deallocate( work_LMloc )

   end subroutine finalize_fields
!----------------------------------------------------------------------------
end module fields
