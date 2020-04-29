module module_fim_cpl_finalize
  use module_fim_cpl_run,only: u_tdcy_phy, v_tdcy_phy, trc_tdcy_phy,       &
                               gfs_u_old, gfs_v_old, gfs_t_old, gfs_q_old, &
                               gfs_oz_old, gfs_cld_old
  implicit none

contains

!*********************************************************************
!       Finish the FIM DYN-PHY coupler component.  
!       T. Henderson            February, 2009
!*********************************************************************
  subroutine cpl_finalize
    ! deallocate coupler state
    deallocate( u_tdcy_phy   )     ! physics forcing of u
    deallocate( v_tdcy_phy   )     ! physics forcing of v
    deallocate( trc_tdcy_phy )     ! physics forcing of tracers
    deallocate( gfs_u_old )
    deallocate( gfs_v_old )
    deallocate( gfs_t_old )
    deallocate( gfs_q_old )
    deallocate( gfs_oz_old )
    deallocate( gfs_cld_old )
    
    return
  end subroutine cpl_finalize
end module module_fim_cpl_finalize
