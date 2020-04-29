module restart_cpl

contains

  subroutine read_restart_cpl (unitno)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read_restart_cpl: read physics fields from the restart file
! SMS doesn't yet properly handle Fortran derived types, so those fields
! need to go through an interface routine (readarr64[t]_r).
!
! !!!!!CRITICAL!!!!! Any changes to fields read in here MUST be made in exactly the same way
! in write_restart_cpl.F90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use module_control,                 only: ntra,ntrb
    use fimnamelist,                    only: nvl
!tgs
    use module_fim_cpl_run,        only: u_tdcy_phy,v_tdcy_phy,trc_tdcy_phy, &
      gfs_u_old, gfs_v_old, gfs_t_old,    &
      gfs_q_old, gfs_oz_old, gfs_cld_old

    implicit none

! Input arguments
    integer, intent(in) :: unitno   ! Unit number to read from

! Local workspace
    integer :: n      ! Loop indices

! cpl state
! These arrays are dimensioned (nip,nvl):
    call readarr64_r (gfs_u_old,    nvl, unitno)
    call readarr64_r (gfs_v_old,    nvl, unitno)
    call readarr64_r (gfs_t_old,    nvl, unitno)
    call readarr64_r (gfs_q_old,    nvl, unitno)
    call readarr64_r (gfs_oz_old,   nvl, unitno)
    call readarr64_r (gfs_cld_old,  nvl, unitno)
!tgs read in physics tendencies
! These arrays are dimensioned (nvl,nip,[other dimensions]):
    call readarr64t_r (u_tdcy_phy, nvl, unitno)
    call readarr64t_r (v_tdcy_phy, nvl, unitno)
    do n=1,ntra+ntrb
      call readarr64t_r (trc_tdcy_phy(:,:,n), nvl, unitno)
    end do

    write (6,*) 'read_restart_cpl: successfully read cpl fields from restart file'

    return

  end subroutine read_restart_cpl

  subroutine write_restart_cpl (unitno)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write_restart_cpl: write physics fields to the restart file
! SMS doesn't yet properly handle Fortran derived types, so those fields
! need to go through an interface routine (writearr64[t]_r).
!
! !!!!!CRITICAL!!!!! Any changes to fields read in here MUST be made in exactly the same way
! in read_restart_cpl.F90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use module_control,                 only: ntra,ntrb
    use fimnamelist,                    only: nvl
!tgs
    use module_fim_cpl_run,        only: u_tdcy_phy,v_tdcy_phy,trc_tdcy_phy, &
      gfs_u_old, gfs_v_old, gfs_t_old,    &
      gfs_q_old, gfs_oz_old, gfs_cld_old

    implicit none

! Input arguments
    integer, intent(in) :: unitno   ! Unit number to write to

! Local workspace
    integer :: n      ! Loop indices

! cpl state
! These arrays are dimensioned (nip,nvl):
    call writearr64_r (gfs_u_old,    nvl, unitno)
    call writearr64_r (gfs_v_old,    nvl, unitno)
    call writearr64_r (gfs_t_old,    nvl, unitno)
    call writearr64_r (gfs_q_old,    nvl, unitno)
    call writearr64_r (gfs_oz_old,   nvl, unitno)
    call writearr64_r (gfs_cld_old,  nvl, unitno)
!tgs Write out physics tendencies
! These arrays are dimensioned (nvl,nip,[other dimensions]):
    call writearr64t_r (u_tdcy_phy, nvl, unitno)
    call writearr64t_r (v_tdcy_phy, nvl, unitno)
    do n=1,ntra+ntrb
      call writearr64t_r (trc_tdcy_phy(:,:,n), nvl, unitno)
    end do

    write (6,*) 'write_restart_cpl: successfully wrote cpl fields to restart file'

    return

  end subroutine write_restart_cpl

end module restart_cpl
