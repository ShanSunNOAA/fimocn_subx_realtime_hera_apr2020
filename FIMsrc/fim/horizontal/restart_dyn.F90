module restart_dyn
    use module_sfc_variables, only: rn2d, rc2d, u10m, v10m, sst_prev, sst_next, fice2d_prev,	&
      fice2d_next, fice2d, runoff
    use module_variables,     only: curr_write_time, nf, of, vof, adbash1, adbash2, adbash3, psrf,	&
      ptdcy, pw2d, pq2d, u_tdcy, v_tdcy, dp_tdcy, dpl_tdcy, tr3d, trdp, 				&
      trc_tdcy, trl_tdcy, us3d, vs3d, ws3d, mp3d, tk3d, dp3d, rh3d, 					&
      relvor, dpinit, pr3d, ex3d, ph3d, sdot, massfx, cumufx,wt_int_rev,k_int_rev
    use module_constants,     only: dpsig, thetac, lat, lon, nprox, proxs, area, cs, sn, sidevec_c,	&
      sidevec_e, sideln, rsideln, rprox_ln, area, rarea, corio,						&
      deg_lat, deg_lon, sigak, sigbk
    USE module_control,       ONLY: dt, nabl, ntra, ntrb, npp, nvlp1, next_date, have_next_sst
    USE fimnamelist,          ONLY:  nvl, pure_sig
    use module_globsum,       only: qmstr, qmstrc, qmstrn, qdtr_set
    use stencilprint

contains

  subroutine read_restart_dyn (unitno)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read required dynamics fields from the restart file. SMS will modify the code to do the
! appropriate single-task reading of the restart file, and subsequent scattering of the
! data to other MPI tasks.
!
! CRITICAL: If you modify this file, you MUST also modify write_restart_dyn.F90 in the
! same way. Otherwise restart will be broken.
!
! readarr32_[i|r] assumes 32-bit data.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    integer, intent(in) :: unitno ! unit number to read from

    integer :: n, t ! indices

    read (unitno, err=90) curr_write_time, nf, of, vof, adbash1, adbash2, adbash3, thetac, dpsig, sigak, sigbk
    read (unitno, err=90) lat, lon, nprox, proxs, area, cs, sn, psrf, ptdcy
    read (unitno, err=90) qmstr, qmstrc, qmstrn, qdtr_set
    read (unitno, err=90) sidevec_c, sidevec_e, sideln, rprox_ln
    read (unitno, err=90) sst_prev, sst_next, fice2d_prev, fice2d_next, fice2d, next_date, have_next_sst
    read (unitno, err=90) rarea, rsideln, corio, deg_lat, deg_lon, rn2d, pw2d, rc2d, u10m, v10m, runoff

! These arrays are dimensioned (nvl,nip[,other dimensions]):
    do n=1,nabl
      call readarr32_r (u_tdcy(:,:,n), nvl, unitno)
      call readarr32_r (v_tdcy(:,:,n), nvl, unitno)
      call readarr32_r (dp_tdcy(:,:,n), nvl, unitno)
      call readarr32_r (dpl_tdcy(:,:,n), nvl, unitno)
    end do

    do t=1,ntra+ntrb
      call readarr32_r (tr3d(:,:,t), nvl, unitno)
      call readarr32_r (trdp(:,:,t), nvl, unitno)
      do n=1,nabl
        call readarr32_r (trc_tdcy(:,:,n,t), nvl, unitno)
        call readarr32_r (trl_tdcy(:,:,n,t), nvl, unitno)
      end do
    end do

    call readarr32_r (us3d, nvl, unitno)
    call readarr32_r (vs3d, nvl, unitno)
    call readarr32_r (ws3d, nvl, unitno)
    call readarr32_r (mp3d, nvl, unitno)
    call readarr32_r (tk3d, nvl, unitno)
    call readarr32_r (dp3d, nvl, unitno)
    call readarr32_r (rh3d, nvl, unitno)
    call readarr32_r (relvor,nvl,unitno)
    call readarr32_r (dpinit, nvl, unitno)

! These arrays are dimensioned (nvlp1,nip):
    call readarr32_r (pr3d, nvlp1, unitno)
    call readarr32_r (ex3d, nvlp1, unitno)
    call readarr32_r (ph3d, nvlp1, unitno)
    call readarr32_r (sdot, nvlp1, unitno)

! These arrays are dimensioned (nvl,npp,nip[,other dimensions]):
! Simplest coding folds nvl*npp into a single dimension. Will need to rewrite massfx
! and cumufx to a transpose if root process can't hold npp 3-d fields in memory
    do n=1,nabl
      call readarr32_r (massfx(:,:,:,n), nvl*npp, unitno)
    end do
    call readarr32_r (cumufx(:,:,:), nvl*npp, unitno)

    write(6,*) 'read_restart_dyn: successfully read dynamics fields from restart file'

    return

90  write(6,*)'read_restart_dyn: Error reading from unit ', unitno, '. Stopping'
    stop

  end subroutine read_restart_dyn

  subroutine write_restart_dyn (unitno)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write required dynamics fields to the restart file. SMS will modify the code to do the
! appropriate single-task writing of the restart file, after gathering the
! data from other MPI tasks.
!
! CRITICAL: If you modify this file, you MUST also modify read_restart_dyn.F90 in the
! same way. Otherwise restart will be broken.
!
! writearr32_[i|r] assumes 32-bit data.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    implicit none

    integer, intent(in) :: unitno ! unit number to write to

    integer :: n, t ! indices

    write (unitno, err=90) curr_write_time, nf, of, vof, adbash1, adbash2, adbash3, thetac, dpsig, sigak, sigbk
    write (unitno, err=90) lat, lon, nprox, proxs, area, cs, sn, psrf, ptdcy
    write (unitno, err=90) qmstr, qmstrc, qmstrn, qdtr_set
    write (unitno, err=90) sidevec_c, sidevec_e, sideln, rprox_ln
    write (unitno, err=90) sst_prev, sst_next, fice2d_prev, fice2d_next, fice2d, next_date, have_next_sst
    write (unitno, err=90) rarea, rsideln, corio, deg_lat, deg_lon, rn2d, pw2d, rc2d, u10m, v10m, runoff

! These arrays are dimensioned (nvl,nip[,other dimensions]):
    do n=1,nabl
      call writearr32_r (u_tdcy(:,:,n), nvl, unitno)
      call writearr32_r (v_tdcy(:,:,n), nvl, unitno)
      call writearr32_r (dp_tdcy(:,:,n), nvl, unitno)
      call writearr32_r (dpl_tdcy(:,:,n), nvl, unitno)
    end do

    do t=1,ntra+ntrb
      call writearr32_r (tr3d(:,:,t), nvl, unitno)
      call writearr32_r (trdp(:,:,t), nvl, unitno)
      do n=1,nabl
        call writearr32_r (trc_tdcy(:,:,n,t), nvl, unitno)
        call writearr32_r (trl_tdcy(:,:,n,t), nvl, unitno)
      end do
    end do

    call writearr32_r (us3d, nvl, unitno)
    call writearr32_r (vs3d, nvl, unitno)
    call writearr32_r (ws3d, nvl, unitno)
    call writearr32_r (mp3d, nvl, unitno)
    call writearr32_r (tk3d, nvl, unitno)
    call writearr32_r (dp3d, nvl, unitno)
    call writearr32_r (rh3d, nvl, unitno)
    call writearr32_r (relvor,nvl, unitno)
    call writearr32_r (dpinit, nvl, unitno)

! These arrays are dimensioned (nvlp1,nip):
    call writearr32_r (pr3d, nvlp1, unitno)
    call writearr32_r (ex3d, nvlp1, unitno)
    call writearr32_r (ph3d, nvlp1, unitno)
    call writearr32_r (sdot, nvlp1, unitno)

! These arrays are dimensioned (nvl,npp,nip[,other dimensions]):
! Simplest coding folds nvl*npp into a single dimension. Will need to rewrite massfx
! and cumufx to a transpose if root process can't hold npp 3-d fields in memory
    do n=1,nabl
      call writearr32_r (massfx(:,:,:,n), nvl*npp, unitno)
    end do
    call writearr32_r (cumufx(:,:,:), nvl*npp, unitno)

    write(6,*) 'write_restart_dyn: successfully wrote dynamics fields to restart file'

    return

90  write(6,*)'write_restart_dyn: Error writing to unit ', unitno, '. Stopping'
    stop

  end subroutine write_restart_dyn

end module restart_dyn
