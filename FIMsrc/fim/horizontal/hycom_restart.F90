module hycom_restart
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read required ocean fields from the restart file. SMS will modify the code
! to do the appropriate single-task reading of the restart file, and subsequent
! scattering of the data to other MPI tasks.
!
! CRITICAL: If you modify this file, you MUST modify write_restart_ocn.F90 in
! the same way. Otherwise restart will be broken.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  use module_variables,only: curr_write_time
  use module_control,  only: npp
  use hycom_variables, only: temp,saln,passv_tr,dp,uclin,vclin,		&
        geop,montg,pbot,gzbot,psikk,spvkk,temice,thkice,covice,		&
        dens,ptrop,pres,dpinit,cumuflx,ocnarea,rivflo,		&
        srfx_ave,pmne_ave,qf2d_ave,hf2d_ave,ssht_ave,			&
        hmixl,saltini,heatini,rhoglb0,temglb0,salglb0,massglb0,		&
        totmass,pcpcor,radcor,masscor,srfxcum,pmnecum
        
  use fimnamelist     ,only: kdm
  use hycom_control   ,only: numtr
  use hycom_constants ,only: nuphill,uphill,dnhill,odepth,wet
  implicit none
contains
!
subroutine read_restart_ocn (unitno)

  integer, intent(in) :: unitno ! unit number to read from

  integer :: n

  read (unitno, err=90) curr_write_time,pbot,gzbot,psikk,spvkk,		&
        temice,covice,thkice,ocnarea,rivflo,nuphill,uphill,dnhill,	&
        srfx_ave,pmne_ave,qf2d_ave,hf2d_ave,ssht_ave,odepth,wet,hmixl,	&
        saltini,heatini,rhoglb0,temglb0,salglb0,massglb0,totmass,	&
        pcpcor,radcor,masscor,srfxcum,pmnecum

!SMS$PARALLEL (dh,i) BEGIN
! exchange needed for bitwise-identical
!SMS$EXCHANGE (odepth,wet)
!SMS$PARALLEL END

  do n=1,numtr
    call readarr32_r (passv_tr(:,:,n), kdm, unitno)
  end do

  do n=1,2
    call readarr32_r (dp      (:,:,n), kdm, unitno)
    call readarr32_r (uclin   (:,:,n), kdm, unitno)
    call readarr32_r (vclin   (:,:,n), kdm, unitno)
  end do

  call readarr32_r (temp,   kdm,   unitno)
  call readarr32_r (saln,   kdm,   unitno)
  call readarr32_r (montg,  kdm,   unitno)
  call readarr32_r (dens,   kdm,   unitno)
  call readarr32_r (dpinit, kdm,   unitno)
  call readarr32_r (pres,   kdm+1, unitno)
  call readarr32_r (ptrop,  3,     unitno)
  call readarr32_r (cumuflx(:,:,:), kdm*npp, unitno)

  write(6,*) 'read_restart_ocn: successfully read ocean fields from restart file'
  return

90 write(6,*) 'read_restart_ocn: Error reading to unit ', unitno, '. Stopping'
  stop
end subroutine read_restart_ocn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Write required ocean fields to the restart file. SMS will modify the code to do the
!! appropriate single-task writing of the restart file, after gathering the
!! data from other MPI tasks.
!
!! CRITICAL: If you modify this file, you MUST also modify read_restart_ocn.F90 in the
!! same way. Otherwise restart will be broken.
!
!! read_restart_ocn and write_restart_ocn belong in a module, but SMS doesn't like multiple
!! subroutines in a file.
!
!! writearr32_r assumes 32-bit data.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_restart_ocn (unitno)

  integer, intent(in) :: unitno ! unit number to write to

  integer :: n

  write (unitno, err=90) curr_write_time,pbot,gzbot,psikk,spvkk,	&
        temice,covice,thkice,ocnarea,rivflo,nuphill,uphill,dnhill,	&
        srfx_ave,pmne_ave,qf2d_ave,hf2d_ave,ssht_ave,odepth,wet,hmixl,	&
        saltini,heatini,rhoglb0,temglb0,salglb0,massglb0,totmass,	&
        pcpcor,radcor,masscor,srfxcum,pmnecum

  do n=1,numtr
    call writearr32_r (passv_tr(:,:,n), kdm, unitno)
  end do

  do n=1,2
    call writearr32_r (   dp   (:,:,n), kdm, unitno)
    call writearr32_r (uclin   (:,:,n), kdm, unitno)
    call writearr32_r (vclin   (:,:,n), kdm, unitno)
  end do

  call writearr32_r (temp,   kdm,   unitno)
  call writearr32_r (saln,   kdm,   unitno)
  call writearr32_r (montg,  kdm,   unitno)	! for ssh
  call writearr32_r (dens,   kdm,   unitno)
  call writearr32_r (dpinit, kdm,   unitno)
  call writearr32_r (pres,   kdm+1, unitno)
  call writearr32_r (ptrop,  3,     unitno)	! for ssh
  call writearr32_r (cumuflx(:,:,:), kdm*npp, unitno)

  write(6,*) 'write_restart_ocn: successfully write ocean fields to restart file'

  return

90 write(6,*)'write_restart_ocn: Error writing to unit ', unitno, '. Stopping'
  stop
end subroutine write_restart_ocn
end module hycom_restart
