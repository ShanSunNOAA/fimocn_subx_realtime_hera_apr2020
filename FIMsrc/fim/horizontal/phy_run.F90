module module_fim_phy_run

  implicit none

#include <gptl.inc>

contains

  subroutine phy_run(its)
!*********************************************************************
!       "Run" method for fim global model physics
!       Alexander E. MacDonald  12/24/2005
!       J. LEE                  12/28/2005
!*********************************************************************

    use module_constants
    use module_control, only: CallPhysics, CallRadiation, itsDFI
    use fimnamelist,    only: nts, UpdateSST, itsStart, digifilt
    use module_physics, only: physics

!  Declare dummy arguments
    integer, intent(in) :: its

!  Declare local variables:
    integer :: ret
    integer :: laststep

    ret = gptlstart ('phy_run')

!...........................................................
! Advance the physics component by one time step unless this 
! is the last (nts+1) iteration.  
! This complexity is required for the NCEP ESMF approach 
! in which single-phase DYN and PHY components alternate 
! execution during each time step.  
!
!TBH:  Restore if statement label once Mark fixes PPP
!TBH skip_last_iteration: if (its <= nts ) then

    if (digifilt) then
      laststep=itsDFI+nts
    else
      laststep=itsStart+nts
    end if
    if (its < laststep) then

!...........................................................
! Do column calculations:
!...........................................................
! Condensational heating parameterizations
!
      if (UpdateSST) then
        call sst_run(its)
      endif
      call physics (its,          &
        CallPhysics,CallRadiation ) ! Timestep interval to call physics,radiation
    endif
!TBH:  Restore if statement label once Mark fixes PPP
!TBH endif skip_last_iteration

    ret = gptlstop ('phy_run')

    return
  end subroutine phy_run

  subroutine sst_run(its)
!*********************************************************************
!       "Run" method for fim global model physics
!       Alexander E. MacDonald  12/24/2005
!       J. LEE                  12/28/2005
!*********************************************************************

    use module_constants
    use module_control,  only: CallSST
    use fimnamelist,     only: nts, itsStart

    implicit none

!  Declare dummy arguments
    integer, intent(in) :: its
    integer             :: mype
!...........................................................
! Advance the physics component by one time step unless this 
! is the last (nts+1) iteration.  
! This complexity is required for the NCEP ESMF approach 
! in which single-phase DYN and PHY components alternate 
! execution during each time step.  
!
!TBH:  Restore if statement label once Mark fixes PPP
!TBH skip_last_iteration: if (its <= nts ) then
!sms$comm_rank(mype)
    if (its < itsStart+nts ) then

!...........................................................
! Do sst update:
!...........................................................
! Condensational heating parameterizations
!
      if (mod(its,CallSST) == 0) call update_sst (its) ! Timestep interval to call update sst (set at once 1 day)

    end if
!TBH:  Restore if statement label once Mark fixes PPP
!TBH endif skip_last_iteration
    return
  end subroutine sst_run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine update_sst (its)
    use module_control,       only: prev_date,next_date,numphr,nip,have_next_sst
    use fimnamelist,          only: yyyymmddhhmm
    use module_sfc_variables, only: ts2d,sst_prev,sst_next,fice2d,fice2d_prev,fice2d_next,slmsk2d,hice2d
    use gfs_physics_internal_state_mod, only: gfs_physics_internal_state, gis_phy
    use module_fim_phy_init, only: sstunit
    USE slint, ONLY: bl_int
    use infnan, only: inf, negint

    implicit none

    integer, intent(in) :: its
! Local variables
    integer,parameter :: im_oc=360,jm_oc=180
    integer :: mype,ipn
    integer :: YEARi,YEARm
    integer :: MONTHi,MONTHm
    integer :: DAYi,DAYm,HH,MM,SS
    integer :: MIDMON

    integer :: ID
    integer :: M1
    integer :: M2
    integer :: MIDM
    integer :: MIDP

    integer :: idat(8),mdat(8),DAYS(12),time_header(4)
    integer :: ierr
    real*8  :: rinc(5)
    real    :: FAC
    real, allocatable :: sst_next_ll(:,:)
    real, allocatable :: fice2d_next_ll(:,:)


    data    DAYS /31,28,31,30,31,30,31,31,30,31,30,31/
!sms$comm_rank(mype)
    READ(yyyymmddhhmm(1:4),   FMT='(I4)') YEARi
    READ(yyyymmddhhmm(5:6),   FMT='(I2)') MONTHi
    READ(yyyymmddhhmm(7:8),   FMT='(I2)') DAYi
    idat = 0
    idat(1) = YEARi
    idat(2) = MONTHi
    idat(3) = DAYi

    rinc(:) = 0
    rinc(2) = (float(its))/float(numphr)
    print*,'In update_sst yeari,monthi,dayi, fhour ',YEARi,MONTHi,DAYi,rinc(2)
    call w3movdat (rinc, idat, mdat)

    YEARm =mdat(1)
    MONTHm=mdat(2)
    DAYm  =mdat(3)
    HH   = mdat(5)
    MM   = mdat(6)
    SS   = mdat(7)

    IF (mod(YEARm,4) == 0) THEN
      DAYS(2) = 29.0
    ENDIF
    IF (YEARm == 1900) THEN
      DAYS(2)=28.0
    ENDIF

    MIDMON = DAYS(MONTHm)/2 + 1

    IF (DAYm < MIDMON) THEN

      M1   = MOD(MONTHm+10,12) + 1
      M2   = MONTHm
      MIDM = DAYS(M1)/2 + 1
      MIDP = DAYS(M1) + MIDMON
      ID = DAYm + DAYS(M1)
      have_next_sst=.false.

    ELSE

      M2   = MOD(MONTHm,12) + 1
      M1   = MONTHm
      MIDM = MIDMON
      MIDP = DAYS(M2)/2 + 1 + DAYS(M1)
      ID = DAYm

    ENDIF

    if (DAYm == MIDMON .AND. .not. have_next_sst) then    
      sst_prev      = sst_next
      fice2d_prev   = fice2d_next
      prev_date     = next_date
      have_next_sst = .true.
!SMS$SERIAL (<sst_next,fice2d_next,OUT> : default=ignore) BEGIN 
      allocate(sst_next_ll(im_oc,jm_oc))
      allocate(fice2d_next_ll(im_oc,jm_oc))

! find matching dates to read SST
      next_date(:)=negint
      do while (M1 /= prev_date(2) .or. M2 /= next_date(2))
        read(sstunit,iostat=ierr) time_header
        if (ierr.ne.0) then
          write (*,'(a)') 'update_sst: error reading time_header'
          stop
        endif
        read(sstunit,iostat=ierr) sst_next_ll
        if (ierr.ne.0) then
          write (*,'(a)') 'update_sst: error reading sst_next_ll'
          stop
        endif
        read(sstunit,iostat=ierr) fice2d_next_ll
        if (ierr.ne.0) then
          write (*,'(a)') 'update_sst: error reading fice2d_next_ll'
          stop
        endif
        next_date=time_header
        print*,'New sst records read',time_header,M1,prev_date(2),M2,next_date(2)
      end do

      CALL bl_int (sst_next_ll(:,:), sst_next)
      CALL bl_int (fice2d_next_ll(:,:), fice2d_next)
      deallocate(sst_next_ll)
      deallocate(fice2d_next_ll)
!SMS$SERIAL END
    end if	!if (DAYm == MIDMON .AND. .not. have_next_sst)
!SMS$SERIAL (<FAC,out> : default=ignore) BEGIN
    FAC = (real(ID -   MIDM)*86400 + HH*3600 + MM*60 + SS)/ &
      (real(MIDP - MIDM)*86400         )
    print*,'In update_sst year,month,day ',YEARm,MONTHm,DAYm,FAC
!SMS$SERIAL END
! replace ts2d over ocean points
!SMS$PARALLEL(dh, ipn) BEGIN
    DO ipn=1,nip
!  need logic to keep ice's temperature the same, only update sst and ice fraction
! update ice fraction 1st
! if there is new ice, set it to -1.8 and hice=0.0
! if ice melts, then set ts2d to sst and hice=0.0
      IF (slmsk2d(ipn) /= 1) THEN
        fice2d(ipn) = fice2d_next(ipn)*(FAC) + fice2d_prev(ipn)*(1.0 - FAC)
        if (fice2d(ipn) > 1) fice2d(ipn)=1.0
        if (fice2d(ipn) < 0) fice2d(ipn)=0.0
      ENDIF
      IF (fice2d(ipn) >= 0.5 .AND. slmsk2d(ipn) == 0) THEN ! freeze open ocean
        slmsk2d(ipn)=2.0
        ts2d(ipn)=271.35
        hice2d(ipn)=0.0
      ENDIF
      IF (fice2d(ipn) < 0.5 .AND. slmsk2d(ipn) == 2) THEN ! melt sea-ice
        slmsk2d(ipn)=0.0
        hice2d(ipn)=0.0
      ENDIF
      IF (slmsk2d(ipn) == 0) THEN
        ts2d(ipn)=sst_next(ipn)*(FAC)+sst_prev(ipn)*(1.0 - FAC)
      ENDIF

      gis_phy%sfc_fld%TSEA  (ipn,1) = ts2d    (ipn)
      gis_phy%sfc_fld%HICE  (ipn,1) = hice2d  (ipn)
      gis_phy%sfc_fld%FICE  (ipn,1) = fice2d  (ipn)
      gis_phy%sfc_fld%SLMSK (ipn,1) = slmsk2d (ipn)
    ENDDO
!SMS$PARALLEL END
    RETURN

  end subroutine update_sst

end module module_fim_phy_run
