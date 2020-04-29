SUBROUTINE add_incr ()

! read increment produced by EnKF  on  gaussian grid
! interpolate to icos, use weights to distrubute vertically and add
! to background

  USE module_control,ONLY: nvlp1,nip,ntra,ntrb
  USE fimnamelist   ,ONLY: glvl,nvl,pure_sig,incr_fname
  USE module_constants,ONLY: lat,lon,sigak,sigbk,thetac,p1000,cp,rd
  USE module_variables, ONLY: us3d,vs3d,pr3d,dp3d,ex3d,tr3d,sdot,ex3d,&
    wt_int_rev,k_int_rev
  USE kinds, ONLY : r_single
  USE units, ONLY: getunit, returnunit
  USE slint, ONLY: slint_init,bilinear_interp,bilinear_interp_uv

  IMPLICIT NONE

!locals

  REAL(r_single), ALLOCATABLE, DIMENSION(:,:) :: latlongauss
  REAL(r_single), ALLOCATABLE, DIMENSION(:,:) :: vargauss
  REAL(r_single), ALLOCATABLE, DIMENSION(:) :: tmp

  INTEGER :: nlons,nlats,nlevs,npts
  REAL :: clip

  REAL, DIMENSION(nvl,nip) :: tmp3d !for temporarry storage of thv increment
  REAL, DIMENSION(nip,2) :: latlonicos
  REAL, DIMENSION(nvl,nip) :: ug,vg
  INTEGER :: unitno,i,j,k,l,ivar,ipn,status,nvars
  logical::error

  clip=TINY(1.)

  unitno=getunit()

  status=0
!sms$serial (default=ignore) begin
  OPEN (unit=unitno,file=TRIM(incr_fname),form="unformatted",iostat=status,status='old')
!sms$serial end
  IF ( status /= 0 ) THEN
    PRINT *,'Failed to open file: ',TRIM(incr_fname)
    STOP
  ENDIF

! Set small defaults for non-root tasks.

  nlons=1
  nlats=1
  nlevs=1
  nvars=1

  status=0
!sms$serial (default=ignore) begin
  READ (unitno,iostat=status) nlons,nlats,nlevs,nvars
!sms$serial end
  IF ( status /= 0 ) THEN
    PRINT *,'Failed to read nlons,nlats,nlevs in file: ',TRIM(incr_fname)
    STOP
  ENDIF

  npts = nlons*nlats

  ALLOCATE(latlongauss(npts,2),vargauss(nlevs,npts))
  ALLOCATE(tmp(npts))

!longitude first

  status=0
!sms$serial (default=ignore) begin
  READ (unitno,iostat=status) tmp
!sms$serial end
  IF ( status /= 0 ) THEN
    PRINT *,'Failed to read lon in file: ',TRIM(incr_fname)
    STOP
  ENDIF

!sms$serial (default=ignore) begin
  latlongauss(:,2)=tmp
!sms$serial end

!latitude second

  status=0
!sms$serial (default=ignore) begin
  READ (unitno,iostat=status) tmp
!sms$serial end
  IF ( status /= 0 ) THEN
    PRINT *,'Failed to read lat in file: ',TRIM(incr_fname)
    STOP
  ENDIF

!sms$serial (<lat,lon,in>:default=ignore) begin
  latlongauss(:,1)=tmp
  latlonicos(:,1)=lat
  latlonicos(:,2)=lon
  CALL slint_init(latlongauss,npts,latlonicos,nip)
!sms$serial end

!for all vars all levels of vargauss need to be available for interp

!u v

!u

  error=.false.
!sms$serial (default=ignore) begin
  do k=1,nlevs
    if (.not.error) then
      read (unitno,iostat=status) tmp
      if (status.eq.0) then
        vargauss(k,:)=tmp
        call bilinear_interp(vargauss(k,:),ug(k,:))
      else
        error=.true.
        exit
      endif
    endif
  enddo
!sms$serial end
  IF (error) THEN
    PRINT *,'Failed to read u in file: ',TRIM(incr_fname)
    STOP
  ENDIF

!sms$serial (<k_int_rev,wt_int_rev,in>,<us3d,inout>:default=ignore) begin
  DO k=1,nlevs
    IF (pure_sig) THEN
      us3d(k,:)=us3d(k,:)+ug(k,:)
    ELSE
      us3d(k,:) = us3d(k,:) + &
        (1.-wt_int_rev(k,k_int_rev(k,:)))*ug(k,k_int_rev(k,:))+&
        wt_int_rev(k,k_int_rev(k,:))*ug(k,k_int_rev(k,:)+1)
    ENDIF
  ENDDO
!sms$serial end

!v

  error=.false.
!sms$serial (default=ignore) begin
  do k=1,nlevs
    if (.not.error) then
      read (unitno,iostat=status) tmp
      if (status.eq.0) then
        vargauss(k,:)=tmp
        call bilinear_interp(vargauss(k,:),ug(k,:))
      else
        error=.true.
        exit
      endif
    endif
  enddo
!sms$serial end
  IF (error) THEN
    PRINT *,'Failed to read v in file: ',TRIM(incr_fname)
    STOP
  ENDIF

!sms$serial (<k_int_rev,wt_int_rev,in>,<vs3d,inout>:default=ignore) begin
  DO k=1,nlevs
    IF ( pure_sig) THEN
      vs3d(k,:)=vs3d(k,:)+ug(k,:)
    ELSE
      vs3d(k,:) = vs3d(k,:) + &
        (1.-wt_int_rev(k,k_int_rev(k,:)))*ug(k,k_int_rev(k,:))+&
        wt_int_rev(k,k_int_rev(k,:))*ug(k,k_int_rev(k,:)+1)
    ENDIF
  ENDDO
!sms$serial end

!thv

  error=.false.
!sms$serial (default=ignore) begin
  do k=1,nlevs
    if (.not.error) then
      read (unitno,iostat=status) tmp
      if (status.eq.0) then
        vargauss(k,:)=tmp
        call bilinear_interp(vargauss(k,:),ug(k,:))
      else
        error=.true.
      endif
    endif
  enddo
!sms$serial end
  IF (error) THEN
    PRINT *,'Failed to read thv in file: ',TRIM(incr_fname)
    STOP
  ENDIF

!save incr for virt potential

!sms$serial (<k_int_rev,wt_int_rev,in>:default=ignore) begin
  DO k=1,nlevs
    IF ( pure_sig) THEN
      tmp3d(k,:)=ug(k,:)
    ELSE
      tmp3d(k,:)=&
        (1.-wt_int_rev(k,k_int_rev(k,:)))*ug(k,k_int_rev(k,:))+&
        wt_int_rev(k,k_int_rev(k,:))*ug(k,k_int_rev(k,:)+1)
    ENDIF
  ENDDO
!sms$serial end

!qv

  error=.false.
!sms$serial (default=ignore) begin
  do k=1,nlevs
    if (.not.error) then
      read (unitno,iostat=status) tmp
      if (status.eq.0) then
        vargauss(k,:)=tmp
        call bilinear_interp(vargauss(k,:),ug(k,:))
      else
        error=.true.
      endif
    endif
  enddo
!sms$serial end
  IF (error) THEN
    PRINT *,'Failed to read qv in file: ',TRIM(incr_fname)
    STOP
  ENDIF

!sms$serial (<k_int_rev,wt_int_rev,in>,<tr3d,inout>:default=ignore) begin
  DO k=1,nlevs
    IF ( pure_sig) THEN
!recalculate pottemp from virtual
      tr3d(k,:,1)=tr3d(k,:,1)+&
        (tmp3d(k,:)-0.6078*tr3d(k,:,1)*ug(k,:))/(1.+0.6078*tr3d(k,:,2))
      tr3d(k,:,2)=tr3d(k,:,2)+ug(k,:)
    ELSE
      tr3d(k,:,1)=tr3d(k,:,1)+&
        (tmp3d(k,:)-0.6078*tr3d(k,:,1)*ug(k,:))/(1.+0.6078*tr3d(k,:,2))
      tr3d(k,:,2) = tr3d(k,:,2) + &
        (1.-wt_int_rev(k,k_int_rev(k,:)))*ug(k,k_int_rev(k,:))+&
        wt_int_rev(k,k_int_rev(k,:))*ug(k,k_int_rev(k,:)+1)
    ENDIF
  ENDDO
!sms$serial end

!there is no increment for clw ivar=3

!ozone

  ivar=4

  error=.false.
!sms$serial (default=ignore) begin
  do k=1,nlevs
    if (.not.error) then
      read (unitno,iostat=status) tmp
      if (status.eq.0) then
        vargauss(k,:)=tmp
        call bilinear_interp(vargauss(k,:),ug(k,:))
      else
        error=.true.
      endif
    endif
  enddo
!sms$serial end
  IF (error) THEN
    PRINT *,'Failed to read ozone in file: ',TRIM(incr_fname)
    STOP
  ENDIF

!sms$serial (<k_int_rev,wt_int_rev,in>,<tr3d,inout>:default=ignore) begin
  DO k=1,nlevs
    IF ( pure_sig) THEN
      tr3d(k,:,ivar)=tr3d(k,:,ivar)+ug(k,:)
    ELSE
      tr3d(k,:,ivar) = tr3d(k,:,ivar) + &
        (1.-wt_int_rev(k,k_int_rev(k,:)))*ug(k,k_int_rev(k,:))+&
        wt_int_rev(k,k_int_rev(k,:))*ug(k,k_int_rev(k,:)+1)
    ENDIF
  ENDDO

  DO i=1,nip
    DO ivar=1,ntra+ntrb
      DO k=1,nlevs
        IF ( tr3d(k,i,ivar) < clip ) tr3d(k,i,ivar)=clip
      ENDDO
    ENDDO
  ENDDO
!sms$serial end

!  DO i=1,nip
!     WHERE( tr3d(:,i,:) < clip ) tr3d(:,i,:)=clip
!  ENDDO

!  WHERE (tr3d < clip) tr3d=clip

!psfc

  status=0
!sms$serial (default=ignore) begin
  read (unitno,iostat=status) tmp
!sms$serial end
  IF ( status /= 0 ) THEN
    PRINT *,'Failed to read psfc in file: ',TRIM(incr_fname)
    STOP
  ENDIF

  k=1

!sms$serial (<dp3d,ex3d,pr3d,inout>:default=ignore) begin
  vargauss(k,:)=tmp
  CALL bilinear_interp(vargauss(k,:),ug(k,:))

  pr3d(k,:)=pr3d(k,:)+ug(k,:)

  DO k=2,nlevs+1
    pr3d(k,:)=sigak(k)+sigbk(k)*pr3d(1,:)
  ENDDO

  IF ( .NOT. pure_sig ) THEN
    DO k=2,nlevs
      pr3d(k,:)=pr3d(k,:)+ug(1,:)*&
        (pr3d(k,:)-pr3d(nlevs+1,:))/(pr3d(1,:)-pr3d(nlevs+1,:))
    ENDDO
  ENDIF

  DO k=1,nlevs
    dp3d(k,:)=pr3d(k,:)-pr3d(k+1,:)
    ex3d(k,:)=cp*(pr3d(k,:)/p1000)**(rd/cp)
  ENDDO

  ex3d(nlevs+1,:)=cp*(pr3d(nlevs+1,:)/p1000)**(rd/cp)
!sms$serial end

  CALL returnunit(unitno)

  DEALLOCATE(vargauss,latlongauss,tmp)

END SUBROUTINE add_incr
