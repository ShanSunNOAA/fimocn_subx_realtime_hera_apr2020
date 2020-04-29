MODULE module_outvar_enkf_spectral

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: outvar_enkf_spectral,cycle_hour,&
    jcap,nlonsgauss,nlatsgauss,idate,fhour,timestr

  INTEGER, PARAMETER :: cycle_hour=6
  INTEGER , DIMENSION(4) :: idate
  CHARACTER(len=6) :: timestr
  REAL :: fhour
  INTEGER :: jcap,nlonsgauss,nlatsgauss

CONTAINS

  SUBROUTINE glvl2gauss(glvl,jcap,lonb,latb,nlons_lin,nlats_lin,nwaves)

    INTEGER, INTENT(in) :: glvl
    INTEGER, INTENT(out) :: jcap,lonb,latb,nlons_lin,nlats_lin,nwaves

!??    latf=(jcap+1)*3/2
!??    lonf=(jcap+1)*3

    IF (glvl <= 6) THEN

      jcap=126
      lonb=384
      latb=190

    ELSE IF (glvl == 7) THEN

      jcap=254
      lonb=768
      latb=384

    ELSE IF (glvl == 8) THEN

      jcap=382
      lonb=1152
      latb=576

    ELSE
      WRITE(6,*)'In module_outvar_enkf.F90 min(glvl)=6 max(glvl)=8'
      WRITE(6,*)'glvl = ',glvl,' outside of range'
      WRITE(6,*)'Stopping'
      CALL FLUSH(6)
      STOP
    END IF

    nlats_lin=jcap+2
    nlons_lin=2*nlats_lin
    nwaves=(jcap+1)*(jcap+2)

  END SUBROUTINE glvl2gauss

  SUBROUTINE outvar_enkf_spectral(time,g3sig,ph2d)

    USE slint, ONLY: slint_init,bilinear_interp,bilinear_interp_uv
    USE sigio_module
    USE module_constants       ,ONLY: lat,lon,rd,cp,p1000,grvity,pi,sigak,sigbk
    USE module_control         ,ONLY: nip,nvlp1,ntra,dt,nvarsig,nvartracersig,nvlsig
    USE fimnamelist            ,ONLY: nvl,FixedGridOrder,ptop,glvl,curve, &
                                      ens_member,yyyymmddhhmm
    USE specmod, ONLY: init_spec_vars,gaulats,sptezv_s,sptez_s
    USE kinds, ONLY: i_kind,r_double,r_kind,r_single
    USE units,          ONLY: getunit, returnunit

!    USE funcphys,               ONLY: fpvsl

! External variable declarations:
!sms$distribute (dh,1) begin
    REAL   ,INTENT(IN)    :: ph2d(nip)
!sms$distribute end
!sms$distribute (dh,2) begin
    REAL   ,INTENT(IN)    :: g3sig(nvlsig+1,nip,nvarsig)
!sms$distribute end

    INTEGER, INTENT(IN)    :: time

!nip=10*(2**glvl)**2+2

! local vars

    TYPE(sigio_head) :: spectral_head
    TYPE(sigio_data) :: spectral_data

    REAL, DIMENSION(nip)        :: topo
    REAL, DIMENSION(nvlsig+1,nip) :: tvirt

    INTEGER :: jcap,nwaves,npts,ntracer,nvar,lonb,latb,&
      idrt,idvc,nlats_lin,nlons_lin
    INTEGER :: i,j,k,ivar,itracer,imo,jmo

    CHARACTER(len=10)  :: canaldate
    CHARACTER(len=250) :: filename
    CHARACTER(len=10)  :: tmpfile

    REAL :: exn,pdryini
    INTEGER :: indx,ierr,big_endian_lun

    REAL(r_single),ALLOCATABLE,DIMENSION(:) :: ug,vg
    REAL,DIMENSION(nip,2)                   :: latlonicos
    REAL,ALLOCATABLE,DIMENSION(:,:)         :: latlongauss
    INTEGER,DIMENSION(2)                    :: iens

    INTEGER :: levs,itrun,iorder,irealf,igen,latf,lonf,latr,lonr,&
      ntrac,icen2,idpp,idsl,idvm,idvt,idrun,idusr,ncldt,ixgr,ivs,&
      nvcoord

    LOGICAL :: tmpexist

    ntracer=nvarsig-nvartracersig+1

    idrt=4 !for gaussian grid

!assign values for the header

    fhour=time
    levs=nvlsig
    itrun=1
    iorder=2
    irealf=1
    igen=82
    ntrac=ntracer
    icen2=0
    idpp=0
    idsl=0
    idvc=2 !for hybrid vertical coord
    idvm=0
    idvt=21
    idrun=0
    idusr=0
    ncldt=1
    ixgr=0
    ivs=198410
    nvcoord=2
    pdryini = 98.24468

    big_endian_lun=getunit(82)
    if (big_endian_lun.lt.0) then
      big_endian_lun=getunit(76)
      if (big_endian_lun.lt.0) then
        write (*,'(a)') &
          'outvar_enkf_spectral: Could not obtain big-endian unit number.'
        stop
      endif
    endif

!sms$serial (<g3sig,ph2d,lat,lon,in>:default=ignore) begin

    iens=(/3,ens_member/)

    canaldate=yyyymmddhhmm(1:10)

    READ(canaldate(1:4),'(i4)') idate(4)
    READ(canaldate(5:6),'(i2)') idate(2)
    READ(canaldate(7:8),'(i2)') idate(3)
    READ(canaldate(9:10),'(i2)') idate(1)

    WRITE(timestr,'(i6.6)')time

    WRITE(6,*)'Writing out spectral enkfio_out file at time=',timestr
    CALL FLUSH(6)

    CALL glvl2gauss(glvl,jcap,lonb,latb,nlons_lin,nlats_lin,nwaves)

    nlonsgauss=lonb !nlons_lin
    nlatsgauss=latb !nlats_lin

!    nlonsgauss=nlons_lin
!    nlatsgauss=nlats_lin

!    latf=latb
!    lonf=lonb
!    latr=latb
!    lonr=lonb

    latf=nlatsgauss
    lonf=nlonsgauss
    latr=nlatsgauss
    lonr=nlonsgauss

    CALL init_spec_vars(nlonsgauss,nlatsgauss,jcap,idrt)

    CALL sigio_alhead(spectral_head,ierr,nvlsig,nvcoord,nvartracersig)

    spectral_head%fhour=fhour
    spectral_head%idate=idate
    spectral_head%jcap=jcap
    spectral_head%levs=levs
    spectral_head%itrun=itrun
    spectral_head%iorder=iorder
    spectral_head%irealf=irealf
    spectral_head%igen=igen
    spectral_head%latf=latf
    spectral_head%lonf=lonf
    spectral_head%latb=latb
    spectral_head%lonb=lonb
    spectral_head%latr=latr
    spectral_head%lonr=lonr
    spectral_head%ntrac=ntrac
    spectral_head%icen2=icen2
    spectral_head%iens=iens
    spectral_head%idpp=idpp
    spectral_head%idsl=idsl
    spectral_head%idvc=idvc
    spectral_head%idvm=idvm
    spectral_head%idvt=idvt
    spectral_head%idrun=idrun
    spectral_head%idusr=idusr
    spectral_head%ncldt=ncldt
    spectral_head%ixgr=ixgr
    spectral_head%ivs=ivs
    spectral_head%nvcoord=nvcoord
    spectral_head%pdryini=pdryini
    spectral_head%vcoord(:,1)=sigak
    spectral_head%vcoord(:,2)=sigbk

!    spectral_head%cfvars(1)=

    CALL sigio_aldata(spectral_head,spectral_data,ierr)

    imo=lonb
    jmo=latb

    npts = nlonsgauss*nlatsgauss

    ALLOCATE(latlongauss(npts,2),ug(npts),vg(npts))

    indx=0

    DO j = 1, nlatsgauss
      DO i = 1, nlonsgauss
        indx=indx+1
!          indx=(j - 1) * nlonsgauss + i
        latlongauss(indx,1)=ASIN(gaulats(j))
        latlongauss(indx,2)=2.*pi*float(i-1)/nlonsgauss
      END DO
    END DO

!    WRITE(cjcap,'(i6)')jcap
!    WRITE(cnlats,'(i6)')nlatsgauss
!    WRITE(cnlons,'(i6)')nlonsgauss

!    cstring=&
!         't'//TRIM(ADJUSTL(cjcap))//'_'//&
!         TRIM(ADJUSTL(cnlons))//'_'//&
!         TRIM(ADJUSTL(cnlats))//'_'//&
!         TRIM(ADJUSTL(canaldatep1))//'_'//&
!         TRIM(ADJUSTL(timestr))//'_'//&
!         cmember

    latlonicos(:,1)=lat
    latlonicos(:,2)=lon

    CALL slint_init(latlonicos,nip,latlongauss,npts)

! zero-initialize arrays

    spectral_data%d   = 0.
    spectral_data%hs  = 0.
    spectral_data%ps  = 0.
    spectral_data%q   = 0.
    spectral_data%t   = 0.
    spectral_data%xgr = 0.
    spectral_data%xss = 0.
    spectral_data%z   = 0.

!surface height

    topo=ph2d/grvity

    CALL bilinear_interp(topo,ug)
    CALL sptez_s(spectral_data%hs(:),ug,-1)

!p - needs to be log(ps)
    ivar=1

    CALL bilinear_interp(g3sig(1,:,ivar),ug)
    ug=LOG(ug/1000._r_kind)  !Pa -> cbar
    CALL sptez_s(spectral_data%ps(:),ug,-1)

!t - needs to be virtual temp
    ivar=3

    DO k=1,nvlsig
      tvirt(k,:)=g3sig(k,:,ivar)*(1.+0.6078*g3sig(k,:,ntracer))
      CALL bilinear_interp(tvirt(k,:),ug)
      CALL sptez_s(spectral_data%t(:,k),ug,-1)
    ENDDO

!u,v
    ivar=4

    DO k=1,nvlsig
      CALL bilinear_interp_uv(g3sig(k,:,ivar),ug,g3sig(k,:,ivar+1),vg)
      CALL sptezv_s(spectral_data%d(:,k),spectral_data%z(:,k),ug,vg,-1)
    ENDDO

!tracers

!order in gfs: qv,ozone,qw
!order in fim: qv,qw,ozone

!qv
    ivar=ntracer

    itracer=1

    DO k=1,nvlsig
      CALL bilinear_interp(g3sig(k,:,ivar),ug)
      CALL sptez_s(spectral_data%q(:,k,itracer),ug,-1)
    ENDDO

!ozone

    ivar=ntracer+2

    itracer=2

    DO k=1,nvlsig
      CALL bilinear_interp(g3sig(k,:,ivar),ug)
      CALL sptez_s(spectral_data%q(:,k,itracer),ug,-1)
    ENDDO

!qw

    ivar=ntracer+1

    itracer=3

    DO k=1,nvlsig
      CALL bilinear_interp(g3sig(k,:,ivar),ug)
      CALL sptez_s(spectral_data%q(:,k,itracer),ug,-1)
    ENDDO

    filename='fim_out_spectral_'//TRIM(timestr)

    CALL sigio_swohdc(big_endian_lun,filename,spectral_head,spectral_data,ierr)

    CALL sigio_axdata(spectral_data,ierr)

    DEALLOCATE(latlongauss,ug,vg)

!sms$serial end

    CALL returnunit(big_endian_lun)

    RETURN

  END SUBROUTINE outvar_enkf_spectral

END MODULE module_outvar_enkf_spectral
