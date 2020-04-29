!-------------------------------------------------------------------------------
! Module : sptsnv
! 
! This module contains subroutines and variables that are used to
! initialize and conduct spectral transforms for scalar and vector fields.
! 
! The main purpose of the creation of this module is to provide a mechanism 
! which only requires one time initialization. The sptez*_mod subroutines 
! are adapted from the similar sptez* subroutines in sp library.  The 
! entirety of the sp library is unchanged and thus could be updated readily 
! once its new version becomes available.
!
! History:
!   May 2014, N. Wang, initial version
!
!-------------------------------------------------------------------------------
MODULE sptsnv
  IMPLICIT NONE
  REAL,ALLOCATABLE,SAVE:: eps(:),epstop(:)
  REAL,ALLOCATABLE,SAVE:: enn1(:)
  REAL,ALLOCATABLE,SAVE:: elonn1(:)
  REAL,ALLOCATABLE,SAVE:: eon(:),eontop(:)
  REAL(8),ALLOCATABLE,SAVE:: afft(:)
  REAL,ALLOCATABLE,SAVE:: clat(:),slat(:),wlat(:)
  REAL,ALLOCATABLE,SAVE:: pln(:,:)
  REAL,ALLOCATABLE,SAVE:: plntop(:,:)
  LOGICAL, SAVE :: initialized = .false.
  
  PRIVATE
  PUBLIC :: init_sp, end_sp, sptez_mod, sptezm_mod, sptezmv_mod
CONTAINS

  SUBROUTINE init_sp (iromb,maxwv,idrt,imax,jmax)
    IMPLICIT NONE
    INTEGER :: iromb, maxwv,idrt,imax,jmax
    INTEGER jb, je

    jb=1
    je=(jmax+1)/2

    WRITE(6,*) 'Allocate arrays for spectral transform pre-calculated data.'

    ALLOCATE(eps((maxwv+1)*((iromb+1)*maxwv+2)/2),epstop(maxwv+1))
    ALLOCATE(enn1((maxwv+1)*((iromb+1)*maxwv+2)/2))
    ALLOCATE(elonn1((maxwv+1)*((iromb+1)*maxwv+2)/2))
    ALLOCATE(afft(50000+4*imax))
    ALLOCATE(eon((maxwv+1)*((iromb+1)*maxwv+2)/2),eontop(maxwv+1))
    ALLOCATE(clat(jb:je),slat(jb:je),wlat(jb:je))
    ALLOCATE(pln((maxwv+1)*((iromb+1)*maxwv+2)/2,jb:je))
    ALLOCATE(plntop(maxwv+1,jb:je))

    afft(:) = 0.
    CALL sptranf0(iromb,maxwv,idrt,imax,jmax,jb,je,eps,epstop,  &
          enn1,elonn1,eon,eontop,afft,clat,slat,wlat,pln,plntop)
    initialized = .true.

  END SUBROUTINE init_sp

  SUBROUTINE end_sp()
    IF (.not. initialized) THEN
      WRITE(6,*) 'end_sp in sptsnv module: Attempt to deallocate arrays that have not been allocated.'
      CALL FLUSH(6)
      RETURN
    ENDIF
    WRITE(6,*) 'Deallocate arrays for spectral transform pre-calculated data.'
    DEALLOCATE(eps,epstop,enn1,elonn1,afft,eon,eontop,pln,plntop)
    DEALLOCATE(clat,slat,wlat)
    initialized = .false.
  END SUBROUTINE end_sp

  SUBROUTINE sptez_mod(iromb,maxwv,idrt,imax,jmax,wave,grid,idir)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: iromb,maxwv,idrt,imax,jmax, idir
    REAL, INTENT(IN) :: wave((maxwv+1)*((iromb+1)*maxwv+2))
    REAL, INTENT(OUT) :: grid(imax,jmax)

    INTEGER :: ip, is, jn, js, kw, mx, kg, jb, je, jc 
    INTEGER :: ncpus


    IF (.not. initialized) THEN
      CALL init_sp(iromb,maxwv,idrt,imax,jmax)
    ENDIF

    mx=(maxwv+1)*((iromb+1)*maxwv+2)/2
    ip=1
    is=1
    jn=imax
    js=-jn 
    kw=2*mx
    kg=imax*jmax
    jb=1
    je=(jmax+1)/2
    jc=ncpus()

    CALL sptranf_s(iromb,maxwv,idrt,imax,jmax,1,ip,is,jn,js,kw,kg,jb,je,jc, &
                 wave,grid(1,1),grid(1,je+1))

  END SUBROUTINE sptez_mod 

  SUBROUTINE sptezm_mod(iromb,maxwv,idrt,imax,jmax,kmax,wave,grid,idir)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: iromb,maxwv,idrt,imax,jmax,kmax,idir
    REAL, INTENT(IN) :: wave((maxwv+1)*((iromb+1)*maxwv+2),kmax)
    REAL, INTENT(OUT) :: grid(imax,jmax,kmax)

    INTEGER :: ip, is, jn, js, kw, mx, kg, jb, je, jc 
    INTEGER :: ncpus


    IF (.not. initialized) THEN
      CALL init_sp(iromb,maxwv,idrt,imax,jmax)
    ENDIF

    mx=(maxwv+1)*((iromb+1)*maxwv+2)/2
    ip=1
    is=1
    jn=imax
    js=-jn 
    kw=2*mx
    kg=imax*jmax
    jb=1
    je=(jmax+1)/2
    jc=ncpus()

    CALL sptranf_s(iromb,maxwv,idrt,imax,jmax,kmax,ip,is,jn,js,kw,kg,jb,je,jc, &
                 wave,grid(1,1,1),grid(1,je+1,1))

  END SUBROUTINE sptezm_mod 

  SUBROUTINE sptezmv_mod(iromb,maxwv,idrt,imax,jmax,kmax,waved,wavez,gridu,gridv,idir)
    INTEGER, INTENT(IN) :: iromb,maxwv,idrt,imax,jmax,kmax,idir
    REAL, INTENT(IN) :: waved((maxwv+1)*((iromb+1)*maxwv+2),kmax)
    REAL, INTENT(IN) :: wavez((maxwv+1)*((iromb+1)*maxwv+2),kmax)
    REAL, INTENT(OUT) ::  gridu(imax,jmax,kmax)
    REAL, INTENT(OUT) :: gridv(imax,jmax,kmax)

    INTEGER :: ip, is, jn, js, kw, mx, kg, jb, je, jc 
    INTEGER :: ncpus

    IF (.not. initialized) THEN
      CALL init_sp(iromb,maxwv,idrt,imax,jmax)
      initialized = .true.
    ENDIF

    mx=(maxwv+1)*((iromb+1)*maxwv+2)/2
    ip=1
    is=1
    jn=imax
    js=-jn
    kw=2*mx
    kg=imax*jmax
    jb=1
    je=(jmax+1)/2
    jc=ncpus()

    CALL sptranfv_s(iromb,maxwv,idrt,imax,jmax,kmax,ip,is,jn,js,kw,kg,jb,je,jc,  &
                  waved,wavez,gridu(1,1,1),gridu(1,je+1,1),gridv(1,1,1),gridv(1,je+1,1))

  END SUBROUTINE sptezmv_mod

!spectral synthesis subroutine (for scalar variable)
  SUBROUTINE sptranf_s(iromb,maxwv,idrt,imax,jmax,kmax,ip,is,jn,js,kw,kg, &
                     jb,je,jc,wave,gridn,grids)
    INTEGER, INTENT(IN) :: iromb,maxwv,idrt,imax,jmax,kmax,ip,is,jn,js
    INTEGER, INTENT(IN) :: kw,kg,jb,je,jc
    REAL, INTENT(IN) :: wave((maxwv+1)*((iromb+1)*maxwv+2)*kmax)
    REAL, INTENT(OUT) :: gridn(*), grids(*)

    REAL(8) afft_lc(50000+4*imax)
    REAL g(imax,2)
    INTEGER :: idir
    INTEGER :: i, k, j, kws, wtop, ijkn, ijks
    INTEGER :: mp

    mp = 0
    afft_lc = afft 
    idir = 1

    DO k=1,kmax
      kws=(k-1)*kw
      wtop=0
      DO j=jb,je
        call sptranf1(iromb,maxwv,idrt,imax,jmax,j,j,  &
                      eps,epstop,enn1,elonn1,eon,eontop,  &
                      afft_lc,clat(j),slat(j),wlat(j),  &
                      pln(1,j),plntop(1,j),mp,  &
                      wave(kws+1),wtop,g,idir)
        IF(ip.eq.1.and.is.eq.1) THEN
          DO i=1,imax
            ijkn=i+(j-jb)*jn+(k-1)*kg
            ijks=i+(je-j)*jn+(k-1)*kg
            gridn(ijkn)=g(i,1)
            grids(ijks)=g(i,2)
          ENDDO
        ELSE
          DO i=1,imax
            ijkn=mod(i+ip-2,imax)*is+(j-jb)*jn+(k-1)*kg+1
            ijks=mod(i+ip-2,imax)*is+(j-jb)*js+(k-1)*kg+1
            gridn(ijkn)=g(i,1)
            grids(ijks)=g(i,2)
          ENDDO
        ENDIF
      ENDDO
    ENDDO

  END SUBROUTINE 

!spectral analysis subroutine (for scalar variable)
  SUBROUTINE sptranf_a(iromb,maxwv,idrt,imax,jmax,kmax,ip,is,jn,js,kw,kg, &
                     jb,je,jc,wave,gridn,grids)
    INTEGER, INTENT(IN) :: iromb,maxwv,idrt,imax,jmax,kmax,ip,is,jn,js
    INTEGER, INTENT(IN) :: kw,kg,jb,je,jc
    REAL, INTENT(OUT) :: wave((maxwv+1)*((iromb+1)*maxwv+2)*kmax)
    REAL, INTENT(IN) :: gridn(*), grids(*)

    REAL(8) afft_lc(50000+4*imax)
    REAL g(imax,2)
    INTEGER :: idir
    INTEGER :: i, k, j, kws, wtop, ijkn, ijks
    INTEGER :: mp

    mp = 0
    afft_lc = afft 
    idir = -1
    wave=0

    DO k=1,kmax
      kws=(k-1)*kw
      wtop=0
      DO j=jb,je
        IF(wlat(j).gt.0.) THEN
          IF(ip.eq.1.and.is.eq.1) THEN
            DO i=1,imax
              ijkn=i+(j-jb)*jn+(k-1)*kg
              ijks=i+(j-jb)*js+(k-1)*kg
              g(i,1)=gridn(ijkn)
              g(i,2)=grids(ijks)
            ENDDO
          ELSE
            DO i=1,imax
              ijkn=mod(i+ip-2,imax)*is+(j-jb)*jn+(k-1)*kg+1
              ijks=mod(i+ip-2,imax)*is+(j-jb)*js+(k-1)*kg+1
              g(i,1)=gridn(ijkn)
              g(i,2)=grids(ijks)
            ENDDO
          ENDIF
          CALL sptranf1(iromb,maxwv,idrt,imax,jmax,j,j,  &
                       eps,epstop,enn1,elonn1,eon,eontop,  &
                       afft_lc,clat(j),slat(j),wlat(j),  &
                       pln(1,j),plntop(1,j),mp,  &
                       wave(kws+1),wtop,g,idir)
        ENDIF
      ENDDO
    ENDDO
  END SUBROUTINE 

!spectral synthesis subroutine (for vectorial variable)
  SUBROUTINE sptranfv_s(iromb,maxwv,idrt,imax,jmax,kmax, ip,is,jn,js,kw,kg,jb,je,jc,   &
                      waved,wavez,gridun,gridus,gridvn,gridvs)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: iromb,maxwv,idrt,imax,jmax,kmax,ip,is,jn,js
    INTEGER, INTENT(IN) :: kw,kg,jb,je,jc
    REAL, INTENT(IN) :: waved((maxwv+1)*((iromb+1)*maxwv+2)*kmax)
    REAL, INTENT(IN) :: wavez((maxwv+1)*((iromb+1)*maxwv+2)*kmax)
    REAL, INTENT(OUT):: gridun(*),gridus(*),gridvn(*),gridvs(*)

    INTEGER mp(2)
    REAL w((maxwv+1)*((iromb+1)*maxwv+2)/2*2,2)
    REAL wtop(2*(maxwv+1),2)
    REAL g(imax,2,2)
    REAL winc((maxwv+1)*((iromb+1)*maxwv+2)/2*2,2)

    REAL(8) afft_lc(50000+4*imax)
    INTEGER :: idir
    INTEGER :: i, k, j, kws, ijkn, ijks, mx

    mx=(maxwv+1)*((iromb+1)*maxwv+2)/2
    mp=1
    afft_lc=afft
    idir = 1

    DO k=1,kmax
      kws=(k-1)*kw
      CALL spdz2uv(iromb,maxwv,enn1,elonn1,eon,eontop,waved(kws+1),wavez(kws+1),  &
                   w(1,1),w(1,2),wtop(1,1),wtop(1,2))
      DO j=jb,je
        CALL sptranf1(iromb,maxwv,idrt,imax,jmax,j,j,eps,epstop,enn1,elonn1,eon,eontop, &
                      afft_lc,clat(j),slat(j),wlat(j),pln(1,j),plntop(1,j),mp, &
                      w(1,1),wtop(1,1),g(1,1,1),idir)
        CALL sptranf1(iromb,maxwv,idrt,imax,jmax,j,j,eps,epstop,enn1,elonn1,eon,eontop, &
                      afft_lc,clat(j),slat(j),wlat(j),pln(1,j),plntop(1,j),mp, &
                      w(1,2),wtop(1,2),g(1,1,2),idir)
        IF(ip.eq.1.and.is.eq.1) THEN
          DO i=1,imax
            ijkn=i+(j-jb)*jn+(k-1)*kg
            ijks=i+(je-j)*jn+(k-1)*kg
            gridun(ijkn)=g(i,1,1)
            gridus(ijks)=g(i,2,1)
            gridvn(ijkn)=g(i,1,2)
            gridvs(ijks)=g(i,2,2)
          ENDDO
        ELSE
          DO i=1,imax
            ijkn=mod(i+ip-2,imax)*is+(j-jb)*jn+(k-1)*kg+1
            ijks=mod(i+ip-2,imax)*is+(j-jb)*js+(k-1)*kg+1
            gridun(ijkn)=g(i,1,1)
            gridus(ijks)=g(i,2,1)
            gridvn(ijkn)=g(i,1,2)
            gridvs(ijks)=g(i,2,2)
          ENDDO
        ENDIF
      ENDDO
    ENDDO

  END SUBROUTINE sptranfv_s

!spectral analysis subroutine (for vectorial variable)
  SUBROUTINE sptranfv_a(iromb,maxwv,idrt,imax,jmax,kmax, ip,is,jn,js,kw,kg,jb,je,jc,   &
                      waved,wavez,gridun,gridus,gridvn,gridvs)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: iromb,maxwv,idrt,imax,jmax,kmax,ip,is,jn,js
    INTEGER, INTENT(IN) :: kw,kg,jb,je,jc
    REAL, INTENT(OUT) :: waved((maxwv+1)*((iromb+1)*maxwv+2)*kmax)
    REAL, INTENT(OUT) :: wavez((maxwv+1)*((iromb+1)*maxwv+2)*kmax)
    REAL, INTENT(IN):: gridun(*),gridus(*),gridvn(*),gridvs(*)

    REAL(8) afft_lc(50000+4*imax)

    INTEGER mp(2)
    REAL w((maxwv+1)*((iromb+1)*maxwv+2)/2*2,2)
    REAL wtop(2*(maxwv+1),2)
    REAL g(imax,2,2)
    REAL winc((maxwv+1)*((iromb+1)*maxwv+2)/2*2,2)

    INTEGER :: idir
    INTEGER :: i, k, j, kws, ijkn, ijks, mx

    mx=(maxwv+1)*((iromb+1)*maxwv+2)/2
    mp=1
    afft_lc=afft
    idir = -1
    waved=0
    wavez=0

    DO k=1,kmax
      kws=(k-1)*kw
      w=0
      wtop=0
      DO j=jb,je
        IF(wlat(j).gt.0.) THEN
          IF(ip.eq.1.and.is.eq.1) THEN
            DO i=1,imax
              ijkn=i+(j-jb)*jn+(k-1)*kg
              ijks=i+(j-jb)*js+(k-1)*kg
              g(i,1,1)=gridun(ijkn)/clat(j)**2
              g(i,2,1)=gridus(ijks)/clat(j)**2
              g(i,1,2)=gridvn(ijkn)/clat(j)**2
              g(i,2,2)=gridvs(ijks)/clat(j)**2
            ENDDO
          ELSE
            DO i=1,imax
              ijkn=mod(i+ip-2,imax)*is+(j-jb)*jn+(k-1)*kg+1
              ijks=mod(i+ip-2,imax)*is+(j-jb)*js+(k-1)*kg+1
              g(i,1,1)=gridun(ijkn)/clat(j)**2
              g(i,2,1)=gridus(ijks)/clat(j)**2
              g(i,1,2)=gridvn(ijkn)/clat(j)**2
              g(i,2,2)=gridvs(ijks)/clat(j)**2
            ENDDO
          ENDIF
          CALL sptranf1(iromb,maxwv,idrt,imax,jmax,j,j,eps,epstop,enn1,elonn1,eon,eontop, &
                        afft_lc,clat(j),slat(j),wlat(j),pln(1,j),plntop(1,j),mp, &
                        w(1,1),wtop(1,1),g(1,1,1),idir)
          CALL sptranf1(iromb,maxwv,idrt,imax,jmax,j,j,eps,epstop,enn1,elonn1,eon,eontop, &
                        afft_lc,clat(j),slat(j),wlat(j),pln(1,j),plntop(1,j),mp, &
                        w(1,2),wtop(1,2),g(1,1,2),idir)
        ENDIF
      ENDDO
      CALL spuv2dz(iromb,maxwv,enn1,elonn1,eon,eontop,w(1,1),w(1,2),wtop(1,1),wtop(1,2), &
                   winc(1,1),winc(1,2))
      waved(kws+1:kws+2*mx)=waved(kws+1:kws+2*mx)+winc(1:2*mx,1)
      wavez(kws+1:kws+2*mx)=wavez(kws+1:kws+2*mx)+winc(1:2*mx,2)
    ENDDO

  END SUBROUTINE sptranfv_a

END MODULE sptsnv


