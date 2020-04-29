module mdul_mkmodisland
contains
subroutine mkmodisland(vtype2d,nip,setO2L,setL2O,nset,no2l,nl2o)

  use libmf
  use netcdf
  use module_control, only: npp
  use fimnamelist,    only: glvl,landsmoothfact,landdatdir,	&
                      landdatfile,landglvldir,niland,njland
  use module_constants, only: lat,lon,nprox,proxs,prox

  ! -- only needed inside fim
  !use global_bounds,      only: set_global_bounds, set_global_bounds_halo, ims, ime, ips, ipe, ihs, ihe

  implicit none

  integer, intent(IN) :: nip,nset
  real, intent(INOUT) :: vtype2d(nip)
  integer,intent(OUT) :: setO2L(nset,2),setL2O(nset,2),no2l,nl2o

  integer i,j,k,l,m,n,ips,ipe,irc,iuglvl,isn,ipn
  integer :: ni,nj,ncid,varid
  integer,allocatable :: iwork2d(:)
  integer verb,diagcomp,verbcomp

  character(24)   :: varname
  character(16)   :: header
  character gpath*120,gfsfcpath*120,spath*120,opath*120,obspath*120
  character cglvl*1,version*2,compdata*3

  real, allocatable :: vtype(:)
  real(kind=4), allocatable :: land(:,:),dum(:,:)
  real gridscale
  real blat,blon,dlat,dlon,dlonkm,dlatkm
  real t1,t2,t3,t4,t5
  real radinf,radinfj,vtypecomp
  real rlat,rlon

  integer,allocatable :: setO2L999(:,:)
  integer,allocatable :: setL2O999(:,:)

  common /grid/ blat,dlat,blon,dlon,radinf,radinfj,cglvl,version

!sms$ignore begin

! -- turn on output here
  !
  verb=0

  ! -- turn on verbose output comp2modis
  !
  verbcomp=1

  ! -- diagcomp writes a gras .obs file in comp2modis
  !
  diagcomp=1

  version='v1'

  allocate (vtype(nip),iwork2d(nip),stat=irc)

  ni=niland
  nj=njland

  if(glvl <= 9) then
     gridscale=15.0*(10-glvl) 
  else
     print*,'EEE invalid glvl must be <= 9'
     stop 'bad glvl'
  endif


  write(cglvl,'(i1)') glvl

  varname="LU_INDEX"

  blon=-179.9750
  blat=-89.9750

  dlon=360.0/(ni)
  dlat=180.0/(nj)

  blon=-180.0+dlon*0.5
  blat=-90.0+dlat*0.5

  dlonkm=dlon*111.0
  dlatkm=dlat*111.0

  if(verb.eq.1) then
     print*,'dlonkm: ',dlonkm,' dlatkm: ',dlatkm
     print*,'blon: ',blon,' blat: ',blat
     print*,'dlon: ',dlon,' dlat: ',dlat,' elon: ',blon+(ni-1)*dlon,' elat: ',blat+(nj-1)*dlat
  endif

  ! from prep/getlvl.F90
  !
  allocate (iwork2d(nip),stat=irc)
  allocate(lat(nip),stat=irc)
  allocate(lon(nip),stat=irc)
  allocate(nprox(nip),stat=irc)
  allocate(prox(npp,nip),stat=irc)
  allocate(proxs(npp,nip),stat=irc)
  
  gpath=trim(landglvldir)//'glvl.dat'
  spath=trim(landdatdir)//'/'//trim(landdatfile)

  if(verb.eq.1) then
     write(*,'(a,a)') '   glvl: ',gpath
     write(*,'(a,a)') ' source: ',spath
     write(*,'(a,a)') '    out: ',opath
  endif

  iuglvl=12

  !--  read in the glvl.dat for lat/lons
  !
  write(6,*)'mkmodisland: trying to open file=',trim(gpath)
  open(iuglvl,file=gpath,form='unformatted',status='old',err=800)
  
  read(iuglvl) header
  read(iuglvl) header

  ! -- only needed in fim
  !call set_global_bounds ()  ! define ims, ime, ips, ipe
  
  ips=1
  ipe=nip

  ! -- read as in fim/horizontal/dyn_init.f90
  !
  read (iuglvl,err=90) lat
  read (iuglvl,err=90) lon
  read (iuglvl,err=90) nprox

  do isn=1,size(proxs,1)
     read (iuglvl,err=90) iwork2d
      do ipn=ips,ipe
        proxs(isn,ipn) = iwork2d(ipn)
      end do
    end do

    do isn=1,size(prox,1)
      read (iuglvl,err=90) iwork2d
      do ipn=ips,ipe
        prox(isn,ipn) = iwork2d(ipn)
      end do
    end do

  close(iuglvl)

  call cpu_time(t2)
  
  ! -- readn the .nc file with lu_type = landuse index
  !
  allocate(land(ni,nj),stat=irc)

  call check( nf90_open(spath, NF90_NOWRITE, ncid) )

  ! Get the varid of the data variable, based on its name.
  call check( nf90_inq_varid(ncid,varname, varid) )

  ! Read the data.
  ! call check( nf90_get_var(ncid, varid, data_in, start=start, count=count) )
  call check( nf90_get_var(ncid, varid, land) )

  if(verb.eq.1) then
     call qprntn(land,varname,1,1,ni,nj,50,6)
  endif

  ! -- allocate arrays to hold indices to change gfs land to ocean and vice versa
  !
  allocate(setO2L999(nset,2),stat=irc)
  allocate(setL2O999(nset,2),stat=irc)

  call cpu_time(t3)

  if(verb.eq.1)  print*,'TTT(mkmodisland) time in .nc read',(t3-t2)

  ! --  analyze to icos grid
  !
  vtype=0.0
  radinf=gridscale*landsmoothfact*0.5
  radinfj=(radinf*km2deglat)/dlat

  do i=1,nip
     rlat=lat(i)*rad2deg
     rlon=lon(i)*rad2deg
     if(rlon > 180.0) rlon=rlon-360.0
     vtypecomp=0.0
     call anlland(vtypecomp,land,rlat,rlon,ni,nj,vtype(i))
     
     !write(*,'("i:",i6,2x,4f12.3,2x,"ib,ie: ",2(i5,1x),"jb,je: ",2(i5,1x))') i,rlat,rlon,ri,rj,ib,ie,jb,je
  enddo

  call cpu_time(t4)
  if(verb.eq.1) print*,'TTT(mkmodisland) time in analysis',(t4-t3)

  call comp2modis(vtype,vtype2d,prox,lat,lon,		&
       setO2L,setL2O,setO2L999,setL2O999,no2l,nl2o,	&
       nset,npp,nip,verbcomp,diagcomp)

  ! replace  vtype2d with the correctd vtype
  !
  do i=1,nip
     vtype2d(i)=vtype(i)
  enddo

  call cpu_time(t5)
  if(verb.eq.1) print*,'TTT(mkmodisland) time in comp2modis',(t5-t4)

  goto 900

800 continue
  print*,'error in open of g?glvl.dat'
  goto 900

810 continue
  print*,'error in open of gfsfc.dat'
  goto 900
  
801 continue
  print*,'error in open of land.dat'
  goto 900
  
802 continue
  print*,'read end'
  goto 900
  
803 continue
  print*,'read err'
  goto 900
  
90 continue
  print*,'read err'
  goto 900
  
900 continue

  print*,'III(mkmodisland) completed sucessfully total cpu time: ',(t5-t2)
  return

!sms$ignore end

end subroutine mkmodisland


  subroutine check(status)
  use netcdf
    implicit none
    integer, intent ( in) :: status
!sms$ignore begin
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
!sms$ignore end
  end subroutine check  



subroutine anlland(z0test,land,rlat,rlon,ni,nj,z0out)

 use libmf

  implicit none

  real   ,intent(IN) :: z0test,rlat,rlon
  real  ,intent(OUT) :: z0out
  integer,intent(IN) :: ni,nj
  real(kind=4),intent(IN) :: land(ni,nj)
  real(kind=4) :: z0s(ni*20),o0s(ni*20)
  real         :: blat,dlat,blon,dlon,radinf,radinfj,ri,rj,rlatfact,	&
                  rjb,rje,rib,rie,z0bar,z0rms,distmin,tlat,tlon,dx,dy,	&
                  tval,odominant,dz0,dist,z0min,radinfi
  integer      :: iz0s(ni*20),n0s(ni*20),ib,ie,jb,je,nz0,ii,jj,		&
                  ndominant,ntot,n0,ncnt,n
  character cglvl*1,version*2

  common /grid/ blat,dlat,blon,dlon,radinf,radinfj,cglvl,version

  integer verb

!sms$ignore begin
  verb=0

  ri=(rlon-blon)/dlon+1.0
  rj=(rlat-blat)/dlat+1.0
  
  rlatfact=cos(rlat*deg2rad)

  radinfi=0.0
  if(rlatfact > 0.0) radinfi=radinfj/rlatfact
  
  rjb=rj-radinfj-1.0
  rje=rj+radinfj+1.0

  if(rjb < 1) rjb=1.0
  if(rje > nj) rje=nj
  
  if(radinfi == 0) then
     rib=1
     rie=ni
  else
     rib=ri-radinfi-1.0
     rie=ri+radinfi+1.0
  endif

  if(rlatfact == 0.0) then
     rib=1.0
     rie=ni/2
  else
     if(rib < 1) rib=1.0
     if(rie > ni) rie=ni
  endif
  
  ib=nint(rib)
  ie=nint(rie)

  jb=nint(rjb)
  je=nint(rje)

  if(ib < 1) ib=1
  if(ie > ni) ie=ni

  if(jb < 1) jb=1
  if(je > nj) je=nj


  z0bar=0.0
  z0rms=0.0
  nz0=0

  distmin=1e20

  do ii=ib,ie
     do jj=jb,je
        tlat=blat+(jj-1)*dlat
        tlon=blon+(ii-1)*dlon
        dy=(tlat-rlat)
        dx=(tlon-rlon)*rlatfact
        dist=sqrt(dx*dx+dy*dy)*deglat2km
        if(dist < distmin) then
           distmin=dist
           z0min=land(ii,jj)
        endif

        if(dist <= radinf) then
           nz0=nz0+1
           z0s(nz0)=land(ii,jj)
        endif

     end do
  end do

!  do n=1,nz0
!     z0bar=z0bar+z0s(n)
!  enddo
 
  if(nz0 == 1) then
     z0bar=z0s(1)
  elseif (nz0 > 1) then
     call indexx(nz0,z0s,iz0s)

     n0=1
     o0s(n0)=z0s(iz0s(n0))

     !print*, '000',iz0s(n0)

     ncnt=1

     do n=2,nz0
        tval=z0s(iz0s(n))
        !print*, '111',n,tval
        if(tval /= o0s(n0)) then
           n0s(n0)=ncnt
           !print*,'qqq',n0,o0s(n0),n0s(n0)
           n0=n0+1
           o0s(n0)=z0s(iz0s(n))
           ncnt=1
        else
           ncnt=ncnt+1
        endif

     enddo

     n0s(n0)=ncnt

     ndominant=-999
     ntot=0
     do n=1,n0
        if(n0s(n) > ndominant) then
           ndominant=n0s(n)
           odominant=o0s(n)
        endif
        ntot=ntot+n0s(n)
     enddo

     !if(rlat >= -90.0 .and. rlat <= 90.0 .and. rlon >= 0.0 .and. rlon <= 360.0)  &
     !if(rlat >= 20.0 .and. rlat <= 90.0 .and. rlon >= 180.0 .and. rlon <= 360.0)  &
     !     write(*,'(a,2x,2(f7.2,1x),2x,3(i5,2x),5x,(i5,1x,f3.0),5x,10(i3,2x,f3.0,1x,i4,4x))') &
     !     'rlat,rlon,nz0,n0=ntot,ndom',rlat,rlon,nz0,n0,ntot,ndominant,odominant,(n,o0s(n),n0s(n),n=1,n0)

     z0bar=odominant

  else
     print*,'in anlland no points! nz0 = 0',rlat,rlon
     stop 'no points'
  endif

  if(verb.eq.1) then
     dz0=z0test-z0bar
     write(*,'("final",2x,2(f7.2,1x),2x,i6,1x,2(f7.2,1x),2x,(2(f7.2,1x)),2x,4(i5,1x))') &
           rlat,rlon,nz0,z0bar,dz0,distmin,z0min,ib,(ie-ib),jb,(je-jb)
  endif

  z0out=z0bar

  return

!sms$ignore end

end subroutine anlland


subroutine comp2modis(modis,gfs,prox,lat,lon,		&
     setO2L,setL2O,setO2L999,setL2O999,no2l,nl2o,	&
     nset,npp,nip,verb,diag)

  use libmf

  implicit none

  character(8) :: stid
  character cglvl*1,version*2
  
  
  integer,intent(IN)    :: nip,npp,diag,verb,nset
  real   ,intent(IN)    :: lat(nip),lon(nip),gfs(nip)
  real   ,intent(INOUT) :: modis(nip)
  integer,intent(OUT)   :: setO2L(nset,2),setL2O(nset,2),		&
                           setO2L999(nset,2),setL2O999(nset,2)
  integer,intent(OUT)   :: no2l,nl2o
  
  real    :: radinf,radinfj,vtypecomp,vtypecorr,ostate,tim
  real    :: rlat,rlon
  real    :: blat,blon,dlat,dlon,dlonkm,dlatkm

  integer :: testType,istate,no2l999,nl2o999,iprox
  integer :: i,j,k,l,m,n,ips,ipe,nlev,nflag,iok

  integer :: prox(npp,nip)
  integer :: nextps(npp),nextps2(npp,npp),nextps3(npp,npp,npp),		&
             nextps4(npp,npp,npp,npp)
  integer :: nfim,nmodis,ntot,iuobsdiff,iModis,vtypeIN
  
  common /grid/ blat,dlat,blon,dlon,radinf,radinfj,cglvl,version

!sms$ignore begin

  no2l=0
  nl2o=0

  ntot=0

  no2l999=0
  nl2o999=0
  
  iuobsdiff=17
  ! -- compiler dependent feature access='stream' put  
!tgs Lahey compiler gives an error in the following 3 lines, for now commented
!out
!  if(diag == 1) then
!     open(iuobsdiff,file='g'//cglvl//'vtype'//version//'.diff.obs',form='unformatted',access='stream',status='unknown',err=812)
!  endif

  do i=1,nip

     rlat=lat(i)*rad2deg
     rlon=lon(i)*rad2deg

     istate=0
     vtypeIN=modis(i)

     ! -- case modis h2o
     !
     iModis=testModisO(modis(i))
     if(iModis == 1) then

        ! -- check if gfs h2o; if not scan for an h2o point...
        !
        iok=testGfsO(gfs(i))
        if(iok == 0) then

           call scanProx(0,i,lat,lon,rlat,rlon,modis,gfs,prox,		&
                nextps,nextps2,nextps3,nextps4,				&
                iprox,istate,nip,npp,verb)
              
           no2l=no2l+1
           ntot=ntot+1
           
           if(istate == -999) then 
              ! -- modis has a lake, but no gfs lake points, set to nearest modis land
              !
              call scanProx(11,i,lat,lon,rlat,rlon,modis,gfs,prox,	&
                   nextps,nextps2,nextps3,nextps4,			&
                   iprox,istate,nip,npp,verb)

              if(istate == -999) then
                 print*,'LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL i: ',i,rlat,rlon,modis(i),gfs(i)
                 stop 'lllll'
              endif

              if(verb.eq.1) write(*,'(a,i10,2x,4(f7.2,2x),2x,a,i5)') 'MMMLLL bbbbbbbbbbbbbbb MMMLLL i,lat,lon:', &
                                    i,rlat,rlon,modis(i),gfs(i),' istate: ',istate
              modis(i)=modis(iprox)
              istate=101
              if(verb.eq.1) write(*,'(a,i10,2x,4(f7.2,2x),2x,a,i5)') 'MMMLLL make modis land MMMLLL i,lat,lon:', &
                                    i,rlat,rlon,modis(i),gfs(i),' istate: ',istate

              nl2o999=nl2o999+1
              setO2L999(nl2o999,1)=i
              setO2L999(nl2o999,2)=iprox

              setO2L(no2l,1)=i
              setO2L(no2l,2)=iprox

           else

              setO2L(no2l,1)=i
              setO2L(no2l,2)=iprox

           endif

           if(verb.eq.1) then
              write(*,'(a,i10,2x,4(f7.2,2x),2x,a,i5)') 'GGGLLL GFS == LAND && MODIS == sea      i,lat,lon:', &
                      i,rlat,rlon,modis(i),gfs(i),' istate: ',istate
              write(*,'(a,i10,2x,4(f7.2,2x),2x,a,i5)') 'GGGLLL GFS == LAND && MODIS == sea  iprox,lat,lon:', &
                      iprox,rlat,rlon,modis(i),gfs(iprox),' istate: ',istate
           endif

           if(istate /= -999) then
           endif

        endif

        ! -- case modis land
        !

     else

        ! -- check if gfs land, if not scan for a nearby land point...
        !

        iok=(testGfsL(gfs(i)))

        if(iok == 0) then

           nl2o=nl2o+1
           ntot=ntot+1

           call scanProx(1,i,lat,lon,rlat,rlon,modis,gfs,prox,    &
                nextps,nextps2,nextps3,nextps4,                   &
                iprox,istate,                           &
                nip,npp,verb)
           
           if(istate == -999) then 

              istate=-9999

              ! -- modis land (islands) but gfs ocean; - set modis to ocean
              !
              if(verb.eq.1) write(*,'(a,i10,2x,4(f7.2,2x),2x,a,i5)') 'MMMOOO bbbbbbbbbbbbbbb MMMOOO i,lat,lon:', &
                                    i,rlat,rlon,modis(i),gfs(i),' istate: ',istate
              modis(i)=17
              istate=201
              if(verb.eq.1) write(*,'(a,i10,2x,4(f7.2,2x),2x,a,i5)') 'MMMOOO make modis land MMMOOO i,lat,lon:', &
                                    i,rlat,rlon,modis(i),gfs(i),' istate: ',istate

              setL2O(nl2o,1)=i
              setL2O(nl2o,2)=i

              no2l999=no2l999+1
              setL2O999(no2l999,1)=i
              setL2O999(no2l999,2)=i

           else

              setL2O(nl2o,1)=i
              setL2O(nl2o,2)=iprox

           endif

           if(verb.eq.1) then
              write(*,'(a,i10,2x,4(f7.2,2x),2x,a,i5)') 'GGGOOO GFS == sea  && MODIS == LAND     i,lat,lon:', &
                      i,rlat,rlon,modis(i),gfs(i),' istate: ',istate
              write(*,'(a,i10,2x,4(f7.2,2x),2x,a,i5)') 'GGGOOO GFS == sea  && MODIS == LAND iprox,lat,lon:', &
                      iprox,rlat,rlon,modis(i),gfs(iprox),' istate: ',istate
           endif

        endif

     endif


     if(diag == 1 .and. istate /= 0) then
        write(stid,'(a1,i7.7)') 's',i
        tim = 0.0 
        nlev = 1 
        nflag = 1

        if(istate /= 0) then
           if(mod(ntot,200) == 0) print*,'conflict point i: ',i,'istate: ',istate
        endif

        ostate=istate*1.0
        vtypeCORR=modis(i)
        write (iuobsdiff) stid,rlat,rlon,tim,nlev,nflag 
        write (iuobsdiff) ostate,gfs(i),modis(i)

     endif

  enddo

  if(diag == 1) then
     nlev = 0 
     write (iuobsdiff) stid,rlat,rlon,tim,nlev,nflag 
     close(iuobsdiff)
  endif

  print*
  write(*,'(a,2(i5,1x),a,2(i5,1x),a,i5)') 'LLL no2l: ',no2l,no2l999,'OOO nl2o: ',nl2o,nl2o999,'NNN ntot: ',ntot
  print*
  go to 900
812 continue
  print*,'error in open of output '//cglvl//'z0.obs'
  goto 900

900 continue

!sms$ignore end

end subroutine comp2modis


subroutine scanProx(testType,i,lat,lon,rlat,rlon,modis,gfs,prox,	&
     nextps,nextps2,nextps3,nextps4,					&
     iprox,istate,nip,npp,verb)

  use libmf,only:rad2deg

  implicit none

  integer,intent(IN) :: nip,npp
  integer,intent(IN) :: testType,i
  real   ,intent(IN) :: rlat,rlon
  real   ,intent(IN) :: gfs(nip),modis(nip),lat(nip),lon(nip)
  integer,intent(IN) :: prox(npp,nip)
  integer :: nextps(npp),nextps2(npp,npp),nextps3(npp,npp,npp),		&
             nextps4(npp,npp,npp,npp),ii,jj,j,k,l,m
  integer :: verb,istate,iset,testO,iprox

  istate=-999

  iset=0
  do j=1,npp
     nextps(j)=prox(j,i)

     if(testType == 0)  testO=testGfsO(gfs(prox(j,i)))
     if(testType == 1)  testO=testGfsL(gfs(prox(j,i)))
     if(testType == 10) testO=testModisO(modis(prox(j,i)))
     if(testType == 11) testO=testModisL(modis(prox(j,i)))

     if(testO == 1) then
        iset=1
        iprox=prox(j,i)
        if(testType == 0)  istate=-1
        if(testType == 1)  istate=1
        if(testType == 10) istate=-11
        if(testType == 11) istate=11
        if(verb.eq.1) write(*,'(a,1x,2(i7,1x),2x,3(f6.1,1x))') 'FF1111',i,iprox,gfs(iprox),modis(i),modis(iprox)
        return
     endif
  end do

  ! -- search in 2nd ring
  !
  if(iset == 0) then

     do j=1,npp
        ii=nextps(j)
        do jj=1,npp
           nextps2(j,jj)=prox(jj,ii)

           if(testType == 0)  testO=testGfsO(gfs(prox(jj,ii)))
           if(testType == 1)  testO=testGfsL(gfs(prox(jj,ii)))
           if(testType == 10) testO=testModisO(modis(prox(jj,ii)))
           if(testType == 11) testO=testModisL(modis(prox(jj,ii)))

           if(testO == 1) then
              iset=1
              iprox=prox(jj,ii)
              if(testType == 0)  istate=-2
              if(testType == 1)  istate=2 
              if(testType == 10) istate=-12 
              if(testType == 11) istate=12 
              if(verb.eq.1) write(*,'(a,1x,2(i7,1x),2x,3(f6.1,1x))') 'FF22222222', &
                                    i,iprox,gfs(iprox),modis(i),modis(iprox)
              return
           endif
        enddo
     enddo

  endif

  ! -- search in 3nd ring
  !
  if(iset == 0) then

     do k=1,npp
        do j=1,npp
           ii=nextps2(k,j)
           do jj=1,npp
              nextps3(k,j,jj)=prox(jj,ii)

              if(testType == 0)  testO=testGfsO(gfs(prox(jj,ii)))
              if(testType == 1)  testO=testGfsL(gfs(prox(jj,ii)))
              if(testType == 10) testO=testModisO(modis(prox(jj,ii)))
              if(testType == 11) testO=testModisL(modis(prox(jj,ii)))

              if(testO == 1) then
                 iset=1
                 iprox=prox(jj,ii)
                 if(testType == 0)  istate=-3
                 if(testType == 1)  istate=3
                 if(testType == 10) istate=-13
                 if(testType == 11) istate=13
                 if(verb.eq.1) write(*,'(a,1x,2(i7,1x),2x,3(f6.1,1x))') 'FF333333333333', &
                                       i,iprox,gfs(iprox),modis(i),modis(iprox)
                 return
              endif
           enddo
        enddo
     enddo

  endif

  ! -- search in 4th ring
  !
  if(iset == 0) then
     do l=1,npp
        do k=1,npp
           do j=1,npp
              ii=nextps3(l,k,j)
              do jj=1,npp
                 nextps4(l,k,j,jj)=prox(jj,ii)

                 if(testType == 0)  testO=testGfsO(gfs(prox(jj,ii)))
                 if(testType == 1)  testO=testGfsL(gfs(prox(jj,ii)))
                 if(testType == 10) testO=testModisO(modis(prox(jj,ii)))
                 if(testType == 11) testO=testModisL(modis(prox(jj,ii)))

                 if(testO == 1) then
                    iset=1
                    iprox=prox(jj,ii)
                    if(testType == 0) istate=-4
                    if(testType == 1) istate=4 
                    if(testType == 10) istate=-14
                    if(testType == 11) istate=14
                    if(verb.eq.1) write(*,'(a,1x,2(i7,1x),2x,3(f6.1,1x))') &
                      'FF4444444444444444',i,iprox,gfs(iprox),modis(i),modis(iprox)
                    return
                 endif
              enddo
           enddo
        enddo
     enddo

  endif

  ! -- 5th ring -- bail if get this far...
  !
  if(iset == 0) then
     do m=1,npp
        do l=1,npp
           do k=1,npp
              do j=1,npp
                 ii=nextps4(m,l,k,j)
                 if(verb .ge. 2) print*,'444 ii',ii,lat(ii)*rad2deg,lon(ii)*rad2deg,modis(ii)
                 do jj=1,npp

                    if(testType == 0)  testO=testGfsO(gfs(prox(jj,ii)))
                    if(testType == 1)  testO=testGfsL(gfs(prox(jj,ii)))
                    if(testType == 10) testO=testModisO(modis(prox(jj,ii)))
                    if(testType == 11) testO=testModisL(modis(prox(jj,ii)))

                    if(testO == 1) then
                       iset=1
                       iprox=prox(jj,ii)
                       if(testType == 0)  istate=-5
                       if(testType == 1)  istate=5
                       if(testType == 10) istate=-15
                       if(testType == 11) istate=15
                       if(verb.eq.1) write(*,'(a,1x,2(i7,1x),2x,3(f6.1,1x))') 'FF55555555555555555555', &
                                             i,iprox,gfs(iprox),modis(i),modis(iprox)
                       return
                    endif
                 enddo
              enddo
           enddo
        enddo
     enddo
  endif

  ! -- pathological point -- can't find type after 4th ring... optional code to check 5th ring
  !
  istate=-999
  if(verb.eq.1) print*, 'BBB bailing after 5th ring........',i,gfs(i),modis(i)
  return


  if(iset == 0) then
     print*,'MMM error need to expand search...one more time..'
     print*,'i= ',i,rlat,rlon
  endif

999 continue
  return

end subroutine scanProx

subroutine deconflictGfsland(gfs,work,setO2L,setL2O,nip,nset,no2l,nl2o)

  implicit none

  integer,intent(IN)   :: nip,nset,no2l,nl2o,setO2L(nset,2),setL2O(nset,2)
  real  ,intent(INOUT) :: gfs(nip)
  real  ,intent(OUT)   :: work(nip)
  integer  :: i,iin,iout

  ! -- store the input array to a work array so we do not overwrite changed values
  !
  work=gfs

  ! -- change water to land
  !
  do i=1,no2l
     iin=setO2L(i,1)
     iout=setO2L(i,2)
     gfs(iin)=work(iout)
  enddo

  ! -- change land to water
  !
  do i=1,nl2o
     iin=setL2O(i,1)
     iout=setL2O(i,2)
     gfs(iin)=work(iout)
  enddo

end subroutine deconflictGfsland


function testModisO(r1)
  integer :: testModisO
  REAL, INTENT(IN) :: r1
  testModisO=0
  if(r1 == 17.0 .or. r1 == 21.0) testModisO=1
END FUNCTION testModisO

function testModisL(r1)
  integer :: testModisL
  REAL, INTENT(IN) :: r1
  testModisL=0
  if(r1 /= 17.0 .and. r1 /= 21.0) testModisL=1
END FUNCTION testModisL

function testGfsO(r1)
  integer :: testGfsO
  REAL, INTENT(IN) :: r1
  testGfsO=0
  if(r1 == 0.0) testGfsO=1
END FUNCTION testGfsO

function testGfsL(r1)
  integer :: testGfsL
  REAL, INTENT(IN) :: r1
  testGfsL=0
  if(r1 /= 0.0) testGfsL=1
END FUNCTION testGfsL
end module mdul_mkmodisland
