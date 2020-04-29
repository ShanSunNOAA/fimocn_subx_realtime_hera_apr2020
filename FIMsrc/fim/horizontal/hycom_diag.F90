module hycom_diag
contains

   subroutine transec(nstep,leap)

   use module_control  ,only: nip,dt
   use fimnamelist     ,only: ArchvTimeUnit, ocn_trnsecdir
   use hycom_constants ,only: wet,ArchvStepOcn
   use hycom_variables ,only: curl_ave,mssfx_ave
   use hycom_sverbal   ,only: sverbal,ocn_flxsum
   implicit none

   integer,intent(IN) :: leap		! leapfrog time slot
   integer,intent(IN) :: nstep		! (atmo) time step

   real,external :: its2time
   integer       :: i
   real          :: time

   time=its2time(nstep)
! --- perform various online diagnostics

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO
   do i=1,nip
    if (wet(i) > 0 ) then
     curl_ave (    i)=curl_ave (    i)/float(ArchvStepOcn)
     mssfx_ave(:,:,i)=mssfx_ave(:,:,i)/float(ArchvStepOcn)
    end if
   end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END

   if (len_trim(ocn_trnsecdir).gt.0) then
! --- diagnose sverdrup transport at specific latitudes
    call sverbal(nstep,curl_ave)

! --- diagnose western bdry current transport (for comparison with sverdrup
! --- transport) and transport through various passages
    call ocn_flxsum(nstep,mssfx_ave)
   end if				! ocn_trnsecdir exists
!
   return
   end subroutine transec


   subroutine glbdia(nstep,leap,what)

   use module_control  ,only: nip,dt
   use module_constants,only: area,rarea,grvity,deg_lat,deg_lon
   use fimnamelist     ,only: ArchvTimeUnit,kdm,diag_intvl,coupled,	&
                              ocn_trnsecdir,do_radcor,do_pcpcor
   use hycom_control   ,only: bclin_frq
   use hycom_constants ,only: wet,epsil,batrop,ArchvStepOcn,spcifh,	&
       onem,rhoice,stdsal,saldif,thref,grvity,fusion
   use hycom_variables ,only: srfx_ave,pmne_ave,montg,dp,temp,saln,	&
       temglb0,salglb0,rhoglb0,curl_ave,mssfx_ave,thkice,covice,dens,	&
       srfxcum,pmnecum,massglb0,totmass,ocnarea,radcor,pcpcor,radtot,	&
       pcptot
   implicit none

   integer,intent(IN) :: leap		! leapfrog time slot
   integer,intent(IN) :: nstep		! (atmo) time step
   character(len=*),intent(IN),optional :: what

   real,external :: its2time
!  real*8   :: rhocol,temcol,salcol
!  integer  :: i,k
   integer  :: i
   real     :: time,elapsd
   real*8   :: tmean,smean,rmean,ttl_nh_icevol,ttl_sh_icevol,		&
               ttl_nh_icecov,ttl_sh_icecov,totice,glbssh,thrmfx,salnfx
   character(len=80) :: string
!  real,parameter :: rst_days = 3650.	! days
   real,parameter :: rst_days = 365.	! days

!SMS$DISTRIBUTE (dh,1) BEGIN
  real icevol_nh(nip),icevol_sh(nip),icecov_nh(nip),icecov_sh(nip)
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE (dh,2) BEGIN
  real work1(kdm,nip),work2(kdm,nip)
!SMS$DISTRIBUTE END
 
   time=its2time(nstep)
!
! --- compare area-integrated surface fluxes to volume-integrated tendencies
!
   if (mod(nstep,diag_intvl).ne.0) return
     icecov_nh(:)=0.
     icevol_nh(:)=0.
     icecov_sh(:)=0.
     icevol_sh(:)=0.
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO
  do i=1,nip
   if (deg_lat(i) > 0.) then
     icecov_nh(i)=covice(i)
     icevol_nh(i)=thkice(i)
   else
     icecov_sh(i)=covice(i)
     icevol_sh(i)=thkice(i)
   end if
  end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END

   call glob2d(icecov_nh,ttl_nh_icecov)
   call glob2d(icecov_sh,ttl_sh_icecov)
   call glob2d(icevol_nh,ttl_nh_icevol)
   call glob2d(icevol_sh,ttl_sh_icevol)
   call glob2d(montg(1,:),glbssh)

   call glob2d(thkice,totice)
   call glob2d(srfxcum,thrmfx)		! W
   call glob2d(pmnecum,salnfx)

   work2(:,:)=dp(:,:,leap)
   call glob3d(temp,work2,tmean)	! deg kg m/s^2
   call glob3d(saln,work2,smean)	! g m/s^2
   call glob3d(dens,work2,rmean)
   work1(:,:)=1.
   call glob3d(work2,work1,totmass)	! N
   totmass=totmass/grvity		! N => kg

   totice=totice*rhoice			! kg
   tmean=tmean/(totmass*grvity)-totice*fusion/(spcifh*totmass)		! deg
   smean=smean/(totmass*grvity)-totice*saldif/totmass			! g/kg
   rmean=rmean/(totmass*grvity)
   salnfx=-salnfx*stdsal/thref		! g/sec

   if (present(what)) then
     string=trim(what)
   else
     string='(glbdia)'
   endif

   if (nstep==0) then                 		! save initial global T/S etc.
     temglb0=tmean
     salglb0=smean
     rhoglb0=rmean
     massglb0=totmass
     print '(a,3f9.4,es11.4)','(glbdia) initl global T,S,sig,mass:',	&
       temglb0,salglb0,rhoglb0,massglb0
     print '(a,2f9.3)',' initl ice extent (Mm^2) in NH/SH=',		&
       ttl_nh_icecov*1.e-12,ttl_sh_icecov*1.e-12
     print '(a,2f9.3)',' initl ice thickness (m) in NH/SH=',		&
       ttl_nh_icevol/max(epsil,ttl_nh_icecov),				&
       ttl_sh_icevol/max(epsil,ttl_sh_icecov)
   else				! compare flux-inferred and actual tendencies
    elapsd=nstep*batrop		! time (sec) elapsed since model start

!   if (mod(elapsd+30.,365.*86400.).le.60.) then		! once a year
     if (coupled) then
      radtot=374.8			! incoming SW in W/m2 from CORE
      pcptot=1.1185/(365.*86400.)	! m/sec from CORE
     end if

     radcor=-(tmean-temglb0)*totmass*spcifh/		&
       (       ocnarea*rst_days*86400.*radtot)
     pcpcor= (smean-salglb0)*totmass*thref/		&
       (stdsal*ocnarea*rst_days*86400.*pcptot)
     write(*,'(a,1x,a,f8.2,a,2f9.3,a,2f9.5)') trim(string),		&
       ArchvTimeUnit,time,' heatflux & freshwater adjustment factor=',	&
         1.+radcor,1.+pcpcor,' T/S gain=',tmean-temglb0,smean-salglb0
    end if

    write(*,104) trim(string),ArchvTimeUnit,time,			&
      '  glob T,S,sig,mass',tmean,smean,rmean,totmass,			&
      '             initl:',temglb0,salglb0,rhoglb0,massglb0
    write(*,103) trim(string),ArchvTimeUnit,time,			&
      '  globl T,S,mass drift',leap,					&
      tmean-temglb0,smean-salglb0,totmass-massglb0
    write(*,103) trim(string),ArchvTimeUnit,time,			&
      '  srf.flx-induced drift (deg,psu)',leap,				&
      thrmfx*batrop/(spcifh*totmass),					&
      salnfx*batrop/(       totmass)

    write(*,102) trim(string),ArchvTimeUnit,time,			&
      '  ice extent (Mm^2) in NH/SH=',ttl_nh_icecov*1.e-12,		&
      ttl_sh_icecov*1.e-12
    write(*,102) trim(string),ArchvTimeUnit,time,			&
      '  ice thickness (m) in NH/SH=',					&
      ttl_nh_icevol/max(epsil,ttl_nh_icecov),				&
      ttl_sh_icevol/max(epsil,ttl_sh_icecov)
    write(*,102) trim(string),ArchvTimeUnit,time,			&
      '  ssh (m) =',glbssh/ocnarea

!   print 100,ArchvTimeUnit,time,					&
!      '  extrapolated temp drift (deg/century):'			&
!      ,watcumu/(spcifh*totmass) * 36500.*86400./elapsd
!   print 100,ArchvTimeUnit,time,					&
!      '  extrapolated saln drift (psu/century):'			&
!      ,-pmecumu*stdsal/(thref*totmass) * 36500.*86400./elapsd
!  end if	! once a year

   return
100 format ("(hycom_glbdia) ",a,f8.1,a,3f9.4)
101 format ("(hycom_glbdia) ",a,f8.1,a,3es11.3)
102 format (a,1x,a,f8.1,a,2f9.3)
103 format (a,1x,a,f8.1,a,i2,3es11.3)
104 format (a,1x,a,f8.1,a,f8.4,2f9.4,es11.4/23x,a,f8.4,2f9.4,es11.4)
  end subroutine glbdia


  subroutine glob3d(field,dp,totout)

! --- compute global dp-weighted integral of 'field'

  use module_control  ,only: nip
  use module_constants,only: area,grvity,inv_perm
  use fimnamelist     ,only: kdm
  use hycom_constants ,only: wet
  use hycom_variables ,only: totmass
  implicit none
!SMS$DISTRIBUTE (dh,2) BEGIN
  real  ,intent(IN)  :: field(kdm,nip),dp(kdm,nip)
!SMS$DISTRIBUTE END
  real*8,intent(OUT) :: totout
  integer i,k,ip
  real colsum

  totout=0.
!ss use SERIAL to guarantee bit-exact results across different MPI task counts
!SMS$SERIAL (<wet,inv_perm,field,dp,area,IN>,<totout,OUT> : default=ignore) BEGIN
  do i=1,nip
    ip=inv_perm(i)
    if (wet(ip) > 0) then
      colsum=0.
      do k=1,kdm
        colsum=colsum+field(k,ip)*dp(k,ip)*area(ip)
      end do
      totout=totout+colsum
    end if
  end do
!SMS$SERIAL END
  return
  end subroutine glob3d


!   real*8 function glob3d(field,dp)
! 
! ! --- compute global average of 'field'
! 
!   use module_control  ,only: nip
!   use module_constants,only: area,grvity
!   use fimnamelist     ,only: kdm
!   use hycom_constants ,only: wet
!   use hycom_variables ,only: totmass
!   implicit none
! !SMS$DISTRIBUTE (dh,2) BEGIN
!   real,intent(IN) :: field(kdm,nip),dp(kdm,nip)
! !SMS$DISTRIBUTE END
!   integer i,k
!   real*8 totl,colsum
! 
!   totl=0.
! !SMS$PARALLEL (dh,i) BEGIN
!   do i=1,nip
!    if (wet(i) > 0) then
!     colsum=0.
!     do k=1,kdm
!      colsum=colsum+field(k,i)*dp(k,i)*area(i)
!     end do
!     totl=totl+colsum
!    end if
!   end do
! !SMS$PARALLEL END
! !SMS$REDUCE(totl,SUM)
!   glob3d=totl/(grvity*totmass)
!   return
!   end function glob3d


  subroutine glob2d(field,totout)

! --- compute global surface integral of 'field'

  use module_control  ,only: nip
  use module_constants,only: area,inv_perm
  use hycom_constants ,only: wet
  implicit none
!SMS$DISTRIBUTE (dh,1) BEGIN
  real  ,intent(IN)  :: field(nip)
!SMS$DISTRIBUTE END
  real*8,intent(OUT) :: totout
  integer i,ip

  totout=0.
!ss use SERIAL to guarantee bit-exact results across different MPI task counts
!SMS$SERIAL (<wet,inv_perm,field,area,IN>,<totout,OUT> : default=ignore) BEGIN
  do i=1,nip
    ip=inv_perm(i)
    if (wet(ip) > 0) totout=totout+field(ip)*area(ip)
  end do
!SMS$SERIAL END
  return
  end subroutine glob2d


! ! real*8 function glob2d(field)
!   real   function glob2d(field)
! 
! ! --- compute global surface integral of 'field'
! 
!   use module_control  ,only: nip
!   use module_constants,only: area
!   use hycom_constants ,only: wet
!   implicit none
! !SMS$DISTRIBUTE (dh,1) BEGIN
!   real,intent(IN) :: field(nip)
! !SMS$DISTRIBUTE END
!   integer i
! ! real*8 totl
!   real   totl
! 
!   totl=0.
! !SMS$PARALLEL (dh,i) BEGIN
!   do i=1,nip
!    if (wet(i) > 0) totl=totl+field(i)*area(i)
!   end do
! !SMS$PARALLEL END
! !SMS$REDUCE(totl,SUM)
!   glob2d=totl
!   return
!   end function glob2d


   subroutine prt_column(nstep,what)

   use module_control  ,only: nip
   use module_constants,only: deg_lat, deg_lon, nprox, prox, pi
   use hycom_variables ,only: temp,saln,dens,dp
   use fimnamelist     ,only: ArchvTimeUnit,itest,kdm
   use hycom_constants ,only: onem

   implicit none
   character,intent(IN) :: what*(*)
   integer,intent(IN) :: nstep
   real           :: time
   real,external  :: its2time

   integer i,k

#include <gptl.inc>

   if (itest.le.0) return
   time=its2time(nstep)

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO
   do i=1,nip
    if (i.eq.itest) then
    write(*,'(i8,1x,a,f6.1,i7,2a/(29x,i3,3f8.3,f9.3))') nstep,		&
      ArchvTimeUnit,time,i,' entering ',trim(what),(k,temp(k,i),	&
      saln(k,i),dens(k,i),dp(k,i,2)/onem,k=1,kdm)
    end if
   end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END

   call flush(6)
   return
   end subroutine prt_column


  subroutine slabsum(nstep,leap,text)

! --- compute integrals of T/S over sets of layers

  use module_control  ,only: nip
  use module_constants,only: area,grvity
  use fimnamelist     ,only: kdm
  use hycom_constants ,only: wet,rhoice,saldif
  use hycom_variables ,only: dp,temp,saln,thkice
  implicit none
  integer  ,intent(IN) :: nstep,leap		! time step, leapfrog time slot
  character,intent(IN) :: text*(*)
  real*8 wgtglb,temglb,salglb,wgtcol,temcol,salcol
  integer :: i,k,k1,k2,nslab=8

  do k1=1,kdm,nslab
   k2=min(kdm,k1+nslab-1)
   wgtglb=0.
   temglb=0.
   salglb=0.
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE(wgtcol,temcol,salcol) REDUCTION(+:wgtglb,temglb,salglb)
   do i=1,nip
    if (wet(i) > 0) then
     wgtcol=0.
     temcol=0.
     salcol=0.
     if (k1.eq.1) salcol=-saldif*thkice(i)*rhoice*grvity*area(i)
     do k=k1,k2
      wgtcol=wgtcol+dp(k,i,leap)          *area(i)
      temcol=temcol+dp(k,i,leap)*temp(k,i)*area(i)
      salcol=salcol+dp(k,i,leap)*saln(k,i)*area(i)
     end do
     wgtglb=wgtglb+wgtcol
     temglb=temglb+temcol
     salglb=salglb+salcol
    end if
   end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END
!SMS$REDUCE(wgtglb,temglb,salglb,SUM)
   print '(i8,2(a,i2),a,a14,2(3x,2f14.8))',nstep,' lyrs ',k1,'-',k2,	&
     '  ',text,temglb/wgtglb,salglb/wgtglb
  end do			! k1 loop

  call ocnglobsum(nstep,leap,text)

  return
  end subroutine slabsum


  subroutine ocnglobsum(nstep,leap,text)

! --- compute integrals of T/S over sets of layers

  use module_control  ,only: nip
  use module_constants,only: area,grvity
  use fimnamelist     ,only: kdm
  use hycom_constants ,only: wet,rhoice,saldif,onecm
  use hycom_variables ,only: dp,temp,saln,thkice,totmass,ocnarea
  implicit none
  integer  ,intent(IN) :: nstep,leap		! time step, leapfrog time slot
  character,intent(IN) :: text*(*)
  real*8 wgtglb,temglb,salglb,wgtcol,temcol,salcol
  integer :: i,k

  wgtglb=0.
  temglb=0.
  salglb=0.
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE(wgtcol,temcol,salcol) REDUCTION(+:wgtglb,temglb,salglb)
  do i=1,nip
   if (wet(i) > 0) then
    wgtcol=0.
    temcol=0.
    salcol=0.
    salcol=-saldif*thkice(i)*rhoice*grvity*area(i)
    do k=1,kdm
     wgtcol=wgtcol+dp(k,i,leap)          *area(i)
     temcol=temcol+dp(k,i,leap)*temp(k,i)*area(i)
     salcol=salcol+dp(k,i,leap)*saln(k,i)*area(i)
    end do
    wgtglb=wgtglb+wgtcol
    temglb=temglb+temcol
    salglb=salglb+salcol
   end if
  end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END
!SMS$REDUCE(wgtglb,temglb,salglb,SUM)
  print '(i8,a20,3x,2f14.8,3x,a,f6.2,a)',nstep,text,temglb/wgtglb,	&
    salglb/wgtglb,'d(pbot)=',(wgtglb-grvity*totmass)/(onecm*ocnarea),' cm'
  return
  end subroutine ocnglobsum

end module hycom_diag
