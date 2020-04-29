module hycom_hystat
use stencilprint
use findmaxmin1
use findmaxmin2
contains

!*********************************************************************
!    hystat
!	Hydrostatic equation for HYCOM
!       R. bleck		December 2009
!*********************************************************************

  subroutine hystat (nstep,	&
    leap,geop,montg,dens,	& ! geopot., montg.pot., density
    spcifv,temp,saln,		& ! specif.volume, temperature, salinity
    dp,pres,ptrop,		& ! lyr.thknss, intfc pressure, barotp.pres.
    psikk,spvkk,pbot,gzbot)	  ! montg.pot, specif.vol, pres, geopot. at bottom

  use module_control  ,only: nip
  use module_constants,only: grvity,perm
  use fimnamelist     ,only: kdm,itest,diag_intvl,janjic_ocn
  use hycom_constants ,only: wet,thref,onem,r_onem,land_spval
  use hycom_control   ,only: kgap,modesplit

  implicit none

! Dimension and type external variables:
  integer,intent (IN) :: nstep			! model time step
  integer,intent (IN) :: leap			! leap frog time slot
!SMS$DISTRIBUTE (dh,1) BEGIN
  real,intent (IN)    :: pbot(nip)		! initial bottom pressure
  real,intent (IN)    :: gzbot(nip)		! bottom geopotential
  real,intent (IN)    :: psikk(nip),spvkk(nip)	! bottom montg.pot, specif.vol
! Declare local variables:
  real       :: oneta(nip),work(nip)
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE (dh,2) BEGIN
  real,intent (INOUT) :: montg (kdm,nip)	! montgomery potential
  real,intent (INOUT) :: geop(kdm+1,nip)	! geopotential
  real,intent (INOUT) :: pres(kdm+1,nip)	! interface pressure
  real,intent (INOUT) :: spcifv(kdm,nip)	! specific vol. (inv.density)
  real,intent (IN)    :: dens  (kdm,nip)	! potential density
  real,intent (IN)    ::  temp (kdm,nip)	! temperature
  real,intent (IN)    :: saln  (kdm,nip)	! salinity
  real,intent (IN)    :: dp    (kdm,nip,2)	! layer thickness
  real,intent (IN)    :: ptrop(3,nip)		! barotropic pressure
!SMS$DISTRIBUTE END
  real       :: q,thstar
  integer    :: i,k
  character  :: string*2,text*32
  logical    :: vrbos	! switch for 'verbose' mode

  write (text,'(a,i8)') '(ocn_hystat)  step',nstep

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (vrbos,thstar,q)
  do i=1,nip		!  global icos loop
   work(i)=land_spval
   if (wet(i) > 0 ) then
    montg(:,i)=0.
    if (modesplit) oneta(i)=0.
    vrbos=i.eq.itest .and. mod(nstep,diag_intvl).eq.0

    do k=1,kdm
     pres(k+1,i)=pres(k,i)+dp(k,i,leap)

! --- use upper interface pressure in converting sigma to sigma-star.
! --- this is to avoid density variations in layers intersected by sea floor

     thstar=dens(k,i)
!!     +kappaf(temp(k,i),saln(k,i),pres(k,i),wgtkap(i))
     q=-thref*thstar				! - dens anomal / ref density
     spcifv(k,i)=thref*q*(1.+q*(1.+q))		! specific volume anomaly
    end do

    if (modesplit) then
     if (janjic_ocn) then
      print *,'janjic scheme not allowed when modesplit=true'
      stop
     end if

! --- store (1+eta) (= p_total/p_prime) in -oneta-
     oneta(i)=1.+ptrop(leap,i)/pres(kdm+1,i)
 
! --- m_prime in lowest layer:
     montg(kdm,i)=psikk(i)+pres(kdm+1,i)*(spcifv(kdm,i)-spvkk(i))	&
                 +ptrop(leap,i)*spcifv(kdm,i)

! --- m_prime in remaining layers:
     do k=kdm-1,1,-1
      montg(k,i)=montg(k+1,i)+pres(k+1,i)*oneta(i)			&
                *(spcifv(k,i)-spcifv(k+1,i))
     end do

    else				! no modesplit
! --- m_total in lowest layer:
     montg(kdm,i)=psikk(i)+(pres(kdm+1,i)-pbot(i))*(thref+spvkk(i))	&
                 +pres(kdm+1,i)*(spcifv(kdm,i)-spvkk(i))
     work(i)=spcifv(kdm,i)-spvkk(i)

! --- m_total in remaining layers:
     do k=kdm-1,1,-1
      montg(k,i)=montg(k+1,i)+pres(k+1,i)*(spcifv(k,i)-spcifv(k+1,i))
     end do

!    if (janjic_ocn) then
      geop(kdm+1,i)=gzbot(i)+(pres(kdm+1,i)-pbot(i))*thref
      do k=kdm,1,-1
! --- omit contribution of reference density to hydrostatic equation
! --- (i.e. the dynamically irrelevant purely pressure-dependent part)
       geop(k,i)=geop(k+1,i)+spcifv(k,i)*(pres(k+1,i)-pres(k,i))
      end do
     end if
!   end if				! modesplit

    if (vrbos) then
!    do k=1,kdm
!     print '(i7,i3,a,4f11.2)',perm(i),k,' uppr/lowr geopot:',		&
!       montg(k,i)-pres(k  ,i)*(thref+spcifv(k,i)),			&
!       geop(k  ,i)-pres(k  ,i)*thref,					&
!       montg(k,i)-pres(k+1,i)*(thref+spcifv(k,i)),			&
!       geop(k+1,i)-pres(k+1,i)*thref
!    end do
     write (*,103) nstep,perm(i),					&
      ' HYSTAT  temp    saln  dnsity   thkns    dpth   montg   geopot',	&
       (k,temp(k,i),saln(k,i),dens(k,i),				&
        dp(k,i,leap)*r_onem,pres(k+1,i)*r_onem,				&
         montg(k,i)/9.8,geop(k,i)/9.8,k=1,kdm)
103  format (i9,i7,a/(i21,3f8.2,2f8.1,f8.3,f8.2))
    end if

   end if			! ocean point
  end do			! horizontal loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END

  if (mod(nstep,diag_intvl).eq.0) then
   do k=1,kdm,kgap
    write (string,'(i2)') k
    call findmxmn2(montg,kdm,nip,k,'(ohystat) montg '//string,wet)
   end do
   if (modesplit) then
    call findmxmn1(oneta,nip,'oneta',wet)
    if (mod(nstep,diag_intvl).eq.0)				&
     call stencl(oneta,1,1000.,trim(text)//'  oneta x 1000')
   end if
  end if

  if (mod(nstep,diag_intvl).eq.0) then
   call stencl(geop,kdm+1,1.,trim(text)//'  geopotential')
   call stencl(montg,kdm,1.,trim(text)//'  Montgomery potential')
   call stencl(spcifv,kdm,10./thref**2,trim(text)//'  specific volume x 10')
   call stencl(work,1,100./thref**2,trim(text)//'  botm spcif vol diff x 100')
   work(:)=geop(1,:)/grvity
   call stencl(work,1,100.,trim(text)//'  sea surface height (cm)')
  end if

  return
  end subroutine hystat
end module hycom_hystat
