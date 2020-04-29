module hycom_convec
contains
!*********************************************************************
!  convec
!    Performs convective adjustment in ocean model
!    R. Bleck        April 2012
! 
!*********************************************************************

  subroutine convec(nstep,leapn,uvel,vvel,dens,dp,pres,temp,saln,passv_tr)

  use module_control  ,only: nip
  use module_constants,only: perm
  use fimnamelist     ,only: kdm,itest,diag_intvl
  use hycom_control   ,only: numtr
  use hycom_constants ,only: wet,epsil,onem,onecm,onemm
  use hycom_sigetc

  implicit none
  integer,intent(IN) :: nstep                   ! model time step
  integer,intent(IN) :: leapn                   ! leapfrog time slot (n='new')
!SMS$DISTRIBUTE (dh,2) BEGIN
  real,intent(INOUT) :: uvel (kdm,nip,2)	! u velocity
  real,intent(INOUT) :: vvel (kdm,nip,2)	! v velocity
  real,intent(INOUT) :: dens (kdm,nip  )	! density
  real,intent(INOUT) :: dp   (kdm,nip,2)	! layer thickness
  real,intent(INOUT) :: pres (kdm+1,nip)	! interface pressure
  real,intent(INOUT) :: temp (kdm,nip)		! temperature
  real,intent(INOUT) :: saln (kdm,nip)		! salinity
  real,intent(INOUT),optional :: passv_tr(kdm,nip,numtr)   ! passive tracer
!SMS$DISTRIBUTE END

  real*8 tem,sal,uuu,vvv
  real q1,q2,sigup,siglo,thet,trc(numtr),homog,dpcol(kdm),dncol(kdm),	&
       ucol(kdm),vcol(kdm),tcol(kdm),scol(kdm),pcol(kdm+1),		&
       trcol(kdm,numtr),						&
       totem,tosal,tndcyt,tndcys		!  col.integrals (diag.use only)
  integer i,k,kn,kp,kbase,kmax
  logical :: event,vrbos
  real,parameter :: acurcy = 2.e-6		! 32-bit 

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (vrbos,pcol,ucol,vcol,tcol,scol,dncol,	&
!$OMP trcol,dpcol,totem,tosal,event,kbase,homog,kmax,sigup,siglo,	&
!$OMP q1,q2,uuu,vvv,tem,sal,thet,trc,tndcyt,tndcys)
  do i=1,nip
   if (wet(i) > 0 ) then
    vrbos=i.eq.itest .and. mod(nstep,diag_intvl).eq.0
    if (vrbos) then
     write (*,103) nstep,perm(i),					&
      '  enter convec: temp   saln   dens  thkns   dpth  u (cm/s)  v',	&
      (k,temp(k,i),saln(k,i),dens(k,i),dp(k,i,leapn)/onem,		&
      pres(k+1,i)/onem,uvel(k,i,leapn),vvel(k,i,leapn),k=1,kdm)
103  format (i8,i8,a/(i29,0p,4f7.2,f7.1,2p,2f7.2))
    end if

! --- extract single column from 3-d mesh

    pcol(1)=pres(1,i)
    do k=1,kdm
     ucol(k)=uvel(k,i,leapn)
     vcol(k)=vvel(k,i,leapn)
     tcol(k)=temp(k,i)
     scol(k)=saln(k,i)
     dncol(k)=dens(k,i)
     if (present(passv_tr)) trcol(k,:)=passv_tr(k,i,:)
     dpcol(k)=dp(k,i,leapn)
     pcol(k+1)=pcol(k)+dpcol(k)
    end do

    totem=0.
    tosal=0.
    do k=1,kdm
     totem=totem+tcol(k)*dpcol(k)
     tosal=tosal+scol(k)*dpcol(k)
    end do

! --- convect variables

 9  event=.false.
    do kbase=1,kdm-1
     homog=dpcol(kbase)
     kmax=kbase
     do k=kbase+1,kdm
      if (dncol(k).ge.dncol(kbase)) exit
      kmax=k
      event=.true.
      sigup=dncol(kbase)
      siglo=dncol(k    )
      q1=max(homog,epsil)
      q2=max(dpcol(k),0.)
      uuu=(q1*ucol(kbase)+q2*ucol(k))/(q1+q2)
      vvv=(q1*vcol(kbase)+q2*vcol(k))/(q1+q2)
      tem=(q1*tcol(kbase)+q2*tcol(k))/(q1+q2)
      sal=(q1*scol(kbase)+q2*scol(k))/(q1+q2)
      thet=sigocn(tem,sal)
      homog=homog+dpcol(k)
      ucol(kbase)=uuu
      vcol(kbase)=vvv
      tcol(kbase)=tem
      scol(kbase)=sal
      dncol(kbase)=thet
      if (present(passv_tr)) then
       trc(:)=(q1*trcol(kbase,:)+q2*trcol(k,:))/(q1+q2)
       trcol(kbase,:)=trc(:)
      end if

      if (vrbos) then
       write (*,100) nstep,perm(i),kbase,k,'    upr,lwr,final dncol:',	&
        sigup,siglo,thet,q2/(q1+q2)
100    format (2i8,2i3,a,3f8.3,f6.2)
      end if
     end do

     if (kmax.gt.kbase) then
      if (vrbos) then
       write (*,'(2i8,2(a,i2))') nstep,perm(i),			&
        ' homogenizing layers k= ',kbase,'...',kmax
      end if

      do kp=kbase+1,kmax
       ucol(kp)=uuu
       vcol(kp)=vvv
       tcol(kp)=tem
       scol(kp)=sal
       dncol(kp)=thet
       if (present(passv_tr)) trcol(kp,:)=trc(:)
      end do
     end if
    end do
    if (event) go to 9

    tndcyt=-totem
    tndcys=-tosal
    totem=10.*pres(kdm+1,i)
    tosal=35.*pres(kdm+1,i)
    do k=1,kdm
     tndcyt=tndcyt+tcol(k)*dpcol(k)
     tndcys=tndcys+scol(k)*dpcol(k)
    end do
    if (abs(tndcyt).gt.acurcy*totem) 					&
      write (*,'(i8,a,2es16.8,es9.1)') perm(i),				&
      '  convec - bad temp.intgl.',totem,tndcyt,tndcyt/totem
    if (abs(tndcys).gt.acurcy*tosal) 					&
      write (*,'(i8,a,2es16.8,es9.1)') perm(i),				&
      '  convec - bad saln.intgl.',tosal,tndcys,tndcys/tosal

! --- put 1-d column back into 3-d grid

    do k=1,kdm
     uvel(k,i,leapn)=ucol(k)
     vvel(k,i,leapn)=vcol(k)
     temp(k,i)=tcol(k)
     saln(k,i)=scol(k)
     dens(k,i)=dncol(k)
     if (present(passv_tr)) passv_tr(i,k,:)=trcol(k,:)
    end do

    if (vrbos) then
     write (*,103) nstep,perm(i),					&
      '   exit convec: temp   saln   dens  thkns   dpth  u (cm/s)  v',	&
      (k,temp(k,i),saln(k,i),dens(k,i),dp(k,i,leapn)/onem,		&
      pres(k+1,i)/onem,uvel(k,i,leapn),vvel(k,i,leapn),k=1,kdm)
    end if
   end if			! ocean point
  end do			! horiz. loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END

  return
  end subroutine convec
end module hycom_convec
