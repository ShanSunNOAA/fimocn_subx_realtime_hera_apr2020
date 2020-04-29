module hycom_diamix
contains
!*********************************************************************
!  diamix
!    diapycnal column mixing using McDougall-Dewar scheme
!
!    R. Bleck   				June 2011
!*********************************************************************

  subroutine diamix(nstep,leapn,dp,dens,pres,temp,saln,uvel,vvel)

  use module_control  ,only: nip
  use module_constants,only: perm
  use fimnamelist     ,only: kdm,itest,diag_intvl,diapyc,diapyn
  use hycom_control   ,only: bclin_frq
  use hycom_constants ,only: wet,epsil,onem,tencm,onecm,onemm,batrop
  use hycom_sigetc
  use hycom_diapfl

  implicit none
  integer,intent(IN) :: nstep			! model time step
  integer,intent(IN) :: leapn			! leapfrog time slot (n='new')
!SMS$DISTRIBUTE  (dh,2) BEGIN
  real,intent(INOUT) :: uvel (kdm,nip,2)	! u velocity
  real,intent(INOUT) :: vvel (kdm,nip,2)	! v velocity
  real,intent(INOUT) :: dens (kdm,nip  )	! density
  real,intent(INOUT) :: dp   (kdm,nip,2)	! layer thickness
  real,intent(INOUT) :: pres (kdm+1,nip)	! interface pressure
  real,intent(INOUT) :: temp (kdm,nip)		! temperature
  real,intent(INOUT) :: saln (kdm,nip)		! salinity
!SMS$DISTRIBUTE END
  logical vrbos,event
  integer i,k,kk,k1,kp,kmax,it,iter,klist(kdm),kplist(kdm)
  real delt
  real*8 dpcol(kdm),pcol(kdm+1),tcol(kdm),scol(kdm),rcol(kdm),		&
       ucol(kdm),vcol(kdm),dpsub(kdm),tsub(kdm),ssub(kdm),rsub(kdm),	&
       pold(kdm+1),colint,colins,cloutt,clouts,colinu,cloutu,		&
       colinv,cloutv,usum,vsum,delp
  real   ,parameter :: acurcy=1.e-11
  character string*4
  integer,parameter :: mlskip=1		! 0/1: include/exclude mixed layer

  if (diapyc.eq.0. .and. diapyn.eq.0.) return
  kk=kdm-mlskip

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (vrbos,tcol,scol,rcol,dpcol,ucol,vcol,	&
!$OMP pcol,pold,colint,colins,colinu,colinv,kmax,klist,iter,event,	&
!$OMP k1,delt,it,kp,dpsub,tsub,ssub,rsub,kplist,cloutu,cloutv,		&
!$OMP usum,vsum,delp,cloutt,clouts)
  do i=1,nip
  if (wet(i) > 0 ) then
   vrbos=i.eq.itest .and. mod(nstep,diag_intvl).eq.0

   if (vrbos) then
    print 104,nstep,perm(i),						&
      'DIAMIX  in:  thkns    temp    saln    dens      u       v'
    print 105,(k,dp(k,i,leapn)/onem,temp(k,i),saln(k,i),dens(k,i),      &
      uvel(k,i,leapn),vvel(k,i,leapn),k=1,kdm)
104  format (i9,i7,3x,a)
105  format (i28,f9.2,5f8.2)
   end if

   do k=1,kk
    tcol (k)=temp (k+mlskip,i)
    scol (k)=saln (k+mlskip,i)
    rcol (k)=dens (k+mlskip,i)
    dpcol(k)=dp   (k+mlskip,i,leapn)
    ucol (k)=uvel (k+mlskip,i,leapn)
    vcol (k)=vvel (k+mlskip,i,leapn)
   end do

   pcol(1)=pres(1,i)
   pold(1)=pcol(1)
   do k=1,kk
    pcol(k+1)=pcol(k)+dpcol(k)
    pold(k+1)=pcol(k+1)
   end do

! --- column integrals are computed for diagnostic purposes only
   colint=0.
   colins=0.
   colinu=0.
   colinv=0.
   do k=1,kk
    colint=colint+tcol(k)*dpcol(k)
    colins=colins+scol(k)*dpcol(k)
    colinu=colinu+ucol(k)*dpcol(k)
    colinv=colinv+vcol(k)*dpcol(k)
   end do

! --- McD-D cannot inflate multiple near-massless layers sandwiched between
! --- a pair of thicker layers. If -n- such layers (n>1) are present, we must
! --- inflate them one-by-one by solving the diffusion eqn repeatedly using a
! --- vertical profile from which the remaining n-1 layers are excluded.

   kmax=0
   do k=1,kk
    if (pcol(k).gt.pcol(kk+1)-onemm) exit
    kmax=k
    if (dpcol(k).gt.tencm) then
     klist(k)=-1		! mass-containing layer
!   else if (rcol(k).lt.dens(ktfoss,i)) then
!    klist(k)=99		! massless layers beneath ML (non-inflatable)
    else
     klist(k)=0			! inflatable massless layer
    end if
   end do
   if (kmax.eq.0) go to 3	! no action - single layer touching bottom

! --- label massless layers such that in a contiguous sequence of such layers
! --- each one gets a different label. this way, they can be inflated one
! --- after the other.

   iter=0
1  iter=iter+1
   event=.false.
   if (mod(iter,2).eq.1) then
! --- top-down scan
    k1=1
    do k=2,kmax
     if (klist(k).ne.99) then
      if (klist(k).eq.0 .and. klist(k1).ne.iter			&
                        .and. klist(k1).ne.0) then
       klist(k)=iter		! include layer -k- in iteration 'iter'
       event=.true.
      end if
      k1=k
     end if
    end do
   else
! --- bottom-up scan
    k1=kmax
    do k=kmax-1,1,-1
     if (klist(k).ne.99) then
      if (klist(k).eq.0 .and. klist(k1).ne.iter			&
                        .and. klist(k1).ne.0) then
       klist(k)=iter		! include layer -k- in iteration 'iter'
       event=.true.
      end if
      k1=k
     end if
    end do
   end if
   if (event) go to 1

   if (vrbos) then
   print '(i9,i7,a/(20i3))',nstep,perm(i),' (diamix) klist:',	&
     (klist(k),k=1,kmax)
   end if

! --- diffusion eqn will be solved 'iter' times. during the first (iter-1)
! --- iterations, we select a different member from each set of consecutive
! --- massless layers. all layers participate in the final iteration.

   delt=batrop*bclin_frq/iter

   if (kmax.gt.2) then
    do it=1,iter
     kp=0
     do k=1,kmax
      if (klist(k).ne.99) then
       if (it.eq.iter .or. klist(k).lt.0 .or. klist(k).eq.it) then
        kp=kp+1
        dpsub(kp)=dpcol(k)
        tsub (kp)=tcol (k)
        ssub (kp)=scol (k)
        rsub (kp)=rcol (k)
        kplist(kp)=k
       end if
      end if
     end do

! --- now do the actual mixing
     if (vrbos) then
      print '(i9,i7,a,i3,a/(20i3))',nstep,perm(i),' pass',it,		&
       '  calling diapfl with the following layers:',(kplist(k),k=1,kp)
     end if

     call diapfl(nstep,delt,kp,dpsub,rsub,tsub,ssub,vrbos,perm(i))

     do k1=1,kp
      k=kplist(k1)
      dpcol(k)=dpsub(k1)
      rcol (k)=rsub (k1)
      tcol (k)=tsub (k1)
      scol (k)=ssub (k1)
     end do

     do k=1,kk
      pcol(k+1)=pcol(k)+dpcol(k)
      if (pcol(k).gt.pcol(kk+1)-onemm) exit
      kmax=k
     end do
     if (kmax.le.2) exit
    end do			! iter
   end if			! kmax > 2

! --- momentum diffusion

   do k=1,kk+1
    pcol(k)=min(pold(kk+1),pcol(k))	! suppress roundoff errors in bottm prs
   end do
   pcol(kk+1)=pold(kk+1)		! suppress roundoff errors in bottm prs
   cloutu=0.
   cloutv=0.
   kp=1
   do k=1,kk
    if (dpcol(k).gt.0.) then
     usum=0.
     vsum=0.
2    delp=min(pcol(k+1),pold(kp+1))-max(pcol(k),pold(kp))
     usum=usum+ucol(kp)*delp
     vsum=vsum+vcol(kp)*delp
     if (pold(kp+1).lt.pcol(k+1)) then
      kp=kp+1
      go to 2
     end if
     uvel(k,i,leapn)=usum/dpcol(k)
     vvel(k,i,leapn)=vsum/dpcol(k)
     cloutu=cloutu+usum
     cloutv=cloutv+vsum
    end if
   end do

! --- check conservation of column integrals
   cloutt=0.
   clouts=0.
   do k=1,kk
    cloutt=cloutt+tcol(k)*dpcol(k)
    clouts=clouts+scol(k)*dpcol(k)
   end do
   if (abs(colint-cloutt).gt.acurcy*10.*pres(kdm+1,i)) then
    print 103,perm(i),'(diamix) T column intgl.error',			&
      colint,cloutt,(cloutt-colint)/(10.*pres(kdm+1,i))
   end if
   if (abs(colins-clouts).gt.acurcy*35.*pres(kdm+1,i)) then
    print 103,perm(i),'(diamix) S column intgl.error',			&
      colins,clouts,(clouts-colins)/(35.*pres(kdm+1,i))
   end if
   if (abs(colinu-cloutu).gt.acurcy*.1*pres(kdm+1,i)) then
    print 103,perm(i),'(diamix) U column intgl.error',			&
      colinu,cloutu,(cloutu-colinu)/(.1*pres(kdm+1,i))
   end if
   if (abs(colinv-cloutv).gt.acurcy*.1*pres(kdm+1,i)) then
    print 103,perm(i),'(diamix) V column intgl.error',			&
      colinv,cloutv,(cloutv-colinv)/(.1*pres(kdm+1,i))
   end if
103 format (i7,3x,a,2es14.6,es9.1)

   do k=1,kk
    temp (k+mlskip,i)      =tcol (k)
    saln (k+mlskip,i)      =scol (k)
    dens (k+mlskip,i)      =rcol (k)
    dp   (k+mlskip,i,leapn)=dpcol(k)
    uvel (k+mlskip,i,leapn)=ucol (k)
    vvel (k+mlskip,i,leapn)=vcol (k)
   end do
   do k=1,kdm
    pres(k+1,i)=pres(k,i)+dp(k,i,leapn)
   end do
 3 continue

   if (vrbos) then
    print 104,nstep,perm(i),						&
      'DIAMIX out:  thkns    temp    saln    dens      u       v'
    print 105,(k,dp(k,i,leapn)/onem,temp(k,i),saln(k,i),		&
      dens(k,i),uvel(k,i,leapn),vvel(k,i,leapn),k=1,kdm)
   end if

  end if			! ocean point
  end do			! horiz. loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END

  return
  end subroutine diamix
end module hycom_diamix
