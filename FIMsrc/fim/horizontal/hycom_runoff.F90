module hycom_runoff
use findmaxmin1
use stencilprint
contains
!*********************************************************************
!  runoff
!   construct watersheds for precip runoff
!
!   R. Bleck					May  2011 
!   R. Bleck	improved watershed detection	Dec  2011 
!   R. Bleck	generalization for curve > 0	May  2014
!*********************************************************************
  subroutine runoff(topo,nuphill,uphill,dnhill)
  use module_control,  only: nip,npp
  use fimnamelist,     only: glvl
  use module_constants,only: nprox,prox,deg_lat,deg_lon,perm,inv_perm
  use hycom_constants ,only: wet
  use fimnamelist     ,only: itest
! use module_dffusn_lev
! use ersatzpipe,only: pipe2
  implicit none
!SMS$DISTRIBUTE  (dh,1) BEGIN
  real,intent(INOUT)  :: topo(nip)	! surface geopotential
  integer,intent(OUT) :: nuphill(nip)	! no.of cells draining into cell -i-
  integer,intent(OUT) :: dnhill(nip)	! downhill ngbor of -i- (local index)
  integer             :: dnhillp(nip)	! downhill ngbor of -i- (global index)
  integer             :: rivmou(nip)	! ultimate outflow point for cell -i-
  real                :: topo1(nip)	! auxiliary topo array
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE  (dh,2) BEGIN
  integer,intent(OUT) :: uphill(npp,nip)! neighboring cells draining into -i-
  integer :: uphillp(npp,nip)! neighboring cells draining into -i- (glob.index)
!SMS$DISTRIBUTE END
  integer i,ix,edg,iter,nxt,nr,totcst,totinl
  integer,allocatable :: shedsz(:),shedno(:)
  character(len=1) :: char1
  real    hgtmin,dfflen(1),addshd,ncount,ncoast,inland,ntot,nshd,nrivr,rnxt,sum
  integer itestp
!
  print *,'entering runoff ...'

  itestp=0
  if (itest.gt.0) then
!SMS$SERIAL (<perm,IN>) BEGIN
   itestp=perm(itest)			! test point (global index)
!SMS$SERIAL END
  end if
!
! --- optional: smooth terrain height
! call findmxmn1(topo,nip,'original topo')
! call stencl(topo,1,1.,'(runoff) original surf.geopot.')
!
! dfflen(1)=sqrt(5.e14/float(nip))*.2		! avg.mesh size / 5
! do iter=1,1
!  print *,'(runoff) smooth topography, dfflen =',dfflen(1)
!  call dffusn_lev(topo,dfflen,1,1,1,.false.)
! end do

! ---elevate points that otherwise won't drain
 
  iter=0
4 ncount=0.
  iter=iter+1
!SMS$PARALLEL (dh,i) BEGIN
  topo1(:)=topo(:)
!SMS$EXCHANGE (topo1)
!$OMP PARALLEL DO PRIVATE(hgtmin) REDUCTION(+:ncount)
  do i=1,nip
   if (wet(i) == 0) then	! we are on land
    hgtmin=99999.
    do edg=1,nprox(i)
     hgtmin=min(hgtmin,topo1(prox(edg,i)))
    end do
    hgtmin=hgtmin+.1		! maintain nonzero slope
    if (topo1(i).lt.hgtmin) then
!     print '(a,i8,2(a,f10.2))','(runoff) raising terrain at',	&
!      perm(i),'  from',topo1(i),'  to', hgtmin
     topo(i)=hgtmin
     ncount=ncount+1.
    end if
   end if			! land point
  end do			! horiz. loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END
!SMS$REDUCE (ncount,SUM)
  if (mod(iter,100).eq.1) print '(i7,a,i7)',int(ncount),	&
    '  points raised, iteration',iter
  if (ncount.gt.0.) go to 4
  call findmxmn1(topo,nip,'altered topo')
  call stencl(topo,1,1.,'(runoff) altered surf.geopot.')
  call flush(6)

! call pipe2(0,topo,1,'topo')
!sms$compare_var(topo, "hycom_runoff.F90 - topo")

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- step 1: find lowest-elevation neighbor of each grid cell

  ntot=0.
  nshd=0.
!SMS$PARALLEL (dh,i) BEGIN
!SMS$EXCHANGE (topo)
!$OMP PARALLEL DO PRIVATE(hgtmin,ix) REDUCTION(+:ntot,nshd)
  do i=1,nip
   dnhill(i)=0
   dnhillp(i)=0
   if (wet(i) == 0) then	! we are on land
    ntot=ntot+1.
    hgtmin=topo(i)
    do edg=1,nprox(i)		! scan neighbors for lowest elevation
     ix=prox(edg,i)
     if (wet(ix) > 0) then	! neighbor -ix- is coastal outflow point
      dnhill(i)=ix		! point -i- drains into -ix-
      dnhillp(i)=perm(ix)
      if (i.eq.itest)							&
        print 101,'water at',perm(i),' (la/lo=',deg_lat(i),deg_lon(i),	&
          ') flows to',perm(ix),' (la/lo=',deg_lat(ix),deg_lon(ix),')'
        print 100,'point',perm(ix),' is coastal point'
      exit
     else if (topo(ix).le.hgtmin) then
      hgtmin=topo(ix)		! found lower neighbor
      dnhill(i)=ix
      dnhillp(i)=perm(ix)
     end if
    end do			! loop through edges

    ix=dnhill(i)
    if (ix.gt.0) then		! flow is directed from -i- to -ix-
     nshd=nshd+1.
     if (i.eq.itest)							&
      print 101,'water at',perm(i),' (la/lo=',deg_lat(i),deg_lon(i),	&
        ') flows to',perm(ix),' (la/lo=',deg_lat(ix),deg_lon(ix),')'
    else			! no lower neighbor exists
     print '(a,i8,a,2f7.2)','(runoff) no outflow from point',		&
       perm(i),'  la/lo=',deg_lat(i),deg_lon(i)
!    dnhill(i)=i		! water stays where it is
!    dnhillp(i)=i		! water stays where it is
     stop '(runoff error)'
    end if
   end if			! -i- is land point
  end do			! horiz. loop
!$OMP END PARALLEL DO
!SMS$EXCHANGE (dnhillp)
!SMS$PARALLEL END
!SMS$REDUCE (ntot,nshd,SUM)

  print 100,'downhill neighbors found for ',int(nshd),' of',int(ntot),' points'
100 format ('(runoff) ',3(a,i8))
  call flush(6)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- step 2: looking upstream, find all points that drain into a given cell

  nshd=0.
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE(ix) REDUCTION(+:nshd)
  do i=1,nip
   nuphill(i)=0				! # of cells draining into -i-
   uphill(:,i)=0
   do edg=1,nprox(i)
    ix=prox(edg,i)
    if (dnhillp(ix).eq.perm(i)) then	! flow is directed from -ix- to -i-
     nuphill(i)=nuphill(i)+1
     uphill(nuphill(i),i)=ix
     nshd=nshd+1.
    end if
   end do			! loop through edges
  end do			! horiz. loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END
!SMS$REDUCE (nshd,SUM)

  print 100,'precip will be collected from',int(nshd),' of',int(ntot),' points'

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- step 3: construct individual watersheds (for diagnostic use)

! --- 3a: identify coastal outflow points and inland depressions

  nrivr=0.
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE(ix) REDUCTION(+:nrivr)
  do i=1,nip
   rivmou(i)=0
   if (wet(i) > 0) then
    do edg=1,nprox(i)
     ix=prox(edg,i)
     if (dnhillp(ix).eq.perm(i)) then
      rivmou(i)=perm(i)			! -i- is coastal outflow point
      exit
     end if
    end do
    if (rivmou(i).gt.0) then
     nrivr=nrivr+1.
     print 101,'point',perm(i),' (la/lo=',			&
       deg_lat(i),deg_lon(i),') is coastal outflow point'
    end if
   else if (dnhillp(i).eq.perm(i)) then
    rivmou(i)=-perm(i)			! -i- is inland depression
    nrivr=nrivr+1.
    print 101,'point',perm(i),' (la/lo=',			&
      deg_lat(i),deg_lon(i),') is inland depression'
   end if
  end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END
!SMS$REDUCE (nrivr,SUM)
  print *,'total of',int(nrivr),'  watersheds'
  call flush(6)

! --- 3b: generate catalog of all watersheds

  allocate(shedno(int(nrivr)),shedsz(int(nrivr)))
  shedno(:)=0
  shedsz(:)=0
!SMS$SERIAL(<perm,rivmou,IN>,<shedno,OUT>:DEFAULT=IGNORE) BEGIN
  nr=0
  do i=1,nip
   if (rivmou(i).ne.0) then
    nr=nr+1
    shedno(nr)=perm(i)
    shedsz(nr)=1
    print 100,'watershed no.',nr,'  terminates at point',perm(i)
   end if
  end do
!SMS$SERIAL END
  call flush(6)

! --- 3c: identify all points belonging to a given watershed

  totcst=0
  totinl=0
3 addshd=0.
  ncoast=0.
  inland=0.

!SMS$EXCHANGE (rivmou)
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE(ix) REDUCTION(+:addshd,ncoast,inland)
  do i=1,nip
   if (wet(i) == 0) then	! we are on land
    ix=dnhill(i)		! water is directed from -i- to -ix-
    if (rivmou(i).eq.0 .and. rivmou(ix).ne.0) then
     rivmou(i)=rivmou(ix)
     addshd=addshd+1.
     if (rivmou(i).gt.0) then
      ncoast=ncoast+1.
     else
      inland=inland+1.
     end if
    end if			! point -i- added to watershed
   end if			! land point
  end do			! horiz. loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END
!SMS$REDUCE (addshd,ncoast,inland,SUM)
  call flush(6)

  print *,'addshd =',int(addshd)
  totcst=totcst+int(ncoast)
  totinl=totinl+int(inland)
  if (addshd.gt.0.) go to 3

  print 100,'found',totcst,' points draining to ocean,',	&
    totinl,' not draining'

  totcst=0
  do nr=1,int(nrivr)
   sum=0.
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO REDUCTION(+:sum)
   do i=1,nip
    if (wet(i) == 0) then
     if (shedno(nr).eq.abs(rivmou(i))) sum=sum+1.
    end if
   end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END
!SMS$REDUCE (sum,SUM)
   shedsz(nr)=sum
   totcst=totcst+sum
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO
   do i=1,nip
    if (i.eq.shedno(nr))					&
     print 101,'watershed draining at',perm(i),			&
      ' (la/lo=',deg_lat(i),deg_lon(i),')    size:',shedsz(nr)
   end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END
  end do
  print *,'total points in watersheds:',totcst
  call flush(6)

! --- convert -rivmou- from integer to real for plotting by icosphere utility

!SMS$SERIAL (<inv_perm,rivmou,IN>) BEGIN
  do i=1,nip
   topo1(i)=0.
   if (rivmou(inv_perm(i)).gt.0) topo1(i)=rivmou(inv_perm(i))
  end do
  open (19,file='watersheds.bin',form='unformatted')
  write (19) nip,topo1
  close (19)
! print '(a/(5(i9,i7)))','-watersheds.bin- array',(i,nint(topo1(i)),i=1,nip)
!SMS$SERIAL END
  call flush(6)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- optional: follow drainage from -itest- to coast
 
  if (itestp.gt.0) then
   nxt=itestp
   print *,'follow drainage from point',nxt,' to coast'
1  nxt=-nxt
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO REDUCTION(max:nxt)
   do i=1,nip
    if (wet(i)==0 .and. perm(i).eq.-nxt) then
     nxt=dnhill(i)
     print 101,'water at',perm(i),' (la/lo=',deg_lat(i),deg_lon(i),	&
       ') flows to',perm(nxt),' (la/lo=',deg_lat(nxt),deg_lon(nxt),')'
101  format (2(a,i8,a,2f7.2),a)
     if (wet(nxt) > 0 .or. i.eq.dnhill(i)) then
      print 100,'point',perm(nxt),' is drainage end point'
      nxt=0
     else
      nxt=perm(nxt)
     end if
    end if
   end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END
  rnxt=nxt
!SMS$REDUCE (rnxt,MAX)
   if (rnxt.gt.0.) go to 1
   call flush(6)
  end if			! itest > 0

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- optional: more diagnostics
 
  topo1(:)=dnhillp(:)
  call stencl(topo1,1,1.,'(hycom_runoff) dnhill array')
! call pipe2(0,topo1,1,'dnhillp')
!sms$compare_var(topo1, "hycom_runoff.F90 - dnhill")
!sms$compare_var(nuphill, "hycom_runoff.F90 - nuphill")

  do edg=1,5
   write (char1,'(i1)') edg
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO
   do i=1,nip
    uphillp(edg,i)=0
    topo1(i)=0
    if (uphill(edg,i).ne.0) then
     uphillp(edg,i)=perm(uphill(edg,i))
     topo1(i)=uphillp(edg,i)
    end if
   end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END
   call stencl(topo1,1,1.,'(hycom_runoff) uphill array edg='//char1)
   topo1(:)=uphillp(edg,:)
!  call pipe2(0,topo1,1,'uphillp'//char1)
!sms$compare_var(topo1, "hycom_runoff.F90 - uphill")
  end do

  print *,'... exiting runoff'
  return
  end subroutine runoff
end module hycom_runoff
