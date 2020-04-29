subroutine lay2lay(kkold,pold,thold,				&
                   varin1,varin2,varin3,varin4,varin5,		&
                   kknew,pnew,thnew,				&
                   varou1,varou2,varou3,varou4,varou5,		&
                   targt,nip)
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!   Layer-to-layer conversion routine, obtained by combining routines
!   lay2lev4 (layer to level) and lin2stp (linear to step fct)
!
!   (Lay2lev based on Bleck, R., 1984: Vertical Coordinate Transformation
!   of Vertically Discretized Atmospheric Fields. Mon. Wea. Rev., 112,
!   2537-2541)
!
!   Specifically, lay2lay ....

!     (1) fits continuous, piecewise linear curves to input stairstep
!         profiles thold(pold), varin1(pold), varin2(pold), ...;
!     (2) creates a new stairstep profile of theta by piecewise
!         integrating over the continous theta curve such that the new
!         steps theta(pnew) match prescribed 'target' values;
!     (3) modifies the new pressure -pnew- where necessary to satisfy
!         minimum layer thickness constraints;
!     (4) integrates the continuous curves for thold, varin1, varin2,...
!         over pressure intervals obtained in (3) to produce new
!         stairstep profiles varout1(pnew), varout2(pnew), ...;
!
!   The routine is presently configured to conserve the height of the
!   air column. It can easily be modified to conserve column thermal
!   energy instead.
!
!   Rainer Bleck	2008
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! --- More details:
!
! --- 'Old' data consisting of 3-dim. arrays of
!
! ---     (a) 'kkold+1' interface pressure (or Exner fct) values (pold)
! ---     (b) 'kkold' layer averages of pot.temperature (thold)
! ---     (b) 'kkold' layer averages of dep.variables (varin1,varin2,..)
!
! --- ...are transformed from piecewise constant (stairstep-type)
! --- functions of -p- to continuous, piecewise linear functions of -p-
! --- preserving the column mean of the input data. This intrinsically
! --- nonunique curve fitting problem is rendered unique by minimizing
! --- zigzags in the curves, i.e., imposing penalties for large 2nd
! --- derivatives. Additional penalties are imposed for large overshoots
! --- in the center portion of each layer. (Sharply raising the latter
! --- penalties will force the algorithm to generate curves resembling
! --- the input stairstep curves.)
!
! --- The resulting line segments are integrated to form layers
! --- constrained to yield prescribed values of pot.temperature (targt).
! --- New values of intfc pressure, theta, and dep.variables are returned
! --- in pnew(:,:,kknew+1), thnew(:,:,kknew), varou1(:,:,kknew), ....
!
  use module_hybgen, only: regrid_1d,remap_1d
  use module_constants,only: p1000,rd,cp,perm
  use fimnamelist, only: yyyymmddhhmm,PrintIpnDiag,ptop
  use units, only: getunit, returnunit
!
  implicit none

  integer,intent(IN) :: nip,kkold,kknew

!sms$distribute(dh,2) begin
  real,intent(IN)  :: pold(kkold+1,nip),thold (kkold,nip),		&
                      varin1(kkold,nip),varin2(kkold,nip),		&
                      varin3(kkold,nip),varin4(kkold,nip),		&
                      varin5(kkold,nip),				&
                      targt (kknew,nip)
  real,intent(OUT) :: pnew(kknew+1,nip),thnew (kknew,nip),		&
                      varou1(kknew,nip),varou2(kknew,nip),		&
                      varou3(kknew,nip),varou4(kknew,nip),		&
                      varou5(kknew,nip)
!sms$distribute end
!
  integer,parameter :: nsize = 321		!  must be 5*kkold+1 or larger
  integer ipn,k,ka,kb,l,n,indx(nsize)
  real matrix(nsize,nsize),pint(nsize),solth(nsize),			&
           solu1(nsize),solu2(nsize),solu3(nsize),solu4(nsize),		&
           solu5(nsize),excol(kknew+1),thcol(kknew),			&
           v1col(kknew),v2col(kknew),v3col(kknew),v4col(kknew),		&
           v5col(kknew),tarcol(kknew),prcol(kknew+1),			&
           exold(kknew+1),exwrk(kknew+1),unusd1(nsize),unusd2(nsize),   &
           unusd3(nsize),oddev,p2ex,ex2p,arg
  real pr_extd(nsize+3),th_extd(nsize+3),tg_extd(kknew+1)
!
  real,parameter :: penlty = 0.		! penalty for midlyr overshoots
  real,parameter :: flag = -.03125	! missing data
  logical realyr,vrbos,fail
  integer :: unitno                     ! unit number to write to
  integer :: ios                        ! iostat return
  integer :: ret

#include <gptl.inc>
!
  ex2p(arg)=p1000*(arg/cp)**(cp/rd)		!  convert Pi => p
!     p2ex(arg)=cp*(arg/p1000)**(rd/cp)		!  convert p  => Pi

  print *,'entering lay2lay...'
  if (nsize.lt.5*kkold+1) stop '(nsize too small in subr.lay2lev)'
!
  n=kkold				!  number of input layers
!
  do k=1,nsize
    solth(k)=0.
    solu1(k)=0.
    solu2(k)=0.
    solu3(k)=0.
    solu4(k)=0.
    solu5(k)=0.
    do l=1,nsize
      matrix(l,k)=0.
    end do
  end do
!
  do k=1,n
! --- upper left quadrant:
    matrix(k,4*k-3)=1.
    matrix(k,4*k-2)=2.
    matrix(k,4*k-1)=2.
    matrix(k,4*k  )=2.
    matrix(k,4*k+1)=1.
! --- lower right quadrant:
    matrix(n+4*k-3,4*n+1+k)=1.
    matrix(n+4*k-2,4*n+1+k)=2.
    matrix(n+4*k-1,4*n+1+k)=2.
    matrix(n+4*k  ,4*n+1+k)=2.
    matrix(n+4*k+1,4*n+1+k)=1.
  end do
!
! --- lower left quadrant:
  do k=2,4*n
    matrix(n+k+1,k-1)=1.
    matrix(n+k-1,k+1)=1.
    matrix(n+k  ,k  )=6.
  end do
!
  do k=2,4*n-1
    matrix(n+k+1,k  )=-4.
    matrix(n+k  ,k+1)=-4.
  end do
!
  matrix(  n+1,    1)=1.
  matrix(5*n+1,4*n+1)=1.
  matrix(5*n  ,4*n+1)=-2.
  matrix(5*n+1,4*n  )=-2.
  matrix(  n+1,    2)=-2.
  matrix(  n+2,    1)=-2.
  matrix(5*n  ,4*n  )=5.
  matrix(  n+2,    2)=5.
!
! --- penalize overshoots at layer midpoints
!
  do k=1,n
    matrix(n+4*k-2,4*k-2)=matrix(n+4*k-2,4*k-2)+penlty*.25
    matrix(n+4*k-1,4*k-1)=matrix(n+4*k-1,4*k-1)+penlty
    matrix(n+4*k  ,4*k  )=matrix(n+4*k  ,4*k  )+penlty*.25
  end do
!
! --- decompose matrix
  ret = gptlstart ('ludcmp')
  call ludcmp(matrix,5*n+1,nsize,indx,oddev)
  ret = gptlstop ('ludcmp')

  unitno = getunit ()
  if (unitno < 0) then
    print*,'lay2lay: getunit failed. Stopping'
    stop
  end if
  
!SMS$PARALLEL (dh,ipn) BEGIN

!
!! Cray modification (thanks to Pete Johnsen):
!!      Make multi-line OMP directive conform to the
!!      OpenMP Application Program Interface Version 2.5 May 2005,
!!      Sect. 2.1.2 Free Source Form Directives, according to which:
!!         Continued directive lines must have an ampersand as the
!!         last nonblank character on the line, prior to any comment
!!         placed inside the directive. Continuation directive lines
!!         can have an ampersand after the directive sentinel with
!!         optional white space before and after the ampersand.
!TBH:  NOTE that this file has not been tested with OpenMP in quite a while 
!TBH:  so the directive below may require adjustment.  

!$OMP PARALLEL DO PRIVATE(vrbos,k,ret,solth,solu1,solu2,solu3,solu4,solu5, &
!$OMP             pint,realyr,tarcol,pr_extd,th_extd,tg_extd,		&
!$OMP             excol,exold,exwrk,thcol,prcol,v1col,v2col,v3col,	&
!$OMP             v4col,v5col,ios, unusd1, unusd2, unusd3)
  do ipn=1,nip
    vrbos=ipn.eq.PrintIpnDiag
!
! --- open file for saving step-by-step details of the interpolation procedure
    if (vrbos) then
      print  '(3a,i8,a)','store data for date = ',		&
           yyyymmddhhmm,', ipn =',perm(ipn),			&
           ' in file "stairstep_details" for offln plotting'
      open (unitno, file='stairstep_details', form='formatted',	&
            status='new', iostat=ios)
      if (ios /= 0) then
        print*,'lay2lay: failed to open unit number ', unitno,' for writing'
      end if
      write (unitno,'(a10,i8,4i5)') yyyymmddhhmm,ipn,kkold,4*n+1,kknew,kknew
    end if
!
    do k=1,kknew
      varou1(k,ipn)=flag
      varou2(k,ipn)=flag
      varou3(k,ipn)=flag
      varou4(k,ipn)=flag
      varou5(k,ipn)=flag
    end do
!
    do k=1,n
      solth(k)=8.*thold (k,ipn)
      solu1(k)=8.*varin1(k,ipn)
      solu2(k)=8.*varin2(k,ipn)
      solu3(k)=8.*varin3(k,ipn)
      solu4(k)=8.*varin4(k,ipn)
      solu5(k)=8.*varin5(k,ipn)
    end do
!
    do k=n+1,nsize
      solth(k)=0.
      solu1(k)=0.
      solu2(k)=0.
      solu3(k)=0.
      solu4(k)=0.
      solu5(k)=0.
    end do
!
! --- penalize overshoots at layer midpoints
!
    do k=1,n
      solth(n+4*k-2)=thold (k,ipn)*penlty*.25
      solth(n+4*k-1)=thold (k,ipn)*penlty
      solth(n+4*k  )=thold (k,ipn)*penlty*.25
!
      solu1(n+4*k-2)=varin1(k,ipn)*penlty*.25
      solu1(n+4*k-1)=varin1(k,ipn)*penlty
      solu1(n+4*k  )=varin1(k,ipn)*penlty*.25
!
      solu2(n+4*k-2)=varin2(k,ipn)*penlty*.25
      solu2(n+4*k-1)=varin2(k,ipn)*penlty
      solu2(n+4*k  )=varin2(k,ipn)*penlty*.25
!
      solu3(n+4*k-2)=varin3(k,ipn)*penlty*.25
      solu3(n+4*k-1)=varin3(k,ipn)*penlty
      solu3(n+4*k  )=varin3(k,ipn)*penlty*.25
!
      solu4(n+4*k-2)=varin4(k,ipn)*penlty*.25
      solu4(n+4*k-1)=varin4(k,ipn)*penlty
      solu4(n+4*k  )=varin4(k,ipn)*penlty*.25
!
      solu5(n+4*k-2)=varin5(k,ipn)*penlty*.25
      solu5(n+4*k-1)=varin5(k,ipn)*penlty
      solu5(n+4*k  )=varin5(k,ipn)*penlty*.25
    end do
!
! --- replace values in massless layers by values from nearest 'real' layer
!
!cc      realyr=.false.
!cc      do 8 k=1,n
!cc      if (pold(ipn,k+1).gt.pold(ipn,k)+.01) realyr=.true.
!cc      if (k.eq.1) go to 8
!cc      if (realyr .and. pold(ipn,k).ge.pold(ipn,k+1)-.01) then
!cc        solth(k)=solth(k-1)
!cc        solu1(k)=solu1(k-1)
!cc        solu2(k)=solu2(k-1)
!cc        solu3(k)=solu3(k-1)
!cc      end if
!cc   8  continue
!
!cc      realyr=.false.
!cc      do 9 k=n,1,-1
!cc      if (pold(ipn,k+1).gt.pold(ipn,k)+.01) realyr=.true.
!cc      if (k.eq.n) go to 9
!cc      if (realyr .and. pold(ipn,k).ge.pold(ipn,k+1)-.01) then
!cc       solth(k)=solth(k+1)
!cc       solu1(k)=solu1(k+1)
!cc       solu2(k)=solu2(k+1)
!cc       solu3(k)=solu3(k+1)
!cc      end if
!cc   9  continue
!
    if (vrbos) then
      write (*,101) perm(ipn),' (lay2lay) thold  inpt',			&
        (pold(k,ipn),.125*solth(k),k=1,n),pold(n+1,ipn)
!!    write (*,101) ipn,'vrbl 1 inpt',(pold(k,ipn),.125*solu1(k),k=1,n)
!!    write (*,101) ipn,'vrbl 2 inpt',(pold(k,ipn),.125*solu2(k),k=1,n)
!!    write (*,101) ipn,'vrbl 3 inpt',(pold(k,ipn),.125*solu3(k),k=1,n)
!!    write (*,101) ipn,'vrbl 4 inpt',(pold(k,ipn),.125*solu4(k),k=1,n)
!!    write (*,101) ipn,'vrbl 5 inpt',(pold(k,ipn),.125*solu5(k),k=1,n)
      write (unitno,*) 'input profile'
      write (unitno,100) (pold(k,ipn),.125*solth(k),k=1,n)
100   format (2es15.7)
    end if
101 format (i10,3x,a,' profile:'/(4(f9.1,es10.3)))
!
    ret = gptlstart ('lubksb')
    call lubksb(matrix,5*n+1,nsize,indx,solth)
    call lubksb(matrix,5*n+1,nsize,indx,solu1)
    call lubksb(matrix,5*n+1,nsize,indx,solu2)
    call lubksb(matrix,5*n+1,nsize,indx,solu3)
    call lubksb(matrix,5*n+1,nsize,indx,solu4)
    call lubksb(matrix,5*n+1,nsize,indx,solu5)
    ret = gptlstop ('lubksb')
!
! --- assign pressure values to the end points of the 4*n line segments
!
    do k=1,n
      pint(4*k-3)=pold(k,ipn      )
    end do
    pint(4*n+1)=pold(kkold+1,ipn)
    do k=1,n
      pint(4*k-2)=.75*pint(4*k-3)+.25*pint(4*k+1)
      pint(4*k-1)=.50*pint(4*k-3)+.50*pint(4*k+1)
      pint(4*k  )=.25*pint(4*k-3)+.75*pint(4*k+1)
    end do
!
    if (vrbos) then
      write (*,101) perm(ipn),' (lay2lay) theta  pcwise lin.',		&
        (pint(k),solth(k),k=1,4*n+1)
!!    write (*,101) ipn,'vrbl 1 pcwise lin.',(pint(k),solu1(k),k=1,4*n+1)
!!    write (*,101) ipn,'vrbl 2 pcwise lin.',(pint(k),solu2(k),k=1,4*n+1)
!!    write (*,101) ipn,'vrbl 3 pcwise lin.',(pint(k),solu3(k),k=1,4*n+1)
!!    write (*,101) ipn,'vrbl 4 pcwise lin.',(pint(k),solu4(k),k=1,4*n+1)
!!    write (*,101) ipn,'vrbl 5 pcwise lin.',(pint(k),solu5(k),k=1,5*n+1)
      write (unitno,*) 'piecewise linear fit'
      write (unitno,100) (pint(k),solth(k),k=1,4*n+1)
    end if

! --- remove superadiabats from piecewise linear profile
    call superad(4*n+1,solth,pint,vrbos,perm(ipn))
!
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
! --- Step 1 completed. Now break -solth- into stairsteps whose 'risers'
! --- are at prescribed -targt- values
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!
    tarcol(:)=targt(:,ipn)
!
! --- extend fitted curve at bottom and top to match range of target values
! --- extend list of target values to match range of fitted curve
!
    pr_extd(2:4*n+2)=pint(1:4*n+1)
    pr_extd(    1)=pint(    1)
    pr_extd(4*n+3)=pint(4*n+1)
!
    th_extd(2:4*n+2)=solth(1:4*n+1)
    th_extd(    1)=min(tarcol(    1),solth(    1))
    th_extd(4*n+3)=max(tarcol(kknew),solth(4*n+1))
!
    tg_extd(1:kknew)=tarcol
    tg_extd(      1)=min(tarcol(    1),solth(    1))
    tg_extd(kknew+1)=max(tarcol(kknew),solth(4*n+1))
!
    ret = gptlstart ('lin2stp')
    call lin2stp(th_extd,pr_extd,4*n+3,tg_extd,excol,kknew,		&
                 vrbos,perm(ipn))
    ret = gptlstop ('lin2stp')
    do k=kknew,2,-1
      excol(k+1)=excol(k)
    end do
    excol(1)=pint(1)
    thcol(:)=tarcol(:)
!
    if (vrbos) then
      write (unitno,*) 'output profile before hybridization'
      write (unitno,100) (excol(k),thcol(k),k=1,kknew)
    end if
!
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
! --- Step 2 completed. Now hybridize the grid
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!
! --- FIM atmosphere ends at pressure -ptop-
    excol(kknew+1)=cp*(ptop/p1000)**(rd/cp)
    do k=kknew,kknew/2,-1
      excol(k)=max(excol(k),excol(k+1))
    end do
!
! --- inflate massless layers to create hybrid-isentropic grid.
! --- it is assumed here that -pint- represents Exner fctn, *not* pressure
!
    prcol(  1)=ex2p(excol(  1))
    do k=1,kknew
      excol(k+1)=min(excol(k),excol(k+1))		! remove superadiabats
      prcol(k+1)=ex2p(excol(k+1))
    end do
    exold(:)=excol(:)
!
! --- for consistency between initialization and time integration,
! --- regrid_1d (part of the FIM grid generator) is used to here
! --- to inflate massless layers. however, vertical regridding of
! --- dependent variables will be done separately by lin2stp.
!
    if (vrbos) then
      print *,'lay2lay calling regrid_1d/remap_1d ...'
    endif
!
    ret = gptlstart ('regrid_1d_fromlay2lay')
    call regrid_1d(0,tarcol,thcol,excol,prcol,vrbos,perm(ipn),ipn)
    ret = gptlstop ('regrid_1d_fromlay2lay')
    exwrk(:)=exold(:)
    unusd1(:) = 0.
    unusd2(:) = 0.
    unusd3(:) = 0.
    ret = gptlstart ('regrid_1d_fromlay2lay')
    call remap_1d(0,tarcol,1,thcol,unusd1,unusd2,exwrk,excol,		&
                  unusd3,prcol,-1,vrbos,perm(ipn),ipn)
    ret = gptlstop ('regrid_1d_fromlay2lay')
!
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
! --- Step 3 completed. Now interpolate input variables to the new grid.
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!
    do k=1,4*n+1
      pint(k)=max(pint(k),excol(kknew+1))
    end do
! 
! --- do full-column regridding for all variables
!
    ret = gptlstart ('lin2stp')
     if (vrbos) print *,'convert vrbl1 from linear to steppy'
    call lin2stp(pint,solu1,4*n+1,excol,v1col,kknew,vrbos,perm(ipn))
     if (vrbos) print *,'convert vrbl2 from linear to steppy'
    call lin2stp(pint,solu2,4*n+1,excol,v2col,kknew,vrbos,perm(ipn))
     if (vrbos) print *,'convert vrbl3 from linear to steppy'
    call lin2stp(pint,solu3,4*n+1,excol,v3col,kknew,vrbos,perm(ipn))
     if (vrbos) print *,'convert vrbl4 from linear to steppy'
    call lin2stp(pint,solu4,4*n+1,excol,v4col,kknew,vrbos,perm(ipn))
     if (vrbos) print *,'convert vrbl5 from linear to steppy'
    call lin2stp(pint,solu5,4*n+1,excol,v5col,kknew,vrbos,perm(ipn))
     if (vrbos) print *,'convert theta from linear to steppy'
    call lin2stp(pint,solth,4*n+1,excol,thcol,kknew,vrbos,perm(ipn))
    ret = gptlstop ('lin2stp')
!
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
! --- Steps 4 completed. Finish up.
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!
    do k=1,kknew
      pnew  (k,ipn)=excol(k)
      thnew (k,ipn)=thcol(k)
      varou1(k,ipn)=v1col(k)
      varou2(k,ipn)=v2col(k)
      varou3(k,ipn)=v3col(k)
      varou4(k,ipn)=v4col(k)
      varou5(k,ipn)=v5col(k)
    end do
    pnew(kknew+1,ipn)=excol(kknew+1)
!
    if (vrbos) then
      write (*,101) perm(ipn),' (lay2lay) theta  outp',		&
        (excol(k),thcol(k),k=1,kknew),excol(kknew+1)
      write (*,101) perm(ipn),' (lay2lay) vrbl 1 outp',		&
        (excol(k),v1col(k),k=1,kknew),excol(kknew+1)
      write (*,101) perm(ipn),' (lay2lay) vrbl 2 outp',		&
        (excol(k),v2col(k),k=1,kknew),excol(kknew+1)
      write (*,101) perm(ipn),' (lay2lay) vrbl 3 outp',		&
        (excol(k),v3col(k),k=1,kknew),excol(kknew+1)
      write (*,101) perm(ipn),' (lay2lay) vrbl 4 outp',		&
        (excol(k),v4col(k),k=1,kknew),excol(kknew+1)
      write (*,101) perm(ipn),' (lay2lay) vrbl 5 outp',		&
        (excol(k),v5col(k),k=1,kknew),excol(kknew+1)
      write (unitno,*) 'output profile after hybridization'
      write (unitno,100) (excol(k),thcol(k),k=1,kknew)
      close (unitno)
    end if
  end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END
!
  call returnunit(unitno)

  print *,'... exiting lay2lay'
  return
end subroutine lay2lay


  subroutine superad(kdm,solth,pint,vrbose,ipn)

! --- remove superadiabats from piecewise linear temperature profile

  integer,intent(IN)    :: kdm			! number of levels
  real   ,intent(INOUT) :: solth(kdm)		! pot.temperature
  real   ,intent(IN)    :: pint(kdm)		! Exner function
  logical,intent(IN)    :: vrbose
  integer,intent(IN)    :: ipn
  real*8  theta(kdm),exnr(kdm),dphi,test,base,thmod
  integer k,ka,kb
  logical vrbos,fail

  vrbos=vrbose
  theta(:)=solth(:)
  exnr(:)=pint(:)

  do k=1,kdm-1
   if (theta(k+1).lt.theta(k)) then			! found superadiabat
    if (vrbos .or. theta(k+1).lt.theta(k)-5.) vrbos=.true.

! --- find set of layers whose combined geopot.thickness can be replicated
! --- by -theta- profile modified to remove unstable stratification.

    base=theta(k)
    dphi=(theta(k)+theta(k+1))*(exnr(k)-exnr(k+1))
    test=2.*base              *(pint(k)-pint(k+1))
    do ka=k+1,kdm-1			! search column above level -k-
     dphi=dphi+(theta(ka)+theta(ka+1))*(pint(ka)-pint(ka+1))
     if (test +(base     +theta(ka+1))*(pint(ka)-pint(ka+1)).gt.dphi	&
         .or. theta(ka+1).lt.base					&
         .or. pint(ka+1).gt.pint(k)-.001) then
      test=test+2.*base*(pint(ka)-pint(ka+1))
      if (ka.eq.kdm-1) print 102,ipn,k,					&
        ' (lay2lay) running out of layers to remove superadiabat:',	&
        theta(k:kdm)
     else				! no need to go beyond level ka+1

! --- assume -theta- to vary linearly between levels -k- and -ka-.
! --- use geopot.thickness constraint to determine theta(ka)

      if (vrbos) print 102,ipn,k,					&
        ' (lay2lay) found superadiabat:',(theta(kb),kb=k,ka+1)
 102  format(i8,i4,a/(10f8.2))

      thmod=(dphi-base*(exnr(k)-exnr(ka))				&
       -theta(ka+1)*(exnr(ka)-exnr(ka+1)))/(exnr(k)-exnr(ka+1))

      fail=.false.
      if (thmod.lt.theta(k)-.01 .or. thmod.gt.theta(ka+1)+.01) then
       if (.not.vrbos) print 102,ipn,k,					&
         ' (lay2lay) found superadiabat:',(theta(kb),kb=k,ka+1)
       print 102,ipn,ka,' (lay2lay) thmod out of range:',thmod
       if (thmod.lt.theta(k)-.5 .or. thmod.gt.theta(ka+1)+.5)		&
         fail=.true.
       vrbos=.true.
      end if

      if (.not.fail) then
       do kb=k+1,ka-1
        wgt=1.
        if (exnr(k).gt.exnr(ka))					&
          wgt=(exnr(kb)-exnr(ka))/(exnr(k)-exnr(ka))
        theta(kb)=theta(k)*wgt+thmod*(1.-wgt)
       end do
       theta(ka)=thmod
       if (vrbos) print 102,ipn,k,					&
         ' (lay2lay) modified profile:',(theta(kb),kb=k,ka+1)

! --- compute new dphi for testing purposes
       if (vrbos) then
        test=0.
        do kb=k,ka
         test=test+(theta(kb)+theta(kb+1))*(exnr(kb)-exnr(kb+1))
        end do
        print '(i8,i4,a,2f11.2)',ipn,k,					&
          ' (lay2lay) old, new dphi:',dphi,test
       end if
      end if
      exit
     end if			! test < dphi
    end do			! ka loop
   end if			! superadiabat
  end do			! k loop

! --- if not all superadiabats got removed, do it now by brute force
  do k=2,kdm-1
   theta(k)=max(theta(k),theta(k-1))
  end do
  do k=kdm-1,1,-1
   theta(k)=min(theta(k),theta(k+1))
  end do
  solth(:)=theta(:)
  return
  end subroutine superad
!
!
SUBROUTINE lubksb(a,n,np,indx,b)
  INTEGER n,np,indx(n)
  REAL a(np,np),b(n)
  INTEGER i,ii,j,ll
  REAL sum
  ii=0
!
  do i=1,n
    ll=indx(i)
    sum=b(ll)
    b(ll)=b(i)
    if (ii.ne.0)then
      do j=ii,i-1
        sum=sum-a(i,j)*b(j)
      end do
    else if (sum.ne.0.) then
      ii=i
    endif
    b(i)=sum
  end do
  do i=n,1,-1
    sum=b(i)
    do j=i+1,n
      sum=sum-a(i,j)*b(j)
    end do
    b(i)=sum/a(i,i)
  end do
  return
END SUBROUTINE lubksb
!  (C) Copr. 1986-92 Numerical Recipes Software 'W3.
!
!
SUBROUTINE ludcmp(a,n,np,indx,d)
  INTEGER n,np,indx(n)
  REAL d,a(np,np),TINY
  PARAMETER (TINY=1.0e-20)
  INTEGER i,imax,j,k
  REAL aamax,dum,sum,vv(n)
  d=1.
  do i=1,n
    aamax=0.
    do j=1,n
      if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
    end do
    if (aamax.eq.0.) pause 'singular matrix in ludcmp'
    vv(i)=1./aamax
  end do
  do j=1,n
    do i=1,j-1
      sum=a(i,j)
      do k=1,i-1
        sum=sum-a(i,k)*a(k,j)
      end do
      a(i,j)=sum
    end do
    aamax=0.
    do i=j,n
      sum=a(i,j)
      do k=1,j-1
        sum=sum-a(i,k)*a(k,j)
      end do
      a(i,j)=sum
      dum=vv(i)*abs(sum)
      if (dum.ge.aamax) then
        imax=i
        aamax=dum
      endif
    end do
    if (j.ne.imax)then
      do k=1,n
        dum=a(imax,k)
        a(imax,k)=a(j,k)
        a(j,k)=dum
      end do
      d=-d
      vv(imax)=vv(j)
    endif
    indx(j)=imax
    if(a(j,j).eq.0.)a(j,j)=TINY
    if(j.ne.n)then
      dum=1./a(j,j)
      do i=j+1,n
        a(i,j)=a(i,j)*dum
      end do
    endif
  end do
  return
END SUBROUTINE ludcmp
!  (C) Copr. 1986-92 Numerical Recipes Software 'W3.
