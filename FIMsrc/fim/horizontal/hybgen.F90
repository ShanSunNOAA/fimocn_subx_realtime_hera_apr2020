module module_hybgen
  use findmaxmin2
  use global_bounds, only: ims, ime, ips, ipe, ihe

#include <gptl.inc>

contains

!***********************************************************************
!     hybgen
!       Adapted from 2nd generation restep-based grid generator in HYCOM
!       R. Bleck     Mar 2006
!       R. Bleck     Mar 2008       arbitrary number of tracers
!       R. Bleck     Jul 2009       inflation of top layers
!       R. Bleck     Nov 2009       merged with hybgen_sig
!       S. Sun       Nov 2009       added option to smooth th-dot
!       S. Sun       Nov 2009       added lateral intfc smoothing option
!***********************************************************************

  subroutine hybgen (its,   &
            targt,          &	! target potential temperature (1-D)
            us3d,vs3d,      &	! horiz.velocity (3-D)
            tr3d,           &	! mass field tracers (ntra 3-D fields)
            sdot,           &	! intfc displacement sdot*(dp/ds)*dt (3-D)
            ex3d,dp3d,pr3d, &	! Exner function, layer thknss, pressure (3-D)
            TimingBarriers)
 
    use module_control  ,only: nvlp1,ntra,kbl,nip,dt
    use fimnamelist     ,only: nvl,PrintIpnDiag,intfc_smooth,pure_sig,miglim
    use module_constants,only: p1000,rd,cp,sigak,sigbk,perm,smoo_coeff
    use module_variables,only: worka,workb
    use dffusn,          only: dffusn_lyr, dffusn_lev

    implicit none

! Type and dimension external variables:

    integer,intent(IN)     :: its			! model time step
    real   ,intent(IN)     :: targt(nvl)	        ! target pot.temp.
!JR Let SMS distribute arrays because they may be exchanged or passed to compare_var
!sms$distribute (dh,2) begin
    real   ,intent(INOUT)  :: tr3d(nvl,nip,ntra)	! mass field tracers
    real   ,intent(INOUT)  :: ex3d(nvlp1,nip)	! Exner function
    real   ,intent(INOUT)  :: dp3d(nvl  ,nip)	! layer thickness
    real   ,intent(INOUT)  :: pr3d(nvlp1,nip)	! pressure
    real   ,intent(OUT)    :: sdot(nvlp1,nip)	! sdot*(dp/ds)
    real   ,intent(INOUT)  :: us3d(nvl  ,nip)	! u velocity
    real   ,intent(INOUT)  :: vs3d(nvl  ,nip)	! v velocity

! Local variables

    real    :: exsmo3d(nvlp1,nip)		! 3d array for smoothing
    real    :: exdif3d(nvl,nip)
!SMS$DISTRIBUTE END
    logical,intent(in)     :: TimingBarriers
    integer :: ipn,k,k1,k2
    logical :: vrbos				! switch for 'verbose' mode
    real    :: trcol(nvl,ntra)			! tracer column vectors
    real    :: valmin
    logical :: thsmoo = .false.			! use smoothed thdot in thcol
    character :: string*20
    integer :: ret                              ! return code from gptl routines

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- if thsmoo=true, smooth diabatic theta increment (thdot) laterally.
! --- the smoothed field is used   o n l y   to determine coordinate movement
! --- (regridding). the regular theta field (tracer 1) is not smoothed.

    ret = gptlstart ('hybgen')

    if (pure_sig) then
      thsmoo       = .false.
      intfc_smooth = 0.
      miglim = -1
    end if

!JR TODO: Need to thread the true branch

    if (thsmoo .and. its > 0) then
      workb(:,:) = tr3d(:,:,1) - worka(:,:)	! worka = theta before physics call
!JR exchange (delp,fld) moved from dffusn_lyr to here
      if (TimingBarriers) then
        ret = gptlstart ('hybgen_barrier')
!SMS$BARRIER
        ret = gptlstop ('hybgen_barrier')
      end if
      ret = gptlstart ('hybgen_exchange')
!SMS$EXCHANGE(dp3d,workb)
      ret = gptlstop  ('hybgen_exchange')

      ret = gptlstart ('dffusn_lyr')
      call dffusn_lyr (workb, dp3d, dt*20.)		! lateral thdot smoothing
      ret = gptlstop  ('dffusn_lyr')

      worka(:,:) = worka(:,:) + workb(:,:)

    else						! thsmoo = false

!sms$ignore begin
!$OMP PARALLEL DO PRIVATE (k) SCHEDULE (runtime)
      do ipn=ips,ipe
        do k=1,nvl
          worka(k,ipn) = tr3d(k,ipn,1)
        end do
      end do
!$OMP END PARALLEL DO
!sms$ignore end

    end if					! thsmoo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
!sms$compare_var(ex3d, "hybgen.F90 - ex3d7 ")

!sms$ignore begin
!$OMP PARALLEL DO PRIVATE (vrbos) SCHEDULE (runtime)
    do ipn=ips,ipe
      vrbos = ipn == PrintIpnDiag

!   if (vrbos .and. pure_sig) then
!    write (6,'(a/(5f14.6))') '(hybgen) sigak array:',sigak
!    write (6,'(a/(5f14.6))') '(hybgen) sigbk array:',sigbk
!   end if

! --- call single-column version of hybgen
  
! --- interface smoothing requires separation of regridding and remappping.
! --- step 1: subr. regrid_1d does the regridding 
! --- step 2: subr. dffusn_lev does the (optional) interface smoothing
! --- step 3: subr. remap_1d does the remapping, i.e., it vertically advects
! ---         all prognostic variables (after re-inflating layers that may
! ---         have become too thin during smoothing)

      exsmo3d(:,ipn) = ex3d(:,ipn)
      call regrid_1d (its, targt, worka(1,ipn), exsmo3d(1,ipn),		&
                      pr3d(1,ipn), vrbos, perm(ipn), ipn)
! --- (on return from regrid_1d, exsmo3d is new, pr3d is old)
    end do
!$OMP END PARALLEL DO
!sms$ignore end

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (intfc_smooth > 0.) then
      if (TimingBarriers) then
        ret = gptlstart ('hybgen_barrier')
!SMS$BARRIER
        ret = gptlstop ('hybgen_barrier')
      end if
!JR exchange (fld) moved from dffusn_lyr to here
      ret = gptlstart ('hybgen_exchange')
!SMS$EXCHANGE(exsmo3d)
      ret = gptlstop  ('hybgen_exchange')

      ret = gptlstart ('dffusn_lev')
      call dffusn_lev (exsmo3d, smoo_coeff, nvlp1, kbl+2, nvl)
      ret = gptlstop  ('dffusn_lev')

! --- check occasionally whether smoothing is generating neg.lyr.thknss
      if (mod(its,100) == 0) then
        do k=nvl/3,2*nvl/3
!SMS$PARALLEL(dh, ipn) BEGIN
          do ipn=1,nip
            exdif3d(k,ipn) = exsmo3d(k,ipn) - exsmo3d(k+1,ipn)
          end do
          valmin = minval(exdif3d(k,:))
!SMS$REDUCE(valmin,MIN)
!SMS$PARALLEL END
          if (valmin < 0.) then
            write (string,'(a,i2)') 'exdif aftr smoo k=',k
            call findmxmn2(exdif3d,nvl,nip,k,string)
          end if
        end do
      end if
    end if
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!sms$ignore begin
!$OMP PARALLEL DO PRIVATE (vrbos,trcol) SCHEDULE (runtime)
    do ipn=ips,ipe
      vrbos = ipn == PrintIpnDiag

      trcol(:,:) = tr3d(:,ipn,:)

! --- (on entry to remap_1d, ex3d,pr3d are old)
      call remap_1d (its, targt, ntra, trcol, us3d(1,ipn), vs3d(1,ipn),	&
                     ex3d(1,ipn), exsmo3d(1,ipn), dp3d(1,ipn),		&
                     pr3d(1,ipn), miglim, vrbos, perm(ipn), ipn)
! --- (on exit from remap_1d, exsmo3d,pr3d,dp3d are new)

      tr3d(:,ipn,:) = trcol(:,:)
      sdot(:,ipn) = exsmo3d(:,ipn) - ex3d(:,ipn)
      ex3d(:,ipn) = exsmo3d(:,ipn)
    end do
!$OMP END PARALLEL DO
!sms$ignore end

!!sms$compare_var(ex3d, "hybgen.F90 - ex3d8 ")

    ret = gptlstop ('hybgen')

    return
  end subroutine hybgen

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Nothing below needs to be processed by SMS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine regrid_1d (its, targt, thcol, excol, prcol, vrbos, glbl, ipn)

  use module_control  ,only: nvlp1,kbl,nip,dt
  use fimnamelist     ,only: nvl,thktop,pure_sig,ptop
  use module_constants,only: p1000,rd,cp,dpsig,grvity,deg_lat,deg_lon,	&
                             sigak,sigbk
  implicit none

  real,intent(IN)    :: prcol(nvlp1)	! must be consistent with initl excol
  real,intent(INOUT) :: excol(nvlp1)	! initl/final exner function
  real,intent(IN)    :: targt(nvl)	! target pot.temp.
  real,intent(INOUT) :: thcol(nvl)	! actual pot. temp.
  integer,intent(IN) :: its,glbl,ipn
  logical,intent(IN) :: vrbos

! Local variables

  integer :: k,k1,k2,iter		! layer/level indices
  real    :: thnew(nvl)			! new pot.temp.
  real    :: exnew(nvlp1)		! new Exner function
  real    :: exwrk(nvlp1)		! intermediate column values
  real    :: heatfx(nvl)		! turbulent vertical heat flux
  real    :: arg,wgt1,wgt2,try,coeff,cnv
  real    :: ex2p,dex2dp,p2ex,dp2dex
  real    :: eqlb_slak
  logical :: event
  real,parameter :: dffudt=.1		! therm.diffu.coeff x time step [m^2]
  real,parameter :: thin=rd/p1000	! approx. 1 Pa in Exner fcn units
  real,parameter :: tolrnce=0.001	! in degrees

  ex2p(arg) = p1000*(arg/cp)**(cp/rd)	!  convert Pi => p
  p2ex(arg) = cp*(arg/p1000)**(rd/cp)	!  convert p  => Pi

!sms$ignore begin

  if (vrbos) print *,'entering restep_1d, ipn=',glbl

  eqlb_slak = .1*float(its)/(10.+float(its))

#ifdef DEBUGPRINT
  if (vrbos) then
    write (6,99) its,glbl,deg_lat(ipn),deg_lon(ipn),'i n p u t  profile:'
    do k2=1,nvl,10
     write (6,100) (prcol(k1),k1=k2,min(nvlp1,k2+10) )
     write (6,102) (excol(k1),k1=k2,min(nvlp1,k2+10) )
     write (6,101) (thcol(k1),k1=k2,min(nvl  ,k2+9 ) )
     write (6,101) (targt(k1),k1=k2,min(nvl  ,k2+9 ) )
     write (6,101)
    end do
  end if
#endif
 99  format ('its,ipn=',i6,i8,'  lat/lon=',2f7.1,'  hybgen  ',a/	&
       '(4-line groups: pressure, exn.fcn, theta, target)')
!100 format (-2p,11f7.1)
!102 format (   11f7.1)
 100 format (-2p,11f7.2)
 102 format (   11f7.2)
 101 format (5x,10f7.2)

  do k=1,nvl
    if (excol(k) < excol(k+1)-thin) then
      write (6,99) its,glbl,deg_lat(ipn),deg_lon(ipn),'i n p u t  profile:'
      do k2=1,nvl,10
        write (6,100) (prcol(k1),k1=k2,min(nvlp1,k2+10) )
        write (6,102) (excol(k1),k1=k2,min(nvlp1,k2+10) )
        write (6,101) (thcol(k1),k1=k2,min(nvl  ,k2+9 ) )
        write (6,101) (targt(k1),k1=k2,min(nvl  ,k2+9 ) )
        write (6,101)
      end do
      write (6,103) 'nonmonotonic Exner function on input at ipn,k =',	&
                    glbl,k,excol(k),excol(k+1)
103   format (a,i8,i4,2f9.2)
!   stop '(error: non-monotonic Exner fcn)'
      excol(k+1) = excol(k)
    end if
  end do

  if (pure_sig) then
!  if (vrbos) then
!   write (6,'(a/(5f14.6))') '(regrid_1d) sigak array:',sigak
!   write (6,'(a/(5f14.6))') '(regrid_1d) sigbk array:',sigbk
!  end if

! --- ><><><><><><><><><><><><><><><><><><>
! ---   restore to sigma coordinate grid
! --- ><><><><><><><><><><><><><><><><><><>

    exnew(    1) = excol(    1)
    exnew(nvlp1) = excol(nvlp1)
    do k=2,nvl
! --- sigak,sigbk define the sigma-p levels used in the GFS model.
      exnew(k) = p2ex(sigak(k) + sigbk(k)*prcol(1))
    end do

#ifdef DEBUGPRINT
    if (vrbos) then
      write (6,98) its,glbl,'o u t p u t     p r e s s u r e'
98    format ('time step',i6,'  ipn=',i8,'   hybgen    ',a)
      do k2=1,nvl,10
        write (6,100) (ex2p(exnew(k1)),k1=k2,min(nvlp1,k2+10) )
      end do
    end if
#endif
  else

! --- ><><><><><><><><><><><><><><><><><><>
! ---    restore to hybrid-isentropic grid
! --- ><><><><><><><><><><><><><><><><><><>

   ! --- eliminate static instabilities

    do k=nvl-1,1,-1
      thcol(k) = min(thcol(k),thcol(k+1))		! non-conservative!
    end do
    exwrk(:) = excol(:)

! --- 'thcol' profile is now stably (or neutrally) stratified.

! --- check whether theta values exceed 'targt' range at bottom of column.
! --- if so, homogenize layers to eliminate this condition

!?     do  k=1,nvl-1
!?      if (thcol(k) >= targt(1)-tolrnce) then
!?       exit				!  no action required
!?      else
!? ! --- average of layers 1...k colder than lowest target
!? 
!?       if (vrbos) then
!?        write (*,107) glbl,'lyrs',1,k+1,' bfore range limiting:',	&
!?         (exwrk(k2),thcol(k2),k2=1,k+1),exwrk(k+2)
!?       end if
!? 
!? ! --- compute pot.temp. obtained by homogenizing layers 1,...,k+1
!?       wgt1=max(thin,exwrk(k  )-exwrk(k+1))
!?       wgt2=max(  0.,exwrk(k+1)-exwrk(k+2))
!?       try=(thcol(k)*wgt1+thcol(k+1)*wgt2)/(wgt1+wgt2)
!?       if (try < targt(1)+tolrnce) then
!? ! --- average of layers 1,...,k+1 still too cold. continue adding layers
!?        exwrk(k+1)=exwrk(1)
!?        do k1=1,k+1
!?         thcol(k1)=try
!?        end do
!? 
!?        if (vrbos) then
!?         write (*,107) glbl,'lyrs',1,k+1,' after range limiting:',	&
!?          (exwrk(k2),thcol(k2),k2=1,k+1),exwrk(k+2)
!?        end if
!? 
!?       else
!? ! --- adding all of layer k+1 is overkill; entrain only part of lyr k+1
!?        exwrk(k+1)=min(exwrk(k  ),max(exwrk(k+2),			&
!?                      (exwrk(k  )*(thcol(k  )-targt(1))		&
!?                  +    exwrk(k+1)*(thcol(k+1)-thcol(k)))		&
!?                  /               (thcol(k+1)-targt(1))))
!?        do k1=1,k
!?         thcol(k1)=targt(1)
!?        end do
!? 
!?        if (vrbos) then
!?         write (*,107) glbl,'lyrs',1,k+1,' after range limiting:',	&
!?          (exwrk(k2),thcol(k2),k2=1,k+1),exwrk(k+2)
!?        end if
!? 
!?        exit				!  range limiting completed
!?       end if
!?      end if
!?     end do

! --- check whether theta values exceed 'targt' range at top of column.
! --- if so, homogenize layers to eliminate this condition

    do k=nvl,2,-1
      if (thcol(k) <= targt(nvl)+tolrnce) then
        exit				!  no action required
      end if
! --- average of layers k...nvl warmer than highest target

#ifdef DEBUGPRINT
      if (vrbos) then
        write (*,107) glbl,'lyrs',k-1,nvl,' bfore range limiting:',	&
             (exwrk(k2),thcol(k2),k2=k-1,nvl),exwrk(nvlp1)
      end if
#endif

! --- compute pot.temp. obtained by homogenizing layers k-1,...,nvl
      wgt1 = max(  0.,exwrk(k-1)-exwrk(k  ))
      wgt2 = max(thin,exwrk(k  )-exwrk(k+1))
      try = (thcol(k-1)*wgt1 + thcol(k)*wgt2)/(wgt1+wgt2)
      if (try > targt(nvl)-tolrnce) then
!---average of layers k-1,...,nvl still too warm. continue adding layers
        exwrk(k) = exwrk(nvlp1)
        do k1=k-1,nvl
          thcol(k1) = try
        end do

#ifdef DEBUGPRINT
        if (vrbos) then
          write (*,107) glbl,'lyrs',k-1,nvl,' after range limiting:',	&
               (exwrk(k2),thcol(k2),k2=k-1,nvl),exwrk(nvlp1)
        end if
#endif

      else
! --- adding all of layer k-1 is overkill; entrain only part of lyr k-1
        exwrk(k) = min(exwrk(k-1),max(exwrk(k+1),			&
                   (exwrk(k  )*(thcol(k  )-thcol(k-1))			&
               +    exwrk(k+1)*(targt(nvl)-thcol(k  )))			&
               /               (targt(nvl)-thcol(k-1))))
        do k1=k,nvl
          thcol(k1) = targt(nvl)
        end do

#ifdef DEBUGPRINT
        if (vrbos) then
          write (*,107) glbl,'lyrs',k-1,nvl,' after range limiting:',	&
               (exwrk(k2),thcol(k2),k2=k-1,nvl),exwrk(nvlp1)
        end if
#endif
        exit				!  range limiting completed
      end if
    end do
107 format (i7,2x,a,i3,'-',i3,a,20(f7.1,f8.3))

! --- now convert column to purely isentropic coordinates, i.e.,
! --- find pressure levels where all pot.temps are on target

    call restp_1d(thcol,exwrk,nvl,thnew,exnew,targt,nvl,vrbos,glbl)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- invoke heat diffusion (McDougall-Dewar) to inflate thin layers

    heatfx(nvl) = 0.
    heatfx(  1) = 0.
    do iter=1,3				! apply the scheme repeatedly
      exwrk(:) = exnew(:)
      do k=2,nvl-1
        heatfx(k) = 0.    
        if (exwrk(k  ) < exwrk(    1)-.01 .and.				&
            exwrk(k+1) > exwrk(nvlp1)+.01) then
          cnv = grvity/targt(k)				! Exner fcn units/meter
          heatfx(k) = dffudt*cnv*cnv*.5*(targt(k+1)-targt(k-1))		&
                      /max(.03*cnv,exwrk(k)-exwrk(k+1))
        end if
      end do

      do k=1,nvl-1
        if (exwrk(k+1) < exwrk(    1)-.01 .and.				&
            exwrk(k+1) > exwrk(nvlp1)+.01) then

          exnew(k+1) = max(exwrk(k+2),min(exwrk(k),			&
                           exwrk(k+1)+(heatfx(k+1)-heatfx(k))		&
                       /(targt(k+1)-targt(k))))
        end if
      end do

      event = .false.
      do k=2,nvlp1
        if (exnew(k) > exnew(k-1)+thin) then
          event=.true.
          print '(a,i8,i4,a,2F8.1)','dp<0 due to heat diffusion at ipn,k=',&
               glbl,k,'  lat/lon =',deg_lat(ipn),deg_lon(ipn)
        end if
      end do				! iter

      if (vrbos .or. event) then
        print '(i8,a,i2,a)',glbl,' heat diffusion, iter',iter,		&
              '  Ex.fcn (3-line groups: old,new,dif x 10^4)'
        do k2=1,nvl,10
          write (6,108) (exwrk(k1),k1=k2,min(nvlp1,k2+9) )
          write (6,108) (exnew(k1),k1=k2,min(nvlp1,k2+9) )
          write (6,109) (int(1.e4*(exnew(k1)-exwrk(k1))),k1=k2,min(nvlp1,k2+9) )
          write (6,108)
        end do
108     format (10f8.2)
109     format (10i8)
      end if
      if (.not.event) exit
    end do			! iter
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! --- suppress computational mode causing gradual depletion of alternate layers

    call equilb(thnew,exnew,nvl,eqlb_slak,vrbos,glbl)

! --- inflate massless layers

    call inflate(prcol(1),exnew,vrbos,glbl)
  end if			!  pure sigma or hybrid-isentropic option

  if (vrbos) then
!  write (6,104) glbl,'  hybgen: old Exner fcn (excol)',excol
!  write (6,104) glbl,'  hybgen: new Exner fcn (exnew)',exnew
!  write (6,104) glbl,'  hybgen: Exner fcn tndcy (displ)',exnew-excol
   write (6,104) glbl,'  hybgen: old prs',(ex2p(excol(k)),k=1,nvlp1)
   write (6,104) glbl,'  hybgen: new prs',(ex2p(exnew(k)),k=1,nvlp1)
   write (6,104) glbl,'  hybgen: prs tndcy',			&
      (ex2p(exnew(k))-ex2p(excol(k)),k=1,nvlp1)
 104 format (i8,a/-2p,(8f9.2))
  end if

  excol(:) = exnew(:)
  if (vrbos) print *,'...exiting regrid_1d'
  return

!sms$ignore end
end subroutine regrid_1d


subroutine remap_1d (its, targt, ntra, trcol, ucol, vcol, exold, excol,	&
                      dpcol, prcol, miglim, vrbos, glbl, ipn)

  use module_control  ,only: nvlp1,kbl,nip,dt
  use fimnamelist     ,only: nvl,thktop,pure_sig,intfc_smooth,slak
  use module_constants,only: p1000,rd,cp,dpsig,deg_lat,deg_lon
  implicit none

  integer,intent(IN) :: ntra		! number of tracer fields
  real,intent(INOUT) :: trcol(nvl,ntra)
  real,intent(INOUT) :: ucol(nvl)
  real,intent(INOUT) :: vcol(nvl)
  real,intent(IN)    :: exold(nvlp1)
  real,intent(INOUT) :: excol(nvlp1)
  real,intent(INOUT) :: dpcol(nvl)
  real,intent(INOUT) :: prcol(nvlp1)
  real,intent(IN)    :: targt(nvl)	! target pot.temp.
  integer,intent(IN) :: its,glbl,ipn
  integer,intent(IN) :: miglim		! intfc migration limit (< 0: no limit)
  logical,intent(IN) :: vrbos

! Type and dimension of local variables:

  integer :: k,k1,k2,n			! layer/level indices
  real    :: trnew(nvl,ntra)		! new tracers (1=pot.temp.)
  real    :: prnew(nvlp1)		! new pres & lyr thickness
  real    :: exnew(nvlp1)		! new Exner fcn
  real    :: unew(nvl),vnew(nvl)	! new velocities
  real    :: exwrk(nvlp1)		! intermediate column values
  real    :: thwrk(nvlp1)		! intermediate column values
  real    :: pk1col(nvlp1)		! vert.coord. used for theta remap
  real    :: pk1new(nvlp1)		! vert.coord. used for theta remap

  real    :: ex2p,arg,colin,clout,colscl,qq
  real    :: dplo,dpup,devilo,deviup,tha,thb
  real    :: kappa(nvl),avgkap(nvl)     ! variables used in kappa diagno
  real    :: exdif(nvlp1),rmsdsp        ! variables used in kappa diagno
  real    :: exsav(nvlp1)               ! variables used in kappa diagno

  logical,parameter :: kappa_diag=.FALSE.	! vertical diffusivity diagno
  real   ,parameter :: small=1.e-6
  real   ,parameter :: acurcy=2.e-6	! for 32-bit word length

  integer,parameter :: cnsv=1		! cnsv=1: conserve pot+intern.energy
! integer,parameter :: cnsv=2		! cnsv=2: conserve column height

  ex2p(arg) = p1000*(arg/cp)**(cp/rd)	! convert Pi => p
  
!sms$ignore begin

  if (vrbos) print *,'entering remap_1d, ipn=',glbl

  if (cnsv == 1) then
!   pk1col(:) = exold(:)*prcol(:)	! old pk1 = p^{kappa+1)
    pk1col(:) = exold(:)*prcol(:)/cp	! old pk1 = p^{kappa+1)
  else if (cnsv == 2) then
    pk1col(:) = exold(:)		! old pk1 = p^k (Exner fcn)
  end if

  exnew(:) = excol(:)

  if (.not.pure_sig) then
    if (intfc_smooth > 0.) then
      call inflate (prcol(1), exnew, vrbos, glbl)
    end if

! --- (optional:) retard restoration to reduce overshooting
    if (its > 0) then
      do k=1,nvlp1
        exnew(k) = slak*exnew(k) + (1. - slak)*exold(k)
      end do
    end if

! --- limit interface migration through the 'miglim' parameter
    if (miglim > 0) then
      exwrk(:)=exnew(:)
      do k=1+miglim,nvlp1-miglim
!       exnew(k)=max(exold(k+miglim),min(exnew(k),exold(k-miglim)))
        exnew(k)=max(.0001*exold(k)+.9999*exold(k+miglim),		&
                 min(.0001*exold(k)+.9999*exold(k-miglim),exnew(k)))
      end do
      do k=2,miglim
        exnew(      k)=min(exold(    1),exnew(      k))
        exnew(nvlp1-k)=max(exold(nvlp1),exnew(nvlp1-k))
      end do
      if (vrbos) print '(i7,a,i3/(i18,2f9.2,f9.3))',glbl,		&
        '  (remap_1d)   exold    exnew     clip    MIGLIM =',		&
        miglim,(k,exold(k),exwrk(k),exnew(k)-exwrk(k),k=1,nvlp1)
    end if
  end if	! if (.not.pure_sig)

! --- update pressure, layer thickness, pressure^(kappa+1)

  do k=nvlp1,1,-1
    prnew(k) = ex2p(exnew(k))
  end do

  do k=1,nvl
    dpcol(k) = prnew(k) - prnew(k+1)
  end do

  if (cnsv == 1) then
!   pk1new(:) = exnew(:)*prnew(:)	! new pk1 = p^(kappa+1)
    pk1new(:) = exnew(:)*prnew(:)/cp	! new pk1 = p^(kappa+1)
  else if (cnsv == 2) then
    pk1new(:) = exnew(:)		! new pk1 = p^k (Exner fcn)
  end if

  prnew (nvlp1) = prcol (nvlp1)		! safeguard against roundoff error
  pk1new(nvlp1) = pk1col(nvlp1)		! safeguard against roundoff error
  prnew (1) = prcol (1)			! safeguard against roundoff error
  pk1new(1) = pk1col(1)			! safeguard against roundoff error

! --- interface movement spawns vertical advection of dependent variables.
! --- we have 3 advection/remapping choices: PCM,PLM,PPM (see subr. pcwise)

  if (vrbos) print *,'now advecting -u- field'
  call pcwise(nvl,prcol,ucol,prnew,unew,miglim,vrbos,glbl)		! u vel
  if (vrbos) print *,'now advecting -v- field'
  call pcwise(nvl,prcol,vcol,prnew,vnew,miglim,vrbos,glbl)		! v vel

  do n=2,ntra			! all tracers except pot.temp
  if (vrbos) print *,'now advecting tracer',n
   call pcwise(nvl,prcol,trcol(1,n),prnew,trnew(1,n),miglim,vrbos,glbl)
  end do

! --- now advect pot.temp (k=1)
  if (vrbos) print *,'now advecting tracer 1 (theta)'
  call pcwise(nvl,pk1col,trcol,pk1new,trnew,miglim,vrbos,glbl)		! theta

  if (.not.pure_sig) then
! --- redistribute theta among neighboring layers to help them stay on target.
! --- this is to counteract a comput.mode associated with vertical advection

    do k=nvl,3,-1
      dplo = max(pk1new(k-1)-pk1new(k  ), pk1new(1)*small)
      dpup = max(pk1new(k  )-pk1new(k+1), pk1new(1)*small)
      tha = trnew(k-1,1)
      thb = trnew(k  ,1)
      devilo = (tha-targt(k-1))*dplo
      deviup = (thb-targt(k  ))*dpup
      if (deviup > 0. .and. devilo < 0.) then
        trnew(k-1,1) = tha + min(deviup,-devilo)/dplo
        trnew(k  ,1) = thb - min(deviup,-devilo)/dpup
      else if (deviup < 0. .and. devilo > 0.) then
        trnew(k-1,1) = tha - min(devilo,-deviup)/dplo
        trnew(k  ,1) = thb + min(devilo,-deviup)/dpup
      end if

      if (vrbos .and. deviup*devilo < 0.) then   
        write (6,'(a,i8,i4,2(3x,a,2f9.4))') 'ipn,k =',glbl,k,		&
             'targ dev''n',tha-targt(k-1),thb-targt(k),'cut to',	&
             trnew(k-1,1)-targt(k-1),trnew(k,1)-targt(k)
      end if
    end do

! --- fill massless cells with data from mass-containing layer below

    if (thktop == 0.) then
      do k=3,nvl
        qq = 1./(dpcol(k)+small)
        ucol(k) = (ucol(k)*dpcol(k) + ucol(k-1)*small)*qq
        vcol(k) = (vcol(k)*dpcol(k) + vcol(k-1)*small)*qq
        trnew(k,2:ntra) = (trnew(k  ,2:ntra)*dpcol(k)			&
                         + trnew(k-1,2:ntra)*small)*qq
      end do
    end if
  end if			! hybrid-isentropic option

! --- column integrals (colin/clout) are for diagnostic purposes only
  colin = 0.
  clout = 0.
  colscl = 0.
  do k=1,nvl
    colin = colin + trcol(k,1)*(pk1col(k)-pk1col(k+1))
    clout = clout + trnew(k,1)*(pk1new(k)-pk1new(k+1))
    colscl = colscl + abs(trcol(k,1)*(pk1col(k)-pk1col(k+1)))
  end do

  if (abs(clout-colin) > max(1.e-33,acurcy*colscl)) then
    write (6,106) glbl,'hybgen - column intgl.error', colin,clout,	&
       (clout-colin)/max(1.e-33,abs(colin))
106 format (i8,3x,a,2es16.8,es9.1)
  end if

  !--------------------------------------------------------------------
  ! --- vertical diffusivity diagnostics (optional)

#ifdef DEBUGPRINT
  rmsdsp = 0.
  do k=1,nvl
    rmsdsp = rmsdsp + (exnew(k+1)-excol(k+1))**2
  end do
  rmsdsp = sqrt(rmsdsp/float(nvl))

  if (kappa_diag .and. rmsdsp > .01) then
    thwrk(:) = trnew(:,1)
    exwrk(:) = exnew(:)
    call restp_1d (thwrk, exwrk, nvl, trnew, exnew, targt, nvl, vrbos, glbl)
    exdif(    1) = 0.
    exdif(nvlp1) = 0.
    do k=2,nvl
      exdif(k) = (exnew(k)-exsav(k))/dt
    end do

    call diagkp (nvl, exsav, exdif, targt, kappa, .false., glbl)
!!!    call diagkp(nvl,exsav,exdif,targt,kappa,.false.,glbl)
    if (minval(kappa) < -1.e-2) then
!!!    if (maxval(kappa) >  1.e-2) then
      call diagkp (nvl, exsav, exdif, targt, kappa, .true., glbl)

      write (6,'(a,2i5,5x,a/12(f6.1))') 'its,ipn =',its,glbl,		&
           'vert.diffusivity (cm^2/s):',(1.e4*kappa(k),k=2,nvl-1)
!!!     print *,'maxval:',maxval(kappa)
      print *,'minval:',minval(kappa)
    end if
    avgkap(:) = avgkap(:) + kappa(:)
  end if			!  diffusivity diagnostics

  if (kappa_diag .and. mod(its,15) == 0) then
    write (6,'(a,i5,5x,a/12(f6.1))') 'its =',its,			&
         'avg.vert.diffusivity (cm^2/s):',(1.e4*avgkap(k)/nip,k=2,nvl-1)
  end if

#endif

  ucol(:) = unew(:)
  vcol(:) = vnew(:)
  trcol(:,:) = trnew(:,:)
  excol(:) = exnew(:)
  prcol(:) = prnew(:)

  do k=1,nvl
    if (excol(k) < excol(k+1)) then
      write (6,99) its,glbl,deg_lat(ipn),deg_lon(ipn),'o u t p u t  profile:'
99    format ('its,ipn=',i6,i8,'  lat/lon=',2f7.1,'  hybgen  ',a/	&
              '(4-line groups: pressure, exn.fcn, theta, target)')
      do k2=1,nvl,10
        write (6,100) (prcol(k1)  ,k1=k2,min(nvlp1,k2+10) )
        write (6,102) (excol(k1)  ,k1=k2,min(nvlp1,k2+10) )
        write (6,101) (trcol(k1,1),k1=k2,min(nvl  ,k2+9 ) )
        write (6,101) (targt(k1)  ,k1=k2,min(nvl  ,k2+9 ) )
        write (6,101)
      end do
!100  format (-2p,11f7.1)
!102  format (   11f7.1)
 100  format (-2p,11f7.2)
 102  format (   11f7.2)
 101  format (5x,10f7.2)
      write (6,103) 'nonmonotonic Exner function on return at ipn,k =',	&
         glbl,k,excol(k),excol(k+1)
103   format (a,i8,i4,2f9.2)
!   stop '(error: non-monotonic Exner fcn)'
      excol(k+1) = excol(k)
    end if
  end do

#ifdef DEBUGPRINT
  if (vrbos) then
    write (6,99) its,glbl,deg_lat(ipn),deg_lon(ipn),'o u t p u t  profile:'
    do k2=1,nvl,10
      write (6,100) (prcol(k1)  ,k1=k2,min(nvlp1,k2+10) )
      write (6,102) (excol(k1)  ,k1=k2,min(nvlp1,k2+10) )
      write (6,101) (trcol(k1,1),k1=k2,min(nvl  ,k2+9 ) )
      write (6,101) (targt(k1)  ,k1=k2,min(nvl  ,k2+9 ) )
      write (6,101)
    end do
  end if
#endif

  if (vrbos) print *,'...exiting remap_1d'
  return

!sms$ignore end
  end subroutine remap_1d


  subroutine restp_1d(thold,pkold,kold,thnew,pknew,targt,knew,vrbos,ipn)
  
! --- convert a stairstep (i.e., piecewise constant) theta profile into a
! --- stairstep profile constrained to have prescribed theta ('targt') steps.

! --- input  variables: thold,pkold,targt,kold,knew,vrbos
! --- output variables: thnew,pknew

  use module_constants, only: cp,rd,p1000
  
  implicit none
  integer,intent(IN)  :: kold,knew,ipn
  real,   intent(IN)  :: thold(kold),pkold(kold+1),targt(knew)
  logical,intent(IN)  :: vrbos
  real,   intent(OUT) :: thnew(knew),pknew(knew+1)

  integer k,ko
  real oldth(kold)
  real cloutt,colint,colscl,pinteg,tha,thb,ex2p,arg
  real, parameter :: acurcy=1.e-6	! for 32-bit word length
  
  ex2p(arg)=p1000*(arg/cp)**(cp/rd)	!  convert Pi => p (mb)
  
!sms$ignore begin

  if (vrbos) then
    write (6,101) ipn,							&
    'restp1 -- input profile:   theta    thknss    press     p^kap',	&
    ex2p(pkold(1)),pkold(1),						&
    (k,thold(k),ex2p(pkold(k))-ex2p(pkold(k+1)),ex2p(pkold(k+1)),	&
    pkold(k+1),k=1,kold)
 101 format (i8,4x,a/54x,f9.1,f11.3/(33x,i3,f9.3,2f9.1,f11.3))
  end if
  
! --- remove theta inversions from input profile
  oldth(kold)=thold(kold)
  do k=kold,2,-1
    oldth(k-1)=min(oldth(k),thold(k-1))
  end do
  
  thnew(:)=targt(:)
  thnew(   1)=min(oldth(1),oldth(kold),targt(   1))
  thnew(knew)=max(oldth(1),oldth(kold),targt(knew))
  pknew(     1)=pkold(     1)
  pknew(knew+1)=pkold(kold+1)
  
! --- column integrals (colin/clout) are computed for diagnostic purposes only
  
  cloutt=0.
  colint=0.
  colscl=0.
  do k=1,kold
    colint=colint+oldth(k)*(pkold(k+1)-pkold(k))
    colscl=colscl+abs(oldth(k)*(pkold(k+1)-pkold(k)))
  end do
  
! --- find interface pknew(k+1) separating layers k and k+1 by requiring
! --- that integral over pk*d(theta) from thnew(k) to thnew(k+1) be preserved.
  
  ko=1
  do k=1,knew-1
    pinteg=0.
    thb=thnew(k)
  5 tha=thb
    thb=min(thnew(k+1),max(thnew(k),oldth(ko)))
    pinteg=pinteg+pkold(ko)*(thb-tha)
    if (oldth(ko) < thnew(k+1)) then
      if (ko < kold) then
        ko=ko+1
        go to 5
      end if
      tha=thb
      thb=thnew(k+1)
      pinteg=pinteg+pkold(kold+1)*(thb-tha)
    end if
    pknew(k+1)=pknew(k)
    if (thnew(k+1) > thnew(k)) pknew(k+1)=pinteg/(thnew(k+1)-thnew(k))
    cloutt=cloutt+thnew(k)*(pknew(k+1)-pknew(k))
  enddo
  
  cloutt=cloutt+thnew(knew)*(pknew(knew+1)-pknew(knew))
  if (abs(cloutt-colint) > max(1.e-33,acurcy*colscl)) then
    write (6,100) ipn,'restp1 - column intgl.error',			&
    colint,cloutt,(cloutt-colint)/max(1.e-33,abs(colint))
 100 format (i8,3x,a,2es16.8,es9.1)
  end if
  
  if (vrbos) then
    write (6,101) ipn,					     		&
    'restp1 -- outpt profile:   theta    thknss    press     p^kap',	&
    ex2p(pknew(1)),pknew(1),						&
    (k,thnew(k),ex2p(pknew(k))-ex2p(pknew(k+1)),ex2p(pknew(k+1)),	&
    pknew(k+1),k=1,knew)
  end if
  
  return

!sms$ignore end
  end subroutine restp_1d


  subroutine inflate(psurf,exner,vrbos,ipn)

! --- <><><><><><<><><><><><><><><><><><><><><><><>
! --- hybridization (inflation of massless layers)
! --- <><><><><><<><><><><><><><><><><><><><><><><>

  use module_constants, only: cp,rd,p1000,dpsig
  use module_control  , only: nvlp1,kbl
  use fimnamelist     , only: nvl,PrintIpnDiag,thktop,ptop
 
  implicit none
  integer, intent(IN)       :: ipn		! icos index
  logical, intent(IN)       :: vrbos
  real   , intent(IN)       :: psurf		! surface pressure
  real   , intent(INOUT)    :: exner(nvlp1)	! Exner function

  real    :: dp0,dp1,dpsum,arg,arg1,exold,q,extop
  real    :: dex2dp,dp2dex,shrink,top,cushn,taper
  integer :: last,k,k1,k2
  real,parameter :: thin=rd/p1000	! approx. 1 Pa in Exner fcn units

! taper(k1,k2)=1./float(k2-k1)
  taper(k1,k2)=1./(1.+(.25*(k2-k1))**2)
! taper(k1,k2)=1./(1.+(.5*(k2-k1))**2)
! dex2dp(arg)=p1000*arg/rd		!  convert d Pi => dp (for p ~ p1000)
  dp2dex(arg)=rd*arg/p1000		!  convert dp => d Pi (for p ~ p1000)

! --- simplified cushion function suitable for nonnegative first argument
  cushn(arg,arg1)=.25*min(arg,2.*arg1)**2/arg1+max(arg,2.*arg1)-arg1
! cushn(arg,arg1)=arg1			! no cushn

!sms$ignore begin

  if (vrbos) print *,'entering inflate, ipn=',ipn

  dpsum=0.
  extop=cp*(ptop/p1000)**(rd/cp)
  top=.2		! sigma layers are shrunk in proportion to (psurf-top) 
  shrink=min(1.,(psurf/p1000-top)/(1.-top))
  last=1

! --- step 1: inflate layers from bottom up

  do k=1,nvl-1

! --- set lower limits for layer thknss (dp0)
! --- and upper limits for upper intfc pressure (dpsum)

   dp0=dp2dex(dpsig(k))*shrink	!  convert to Exner fcn, shrink abv mountains
   dpsum=dpsum+dp0		!  minimum cumulative distance to ground
   exold=exner(k+1)
!  if (exner(k)-dp0 < extop) then
!   print '(2a,i8,i4,a,f9.2)','warning: stack of min.thknss layers',	&
!    ' extends past top of atmosphere at ipn,k=',ipn,k,'  extop=',extop
!  end if

   if (k <= kbl) then
! --- maintain -kbl- fixed-depth layers near surface
    dp1=dp0
    last=k+1
    exner(k+1)=exner(k)-dp1

   else				!  k > kbl
    if (k > last) then		!  out of sigma domain -> reduce dp0
!     dp0=dp0/min(20.,1.25*2.**(k-last))
     dp0=dp0/min(10.,1.25*2.**(k-last))
!!   dp0=dp0/min(100.,12.5*2.**(k-last))
!!   dp0=dp0/min(100.,12.5*float(k-last)**2)
!!   dp0=dp0/min(10. ,1.25*float(k-last)**2)
! next line sets 1 hPa as minimum thickness in isentropically resolved layers
!    26 Oct 2010 - Rainer, Stan
     dp0 = max(dp0,dp2dex(100.))
    end if
    dp1=cushn(max(0.,exner(k)-min(exold,exner(1)-dpsum)),dp0)
    if (exner(k)-dp1 < exold-thin) then	!  inflation required
     if (k == last) last=k+1			!  still in sigma domain
     exner(k+1)=max(extop,exner(k)-dp1)
    end if
   end if			!  k < or > kbl

   if (vrbos .and. exold.ne.exner(k+1)) then
    write (6,105) ipn,k+1,exold,exner(k+1),			&
     max(0.,exner(k)-min(exold,exner(1)-dpsum)),dp0,dp1
 105 format (i8,'  k=',i3,'  exn.',f8.2,'  =>',f8.2,            &
      '  arg1/2,cush =',3f8.2)
   end if

  end do			!  k loop

! --- step 2: inflate layers from top down

  if (thktop > 0.) then
   do k=nvl,2,-1
    q=exner(k)/cp
    q=q*q*sqrt(q)		!  (p/p0)**(5/7) = (Exn/cp)**(5/2)
    dp0=dp2dex(thktop)*taper(k,nvl)/q
    exold=exner(k)
    dp1=cushn(max(0.,exold-exner(k+1)),dp0)
    if (exner(k+1)+dp1 > exold+thin) then
     exner(k)=exner(k+1)+dp1

     if (vrbos) then
      write (6,105) ipn,k,exold,exner(k),max(0.,exold-exner(k+1)),dp0,dp1
     end if
    
    end if
    if (exner(k) < exner(k-1)-dp0) exit
   end do			!  k loop
  end if			!  thktop > 0

  if (vrbos) print *,'...exiting inflate'
  return

!sms$ignore end
  end subroutine inflate


  subroutine pcwise(kk,xold,yold,xnew,ynew,miglim,vrbos,ipn)

!-- 1-dim PCM/PLM/PPM transport routine, extracted from HYCOM's fct3d.f.
!-- if miglim.ne.1, |xold-xnew| can exceed cell width (i.e. no CFL constraints)

!-- there are 3 advection/remapping choices:
!-- remap_optn=1			! PCM (piecewise constant)
!-- remap_optn=2			! PLM (piecewise linear)
!-- remap_optn=3			! PPM (piecewise parabolic)

  use fimnamelist,only: remap_optn
  implicit none
  integer  ,intent(IN)  :: kk		! vert.dimension
  integer  ,intent(IN)  :: miglim	! intfc migration limit (<0: no limit)
  logical  ,intent(IN)  :: vrbos	!  switch for 'verbose' mode
  integer  ,intent(IN)  :: ipn		! grid cell number

!-- xold/new	- old/new cell boundaries
!-- fldold/new	- mixing ratio of dep.variable before and after transport

    real   ,intent(IN)  :: xold(kk+1),xnew(kk+1),yold(kk)
    real   ,intent(OUT) :: ynew(kk)

    real*8 zold(kk+1),znew(kk+1),delz(kk+1),fco(kk),fcn(kk),		&
         vertfx(kk+1),vertdv(kk),fldold(kk),fldnew(kk)
    real*8 a(kk),b(kk),c(kk),dx,fcdx,yl,yr,colscl
    real*8 amount,bfore,after,slab,dslab
    integer k,lyr
    real*8, parameter :: athird=1./3.
    real*8, parameter :: small=1.e-22
    real*8, parameter :: acurcy=1.e-11

!sms$ignore begin

!-- make sure -x- increases with -k-. change sign of -xold/new- if necessary
    if (xold(1) < xold(kk+1)) then
      zold(:)=xold(:)
      znew(:)=xnew(:)
    else
      zold(:)=-xold(:)
      znew(:)=-xnew(:)
    end if
    delz(:)=znew(:)-zold(:)
    fldold(:)=yold(:)

    if (vrbos) then
     write (*,100) ipn,'(pcwise) input profile:',			&
      'old_indep   old_binsz  lowr.displ  uppr.displ   old_depnd',	&
      (k,zold(k),zold(k+1)-zold(k),delz(k),delz(k+1),fldold(k),k=1,kk),	&
      kk+1,zold(kk+1)
    end if
 100 format (i8,3x,a/9x,a/(i6,f12.1,3f12.2,es12.3))

!-- deduce old and new cell width from -zold,znew-
    do 15 k=1,kk
    fco(k)=max(dble(0.),zold(k+1)-zold(k))
15  fcn(k)=max(dble(0.),znew(k+1)-znew(k))

    bfore=0.
    colscl=0.
    do k=1,kk
      bfore=bfore+fldold(k)*fco(k)
      colscl=colscl+abs(fldold(k)*fco(k))
    end do
    fldnew(:)=fldold(:)

!-- start by filling zero-width cells with data from neighboring cells

    do 17 k=kk-1,1,-1
17  fldnew(k)=(fldnew(k)*fco(k)+fldnew(k+1)*small)			&
             /(          fco(k)+            small)
    do 18 k=2,kk
18  fldnew(k)=(fldnew(k)*fco(k)+fldnew(k-1)*small)			&
             /(          fco(k)+            small)

!-- fit 0th, 1st, or 2nd deg. polynomial to -fldnew- in each cell
    a(1 )=fldnew(1 )
    b(1 )=0.
    c(1 )=0.
    a(kk)=fldnew(kk)
    b(kk)=0.
    c(kk)=0.

    do 16 k=2,kk-1

    if (remap_optn.eq.1) then
!-- piecewise constant method:
      a(k)=fldnew(k)
      b(k)=0.
      c(k)=0.

    else if (remap_optn.eq.2) then
!-- piecewise linear method:
!-- fit linear function a+bx to tracer in each cell (-.5 < x < +.5)
      a(k)=fldnew(k)
      b(k)=0.
      if (fldnew(k) <= min(fldnew(k-1),fldnew(k+1)) .or.		&
          fldnew(k) >= max(fldnew(k-1),fldnew(k+1))) then
        b(k)=0.
      else if ((fldnew(k+1)-fldnew(k-1))*(fldnew(k-1)+fldnew(k+1)	&
        -2.*fldnew(k)) > 0.) then
        b(k)=fldnew(k)-fldnew(k-1)
      else
        b(k)=fldnew(k+1)-fldnew(k)
      end if
      c(k)=0.

    else if (remap_optn.eq.3) then
!-- piecewise parabolic method:
!-- construct parabola  a+bx+cx^2  whose integral over [-.5,+.5] equals
!-- fldnew(k) and which passes though points yl,yr at [-.5,+.5] resp.
!!    yl=.5*(fldnew(k-1)+fldnew(k))
!!    yr=.5*(fldnew(k+1)+fldnew(k))
      yl=(max(small,fco(k-1))*fldnew(k)+max(small,fco(k))*fldnew(k-1))/	&
         (max(small,fco(k-1))          +max(small,fco(k)))
      yr=(max(small,fco(k+1))*fldnew(k)+max(small,fco(k))*fldnew(k+1))/	&
         (max(small,fco(k+1))          +max(small,fco(k)))
      a(k)=1.5*fldnew(k)-.25*(yl+yr)
      b(k)=yr-yl
      c(k)=6.*(.5*(yl+yr)-fldnew(k))
      if (abs(yr-yl) < 6.*abs(.5*(yl+yr)-fldnew(k))) then
!-- apex of parabola lies inside interval [-.5,+.5], implying an over-
!-- or undershoot situation. change curve to prevent over/undershoots.
        if (abs(yr-yl) > 2.*abs(.5*(yl+yr)-fldnew(k))) then
!-- put apex of parabola on edge of interval [-.5,+.5]
          if ((yr-yl)*(.5*(yl+yr)-fldnew(k)) > 0.) then
!-- apex at x=-.5
            a(k)=.25*(3.*fldnew(k)+yl)
            c(k)=3.*(fldnew(k)-yl)
            b(k)=c(k)
          else
!-- apex at x=+.5
            a(k)=.25*(3.*fldnew(k)+yr)
            c(k)=3.*(fldnew(k)-yr)
            b(k)=-c(k)
          end if
        else			!  -1/6 < x < +1/6
!-- moving apex won't help. replace parabola by constant.
          a(k)=fldnew(k)
          b(k)=0.
          c(k)=0.
        end if
      end if
    else
      print *,'subr.pcwise -- illegal remap_optn=',remap_optn
      stop '(pcwise)'
    end if
16  continue

!-- get flux by summing -fldnew- over upstream slab of thickness -delz-

    if (miglim.eq.1) then
      do 27 k=2,kk
      slab=0.
      amount=0.
      vertfx(k)=0.
      if (delz(k) > 0.) then			! interface moves in +k dir.
        lyr=k
        if (fco(lyr) > 0.) then
          dslab=min(slab+fco(lyr), delz(k))	&
               -min(slab         , delz(k))
          dx=dslab/fco(lyr)
          fcdx=a(lyr)				&
              +b(lyr)*.5*(dx-1.)		& !  not needed in pcm
              +c(lyr)*(.25-dx*(.5-dx*athird))	  !  not needed in pcm,plm
          amount=amount+fcdx*dslab
          slab=slab+dslab
        end if
      else if (delz(k) < 0.) then		! interface moves in -k dir.
        lyr=k-1
        if (fco(lyr) > 0.) then
          dslab=min(slab+fco(lyr),-delz(k))	&
               -min(slab         ,-delz(k))
          dx=dslab/fco(lyr)
          fcdx=a(lyr)				&
              +b(lyr)*.5*(1.-dx)		& !  not needed in pcm
              +c(lyr)*(.25-dx*(.5-dx*athird))	  !  not needed in pcm,plm
          amount=amount+fcdx*dslab
          slab=slab+dslab
        end if
      end if
      if (slab.ne.0.) vertfx(k)=-delz(k)*amount/slab
27    continue

    else				! miglim .ne. 1
      do 22 k=2,kk
      slab=0.
      amount=0.
      vertfx(k)=0.
      if (delz(k) > 0.) then			! interface moves in +k dir.
        lyr=k-1
24      lyr=lyr+1
        if (slab >= delz(k)) goto 23
        if (fco(lyr) > 0.) then
          dslab=min(slab+fco(lyr), delz(k))	&
               -min(slab         , delz(k))
          dx=dslab/fco(lyr)
          fcdx=a(lyr)				&
              +b(lyr)*.5*(dx-1.)		& !  not needed in pcm
              +c(lyr)*(.25-dx*(.5-dx*athird))	  !  not needed in pcm,plm
          amount=amount+fcdx*dslab
          slab=slab+dslab
        end if
        if (lyr < kk) go to 24
      else if (delz(k) < 0.) then		! interface moves in -k dir.
        lyr=k
25      lyr=lyr-1
        if (slab >= -delz(k)) goto 23
        if (fco(lyr) > 0.) then
          dslab=min(slab+fco(lyr),-delz(k))	&
               -min(slab         ,-delz(k))
          dx=dslab/fco(lyr)
          fcdx=a(lyr)				&
              +b(lyr)*.5*(1.-dx)		& !  not needed in pcm
              +c(lyr)*(.25-dx*(.5-dx*athird))	  !  not needed in pcm,plm
          amount=amount+fcdx*dslab
          slab=slab+dslab
        end if
        if (lyr > 1) go to 25
      end if
23    if (slab.ne.0.) vertfx(k)=-delz(k)*amount/slab
22    continue
    end if

    vertfx(   1)=0.			!  don't allow flux through lower bdry
    vertfx(kk+1)=0.			!  don't allow flux through upper bdry
    do 26 k=1,kk
26  vertdv(k)=vertfx(k+1)-vertfx(k)

    if (vrbos) write (*,'(a/(i3,4es12.3))')				&
     'pcwise:   flux  flx.div/thk    old_thk     new_thk',		&
      (k,vertfx(k),vertdv(k)/max(small,fcn(k)),fco(k),fcn(k),k=1,kk),	&
       kk+1,vertfx(kk+1)

    after=0.
    do 4 k=1,kk
    amount=fldnew(k)*fco(k)-vertdv(k)
    fldnew(k)=(amount+fldnew(k)*small)/(fcn(k)+small)

!   if (miglim.eq.1 .and. k.gt.1 .and. k.lt.kk) then
!    if (fldnew(k).lt.min(fldold(k-1),fldold(k),fldold(k+1))		&
!      -max(1.,abs(fldnew(k)))*acurcy .or.				&
!        fldnew(k).gt.max(fldold(k-1),fldold(k),fldold(k+1))		&
!      +max(1.,abs(fldnew(k)))*acurcy)					&
!    print '(i8,i4,a,f11.2,a,3es11.3)',ipn,k,' (pcwise)',		&
!      fldnew(k),' outside range',fldnew(k)-fldold(k-1),		&
!      fldnew(k)-fldold(k),fldnew(k)-fldold(k+1)
!   end if

    after=after+amount
 4  continue

    if (abs(bfore-after) > max(1.d-33,acurcy*colscl)) then
      write (*,104) ipn,'pcwise - column intgl.error',bfore,after,	&
      (after-bfore)/max(1.d-33,abs(bfore))
    end if
 104 format (i8,3x,a,2es16.8,es9.1)

    if (vrbos) then
     write (*,100) ipn,'(pcwise) output profile:',			&
      'new_indep   new_binsz  lowr.displ  uppr.displ   new_depnd',	&
      (k,znew(k),znew(k+1)-znew(k),delz(k),delz(k+1),fldnew(k),k=1,kk),	&
      kk+1,znew(kk+1)
    end if

    ynew(:)=fldnew(:)
    return

!sms$ignore end
  end subroutine pcwise


   subroutine equilb(thet,pres,kk,slak,vrbos,ipn)

! --- expand thin layers at the expense of thick layers above and below.
! --- do this without changing theta

   use module_constants,only: rd,p1000

   implicit none
   integer,intent(IN) :: kk,ipn		! no. of layers, test point location
   logical,intent(IN) :: vrbos		! if true, print results at test point
   real,intent(IN)    :: thet(kk)	! pot.temp. in layers
   real,intent(INOUT) :: pres(kk+1)	! Exner fcn on interfaces
   real,intent(IN)    :: slak		! retardation coefficient
   integer k,k1,k2,ncount
   real dp1,dp2,dp3,dp4,dp5,th1,th2,th3,th4,th5,dis3,dis4,	&
        ratio,goal,pnew(kk+1)
   logical event
   real,parameter :: thin=rd/p1000	! approx. 1 Pa in Exner fcn units

!sms$ignore begin

   if (vrbos) then
    write (6,99) ipn,'  equilb input profile:'
    do k2=1,kk,10
     write (6,100) (pres(k1),k1=k2,min(kk+1,k2+10) )
     write (6,101) (thet(k1),k1=k2,min(kk  ,k2+9 ) )
     write (6,100)
    end do
 99  format ('ipn=',i8,a,i9,a)
 100 format (11f7.1)
 102 format (11i7)
 101 format (5x,10f7.2)
   end if

! --- scenario 1: sequence of 5 thin-thick-thin-thick-thin layers

   pnew(:)=pres(:)
   ncount=0
   do 1 k=3,kk-2
    if (pnew(k-2) > pnew(1)-thin) go to 1
    dp1=pnew(k-2)-pnew(k-1)
    dp2=pnew(k-1)-pnew(k  )
    dp3=pnew(k  )-pnew(k+1)
    dp4=pnew(k+1)-pnew(k+2)
    dp5=pnew(k+2)-pnew(k+3)
    th1=thet(k-2)
    th2=thet(k-1)
    th3=thet(k  )
    th4=thet(k+1)
    th5=thet(k+2)
! --- look for small dp1,dp3,dp5 in combination with large dp2,dp4
    if (dp2 > dp1 .and. dp4 > dp5) then
     goal=.5*(dp3+min(dp2,dp4))		!  desired thknss of lyr 3
     if (dp2 > goal .and. dp4 > goal) then
! --- thin-thick-thin-thick-thin combination found -> inflate lyr 3
      dis3=min(dp2-goal,goal-dp3) * slak
      dis4=min(dp4-goal,goal-dp3) * slak
      if (th3 > th2 .and. th4 > th3) then
! --- theta conservation requires  dis3*(th3-th2)+dis4*(th4-th3)=0
       ratio=(th4-th3)/(th3-th2)
       if (ratio > 1.) then
        dis4=dis3/ratio
       else
        dis3=dis4*ratio
       end if
      end if
! --- ready to expand middle layer
      pnew(k  )=pnew(k  )+dis3
      pnew(k+1)=pnew(k+1)-dis4

      if (vrbos) then
       write (6,'(a,5f8.2)') 'thknss quintuplet',dp1,dp2,		&
        dp3,dp4,dp5, '          becomes',dp1,pnew(k-1)-pnew(k),		&
         pnew(k)-pnew(k+1),pnew(k+1)-pnew(k+2),dp5,			&
          'orig dis3,dis4, mod dis3,dis4,ratio =',			&
           min(dp2-goal,goal-dp3),min(dp4-goal,goal-dp3),		&
            dis3,dis4,ratio
      end if
      ncount=ncount+1
     end if
    end if
 1 continue

! --- scenario 2: sequence of 3 thin-thick-thin layers

!  do 2 k=2,kk-1
!   if (pnew(k) > pnew(1)-thin) go to 2
!   dp2=pnew(k-1)-pnew(k  )
!   dp3=pnew(k  )-pnew(k+1)
!   dp4=pnew(k+1)-pnew(k+2)
!   th2=thet(k-1)
!   th3=thet(k  )
!   th4=thet(k+1)
! --- look for small dp2,dp4 in combination with large dp3
!   if (dp3 > dp2 .and. dp3 > dp4) then
!    goal=.5*(dp3+max(dp2,dp4))		!  desired thknss of lyr 3
!    if (dp2 < goal .and. dp4 < goal) then
! --- thin-thick-thin combination found -> deflate lyr 3
!     dis3=min(goal-dp2,dp3-goal) * slak
!     dis4=min(goal-dp4,dp3-goal) * slak
!     if (th3 > th2 .and. th4 > th3) then
! --- theta conservation requires  dis3*(th3-th2)+dis4*(th4-th3)=0
!      ratio=(th4-th3)/(th3-th2)
!      if (ratio > 1.) then
!       dis4=dis3/ratio
!      else
!       dis3=dis4*ratio
!      end if
!     end if
! --- ready to shrink middle layer
!     pnew(k  )=pnew(k  )-dis3
!     pnew(k+1)=pnew(k+1)+dis4

!     if (vrbos) then
!      write (6,'(a,i3,a,3f9.3/a,3f9.3)') 'k=',k,' thknss triple',	&
!       dp2,dp3,dp4,'            becomes',pnew(k-1)-pnew(k),		&
!        pnew(k)-pnew(k+1),pnew(k+1)-pnew(k+2)
!      write (6,'(a,3es11.3)') 'pres.chgs.,ratio=',-dis3,dis4,ratio
!     end if

!     ncount=ncount+1
!    end if
!   end if
!2 continue

   event = ncount>kk/3			!  find interesting cases

   do k=1,kk
    if (pnew(k+1) > pnew(k)+thin) then
     event=.true.
     write (*,'(a,i3)')						&
     'error: nonmonotonic pressure on return from equilb, k=',k
    end if
   end do

   if (event .or. vrbos) then
    write (6,99) ipn,'  equilb input profile:'
    do k2=1,kk,10
     write (6,100) (pres(k1),k1=k2,min(kk+1,k2+10) )
     write (6,101) (thet(k1),k1=k2,min(kk  ,k2+9 ) )
     write (6,100)
    end do
    write (6,99) ipn,'  equilb output profile:',ncount,' inflations'
    do k2=1,kk,10
     write (6,100) (pnew(k1),k1=k2,min(kk+1,k2+10) )
     write (6,102) (int(1000.*(pnew(k1)-pres(k1))),k1=k2,min(kk+1,k2+10) )
     write (6,101) (thet(k1),k1=k2,min(kk  ,k2+9 ) )
     write (6,100)
    end do
   end if

   pres(:)=pnew(:)
   return

!sms$ignore end
end subroutine equilb
   

subroutine diagkp(kk,z,zdot,theta,kappa,vrbos,ipn)

! --- diagnose vertical mixing coefficient resulting from manipulation
! --- of the air column (such as vertical advection, regridding,...)

! --- approach: transform diffusion eqn for theta
! ---  (d theta/dt)_z = (d/dz) [kappa d theta/dz]
! --- into
! ---  (dz/dt)_theta = -(d/d theta) [kappa d theta/dz]
! --- ('d' = 'partial'). represent kappa in layers as averages of kappa
! --- on interfaces. solve tridiagonal system for interface kappas.
! --- kappa values returned are   l a y e r   averages.
! --- kappa is assumed zero in first and last layer (ref. McDougall & Dewar).

! --- input variables:
! --- z(kk+1)    -  interface depths (meters, massless layers allowed)
! --- zdot(kk+1) -  rate of interface movement (m/s)
! --- theta(kk)  -  buoyancy (independent variable!)

! --- output: kappa (m^2/s) in layers whose thickness exceeds 'thresh'

   implicit none
   integer i,j,k,kk,kp
   real z(kk+1),theta(kk),zdot(kk+1),kappa(kk),zmid(kk),	&
        thmid(kk),zdotm(kk),a(kk+1,kk+1),a1(kk+1,kk+1),		&
        b(kk+1),b1(kk+1),d,dpidz
   integer indx(kk+1),klist(kk),ipn
   logical vrbos
   real,parameter :: thresh=1.e-3

!sms$ignore begin

! --- weed out massless layers
   kp=0
   do 5 k=1,kk
   klist(k)=0
   if (z(k+1) < z(k)-thresh) then	!  atmospheric case: z = Exn.fcn, decreasing with k
     kp=kp+1
     klist(k)=kp
     zmid(kp)=.5*(z(k)+z(k+1))		! mid layer depth
     thmid(kp)=theta(k)
     zdotm(kp)=zdot(k)
   end if
5  continue

   if (vrbos) then
     write (*,'(a,i8,a/(5(f9.1,f6.1)))') 'ipn =',ipn,			&
      '  input profile:',(z(k),theta(k),k=1,kk),z(kk+1)
     write (*,'(2(i5,a))') kp,' non-massless layers  =>',kp-3,' unknowns'
   end if

   a=0.
   do 10 k=3,kp-1
   a(k,k+1)= (thmid(k+1)-thmid(k-1))/(zmid(k+1)-zmid(k-1))
   a(k,k-1)=-(thmid(k  )-thmid(k-2))/(zmid(k  )-zmid(k-2))
   a(k,k  )=a(k,k+1)+a(k,k-1)
10 continue

   do 11 k=3,kp-1
11 b(k)=2.*zdotm(k)*(thmid(k-1)-thmid(k))

   if (vrbos) then
     write(*,*) 'ipn =',ipn,'  input matrix * 10^3'
     do k=3,kp-1
       write(*,'(3p,20f8.2)') (a(k,j),j=3,kp-1)
     end do
   end if

   a1=a
   b1=b

   if (vrbos) write(*,*) 'ipn =',ipn,'  input rhs * 10^3'
   if (vrbos) write(*,'(3p,20f8.2)')  (b(j),j=3,kp-1)

   call ludcmp(a(3,3),kp-3,kk+1,indx,d)
   call lubksb(a(3,3),kp-3,kk+1,indx,b(3))

   if (vrbos) then
! --- did the equation solver return a credible solution?
     do 22 k=3,kp-1
     b1(k)=0.
     do 22 i=3,kp-1
22   b1(k)=b1(k)+a1(k,i)*b(i)
     write(*,*) 'ipn =',ipn,'  matrix * kappa * 10^3'
     write(*,'(3p,20f8.2)')  (b1(j),j=3,kp-1)

     do k=2,kk-1
       j=klist(k)
       if (j > 1 .and. j < kp) then
         dpidz=9.806/thmid(j)			!  exn.fcn => z conversion
         write(*,'(a,i4,6(a,f8.2))')	 				&
!!       write(*,'(a,i4,2(a,f8.1),a,es11.2,3(a,f8.1))')	 		&
         'k=',k,' z_uppr=',z(k),' zmid=',zmid(j),			&
         ' zdot(x10^3)=',1.e3*zdot(k),					&
         ' k_uppr(cm^2/s)=',1.e4*b(j)*dpidz**2,				&
         ' k_mid=',.5e4*(b(j)+b(j+1))*dpidz**2
       end if
     end do
   end if

! --- set diffusivity=0 in 1st and last layer
   b(   1)=0.
   b(   2)=0.
   b(kp  )=0.
   b(kp+1)=0.
   kappa=0.
   do k=2,kk-1
     dpidz=9.806/theta(k)			!  atmosphere
     j=klist(k)
     if (j > 1 .and. j < kp) kappa(k)=.5*(b(j)+b(j+1))*dpidz**2
   end do
   return

!sms$ignore end
 end subroutine diagkp


!! SUBROUTINE ludcmp(a,n,np,indx,d)
!!       INTEGER n,np,indx(n),NMAX
!!       REAL d,a(np,np),TINY
!!       PARAMETER (NMAX=500,TINY=1.0e-20)
!!       INTEGER i,imax,j,k
!!       REAL aamax,dum,sum,vv(NMAX)
!!       d=1.
!! 
!!       do 12 i=1,n
!!         aamax=0.
!!         do 11 j=1,n
!!           if (abs(a(i,j)) > aamax) aamax=abs(a(i,j))
!! 11      continue
!!         if (aamax == 0.) pause 'singular matrix in ludcmp'
!!         vv(i)=1./aamax
!! 12    continue
!!       do 19 j=1,n
!!         do 14 i=1,j-1
!!           sum=a(i,j)
!!           do 13 k=1,i-1
!!             sum=sum-a(i,k)*a(k,j)
!! 13        continue
!!           a(i,j)=sum
!! 14      continue
!!         aamax=0.
!!         do 16 i=j,n
!!           sum=a(i,j)
!!           do 15 k=1,j-1
!!             sum=sum-a(i,k)*a(k,j)
!! 15        continue
!!           a(i,j)=sum
!!           dum=vv(i)*abs(sum)
!!           if (dum >= aamax) then
!!             imax=i
!!             aamax=dum
!!           endif
!! 16      continue
!!         if (j.ne.imax)then
!!           do 17 k=1,n
!!             dum=a(imax,k)
!!             a(imax,k)=a(j,k)
!!             a(j,k)=dum
!! 17        continue
!!           d=-d
!!           vv(imax)=vv(j)
!!         endif
!!         indx(j)=imax
!!         if(a(j,j) == 0.)a(j,j)=TINY
!!         if(j.ne.n)then
!!           dum=1./a(j,j)
!!           do 18 i=j+1,n
!!             a(i,j)=a(i,j)*dum
!! 18        continue
!!         endif
!! 19    continue
!!       return
!! end subroutine ludcmp
!! !  (C) Copr. 1986-92 Numerical Recipes Software 'W3.
!! 
!! 
!! SUBROUTINE lubksb(a,n,np,indx,b)
!!       INTEGER n,np,indx(n)
!!       REAL a(np,np),b(n)
!!       INTEGER i,ii,j,ll
!!       REAL sum
!!       ii=0
!!       do 12 i=1,n
!!         ll=indx(i)
!!         sum=b(ll)
!!         b(ll)=b(i)
!!         if (ii.ne.0)then
!!           do 11 j=ii,i-1
!!             sum=sum-a(i,j)*b(j)
!! 11        continue
!!         else if (sum.ne.0.) then
!!           ii=i
!!         endif
!!         b(i)=sum
!! 12    continue
!!       do 14 i=n,1,-1
!!         sum=b(i)
!!         do 13 j=i+1,n
!!           sum=sum-a(i,j)*b(j)
!! 13      continue
!!         b(i)=sum/a(i,i)
!! 14    continue
!!       return
!! end subroutine lubksb
!  (C) Copr. 1986-92 Numerical Recipes Software 'W3.

end module module_hybgen
