module mdul_transects

!            Instructions for constructing transects in FIM
!            ==============================================
! 
! The centerpiece of the transport diagnostics package is an offline utility
! that constructs strings ("chains") of icosahedral cell edges approximating
! user-selected transect lines.
! 
! Sources for the chain finder utility and for plotting the resulting
! transects are in the fim/ihycomproc/bsncrop branch of svn.
! 
! Step-by-step instructions:
! 
! 1. Compile lalosecgen.F90 by issuing the command 'make lalosecgen'. The
! routine requires a grid configuration file glvl<G>.dat which can be
! found under the name glvl.dat (without the <G>) in the 'prep'
! subdirectory of a FIM run using curve=0 and the desired G level.
! 
! 2. In the namelist file lalobounds.nml, specify glvl and the end points
! as well as identifiers for up to 20 transects. Zonal transects are
! defined by a latitude 'reflat' and 2 longitude values 'longlo,longhi'.
! If longlo=longhi, the resulting transect will span the whole latitude
! circle. 
!
! 3. Execute lalosecgen.
! 
! 4. For each specified transect, lalosecgen generates 2 output files with
! names like
! 
!     G5_65N_latcirc.dat
!     G5_65N_latcirc.bin
! 
! In this particular example, G5 indicates glvl=5 and 'latcirc' is the
! arbitrarily chosen identifier mentioned earlier. The .dat file serves as
! input to FIM while the .bin file serves as input to 'chainplot', a
! plotting utility based on Ning Wang's icosphere.
! 
! 5. Cell edge information is given in G5_65N_latcirc.dat in the following
! form: 
!
!   145 pts
!         14 1 clockw  64.230  10.000        14 6 clockw  64.230  10.000
!         45 4 countr  65.276  14.717        45 5 countr  65.276  14.717
!         77 1 clockw  64.203  19.089       108 2 clockw  64.959  24.072
!        108 1 clockw  64.959  24.072       108 6 clockw  64.959  24.072
!        139 4 countr  65.546  29.308       139 5 countr  65.546  29.308
!        171 1 clockw  64.039  32.981       170 5 countr  65.966  34.758
!        202 1 clockw  64.305  38.133       201 5 countr  66.170  40.319
!        233 1 clockw  64.387  43.347       264 2 clockw  64.387  48.653
!       ....
!       ....
!       8517 6 countr  65.546 350.692      8548 2 clockw  64.959 355.928
!       8548 1 clockw  64.959 355.928      8579 2 clockw  64.203   0.911
!       8579 1 clockw  64.203   0.911      8578 4 countr  65.276   5.283
!       8578 5 countr  65.276   5.283
! 
! 
! In this example, the first (westernmost) member of the chain is edge #1 of
! icos cell 14, the 2nd one is edge #6 of cell 14, and the last one is edge
! #5 of cell 8578 (which by design is adjacent to edge #1 of cell 14). The
! numbers  64.230  10.000  indicate latitude and longitude of cell 14, and
! so forth. "clockw" and "countr" indicate whether the edges in a given cell
! are passed through in clockwise or counterclockwise direction. (Each cell
! contributes 2 edges to the chain, but the 1st edge coincides with the
! 2nd edge of the previous cell, hence is not listed to avoid duplication.)
! This information is needed to determine whether the mass flux through a
! given edge is to be added to, or subtracted from, the transport total.
! 
! 6. List the transects to be used by FIM in a namelist file called
! transects_atm_G<glvl>.nml. This file must reside in the same directory as
! the .dat files themselves. Specify the path to this directory in
! 'atm_trnsecdir' in FIMnamelist.
! 
! 
!                computing mass fluxes across transects
!                ======================================
! 
! The following additions to the FIM source code are required to compute
! time-integrated mass fluxes across transects:
! 
! 1. Add the namelist variable atm_trnsecdir to  cntl/fimnamelist.F90
! (2 places!). Specify atm_trnsecdir in FIMnamelist.
!
! 2. In cntl/module_control, add the lines
!       integer, parameter :: numsecs=10         ! max.# of transects
!       character*40 :: thruflsecs(numsecs)      ! transects for thruflow diag
!
! 3. In fim/horizontal, add atm_transects.o to FIM_HORIZONTAL_OBJS and import
! the source code atm_transects.F90
!
! 4. Insert the following in module_constants.F90 (these are nondistributed
! arrays):
!   integer    ,allocatable :: lgthset(:),cellset(:,:),edgeset(:,:)
!   character*6,allocatable :: sense(:,:)
!   real       ,allocatable :: ctrlat(:,,:),ctrlon(:,:),transp(:,:,:)
! 
! 5. Add the variable flxavg(:,:,:) in module_variables.F90 and allocate
! it in dyn_alloc.F90 the way cumuflx is allocated.
! 
! 6. Add the line 'use mdul_transects,only: trnsec_init' in dyn_init.F90
! and insert a call to trnsec_init (no arguments).  This will read in the
! desired transect(s) from the directory specified by the FIMnamelist
! variable atm_trnsecdir.
! 
! 7. In dyn_run.F90, include flxavg in 'use module_variables' and add the
! line 'use mdul_transects,only: trnsec1'. Then add a call to trnsec1 near
! the end of dyn_run:
!
!   call trnsec1(its,                &
!        nf,of,vof,                  & ! time slots for Adams Bashforth
!        adbash1,adbash2,adbash3,    & ! Adams Bashforth time dif. weights
!        flxavg,massfx)              ! time-integrated and instant. massfx
!
! This will accumulate mass fluxes exactly as subr.transp1 does, but over a
! longer time period.
! 
! 8. In output.F90, include flxavg in 'use module variables' and add the
! line 'use mdul_transects,only: atm_flxsum'. Insert a call to flxsum, to
! be executed at the end of the run:
!
!     if (its.ge.itsStart-1+nts) call atm_flxsum(itsStart-1,nts,flxavg)
!
! This will compute the mass flux across each transect, layer by layer, and
! averaged over TotalTime, the period covered by the model run. Results are
! written in ASCII form to stdout. If restarts are required, flux results
! from consecutive restart runs can be averaged off-line to obtain fluxes
! valid for the overall integration period.
!
! 9. The variable 'trform' provides the option to transform mass fluxes
! from the native vertical grid to either theta or sigma coordinates.
!
! 10. A utility atm_trnsecplot.F90 is available for extracting the pertinent
! information from stdout and plotting the layer mass fluxes.


contains
  subroutine trnsec_init

! --- read chains of edge indices needed in mass flux diagnostics.
! --- to be called by dyn_init.

  use module_control  ,only: glvl,nvl,numsecs,thruflsecs
  use fimnamelist     ,only: atm_trnsecdir
  use module_constants,only: lgthset,cellset,edgeset,sense,		&
                             ctrlat,ctrlon,transp,inv_perm
  use module_variables,only: flxavg
  use units           ,only: getunit, returnunit

  implicit none
  integer n,maxlen,lun1,lun2,len,ios
  integer,allocatable :: tmp(:)
  character*120 flnm1,flnm2
  character*2 chrglvl

  print *,'entering trnsec_init...'
  if (len_trim(atm_trnsecdir).gt.0) then
   thruflsecs(:)=' '
   write (chrglvl,'(a,i1)') 'G',glvl
   flnm1=trim(atm_trnsecdir)//'transects_atm_'//chrglvl//'.nml'
   print *,'get atm.transects for online transport diagnostics from ',	&
     trim(flnm1)
  else
   print *,'no atm.transects found for online transport diagnostics'
  print *,'... exiting trnsec_init'
   return
  end if                        ! atm_trnsecdir exists

  print *,'initializing online mass transport diagnostics...'
  allocate(lgthset(numsecs))
  lgthset(:)=0
  lun1 = getunit ()
  lun2 = getunit ()

! --- determine space requirements for storing edge chain info
  open (lun1,file=flnm1,form='formatted',action='read',iostat=ios)
  do n=1,numsecs			! look for thruflow transects
   read (lun1,'(a40)',end=2) thruflsecs(n)
   if (index(thruflsecs(n),chrglvl).eq.0) then
    print *,'cannot find ',chrglvl,' in file name ',			&
      trim(thruflsecs(n)),' ... proceeding with',n-1,' transects'
    exit
   else
    flnm2=trim(atm_trnsecdir)//trim(thruflsecs(n))
    print *,'checking ',trim(flnm2),' for length of transect'
    open (lun2,file=flnm2,form='formatted',action='read',iostat=ios)
    read (lun2,'(i5)') lgthset(n)
    close (lun2)
   end if
  end do
2 close (lun1)

  maxlen=maxval(lgthset(:))
  print *,'max length of atm.transects:',maxlen
  allocate (edgeset(maxlen,numsecs),cellset(maxlen,numsecs),		&
            sense(maxlen,numsecs),ctrlat(maxlen,numsecs),		&
            ctrlon(maxlen,numsecs),transp(nvl,maxlen,numsecs),		&
            tmp(maxlen))
  edgeset(:,:)=0
  cellset(:,:)=0
  ctrlat (:,:)=0.
  ctrlon (:,:)=0.
  transp (:,:,:)=0.
  sense  (:,:)=' '

! --- read the chosen transects

  do n=1,numsecs				! loop over transects
   if (thruflsecs(n)(1:1).ne.' ') then
    flnm2=trim(atm_trnsecdir)//thruflsecs(n)
    print '(/2a)','read cell/edge info for atm.thruflow transect ',	&
     trim(flnm2)
    open (lun2,file=flnm2,form='formatted',action='read',iostat=ios)
!SMS$SERIAL BEGIN
    read (lun2,101) (cellset(len,n),edgeset(len,n),			&
      sense(len,n),ctrlat(len,n),ctrlon(len,n),len=1,lgthset(n))
!SMS$SERIAL END
101 format (/2(i10,i2,1x,a6,2f8.3))
    close (lun2)

    print '(/2a)','(trnsec_init) atm.transects (GLOBAL indices) for ',	&
      trim(thruflsecs(n))
    print 101,(cellset(len,n),edgeset(len,n),sense(len,n),		&
     ctrlat(len,n),ctrlon(len,n),len=1,lgthset(n))

!SMS$SERIAL (<inv_perm,IN>) BEGIN
    do len=1,lgthset(n)
     tmp(len)=inv_perm(cellset(len,n))
!    print '(2(a,i8))','map',cellset(len,n),' onto local cell',		&
!      tmp(len)
    end do
    do len=1,lgthset(n)
     cellset(len,n)=tmp(len)
    end do
!SMS$SERIAL END

    print '(/2a)','(trnsec_init) atm.transects (LOCAL indices) for ',	&
      trim(thruflsecs(n))
    print 101,(cellset(len,n),edgeset(len,n),sense(len,n),		&
     ctrlat(len,n),ctrlon(len,n),len=1,lgthset(n))
   end if
  end do

  call flush(6)
  call returnunit (lun1)
  call returnunit (lun2)
! --- initialize mass flux time integral
  flxavg(:,:,:)=0.
  print *,'...exiting trnsec_init'
  return
  end subroutine trnsec_init


   subroutine atm_flxsum(nstart,nstep,massflx)

! --- sum up mass fluxes normal to a specified line of icos edges.
! --- to be called at end of run.

   use module_control  ,only: nip,npp,nvl,numsecs,			&
                              thruflsecs,ArchvTimeUnit,dt
   use fimnamelist     ,only: atm_trnsecdir
   use module_constants,only: grvity,lgthset,cellset,edgeset,sense,	&
                              ctrlat,ctrlon,transp,prox,nprox,perm
   implicit none
   integer,intent(IN) :: nstart,nstep
!SMS$DISTRIBUTE (dh,3) BEGIN
   real   ,intent(IN) :: massflx(nvl,npp,nip)	! time-integrated mass flux
!SMS$DISTRIBUTE END
   integer            :: n,i,ix,len,k
   real               :: sum0
   real               :: valu, zero=38.		! linout-related variables
   logical            :: vrbos=.false.
   real,  external    :: its2time

   call flush(6)
!  print *,'entering atm_flxsum, nstep =',nstart+nstep
   print *,'-------atm.flxsum--------atm.flxsum--------atm.flxsum--------atm.flxsum-------'
   do n=1,numsecs
    if (thruflsecs(n)(1:1).ne.' ') then
     print '(2a)','starting atmo.thruflow calculation for ',		&
       trim(thruflsecs(n))
     if (vrbos) then
      print '(2a)',trim(thruflsecs(n)),' -- cell/edge info:'
      print 100,(cellset(len,n),edgeset(len,n),sense(len,n),		&
       ctrlat(len,n),ctrlon(len,n),len=1,lgthset(n))
100   format (2(i10,i2,1x,a6,2f8.2))
     end if

     do k=1,nvl
      transp(k,:,n)=0.
      do len=1,lgthset(n)
       sum0=0.
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE(ix) REDUCTION(+:sum0)
       do i=1,nip
        if (cellset(len,n).eq.i) then
         ix=prox(edgeset(len,n),i)
! --- are we circling the cell in clock- or counterclockwise direction?
         sum0=massflx(k,edgeset(len,n),cellset(len,n))		! N
         sum0=sum0/(dt*float(nstep)*grvity*1.e6)		! Ktons/s
         if (sense(len,n).eq.'clockw') then
          if (vrbos) print  101,k,perm(i),'  exports',sum0,		&
            ' t/s across edge',edgeset(len,n),'    to',perm(ix)
         else if (sense(len,n).eq.'countr') then
          if (vrbos) print  101,k,perm(i),'  imports',sum0,		&
            ' t/s across edge',edgeset(len,n),'  from',perm(ix)
          sum0=-sum0
         else
          stop '(wrong characters in sense array)'
         end if
101      format('(flxsum)  k=',i2,'  cell',i8,a,f8.1,a,i2,a,i8)
        end if			! cell i belongs to transect
       end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END
!SMS$REDUCE (sum0,SUM)
       transp(k,len,n)=sum0
      end do			! loop over cells in transect

      if (vrbos) print 103,'(flxsum)  k=',k,				&
        'transport contribution from indiv.cell edges:'			&
        ,(cellset(len,n),edgeset(len,n),transp(k,len,n),		&
        len=1,lgthset(n))
103   format (a,i2,2x,a/(3(i11,i3,f11.1)))
     end do			! vertical loop

! --- summation over transect
     do len=2,lgthset(n)
      do k=1,nvl
       transp(k,1,n)=transp(k,1,n)+transp(k,len,n)
      end do
     end do

     do k=1,nvl
      print 102,ArchvTimeUnit,its2time(nstart+nstep),'  k=',k,		&
        trim(thruflsecs(n)),'thruflow (Ktons/s):',transp(k,1,n)
102   format (a,f11.1,a,i3,2(2x,a),f11.1)
     end do			! vertical loop

     call flush(6)
     do k=nvl,1,-1
      call linout(zero,'o',-1)
      valu=transp(k,1,n)*2.e-3 + zero
      call linout(valu,'X', 1,k)
     end do
     call flush(6)

! --- summation over layers
     do k=2,nvl
      transp(1,1,n)=transp(1,1,n)+transp(k,1,n)
     end do			! vertical loop
     print 107,ArchvTimeUnit,its2time(nstart+nstep),			&
       trim(thruflsecs(n)),'total thruflow (Ktons/s)',transp(1,1,n)
107  format (a,f11.1,2(2x,a),f11.1)
    end if
   end do			! loop over transects

   print *,'-------atm.flxsum--------atm.flxsum--------atm.flxsum--------atm.flxsum-------'
!  print *,'...exiting atm_flxsum'
   call flush(6)
   return
   end subroutine atm_flxsum


  subroutine trnsec1(its,nf,of,vof,adbash1,adbash2,adbash3,		&
                     flxavg,massfx)

! --- build up mass flux time integral. to be called each time step.

  use module_control  ,only: nvl,npp,nip
  use fimnamelist     ,only: ptop,PrintIpnDiag
  use module_constants,only: nprox,prox,thetac,perm,inv_perm,p1000
  use module_variables,only: tr3d,dp3d

  implicit none
  integer,intent(IN)  :: its		! model time step
  integer,intent(IN)  :: nf,of,vof	! time slots: new,old,very old
  real   ,intent(IN)  :: adbash1,adbash2,adbash3

!sms$distribute (dh,3) begin
  real,intent(IN)     :: massfx(nvl,npp,nip,3)
  real,intent(INOUT)  :: flxavg(nvl,npp,nip)
!sms$distribute end
  integer    :: ico		! Index for icos point number
  integer    :: edg,k,n
  logical    :: vrbos
  real*8     :: thold(nvl),thnew(nvl),prold(nvl+1),prnew(nvl+1),	&
                velo(nvl),flux(nvl),znew(nvl+1),w1,w2
  real       :: vrbl4
! integer,parameter :: trform=0		! do nothing
  integer,parameter :: trform=1		! transform fluxes to pure theta coord.
! integer,parameter :: trform=2		! transform fluxes to pure sigma coord.

!SMS$PARALLEL (dh,ico) BEGIN
!$OMP PARALLEL DO PRIVATE(flux,prold,thold,velo,w1,w2,prnew,thnew,vrbos)
  do ico=1,nip			! horizontal loop
   vrbos=ico.eq.PrintIpnDiag
   do edg=1,nprox(ico)
    do k=1,nvl
     flux(k)=adbash1*massfx(k,edg,ico, nf)				&
            +adbash2*massfx(k,edg,ico, of)				&
            +adbash3*massfx(k,edg,ico,vof)
    end do

! --- transform fluxes to pure isentropic or pure sigma space before processing
    if (trform.gt.0) then
     prold(nvl+1)=ptop
     do k=nvl,1,-1
      prold(k)=prold(k+1)+.5*(dp3d(k,ico)+dp3d(k,prox(edg,ico)))
      thold(k)=.5*(tr3d(k,ico,1)+tr3d(k,prox(edg,ico),1))
      velo(k)=flux(k)/(.5*(dp3d(k,ico)+dp3d(k,prox(edg,ico))))	! flux => velo
     end do
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- theta option: find pressure on isentropic layer interfaces
     if (trform.eq.1) then
      call rstep1d(thold,prold,nvl,thnew,prnew,thetac,nvl,vrbos,perm(ico))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- sigma option: find pressure on terrain-following surfaces
     else if (trform.eq.2) then
      do k=1,nvl+1
       w1=float(nvl+1-k)/float(nvl)
       w2=w1*w1
       w1=w1*sqrt(w1)				! reduce dp as p => 0
       prnew(k)=p1000*w1+prold(nvl+1)*(1.-w1) + (prold(1)-p1000)*w2
       if (k.gt.1 .and. prnew(k).ge.prnew(k-1)) then
        print *,'(trnsec1) dp<0 at ipn =',perm(ico)
        stop '(error trnsec1)'
       end if
      end do

      if (vrbos) then
       print '(a,i7/-2p,(10f8.2))',		&
         '(trnsc1) old pressure levels at ico=',ico,(prold(k),k=1,nvl+1)
       print '(a,i7/-2p,(10f8.2))',		&
         '(trnsc1) new pressure levels at ico=',ico,(prnew(k),k=1,nvl+1)
       print '(a,i7/-2p,(10f8.2))',		&
         '(trnsc1) thickness levels at ico=',ico,(prnew(k)-prnew(k+1),k=1,nvl)
       call flush(6)
       do k=1,nvl+1
        vrbl4=75.*prnew(k)/p1000
        call linout(vrbl4,'X',-1)
        call linout(vrbl4,'X', 1,k)
       end do
       call flush(6)
      end if
     end if			! trform = 1 or 2
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! --- interpolate fluxes to new layers. note: input=velocity, output=flux
     call remap1d(nvl,prold,velo,prnew,flux,vrbos,perm(ico))
    end if			! trform > 0

    do k=1,nvl
     flxavg(k,edg,ico)=flxavg(k,edg,ico)+flux(k)
    end do
   end do			! loop over edges
  end do			! horizontal loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END

  return
  end subroutine trnsec1


  subroutine rstep1d(thold,pold,kold,thnew,pnew,targt,knew,vrbos,ipn)
  
! --- convert a stairstep (i.e., piecewise constant) theta profile into a
! --- stairstep profile constrained to have prescribed theta ('targt') steps.

! --- input  variables: thold,pold,targt,kold,knew,vrbos
! --- output variables: thnew,pknew

  implicit none
  integer,intent(IN)  :: kold,knew,ipn
  real   ,intent(IN)  :: targt(knew)
  real*8 ,intent(IN)  :: thold(kold),pold(kold+1)
  logical,intent(IN)  :: vrbos		! if true, print diagnostics at ipn
  real*8 ,intent(OUT) :: thnew(knew),pnew(knew+1)

  integer k,ko
  real*8 oldth(kold)
  real*8 cloutt,colint,colscl,pinteg,tha,thb
! real,parameter :: acurcy=1.e-6	! for 32-bit arithmetic
  real,parameter :: acurcy=1.e-12	! for 64-bit arithmetic
  
  if (vrbos) write (6,101) ipn,						&
     'rstep1d -- input profile:   theta    thknss     press',		&
     pold(1),(k,thold(k),pold(k)-pold(k+1),pold(k+1),k=1,kold)
101 format (i8,4x,a/56x,f9.1/(34x,i3,f9.3,f9.1,f10.1))
  
! --- remove theta inversions from input profile
  oldth(kold)=thold(kold)
  do k=kold-1,1,-1
    oldth(k)=min(oldth(k+1),thold(k))
  end do
  
  thnew(:)=targt(:)
  thnew(   1)=min(oldth(1),oldth(kold),targt(   1))
  thnew(knew)=max(oldth(1),oldth(kold),targt(knew))
  pnew(     1)=pold(     1)
  pnew(knew+1)=pold(kold+1)
  
! --- column integrals (colin/clout) are computed for diagnostic purposes only
  
  cloutt=0.
  colint=0.
  colscl=0.
  do k=1,kold
    colint=colint+oldth(k)*(pold(k+1)-pold(k))
    colscl=colscl+abs(oldth(k)*(pold(k+1)-pold(k)))
  end do
  
! --- find interface pnew(k+1) separating layers k and k+1 by requiring
! --- that integral over p*d(theta) from thnew(k) to thnew(k+1) be preserved.
  
  ko=1
  do k=1,knew-1
    pinteg=0.
    thb=thnew(k)
 5  tha=thb
    thb=min(thnew(k+1),max(thnew(k),oldth(ko)))
    pinteg=pinteg+pold(ko)*(thb-tha)
    if (oldth(ko) < thnew(k+1)) then
      if (ko < kold) then
        ko=ko+1
        go to 5
      end if
      tha=thb
      thb=thnew(k+1)
      pinteg=pinteg+pold(kold+1)*(thb-tha)
    end if
    pnew(k+1)=pnew(k)
    if (thnew(k+1) > thnew(k)) pnew(k+1)=pinteg/(thnew(k+1)-thnew(k))
    cloutt=cloutt+thnew(k)*(pnew(k+1)-pnew(k))
  enddo
  
  cloutt=cloutt+thnew(knew)*(pnew(knew+1)-pnew(knew))
  if (abs(cloutt-colint) > max(1.e-33,acurcy*colscl)) then
    write (6,100) ipn,'rstep1d - column intgl.error',			&
    colint,cloutt,(cloutt-colint)/max(1.e-33,abs(colint))
100 format (i8,3x,a,2es16.8,es9.1)
  end if
  
  if (vrbos) write (6,101) ipn,					    	&
     'rstep1d -- outpt profile:   theta    thknss     press',		&
     pnew(1),(k,thnew(k),pnew(k)-pnew(k+1),pnew(k+1),k=1,knew)
  
  return
  end subroutine rstep1d


  subroutine remap1d(kk,xold,yold,xnew,ynew,vrbos,ipn)

!-- 1-dim PCM/PLM/PPM transport routine, extracted from HYCOM's fct3d.f.
!-- i m p o r t a n t:  this version returns ynew*dp, not ynew

!-- there are 3 advection/remapping choices:
!-- remap_optn=1			! PCM (piecewise constant)
!-- remap_optn=2			! PLM (piecewise linear)
!-- remap_optn=3			! PPM (piecewise parabolic)

    implicit none
    integer  ,intent(IN)  :: kk		! vert.dimension
    logical  ,intent(IN)  :: vrbos	! if true, print diagnostics at ipn
    integer  ,intent(IN)  :: ipn	! grid cell number

!-- xold/new	- old/new cell boundaries
!-- fldold/new	- mixing ratio of dep.variable before and after transport

    real*8 ,intent(IN)  :: xold(kk+1),xnew(kk+1),yold(kk)
    real*8 ,intent(OUT) :: ynew(kk)

    real*8 zold(kk+1),znew(kk+1),delz(kk+1),fco(kk),fcn(kk),		&
         vertfx(kk+1),vertdv(kk),fldold(kk),fldnew(kk)
    real*8 a(kk),b(kk),c(kk),dx,fcdx,yl,yr,colscl
    real*8 amount,bfore,after,slab,dslab
    integer k,lyr
    real*8, parameter :: athird=1./3.
    real*8, parameter :: small=1.e-22
    real*8, parameter :: acurcy=1.e-11
    integer,parameter :: remap_optn=1

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

    if (vrbos) write (*,100) ipn,'(remap1d) input profile:',		&
      'old_indep   old_binsz  lowr.displ  uppr.displ   old_depnd',	&
      (k,zold(k),zold(k+1)-zold(k),delz(k),delz(k+1),fldold(k),k=1,kk),	&
      kk+1,zold(kk+1)
100 format (i8,3x,a/9x,a/(i6,f12.1,es12.3,2f12.1,es12.3))

!-- deduce old and new cell width from -zold,znew-
    do 15 k=1,kk
    fco(k)=max(0.,zold(k+1)-zold(k))
 15 fcn(k)=max(0.,znew(k+1)-znew(k))

    bfore=0.
    colscl=0.
    do k=1,kk
      bfore=bfore+fldold(k)*fco(k)
      colscl=colscl+abs(fldold(k)*fco(k))
    end do
    fldnew(:)=fldold(:)

!-- start by filling zero-width cells with data from neighboring cells

!    do 17 k=kk-1,1,-1
!17  fldnew(k)=(fldnew(k)*fco(k)+fldnew(k+1)*small)			&
!             /(          fco(k)+            small)
!    do 18 k=2,kk
!18  fldnew(k)=(fldnew(k)*fco(k)+fldnew(k-1)*small)			&
!             /(          fco(k)+            small)

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
      print *,'subr.remap1d -- illegal remap_optn=',remap_optn
      stop '(remap1d)'
    end if
16  continue

!-- get flux by summing -fldnew- over upstream slab of thickness -delz-

    do 22 k=2,kk
    slab=0.
    amount=0.
    vertfx(k)=0.
    if (delz(k) > 0.) then			! interface moves in +k dir.
      lyr=k-1
24    lyr=lyr+1
      if (slab >= delz(k)) goto 23
      if (fco(lyr) > 0.) then
        dslab=min(slab+fco(lyr), delz(k))	&
             -min(slab         , delz(k))
        dx=dslab/fco(lyr)
        fcdx=a(lyr)				&
            +b(lyr)*.5*(dx-1.)			& !  not needed in pcm
            +c(lyr)*(.25-dx*(.5-dx*athird))	  !  not needed in pcm,plm
        amount=amount+fcdx*dslab
        slab=slab+dslab
      end if
      if (lyr < kk) go to 24
    else if (delz(k) < 0.) then		! interface moves in -k dir.
      lyr=k
25    lyr=lyr-1
      if (slab >= -delz(k)) goto 23
      if (fco(lyr) > 0.) then
        dslab=min(slab+fco(lyr),-delz(k))	&
             -min(slab         ,-delz(k))
        dx=dslab/fco(lyr)
        fcdx=a(lyr)				&
            +b(lyr)*.5*(1.-dx)			& !  not needed in pcm
            +c(lyr)*(.25-dx*(.5-dx*athird))	  !  not needed in pcm,plm
        amount=amount+fcdx*dslab
        slab=slab+dslab
      end if
      if (lyr > 1) go to 25
    end if
23  if (slab.ne.0.) vertfx(k)=-delz(k)*amount/slab
22  continue

    vertfx(   1)=0.			!  don't allow flux through lower bdry
    vertfx(kk+1)=0.			!  don't allow flux through upper bdry
    do 26 k=1,kk
26  vertdv(k)=vertfx(k+1)-vertfx(k)

    if (vrbos) write (*,'(a/(i3,4es12.3))')				&
     'remap1d:  flux  flx.div/thk    old_thk     new_thk',		&
      (k,vertfx(k),vertdv(k)/max(small,fcn(k)),fco(k),fcn(k),k=1,kk),	&
       kk+1,vertfx(kk+1)

    after=0.
    do 4 k=1,kk
    amount=fldnew(k)*fco(k)-vertdv(k)
!   fldnew(k)=(amount+fldnew(k)*small)/(fcn(k)+small)
    fldnew(k)= amount			! return dependent variable x thickness

!   if (miglim.eq.1 .and. k.gt.1 .and. k.lt.kk) then
!    if (fldnew(k).lt.min(fldold(k-1),fldold(k),fldold(k+1))		&
!      -max(1.,abs(fldnew(k)))*acurcy .or.				&
!        fldnew(k).gt.max(fldold(k-1),fldold(k),fldold(k+1))		&
!      +max(1.,abs(fldnew(k)))*acurcy)					&
!    print '(i8,i4,a,f11.2,a,3es11.3)',ipn,k,' (remap1d)',		&
!      fldnew(k),' outside range',fldnew(k)-fldold(k-1),		&
!      fldnew(k)-fldold(k),fldnew(k)-fldold(k+1)
!   end if

    after=after+amount
 4  continue

    if (abs(bfore-after) > max(1.e-33,acurcy*colscl)) then
      write (*,104) ipn,'remap1d - column intgl.error',bfore,after,	&
      (after-bfore)/max(1.e-33,abs(bfore))
    end if
104 format (i8,3x,a,2es16.8,es9.1)

    if (vrbos) then
      write (*,100) ipn,'(remap1d) output profile:',			&
      'new_indep   new_binsz  lowr.displ  uppr.displ   new_depnd',	&
      (k,znew(k),znew(k+1)-znew(k),delz(k),delz(k+1),fldnew(k),		&
      k=1,kk),kk+1,znew(kk+1)
    end if

    ynew(:)=fldnew(:)
    return
  end subroutine remap1d


   subroutine linout(value,char,index,labl)
   implicit none
   real       ,intent(IN) :: value
   integer    ,intent(IN) :: index
   integer    ,intent(IN),optional :: labl
   character*1,intent(IN) :: char
   integer, parameter :: length=77
   character(len=length) line
   integer l,n
   save line

! --- replace n-th element of array 'line' by character 'char', where
! --- n = 'value' modulo 'length'
! --- index < 0   -- initialize 'line' by blanks before adding 'char'
! --- index > 0   -- output 'line' after adding 'char'
! --- 'labl', if given, is a 1-3 digit number printed at the beginning of 'line'

   if (index.lt.0) line=' '
   n=int(mod(value+float(length-1),float(length)))+1
   line(n:n)=char

   if (index.gt.0) then
     if (present(labl)) then
       write (*,'(i3,a)') labl,line
     else
       write (*,'("=+",a)') line
    end if
   end if
   return
   end subroutine linout
end module mdul_transects
