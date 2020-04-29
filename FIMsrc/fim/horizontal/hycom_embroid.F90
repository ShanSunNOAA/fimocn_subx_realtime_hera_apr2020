module hycom_embroid
use stencilprint
contains
!*********************************************************************
!  embroid
!
!  extrapolates land surface properties to 'new' grid cells where
!  these properties are not yet defined
!
!   R. Bleck                    Nov  2010
!*********************************************************************
   subroutine embroid
   use module_control,  only: nip
   use fimnamelist,     only: PrintIpnDiag,stencl_step,itest
   use module_constants,only: nprox,prox,deg_lat,deg_lon,perm
   use module_variables,only: ph3d
   use hycom_constants, only: wet
   use gfs_physics_internal_state_mod, only:gis_phy
!
   implicit none
   integer,parameter :: lsoil=4		! number of soil layers
   integer k,atm_step,atm_Diag
!SMS$DISTRIBUTE (dh,1) BEGIN
   real    oldmsk(nip)
!
! --- the following variables are patterned after parameters appearing
! --- in calls to gbphys and gloopr in do_physics_one_step.F90.
! --- caution: this list of surface parameters may not be complete
!
   real	alnsf1(nip),		&	! albedo for near-IR scattered
        alnwf1(nip),		&	! albedo for near-IR beam
        alvsf1(nip),		&	! albedo for visible scattered
        alvwf1(nip),		&	!  albedo for visible beam
        canopy1(nip),		&	! canopy water (m)
        facsf1(nip),		&	! frac'l cov. w.strong cosz dependency
        facwf1(nip),		&	! frac'l cov. w.weak cosz dependency
        shdmax1(nip),		&	! max.frac'l cov. of green vegetation
        shdmin1(nip),		&	! min.frac'l cov. of green vegetation
        sheleg1(nip),		&	! snow depth equivalent (mm)
        slmsk1(nip),		&	! sea-land mask (0-sea,1-land,2-sea ice)
        slope1(nip),		&	! slope class (from namelist_soilveg)
!       sncovr1(nip),		&	! fractional snow cover
        snoalb1(nip),		&	! upper bound for albedo of deep snow
        snwdph1(nip),		&	! snow depth (mm)
        stype1(nip),		&	! soil type (integer 1-9)
        tg31(nip),		&	! deep soil temperature (K)
        vfrac1(nip),		&	! vegetation fraction
        vtype1(nip),		&	! vegetation type (integer 1-13)
        zorl1(nip) 		 	! roughness in cm
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE (dh,2) BEGIN
   integer pprox(6,nip)
!
! --- the following variables are patterned after parameters appearing
! --- in calls to gbphys and gloopr in do_physics_one_step.F90.
! --- caution: this list of surface parameters may not be complete
!
   real slc1(lsoil,nip),	&	! liquid soil moisture
        smc1(lsoil,nip),	&	! soil volumetric water content
        stc1(lsoil,nip)			! soil temperature (K)
!SMS$DISTRIBUTE END
   real ovrn,switch
   integer i,n,n1,n2,newpt,newtot,npass

   print *,'entering embroid ...'
   pprox(:,:)=0		! need this to pass lahey
   newtot=0
   npass=0

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO
   do i=1,nip
    alnsf1 (i)=gis_phy%sfc_fld%alnsf (i,1)
    alnwf1 (i)=gis_phy%sfc_fld%alnwf (i,1)
    alvsf1 (i)=gis_phy%sfc_fld%alvwf (i,1)
    alvwf1 (i)=gis_phy%sfc_fld%alvsf (i,1)
    canopy1(i)=gis_phy%sfc_fld%canopy(i,1)
    facsf1 (i)=gis_phy%sfc_fld%facsf (i,1)
    facwf1 (i)=gis_phy%sfc_fld%facwf (i,1)
    shdmax1(i)=gis_phy%sfc_fld%shdmax(i,1)
    shdmin1(i)=gis_phy%sfc_fld%shdmin(i,1)
    sheleg1(i)=gis_phy%sfc_fld%sheleg(i,1)
    slmsk1 (i)=gis_phy%sfc_fld%slmsk (i,1)
    slope1 (i)=gis_phy%sfc_fld%slope (i,1)
!   sncovr1(i)=gis_phy%sfc_fld%sncovr(i,1)
    snoalb1(i)=gis_phy%sfc_fld%snoalb(i,1)
    snwdph1(i)=gis_phy%sfc_fld%snwdph(i,1)
    stype1 (i)=gis_phy%sfc_fld%stype (i,1)
    tg31   (i)=gis_phy%sfc_fld%tg3   (i,1)
    vfrac1 (i)=gis_phy%sfc_fld%vfrac (i,1)
    vtype1 (i)=gis_phy%sfc_fld%vtype (i,1)
    zorl1  (i)=gis_phy%sfc_fld%zorl  (i,1)
    slc1 (:,i)=gis_phy%sfc_fld%slc (:,i,1)
    smc1 (:,i)=gis_phy%sfc_fld%smc (:,i,1)
    stc1 (:,i)=gis_phy%sfc_fld%stc (:,i,1)

    if (i.eq.itest) then
     print 103,perm(i),i,'  wet/slmsk/stype on entry to embroid:',	&
       wet(i),slmsk1(i),stype1(i)
     call flush(6)
103  format (i7,' (lcl',i8,')',a,i2,2f5.1)
    end if
   end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END

   atm_step=stencl_step
   stencl_step=1
   atm_Diag=PrintIpnDiag
   PrintIpnDiag=itest

   call stencl(slmsk1,1,100.,'(ocn_embroid) slmsk1 x 100')
   call stencl(stype1,1,100.,'(ocn_embroid) soil type x 100')
!  call stencl(sncovr1,1,1.,'(ocn_embroid) snow cover')

   PrintIpnDiag=atm_Diag
   stencl_step=atm_step
   
1  do n=6,0,-1
    npass=npass+1
    print '(a,i5,2(a,i10))','pass',npass,				&
     ' -- looking for uninitialized points with',n,' neighbors'

! --- look for uninitialized cells with -n- initialized neighbors
    newpt=0
!SMS$PARALLEL (dh,i) BEGIN
    do i=1,nip
     oldmsk(i)=slmsk1(i)
    end do

!SMS$EXCHANGE(alnsf1,alnwf1,alvsf1,alvwf1,canopy1,facsf1,facwf1,shdmax1,shdmin1,sheleg1,slope1,snoalb1,snwdph1,stype1,tg31,vfrac1,vtype1,zorl1,slc1,smc1,stc1,oldmsk)

!JR Have not threaded this loop because 1) handling newpt and newtot is tricky and 2) it appears
!JR that this routine is only invoked at initialization?

    do i=1,nip

! --- test for missing soil type at land points
     if (slmsk1(i).eq.1. .and.  stype1(i).eq.0.) then
      print 102,npass,'error -- zero soil type at land point',perm(i),i
      call flush(6)
102   format ('pass',i5,1x,a,i8,' (lcl',i8,')',a,l2,2f5.1)
      stype1(i)=1.
     end if
      
     if (wet(i)>0 .and. oldmsk(i)==1.) then

! --- ocean model wants this land cell to be ocean
      switch=0.
      if (i.eq.itest) then
       print 101,npass,'turn land into OCEAN at',perm(i),i,		&
        '  lat/lon',deg_lat(i),deg_lon(i)
101   format ('pass',i5,1x,a,i8,' (lcl',i8,')',a,2f8.2)
       call flush(6)
      end if
      slmsk1(i)=0.		! convert to ocean
      go to 2
     else if (wet(i)==0 .and. oldmsk(i).ne.1. .and. oldmsk(i).ne.2.) then

! --- ocean model wants this ocean cell to be land
      if (ph3d(1,i).gt.0.) then
       print '(a,i8," (lcl",i8,")",f11.2)',				&
         'error: nonzero terrain height at ocean point',		&
         perm(i),i,ph3d(1,i)
       call flush(6)
      end if
      switch=1.
      if (i.eq.itest) then
       print 101,npass,'uninitialized LAND  at',perm(i),i,		&
        '  lat/lon',deg_lat(i),deg_lon(i)
       call flush(6)
      end if

     else
      go to 2
     end if

! --- check availability of neighbors
     n2=0
     do n1=1,nprox(i)
      if (abs(oldmsk(prox(n1,i))-switch) .lt. .001) then
       n2=n2+1
       pprox(n2,i)=prox(n1,i)
      end if
     end do

! --- are there -n- or more neighbors?
     if (n2.ge.n) then
      if (n2.gt.0) then
       print 104,npass,'available neighbors at',perm(i),i,		&
        (n1,perm(pprox(n1,i)),perm(pprox(n1,i)),n1=1,n2)
104    format ('pass',i5,1x,a,i8,' (lcl',i8,')'/(i32,i8,' (lcl',i8,')'))
       call flush(6)

! --- deduce parameter values in uninitialized cell from neighboring cells
       alnsf1 (i)=0.
       alnwf1 (i)=0.
       alvsf1 (i)=0.
       alvwf1 (i)=0.
       canopy1(i)=0.
       facsf1 (i)=0.
       facwf1 (i)=0.
       shdmax1(i)=0.
       shdmin1(i)=0.
       sheleg1(i)=0.
       slope1 (i)=0.
!      sncovr1(i)=0.
       snoalb1(i)=0.
       snwdph1(i)=0.
       stype1 (i)=0.
       tg31   (i)=0.
       vfrac1 (i)=0.
       vtype1 (i)=0.
       zorl1  (i)=0.
       slc1(:,i) =0.
       smc1(:,i) =0.
       stc1(:,i) =0.
99     format (i7,' (lcl',i8,')   ',a,6f8.2)
       print 99,perm(i),i,'alnsf ',(alnsf1 (pprox(n1,i)),n1=1,n2)
       print 99,perm(i),i,'alnwf ',(alnwf1 (pprox(n1,i)),n1=1,n2)
       print 99,perm(i),i,'alvsf ',(alvsf1 (pprox(n1,i)),n1=1,n2)
       print 99,perm(i),i,'alvwf ',(alvwf1 (pprox(n1,i)),n1=1,n2)
       print 99,perm(i),i,'canopy',(canopy1(pprox(n1,i)),n1=1,n2)
       print 99,perm(i),i,'facsf ',(facsf1 (pprox(n1,i)),n1=1,n2)
       print 99,perm(i),i,'facwf ',(facwf1 (pprox(n1,i)),n1=1,n2)
       print 99,perm(i),i,'shdmax',(shdmax1(pprox(n1,i)),n1=1,n2)
       print 99,perm(i),i,'shdmin',(shdmin1(pprox(n1,i)),n1=1,n2)
       print 99,perm(i),i,'sheleg',(sheleg1(pprox(n1,i)),n1=1,n2)
!      print 99,perm(i),i,'slope ',(slope1 (pprox(n1,i)),n1=1,n2)
!      print 99,perm(i),i,'sncovr',(sncovr1(pprox(n1,i)),n1=1,n2)
!      print 99,perm(i),i,'snoalb',(snoalb1(pprox(n1,i)),n1=1,n2)
!      print 99,perm(i),i,'snwdph',(snwdph1(pprox(n1,i)),n1=1,n2)
!      print 99,perm(i),i,'stype ',(stype1 (pprox(n1,i)),n1=1,n2)
!      print 99,perm(i),i,'tg3   ',(tg31   (pprox(n1,i)),n1=1,n2)
!      print 99,perm(i),i,'vfrac ',(vfrac1 (pprox(n1,i)),n1=1,n2)
!      print 99,perm(i),i,'vtype ',(vtype1 (pprox(n1,i)),n1=1,n2)
!      print 99,perm(i),i,'zorl  ',(zorl1  (pprox(n1,i)),n1=1,n2)
!      print 99,perm(i),i,'slc   ',(slc1 (1,pprox(n1,i)),n1=1,n2)
!      print 99,perm(i),i,'smc   ',(smc1 (1,pprox(n1,i)),n1=1,n2)
!      print 99,perm(i),i,'stc   ',(stc1 (1,pprox(n1,i)),n1=1,n2)
       call flush(6)
       do n1=1,n2
        alnsf1 (i)=alnsf1 (i)+alnsf1 (pprox(n1,i))
        alnwf1 (i)=alnwf1 (i)+alnwf1 (pprox(n1,i))
        alvsf1 (i)=alvsf1 (i)+alvwf1 (pprox(n1,i))
        alvwf1 (i)=alvwf1 (i)+alvsf1 (pprox(n1,i))
        canopy1(i)=canopy1(i)+canopy1(pprox(n1,i))
        facsf1 (i)=facsf1 (i)+facsf1 (pprox(n1,i))
        facwf1 (i)=facwf1 (i)+facwf1 (pprox(n1,i))
        shdmax1(i)=shdmax1(i)+shdmax1(pprox(n1,i))
        shdmin1(i)=shdmin1(i)+shdmin1(pprox(n1,i))
        sheleg1(i)=sheleg1(i)+sheleg1(pprox(n1,i))
!       slope1 (i)=slope1 (i)+slope1 (pprox(n1,i))
        slope1 (i)=           slope1 (pprox(n1,i))
!       sncovr1(i)=sncovr1(i)+sncovr1(pprox(n1,i))
        snoalb1(i)=snoalb1(i)+snoalb1(pprox(n1,i))
        snwdph1(i)=snwdph1(i)+snwdph1(pprox(n1,i))
!       stype1 (i)=stype1 (i)+stype1 (pprox(n1,i))
        stype1 (i)=           stype1 (pprox(n1,i))
        tg31   (i)=tg31   (i)+tg31   (pprox(n1,i))
        vfrac1 (i)=vfrac1 (i)+vfrac1 (pprox(n1,i))
!       vtype1 (i)=vtype1 (i)+vtype1 (pprox(n1,i))
        vtype1 (i)=           vtype1 (pprox(n1,i))
        zorl1  (i)=zorl1  (i)+zorl1  (pprox(n1,i))
        slc1(:,i) =slc1(:,i) +slc1 (:,pprox(n1,i))
        smc1(:,i) =smc1(:,i) +smc1 (:,pprox(n1,i))
        stc1(:,i) =stc1(:,i) +stc1 (:,pprox(n1,i))
       end do
       ovrn=1./n2
       alnsf1 (i)=alnsf1 (i)*ovrn
       alnwf1 (i)=alnwf1 (i)*ovrn
       alvsf1 (i)=alvwf1 (i)*ovrn
       alvwf1 (i)=alvsf1 (i)*ovrn
       canopy1(i)=canopy1(i)*ovrn
       facsf1 (i)=facsf1 (i)*ovrn
       facwf1 (i)=facwf1 (i)*ovrn
       shdmax1(i)=shdmax1(i)*ovrn
       shdmin1(i)=shdmin1(i)*ovrn
       sheleg1(i)=sheleg1(i)*ovrn
!      slope1 (i)=slope1 (i)*ovrn
!      sncovr1(i)=sncovr1(i)*ovrn
       snoalb1(i)=snoalb1(i)*ovrn
       snwdph1(i)=snwdph1(i)*ovrn
!      stype1 (i)=stype1 (i)*ovrn
       tg31   (i)=tg31   (i)*ovrn
       vfrac1 (i)=vfrac1 (i)*ovrn
!      vtype1 (i)=vtype1 (i)*ovrn
       zorl1  (i)=zorl1  (i)*ovrn
       slc1 (:,i)=slc1 (:,i)*ovrn
       smc1 (:,i)=smc1 (:,i)*ovrn
       stc1 (:,i)=stc1 (:,i)*ovrn
       if (switch.eq.1.) then
        print 100,npass,'LAND  parms inferred at',perm(i),i,		&
         '  lat/lon',deg_lat(i),deg_lon(i)
       end if
       if (switch.ne.1.) then
        print 100,npass,'OCEAN parms interpolated to',perm(i),i,	&
         '  lat/lon',deg_lat(i),deg_lon(i)
 100    format ('pass',i5,1x,a,i8,' (lcl',i8,')',a,2f8.2)
       end if
!      print 99,perm(i),i,'alnsf ',alnsf1 (i)
!      print 99,perm(i),i,'alnwf ',alnwf1 (i)
!      print 99,perm(i),i,'alvsf ',alvsf1 (i)
!      print 99,perm(i),i,'alvwf ',alvwf1 (i)
!      print 99,perm(i),i,'canopy',canopy1(i)
!      print 99,perm(i),i,'facsf ',facsf1 (i)
!      print 99,perm(i),i,'facwf ',facwf1 (i)
!      print 99,perm(i),i,'shdmax',shdmax1(i)
!      print 99,perm(i),i,'shdmin',shdmin1(i)
!      print 99,perm(i),i,'sheleg',sheleg1(i)
!      print 99,perm(i),i,'slope ',slope1 (i)
!      print 99,perm(i),i,'sncovr',sncovr1(i)
!      print 99,perm(i),i,'snoalb',snoalb1(i)
!      print 99,perm(i),i,'snwdph',snwdph1(i)
!      print 99,perm(i),i,'stype ',stype1 (i)
!      print 99,perm(i),i,'tg3   ',tg31   (i)
!      print 99,perm(i),i,'vfrac ',vfrac1 (i)
!      print 99,perm(i),i,'vtype ',vtype1 (i)
!      print 99,perm(i),i,'zorl  ',zorl1  (i)
!      print 99,perm(i),i,'slc   ',slc1 (1,i)
!      print 99,perm(i),i,'smc   ',smc1 (1,i)
!      print 99,perm(i),i,'stc   ',stc1 (1,i)
       call flush(6)
       newpt=newpt+1
       slmsk1(i)=switch

      else			! n2=0

! --- no neighbors with desired properties
! --- most likely an isolated single-point island => submerge
!!!    wet(i)=1

       print 100,npass,'no neighbors to infer properties from at ',	&
         perm(i),i,'  lat/lon',deg_lat(i),deg_lon(i)
!!!    print 102,npass,'change',perm(i),i,'  into ocean point'
       call flush(6)

      end if			! n2 > or = 0
     end if			! desired number of neighbors exist/don't exist
2    continue
    end do			! horiz.loop
!SMS$PARALLEL END
!SMS$REDUCE(newpt,SUM)

    if (newpt.gt.0) then
     print *,newpt,'  new points initialized. update halo and begin new search'
     newtot=newtot+newpt
     call flush(6)
     go to 1
    end if			! newpt > 0
    print '(a,i3,a)','no points found having',n,' neighbors'
   end do			! loop over number of neighbors

!SMS$PARALLEL (dh,i) BEGIN
     do i=1,nip
      gis_phy%sfc_fld%alnsf (i,1)=alnsf1 (i)
      gis_phy%sfc_fld%alnwf (i,1)=alnwf1 (i)
      gis_phy%sfc_fld%alvwf (i,1)=alvsf1 (i)
      gis_phy%sfc_fld%alvsf (i,1)=alvwf1 (i)
      gis_phy%sfc_fld%canopy(i,1)=canopy1(i)
      gis_phy%sfc_fld%facsf (i,1)=facsf1 (i)
      gis_phy%sfc_fld%facwf (i,1)=facwf1 (i)
      gis_phy%sfc_fld%shdmax(i,1)=shdmax1(i)
      gis_phy%sfc_fld%shdmin(i,1)=shdmin1(i)
      gis_phy%sfc_fld%sheleg(i,1)=sheleg1(i)
      gis_phy%sfc_fld%slope (i,1)=slope1 (i)
      gis_phy%sfc_fld%slmsk (i,1)=slmsk1 (i)
!     gis_phy%sfc_fld%sncovr(i,1)=sncovr1(i)
      gis_phy%sfc_fld%snoalb(i,1)=snoalb1(i)
      gis_phy%sfc_fld%snwdph(i,1)=snwdph1(i)
      gis_phy%sfc_fld%stype (i,1)=stype1 (i)
      gis_phy%sfc_fld%tg3   (i,1)=tg31   (i)
      gis_phy%sfc_fld%vfrac (i,1)=vfrac1 (i)
      gis_phy%sfc_fld%vtype (i,1)=vtype1 (i)
      gis_phy%sfc_fld%zorl  (i,1)=zorl1  (i)
      gis_phy%sfc_fld%slc (:,i,1)=slc1 (:,i)
      gis_phy%sfc_fld%smc (:,i,1)=smc1 (:,i)
      gis_phy%sfc_fld%stc (:,i,1)=stc1 (:,i)

      if (i.eq.itest) then
       print 103,perm(i),i,'  wet/slmsk/stype on exit from embroid:',	&
         wet(i),slmsk1(i),stype1(i)
       call flush(6)
      end if
     end do
!SMS$EXCHANGE (wet)
!SMS$PARALLEL END

   print *,'... exiting embroid, total number of converted points:',newtot
   return
   end subroutine embroid
end module hycom_embroid
