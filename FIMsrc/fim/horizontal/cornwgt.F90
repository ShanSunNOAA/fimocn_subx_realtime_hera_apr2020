!***********************************************************************
!     cornwgt
!     construct weights for interpolation of cell values to cell corners
!       R. Bleck     Aug 2012
!***********************************************************************

module module_cornwgt

  implicit none

  private locplt,sinangl

CONTAINS

  subroutine rdcnfig ! read vertex coords from configuration file

    use module_constants,only : trlpt,ldgpt
    use module_control,  only : nip,npp
    use units,           only : getunit, returnunit

    character(len=16)           :: header
    character(len=24),parameter :: fn = 'icos_grid_info_level.dat' ! input file
    integer                     :: ierr                            ! status var
    integer                     :: edg                             ! edge index
    integer                     :: ilr
    integer                     :: ipn
    integer                     :: ixy
    integer                     :: unitno                          ! io unit
    real,allocatable            :: corners(:,:,:,:)

!SMS$SERIAL (<trlpt,ldgpt,OUT>:default=ignore) BEGIN

    allocate(corners(npp,2,2,nip),stat=ierr)
    if (ierr.ne.0) stop 'Could not allocate corners'

    unitno = getunit ()
    if (unitno < 0) then
      print*,'rdcnfig: getunit failed. Stopping'
      stop
    end if

    print *,'(cornwgt) get vertex info from file  icos_grid_info_level.dat'

    open (unitno, file=fn, form='unformatted', action='read', status='old', iostat=ierr)
    if (ierr /= 0) then
      print *,'(cornwgt) iostat =', ierr, ' attempting to open file=', fn
      stop
    end if

    read (unitno, iostat=ierr) header
    if (ierr /= 0) then
      print *,'(cornwgt) iostat =',ierr,' after reading header 1. exiting.'
      stop
    end if
    print *,header

    read (unitno, iostat=ierr) header
    if (ierr /= 0) then
      print *,'(cornwgt) iostat =',ierr,' after reading header 2. exiting.'
      stop
    end if
    print *,header

    read (unitno,iostat=ierr) ! read past lat,lon
    if (ierr /= 0) then
      print *,'(cornwgt) iostat =',ierr,' after reading lat,lon. exiting.'
      stop
    end if

    do edg=1,npp
      read (unitno, iostat=ierr) ! read past prox
      if (ierr /= 0) then
        print *,'(cornwgt) iostat =',ierr,' while reading prox. exiting.'
        stop
      end if
    end do

    read (unitno,iostat=ierr) ! read past nprox
    if (ierr /= 0) then
      print *,'(cornwgt) iostat =',ierr,' after reading nprox. exiting.'
      stop
    end if

    do ixy = 1,2
      do ilr = 1,2
        do edg = 1,npp
          READ (unitno,iostat=ierr) (corners(edg,ilr,ixy,ipn),ipn=1,nip)
          if (ierr /= 0) then
            print *,'(cornwgt) iostat =',ierr,' while reading corner info. exiting.'
            stop
          end if
        enddo
      enddo
    enddo

    close (unitno)

    call returnunit (unitno)

!$OMP PARALLEL DO PRIVATE (ixy,edg)
    do ipn=1,nip
      do ixy=1,2
        do edg=1,npp
! --- going counterclockwise, edges start at location -trlpt-, end at -ldgpt-.
! --- hence, ldgpt of edge n agrees with trlpt of edge n+1.
          trlpt(edg,ixy,ipn) = corners(edg,1,ixy,ipn)
          ldgpt(edg,ixy,ipn) = corners(edg,2,ixy,ipn)
        end do
      end do
    end do
!$OMP END PARALLEL DO

    deallocate(corners)

!SMS$SERIAL END

!SMS$EXCHANGE(trlpt,ldgpt)

  end subroutine rdcnfig


  subroutine cornwgt      ! compute interpolation weights
    use global_bounds,   only: ims, ime, ips, ipe, ihe
    use module_control,  only: nip,npp
    use fimnamelist,     only: PrintIpnDiag,eqwgt
    use module_constants,only: lat,lon,deg_lat,deg_lon,prox,nprox,  &
      perm,idxtrl,idxldg,wgtldg,wgttrl,  &
      nedge,permedge,trlpt,ldgpt
    use findmaxmin2

    real :: ang1tr(npp,ims:ime),ang2tr(npp,ims:ime),ang3tr(npp,ims:ime)
    real :: ang1ld(npp,ims:ime),ang2ld(npp,ims:ime),ang3ld(npp,ims:ime)

    real*8 :: vec0(3),vec1(3),vec2(3),vec3(3),vec1x2,vec2x3,vec3x1,  &
      sumvec,angl1,angl2,angl3,sumang
    real   :: tmptr(npp),tmpld(npp)
    integer ipn,ixy,edgcount,edg,edg1,i1,i2,i3l,i3t
    logical vrbos
    character(len=4) string
    real   ,parameter :: rad2dg=57.2957795, athird=1./3.
    integer :: n

    print *,'entering cornwgt...'

!sms$ignore begin

! Halo comp means horizontal loop indices are ips,ihe
!$OMP PARALLEL DO PRIVATE (vrbos,edgcount,edg,ipn,ixy,edg1,n,    &
!$OMP&                     i1,i2,i3l,i3t,vec0,vec1,vec2,vec3,vec1x2,  &
!$OMP&                     vec2x3,vec3x1,sumvec,angl1,angl2,angl3,  &
!$OMP&                     sumang,tmptr,tmpld) SCHEDULE (runtime)
    do ipn=ips,ihe
      vrbos = (ipn.le.ipe .and.      ipn .eq.PrintIpnDiag) .or.  &
        (ipn.gt.ipe .and. perm(ipn).eq.PrintIpnDiag)
      if (vrbos) print '(2(a,i8)/2(a,2i8))','now operating on ipnglb=',  &
        perm(ipn),'  ipnloc=',ipn,'index range',ips,ipe,    &
        '     halo range',ipe+1,ihe
      idxtrl(:,:,ipn)=-999
      idxldg(:,:,ipn)=-999
      wgttrl(:,:,ipn)=1.e33
      wgtldg(:,:,ipn)=1.e33

! --- compensate for an apparent edge index shift in the 'corners' array
      do ixy=1,2
        tmptr(:) = trlpt(:,ixy,ipn)
        tmpld(:) = ldgpt(:,ixy,ipn)
        do edg=1,nprox(ipn)
          edg1 = mod(edg,nprox(ipn))+1
          trlpt(edg1,ixy,ipn) = tmptr(edg)
          ldgpt(edg1,ixy,ipn) = tmpld(edg)
        end do
      end do

      if (vrbos) then
        print '(2(a,i8),a,i2,a)','ipnglb=',perm(ipn),      &
          ' (ipnloc=',ipn,') has',nedge(ipn),' neighbors:'
        print '(a,6i9)','glb:',(perm(prox(permedge(edgcount,ipn),ipn)),  &
          edgcount=1,nedge(ipn))
        print '(a,6i9)','loc:',     (prox(permedge(edgcount,ipn),ipn),  &
          edgcount=1,nedge(ipn))
        print 100,'trailing/leading corner points at ipn=',perm(ipn),  &
          '  lat/lon=',deg_lat(ipn),deg_lon(ipn),      &
          (edg,(trlpt(edg,ixy,ipn)*rad2dg,ixy=1,2),      &
          (ldgpt(edg,ixy,ipn)*rad2dg,ixy=1,2),edg=1,nprox(ipn))
100     format (a,i8,a,2f8.2/(i3,2(f11.2,f7.2)))
      end if

! --- find the three cells meeting at the 2 corner points located at
! --- lat,lon=trlpt(edg,1,ipn),trlpt(edg,2,ipn) and
! --- lat,lon=ldgpt(edg,1,ipn),ldgpt(edg,2,ipn)

      i1=ipn
      do edgcount=1,nedge(ipn)
        edg=permedge(edgcount,ipn)
        i2 =prox(edg,ipn)
        i3l=prox(mod(edg             ,nprox(ipn))+1,ipn)
        i3t=prox(mod(edg-2+nprox(ipn),nprox(ipn))+1,ipn)

        if (vrbos) then
          print 101,'edg=',edg,'-trl corner point',      &
            (trlpt(edg,ixy,ipn)*rad2dg,ixy=1,2),'  surrounded by cells...',&
            perm( i1),'  lat/lon=',deg_lat( i1),deg_lon( i1),  &
            perm( i2),'  lat/lon=',deg_lat( i2),deg_lon( i2),  &
            perm(i3t),'  lat/lon=',deg_lat(i3t),deg_lon(i3t)
          print 101,'edg=',edg,'-ldg corner point',      &
            (ldgpt(edg,ixy,ipn)*rad2dg,ixy=1,2),'  surrounded by cells...',&
            perm( i1),'  lat/lon=',deg_lat (i1),deg_lon( i1),  &
            perm( i2),'  lat/lon=',deg_lat( i2),deg_lon( i2),  &
            perm(i3l),'  lat/lon=',deg_lat(i3l),deg_lon(i3l)
101       format (/a,i2,a,2f8.2,a/(i13,a,2f8.2))
        end if

! --- define vectors pointing at TRAILING icos corner point (vec0)
! -- and at the 3 cells meeting at that point (vec1...vec3)

        vec0(1)=cos(dble(trlpt(edg,1,ipn)))*cos(dble(trlpt(edg,2,ipn)))
        vec0(2)=cos(dble(trlpt(edg,1,ipn)))*sin(dble(trlpt(edg,2,ipn)))
        vec0(3)=sin(dble(trlpt(edg,1,ipn)))

        vec1(1)=cos(dble(lat(i1)))*cos(dble(lon(i1)))
        vec1(2)=cos(dble(lat(i1)))*sin(dble(lon(i1)))
        vec1(3)=sin(dble(lat(i1)))

        vec2(1)=cos(dble(lat(i2)))*cos(dble(lon(i2)))
        vec2(2)=cos(dble(lat(i2)))*sin(dble(lon(i2)))
        vec2(3)=sin(dble(lat(i2)))

        vec3(1)=cos(dble(lat(i3t)))*cos(dble(lon(i3t)))
        vec3(2)=cos(dble(lat(i3t)))*sin(dble(lon(i3t)))
        vec3(3)=sin(dble(lat(i3t)))

        do n=1,3
          vec1(n)=vec1(n)-vec0(n)
          vec2(n)=vec2(n)-vec0(n)
          vec3(n)=vec3(n)-vec0(n)
        end do

! --- form cross products of 3 vector pairs to get sines of enclosed angles

        vec1x2=sinangl(vec1,vec2)
        vec2x3=sinangl(vec2,vec3)
        vec3x1=sinangl(vec3,vec1)

! --- if interpolation is based on plane surface fitted to the 3 icos points,
! --- sine of angle becomes interpolation weight for opposing cell

        if (eqwgt) then
          sumvec=1.
          wgttrl(edg,1,ipn)=athird
          wgttrl(edg,2,ipn)=athird
          wgttrl(edg,3,ipn)=athird
        else
          sumvec=vec1x2+vec2x3+vec3x1
          wgttrl(edg,1,ipn)=vec2x3/sumvec
          wgttrl(edg,2,ipn)=vec3x1/sumvec
          wgttrl(edg,3,ipn)=vec1x2/sumvec
        end if

        angl1=180.-asin(vec1x2)*rad2dg
        angl2=180.-asin(vec2x3)*rad2dg
        angl3=180.-asin(vec3x1)*rad2dg
        sumang=angl1+angl2+angl3

        ang1tr(edg,ipn)=angl1
        ang2tr(edg,ipn)=angl2
        ang3tr(edg,ipn)=angl3

        idxtrl(edg,1,ipn)=i1
        idxtrl(edg,2,ipn)=i2
        idxtrl(edg,3,ipn)=i3t

        if (vrbos) then
          print 102,'vectors, edge',edg,'-trl:',(vec1(n),vec2(n),vec3(n),n=1,3)
          print 103,'angles,  edge',edg,'-trl:',angl1,angl2,angl3,sumang
          print 104,'weights, edge',edg,'-trl:',perm(i1),vec2x3/sumvec,  &
            perm(i2),vec3x1/sumvec,perm(i3t),vec1x2/sumvec
102       format (/a,i2,a,3es14.4,2(/20x,3es14.4))
103       format (a,i2,a,4f11.3)
104       format (a,i2,a/(i9,f11.5))
        end if

! --- define vectors pointing at LEADING icos corner point (vec0)
! --- and at the 3 cells meeting at that point (vec1...vec3)

        vec0(1)=cos(dble(ldgpt(edg,1,ipn)))*cos(dble(ldgpt(edg,2,ipn)))
        vec0(2)=cos(dble(ldgpt(edg,1,ipn)))*sin(dble(ldgpt(edg,2,ipn)))
        vec0(3)=sin(dble(ldgpt(edg,1,ipn)))

        vec1(1)=cos(dble(lat(i1)))*cos(dble(lon(i1)))
        vec1(2)=cos(dble(lat(i1)))*sin(dble(lon(i1)))
        vec1(3)=sin(dble(lat(i1)))

        vec2(1)=cos(dble(lat(i2)))*cos(dble(lon(i2)))
        vec2(2)=cos(dble(lat(i2)))*sin(dble(lon(i2)))
        vec2(3)=sin(dble(lat(i2)))

        vec3(1)=cos(dble(lat(i3l)))*cos(dble(lon(i3l)))
        vec3(2)=cos(dble(lat(i3l)))*sin(dble(lon(i3l)))
        vec3(3)=sin(dble(lat(i3l)))

        do n=1,3
          vec1(n)=vec1(n)-vec0(n)
          vec2(n)=vec2(n)-vec0(n)
          vec3(n)=vec3(n)-vec0(n)
        end do

! --- form cross products of 3 vector pairs to get sines of enclosed angles

        vec1x2=sinangl(vec1,vec2)
        vec2x3=sinangl(vec2,vec3)
        vec3x1=sinangl(vec3,vec1)

! --- if interpolation is based on plane surface fitted to the 3 icos points,
! --- sine of angle becomes interpolation weight for opposing cell

        if (eqwgt) then
          sumvec=1.
          wgtldg(edg,1,ipn)=athird
          wgtldg(edg,2,ipn)=athird
          wgtldg(edg,3,ipn)=athird
        else
          sumvec=vec1x2+vec2x3+vec3x1
          wgtldg(edg,1,ipn)=vec2x3/sumvec
          wgtldg(edg,2,ipn)=vec3x1/sumvec
          wgtldg(edg,3,ipn)=vec1x2/sumvec
        end if

        angl1=180.-asin(vec2x3)*rad2dg
        angl2=180.-asin(vec3x1)*rad2dg
        angl3=180.-asin(vec1x2)*rad2dg
        sumang=angl1+angl2+angl3

        ang1ld(edg,ipn)=angl1
        ang2ld(edg,ipn)=angl2
        ang3ld(edg,ipn)=angl3

        idxldg(edg,1,ipn)=i1
        idxldg(edg,2,ipn)=i2
        idxldg(edg,3,ipn)=i3l

        if (vrbos) then
          print 102,'vectors, edge',edg,'-ldg:',(vec1(n),vec2(n),vec3(n),n=1,3)
          print 103,'angles,  edge',edg,'-ldg:',angl1,angl2,angl3,sumang
          print 104,'weights, edge',edg,'-ldg:',perm(i1),vec2x3/sumvec,  &
            perm(i2),vec3x1/sumvec,perm(i3l),vec1x2/sumvec
        end if
      end do      ! loop over edges

      if (vrbos) then
        do edgcount=1,nedge(ipn)
          edg=permedge(edgcount,ipn)
          print '(/a,i8,a/a)','ipn=',perm(ipn),'  cornwgt output:',  &
            '        trl_index:glb/loc   weight        ldg_index:glb/loc   weight'
          print 105,(edg,perm(idxtrl(edg,n,ipn)),idxtrl(edg,n,ipn),  &
            wgttrl(edg,n,ipn),          &
            perm(idxldg(edg,n,ipn)),idxldg(edg,n,ipn),    &
            wgtldg(edg,n,ipn),n=1,3)
105       format (3('edg=',i1,i11,i8,f10.5,i16,i8,f10.5/))
! --- make line printer plots of points/weights in proper geographic context
          print '(a,i8,i5,a)','next plot: weights for ipn,edg=',  &
            perm(ipn),edg,'-trl'
          call locplt(trlpt(edg,1,ipn),trlpt(edg,2,ipn),    &
            lat(idxtrl(edg,1,ipn)),lon(idxtrl(edg,1,ipn)),    &
            perm(idxtrl(edg,1,ipn)),    wgttrl(edg,1,ipn),    &
            lat(idxtrl(edg,2,ipn)),lon(idxtrl(edg,2,ipn)),    &
            perm(idxtrl(edg,2,ipn)),    wgttrl(edg,2,ipn),    &
            lat(idxtrl(edg,3,ipn)),lon(idxtrl(edg,3,ipn)),    &
            perm(idxtrl(edg,3,ipn)),    wgttrl(edg,3,ipn))
          print '(a,i8,i5,a)','next plot: weights for ipn,edg=',  &
            perm(ipn),edg,'-ldg'
          call locplt(ldgpt(edg,1,ipn),ldgpt(edg,2,ipn),    &
            lat(idxldg(edg,1,ipn)),lon(idxldg(edg,1,ipn)),    &
            perm(idxldg(edg,1,ipn)),    wgtldg(edg,1,ipn),    &
            lat(idxldg(edg,2,ipn)),lon(idxldg(edg,2,ipn)),    &
            perm(idxldg(edg,2,ipn)),    wgtldg(edg,2,ipn),    &
            lat(idxldg(edg,3,ipn)),lon(idxldg(edg,3,ipn)),    &
            perm(idxldg(edg,3,ipn)),    wgtldg(edg,3,ipn))
        end do
      end if
    end do      ! horiz. loop

#ifdef DEBUGPRINT
    do edg=1,5
      write (string,'(a,i1)') 'edg',edg
      call findmxmn2(ang1tr,npp,nip,edg,'(cornwgt) ang1tr '//string)
      call findmxmn2(ang2tr,npp,nip,edg,'(cornwgt) ang2tr '//string)
      call findmxmn2(ang3tr,npp,nip,edg,'(cornwgt) ang3tr '//string)
      call findmxmn2(ang1ld,npp,nip,edg,'(cornwgt) ang1ld '//string)
      call findmxmn2(ang2ld,npp,nip,edg,'(cornwgt) ang2ld '//string)
      call findmxmn2(ang3ld,npp,nip,edg,'(cornwgt) ang3ld '//string)
    end do
#endif

!sms$ignore end

    print *,'...exiting cornwgt'
    return
  end subroutine cornwgt


  subroutine locplt(lat0,lon0,lat1,lon1,ipn1,wgt1,  &
    lat2,lon2,ipn2,wgt2,  &
    lat3,lon3,ipn3,wgt3)

! --- given a 'target' point in lat/lon space and 3 nearby points, plot
! --- the locations of the 3 nearby points in a projection centered on
! --- the target point. also plotted are interpolations weights associated
! --- with each point.

! --- index 0 = target point, indices 1,2,3 = surrounding points
    real,   intent(IN) :: lat0,lon0,lat1,lon1,lat2,lon2,lat3,lon3
    real,   intent(IN) :: wgt1,wgt2,wgt3    ! interpolation weights
    integer,intent(IN) :: ipn1,ipn2,ipn3    ! icos point indices
    real*8 vec0(3),vec1(3),vec2(3),vec3(3),veci(3),vecj(3),    &
      xx(3),yy(3),rnorm
    real range
    integer idim,jdim,i,j,k,n
    real,parameter :: ratio=2.0    ! line spacing / character width
    real,parameter :: rad2dg=57.2957795
    character(len=1),allocatable :: char(:,:)
    character(len=6) string1,string2

!sms$ignore begin

!  print '(a/(f9.3,3x,3f9.3))','entering locplt with points 0...4:',  &
!    lat0*rad2dg,lat1*rad2dg,lat2*rad2dg,lat3*rad2dg,      &
!    lon0*rad2dg,lon1*rad2dg,lon2*rad2dg,lon3*rad2dg
!  print '(2a/3i9/(3f9.4))','(locplt) lat/lon increments',    &
!    ' rel.to point 0:',ipn1,ipn2,ipn3          &
!    ,(lat1-lat0)*rad2dg,(lat2-lat0)*rad2dg,(lat3-lat0)*rad2dg    &
!    ,(lon1-lon0)*rad2dg,(lon2-lon0)*rad2dg,(lon3-lon0)*rad2dg

! --- vec0: target point, to become origin of projection:
    vec0(1)=cos(dble(lat0))*cos(dble(lon0))
    vec0(2)=cos(dble(lat0))*sin(dble(lon0))
    vec0(3)=sin(dble(lat0))

! --- vec1,vec2,vec3: surrounding points:
    vec1(1)=cos(dble(lat1))*cos(dble(lon1))
    vec1(2)=cos(dble(lat1))*sin(dble(lon1))
    vec1(3)=sin(dble(lat1))

    vec2(1)=cos(dble(lat2))*cos(dble(lon2))
    vec2(2)=cos(dble(lat2))*sin(dble(lon2))
    vec2(3)=sin(dble(lat2))

    vec3(1)=cos(dble(lat3))*cos(dble(lon3))
    vec3(2)=cos(dble(lat3))*sin(dble(lon3))
    vec3(3)=sin(dble(lat3))

!  print '(a,4f9.3/(13x,4f9.3))','orig.vectors:',      &
!      (vec0(n),vec1(n),vec2(n),vec3(n),n=1,3)
!  print '(a,3es12.4)','radial distances to target point:',    &
!      sqrt(dot_product(vec1-vec0,vec1-vec0)),        &
!      sqrt(dot_product(vec2-vec0,vec2-vec0)),        &
!      sqrt(dot_product(vec3-vec0,vec3-vec0))

! --- construct vectors veci,vecj spanning the plane normal to vec0
    rnorm=vec0(1)**2+vec0(2)**2
    rnorm=1./sqrt(rnorm*(rnorm+vec0(3)**2))
    veci(1)=(vec0(1)*vec0(3))*rnorm
    veci(2)=(vec0(2)*vec0(3))*rnorm
    veci(3)=-(vec0(1)**2+vec0(2)**2)*rnorm

    rnorm=1./sqrt(vec0(1)**2+vec0(2)**2)
    vecj(1)=-vec0(2)*rnorm
    vecj(2)= vec0(1)*rnorm
    vecj(3)=0.

!  print '(a,2f9.3/(14x,2f9.3))','trans.vectors:',      &
!      (veci(n),vecj(n),n=1,3)
!  print '(a,2f11.8,es11.3)','check normalization and orthogonality:',  &
!   dot_product(veci,veci),dot_product(vecj,vecj),dot_product(veci,vecj)

! --- represent vec1...vec3 in new coordinate system
    xx(1)=dot_product(vec1,veci)
    yy(1)=dot_product(vec1,vecj)
    xx(2)=dot_product(vec2,veci)
    yy(2)=dot_product(vec2,vecj)
    xx(3)=dot_product(vec3,veci)
    yy(3)=dot_product(vec3,vecj)

!  print '(a,3es12.4)','radial distances to target point:',    &
!      (sqrt(xx(n)**2+yy(n)**2),n=1,3)

!  print 100,'lat,lon=',lat1*rad2dg,lon1*rad2dg,' => x,y=',xx(1),yy(1)
!  print 100,'lat,lon=',lat2*rad2dg,lon2*rad2dg,' => x,y=',xx(2),yy(2)
!  print 100,'lat,lon=',lat3*rad2dg,lon3*rad2dg,' => x,y=',xx(3),yy(3)
100 format (a,2f8.2,a,2f9.4)

! --- plot results
    range=max(abs(xx(1)),abs(xx(2)),abs(xx(3)),        &
      abs(yy(1)),abs(yy(2)),abs(yy(3)))*1.2
    jdim=75
    idim=jdim/ratio
!  print '(a,f9.4,2i5)','range,idim,jdim=',range,idim,jdim
    allocate (char(idim,jdim))

    char(:,:)=' '
    do n=1,3
      i=.5*(1.+xx(n)/range)*float(idim)+1.
      j=.5*(1.+yy(n)/range)*float(jdim)+1.
!   print '(a,2f9.4,2i5)','x,y,i,j=',xx(n),yy(n),i,j
!   write (char(i,j),'(i1)') n
      char(i,j)='X'
      if      (n.eq.1) then
        write (string1,'(i6)') ipn1
        write (string2,'(f6.4)') wgt1
      else if (n.eq.2) then
        write (string1,'(i6)') ipn2
        write (string2,'(f6.4)') wgt2
      else if (n.eq.3) then
        write (string1,'(i6)') ipn3
        write (string2,'(f6.4)') wgt3
      end if
      do k=1,6
        char(i-1,j-5+k)=string1(k:k)
        char(i+1,j-5+k)=string2(k:k)
      end do
    end do
    i=.5*float(idim)+1.
    j=.5*float(jdim)+1.
    char(i,j)='X'
!  print '(a,2f9.4,2i5)','x,y,i,j=',0.,0.,i,j

    do j=1,jdim
      char(   1,j)='-'
      char(idim,j)='-'
    end do
    do i=1,idim
      char(i,   1)='|'
      char(i,jdim)='|'
    end do
    do i=1,idim
      print '(99a1)',(char(i,j),j=1,jdim)
    end do
    call flush(6)
    return

!sms$ignore end

  end subroutine locplt

  real*8 function sinangl(vec1,vec2)
! --- compute the sine of the angle enclosed by 2 vectors
    real*8,intent(IN) :: vec1(3),vec2(3)
    sinangl=sqrt(((vec1(2)*vec2(3) - vec2(2)*vec1(3))**2    &
      +(vec1(3)*vec2(1) - vec2(3)*vec1(1))**2    &
      +(vec1(1)*vec2(2) - vec2(1)*vec1(2))**2)    &
      /((vec1(1)**2 + vec1(2)**2 + vec1(3)**2)    &
      *(vec2(1)**2 + vec2(2)**2 + vec2(3)**2)))
    return
  end function sinangl
end module module_cornwgt
