!*********************************************************************
!       sfc2gg
!       sfc-to-global program for fim global model
!       - J-W Bao, Jin Lee, Ning Wang   2007
!           - initial version
!       - Stan Benjamin                 May 2008
!           - correction to change to nearest-neighbor interpolation
!             for land-surface 2-d variables instead of
!             previous bilinear-interpolation, which had caused
!             many land points to be erroneously set as water points
!       - Ning Wang                     Aug 2008
!           - Replaced a call to subroutine bilinear_init() with  
!             a call to nn_init(), to work with the new version of  
!             slint library. 
!       - Ning Wang
!           - Replaced computation of nip with a more general code.
!       - Mike Fiorino 20130207
!           - add alt_land option
!
!*********************************************************************

program sfc2gg
  use sfcio_module
  use module_control, only : control,nip
  use fimnamelist   , only : yyyymmddhhmm,glvl,curve,mtnvar_file,gfsltln_file	&
			    ,atm_ic
  implicit none

  integer(sfcio_intkind),parameter:: lusfc=11,luggg=51,luctl=52
  integer(sfcio_intkind)          :: irets
  integer                         :: iret,n
  CHARACTER(len=9 )               :: jdate 
  CHARACTER(len=2 )               :: hh
  CHARACTER(len=80)               :: sfcanlFile
  type(sfcio_head),allocatable    :: head(:)
  type(sfcio_data)                :: data

  call control()
    
  if (atm_ic == 1) then                     ! GFS
    call GetJdate(yyyymmddhhmm,jdate)       ! Julian date conversion
    hh = yyyymmddhhmm(9:10)
    sfcanlFile = jdate // ".gfs.t" // hh // "z.sfcanl"
  elseif (atm_ic == 2) then                 ! CFSR
    sfcanlFile = "sfcanl.gdas." // yyyymmddhhmm(1:10)
  end if
  print *,' get sfcanl data from ',trim(sfcanlFile)

  open(luctl,file='sfcanl.ctl',status='replace',iostat=iret)
  if (iret.ne.0) then
! fail here with informative message?
  endif
  allocate(head(1))
  do n=1,1
    print *,'calling sfcio_srohdc ...'
    call sfcio_srohdc(lusfc,sfcanlFile,head(n),data,irets)
    if(head(n)%latb.ne.head(1)%latb.or.&
      head(n)%lonb.ne.head(1)%lonb.or.&
      head(n)%lsoil.ne.head(1)%lsoil.or.&
      head(n)%ivs.ne.head(1)%ivs) then
      call errmsg('sfc2gg: incompatible data in file sfcanl.ctl')
!      call errexit(2)		! errexit doesn't call MPI_ABORT -> use STOP
      STOP
    endif
    print *,'calling sfc2gg1 ...'
    call sfc2gg1(luggg,head(n),data,glvl,nip,curve,mtnvar_file,gfsltln_file)
    print *,'calling sfcio_axdata'
    call sfcio_axdata(data,irets)
  enddo
  print *,'calling sfc2gg2 ...'
  call sfc2gg2(luctl,1,head,'sfcanl.ieee')

  print *,'... exiting sfc2gg'

contains

  subroutine eusage
    implicit none
    call errmsg('Usage: sfc2gg sfcfile(s) gggfile ctlfile')
  end subroutine eusage

end program sfc2gg

subroutine sfc2gg1(luggg,head,data,glvl,nip,curve,mtnvar_file,gfsltln_file)

  use headers,only: testcurveheader,testglvlheader, &
                    writecurveheader,writeglvlheader
  use sfcio_module
  use fimnamelist,only: alt_land
  use mdul_mkmodisland
  USE slint, ONLY: nn_init, nn_int
  use ssfc2icos_mtnvar,only: read_mtnvar

  implicit none
! integer, parameter :: imax=1152
! integer, parameter :: jmax= 576

  integer,allocatable :: setO2L(:,:)
  integer,allocatable :: setL2O(:,:)

  integer         :: imax
  integer         :: jmax
  integer         ,intent(in) :: luggg
  type(sfcio_head),intent(in) :: head
  type(sfcio_data),intent(in) :: data
  integer         ,intent(in) :: glvl
  integer         ,intent(in) :: nip
  integer         ,intent(in) :: curve

  integer l
  integer ipn,irc,nset,no2l,nl2o
!
  CHARACTER(len=80)   :: grid_file2,mtnvar_file, gfsltln_file
  CHARACTER(len=80)   :: gfsltln_file_surface, gfsltln_file_mtnvar
  real(4) icos2d(nip)

  real :: st3d(4,nip)        ! soil temperature
  real :: sm3d(4,nip)        ! soil moisture
  real :: slc3d(4,nip)       ! liquid soil moisture

  real :: ts2d(nip)          ! skin temperature
  real :: sheleg2d(nip)
  real :: tg32d(nip)
  real :: zorl2d(nip)
  real :: cv2d(nip)
  real :: cvb2d(nip)
  real :: cvt2d(nip)
  real :: alvsf2d(nip)
  real :: alvwf2d(nip)
  real :: alnsf2d(nip)
  real :: alnwf2d(nip)
  real :: slmsk2d(nip)
  real :: vfrac2d(nip)
  real :: canopy2d(nip)
  real :: f10m2d(nip)
  real :: t2m2d(nip)
  real :: q2m2d(nip)
  real :: vtype2d(nip)
  real :: stype2d(nip)
  real :: facsf2d(nip)
  real :: facwf2d(nip)
  real :: uustar2d(nip)
  real :: ffmm2d(nip)
  real :: ffhh2d(nip)
  real :: work2d(nip)
  real :: hice2d(nip)
  real :: fice2d(nip)
  real     tprcp2d(nip)
  real     srflag2d(nip)
  real     snwdph2d(nip)
  real     slc2d(nip)
  real     shdmin2d(nip)
  real     shdmax2d(nip)
  real     slope2d(nip)
  real     snoalb2d(nip)

  integer :: mvar = 14
  real ,allocatable:: mdrag3d(:,:,:), mdrag(:,:)

  real,    allocatable :: stemp2d(:),vtemp2d(:)

  integer idx

  grid_file2 = "glvl.dat"
  OPEN(66,file=grid_file2,status='old',form='unformatted')
  call TestGlvlHeader (66,grid_file2,'sfc2gg1',glvl)
  call TestCurveHeader(66,grid_file2,'sfc2gg1',curve)

  imax=head%lonb
  jmax=head%latb
  IF (imax == 1152 .AND. jmax == 576) THEN
    gfsltln_file = "gfsltln_t382.dat"
  ELSE IF (imax == 1760 .AND. jmax == 880) THEN
    gfsltln_file = "gfsltln_t574.dat"
  ELSE IF (imax == 2640 .AND. jmax == 1320) THEN
    gfsltln_file = "gfsltln_t878.dat"
  ELSE IF (imax == 3072 .AND. jmax == 1536) THEN
    gfsltln_file = "gfsltln_t1534.dat"
  END IF 
  gfsltln_file_surface = gfsltln_file

  PRINT*, 'gfsltln_file for .sfcanl data:', gfsltln_file

  CALL nn_init(gfsltln_file_surface, imax*jmax, 66, nip)
  CLOSE(66)

!
  CALL nn_int (data%tsea, ts2d)
!
  do l=1,head%lsoil
    CALL nn_int (data%smc(:,:,l), icos2d)
    do ipn=1,nip
      sm3d(l,ipn)=icos2d(ipn)
    end do
  enddo
  do l=1,head%lsoil
    CALL nn_int (data%slc(:,:,l), icos2d)
    do ipn=1,nip
      slc3d(l,ipn)=icos2d(ipn)
    end do
  enddo
! do l=1,head%lsoil
!   CALL nn_int (data%smc(:,:,l), icos2d)
!   do ipn=1,nip
!   slc3d(l,ipn)=icos2d(ipn)
!   end do
! enddo
  CALL nn_int (data%sheleg, sheleg2d)
  do l=1,head%lsoil
    CALL nn_int (data%stc(:,:,l), icos2d)
    do ipn=1,nip
      st3d(l,ipn)=icos2d(ipn)
    end do
  enddo
  CALL nn_int (data%tg3, tg32d)
  CALL nn_int (data%zorl, zorl2d)
! CALL bl_int (data%cv, cv2d)
! CALL bl_int (data%cvb, cvb2d)
! CALL bl_int (data%cvt, cvt2d)
  cv2d=0.   !zero out cv2d suggested by Bao
  cvb2d=0.  !zero out cv2d suggested by Bao
  cvt2d=0.  !zero out cv2d suggested by Bao
  CALL nn_int (data%alvsf, alvsf2d)
  CALL nn_int (data%alvwf, alvwf2d)
  CALL nn_int (data%alnsf, alnsf2d)
  CALL nn_int (data%alnwf, alnwf2d)
  CALL nn_int (data%slmsk, slmsk2d)
  CALL nn_int (data%vfrac, vfrac2d)
  CALL nn_int (data%canopy,canopy2d)
  CALL nn_int (data%f10m , f10m2d)
  CALL nn_int (data%t2m , t2m2d)
  CALL nn_int (data%q2m , q2m2d)
  CALL nn_int (data%vtype, vtype2d)
  CALL nn_int (data%stype, stype2d)
  CALL nn_int (data%facsf, facsf2d)
  CALL nn_int (data%facwf, facwf2d)
  CALL nn_int (data%uustar,uustar2d)
  CALL nn_int (data%ffmm , ffmm2d)
  CALL nn_int (data%ffhh , ffhh2d)
  CALL nn_int (data%hice , hice2d)
  CALL nn_int (data%fice , fice2d)
  CALL nn_int (data%tprcp, tprcp2d)
  CALL nn_int (data%srflag, srflag2d)
  CALL nn_int (data%snwdph, snwdph2d)
  CALL nn_int (data%slc, slc2d)
  CALL nn_int (data%shdmin, shdmin2d)
  CALL nn_int (data%shdmax, shdmax2d)
  CALL nn_int (data%slope, slope2d)
  CALL nn_int (data%snoalb, snoalb2d)

  IF (mtnvar_file == 'global_mtnvar.t382') THEN
    gfsltln_file_mtnvar = 'gfsltln_t382.dat'
    imax = 1152; jmax= 576
  ELSE IF (mtnvar_file == 'global_mtnvar.t574') THEN
    gfsltln_file_mtnvar = 'gfsltln_t574.dat'
    imax = 1760; jmax= 880
  ELSE IF (mtnvar_file == 'global_mtnvar.t878') THEN
    gfsltln_file_mtnvar = 'gfsltln_t878.dat'
    imax = 2640; jmax= 1320
  ELSE IF (mtnvar_file == 'global_mtnvar.t1534') THEN
    gfsltln_file_mtnvar = 'gfsltln_t1534.dat'
    imax = 3072; jmax= 1536
  ENDIF
  PRINT*, 'gfsltln_file for mtnvar data:', gfsltln_file_mtnvar
  
  IF (gfsltln_file_mtnvar /= gfsltln_file_surface) THEN
    OPEN(66,file=grid_file2,status='old',form='unformatted')
    call TestGlvlHeader (66,grid_file2,'sfc2gg1',glvl)
    call TestCurveHeader(66,grid_file2,'sfc2gg1',curve)
    PRINT*, 'gfsltln_file for mtnvar data:', gfsltln_file_mtnvar
    CALL nn_init(gfsltln_file_mtnvar, imax*jmax, 66, nip)
    CLOSE(66)
  ENDIF

  allocate(mdrag3d(imax,jmax,mvar))
  allocate(mdrag(mvar,nip))
  write(6,*)'calling read_mtnvar'
  call read_mtnvar(mdrag3d,imax,jmax,mvar,mtnvar_file)
  write(6,*)'Done calling read_mtnvar'
  do l=1,mvar
    CALL nn_int (mdrag3d(:,:,l), icos2d)
    do ipn=1,nip
      mdrag(l,ipn)=icos2d(ipn)
    end do
  end do
  deallocate(mdrag3d)

  if (alt_land) then

!SMS$SERIAL (<vtype2d,INOUT> : default=ignore)  BEGIN

!mf 20130207 -- replace the gfs vtype2d with an analysis of modis
!mf 20130315 -- new version that calculates arrays to change gfs points to modis 
!            -- allocate arrays with indices to deconflict modis and gfs land vars
!
    nset=nip/50

    allocate(setO2L(nset,2), stat=irc)
    allocate(setL2O(nset,2), stat=irc)
    allocate(stemp2d(nip),   stat=irc)  ! temp arrays for deconfliction
    allocate(vtemp2d(nip),   stat=irc) 

!ret = gptlstart ('mktopo') -- m

! -- analayze modis land veg type to fim grid; and find arrays to deconflict
!
    write(6,*)'Calling mkmodisland'
    call mkmodisland(vtype2d,nip,setO2L,setL2O,nset,no2l,nl2o)
    write(6,*)'Done calling mkmodisland'

! -- deconflict 2d arrays
!
    call deconflictGfsland(tg32d,     stemp2d,setO2L,setL2O,nip,nset,no2l,nl2o)
    call deconflictGfsland(zorl2d,    stemp2d,setO2L,setL2O,nip,nset,no2l,nl2o)
    call deconflictGfsland(alvsf2d,   stemp2d,setO2L,setL2O,nip,nset,no2l,nl2o)
    call deconflictGfsland(alvwf2d,   stemp2d,setO2L,setL2O,nip,nset,no2l,nl2o)
    call deconflictGfsland(alnsf2d,   stemp2d,setO2L,setL2O,nip,nset,no2l,nl2o)
    call deconflictGfsland(alnwf2d,   stemp2d,setO2L,setL2O,nip,nset,no2l,nl2o)
    call deconflictGfsland(slmsk2d,   stemp2d,setO2L,setL2O,nip,nset,no2l,nl2o)
    call deconflictGfsland(vfrac2d,   stemp2d,setO2L,setL2O,nip,nset,no2l,nl2o)
    call deconflictGfsland(canopy2d,  stemp2d,setO2L,setL2O,nip,nset,no2l,nl2o)
    call deconflictGfsland(f10m2d,    stemp2d,setO2L,setL2O,nip,nset,no2l,nl2o)
    call deconflictGfsland(t2m2d,     stemp2d,setO2L,setL2O,nip,nset,no2l,nl2o)
    call deconflictGfsland(q2m2d,     stemp2d,setO2L,setL2O,nip,nset,no2l,nl2o)
!    call deconflictGfsland(vtype2d,   stemp2d,setO2L,setL2O,nip,nset,no2l,nl2o)
    call deconflictGfsland(stype2d,   stemp2d,setO2L,setL2O,nip,nset,no2l,nl2o)
    call deconflictGfsland(facsf2d,   stemp2d,setO2L,setL2O,nip,nset,no2l,nl2o)
    call deconflictGfsland(facwf2d,   stemp2d,setO2L,setL2O,nip,nset,no2l,nl2o)
    call deconflictGfsland(uustar2d,  stemp2d,setO2L,setL2O,nip,nset,no2l,nl2o)
    call deconflictGfsland(ffmm2d,    stemp2d,setO2L,setL2O,nip,nset,no2l,nl2o)
    call deconflictGfsland(ffhh2d,    stemp2d,setO2L,setL2O,nip,nset,no2l,nl2o)
    call deconflictGfsland(hice2d,    stemp2d,setO2L,setL2O,nip,nset,no2l,nl2o)
    call deconflictGfsland(fice2d,    stemp2d,setO2L,setL2O,nip,nset,no2l,nl2o)
    call deconflictGfsland(tprcp2d,   stemp2d,setO2L,setL2O,nip,nset,no2l,nl2o)
    call deconflictGfsland(srflag2d,  stemp2d,setO2L,setL2O,nip,nset,no2l,nl2o)
    call deconflictGfsland(snwdph2d,  stemp2d,setO2L,setL2O,nip,nset,no2l,nl2o)
    call deconflictGfsland(slc2d,     stemp2d,setO2L,setL2O,nip,nset,no2l,nl2o)
    call deconflictGfsland(shdmin2d,  stemp2d,setO2L,setL2O,nip,nset,no2l,nl2o)
    call deconflictGfsland(shdmax2d,  stemp2d,setO2L,setL2O,nip,nset,no2l,nl2o)
    call deconflictGfsland(slope2d,   stemp2d,setO2L,setL2O,nip,nset,no2l,nl2o)
    call deconflictGfsland(snoalb2d,  stemp2d,setO2L,setL2O,nip,nset,no2l,nl2o)

! -- deconflict 3d arrays
!
    do idx = 1,SIZE(st3d,1)
      do ipn=1,nip
        work2d(ipn) = st3d(idx,ipn)
      enddo
      call deconflictGfsland(work2d,  stemp2d,setO2L,setL2O,nip,nset,no2l,nl2o)
      do ipn=1,nip
        st3d(idx,ipn)=work2d(ipn)
      enddo
    enddo

    do idx = 1,SIZE(sm3d,1)
      do ipn=1,nip
        work2d(ipn) = sm3d(idx,ipn)
      enddo
      call deconflictGfsland(work2d,  stemp2d,setO2L,setL2O,nip,nset,no2l,nl2o)
      do ipn=1,nip
        sm3d(idx,ipn)=work2d(ipn)
      enddo
    enddo

    do idx = 1,SIZE(slc3d,1)
      do ipn=1,nip
        work2d(ipn) = slc3d(idx,ipn)
      enddo
      call deconflictGfsland(work2d,  stemp2d,setO2L,setL2O,nip,nset,no2l,nl2o)
      do ipn=1,nip
        slc3d(idx,ipn)=work2d(ipn)
      enddo
    enddo

    do idx = 1,SIZE(mdrag,1)
      do ipn=1,nip
        work2d(ipn) = mdrag(idx,ipn)
      enddo
      call deconflictGfsland(work2d,  stemp2d,setO2L,setL2O,nip,nset,no2l,nl2o)
      do ipn=1,nip
        mdrag(idx,ipn)=work2d(ipn)
      enddo
    enddo

!ret = gptlstop ('mktopo')
!SMS$SERIAL END
  endif

! --- open output file
  open (10,file="gfsfc.dat",form='unformatted')
  call WriteGlvlHeader (10,glvl )
  call WriteCurveHeader(10,curve)
  do idx = 1,SIZE(st3d,1)
    do ipn=1,nip
      work2d(ipn) = st3d(idx,ipn)
    enddo
    write(10) work2d
  enddo
  do idx = 1,SIZE(sm3d,1)
    do ipn=1,nip
      work2d(ipn) = sm3d(idx,ipn)
    enddo
    write(10) work2d
  enddo
  do idx = 1,SIZE(slc3d,1)
    do ipn=1,nip
      work2d(ipn) = slc3d(idx,ipn)
    enddo
    write(10) work2d
  enddo
  write(10) ts2d
  write(10) sheleg2d
  write(10) tg32d
  write(10) zorl2d
  write(10) cv2d
  write(10) cvb2d
  write(10) cvt2d
  write(10) alvsf2d
  write(10) alvwf2d
  write(10) alnsf2d
  write(10) alnwf2d
  write(10) slmsk2d
  write(10) vfrac2d
  write(10) canopy2d
  write(10) f10m2d
  write(10) t2m2d
  write(10) q2m2d
  write(10) vtype2d
  write(10) stype2d
  write(10) facsf2d
  write(10) facwf2d
  write(10) uustar2d
  write(10) ffmm2d
  write(10) ffhh2d
  write(10) hice2d
  write(10) fice2d
  write(10) tprcp2d
  write(10) srflag2d
  write(10) snwdph2d
  write(10) slc2d
  write(10) shdmin2d
  write(10) shdmax2d
  write(10) slope2d
  write(10) snoalb2d
  do idx = 1,SIZE(mdrag,1)
    do ipn=1,nip
      work2d(ipn) = mdrag(idx,ipn)
    enddo
    write(10) work2d
  enddo
  close(10) 
!
end subroutine sfc2gg1

subroutine sfc2gg2(luctl,nsfc,head,cfggg)
  use sfcio_module
  implicit none
  integer,intent(in):: luctl,nsfc
  type(sfcio_head),intent(in):: head(nsfc)
  character*(*) cfggg
  real(4),allocatable:: slat(:),wlat(:)
  integer idat(8),jdat(8),jhr
  real(4) rinc(5),rincin(5)
  integer idatin(8)
  character*10 cdat(8)
!  call w3movdat((/0.,head(1)%fhour,0.,0.,0./),&
!                (/head(1)%idate(4),head(1)%idate(2),head(1)%idate(3),0,&
!                  head(1)%idate(1),0,0,0/),idat)
  rincin    = 0.
  rincin(2) = head(1)%fhour
  idatin    = 0
  idatin(1) = head(1)%idate(4)
  idatin(2) = head(1)%idate(2)
  idatin(3) = head(1)%idate(3)
  idatin(5) = head(1)%idate(1)
  call w3movdat(rincin,idatin,idat)
  call w3pradat(idat,cdat)
  if(nsfc.gt.1) then
!    call w3movdat((/0.,head(2)%fhour,0.,0.,0./),&
!                  (/head(2)%idate(4),head(2)%idate(2),head(2)%idate(3),0,&
!                    head(2)%idate(1),0,0,0/),jdat)
    rincin(2) = head(2)%fhour
    idatin(1) = head(2)%idate(4)
    idatin(2) = head(2)%idate(2)
    idatin(3) = head(2)%idate(3)
    idatin(5) = head(2)%idate(1)
    call w3movdat(rincin,idatin,idat)
    call w3difdat(jdat,idat,2,rinc)
    jhr=nint(rinc(2))
  else
    jhr=12
  endif
  if(cfggg(1:1).eq.'/') then
    write(luctl,'("dset ",a)') cfggg
  else
    write(luctl,'("dset ^",a)') cfggg
  endif
  write(luctl,'("options yrev sequential")')
  write(luctl,'("undef -9.99E+33")')
  write(luctl,'("title sfc2gg")')
  write(luctl,'("xdef",i6," linear",2f12.6)') head(1)%lonb,0.d0,360.d0/head(1)%lonb
  allocate(slat(head(1)%latb),wlat(head(1)%latb))
  call splat(4,head(1)%latb,slat,wlat)
  write(luctl,'("ydef",i6," levels")') head(1)%latb
  write(luctl,'(5f12.6)') 180.d0/acos(-1.d0)*asin(dble(slat(head(1)%latb:1:-1)))
  write(luctl,'("zdef",i6," levels")') head(1)%lsoil
  write(luctl,'(5f12.6)') head(1)%zsoil
  write(luctl,'("tdef",i6," linear ",i2.2,"Z",i2.2,a3,i4.4,1x,i6,"hr")')&
    nsfc,idat(5),idat(3),cdat(2)(1:3),idat(1),jhr
  write(luctl,'("vars",i6)') 35
  write(luctl,'("tsea    ",i3," 99 surface temperature (K)")') 1
  write(luctl,'("smc     ",i3," 99 soil volumetric water content ()")') head(1)%lsoil
  write(luctl,'("sheleg  ",i3," 99 snow depth (m)")') 1
  write(luctl,'("stc     ",i3," 99 soil temperature (K)")') head(1)%lsoil
  write(luctl,'("tg3     ",i3," 99 deep soil temperature (K)")') 1
  write(luctl,'("zorl    ",i3," 99 roughness (cm)")') 1
  write(luctl,'("cv      ",i3," 99 convective cloud cover ()")') 1
  write(luctl,'("cvb     ",i3," 99 convective cloud bottom (kPa)")') 1
  write(luctl,'("cvt     ",i3," 99 convective cloud top (kPa)")') 1
  write(luctl,'("alvsf   ",i3," 99 albedo for visible scattered ()")') 1
  write(luctl,'("alvwf   ",i3," 99 albedo for visible beam ()")') 1
  write(luctl,'("alnsf   ",i3," 99 albedo for near-IR scattered ()")') 1
  write(luctl,'("alnwf   ",i3," 99 albedo for near-IR beam ()")') 1
  write(luctl,'("slmsk   ",i3," 99 sea-land-ice mask (0-sea, 1-land, 2-ice)")') 1
  write(luctl,'("vfrac   ",i3," 99 vegetation fraction ()")') 1
  write(luctl,'("canopy  ",i3," 99 canopy water (m)")') 1
  write(luctl,'("f10m    ",i3," 99 10-meter wind speed over lowest model wind speed ()")') 1
  write(luctl,'("vtype   ",i3," 99 vegetation type (integer 1-13)")') 1
  write(luctl,'("stype   ",i3," 99 soil type (integer 1-9)")') 1
  write(luctl,'("facsf   ",i3," 99 ???")') 1
  write(luctl,'("facwf   ",i3," 99 ???")') 1
  write(luctl,'("uustar  ",i3," 99 ???")') 1
  write(luctl,'("ffmm    ",i3," 99 ???")') 1
  write(luctl,'("ffhh    ",i3," 99 ???")') 1
  write(luctl,'("hice    ",i3," 99 ???")') 1
  write(luctl,'("fice    ",i3," 99 ???")') 1
  write(luctl,'("tprcp   ",i3," 99 ???")') 1
  write(luctl,'("srflag  ",i3," 99 ???")') 1
  write(luctl,'("snwdph  ",i3," 99 ???")') 1
  write(luctl,'("slc     ",i3," 99 ???")') head(1)%lsoil
  write(luctl,'("shdmin  ",i3," 99 ???")') 1
  write(luctl,'("shdmax  ",i3," 99 ???")') 1
  write(luctl,'("slope   ",i3," 99 ???")') 1
  write(luctl,'("snoalb  ",i3," 99 ???")') 1
  write(luctl,'("orog    ",i3," 99 orography (m)")') 1
  write(luctl,'("endvars")')
end subroutine sfc2gg2
