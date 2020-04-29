module hycom_thermf

! --- - - - - - - - - - - - - - - - - - - -
! --- manage high-freq. CORE forcing data
! --- - - - - - - - - - - - - - - - - - - -

  use module_control,  only: nip,dt
  use fimnamelist,     only: glvl,ArchvTimeUnit,itest,diag_intvl,	&
      ann_core,inicondir,do_radcor,do_pcpcor,num_mthly_fields,		&
      num_daily_fields,num_6hrly_fields
  use module_constants,only: area,deg_lat,deg_lon,perm
  use hycom_constants, only: wet,thref,airdns,evaplh,csubp,KelvC,	&
      odepth,rhoice,spcifh,batrop,huuge,numflds,flnm_in,varname,ncid1,	&
      ixuwi,ixvwi,ixtem,ixvap,ixlng,ixsho,ixpcp,ixsss,ixrff,land_spval
  use hycom_control,   only: bclin_frq
  use hycom_variables, only: wspd,taux,tauy,ustar,temp,saln,dp,		&
      covice,temice,thkice,salnow,srforc,airtmn,pmne,prcp,salnow,	&
      uwnd,vwnd,airt,vpmx,srfx,ocnarea,radcor,pcpcor,radtot,pcptot,	&
      wm0,wm1,wm2,wm3,wd0,wd1,wd2,wd3,ws0,ws1,ws2,ws3,        		&
      lm0,lm1,lm2,lm3,ld0,ld1,ld2,ld3,ls0,ls1,ls2,ls3
  use module_sfc_variables,only: us2d,hf2d,qf2d,slmsk2d,rsds,rsus,rlds,rlus
  use stencilprint
  use findmaxmin1
  use findmaxmin2
! use ersatzpipe
  implicit none

  integer :: nfld,nrec,n,i
  character chrglvl*2,flnm*80

  contains

  subroutine opnfor

! --- - - - - - - - - - - - - - - - - - - -
! --- manage high-freq. CORE forcing data
! --- - - - - - - - - - - - - - - - - - - -

  use netcdf

!SMS$DISTRIBUTE(dh,1) BEGIN
  real data1(nip),data2(nip)
!SMS$DISTRIBUTE END
  real srial(nip)
  character flnm*120,chrglvl*2,string*16
  integer:: id_tdm,tz,nfld,nrec,i
  integer :: start(2)	! starting indices for the variable in the netcdf file
  integer :: kount(2)	! lengths for the variable in the netcdf file

  varname(ixuwi)='U_10_MOD'		; flnm_in(ixuwi)='u_10'
  varname(ixvwi)='V_10_MOD'		; flnm_in(ixvwi)='v_10'
  varname(ixtem)='T_10_MOD'		; flnm_in(ixtem)='t_10'
  varname(ixvap)='Q_10_MOD'		; flnm_in(ixvap)='q_10'
  varname(ixlng)='LWDN_MOD'		; flnm_in(ixlng)='lwdn'
  varname(ixsho)='SWDN_MOD'		; flnm_in(ixsho)='swdn'
  varname(ixpcp)='PRECIP'		; flnm_in(ixpcp)='precip'
  varname(ixsss)='SALT'			; flnm_in(ixsss)='sss'
  varname(ixrff)='Foxx_o_roff'		; flnm_in(ixrff)='runoff'

! -- initialize
   lm0=-99
   ld0=-99
   ls0=-99
   srforc(:,:,:)=-999999.

!SMS$SERIAL (<wet,area,perm,IN>,<rsds,rlds,prcp,uwnd,vwnd,airt,vpmx,airtmn,salnow,INOUT>,<data1,data2,radtot,pcptot,OUT> : default=ignore) BEGIN
  kount = (/nip,1/)
  write (chrglvl,'(a,i1)') 'g',glvl

  do nfld=1,numflds
   if (num_mthly_fields==12) then
    flnm=trim(inicondir)//'input_'//chrglvl//'/'//trim(flnm_in(nfld))//'_core1_'//chrglvl//'.nc'
   elseif (num_mthly_fields==12*60) then
    flnm=trim(inicondir)//'input_'//chrglvl//'/'//trim(flnm_in(nfld))//'_core2_'//chrglvl//'.nc'
   else 
    stop '(opnfor) illegal num_mthly_fields'
   end if
   write (*,'(a,i3,2a)') 'opening forcing field',nfld,':  ',trim(flnm)
   call errhandl (nf90_open (path=trim(flnm), mode=NF90_NOWRITE, ncid=ncid1(nfld)))   ! netcdf3
   write (*,'(2a)') 'checking contents of ',trim(flnm_in(nfld))
   call errhandl (nf90_inq_dimid (ncid1(nfld), 'time',id_tdm))
   call errhandl (nf90_inquire_dimension (ncid1(nfld), id_tdm, len=tz))
   write (*,'(a,i7,a)') 'file contains',tz,'  records'

   if (ann_core) then		! if ann_core
      start = (/1,tz/)
      call readcdf_1d(ncid1(nfld),varname(nfld),srial,nip,start,kount)
      data1(:)=land_spval
      do i=1,nip
        if (wet(i)>0) data1(i)=srial(perm(i))
      end do
        if (nfld==ixpcp)   prcp(:)=data1(:)
        if (nfld==ixsss) salnow(:)=data1(:)
        if (nfld==ixsho)   rsds(:)=data1(:)
        if (nfld==ixlng)   rlds(:)=data1(:)
        if (nfld==ixuwi)   uwnd(:)=data1(:)
        if (nfld==ixvwi)   vwnd(:)=data1(:)
        if (nfld==ixtem)   airt(:)=data1(:)
        if (nfld==ixvap)   vpmx(:)=data1(:)
   end if	! if ann_core
  end do	! do nfld

!  do k=1,366,73
!   write (string,'(a,i3)') 'dy=',k
!   call findmxmn2(swdn,  366,nip,k,'(forcing) swdn '//string,wet)
!   call findmxmn2(lwdn,  366,nip,k,'(forcing) lwdn '//string,wet)
!   call findmxmn2(airtem,366,nip,k,'(forcing) airt '//string,wet)
!   call findmxmn2(vpmix ,366,nip,k,'(forcing) vpmx '//string,wet)
!   call findmxmn2(uwind ,366,nip,k,'(forcing) uwnd '//string,wet)
!   call findmxmn2(vwind ,366,nip,k,'(forcing) vwnd '//string,wet)
!   print *
!  enddo
! if (ann_core) then		! if ann_core
!   call findmxmn1(prcp,nip,'(ann forcing) prcp ',wet)
!   call findmxmn1(salnow, nip,'(ann forcing) sss ',wet)
!   call findmxmn1(uwnd,nip,'(ann forcing) uwnd ',wet)
!   call findmxmn1(vwnd,nip,'(ann forcing) vwnd ',wet)
!   call findmxmn1(rsds,nip,'(ann forcing) rsds ',wet)
!   call findmxmn1(rlds,nip,'(ann forcing) rlds ',wet)
!   call findmxmn1(airt,nip,'(ann forcing) airt ',wet)
!   call findmxmn1(vpmx,nip,'(ann forcing) vpmx ',wet)
! end if	! if ann_core

! --- compute minimum air temperature for ice initialization

  nfld=ixtem
  do nrec=1,1460
   start = (/1,nrec/)
   call readcdf_1d(ncid1(nfld),varname(nfld),srial,nip,start,kount)
   do i=1,nip
    if (wet(i)==1) then
      data1(i)=srial(perm(i))-KelvC	! deg C
      airtmn(i)=min(airtmn(i),data1(i))	! deg C
    end if
   end do		! do i
  end do		! do nrec

! --- compute global integrals of precip and long+shortwave radiation

  pcptot=0.
  radtot=0.
  do nrec=1,365		! average 1 yr of daily data
  do nfld=1,numflds
   if (nfld==ixsho .or. nfld==ixlng) then
     start = (/1,nrec/)
     call readcdf_1d(ncid1(nfld),varname(nfld),srial,nip,start,kount)

!    if (mod(nrec,5).eq.0) then
      if (nfld==ixlng) then
       do i=1,nip
        data1(i)=srial(perm(i))
       end do
!      call findmxmn1(data1,nip,'longwv',wet)
      else if (nfld==ixsho) then
       do i=1,nip
        data2(i)=srial(perm(i))
       end do
!      call findmxmn1(data2,nip,'shortwv',wet)
      end if	!nfld
!    end if
   end if	! ixsho,ixlng
  end do	! do nfld

!$OMP PARALLEL DO REDUCTION(+:radtot)
     do 44 i=1,nip
     if (wet(i)==1) radtot=radtot+(data1(i)+data2(i))*area(i)
 44  continue
!$OMP END PARALLEL DO

  end do	! do nrec

  nfld=ixpcp
  do nrec=1,12
   start = (/1,nrec/)
   call readcdf_1d(ncid1(nfld),varname(nfld),srial,nip,start,kount)
   do i=1,nip
    data1(i)=srial(perm(i))
   end do
!  call findmxmn1(srial,nip,'precip',wet)

!$OMP PARALLEL DO REDUCTION(+:pcptot)
   do 45 i=1,nip
   if (wet(i)==1) pcptot=pcptot+data1(i)*area(i)	! data1 = mm/sec
45 continue
!$OMP END PARALLEL DO
  end do	! nrec=1,12
  radtot=radtot/(365.*ocnarea)		! hardwired: 1 yr of daily fields
  pcptot=pcptot*1.e-3/(12.*ocnarea)	! hardwired: 1 yr of monthly fields, m/sec

  write (*,'(2(a,f9.1))') 'avg. incoming rad (W/m^2):',radtot,		&
    ' avg. precip (mm/yr):',pcptot*1.e3*365.*86400.
!SMS$SERIAL END

   call findmxmn1(airtmn,nip,'airtmn',wet)
!  call stencl(airtmn,1,1.,'(ocn thermf) airtmn')
!!sms$compare_var(airtmn, "airtmn ")

  return
  end subroutine opnfor

  subroutine rdforf(nstep)
  use netcdf
  use hycom_dffusn,only: dffusn_lev
  use hycom_variables ,only: mth0,mth1,mth2,mth3,day0,day1,day2,day3,	&
			     six0,six1,six2,six3

  integer,intent(IN) :: nstep
  real x,x1,day,dfflen
!SMS$DISTRIBUTE(dh,1) BEGIN
  real data1(nip),data2(nip),data3(nip),data4(nip)
!SMS$DISTRIBUTE END
  real srial(nip)
  integer n,offset,mo,dy,sx,i
  logical vrbos,abort
  character flnm*120,chrglvl*2,errmsg(4)*28
  real*8 totl
  integer :: start(2)	! starting indices for the variable in the netcdf file
  integer :: kount(2)	! lengths for the variable in the netcdf file

  vrbos=mod(nstep-1,diag_intvl)==0
  if (vrbos) then
    write(*,*) 'entering rdforf, nstep =',nstep
  end if

  write (chrglvl,'(a,i1)') 'g',glvl
  day=float(nstep)*dt/86400.
  x=.5+day*12./365.		! indiv.fields are assigned to middle of month
  mo=x

  if (lm0.lt.0) then		! first call

! --- read   a n n u a l  forcing fields
!SMS$SERIAL (<perm,IN>, <data1,OUT> : default=ignore) BEGIN
    flnm=trim(inicondir)//'input_'//chrglvl//'/'//trim(flnm_in(ixrff))//'_core1_'//chrglvl//'.nc'
    write (*,*) 'opening timely forcing field ',trim(flnm)
    call errhandl (nf90_open (path=trim(flnm), mode=NF90_NOWRITE, ncid=ncid1(ixrff)))   ! netcdf3
    start = (/1,1/)
    kount = (/nip,1/)
    call readcdf_1d(ncid1(ixrff),varname(ixrff),srial,nip,start,kount)	! runoff
    write (*,*) 'done with reading ',trim(flnm)
    data1(:)=srial(perm(:))
!SMS$SERIAL END

! --- spread runoff to offshore points to avoid excessive freshwater buildup
    call findmxmn1(data1,nip,'runoff bef smoo',wet)
    call stencl(data1,1,1.e-3,'runoff (mSv) bef smoo')
    totl=0.
!SMS$PARALLEL (dh,i) BEGIN
    do i=1,nip
      if (wet(i) > 0) totl=totl+data1(i)
    end do
!SMS$PARALLEL END
!SMS$REDUCE(totl,SUM)
    print '(a,-6p,f7.3)','global river runoff (Sv) bfore smoothing:',totl

    dfflen=sqrt(5.e13/float(nip))	! SQRT[earth surface/(10*cell count)]
    n=2*(glvl-3)
    print '(a,i3,a)','(rdforf) smooth runoff field',n,' times'
    do i=1,n
      call dffusn_lev(data1,1,1,1,dfflen,.false.)
    end do
    totl=0.
!SMS$PARALLEL (dh,i) BEGIN
    do i=1,nip
      if (wet(i) > 0) totl=totl+data1(i)
    end do
!SMS$PARALLEL END
!SMS$REDUCE(totl,SUM)
    print '(a,-6p,f7.3)','global river runoff (Sv) after smoothing:',totl
!SMS$PARALLEL (dh,i) BEGIN
    do i=1,nip
      if (wet(i) > 0) srforc(1,i,ixrff)=data1(i)/area(i)		! m^3/s -> m/s
    end do
!SMS$PARALLEL END

! --- construct indices for   m o n t h l y   forcing fields

! --- read   m o n t h l y   forcing fields

! --- forcing fields valid at times mo-1,..,mo+2 go into slots lm0,..,lm3

! --- initial staging of forcing function quadruplet
    lm1=mod(mo+3,4)+1
    lm2=mod(lm1,4)+1
    lm3=mod(lm2,4)+1
    lm0=mod(lm3,4)+1

    do nfld=ixpcp,ixsss
      if (nfld==ixpcp) then
        mth0=mod(mo+num_mthly_fields-2,num_mthly_fields)+1
        mth1=mod(mo+num_mthly_fields-1,num_mthly_fields)+1
        mth2=mod(mo+num_mthly_fields  ,num_mthly_fields)+1
        mth3=mod(mo+num_mthly_fields+1,num_mthly_fields)+1
      else if (nfld==ixsss) then
! --- indices for CYCLED monthly climatology field SSS:
        mth0=mod(mo+10,12)+1
        mth1=mod(mo+11,12)+1
        mth2=mod(mo+12,12)+1
        mth3=mod(mo+13,12)+1
      end if
      write (*,100) '0month',mth0,lm0
      write (*,100) '0month',mth1,lm1
      write (*,100) '0month',mth2,lm2
      write (*,100) '0month',mth3,lm3
 100  format (a,i6,':  =======>   slot',i3)
 101  format (a,i6,':  already in slot',i3)
 102  format (a,f10.3,f8.3,i6,a,5i6,' offset')

!SMS$SERIAL (<perm,IN>, <data1,data2,data3,data4,OUT> : default=ignore) BEGIN
      kount = (/nip,1/)
      start = (/1,mth0/)
      call readcdf_1d(ncid1(nfld),varname(nfld),srial,nip,start,kount)
      data1(:)=srial(perm(:))

      start = (/1,mth1/)
      call readcdf_1d(ncid1(nfld),varname(nfld),srial,nip,start,kount)
      data2(:)=srial(perm(:))

      start = (/1,mth2/)
      call readcdf_1d(ncid1(nfld),varname(nfld),srial,nip,start,kount)
      data3(:)=srial(perm(:))

      start = (/1,mth3/)
      call readcdf_1d(ncid1(nfld),varname(nfld),srial,nip,start,kount)
      data4(:)=srial(perm(:))
!SMS$SERIAL END

!SMS$PARALLEL(dh,i) BEGIN
!$OMP PARALLEL DO
      do i=1,nip
        if (nfld==ixpcp) then
! --- precip is given in mm/s. convert to m/s & add river runoff (converted to m/s) to precip
          srforc(lm0,i,nfld)=data1(i)*1.e-3+srforc(1,i,ixrff)
          srforc(lm1,i,nfld)=data2(i)*1.e-3+srforc(1,i,ixrff)
          srforc(lm2,i,nfld)=data3(i)*1.e-3+srforc(1,i,ixrff)
          srforc(lm3,i,nfld)=data4(i)*1.e-3+srforc(1,i,ixrff)
        else if (nfld==ixsss) then
          srforc(lm0,i,nfld)=data1(i)
          srforc(lm1,i,nfld)=data2(i)
          srforc(lm2,i,nfld)=data3(i)
          srforc(lm3,i,nfld)=data4(i)
        end if
      end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END

      if (vrbos) then
        if (nfld==ixpcp) then
          write(*,102)'day,xmon,mo:',day,x,mo,' 0pcp mth0123:',		&
          mth0,mth1,mth2,mth3
        else if (nfld==ixsss) then
          write(*,102)'day,xmon,mo:',day,x,mo,' 0SSS mth0123:',		&
          mth0,mth1,mth2,mth3
        end if
        write(*,*)'0reading monthly field: ',varname(nfld),'=> ',trim(flnm_in(nfld))
      end if
    end do ! nfld=ixpcp,ixsss 

!SMS$PARALLEL(dh,i) BEGIN
      do i=1,nip
        if (vrbos .and. i==itest) write(*,'(2i5,a,2es10.2)') i,lm0,	&
          ' lm0 readf monthly prcp/SSS:',srforc(lm0,i,ixpcp),		&
          srforc(lm0,i,ixsss)
      end do
!SMS$PARALLEL END

  else if (mod(mo+3,4)+1.ne.lm1) then

! --- time marches on: swap out one member of forcing function quadruplet 
    lm1=mod(mo+3,4)+1
    lm2=mod(lm1,4)+1
    lm3=mod(lm2,4)+1
    lm0=mod(lm3,4)+1
    do nfld=ixpcp,ixsss 
      if (nfld==ixpcp) then
        mth0=mod(mo+num_mthly_fields-2,num_mthly_fields)+1
        mth1=mod(mo+num_mthly_fields-1,num_mthly_fields)+1
        mth2=mod(mo+num_mthly_fields  ,num_mthly_fields)+1
        mth3=mod(mo+num_mthly_fields+1,num_mthly_fields)+1
        if (vrbos) then
          write(*,102)'day,xmon,mo:',day,x,mo,' prcp mth0123:',		&
      		mth0,mth1,mth2,mth3
          write(*,'(i9,a,4i6)') int(day),' time interpolation: moprcp',	&
                mth0,mth1,mth2,mth3
        end if
        write (*,101) 'prcp month',mth0,lm0
        write (*,101) 'prcp month',mth1,lm1
        write (*,101) 'prcp month',mth2,lm2
        write (*,100) 'prcp month',mth3,lm3
      elseif (nfld==ixsss) then
        mth0=mod(mo+10,12)+1
        mth1=mod(mo+11,12)+1
        mth2=mod(mo+12,12)+1
        mth3=mod(mo+13,12)+1
        if (vrbos) then
          write(*,102)'day,xmon,mo:',day,x,mo,'  SSS mth0123:',		&
      		mth0,mth1,mth2,mth3
        end if
        write (*,101) ' SSS month',mth0,lm0
        write (*,101) ' SSS month',mth1,lm1
        write (*,101) ' SSS month',mth2,lm2
        write (*,100) ' SSS month',mth3,lm3
      end if
!SMS$SERIAL (<perm,IN>, <data1,OUT> : default=ignore) BEGIN
      kount = (/nip,1/)
      start = (/1,mth3/)
      call readcdf_1d(ncid1(nfld),varname(nfld),srial,nip,start,kount)
      data1(:)=srial(perm(:))
!SMS$SERIAL END

!SMS$PARALLEL(dh,i) BEGIN
      if (nfld==ixpcp) then
! --- precip is given in mm/s. convert to m/s & add river runoff (converted to m/s) to precip
        srforc(lm3,:,nfld)=data1(:)*1.e-3+srforc(1,:,ixrff)
      elseif (nfld==ixsss) then
        srforc(lm3,:,nfld)=data1(:)
      end if
!SMS$PARALLEL END
    end do  ! nfld=ixpcp,ixsss

!SMS$PARALLEL(dh,i) BEGIN
      do i=1,nip
        if (vrbos .and. i==itest) write(*,'(2i5,a,2es10.2)') i,lm3,	&
          ' lm3 readf monthly prcp/SSS:',srforc(lm3,i,ixpcp),		&
          srforc(lm3,i,ixsss)
      end do
!SMS$PARALLEL END

  end if	! lm0<0 or not

  x=mod(x,1.)
  x1=1.-x
! --- quasi-hermite interpolation:
  wm1=x1*(1.+x *(1.-1.5*x ))
  wm2=x *(1.+x1*(1.-1.5*x1))
  wm0=-.5*x *x1*x1
  wm3=-.5*x1*x *x

  if (vrbos) write(*,107) int(day*86400./dt),				&
    ' time interpolation: moncyc',mth0,mth1,mth2,mth3,' slots',		&
    lm0,lm1,lm2,lm3,' weights',wm0,wm1,wm2,wm3

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! --- construct indices for   d a i l y   forcing fields

  x=.5+mod(day,0.+num_daily_fields)	! indiv.fields are assigned to mid-day
  dy=int(mod(x+num_daily_fields-1.,0.+num_daily_fields)+1.)
! --- increment 'offset' each Jan.1 to avoid assigning the same time slot
! --- to days 365 and 1 [a problem created by mod(365,4) = mod(1,4)]
  offset=(day+num_daily_fields-.5)/num_daily_fields
  day0=mod(dy+num_daily_fields-2,num_daily_fields)+1
  day1=mod(dy+num_daily_fields-1,num_daily_fields)+1
  day2=mod(dy+num_daily_fields  ,num_daily_fields)+1
  day3=mod(dy+num_daily_fields+1,num_daily_fields)+1

  if (ld0.lt.0) then		! first call

! --- read   d a i l y   forcing fields

! --- initial staging of forcing function quadruplet
    ld1=mod(dy+3+offset,4)+1
    ld2=mod(ld1,4)+1
    ld3=mod(ld2,4)+1
    ld0=mod(ld3,4)+1
    if (vrbos) write (*,102) 'day,x,dy:',day,x,dy,' 0day0123:',		&
      day0,day1,day2,day3,offset

! --- forcing fields valid at times dy-1,..,dy+2 go into slots ld0,..,ld3

    write (*,100) '0day ',day0,ld0
    write (*,100) '0day ',day1,ld1
    write (*,100) '0day ',day2,ld2
    write (*,100) '0day ',day3,ld3

    do nfld=ixlng,ixsho
      if (vrbos) write (*,*)'reading daily field: ',trim(flnm_in(nfld))
!SMS$SERIAL (<perm,IN>, <data1,data2,data3,data4,OUT> : default=ignore) BEGIN
      kount = (/nip,1/)
      start = (/1,day0/)
      call readcdf_1d(ncid1(nfld),varname(nfld),srial,nip,start,kount)
      data1(:)=srial(perm(:))

      start = (/1,day1/)
      call readcdf_1d(ncid1(nfld),varname(nfld),srial,nip,start,kount)
      data2(:)=srial(perm(:))

      start = (/1,day2/)
      call readcdf_1d(ncid1(nfld),varname(nfld),srial,nip,start,kount)
      data3(:)=srial(perm(:))

      start = (/1,day3/)
      call readcdf_1d(ncid1(nfld),varname(nfld),srial,nip,start,kount)
      data4(:)=srial(perm(:))
!SMS$SERIAL END

!SMS$PARALLEL(dh,i) BEGIN
      srforc(ld0,:,nfld)=data1(:)
      srforc(ld1,:,nfld)=data2(:)
      srforc(ld2,:,nfld)=data3(:)
      srforc(ld3,:,nfld)=data4(:)
!SMS$PARALLEL END

    end do    ! nfld=ixlng,ixsho

!SMS$PARALLEL(dh,i) BEGIN
      do i=1,nip
        if (vrbos .and. i==itest) write (*,'(2i5,a,2es10.2)') i,ld0,	&
          ' ld0 readf daily:   longwv    shortw',			&
          srforc(ld0,i,ixlng),srforc(ld0,i,ixsho)
      end do
!SMS$PARALLEL END

  else if (mod(dy+3+offset,4)+1.ne.ld1) then

! --- time marches on: swap out one member of forcing function quadruplet
    ld1=mod(dy+3+offset,4)+1
    ld2=mod(ld1,4)+1
    ld3=mod(ld2,4)+1
    ld0=mod(ld3,4)+1

    write (*,101) 'day',day0,ld0
    write (*,101) 'day',day1,ld1
    write (*,101) 'day',day2,ld2
    write (*,100) 'day',day3,ld3

    do nfld=ixlng,ixsho
!SMS$SERIAL (<perm,IN>, <data1,OUT> : default=ignore) BEGIN
      kount = (/nip,1/)
      start = (/1,day3/)
      call readcdf_1d(ncid1(nfld),varname(nfld),srial,nip,start,kount)
      data1(:)=srial(perm(:))
!SMS$SERIAL END

!SMS$PARALLEL(dh,i) BEGIN
!JR Could thread this but maybe not enough work
      srforc(ld3,:,nfld)=data1(:)
!SMS$PARALLEL END
    end do    ! nfld=ixlng,ixsho

!SMS$PARALLEL(dh,i) BEGIN
    do i=1,nip
      if (vrbos .and. i==itest) write (*,'(2i5,a,2es10.2)') i,ld3,	&
        ' ld3 readf daily:   longwv    shortw',				&
         srforc(ld3,i,ixlng),srforc(ld3,i,ixsho)
    end do
!SMS$PARALLEL END
  end if 	! ld0<0 or not

  x=mod(x,1.)
  x1=1.-x
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- quasi-hermite interpolation:
!cc      wd1=x1*(1.+x *(1.-1.5*x ))
!cc      wd2=x *(1.+x1*(1.-1.5*x1))
!cc      wd0=-.5*x *x1*x1
!cc      wd3=-.5*x1*x *x
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- linear interpolation, time-smoothed:
  wd0=.25*x1
  wd1=.50*x1+.25*x
  wd2=.25*x1+.50*x
  wd3=       .25*x

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if (vrbos) write(*,107) int(day*86400./dt),' time interpolation:   days',	&
    day0,day1,day2,day3,' slots',ld0,ld1,ld2,ld3,' weights',wd0,wd1,wd2,wd3

! --- construct indices for   6 - h o u r l y   forcing fields

  x=.5+mod(4.*day,num_6hrly_fields+0.)	! fields are assigned to middle of 6-hr intvl
  sx=int(mod(x+num_6hrly_fields-1.,num_6hrly_fields+0.)+1.)
  six0=mod(sx+num_6hrly_fields-2,  num_6hrly_fields)+1
  six1=mod(sx+num_6hrly_fields-1,  num_6hrly_fields)+1
  six2=mod(sx+num_6hrly_fields  ,  num_6hrly_fields)+1
  six3=mod(sx+num_6hrly_fields+1,  num_6hrly_fields)+1

! --- read   6 - h o u r l y   forcing fields

  if (ls0.lt.0) then		! first call

! --- forcing fields valid at times sx-1,..,sx+2 go into slots ls0,..,ls3
! --- initial staging of forcing function quadruplet
    ls1=mod(sx+3,4)+1
    ls2=mod(ls1,4)+1
    ls3=mod(ls2,4)+1
    ls0=mod(ls3,4)+1

    if (vrbos) write (*,102) 'day,x,sx:',day,x,sx,' 0six0123:',six0,six1,six2,six3

    write (*,100) '0hour ',six0,ls0
    write (*,100) '0hour ',six1,ls1
    write (*,100) '0hour ',six2,ls2
    write (*,100) '0hour ',six3,ls3

    do nfld=ixuwi,ixvap
      if(vrbos) write (*,*)'reading field: ',trim(flnm_in(nfld))
!SMS$SERIAL (<perm,IN>, <data1,data2,data3,data4,OUT> : default=ignore) BEGIN
      kount = (/nip,1/)
      start = (/1,six0/)
      call readcdf_1d(ncid1(nfld),varname(nfld),srial,nip,start,kount)
      data1(:)=srial(perm(:))

      start = (/1,six1/)
      call readcdf_1d(ncid1(nfld),varname(nfld),srial,nip,start,kount)
      data2(:)=srial(perm(:))

      start = (/1,six2/)
      call readcdf_1d(ncid1(nfld),varname(nfld),srial,nip,start,kount)
      data3(:)=srial(perm(:))

      start = (/1,six3/)
      call readcdf_1d(ncid1(nfld),varname(nfld),srial,nip,start,kount)
      data4(:)=srial(perm(:))
!SMS$SERIAL END

!SMS$PARALLEL(dh,i) BEGIN
!JR Could thread this but maybe not enough work
        if (nfld==ixtem) then
          srforc(ls0,:,nfld)=data1(:)-KelvC
          srforc(ls1,:,nfld)=data2(:)-KelvC
          srforc(ls2,:,nfld)=data3(:)-KelvC
          srforc(ls3,:,nfld)=data4(:)-KelvC
        else
          srforc(ls0,:,nfld)=data1(:)
          srforc(ls1,:,nfld)=data2(:)
          srforc(ls2,:,nfld)=data3(:)
          srforc(ls3,:,nfld)=data4(:)
        end if
!SMS$PARALLEL END
    end do 	! nfld=ixuwi,ixvap

!SMS$PARALLEL(dh,i) BEGIN
    do i=1,nip
      if (vrbos .and. i==itest) then
        write (*,'(2i5,a,4es10.2)') i,ls0,' ls0 readf 6hrly: u  v  t  q'&
          ,srforc(ls0,i,ixuwi),srforc(ls0,i,ixvwi)			&
          ,srforc(ls0,i,ixtem),srforc(ls0,i,ixvap)
      endif
    end do
!SMS$PARALLEL END

  else if (mod(sx+3,4)+1.ne.ls1) then
! --- time marches on: swap out one member of forcing function quadruplet
    ls1=mod(sx+3,4)+1
    ls2=mod(ls1,4)+1
    ls3=mod(ls2,4)+1
    ls0=mod(ls3,4)+1

    write (*,101) 'hour',six0,ls0
    write (*,101) 'hour',six1,ls1
    write (*,101) 'hour',six2,ls2
    write (*,100) 'hour',six3,ls3

    do nfld=ixuwi,ixvap
!SMS$SERIAL (<perm,IN>, <data1,data2,data3,data4,OUT> : default=ignore) BEGIN
      kount = (/nip,1/)
      start = (/1,six3/)
      call readcdf_1d(ncid1(nfld),varname(nfld),srial,nip,start,kount)
      data1(:)=srial(perm(:))
!SMS$SERIAL END

!SMS$PARALLEL(dh,i) BEGIN
!JR Could thread this but maybe not enough work
        if (nfld==ixtem) then
          srforc(ls3,:,nfld)=data1(:)-KelvC
        else
          srforc(ls3,:,nfld)=data1(:)
        end if
!SMS$PARALLEL END
    end do 	! nfld=ixuwi,ixvap

!SMS$PARALLEL(dh,i) BEGIN
    do i=1,nip
      if (vrbos .and. i==itest) then
        write (*,'(2i5,a,4es10.2)') i,ls3,' ls3 readf 6hrly: u  v  t  q'&
          ,srforc(ls3,i,ixuwi),srforc(ls3,i,ixvwi)			&
          ,srforc(ls3,i,ixtem),srforc(ls3,i,ixvap)
      endif
    end do
!SMS$PARALLEL END
  end if	! ls0<0 or not

  x=mod(x,1.)
  x1=1.-x
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- quasi-hermite interpolation:
!cc      ws1=x1*(1.+x *(1.-1.5*x ))
!cc      ws2=x *(1.+x1*(1.-1.5*x1))
!cc      ws0=-.5*x *x1*x1
!cc      ws3=-.5*x1*x *x
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- linear interpolation, time-smoothed:
  ws0=.25*x1
  ws1=.50*x1+.25*x
  ws2=.25*x1+.50*x
  ws3=       .25*x
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if (vrbos) then
    write (*,102) 'day,x,sx:',day,x,sx,'  six0123:',six0,six1,six2,six3
    write(*,107) int(day*86400./dt),' time interpolation:  6hrly',	&
    six0,six1,six2,six3,' slots',ls0,ls1,ls2,ls3,' weights',		&
    ws0,ws1,ws2,ws3
  end if
!      write (*,105) itest,(trim(flnm_in(n)),n=1,7 )
!      write (*,106) (   srforc(lm0,itest,n),n=1,7 )
!      write (*,105) itest,(trim(flnm_in(n)),n=8,numflds)
!      write (*,106) (   srforc(lm0,itest,n),n=8,numflds)
 105   format (i6,7(4x,a6))
 106   format (10x,7es10.2)
!107   format (i9,a,4i6,a,4i3/9x,a,4f8.3)
 107   format (i9,a,4i6,a,4i2,a,4f8.3)

  abort=.false.
!SMS$PARALLEL(dh,i) BEGIN
  do i=1,nip
    if (wet(i)==1) then
    do n=1,4	! 6hr frq variables
      if (abs(srforc(3,i,n)).gt.100.) then
        write(*,'(a,i9,i2,es10.2)')' wrong srforc(airtem) i=',i,n,srforc(3,i,n)
        errmsg(1)='wrong fields in hyc_6hr'
        abort=.true.
       end if
    end do

    do n=5,6	! daily 2 radiations
      if (abs(srforc(3,i,n)).gt.1000.) then
        write(*,'(a,i9,i2,es10.2)')' wrong srforc(rad) i=',i,n,srforc(3,i,n)
        errmsg(2)='wrong fields in hyc_rad'
        abort=.true.
      end if
    end do

    if (abs(srforc(3,i,7)).gt.1.) then
      write(*,'(a,i9,es10.2)')' wrong srforc(prec) i=',i,srforc(3,i,7)
      errmsg(3)='wrong fields in hyc_precip'
      abort=.true.
    end if

    if (abs(srforc(3,i,8)).gt.50.) then
      write(*,'(a,i9,es10.2)')' wrong srforc(sss) i=',i,srforc(3,i,8)
      errmsg(4)='wrong fields in hyc_sss'
      abort=.true.
    end if

    end if	! wet
  end do	! i
!SMS$PARALLEL END

  if (abort) then
    print '(a28)',errmsg
    stop '(input error)'
  end if

!sms$compare_var(srforc(1,:,ixrff), "runoff")
!sms$compare_var(srforc(1,:,ixtem), "airtemp")
!sms$compare_var(srforc(1,:,ixpcp), "precip")
!sms$compare_var(srforc(1,:,ixsho), "shortw")
!sms$compare_var(srforc(1,:,ixuwi), "u-wind")

!SMS$PARALLEL (dh,i) BEGIN
  do i=1,nip
    if (vrbos .and. i==itest) then 
      write(*,*) 'chk srforc at i=',i,perm(i)
      write(*,'(a)')'       U         V    air_temp   rel.hum     longw    shortw    precip  srf.saln    runoff'
      write(*,'(9es10.2)') ((srforc(mo,i,n),n=1,9),mo=1,4)
    end if
  end do
!SMS$PARALLEL END

  vrbos=mod(nstep-1,diag_intvl)==0
  if (vrbos) print *,'... exiting rdforf'
  return
  end subroutine rdforf

!*********************************************************************
!  thermf_core
!    Compute surface fluxes from CORE data
!    (note: surface fluxes in coupled runs are computed in hycom_run)
!    Shan Sun
!*********************************************************************

  subroutine thermf_core(nstep)  
  use hycom_sigetc    ,only: qsatur
  use hycom_constants ,only: saldif,fusion,stdsal,grvity
  use hycom_variables ,only: heatini,saltini,massglb0
  use hycom_diag      ,only: glob2d,glob3d

! Dimension and type external variables:
  integer,intent (IN) :: nstep		! model time step, starts at 1

  real :: evapw,snsibl,snsibw,snsibi,c_d,c_l,c_s,exchng
  logical   :: vrbos,abort
  real,parameter :: albw = 0.066		! open water albedo
  real,parameter :: albi = 0.80			! ice albedo
  real           :: day,hrs,sunris,sunset,clock,ratio,sqrtdns
  real*8         :: totice,totlt,totls
  real,external :: its2time

#if ( defined NEED_SINDCOSD )
!JR Define required statement functions for situations where
!JR compiler doesn't support them (e.g. gfortran)
!JR use same constant for pi that is used in function sunlit below
   real :: val, sind
   sind(val) = sin(val*3.14159265359/180.)
#endif

! --- --------------------------------
! --- thermal forcing of ocean surface
! --- --------------------------------

  vrbos=mod(nstep-1,diag_intvl)==0
  if (vrbos) print *,'entering thermf_core, nstep =',nstep

  sqrtdns=sqrt(airdns*.001)	! lahey compiler doesn't allow sqrt in initlzn
  call rdforf(nstep)
  day=float(nstep)*dt/86400.

!!sms$compare_var(srforc(lm0,:,ixpcp), "precip0")
!!sms$compare_var(srforc(lm1,:,ixpcp), "precip1")
!!sms$compare_var(srforc(lm2,:,ixpcp), "precip2")
!!sms$compare_var(srforc(lm3,:,ixpcp), "precip3")

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- optional: weak salin./temp. restoring via corrective surface fluxes
  if (nstep.eq.1 .or. mod(day+.0001,30.).lt..0002) then		!  once a month
! if (nstep.eq.1 .or. mod(day+.0001,1.).lt..0002) then		!  once a day
   call glob2d(thkice,totice)
   totice=totice*rhoice					! kg
   call glob3d(temp,dp,totlt)
   totlt=totlt*spcifh/grvity-totice*fusion		! J
   call glob3d(saln,dp,totls)
   totls=totls/grvity-totice*saldif			! g
! --- initialize total heat/salt content
   if (heatini.le.0.) heatini=totlt
   if (saltini.le.0.) saltini=totls
  end if

  if (nstep.gt.1 .and. mod(day+.0001,30.).lt..0002) then        !  once a month
!ss pcpcor/radcor calculation is moved to hycom_diag.F90
!ss! if (nstep.gt.1 .and. mod(day+.0001,1.).lt..0002) then	!  once a day
!ss   pcpcor=(totls-saltini)/(ocnarea*7300.*86400.)	&	!  g/m^2/sec
!ss     *thref/stdsal					&	!  => m/sec
!ss     /pcptot							!  => rel.units
!ss   write (*,101) nstep, '  overall salt gain (psu):',	&
!ss      (totls-saltini)/massglb0,'  => precip''tn factor',(1.d0+pcpcor)
!ss!  radcor=(heatini-totlt)/(ocnarea*7300.*86400.)	&	!  W/m^2
!ss!    /radtot							!  => rel.units
!ss   radcor=0.
!ss   write (*,101) nstep, '  overall temp gain (deg):',	&
!ss     (totlt-heatini)/(spcifh*massglb0),'  => radiation factor',(1.d0+radcor)
!ss 101 format (i9,2(a,f9.5))
   write (*,'(i9,1x,a,f8.2,a,2f7.3)') nstep,ArchvTimeUnit,its2time(nstep),	&
	'  hycom_thermf adjustment factor=',1.+radcor,1.+pcpcor
  end if
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  abort=.false.
!SMS$PARALLEL(dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (vrbos,hrs,sunris,sunset,clock,ratio,		&
!$OMP c_d,c_l,c_s,exchng,evapw,snsibw,snsibi)
  do i=1,nip
   if (wet(i) == 1) then
    vrbos=i.eq.itest .and. mod(nstep-1,diag_intvl)==0

    prcp(i) = srforc(lm0,i,ixpcp)*wm0+srforc(lm1,i,ixpcp)*wm1		&
             +srforc(lm2,i,ixpcp)*wm2+srforc(lm3,i,ixpcp)*wm3
    salnow(i)=srforc(lm0,i,ixsss)*wm0+srforc(lm1,i,ixsss)*wm1		&
             +srforc(lm2,i,ixsss)*wm2+srforc(lm3,i,ixsss)*wm3
    rsds(i) = srforc(ld0,i,ixsho)*wd0+srforc(ld1,i,ixsho)*wd1		&
             +srforc(ld2,i,ixsho)*wd2+srforc(ld3,i,ixsho)*wd3
    rlds(i) = srforc(ld0,i,ixlng)*wd0+srforc(ld1,i,ixlng)*wd1		&
             +srforc(ld2,i,ixlng)*wd2+srforc(ld3,i,ixlng)*wd3
    uwnd(i) = srforc(ls0,i,ixuwi)*ws0+srforc(ls1,i,ixuwi)*ws1		&
             +srforc(ls2,i,ixuwi)*ws2+srforc(ls3,i,ixuwi)*ws3
    vwnd(i) = srforc(ls0,i,ixvwi)*ws0+srforc(ls1,i,ixvwi)*ws1		&
             +srforc(ls2,i,ixvwi)*ws2+srforc(ls3,i,ixvwi)*ws3
    airt(i) = srforc(ls0,i,ixtem)*ws0+srforc(ls1,i,ixtem)*ws1		&
             +srforc(ls2,i,ixtem)*ws2+srforc(ls3,i,ixtem)*ws3
    vpmx(i) = srforc(ls0,i,ixvap)*ws0+srforc(ls1,i,ixvap)*ws1		&
             +srforc(ls2,i,ixvap)*ws2+srforc(ls3,i,ixvap)*ws3
! -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
! --- optional: superimpose day-night cycle on short-wave radiation:
    hrs=sunlit(mod(day+270.,360.),deg_lat(i))
    hrs=max(4.,min(hrs,20.))
! --- not strictly necessary, but a nice touch: use local time
    sunris=24.-deg_lon(i)/15.-.5*hrs
    sunset=24.-deg_lon(i)/15.+.5*hrs
    clock=mod(day,1.)*24.
    if ((clock    .gt.sunris.and.clock    .lt.sunset)  &
      .or.(clock-24..gt.sunris.and.clock-24..lt.sunset)  &
      .or.(clock+24..gt.sunris.and.clock+24..lt.sunset)) then
      ratio=24./hrs
    else
      ratio=0.
    end if
! --- assume sinusoidal variation of insolation during daylight hours
    clock=mod(clock-sunris+24.,24.)
    rsds(i)=rsds(i)*ratio*sind(180.*clock/hrs)*1.5707963268
! -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
    if (do_pcpcor) prcp(i)=prcp(i)*(1.d0+pcpcor)
    if (do_radcor) rsds(i)=rsds(i)*(1.d0+radcor)

    rsus(i)=-rsds(i)*(covice(i)*albi+(1.-covice(i))*albw)	! positive down
    rlus(i)= - 5.67e-8*(temp(1,i)+KelvC)**4*(1.-.05*covice(i))

! --- drag coefficient:
! --- from Kara,Rochford,Hurlburt, B.L.Met., 103, 439-458 (2002)
! --- paper (2005)
    wspd(i)=max(.5,sqrt(uwnd(i)**2+vwnd(i)**2))
    call thermf_cl(wspd(i),vpmx(i),airt(i),temp(1,i),c_l)
    c_d = cd_coare(wspd(i),vpmx(i),airt(i),temp(1,i))
    c_d=max(0.,c_d)
    c_l=max(0.,c_l)
    c_s=.95*c_l
    exchng=airdns*wspd(i)

! --- get oceanic ustar and ekman depth (needed in mixed layer)
    ustar(i)=sqrt(thref*c_d*airdns)*wspd(i)*(1.-.9*covice(i))
    us2d(i)=ustar(i)/sqrtdns
!   hekman(i,j)=ustar(i,j)*(cekman*4.0)/			&
!              (abs(corio(i,j  ))+abs(corio(i+1,j  ))+		&
!               abs(corio(i,jb ))+abs(corio(i+1,jb )))

! --- evapw = latent heat flux over water (W/m^2) in +p direction
! --- snsibw,snsibi = sensible heat flux over water and ice (W/m^2) in +p direc
! --- (extremely low airt over open water is unrealistic -> set evapw=0)

    evapw=c_l*exchng*evaplh*(vpmx(i)-qsatur(temp(1,i)))	! pos.down W/m^2
    snsibw=c_s*exchng*csubp*(airt(i)-temp(1,i))		! pos.down W/m^2
    snsibi=c_s*exchng*csubp*(airt(i)-temice(i))			! pos.dow, W/m^2

    hf2d(i) = snsibw*(1.-covice(i))+snsibi*covice(i)
    qf2d(i) = evapw *(1.-covice(i))
    srfx(i)=rsds(i)+rsus(i)+rlus(i)+rlds(i)+hf2d(i)+qf2d(i)	! pos.down W/m^2
    taux(i)=airdns*c_d*wspd(i)*uwnd(i)*(1.-.9*covice(i))
    tauy(i)=airdns*c_d*wspd(i)*vwnd(i)*(1.-.9*covice(i))

! --- p-minus-e (m/s): prcp (m/s), qf2d (W/m^2)
    pmne(i)=prcp(i) + qf2d(i)*thref/evaplh			! pos.down m/s

    if (vrbos) then
     write (*,100) '(thermf)  i=',perm(i),'  lat/lon=',			&
      deg_lat(i),deg_lon(i),'  day',day,				&
      'heatflx',srfx(i),'shrtw_dn',rsds(i),'longw_dn',rlds(i),		&
      'snsibl',hf2d(i),'latent',qf2d(i),'PminusE',pmne(i),		&
      'SST',temp(1,i),'airtem',airt(i),'vapres',vpmx(i),		&
      'satvapr',qsatur(temp(1,i)),'senswatr',snsibw,			&
      'sensice',snsibi,'dragcf',c_d,'uwind',uwnd(i),			&
      'vwind',vwnd(i),'taux',taux(i),'tauy',tauy(i),'covi',covice(i)
 100 format(a,i8,a,2f8.2,a,f8.2/(5(a9,"=",es10.3)))
    end if

    if (max(abs(temp(1,i)),abs(airt(i))).gt.100.) then
      write(*,'(a,2i8,2es10.2)')  'wrong temp/airt =',nstep,perm(i),	&
       temp(1,i),airt(i)
      abort=.true.
    end if

   end if			! ocean point
  end do			! horizontal loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END

  if (abort) then
   stop '(error)'
  end if

  if (nstep.eq.1.or.mod(nstep,diag_intvl).eq.0) then
!  call findmxmn1(srfx,nip,'srfx',wet)
!  call findmxmn1(rsds,nip,'rsds',wet)
!  call findmxmn1(rlds,nip,'rlds',wet)
!  call findmxmn1(hf2d,nip,'hf2d',wet)
!  call findmxmn1(qf2d,nip,'qf2d',wet)
!  call findmxmn1(pmne,nip,'pmne',wet)
!  call findmxmn1(us2d,nip,'ustar',wet)

   call findmxmn1(srfx,nip,'srfx')
   call findmxmn1(rsds,nip,'rsds')
   call findmxmn1(rlds,nip,'rlds')
   call findmxmn1(hf2d,nip,'hf2d')
   call findmxmn1(qf2d,nip,'qf2d')
   call findmxmn1(pmne,nip,'pmne')
   call findmxmn1(us2d,nip,'ustar')
  end if

  vrbos=mod(nstep-1,diag_intvl)==0

  if (vrbos) print *,'... exiting thermf_core'
  return
  end subroutine thermf_core

  real function sunlit(day_of_year,latitude)

! --- daylight hours as a function of latitude and time of year.

! --- approach: find angle between the unit vector pointing toward the sun
! --- (a function of time of year) and the unit vector marking the position
! --- on earth (a function of latitude and time of day). sunrise and sunset
! --- occur when this angle is zero.

! --- the unit vector pointing at the sun is S = [sqrt(1-q*q) ; 0 ; q]
! --- where q=sin(23.5)*sin(day) (1yr=360days). the position on earth is
! --- P = [cos(long)*cos(lat) ; sin(long)*cos(lat) ; sin(lat)]
! --- where long = longitude relative to noontime sun.

! --- setting lambda=0 yields P*S = cos(noontime zenith angle).
! --- setting P*S = 0 yields an expression for the sunrise/sunset times.

  real,intent(IN) :: day_of_year	! days from spring
  real,intent(IN) :: latitude		! degrees
  real q,cos_riseset,sind,tand,acosd
!JR Should use ifdef NEED_SINDCOSD instead here to be consistent
  real,parameter :: ang2rad=3.14159265359/180.
  real,parameter :: rad2ang=180./3.14159265359
  sind(q)=sin(q*ang2rad)
  tand(q)=tan(q*ang2rad)
  acosd(q)=rad2ang*acos(q)

  q=sind(23.5)*sind(day_of_year)
  cos_riseset=-tand(latitude)*q/sqrt(1.-q*q)
  if (cos_riseset.gt.1.) then
      sunlit=0.				!  24-hr night
  else if (cos_riseset.lt.-1.) then
      sunlit=24.			!  24-hr day
  else
      sunlit=acosd(cos_riseset)/7.5
  end if
  return
  end function sunlit

  real function cd_coare(va_in,vpmx0,airt0,sst)
  use hycom_sigetc, only: qsatur

! --- va_in = wind speed (m/s)
! --- vpmx0  = water vapor mixing ratio (kg/kg)
! --- airt0  = air temperature (C)
! --- sst   = sea surface temperature (C)

  real va_in,vpmx0,airt0,sst

! ---              Ta-Ts
! ---           ==============
! ---   STABLE:  7    to  0.75 degC
! ---  NEUTRAL:  0.75 to -0.75 degC
! --- UNSTABLE: -0.75 to -8    degC

! ---              Va
! ---           ==============
! ---   Low:     1    to   5   m/s
! ---   High:    5    to  40   m/s

      real    tamts,q,qva,va

      real, parameter :: vamin=  1.0,  vamax=40.0
      real, parameter :: tdmin= -8.0,  tdmax= 7.0
      real, parameter :: tzero=KelvC

      real, parameter :: 					&
        as0_00=-0.06695,   as0_10= 0.09966,  as0_20=-0.02477,	&
        as0_01= 0.3133,    as0_11=-2.116,    as0_21= 0.2726,	&
        as0_02=-0.001473,  as0_12= 4.626,    as0_22=-0.5558, 	&
        as0_03=-0.004056,  as0_13=-2.680,    as0_23= 0.3139 

      real, parameter ::					&
        as5_00= 0.55815,   as5_10=-0.005593, as5_20= 0.0006024,	&
        as5_01= 0.08174,   as5_11= 0.2096,   as5_21=-0.02629,	&
        as5_02=-0.0004472, as5_12=-8.634,    as5_22= 0.2121,	&
        as5_03= 2.666e-6,  as5_13= 18.63,    as5_23= 0.7755

      real, parameter ::					&
        au0_00= 1.891,     au0_10=-0.006304, au0_20= 0.0004406,	&
        au0_01=-0.7182,    au0_11=-0.3028,   au0_21=-0.01769,	&
        au0_02= 0.1975,    au0_12= 0.3120,   au0_22= 0.01303,	&
        au0_03=-0.01790,   au0_13=-0.1210,   au0_23=-0.003394

      real, parameter ::					&
        au5_00= 0.6497,    au5_10= 0.003827, au5_20=-4.83e-5,	&
        au5_01= 0.06993,   au5_11=-0.2756,   au5_21= 0.007710,	&
        au5_02= 3.541e-5,  au5_12=-1.091,    au5_22=-0.2555,	&
        au5_03=-3.428e-6,  au5_13= 4.946,    au5_23= 0.7654

      real, parameter ::					&
        an0_00= 1.057,     an5_00= 0.6825,			&
        an0_01=-0.06949,   an5_01= 0.06945,			&
        an0_02= 0.01271,   an5_02=-0.0001029

      real, parameter ::					&
        ap0_10= as0_00 + as0_10*0.75 + as0_20*0.75**2,		&
        ap0_11=          as0_11*0.75 + as0_21*0.75**2,		&
        ap0_12=          as0_12*0.75 + as0_22*0.75**2,		&
        ap0_13=          as0_13*0.75 + as0_23*0.75**2

      real, parameter ::					&
        ap5_10= as5_00 + as5_10*0.75 + as5_20*0.75**2,		&
        ap5_11=          as5_11*0.75 + as5_21*0.75**2,		&
        ap5_12=          as5_12*0.75 + as5_22*0.75**2,		&
        ap5_13=          as5_13*0.75 + as5_23*0.75**2

      real, parameter ::					&
        am0_10= au0_00 - au0_10*0.75 + au0_20*0.75**2,		&
        am0_11=        - au0_11*0.75 + au0_21*0.75**2,		&
        am0_12=        - au0_12*0.75 + au0_22*0.75**2,		&
        am0_13=        - au0_13*0.75 + au0_23*0.75**2

      real, parameter ::					&
        am5_10= au5_00 - au5_10*0.75 + au5_20*0.75**2,		&
        am5_11=        - au5_11*0.75 + au5_21*0.75**2,		&
        am5_12=        - au5_12*0.75 + au5_22*0.75**2,		&
        am5_13=        - au5_13*0.75 + au5_23*0.75**2

! --- saturation specific humidity (lowe, j.appl.met., 16, 100-103, 1976)

          tamts = airt0-sst - 0.61*(airt0+tzero)*(qsatur(airt0)-vpmx0)
          tamts = min( tdmax, max( tdmin, tamts ) )
          va    = max(vamin,min(vamax,va_in))
          qva   = 1.0/va
          if     (va.le.5.0) then
            if     (tamts.ge. 0.75) then
              cd_coare = 						&
                 (as0_00 + as0_01* va + as0_02* va**2 + as0_03* va**3)	&
               + (as0_10 + as0_11*qva + as0_12*qva**2 + as0_13*qva**3)	&
                 *tamts							&
               + (as0_20 + as0_21*qva + as0_22*qva**2 + as0_23*qva**3)	&
                 *tamts**2 
            elseif (tamts.le.-0.75) then
              cd_coare = 						&
                 (au0_00 + au0_01* va + au0_02* va**2 + au0_03* va**3)	&
               + (au0_10 + au0_11*qva + au0_12*qva**2 + au0_13*qva**3)	&
                 *tamts							&
               + (au0_20 + au0_21*qva + au0_22*qva**2 + au0_23*qva**3)	&
                 *tamts**2 
            elseif (tamts.ge. -0.098)  then
              q =  (tamts+0.098)/0.848  !linear between  0.75 and -0.098
              cd_coare = q*						&
              (  (         as0_01* va + as0_02* va**2 + as0_03* va**3)	&
               + (ap0_10 + ap0_11*qva + ap0_12*qva**2 + ap0_13*qva**3)	&
              ) + (1.0-q)*						&
                 (an0_00 + an0_01* va + an0_02* va**2)
            else
              q = (-tamts-0.098)/0.652  !linear between -0.75 and -0.098
              cd_coare = q*						&
              (  (         au0_01* va + au0_02* va**2 + au0_03* va**3)	&
               + (am0_10 + am0_11*qva + am0_12*qva**2 + am0_13*qva**3)	&
              ) + (1.0-q)*						&
                 (an0_00 + an0_01* va + an0_02* va**2)
            endif !tamts
          else !va>5
            if     (tamts.ge. 0.75) then
              cd_coare = 						&
                 (as5_00 + as5_01* va + as5_02* va**2 + as5_03* va**3)	&
               + (as5_10 + as5_11*qva + as5_12*qva**2 + as5_13*qva**3)	&
                 *tamts							&
               + (as5_20 + as5_21*qva + as5_22*qva**2 + as5_23*qva**3)	&
                 *tamts**2 
            elseif (tamts.le.-0.75) then
              cd_coare = 						&
                 (au5_00 + au5_01* va + au5_02* va**2 + au5_03* va**3)	&
               + (au5_10 + au5_11*qva + au5_12*qva**2 + au5_13*qva**3)	&
                 *tamts							&
               + (au5_20 + au5_21*qva + au5_22*qva**2 + au5_23*qva**3)	&
                 *tamts**2 
            elseif (tamts.ge. -0.098)  then
              q =  (tamts+0.098)/0.848  !linear between  0.75 and -0.098
              cd_coare = q*						&
              (  (         as5_01* va + as5_02* va**2 + as5_03* va**3)	&
               + (ap5_10 + ap5_11*qva + ap5_12*qva**2 + ap5_13*qva**3)	&
              ) + (1.0-q)*						&
                 (an5_00 + an5_01* va + an5_02* va**2)
            else
              q = (-tamts-0.098)/0.652  !linear between -0.75 and -0.098
              cd_coare = q*						&
              (  (         au5_01* va + au5_02* va**2 + au5_03* va**3)	&
               + (am5_10 + am5_11*qva + am5_12*qva**2 + am5_13*qva**3)	&
              ) + (1.0-q)*						&
                 (an5_00 + an5_01* va + an5_02* va**2)
            endif !tamts
          endif !va

   	cd_coare=cd_coare*1.e-3		!!! added by RB
   return
   end function cd_coare

      subroutine thermf_cl(wind,vpmx0,airt0,sst, clh)
      use hycom_sigetc, only: qsatur

! --- wind = wind speed (m/s)
! --- vpmx0 = water vapor mixing ratio (kg/kg)
! --- airt0 = air temperature (C)
! --- sst  = sea temperature (C)

      real,intent(IN)  :: wind,vpmx0,airt0,sst
      real,intent(OUT) :: clh

! --- clh    = coefficent of latent heat exchange
! --- qs     = saturation specific humidity at sst
! --- csh    = coefficent of sensible heat exchange
! --- evap   = evaporation         (w/m**2) into ocean
! --- snsibl = sensible heat flux  (w/m**2) into ocean

! --- Latent and Sensible heat fluxes from an approximation
! --- to the COARE 3.0 bulk algorithm (Fairall et al. 2003).

      real    csh,rair,slat,ssen,tdif,tamts,q,qva,va

! --- 'csubp'  = specific heat of air at constant pressure (j/kg/deg)
! --- 'evaplh' = latent heat of evaporation (j/kg)
! --- 'pairc'  = air pressure (mb) * 100
! --- 'rgas'   = gas constant (j/kg/k)
! --- 'tzero'  = celsius to kelvin temperature offset

      real       csubplcl,evaplhlcl
      real       pairc,rgas,tzero
      parameter (csubplcl =1005.7,				&
                 evaplhlcl=2.47e6)
      parameter (pairc=1013.0*100.0,				&
                 rgas =287.1,   tzero=KelvC)

! --- 'vamin'  = minimum allowed wind speed (for cl)
! --- 'vamax'  = maximum allowed wind speed (for cl)
! --- 'tdmin'  = minimum allowed Ta-Ts      (for cl)
! --- 'tdmax'  = maximum allowed Ta-Ts      (for cl)

! --- 'as0_XX' =  stable Ta-Ts  polynominal coefficients, va<=5m/s
! --- 'as5_XX' =  stable Ta-Ts  polynominal coefficients, va>=5m/s
! --- 'au0_XX' = unstable Ta-Ts polynominal coefficients, va<=5m/s
! --- 'au5_XX' = unstable Ta-Ts polynominal coefficients, va>=5m/s
! --- 'an0_XX' =  neutral Ta-Ts polynominal coefficients, va<=5m/s
! --- 'an5_XX' =  neutral Ta-Ts polynominal coefficients, va>=5m/s
! --- 'ap0_XX' =    +0.75 Ta-Ts polynominal coefficients, va<=5m/s
! --- 'ap5_XX' =    +0.75 Ta-Ts polynominal coefficients, va>=5m/s
! --- 'am0_XX' =    -0.75 Ta-Ts polynominal coefficients, va<=5m/s
! --- 'am5_XX' =    -0.75 Ta-Ts polynominal coefficients, va>=5m/s

      real, parameter :: vamin= 1.2, vamax=40.0
      real, parameter :: tdmin=-8.0, tdmax= 7.0

      real, parameter ::						&
        as0_00=-2.925e-4,   as0_10= 7.272e-5,  as0_20=-6.948e-6,	&
        as0_01= 5.498e-4,   as0_11=-1.740e-4,  as0_21= 1.637e-5,	&
        as0_02=-5.544e-5,   as0_12= 2.489e-5,  as0_22=-2.618e-6
      real, parameter ::						&
        as5_00= 1.023e-3,   as5_10=-2.672e-6,  as5_20= 1.546e-6,	&
        as5_01= 9.657e-6,   as5_11= 2.103e-4,  as5_21=-6.228e-5,	&
        as5_02=-2.281e-8,   as5_12=-5.329e-3,  as5_22= 5.094e-4
      real, parameter ::						&
        au0_00= 2.077e-3,   au0_10=-2.899e-4,  au0_20=-1.954e-5,	&
        au0_01=-3.933e-4,   au0_11= 7.350e-5,  au0_21= 5.483e-6,	&
        au0_02= 3.971e-5,   au0_12=-6.267e-6,  au0_22=-4.867e-7
      real, parameter ::						&
        au5_00= 1.074e-3,   au5_10= 6.912e-6,  au5_20= 1.849e-7,	&
        au5_01= 5.579e-6,   au5_11=-2.244e-4,  au5_21=-2.167e-6,	&
        au5_02= 5.263e-8,   au5_12=-1.027e-3,  au5_22=-1.010e-4
      real, parameter ::						&
        an0_00= 1.14086e-3, an5_00= 1.073e-3,				&
        an0_01=-3.120e-6,   an5_01= 5.531e-6,				&
        an0_02=-9.300e-7,   an5_02= 5.433e-8
      real, parameter ::						&
        ap0_00= as0_00 + as0_10*0.75 + as0_20*0.75**2,			&
        ap0_01= as0_01 + as0_11*0.75 + as0_21*0.75**2,			&
        ap0_02= as0_02 + as0_12*0.75 + as0_22*0.75**2
      real, parameter ::						&
        ap5_00= as5_00 + as5_10*0.75 + as5_20*0.75**2,			&
        ap5_01= as5_01,							&
        ap5_02= as5_02,							&
        ap5_11=          as5_11*0.75 + as5_21*0.75**2,			&
        ap5_12=          as5_12*0.75 + as5_22*0.75**2
      real, parameter ::						&
        am0_00= au0_00 - au0_10*0.75 + au0_20*0.75**2,			&
        am0_01= au0_01 - au0_11*0.75 + au0_21*0.75**2,			&
        am0_02= au0_02 - au0_12*0.75 + au0_22*0.75**2
      real, parameter ::						&
        am5_00= au5_00 - au5_10*0.75 + au5_20*0.75**2,			&
        am5_01= au5_01,							&
        am5_02= au5_02,							&
        am5_11=        - au5_11*0.75 + au5_21*0.75**2,			&
        am5_12=        - au5_12*0.75 + au5_22*0.75**2

! --- saturation specific humidity (lowe, j.appl.met., 16, 100-103, 1976)
! ---   Similar to flxflg.eq.2, but with Cl based on an approximation
! ---   to values from the COARE 3.0 algorithm (Fairall et al. 2003), 
! ---   for Cl over the global ocean in the range 1m/s <= Va <= 40m/s
! ---   and -8degC <= Ta-Ts <= 7degC, that is quadratic in Ta-Ts and
! ---   quadratic in either Va or 1/Va (Kara and Hurlburt, 2004).

! ---   Fairall, C. W., E. F. Bradley, J. E. Hare, A. A. Grachev, and J. B.
! ---   Edson, 2003:  Bulk parameterization of air-sea fluxes:  Updates 
! ---   and verification for the COARE algorithm.  J. Climate, 16, 571-591.

! ---   Kara, K. Birol, H.E. Hurlburt, 2004: A Note on the Stability-Dependent
! ---   Exchange Coefficients of Air-Sea Fluxes for Use in General Circulation
! ---   Models. submitted to J. Atmos. Oceanic Technol. on January 14, 2004. 
! ---   preprint at http://hycom.rsmas.miami.edu/publications.html

        rair = pairc / (rgas * ( tzero + airt0 ))
        slat = evaplhlcl*rair
        ssen = csubplcl *rair

        tdif  = sst - airt0
        tamts = -tdif - 0.61*(airt0+tzero)*(qsatur(airt0)-vpmx0)
        tamts = min( tdmax, max( tdmin, tamts ) )
        va    = min( vamax, max( vamin,  wind ) )
        if     (va.le.5.0) then
          if     (tamts.gt. 0.75) then !stable
            clh =   (as0_00 + as0_01* va + as0_02* va**2)		&
                  + (as0_10 + as0_11* va + as0_12* va**2)*tamts		&
                  + (as0_20 + as0_21* va + as0_22* va**2)*tamts**2
          elseif (tamts.lt.-0.75) then !unstable
            clh =   (au0_00 + au0_01* va + au0_02* va**2)		&
                  + (au0_10 + au0_11* va + au0_12* va**2)*tamts		&
                  + (au0_20 + au0_21* va + au0_22* va**2)*tamts**2
          elseif (tamts.ge.-0.098)  then
            q = (tamts-0.75)/0.848  !linear between  0.75 and -0.098
            q = q**2  !favor  0.75
            clh = (1.0-q)*(ap0_00 + ap0_01* va + ap0_02* va**2)		&
                     + q *(an0_00 + an0_01* va + an0_02* va**2)
          else
            q = (tamts+0.75)/0.652  !linear between -0.75 and -0.098
            q = q**2  !favor -0.75
            clh = (1.0-q)*(am0_00 + am0_01* va + am0_02* va**2)		&
                     + q *(an0_00 + an0_01* va + an0_02* va**2)
          endif !tamts
        else !va>5
          qva = 1.0/va
          if     (tamts.gt. 0.75) then !stable
            clh =   (as5_00 + as5_01* va + as5_02* va**2)		&
                  + (as5_10 + as5_11*qva + as5_12*qva**2)*tamts		&
                  + (as5_20 + as5_21*qva + as5_22*qva**2)*tamts**2
          elseif (tamts.lt.-0.75) then !unstable
            clh =   (au5_00 + au5_01* va + au5_02* va**2)		&
                  + (au5_10 + au5_11*qva + au5_12*qva**2)*tamts		&
                  + (au5_20 + au5_21*qva + au5_22*qva**2)*tamts**2
          elseif (tamts.ge.-0.098)  then
            q = (tamts-0.75)/0.848  !linear between  0.75 and -0.098
            q = q**2  !favor  0.75
            clh = (1.0-q)*(ap5_00 + ap5_01* va + ap5_02* va**2		&
                                  + ap5_11*qva + ap5_12*qva**2)		&
                     + q *(an5_00 + an5_01* va + an5_02* va**2)
          else
            q = (tamts+0.75)/0.652  !linear between -0.75 and -0.098
            q = q**2  !favor -0.75
            clh = (1.0-q)*(am5_00 + am5_01* va + am5_02* va**2		&
                                  + am5_11*qva + am5_12*qva**2)		&
                     + q *(an5_00 + an5_01* va + an5_02* va**2)
          endif !tamts
        endif !va
!       csh  = 0.9554*clh

! ---   qs     = saturation specific humidity at sst
! ---   evap   = evaporation         (w/m**2) into ocean
! ---   snsibl = sensible heat flux  (w/m**2) into ocean
!       qs     = qsatur(sst)
!       evap   = -slat*clh*wind*(qs-vpmx0)
!       snsibl = -ssen*csh*wind* tdif
        return
        end subroutine thermf_cl


  subroutine errhandl (ret)
  use netcdf

  integer, intent(in) :: ret

  if (ret /= NF90_NOERR) then
    write(6,*) nf90_strerror (ret)
    stop 999
  end if
  return
  end subroutine errhandl

   subroutine readcdf_1d(ncid,varname,field,nip,start,kount)

! --- extract a 2-D field from a previously opened netcdf file

   use netcdf
   integer      ,intent(IN)   :: ncid,nip,start(2),kount(2)
   character*(*),intent(IN)   :: varname
   real       ,intent(OUT)    :: field(nip)
   integer :: tz,id_fld

   call errhandl (nf90_inq_varid (ncid, trim(varname), id_fld))
   call errhandl (nf90_get_var (ncid, id_fld, field, start, kount))

   return
   end subroutine readcdf_1d

      subroutine readcdf_2d(ncid,varname,field,idm,jdm,start,count,check)
!
! --- extract a 2-D field (real) from a previously opened netcdf file
!
      use netcdf
      implicit none
      integer        ,intent(IN)      :: ncid,idm,jdm,start(3),count(3)
      character*(*),intent(IN)      :: varname
      logical,optional,intent(IN)  :: check
      real           ,intent(OUT)     :: field(idm,jdm)
      integer :: iz,jz,tz,id_idm,id_jdm,id_tdm,id_fld

      if (present(check)) then
       if (check) then
! --- verify that dimensions of archived field agree with specified dimensions
        call errhandl (nf90_inq_dimid (ncid, 'idm', id_idm))
        call errhandl (nf90_inq_dimid (ncid, 'jdm', id_jdm))
        call errhandl (nf90_inq_dimid (ncid, 'time',id_tdm))

        call errhandl (nf90_inquire_dimension (ncid, id_idm, len=iz))
        call errhandl (nf90_inquire_dimension (ncid, id_jdm, len=jz))
        call errhandl (nf90_inquire_dimension (ncid, id_tdm, len=tz))

        if (iz.ne.idm .or. jz.ne.jdm) then
          print '(2(a,2i5))','(readcdf) array dimensions are',iz,jz,          &
            ', not the epected',idm,jdm
          stop '(readcdf error)'
        end if

        if (start(3)+count(3)-1.gt.tz) then
          print '(2(a,i5))','(readcdf) requested record',                         &
            start(3)+count(3)-1,'  does not exist. records in file:',tz
          stop '(readcdf error)'
        end if
       end if
      end if

      call errhandl (nf90_inq_varid (ncid, trim(varname), id_fld))
      call errhandl (nf90_get_var (ncid, id_fld, field, start=start, count=count))
      return
      end subroutine readcdf_2d


      subroutine readcdf_3d(ncid,varname,field,idm,jdm,kdm,start,count,check)
!
! --- extract a 3-D field (real) from a previously opened netcdf file
!
      use netcdf
      implicit none
      integer        ,intent(IN)      :: ncid,idm,jdm,kdm,start(4),count(4)
      character*(*),intent(IN)      :: varname
      logical,optional,intent(IN)  :: check
      real           ,intent(OUT)     :: field(idm,jdm,kdm)
      integer :: iz,jz,kz,tz,id_idm,id_jdm,id_kdm,id_tdm,id_fld
   
      if (present(check)) then
       if (check) then
! --- verify that dimensions of archived field agree with specified dimensions
        call errhandl (nf90_inq_dimid (ncid, 'idm', id_idm))
        call errhandl (nf90_inq_dimid (ncid, 'jdm', id_jdm))
        call errhandl (nf90_inq_dimid (ncid, 'kdm', id_kdm))
        call errhandl (nf90_inq_dimid (ncid, 'time',id_tdm))
   
        call errhandl (nf90_inquire_dimension (ncid, id_idm, len=iz))
        call errhandl (nf90_inquire_dimension (ncid, id_jdm, len=jz))
        call errhandl (nf90_inquire_dimension (ncid, id_kdm, len=kz))
        call errhandl (nf90_inquire_dimension (ncid, id_tdm, len=tz))

        if (iz.ne.idm .or. jz.ne.jdm .or. kz.ne.kdm) then
          print '(2(a,3i5))','(readcdf) array dimensions are',iz,jz,kz,          &
            ', not the epected',idm,jdm,kdm
          stop '(readcdf error)'
        end if

        if (start(4)+count(4)-1.gt.tz) then
          print '(2(a,i5))','(readcdf) requested record',                         &
            start(4)+count(4)-1,'  does not exist. records in file:',tz
          stop '(readcdf error)'
        end if
       end if
      end if

      call errhandl (nf90_inq_varid (ncid, trim(varname), id_fld))
      call errhandl (nf90_get_var (ncid, id_fld, field, start=start, count=count))
      return
      end subroutine readcdf_3d

end module hycom_thermf
