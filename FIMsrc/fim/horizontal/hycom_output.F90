module hycom_output
use findmaxmin1
use findmaxmin2
use findmaxmin3
contains
!*********************************************************************
!       hycom_output
!       Output program for ocean model
!       R.Bleck                  Nov. 2009
!*********************************************************************

  subroutine output(nstep,leap)

  use module_constants
  use module_control,        only: nip, npp, filename_len, hrs_in_month, dt, ArchvStep0
  use fimnamelist,           only: FixedGridOrder,ArchvTimeUnit,	&
                                   yyyymmddhhmm,kdm,ocnonly,coupled
  use module_core_setup,     only: use_write_tasks
  use restart,               only: write_restart
  use module_sfc_variables,only:        &		! air-sea coupling data
     hf2d	&! sensible heat flux (W/m^2)
    ,qf2d	&! latent heat flux (W/m^2)
    ,rsds	&! downward shortwave radiation (W/m^2)
    ,rsus	&! upward shortwave radiation (W/m^2)
    ,rlds	&! downward longwave radiation (W/m^2)
    ,rlus	&! upward longwave radiation (W/m^2)
    ,rain	&! precip in last time step (m/timestep)
    ,t2m2d	&! 2m temperature (K)
    ,q2m2d	 ! 2m vapor mixing ratio (kg/kg)

  use hycom_control,  only: numtr,modesplit,diag_sst
  use hycom_constants,only: wet,ArchvStepOcn,batrop,thref,onem,odepth,	&
                            land_spval,sstproc,nsstpro,KelvC
  use hycom_variables
  use icosio,		     only: icosio_out
  use module_header,         only: headero
! use ersatzpipe

  implicit none
! External variable declarations:
  integer,intent(IN)    :: nstep
  integer,intent(IN)    :: leap				! leap frog time slot
! Local variables:
!SMS$DISTRIBUTE (dh,1) BEGIN
  real      :: aux(nip)
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE (dh,2) BEGIN
  real      :: edgflx(kdm,nip),utotn(kdm,nip),vtotn(kdm,nip),work(kdm,nip)
!SMS$DISTRIBUTE END
  integer       :: i,k,edg,type,time
  real,external :: its2time
  character(len=filename_len),external :: filename,flexflnm
  character char12*12,char4*4
  logical, parameter:: ave_intvl=.true.	! T: write fluxes averaged over "ArchvIntvl"

  time=nint(its2time(nstep))

  utotn(:,:)=0.
  vtotn(:,:)=0.
   work(:,:)=0.
! --- build up time averages for archiving
!SMS$PARALLEL(dh,<i>) BEGIN
!$OMP PARALLEL DO
  do i=1,nip
    if (wet(i) > 0 ) then
      srfx_ave  (i)=srfx_ave  (i)+srfx  (i)
      pmne_ave  (i)=pmne_ave  (i)+pmne  (i)
      prcp_ave  (i)=prcp_ave  (i)+rain  (i)/dt
      swdn_ave  (i)=swdn_ave  (i)+rsds  (i)
      lwdn_ave  (i)=lwdn_ave  (i)+rlds  (i)
       rad_ave  (i)= rad_ave  (i)+rsds(i)+rlds(i)+rsus(i)+rlus(i)
      wspd_ave  (i)=wspd_ave  (i)+wspd  (i)
      taux_ave  (i)=taux_ave  (i)+taux  (i)
      tauy_ave  (i)=tauy_ave  (i)+tauy  (i)
      hf2d_ave  (i)=hf2d_ave  (i)+hf2d  (i)
      qf2d_ave  (i)=qf2d_ave  (i)+qf2d  (i)

      ssht_ave (i)=ssht_ave   (i)+  ssht(i)
      hmixl_ave(i)=hmixl_ave  (i)+ hmixl(i)
      temice_ave(i)=temice_ave(i)+temice(i)
      covice_ave(i)=covice_ave(i)+covice(i)
      thkice_ave(i)=thkice_ave(i)+thkice(i)
      if (coupled) then
        airt_ave (i)=airt_ave (i)+t2m2d(i) - KelvC
        vpmx_ave (i)=vpmx_ave (i)+q2m2d(i)
      end if
      do k=1,kdm
        utotn(k,i)=uclin(k,i,leap)+utrop(leap,i)
        vtotn(k,i)=vclin(k,i,leap)+vtrop(leap,i)
        temp_ave(k,i)=temp_ave(k,i)+temp(k,i)
        saln_ave(k,i)=saln_ave(k,i)+saln(k,i)
        dens_ave(k,i)=dens_ave(k,i)+dens(k,i)
        uvel_ave(k,i)=uvel_ave(k,i)+utotn(k,i)
        vvel_ave(k,i)=vvel_ave(k,i)+vtotn(k,i)
            work(k,i)=  dp    (k,i,leap)/onem	! => m
          dp_ave(k,i)=  dp_ave(k,i)+work(k,i)
      end do
    end if			! wet(i)>0
  enddo
!$OMP END PARALLEL DO
!SMS$PARALLEL END

    if (mod(nstep-ArchvStep0,ArchvStepOcn).ne.0) return
    print *,'now archiving ocean fields at',time,ArchvTimeUnit,		&
    ',   time step',nstep

    if (diag_sst) print '(a/(i5,3x,a))',				&
      'SST tendency contributions in -dsst- array:',			&
      (i,sstproc(i),i=1,nsstpro)
!SMS$PARALLEL(dh,<i>) BEGIN
!$OMP PARALLEL DO
  do i=1,nip
    if (wet(i) > 0 .and. nstep > 0) then	! no ave needed to archive IC
      srfx_ave (i)=srfx_ave (i)/float(ArchvStepOcn)
      pmne_ave (i)=pmne_ave (i)/float(ArchvStepOcn)
      prcp_ave (i)=prcp_ave (i)/float(ArchvStepOcn)
      swdn_ave (i)=swdn_ave (i)/float(ArchvStepOcn)
      lwdn_ave (i)=lwdn_ave (i)/float(ArchvStepOcn)
       rad_ave (i)= rad_ave (i)/float(ArchvStepOcn)
      wspd_ave (i)=wspd_ave (i)/float(ArchvStepOcn)
      taux_ave (i)=taux_ave (i)/float(ArchvStepOcn)
      tauy_ave (i)=tauy_ave (i)/float(ArchvStepOcn)
      qf2d_ave (i)=qf2d_ave (i)/float(ArchvStepOcn)
      hf2d_ave (i)=hf2d_ave (i)/float(ArchvStepOcn)

      ssht_ave (i)=ssht_ave (i)/float(ArchvStepOcn)
      hmixl_ave(i)=hmixl_ave(i)/float(ArchvStepOcn)
      temice_ave(i)=temice_ave(i)/float(ArchvStepOcn)
      covice_ave(i)=covice_ave(i)/float(ArchvStepOcn)
      thkice_ave(i)=thkice_ave(i)/float(ArchvStepOcn)
      do k=1,kdm
        temp_ave(k,i)=temp_ave(k,i)/float(ArchvStepOcn)
        saln_ave(k,i)=saln_ave(k,i)/float(ArchvStepOcn)
        dens_ave(k,i)=dens_ave(k,i)/float(ArchvStepOcn)
        uvel_ave(k,i)=uvel_ave(k,i)/float(ArchvStepOcn)
        vvel_ave(k,i)=vvel_ave(k,i)/float(ArchvStepOcn)
          dp_ave(k,i)=  dp_ave(k,i)/float(ArchvStepOcn)
      end do
      if (coupled) then
        airt_ave (i)=airt_ave (i)/float(ArchvStepOcn)
        vpmx_ave (i)=vpmx_ave (i)/float(ArchvStepOcn)
      end if
    end if	! if (wet(i) > 0 .and. nstep > 0) then
  enddo
!$OMP END PARALLEL DO
!SMS$PARALLEL END

! call findmxmn2(utrop,3,nip,1,'(outp) utrop',wet)
! call findmxmn2(vtrop,3,nip,1,'(outp) vtrop',wet)
! call findmxmn2(ptrop,3,nip,1,'(outp) ptrop',wet)
  call findmxmn1(ssht ,nip,'(outp) ssh ',wet)
  call findmxmn1(curl,nip,'(outp) curl',wet)
  call findmxmn1(curl_ave,nip,'(outp) curl_ave',wet)
  call findmxmn1(ssht_ave,nip,'(outp) ssht_ave',wet)
  call findmxmn1(srfx_ave,nip,'(outp) srfx_ave',wet)
  call findmxmn1(pmne_ave,nip,'(outp) pmne_ave',wet)
  call findmxmn1(prcp_ave,nip,'(outp) prcp_ave',wet)
  call findmxmn1(swdn_ave,nip,'(outp) swdn_ave',wet)
  call findmxmn1(lwdn_ave,nip,'(outp) lwdn_ave',wet)
  call findmxmn1(wspd_ave,nip,'(outp) wspd_ave',wet)
  call findmxmn1(taux_ave,nip,'(outp) taux_ave',wet)
  call findmxmn1(tauy_ave,nip,'(outp) tauy_ave',wet)
  call findmxmn1(hf2d_ave,nip,'(outp) hf2d_ave',wet)
  call findmxmn1(qf2d_ave,nip,'(outp) qf2d_ave',wet)
  call findmxmn1(temice_ave,nip,'(outp) temice_ave',wet)
  call findmxmn1(covice_ave,nip,'(outp) covice_ave',wet)
  call findmxmn1(thkice_ave,nip,'(outp) thkice_ave',wet)
  if (coupled) then
    call findmxmn1(airt_ave,nip,'(outp) airt_ave',wet)
    call findmxmn1(vpmx_ave,nip,'(outp) vpmx_ave',wet)
  end if

  do k=1,kdm
   write (char4,'(a,i2)') 'k=',k
!  call findmxmn3(dp   ,   kdm,nip,2,k,1,'(outp) thkns '//char4,wet)
   call findmxmn2(work ,   kdm,nip,  k,  '(outp) thkns '//char4,wet)
   call findmxmn2(montg,   kdm,nip,  k,  '(outp) montg '//char4,wet)
   call findmxmn2(utotn,   kdm,nip,  k,  '(outp) utotn '//char4,wet)
   call findmxmn2(vtotn,   kdm,nip,  k,  '(outp) vtotn '//char4,wet)
   call findmxmn2(temp ,   kdm,nip,  k,  '(outp) temp  '//char4,wet)
   call findmxmn2(saln ,   kdm,nip,  k,  '(outp) saln  '//char4,wet)
   call findmxmn2(dens ,   kdm,nip,  k,  '(outp) dens  '//char4,wet)
   call findmxmn2(temp_ave,kdm,nip,  k,  '(outp) tave  '//char4,wet)
   call findmxmn2(saln_ave,kdm,nip,  k,  '(outp) save  '//char4,wet)
   call findmxmn2(dens_ave,kdm,nip,  k,  '(outp) dpav  '//char4,wet)
  enddo				! k loop

  do type=1,numtr
    write (char12,'(a,i2)') '(outp) trc',type
    call findmxmn3(passv_tr,kdm,nip,numtr,  1,type,char12//' k=kdm',wet)
    call findmxmn3(passv_tr,kdm,nip,numtr,kdm,type,char12//' k = 1',wet)
    write (char12,'(a,i2.2)') 'ocn_out_tr',type
    call icosio_out (nstep,time,char12,passv_tr(:,:,type),		&
              flexflnm(char12,nstep),headero(char12(9:12),kdm,nstep))
  end do			! tracer type

! --- 3-D fields, instantaneous:
  if (.not. ave_intvl) then
    call icosio_out (nstep,time,'temp',temp,			&
              flexflnm('ocn_out_temp',nstep),headero('temp',kdm,nstep))
    call icosio_out (nstep,time,'saln',saln,			&
              flexflnm('ocn_out_saln',nstep),headero('saln',kdm,nstep))
    call icosio_out (nstep,time,'dens',dens,			&
              flexflnm('ocn_out_dens',nstep),headero('dens',kdm,nstep))
    call icosio_out (nstep,time,'uvel',utotn,			&
              flexflnm('ocn_out_uvel',nstep),headero('uvel',kdm,nstep))
    call icosio_out (nstep,time,'vvel',vtotn,			&
              flexflnm('ocn_out_vvel',nstep),headero('vvel',kdm,nstep))
    call icosio_out (nstep,time,'thik',work,			&
              flexflnm('ocn_out_thik',nstep),headero('thik',kdm,nstep))
    call icosio_out (nstep,time,'mont',montg,			&
              flexflnm('ocn_out_mont',nstep),headero('mont',kdm,nstep))
  end if

! --- 3-D fields, averaged over ArchvTimeUnit
  call icosio_out (nstep,time,'tave',temp_ave,			&
              flexflnm('ocn_out_tave',nstep),headero('tave',kdm,nstep))
  call icosio_out (nstep,time,'save',saln_ave,			&
              flexflnm('ocn_out_save',nstep),headero('save',kdm,nstep))
  call icosio_out (nstep,time,'rave',dens_ave,			&
              flexflnm('ocn_out_rave',nstep),headero('rave',kdm,nstep))
  call icosio_out (nstep,time,'uave',uvel_ave,			&
              flexflnm('ocn_out_uave',nstep),headero('uave',kdm,nstep))
  call icosio_out (nstep,time,'vave',vvel_ave,			&
              flexflnm('ocn_out_vave',nstep),headero('vave',kdm,nstep))
  call icosio_out (nstep,time,'dpav',dp_ave,			&
              flexflnm('ocn_out_dpav',nstep),headero('dpav',kdm,nstep))
  if (diag_sst)								&
  call icosio_out (nstep,time,'dsst',sstndcy,				&
              flexflnm('ocn_out_dsst',nstep),headero('dsst',kdm,nstep))

! --- monthly mass fluxes through individual edges
! do edg=1,npp
!  edgflx(:,:)=mssfx_ave(:,edg,:)
!  write (char4,'(i1,a)') edg,'flx'
!  call icosio_out (nstep,time,char4,edgflx,			&
!             flexflnm('ocn_out_'//char4,nstep),headero('char4',kdm,nstep))
! end do

  call icosio_out (nstep,time,'topo',odepth,  			&
              flexflnm('ocn_out_2D__',nstep),headero('topo',1,nstep))
! --- 2-D fields, instantaneous:
  if (.not. ave_intvl) then
    call icosio_out (nstep,time,'srfh',ssht,				&
              flexflnm('ocn_out_2D__',nstep),headero('srfh',1,nstep))
    call icosio_out (nstep,time,'srfx',srfx,  				&
              flexflnm('ocn_out_2D__',nstep),headero('srfx',1,nstep))
    call icosio_out (nstep,time,'pmne',pmne,  				&
              flexflnm('ocn_out_2D__',nstep),headero('pmne',1,nstep))
    call icosio_out (nstep,time,'ustr',taux,  		        	&
              flexflnm('ocn_out_2D__',nstep),headero('ustr',1,nstep))
    call icosio_out (nstep,time,'vstr',tauy,  		        	&
              flexflnm('ocn_out_2D__',nstep),headero('vstr',1,nstep))
    call icosio_out (nstep,time,'curl',curl,  				&
              flexflnm('ocn_out_2D__',nstep),headero('curl',1,nstep))
    call icosio_out (nstep,time,'tice',temice,  		        &
              flexflnm('ocn_out_2D__',nstep),headero('tice',1,nstep))
    call icosio_out (nstep,time,'cice',covice*100,  		        &
              flexflnm('ocn_out_2D__',nstep),headero('cice',1,nstep))
    call icosio_out (nstep,time,'hice',thkice,  	             	&
              flexflnm('ocn_out_2D__',nstep),headero('hice',1,nstep))
  end if

! --- 2-D fields averaged over ArchvTimeUnit
  call icosio_out (nstep,time,'shav',ssht_ave,  			&
              flexflnm('ocn_out_2D__',nstep),headero('shav',1,nstep))
  call icosio_out (nstep,time,'hfav',srfx_ave,  			&
              flexflnm('ocn_out_2D__',nstep),headero('hfav',1,nstep))
  call icosio_out (nstep,time,'peav',pmne_ave,  			&
              flexflnm('ocn_out_2D__',nstep),headero('peav',1,nstep))
  call icosio_out (nstep,time,'prcp',prcp_ave,  			&
              flexflnm('ocn_out_2D__',nstep),headero('prcp',1,nstep))
  call icosio_out (nstep,time,'swdn',swdn_ave,  			&
              flexflnm('ocn_out_2D__',nstep),headero('swdn',1,nstep))
  call icosio_out (nstep,time,'lwdn',lwdn_ave,  			&
              flexflnm('ocn_out_2D__',nstep),headero('lwdn',1,nstep))
  call icosio_out (nstep,time,'rads', rad_ave,				&
              flexflnm('ocn_out_2D__',nstep),headero('rads',1,nstep))
  call icosio_out (nstep,time,'wspd',wspd_ave,  			&
              flexflnm('ocn_out_2D__',nstep),headero('wspd',1,nstep))
  call icosio_out (nstep,time,'taux',taux_ave,  			&
              flexflnm('ocn_out_2D__',nstep),headero('taux',1,nstep))
  call icosio_out (nstep,time,'tauy',tauy_ave,  			&
              flexflnm('ocn_out_2D__',nstep),headero('tauy',1,nstep))
  call icosio_out (nstep,time,'crla',curl_ave,  			&
              flexflnm('ocn_out_2D__',nstep),headero('crla',1,nstep))
  call icosio_out (nstep,time,'hmxl',hmixl_ave,				&
              flexflnm('ocn_out_2D__',nstep),headero('hmxl',1,nstep))
  call icosio_out (nstep,time,'hfss',hf2d_ave,  			&
              flexflnm('ocn_out_2D__',nstep),headero('hfss',1,nstep))
  call icosio_out (nstep,time,'hfls',qf2d_ave,  			&
              flexflnm('ocn_out_2D__',nstep),headero('hfls',1,nstep))
  call icosio_out (nstep,time,'tiav',temice_ave,			&
              flexflnm('ocn_out_2D__',nstep),headero('tiav',1,nstep))
  call icosio_out (nstep,time,'ciav',covice_ave*100,			&
              flexflnm('ocn_out_2D__',nstep),headero('ciav',1,nstep))
  call icosio_out (nstep,time,'hiav',thkice_ave,			&
              flexflnm('ocn_out_2D__',nstep),headero('hiav',1,nstep))
  if (coupled) then
    call icosio_out (nstep,time,'airt',airt_ave,  			&
              flexflnm('ocn_out_2D__',nstep),headero('airt',1,nstep))
    call icosio_out (nstep,time,'vpmx',vpmx_ave,  			&
              flexflnm('ocn_out_2D__',nstep),headero('vpmx',1,nstep))
  end if

  if (modesplit) then
   aux(:)=utrop(1,:)
   call icosio_out (nstep,time,'utrp',aux,  				&
              flexflnm('ocn_out_2D__',nstep),headero('utrp',1,nstep))
   aux(:)=vtrop(1,:)
   call icosio_out (nstep,time,'vtrp',aux,  				&
              flexflnm('ocn_out_2D__',nstep),headero('vtrp',1,nstep))
   aux(:)=ptrop(1,:)
   call icosio_out (nstep,time,'ptrp',aux,  				&
              flexflnm('ocn_out_2D__',nstep),headero('ptrp',1,nstep))
  end if
  print *,'... archiving completed at  ',time,ArchvTimeUnit,		&
    ',   time step',nstep

!SMS$PARALLEL(dh,i) BEGIN
!$OMP PARALLEL DO
  do i=1,nip
    if (wet(i) > 0 ) then
      curl_ave (i)=0.
      ssht_ave (i)=0.
      srfx_ave (i)=0.
      pmne_ave (i)=0.
      hf2d_ave (i)=0.
      qf2d_ave (i)=0.
      airt_ave (i)=0.
      vpmx_ave (i)=0.
      prcp_ave (i)=0.
      swdn_ave (i)=0.
      lwdn_ave (i)=0.
       rad_ave (i)=0.
      wspd_ave (i)=0.
      taux_ave (i)=0.
      tauy_ave (i)=0.
      hmixl_ave(i)=0.
      mssfx_ave(:,:,i)=0.
      temice_ave(i)=0.
      covice_ave(i)=0.
      thkice_ave(i)=0.
      sstndcy (:,i)=0.
      do k=1,kdm
        temp_ave(k,i)=0.
        saln_ave(k,i)=0.
        dens_ave(k,i)=0.
        uvel_ave(k,i)=0.
        vvel_ave(k,i)=0.
          dp_ave(k,i)=0.
      end do
    end if
  end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END

! if (nstep.le.1 .or. coupled) then
!   print *,'outputting CORE- or CORE-like monthly forcing fields...'
!   call icosio_out (0,0,'swdn',swdn_12mo,12,				&
!              flexflnm('CORE_mthly_swdn_',0),headero('swdn',12,0))
!   call icosio_out (0,0,'lwdn',lwdn_12mo,12,				&
!              flexflnm('CORE_mthly_lwdn_',0),headero('lwdn',12,0))
!   call icosio_out (0,0,'airt',airt_12mo,12,				&
!              flexflnm('CORE_mthly_airt_',0),headero('airt',12,0))
!   call icosio_out (0,0,'vpmx',vpmx_12mo,12,				&
!              flexflnm('CORE_mthly_vpmx_',0),headero('vpmx',12,0))
!   call icosio_out (0,0,'wspd',wspd_12mo,12,				&
!              flexflnm('CORE_mthly_wspd_',0),headero('wspd',12,0))
!   call icosio_out (0,0,'uwnd',uwnd_12mo,12,				&
!              flexflnm('CORE_mthly_uwnd_',0),headero('uwnd',12,0))
!   call icosio_out (0,0,'vwnd',vwnd_12mo,12,				&
!              flexflnm('CORE_mthly_vwnd_',0),headero('vwnd',12,0))
!   call icosio_out (0,0,'prcp',prcp_12mo,12,				&
!              flexflnm('CORE_mthly_prcp_',0),headero('prcp',12,0))
!   print *,'monthly-avgd CORE-type forcing fields saved in files named CORE_mthly_...'
! end if
   
  return
  end subroutine output
end module hycom_output
