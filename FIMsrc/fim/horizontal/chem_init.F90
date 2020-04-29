module module_fim_chem_init

  implicit none

#include <gptl.inc>

contains

  subroutine chem_init
!*********************************************************************
!       Loads the initial variables and constants for the chemsitry 
!       component.  
!       Tom Henderson           November, 2008
!*********************************************************************

!TODO:  move CallChemistry to a chem module
    use module_control,only: nip,numphr,ntra,ntrb
    use fimnamelist   ,only: nvl,glvl,curve,readrestart
    use module_constants, only: perm
    use module_wrf_control,only: ids,ide,jds,jde,kds,kde, &
         ims,ime,jms,jme,kms,kme, &
         its,ite,jts,jte,kts,kte, &
         num_moist,CallChemistry, &
         CallBiom,num_chem,nvl_gocart
    use module_chem_constants,only: p_gocart
    use module_chem_variables,only: ero1,ero2,ero3,rcav,dm0,pm25,p10,emiss_ab,   &
         emiss_oc,emiss_bc,emiss_sulf,emiss_pm25,emiss_pm10,oh_backgd,sandfrac,clayfrac,      &
         emiss_ash_mass,emiss_ash_height,emiss_ash_dt,&
         emiss_tr_height,emiss_tr_mass,emiss_tr_dt,   &
         h2o2_backgd,no3_backgd,emiss_abu,plumestuff
    use module_variables,only: tr3d,trdp,dp3d
    use module_chem_alloc,only: chem_alloc,chem_alloc2
    use module_wrf_share ,only: wrf_set_array_bounds

    !TODO:  clean up duplication in these modules and add "only"
    ! for chemistry and WRF physics namelists:
    use module_chem_namelist_defaults
    ! contains config_flags
    use fimnamelist
    ! TBH:  Ignore these so PPP doesn't have to translate them
    use module_species_decs
    use module_set_wrfphys
    use units, only: getunit, returnunit

    ! Local variables
    integer       :: ipn	! Index for icos point number
    integer       :: ivl	! Index for vertical level
    integer       :: nv,k,idx,nv_g,ierr,ibegin
    logical       :: debug ! control debug prints
    character(64) :: filename
    character(20) :: dum
    character(12) :: dum2
    character(80)            :: header(10)
    real :: maxv
    real :: serial1(nip),serial2(nvl_gocart,nip)
!SMS$DISTRIBUTE (dh,2) BEGIN
    real, allocatable :: dummy(:,:)
!SMS$DISTRIBUTE END
    integer :: unitno
    integer :: ret
    real*8 :: tchem_init
!SMS$insert integer :: mype
!sms$comm_rank(mype)

    ret = gptlstart ('chem_init')

!SMS$PARALLEL(dh, ipn) BEGIN

!TODO:  Remove implicit dependence on ReadRestart being initilialized in 
!TODO:  dyn_init().  

    debug=.false.
    ibegin=ntra+num_moist-3
    if(num_moist.gt.3)ibegin=ntra+num_moist-2

    if (readrestart .and. chem_opt /= 0) then
      write(6,*) 'chem_init: restarting with chemistry enabled not yet supported. Stopping'
      call flush(6)
      stop
    end if

    if (debug) print *,'DEBUG chem_init:  begin, ntra,num_moist,num_chem = ',ntra,num_moist,num_chem

    ash_height =-999.
    ash_mass = -999.

    unitno = getunit ()
    if (unitno < 0) then
      print*,'chem_init: getunit failed: stopping'
      stop
    end if

    if (chem_opt.eq.0) then
      ret = gptlstop ('chem_init')
      return
    end if

    call wrf_set_array_bounds(nvl,nip, &
         ids,ide,jds,jde,kds,kde, &
         ims,ime,jms,jme,kms,kme, &
         its,ite,jts,jte,kts,kte)

    config_flags%chem_opt = chem_opt
    config_flags%chem_in_opt = chem_in_opt
    config_flags%dust_opt = dust_opt
    config_flags%dmsemis_opt = dmsemis_opt
    config_flags%seas_opt = seas_opt
    config_flags%biomass_burn_opt = biomass_burn_opt

    ! Always allocate arrays that are passed through argument lists...  
    print *,'DEBUG:  chem_init():  calling chem_alloc(',chem_opt,')'
    call chem_alloc(chem_opt,aer_ra_feedback)
    print *,'DEBUG:  chem_init():  back from chem_alloc()'

    if(chem_opt > 0) then
      call set_species

      ! read chem data
      if(chem_opt >= 300 .and. chem_opt < 500) then
!SMS$SERIAL (<perm,IN>,<oh_backgd,h2o2_backgd,no3_backgd,p_gocart,emiss_ab,dm0,ero1,ero2,ero3,sandfrac,clayfrac,OUT>:default=ignore) BEGIN
        write(6,*)'reading gocart background fields'
        call flush(6)
        filename = "oh.dat"
        open (unit=unitno, file=filename, form="unformatted", action='read')
        read (unitno) p_gocart
        read (unitno) serial2
        write(6,*)'minimum oh-value ',minval(serial2(:,15)),maxval(serial2(:,15))
        call flush(6)
        close (unitno)
        do ipn=1,nip
          oh_backgd(:,ipn) = serial2(:,perm(ipn))
        end do

        filename = "h2o2.dat"
        open (unit=unitno, file=filename, form="unformatted", action='read')
        read (unitno) p_gocart
        read (unitno) serial2
        write(6,*)'minimum h2o2-value ',minval(serial2(:,15)),maxval(serial2(:,15))
        call flush(6)
        close(unitno)
        do ipn=1,nip
          h2o2_backgd(:,ipn) = serial2(:,perm(ipn))
        end do

        filename = "no3.dat"
        open (unit=unitno, file=filename, form="unformatted", action='read')
        read (unitno) p_gocart
        read (unitno) serial2
        write(6,*)'p_gocart',p_gocart
        write(6,*)'minimum no3-value ',minval(serial2(:,15)),maxval(serial2(:,15))
        call flush(6)
        close(unitno)
        write(6,*)'reading chemistry emissions files for gocart'
        call flush(6)
        do ipn=1,nip
          no3_backgd(:,ipn) = serial2(:,perm(ipn))
        end do

        filename = "e_bc.dat"
        open (unit=unitno, file=filename, form="unformatted", action='read')
        read (unitno) serial1
        close (unitno)
        maxv=maxval(serial1)
        write(6,*)'maxv on input for bc_ant = ',maxv
        call flush(6)
        do ipn=1,nip
          emiss_ab(ipn,p_e_bc) = serial1(perm(ipn))
        end do

        filename = "e_oc.dat"
        open (unit=unitno, file=filename, form="unformatted", action='read')
        read (unitno) serial1
        close (unitno)
        maxv=maxval(serial1)
        write(6,*)'maxv on input for oc_ant = ',maxv
        call flush(6)
        do ipn=1,nip
          emiss_ab(ipn,p_e_oc) = serial1(perm(ipn))
        end do

        filename = "e_sulf.dat"
        open (unit=unitno, file=filename, form="unformatted", action='read')
        read (unitno) serial1
        close(unitno)
        maxv=maxval(serial1)
        write(6,*)'maxv on input for sulf_ant = ',maxv
        call flush(6)
        do ipn=1,nip
          emiss_ab(ipn,p_e_sulf) = serial1(perm(ipn))
        end do

        filename = "e_pm_25.dat"
        open (unit=unitno, file=filename, form="unformatted", action='read')
        read (unitno) serial1
        close(unitno)
        maxv=maxval(serial1)
        write(6,*)'maxv on input for pm25_ant = ',maxv
        call flush(6)
        do ipn=1,nip
          emiss_ab(ipn,p_e_pm_25) = serial1(perm(ipn))
        end do

        filename = "e_pm_10.dat"
        open (unit=unitno, file=filename, form="unformatted", action='read')
        read (unitno) serial1
        close(unitno)
        maxv=maxval(serial1)
        write(6,*)'maxv on input for pm10_ant = ',maxv
        call flush(6)
        do ipn=1,nip
          emiss_ab(ipn,p_e_pm_10) = serial1(perm(ipn))
        end do

        filename = "dm0.dat"
        open (unit=unitno, file=filename, form='unformatted', action='read')
        read (unitno) serial1
        close (unitno)
        do ipn=1,nip
          dm0(ipn)=serial1(perm(ipn))
        end do

        filename = "erod1.dat"
        open (unit=unitno, file=filename, form='unformatted', action='read')
        read (unitno) serial1
        close (unitno)
        do ipn=1,nip
          ero1(ipn) = serial1(perm(ipn))
        end do

        filename = "erod2.dat"
        open (unit=unitno, file=filename, form='unformatted', action='read')
        read (unitno) serial1
        close (unitno)
        do ipn=1,nip
          ero2(ipn) = serial1(perm(ipn))
        end do

        filename = "erod3.dat"
        open (unit=unitno, file=filename, form='unformatted', action='read')
        read (unitno) serial1
        close (unitno)
        do ipn=1,nip
          ero3(ipn) = serial1(perm(ipn))
        end do
!
! AFWA dust option
!
        if(dust_opt == 3) then
          filename = "sand.dat"
          open (unit=unitno, file=filename, form='unformatted', action='read')
          read (unitno) serial1
          close (unitno)
          do ipn=1,nip
            sandfrac(ipn) = serial1(perm(ipn))
          end do

          filename = "clay.dat"
          open (unit=unitno, file=filename, form='unformatted', action='read')
          read (unitno) serial1
          close (unitno)
          do ipn=1,nip
            clayfrac(ipn) = serial1(perm(ipn))
          end do
        endif

        filename = "e_so2.dat"
        open (unit=unitno, file=filename, form='unformatted', action='read')
        read (unitno) serial1
        close (unitno)
        maxv=maxval(serial1)
        write(6,*)'maxv on input for so2_ant = ',maxv
        call flush(6)
        do ipn=1,nip
          emiss_ab(ipn,p_e_so2) = serial1(perm(ipn))
        enddo
!SMS$SERIAL END
      endif ! chem_opt >= 300 and chem_opt < 500
!
!!!!!!! biomassburning next, lots of arrays!
!
      if(biomass_burn_opt > 0 ) then
!SMS$SERIAL (<perm,IN>,<emiss_abu,plumestuff,OUT>:default=ignore) BEGIN
        filename = "ebu_oc.dat"
        open (unit=unitno, file=filename, form='unformatted', action='read')
        read (unitno) serial1
        close (unitno)
        maxv=maxval(serial1(:))
        write(6,*)'maxv on input2 for oc_biom = ',maxv
        call flush(6)
        do ipn=1,nip
	  emiss_abu(ipn,p_e_oc) = serial1(perm(ipn))
        end do

        filename = "ebu_bc.dat"
        open (unit=unitno, file=filename, form='unformatted', action='read')
        read (unitno) serial1
        close (unitno)
        maxv=maxval(serial1(:))
        write(6,*)'maxv on input2 for bc_bion = ',maxv
        call flush(6)
        do ipn=1,nip
          emiss_abu(ipn,p_e_bc) = serial1(perm(ipn))
        end do

        filename = "ebu_so2.dat"
        open (unit=unitno, file=filename, form='unformatted', action='read')
        read (unitno) serial1
        close (unitno)
        maxv=maxval(serial1(:))
        write(6,*)'maxv on input2 for bso2 = ',maxv
        call flush(6)
        do ipn=1,nip
          emiss_abu(ipn,p_e_so2) = serial1(perm(ipn))
        end do

        filename = "ebu_sulf.dat"
        open (unit=unitno, file=filename, form='unformatted', action='read')
        read (unitno) serial1
        close (unitno)
        maxv=maxval(serial1(:))
        write(6,*)'maxv on input2 for bsulf = ',maxv
        call flush(6)
        do ipn=1,nip
          emiss_abu(ipn,p_e_sulf) = serial1(perm(ipn))
        end do

        filename = "ebu_pm25.dat"
        open (unit=unitno, file=filename, form='unformatted', action='read')
        read (unitno) serial1
        close (unitno)
        maxv=maxval(serial1(:))
        write(6,*)'maxv on input2 for bpm25 = ',maxv
        call flush(6)
        do ipn=1,nip
          emiss_abu(ipn,p_e_pm_25) = serial1(perm(ipn))
        end do

        filename = "ebu_pm10.dat"
        open (unit=unitno, file=filename, form='unformatted', action='read')
        read (unitno) serial1
        close (unitno)
        maxv=maxval(serial1(:))
        write(6,*)'maxv on input2 for bpm10 = ',maxv
        call flush(6)
        do ipn=1,nip
          emiss_abu(ipn,p_e_pm_10) = serial1(perm(ipn))
        end do

        filename = "plumestuff.dat"
        open (unit=unitno, file=filename, form='unformatted', action='read')
        do k=1,8
          read (unitno) serial1
          do ipn=1,nip
            plumestuff(ipn,k) = serial1(perm(ipn))
          enddo
          maxv=maxval(plumestuff(:,k))
          write(6,*)'maxv on input2 for plumestuff(k) = ',k,maxv
          call flush(6)
        enddo
        close (unitno)
!SMS$SERIAL END
      endif ! biomass_burn
!
! read volcanic stuff if necessary
!
      if((chem_opt == 316 .or. chem_opt == 317 .or. chem_opt == 502 .or. chem_opt == 300)) then
        if(ash_mass  >= -900. )then
!SMS$SERIAL (<perm,IN>,<emiss_ash_mass,emiss_ash_height,emiss_ash_dt,OUT>:default=ignore) BEGIN
        filename = "volcanic.dat"
        open (unit=unitno, file=filename, form='unformatted', action='read')
        read (unitno) nv_g
        print *,nv_g
        read (unitno) dum
        print *,dum
        read (unitno) dum2
        print *,dum2
        read (unitno) serial1
        do ipn=1,nip
          emiss_ash_mass(ipn) = serial1(perm(ipn))
        end do
        read (unitno) serial1
        do ipn=1,nip
          emiss_ash_height(ipn) = serial1(perm(ipn))
        end do
        read (unitno) serial1
        do ipn=1,nip
          emiss_ash_dt(ipn) = serial1(perm(ipn))
        end do
        close (unitno)
!SMS$SERIAL END

        maxv=maxval(emiss_ash_mass(:))
        write(6,*)'maxv on emissions input for ashmass = ',maxv
        maxv=maxval(emiss_ash_height(:))
        write(6,*)'maxv on emissions input for ashheight = ',maxv
        maxv=maxval(emiss_ash_dt(:))
        write(6,*)'maxv on emissions input for ashdt = ',maxv
        endif
        if(ash_mass.gt.-100)then
          write(0,*)'using namelist value for ash_mass'
          do ipn=1,nip
            if(emiss_ash_mass(ipn).le.0.)cycle
!
! overwrite ash_mass if nameist value exists
!
            emiss_ash_mass(ipn)=ash_mass
          enddo
        endif
!         maxv=maxval(emiss_ash_mass(:))
!         write(6,*)'maxv on emissions input for ashmass = ',maxv
!         call flush(6)
!         read (unitno) emiss_bc
!           emiss_ash_height(ipn)=emiss_bc(ipn)
!
! overwrite ash_height if nameist value exists
!
        if(ash_height .gt. 0.) then
          write(0,*)'using namelist value for ash_height'
          do ipn=1,nip
            if(emiss_ash_mass(ipn).le.0.)cycle
            emiss_ash_height(ipn)=ash_height
            if(emiss_ash_height(ipn) .lt. 1.) emiss_ash_dt(ipn)=0.
          enddo
        endif
            if(ash_height .lt. -990.) then
              write(0,*)'resetting all ash variables to zero'
              do ipn=1,nip
                 emiss_ash_height(ipn)=0.
                 emiss_ash_mass(ipn)=0.
                 emiss_ash_dt(ipn)=-10.
              enddo
            endif
!           if(ash_height .lt. -1.) emiss_ash_height(ipn)=0.
!         call flush(6)
!         read (unitno) emiss_bc
!         do ipn=1,nip
!           emiss_ash_dt(ipn)=emiss_bc(ipn)
!         enddo
!         maxv=maxval(emiss_ash_dt(:))
!         write(6,*)'maxv on emissions input for duration = ',maxv
!         maxv=maxval(emiss_ash_mass(:))
!         write(6,*)'maxv on emissions input for ashmass = ',maxv
!         maxv=maxval(emiss_ash_height(:))
!         write(6,*)'maxv on emissions input for ashheight = ',maxv
!         maxv=maxval(emiss_ash_dt(:))
!         write(6,*)'maxv on emissions input for ashdt = ',maxv
!         call flush(6)
      endif ! (chem_opt.eq.316.or.chem_opt.eq.317) 
!     do k=1,8
!     do ipn=1,nip
!         plumestuff(ipn,p_e_oc)=plumes(ipn,8)
!     enddo
!     enddo

! Initialize chem arrays
      do nv=ntra+1,ntra+ntrb
        do ipn=1,nip
          do ivl=1,nvl
            tr3d(ivl,ipn,nv) = 1.e-16
            if(chem_opt == 501 .or. chem_opt == 502) tr3d(ivl,ipn,nv) = 0
            if(chem_opt == 500 ) tr3d(ivl,ipn,nv) = 390.
          enddo
        enddo
      enddo
!
! cloud water done twice ????
!
!       do ipn=1,nip
!         do ivl=1,nvl
!           tr3d(ivl,ipn,5) = tr3d(ivl,ipn,3)
!           trdp(ivl,ipn,5) = tr3d(ivl,ipn,3)*dp3d(ivl,ipn)
!         enddo
!       enddo
      do ipn=1,nip
        do ivl=1,nvl
          pm25(ivl,ipn) = 0.
          p10(ivl,ipn) = 0.
        enddo
      enddo
      do ipn=1,nip
        rcav(ipn) = 0.
      enddo

      call chem_alloc2(chem_opt,aer_ra_feedback,bio_emiss_opt,biomass_burn_opt,kemit)

      CallChemistry = max(1,numphr*(int(Chemdt+.01)*60)/3600)
      Callbiom      = max(1,numphr*(int(PLUMERISEFIRE_FRQ+.01)*60)/3600)

      if(chem_in_opt == 1 ) then
        call flush(6)
        allocate(  dummy(nvl,nip))     ! if chem_in_opt=1, read old chem
        call flush(6)
!
! read previous volcanic ash forecast
!
        if(chem_opt == 502 ) then
          open (unit=unitno, file='vash1.in', form="unformatted", action='read')
          read (unitno) header
          read (unitno) dummy
          close (unitno)
          maxv=maxval(dummy(:,:))
          write(6,*)'maxv on input2 for vash1 = ',maxv
          call flush(6)
          do ipn=1,nip
            do ivl=1,nvl
              tr3d(ivl,ipn,ibegin+p_vash_1) = dummy(ivl,ipn)
              trdp(ivl,ipn,ibegin+p_vash_1) = dummy(ivl,ipn)*dp3d(ivl,ipn)
            enddo
          enddo
          open (unit=unitno, file='vash2.in', form="unformatted", action='read')
          read (unitno) header
          read (unitno) dummy
          close (unitno)
          maxv=maxval(dummy(:,:))
          write(6,*)'maxv on input2 for vash2 = ',maxv
          call flush(6)
          do ipn=1,nip
            do ivl=1,nvl
              tr3d(ivl,ipn,ibegin+p_vash_2) = dummy(ivl,ipn)
              trdp(ivl,ipn,ibegin+p_vash_2) = dummy(ivl,ipn)*dp3d(ivl,ipn)
            enddo
          enddo
          open (unit=unitno, file='vash3.in', form="unformatted", action='read')
          read (unitno) header
          read (unitno) dummy
          close (unitno)
          maxv=maxval(dummy(:,:))
          write(6,*)'maxv on input2 for vash3 = ',maxv
          call flush(6)
          do ipn=1,nip
            do ivl=1,nvl
              tr3d(ivl,ipn,ibegin+p_vash_3) = dummy(ivl,ipn)
              trdp(ivl,ipn,ibegin+p_vash_3) = dummy(ivl,ipn)*dp3d(ivl,ipn)
            enddo
          enddo
          open (unit=unitno, file='vash4.in', form="unformatted", action='read')
          read (unitno) header
          read (unitno) dummy
          close (unitno)
          maxv=maxval(dummy(:,:))
          write(6,*)'maxv on input2 for vash4 = ',maxv
          call flush(6)
          do ipn=1,nip
            do ivl=1,nvl
              tr3d(ivl,ipn,ibegin+p_vash_4) = dummy(ivl,ipn)
              trdp(ivl,ipn,ibegin+p_vash_4) = dummy(ivl,ipn)*dp3d(ivl,ipn)
            enddo
          enddo
        endif  ! chem_opt=502
        !
        ! this shpould just be a loop that automatically reads the stuff
        ! next is for volcanic ash (4 or 10 bins) + gocart
        if((chem_opt == 316 .or. chem_opt == 317) .and. ash_mass .ne. 0.) then
          open (unit=unitno, file='vash1.in', form="unformatted", action='read')
          read (unitno) header
          read (unitno) dummy
          close (unitno)
          maxv=maxval(dummy(:,:))
          write(6,*)'maxv on input2 for vash1 = ',maxv
          call flush(6)
          do ipn=1,nip
            do ivl=1,nvl
              tr3d(ivl,ipn,ibegin+p_vash_1) = dummy(ivl,ipn)
              trdp(ivl,ipn,ibegin+p_vash_1) = dummy(ivl,ipn)*dp3d(ivl,ipn)
            enddo
          enddo
          open (unit=unitno, file='vash2.in', form="unformatted", action='read')
          read (unitno) header
          read (unitno) dummy
          close (unitno)
          maxv=maxval(dummy(:,:))
          write(6,*)'maxv on input2 for vash2 = ',maxv
          call flush(6)
          do ipn=1,nip
            do ivl=1,nvl
              tr3d(ivl,ipn,ibegin+p_vash_2) = dummy(ivl,ipn)
              trdp(ivl,ipn,ibegin+p_vash_2) = dummy(ivl,ipn)*dp3d(ivl,ipn)
            enddo
          enddo
          open (unit=unitno, file='vash3.in', form="unformatted", action='read')
          read (unitno) header
          read (unitno) dummy
          close (unitno)
          maxv=maxval(dummy(:,:))
          write(6,*)'maxv on input2 for vash3 = ',maxv
          call flush(6)
          do ipn=1,nip
            do ivl=1,nvl
              tr3d(ivl,ipn,ibegin+p_vash_3) = dummy(ivl,ipn)
              trdp(ivl,ipn,ibegin+p_vash_3) = dummy(ivl,ipn)*dp3d(ivl,ipn)
            enddo
          enddo
          open (unit=unitno, file='vash4.in', form="unformatted", action='read')
          read (unitno) header
          read (unitno) dummy
          close (unitno)
          maxv=maxval(dummy(:,:))
          write(6,*)'maxv on input2 for vash4 = ',maxv
          call flush(6)
          do ipn=1,nip
            do ivl=1,nvl
              tr3d(ivl,ipn,ibegin+p_vash_4) = dummy(ivl,ipn)
              trdp(ivl,ipn,ibegin+p_vash_4) = dummy(ivl,ipn)*dp3d(ivl,ipn)
            enddo
          enddo
        endif

        if(chem_opt == 316 ) then
          open (unit=unitno, file='vash5.in', form="unformatted", action='read')
          read (unitno) header
          read (unitno) dummy
          close (unitno)
          maxv=maxval(dummy(:,:))
          write(6,*)'maxv on input2 for vash5 = ',maxv
          call flush(6)
          do ipn=1,nip
            do ivl=1,nvl
              tr3d(ivl,ipn,ibegin+p_vash_5) = dummy(ivl,ipn)
              trdp(ivl,ipn,ibegin+p_vash_5) = dummy(ivl,ipn)*dp3d(ivl,ipn)
            enddo
          enddo
          open (unit=unitno, file='vash6.in', form="unformatted", action='read')
          read (unitno) header
          read (unitno) dummy
          close (unitno)
          maxv=maxval(dummy(:,:))
          write(6,*)'maxv on input2 for vash6 = ',maxv
          call flush(6)
          do ipn=1,nip
            do ivl=1,nvl
              tr3d(ivl,ipn,ibegin+p_vash_6) = dummy(ivl,ipn)
              trdp(ivl,ipn,ibegin+p_vash_6) = dummy(ivl,ipn)*dp3d(ivl,ipn)
            enddo
          enddo
          open (unit=unitno, file='vash7.in', form="unformatted", action='read')
          read (unitno) header
          read (unitno) dummy
          close (unitno)
          maxv=maxval(dummy(:,:))
          write(6,*)'maxv on input2 for vash7 = ',maxv
          call flush(6)
          do ipn=1,nip
            do ivl=1,nvl
              tr3d(ivl,ipn,ibegin+p_vash_7) = dummy(ivl,ipn)
              trdp(ivl,ipn,ibegin+p_vash_7) = dummy(ivl,ipn)*dp3d(ivl,ipn)
            enddo
          enddo
          open (unit=unitno, file='vash8.in', form="unformatted", action='read')
          read (unitno) header
          read (unitno) dummy
          close (unitno)
          maxv=maxval(dummy(:,:))
          write(6,*)'maxv on input2 for vash8 = ',maxv
          call flush(6)
          do ipn=1,nip
            do ivl=1,nvl
              tr3d(ivl,ipn,ibegin+p_vash_8) = dummy(ivl,ipn)
              trdp(ivl,ipn,ibegin+p_vash_8) = dummy(ivl,ipn)*dp3d(ivl,ipn)
            enddo
          enddo
          open (unit=unitno, file='vash9.in', form="unformatted", action='read')
          read (unitno) header
          read (unitno) dummy
          close (unitno)
          maxv=maxval(dummy(:,:))
          write(6,*)'maxv on input2 for vash9 = ',maxv
          call flush(6)
          do ipn=1,nip
            do ivl=1,nvl
              tr3d(ivl,ipn,ibegin+p_vash_9) = dummy(ivl,ipn)
              trdp(ivl,ipn,ibegin+p_vash_9) = dummy(ivl,ipn)*dp3d(ivl,ipn)
            enddo
          enddo
          open (unit=unitno, file='vash10.in', form="unformatted", action='read')
          read (unitno) header
          read (unitno) dummy
          close (unitno)
          maxv=maxval(dummy(:,:))
          write(6,*)'maxv on input2 for vash10 = ',maxv
          call flush(6)
          do ipn=1,nip
            do ivl=1,nvl
              tr3d(ivl,ipn,ibegin+p_vash_10) = dummy(ivl,ipn)
              trdp(ivl,ipn,ibegin+p_vash_10) = dummy(ivl,ipn)*dp3d(ivl,ipn)
            enddo
          enddo
        endif             !  volcanoes
!
! GOCART options
!
        if(chem_opt >= 300 .and. chem_opt < 500)then
          open (unit=unitno, file='so2.in', form="unformatted", action='read')
          read (unitno) header
          read (unitno) dummy
          close (unitno)
          maxv=maxval(dummy(:,:))
          write(6,*)'maxv on input2 for so2 = ',maxv
          call flush(6)
          do ipn=1,nip
            do ivl=1,nvl
              tr3d(ivl,ipn,ibegin+p_so2) = dummy(ivl,ipn)
              trdp(ivl,ipn,ibegin+p_so2) = dummy(ivl,ipn)*dp3d(ivl,ipn)
            enddo
          enddo
          open (unit=unitno, file='sulf.in', form="unformatted", action='read')
          read (unitno) header
          read (unitno) dummy
          close (unitno)
          maxv=maxval(dummy(:,:))
          write(6,*)'maxv on input2 for sulf = ',maxv
          write(6,*)'p_so2,p_sulf = ',p_so2,p_sulf
          call flush(6)
          do ipn=1,nip
            do ivl=1,nvl
              tr3d(ivl,ipn,ibegin+p_sulf) = dummy(ivl,ipn)
              trdp(ivl,ipn,ibegin+p_sulf) = dummy(ivl,ipn)*dp3d(ivl,ipn)
            enddo
          enddo
          open (unit=unitno, file='dms.in', form="unformatted", action='read')
          read (unitno) header
          read (unitno) dummy
          close (unitno)
          maxv=maxval(dummy(:,:))
          write(6,*)'maxv on input2 for dms = ',maxv
          call flush(6)
          do ipn=1,nip
            do ivl=1,nvl
              tr3d(ivl,ipn,ibegin+p_dms) = dummy(ivl,ipn)
              trdp(ivl,ipn,ibegin+p_dms) = dummy(ivl,ipn)*dp3d(ivl,ipn)
            enddo
          enddo
          open (unit=unitno, file='msa.in', form="unformatted", action='read')
          read (unitno) header
          read (unitno) dummy
          close (unitno)
          maxv=maxval(dummy(:,:))
          write(6,*)'maxv on input2 for msa = ',maxv
          do ipn=1,nip
            do ivl=1,nvl
              tr3d(ivl,ipn,ibegin+p_msa) = dummy(ivl,ipn)
              trdp(ivl,ipn,ibegin+p_msa) = dummy(ivl,ipn)*dp3d(ivl,ipn)
            enddo
          enddo
          open (unit=unitno, file='p10.in', form="unformatted", action='read')
          read (unitno) header
          read (unitno) dummy
          close (unitno)
          maxv=maxval(dummy(:,:))
          write(6,*)'maxv on input2 for p10 = ',maxv
          call flush(6)
          do ipn=1,nip
            do ivl=1,nvl
              tr3d(ivl,ipn,ibegin+p_p10) = dummy(ivl,ipn)
              trdp(ivl,ipn,ibegin+p_p10) = dummy(ivl,ipn)*dp3d(ivl,ipn)
            enddo
          enddo
          open (unit=unitno, file='p25.in', form="unformatted", action='read')
          read (unitno) header
          read (unitno) dummy
          close (unitno)
          maxv=maxval(dummy(:,:))
          write(6,*)'maxv on input2 for p25 = ',maxv
          call flush(6)
          do ipn=1,nip
            do ivl=1,nvl
              tr3d(ivl,ipn,ibegin+p_p25) = dummy(ivl,ipn)
              trdp(ivl,ipn,ibegin+p_p25) = dummy(ivl,ipn)*dp3d(ivl,ipn)
            enddo
          enddo
          open (unit=unitno, file='bc1.in', form="unformatted", action='read')
          read (unitno) header
          read (unitno) dummy
          close (unitno)
          maxv=maxval(dummy(:,:))
          write(6,*)'maxv on input2 for bc1 = ',p_bc1,p_bc2,p_p25,maxv
          call flush(6)
          do ipn=1,nip
            do ivl=1,nvl
              tr3d(ivl,ipn,ibegin+p_bc1) = dummy(ivl,ipn)
              trdp(ivl,ipn,ibegin+p_bc1) = dummy(ivl,ipn)*dp3d(ivl,ipn)
            enddo
          enddo
          open (unit=unitno, file='bc2.in', form="unformatted", action='read')
          read (unitno) header
          read (unitno) dummy
          close (unitno)
          maxv=maxval(dummy(:,:))
          write(6,*)'maxv on input2 for bc2 = ',maxv
          call flush(6)
          do ipn=1,nip
            do ivl=1,nvl
              tr3d(ivl,ipn,ibegin+p_bc2) = dummy(ivl,ipn)
              trdp(ivl,ipn,ibegin+p_bc2) = dummy(ivl,ipn)*dp3d(ivl,ipn)
            enddo
          enddo
          open (unit=unitno, file='oc1.in', form="unformatted", action='read')
          read (unitno) header
          read (unitno) dummy
          close (unitno)
          maxv=maxval(dummy(:,:))
          write(6,*)'maxv on input2 for oc1 = ',maxv
          call flush(6)
          do ipn=1,nip
            do ivl=1,nvl
              tr3d(ivl,ipn,ibegin+p_oc1) = dummy(ivl,ipn)
              trdp(ivl,ipn,ibegin+p_oc1) = dummy(ivl,ipn)*dp3d(ivl,ipn)
            enddo
          enddo
          open (unit=unitno, file='oc2.in', form="unformatted", action='read')
          read (unitno) header
          read (unitno) dummy
          close (unitno)
          maxv=maxval(dummy(:,:))
          write(6,*)'maxv on input2 for oc2 = ',maxv
          call flush(6)
          do ipn=1,nip
            do ivl=1,nvl
              tr3d(ivl,ipn,ibegin+p_oc2) = dummy(ivl,ipn)
              trdp(ivl,ipn,ibegin+p_oc2) = dummy(ivl,ipn)*dp3d(ivl,ipn)
            enddo
          enddo
          open (unit=unitno, file='dust1.in', form="unformatted", action='read')
          read (unitno) header
          read (unitno) dummy
          close (unitno)
          maxv=maxval(dummy(:,:))
          write(6,*)'maxv on input2 for dust1 = ',maxv
          call flush(6)
          do ipn=1,nip
            do ivl=1,nvl
              tr3d(ivl,ipn,ibegin+p_dust_1) = dummy(ivl,ipn)
              trdp(ivl,ipn,ibegin+p_dust_1) = dummy(ivl,ipn)*dp3d(ivl,ipn)
            enddo
          enddo
!
! extra arrays for pure GOCART
!
          if(chem_opt == 300 ) then
            open (unit=unitno, file='dust2.in', form="unformatted", action='read')
            read (unitno) header
            read (unitno) dummy
            close (unitno)
            maxv=maxval(dummy(:,:))
            write(6,*)'maxv on input2 for dust2 = ',maxv
            call flush(6)
            do ipn=1,nip
              do ivl=1,nvl
                tr3d(ivl,ipn,ibegin+p_dust_2) = dummy(ivl,ipn)
                trdp(ivl,ipn,ibegin+p_dust_2) = dummy(ivl,ipn)*dp3d(ivl,ipn)
              enddo
            enddo
            open (unit=unitno, file='dust3.in', form="unformatted", action='read')
            read (unitno) header
            read (unitno) dummy
            close (unitno)
            maxv=maxval(dummy(:,:))
            write(6,*)'maxv on input2 for dust3 = ',maxv
            call flush(6)
            do ipn=1,nip
              do ivl=1,nvl
                tr3d(ivl,ipn,ibegin+p_dust_3) = dummy(ivl,ipn)
                trdp(ivl,ipn,ibegin+p_dust_3) = dummy(ivl,ipn)*dp3d(ivl,ipn)
              enddo
            enddo
            open (unit=unitno, file='dust4.in', form="unformatted", action='read')
            read (unitno) header
            read (unitno) dummy
            close (unitno)
            maxv=maxval(dummy(:,:))
            write(6,*)'maxv on input2 for dust4 = ',maxv
            call flush(6)
            do ipn=1,nip
              do ivl=1,nvl
                tr3d(ivl,ipn,ibegin+p_dust_4) = dummy(ivl,ipn)
                trdp(ivl,ipn,ibegin+p_dust_4) = dummy(ivl,ipn)*dp3d(ivl,ipn)
              enddo
            enddo
            open (unit=unitno, file='dust5.in', form="unformatted", action='read')
            read (unitno) header
            read (unitno) dummy
            close (unitno)
            maxv=maxval(dummy(:,:))
            write(6,*)'maxv on input2 for dust5 = ',maxv
            call flush(6)
            do ipn=1,nip
              do ivl=1,nvl
                tr3d(ivl,ipn,ibegin+p_dust_5) = dummy(ivl,ipn)
                trdp(ivl,ipn,ibegin+p_dust_5) = dummy(ivl,ipn)*dp3d(ivl,ipn)
              enddo
            enddo
          endif   ! chem_opt = 300, dust2

          if( seas_opt == 1 )then
            open (unit=unitno, file='seas1.in', form="unformatted", action='read')
            read (unitno) header
            read (unitno) dummy
            close (unitno)
            maxv=maxval(dummy(:,:))
            write(6,*)'maxv on input2 for seas1 = ',maxv
            call flush(6)
            do ipn=1,nip
              do ivl=1,nvl
                tr3d(ivl,ipn,ibegin+p_seas_1) = dummy(ivl,ipn)
                trdp(ivl,ipn,ibegin+p_seas_1) = dummy(ivl,ipn)*dp3d(ivl,ipn)
              enddo
            enddo
            open (unit=unitno, file='seas2.in', form="unformatted", action='read')
            read (unitno) header
            read (unitno) dummy
            close (unitno)
            maxv=maxval(dummy(:,:))
            write(6,*)'maxv on input2 for seas2 = ',maxv
            call flush(6)
            do ipn=1,nip
              do ivl=1,nvl
                tr3d(ivl,ipn,ibegin+p_seas_2) = dummy(ivl,ipn)
                trdp(ivl,ipn,ibegin+p_seas_2) = dummy(ivl,ipn)*dp3d(ivl,ipn)
              enddo
            enddo
            open (unit=unitno, file='seas3.in', form="unformatted", action='read')
            read (unitno) header
            read (unitno) dummy
            close (unitno)
            maxv=maxval(dummy(:,:))
            write(6,*)'maxv on input2 for seas3 = ',maxv
            call flush(6)
            do ipn=1,nip
              do ivl=1,nvl
                tr3d(ivl,ipn,ibegin+p_seas_3) = dummy(ivl,ipn)
                trdp(ivl,ipn,ibegin+p_seas_3) = dummy(ivl,ipn)*dp3d(ivl,ipn)
              enddo
            enddo
            open (unit=unitno, file='seas4.in', form="unformatted", action='read')
            read (unitno) header
            read (unitno) dummy
            close (unitno)
            maxv=maxval(dummy(:,:))
            write(6,*)'maxv on input2 for seas4 = ',maxv
            call flush(6)
            do ipn=1,nip
              do ivl=1,nvl
                tr3d(ivl,ipn,ibegin+p_seas_4) = dummy(ivl,ipn)
                trdp(ivl,ipn,ibegin+p_seas_4) = dummy(ivl,ipn)*dp3d(ivl,ipn)
              enddo
            enddo
          endif
        endif   ! chem_opt >= 300 .and. chem_opt < 500
!
! tracer dispersion
!
        if( chem_opt == 500 ) then
          open (unit=unitno, file='tr1.in', form="unformatted", action='read')
          read (unitno) header
          read (unitno) dummy
          close (unitno)
          maxv=maxval(dummy(:,:))
          write(6,*)'maxv on input2 for tr1 = ',maxv
          call flush(6)
          do ipn=1,nip
            do ivl=1,nvl
              tr3d(ivl,ipn,ibegin+p_tr1) = dummy(ivl,ipn)
              trdp(ivl,ipn,ibegin+p_tr1) = dummy(ivl,ipn)*dp3d(ivl,ipn)
            enddo
          enddo
          open (unit=unitno, file='tr2.in', form="unformatted", action='read')
          read (unitno) header
          read (unitno) dummy
          close (unitno)
          maxv=maxval(dummy(:,:))
          write(6,*)'maxv on input2 for tr2 = ',maxv
          call flush(6)
          do ipn=1,nip
            do ivl=1,nvl
              tr3d(ivl,ipn,ibegin+p_tr2) = dummy(ivl,ipn)
              trdp(ivl,ipn,ibegin+p_tr2) = dummy(ivl,ipn)*dp3d(ivl,ipn)
            enddo
          enddo
        endif ! chem_opt=500
        deallocate(  dummy)            ! if chem_in_opt=1
      endif   ! if (chem_in_opt == 1)
    endif     ! if (chem_opt > 0)
!
! just like for volcanoes, later be used in chem_prep_fim
!
    if(chem_opt == 501 ) then
      write(6,*)'chem_opt = 500, itry to print emissions ',tr_mass,tr_height,nip
      do ipn=1,nip
        if(ipn .eq. 143724 ) then
          write(6,*)'emissions for point ',ipn
          emiss_tr_mass(ipn) = tr_mass
          emiss_tr_height(ipn) =  tr_height
          emiss_tr_dt(ipn) = 120000000
        endif
      enddo
    endif  ! if (chem_opt == 501)

    call returnunit (unitno)

!SMS$PARALLEL END

    ret = gptlstop ('chem_init')
    ret = gptlget_wallclock ('chem_init', 0, tchem_init)  ! The "0" is thread number
    print"(' CHEMISTRY INIT time:',F10.0)", tchem_init

    return

70  write(6,*)'chem_init: error opening unit=', unitno, '. Stopping'
    stop
90  write(6,*)'chem_init: error reading from unit=', unitno, '. Stopping'
    stop
  end subroutine chem_init
end module module_fim_chem_init

