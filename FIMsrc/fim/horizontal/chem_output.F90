module module_chem_output

contains

  subroutine chem_output(its,nts,aod2d,exch,p10,pm25,pr3d,tk3d,tr,trfall,&
    phys2dwrf,tr1_tavg)

    use module_constants,only:rd
    use module_control,  only:ArchvStep,dt,nip,ntra,ntrb,nvlp1
    use fimnamelist,     only:ArchvTimeUnit,nvl
    use icosio_wrapper,  only: maybe_write
    use module_initial_chem_namelists
    use module_wrf_control,only:num_chem,num_moist,nvl_gocart,nvl_gocart

    integer,intent(in)::its,nts
!sms$distribute (dh,1) begin
    real,intent(in)::aod2d(nip),trfall(nip,num_chem)
    real::intaer(nip),intash(nip),intbc(nip),intdust(nip),intoc(nip),&
      intsulf(nip)
    real,intent(in)::phys2dwrf(:,:) ! (nip,:)
!sms$distribute end
!sms$distribute(dh,2) begin
    real,intent(inout):: tr1_tavg(nvl,nip)
    real,intent(in)::exch(nvl,nip),p10(nvl,nip),pm25(nvl,nip),pr3d(nvlp1,nip),tk3d(nvl,nip),&
      tr(nvl,nip,ntra+ntrb)
    real::d1st(nvl,nip),d2st(nvl,nip),d3st(nvl,nip),d4st(nvl,nip),d5st(nvl,nip),&
      dms1(nvl,nip),qcct(nvl,nip),qict(nvl,nip),qrct(nvl,nip),qsct(nvl,nip),&
      rho_phys(nvl,nip),sea1(nvl,nip),sea2(nvl,nip),sea3(nvl,nip),sea4(nvl,nip),&
      trco(nvl,nip)
!sms$distribute end
    integer::ichem_start,imoist_start,j,k
    real::dpsum

    if (mod(its,ArchvStep)==0.or.(its==nts.and.ArchvTimeUnit.eq.'ts')) then

      ichem_start=ntra
!     imoist_start=4-1

      if ((.not.mp_physics==0).or.(.not.cu_physics==0)) then
        trco(:,:) = tr(:,:,5)
        trco(1,:) = phys2dwrf(:,1)
        trco(2,:) = phys2dwrf(:,2)
        trco(3,:) = phys2dwrf(:,3)
        trco(4,:) = phys2dwrf(:,4)
        trco(5,:) = phys2dwrf(:,5)
        trco(6,:) = phys2dwrf(:,6)
      endif
      if (chem_opt == 300) then
        d3st(:,:) = tr(:,:,ichem_start+p_dust_3)
        d4st(:,:) = tr(:,:,ichem_start+p_dust_4)
        d5st(:,:) = tr(:,:,ichem_start+p_dust_5)
        sea3(:,:) = tr(:,:,ichem_start+p_seas_3)
        sea4(:,:) = tr(:,:,ichem_start+p_seas_4)
      endif
      if (chem_opt == 500) then
        d1st(:,:) = tr(:,:,ichem_start+p_tr1)
        d2st(:,:) = tr(:,:,ichem_start+p_tr2)
      endif
      if (chem_opt >= 300 .and. chem_opt < 500) then
        dms1(:,:) = tr(:,:,ichem_start+p_dms)
        d1st(:,:) = tr(:,:,ichem_start+p_dust_1)
        d2st(:,:) = tr(:,:,ichem_start+p_dust_2)
        sea1(:,:) = tr(:,:,ichem_start+p_seas_1)
        sea2(:,:) = tr(:,:,ichem_start+p_seas_2)
      endif
!     if (mp_physics == 2) then
!       qcct(:,:) = tr(:,:,imoist_start+p_qc)
!       qrct(:,:) = tr(:,:,imoist_start+p_qr)
!       qict(:,:) = tr(:,:,imoist_start+p_qi)
!       qsct(:,:) = tr(:,:,imoist_start+p_qs)
!     endif

      if (chem_opt >= 300 .and. chem_opt < 500) then
!SMS$PARALLEL(dh, j) BEGIN
        do j=1,nip
          dpsum=0.
          intash(j)=0.
          intaer(j)=0.
          intbc(j)=0.
          intoc(j)=0.
          intsulf(j)=0.
          intdust(j)=0.
          do k=1,nvl
            dpsum=dpsum+(pr3d(k,j)-pr3d(k+1,j))
            rho_phys(k,j)=.5*(pr3d(k,j)+pr3d(k+1,j))&
              /(RD*tk3d(k,j)) !*(1.+.608*qv3d(k,j))
            intaer(j)=intaer(j)+tr(k,j,ichem_start+p_p25)*(pr3d(k,j)-pr3d(k+1,j))*rho_phys(k,j)
            intbc(j)=intbc(j)+(tr(k,j,ichem_start+p_bc1)&
              +tr(k,j,ichem_start+p_bc2))*(pr3d(k,j)-pr3d(k+1,j))*rho_phys(k,j)
            intoc(j)=intoc(j)+(tr(k,j,ichem_start+p_oc1)&
              +tr(k,j,ichem_start+p_oc2))*(pr3d(k,j)-pr3d(k+1,j))*rho_phys(k,j)
            intsulf(j)=intsulf(j)+tr(k,j,ichem_start+p_sulf)*(pr3d(k,j)&
              -pr3d(k+1,j))*rho_phys(k,j)
            if (chem_opt == 300)intdust(j)=intdust(j)+(d1st(k,j)&
              +.286*d2st(k,j))*(pr3d(k,j)-pr3d(k+1,j))*rho_phys(k,j)
            if (chem_opt==304.or.chem_opt==316.or.chem_opt==317) then
              intdust(j)=intdust(j)+d1st(k,j)*(pr3d(k,j)-pr3d(k+1,j))&
                *rho_phys(k,j)
            endif
            if (chem_opt == 316) then
              intash(j)=intash(j)+(tr(k,j,ichem_start+p_vash_1)       &
                + tr(k,j,ichem_start+p_vash_2)       &
                + tr(k,j,ichem_start+p_vash_3)       &
                + tr(k,j,ichem_start+p_vash_4)       &
                + tr(k,j,ichem_start+p_vash_5)       &
                + tr(k,j,ichem_start+p_vash_6)       &
                + tr(k,j,ichem_start+p_vash_7)       &
                + tr(k,j,ichem_start+p_vash_8)       &
                + tr(k,j,ichem_start+p_vash_9)       &
                + tr(k,j,ichem_start+p_vash_10))     &
                *(pr3d(k,j)-pr3d(k+1,j))*rho_phys(k,j)
            endif
            if (chem_opt == 317 ) then
              intash(j)=intash(j)+(tr(k,j,ichem_start+p_vash_1)       &
                + tr(k,j,ichem_start+p_vash_2)       &
                + tr(k,j,ichem_start+p_vash_3)       &
                + tr(k,j,ichem_start+p_vash_4))      &
                *(pr3d(k,j)-pr3d(k+1,j))*rho_phys(k,j)
            endif
          enddo
          if (chem_opt == 316 .or. chem_opt == 317 ) intash(j)=intash(j)/dpsum
          intaer(j)=intaer(j)/dpsum
          intbc(j)=intbc(j)/dpsum
          intoc(j)=intoc(j)/dpsum
          intsulf(j)=intsulf(j)/dpsum
          intdust(j)=intdust(j)/dpsum
        enddo
!SMS$PARALLEL END
      endif ! chem_opt >= 300 .and. chem_opt < 500
      if (chem_opt == 502) then
!SMS$PARALLEL(dh, j) BEGIN
        do j=1,nip
          dpsum=0.
          intash(j)=0.
            do k=1,nvl
              dpsum=dpsum+(pr3d(k,j)-pr3d(k+1,j))
              rho_phys(k,j)=.5*(pr3d(k,j)+pr3d(k+1,j))&
              /(RD*tk3d(k,j)) !*(1.+.608*qv3d(k,j))
              intash(j)=intash(j)+(tr(k,j,ichem_start+p_vash_1)       &
                + tr(k,j,ichem_start+p_vash_2)       &
                + tr(k,j,ichem_start+p_vash_3)       &
                + tr(k,j,ichem_start+p_vash_4))      &
                *(pr3d(k,j)-pr3d(k+1,j))*rho_phys(k,j)
          enddo
          intash(j)=intash(j)/dpsum
        enddo
!SMS$PARALLEL END
      endif ! chem_opt = 502

!     if (mp_physics.eq.2) then
!       call maybe_write(its,'qcct',qcct,nvl)
!       call maybe_write(its,'qrct',qrct,nvl)
!       call maybe_write(its,'qict',qict,nvl)
!       call maybe_write(its,'qsct',qsct,nvl)
!     endif
      if ((mp_physics.ne.0).or.(cu_physics.ne.0)) then
        call maybe_write(its,'trco',trco,nvl)
      endif
      if (chem_opt == 500) then
        if(its.gt.1)tr1_tavg(:,:) = tr1_tavg(:,:)/float(ArchvStep-1)
        call maybe_write(its,'c13D',d1st,nvl)
        call maybe_write(its,'c23D',tr1_tavg,nvl)
        tr1_tavg(:,:) = 0.
      endif 
      if (chem_opt.ge.300 .and. chem_opt.lt.500) then
        call maybe_write(its,'ex3D',exch,nvl)
        call maybe_write(its,'pm25',pm25,nvl)
        call maybe_write(its,'pm10',p10,nvl)
        call maybe_write(its,'dms1',dms1,nvl)
        call maybe_write(its,'d1st',d1st,nvl)
        call maybe_write(its,'d2st',d2st,nvl)
        call maybe_write(its,'s1ea',sea1,nvl)
        call maybe_write(its,'s2ea',sea2,nvl)
        if (chem_opt.eq.300) then
          call maybe_write(its,'s3ea',sea3,nvl)
          call maybe_write(its,'s4ea',sea4,nvl)
          call maybe_write(its,'d3st',d3st,nvl)
          call maybe_write(its,'d4st',d4st,nvl)
          call maybe_write(its,'d5st',d5st,nvl)
        endif

        dms1(:,:) = tr(:,:,ichem_start+p_bc1)
        call maybe_write(its,'pbc1',dms1,nvl)
        dms1(:,:) = tr(:,:,ichem_start+p_bc2)
        call maybe_write(its,'pbc2',dms1,nvl)
        dms1(:,:) = tr(:,:,ichem_start+p_oc1)
        call maybe_write(its,'obc1',dms1,nvl)
        dms1(:,:) = tr(:,:,ichem_start+p_oc2)
        call maybe_write(its,'obc2',dms1,nvl)
        dms1(:,:) = tr(:,:,ichem_start+p_sulf)
        call maybe_write(its,'sulf',dms1,nvl)
        dms1(:,:) = tr(:,:,ichem_start+p_so2)
        call maybe_write(its,'pso2',dms1,nvl)
        dms1(:,:) = tr(:,:,ichem_start+p_msa)
        call maybe_write(its,'pmsa',dms1,nvl)
        dms1(:,:) = tr(:,:,ichem_start+p_p25)
        call maybe_write(its,'pp25',dms1,nvl)
        dms1(:,:) = tr(:,:,ichem_start+p_p10)
        call maybe_write(its,'pp10',dms1,nvl)
        if (chem_opt.eq.316.or.chem_opt.eq.317) then
          print *,'p_vash_1,p_vash_4 = ',p_vash_1,p_vash_4
          dms1(:,:) = tr(:,:,ichem_start+p_vash_1)
          call maybe_write(its,'ash1',dms1,nvl)
          dms1(:,:) = tr(:,:,ichem_start+p_vash_2)
          call maybe_write(its,'ash2',dms1,nvl)
          dms1(:,:) = tr(:,:,ichem_start+p_vash_3)
          call maybe_write(its,'ash3',dms1,nvl)
          dms1(:,:) = tr(:,:,ichem_start+p_vash_4)
          call maybe_write(its,'ash4',dms1,nvl)
          if (chem_opt.eq.316) then
            dms1(:,:) = tr(:,:,ichem_start+p_vash_5)
            call maybe_write(its,'ash5',dms1,nvl)
            dms1(:,:) = tr(:,:,ichem_start+p_vash_6)
            call maybe_write(its,'ash6',dms1,nvl)
            dms1(:,:) = tr(:,:,ichem_start+p_vash_7)
            call maybe_write(its,'ash7',dms1,nvl)
            dms1(:,:) = tr(:,:,ichem_start+p_vash_8)
            call maybe_write(its,'ash8',dms1,nvl)
            dms1(:,:) = tr(:,:,ichem_start+p_vash_9)
            call maybe_write(its,'ash9',dms1,nvl)
            dms1(:,:) = tr(:,:,ichem_start+p_vash_10)
            call maybe_write(its,'ash0',dms1,nvl)
          endif !chem_opt=316
        endif !chem_opt=316 or chem_opt=317
      endif !chem_opt.ge.300 .and. chem_opt.lt.500
      if (chem_opt == 500) then
        call maybe_write(its,'fl2D',intaer,1, twodfile=.true.)
      endif
! output for volcanic ash only
      if (chem_opt.eq.502) then
          call maybe_write(its,'iash',intash,1, twodfile=.true.)
          print *,'p_vash_1,p_vash_4 = ',p_vash_1,p_vash_4
          dms1(:,:) = tr(:,:,ichem_start+p_vash_1)
          call maybe_write(its,'ash1',dms1,nvl)
          dms1(:,:) = tr(:,:,ichem_start+p_vash_2)
          call maybe_write(its,'ash2',dms1,nvl)
          dms1(:,:) = tr(:,:,ichem_start+p_vash_3)
          call maybe_write(its,'ash3',dms1,nvl)
          dms1(:,:) = tr(:,:,ichem_start+p_vash_4)
          call maybe_write(its,'ash4',dms1,nvl)
      endif
      if (chem_opt.ge.300 .and. chem_opt.lt.500) then
        call maybe_write(its,'ia2D',intaer,1, twodfile=.true.)
        call maybe_write(its,'ib2D',intbc,1, twodfile=.true.)
        call maybe_write(its,'io2D',intoc,1, twodfile=.true.)
        call maybe_write(its,'is2D',intsulf,1, twodfile=.true.)
        call maybe_write(its,'id2D',intdust,1, twodfile=.true.)
        call maybe_write(its,'ao2D',aod2d,1, twodfile=.true.)
        if (chem_opt.eq.316.or.chem_opt.eq.317) then
          call maybe_write(its,'iash',intash,1, twodfile=.true.)
        endif
      endif !chem_opt.ge.300 .and. chem_opt.lt.500
    endif

  end subroutine chem_output

end module module_chem_output
