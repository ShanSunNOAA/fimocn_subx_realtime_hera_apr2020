module module_wrf_output

contains

  subroutine wrf_output(its,nts,pr3d,tk3d,tr,phys2dwrf)

    use module_constants,only:rd
    use module_control,only:ArchvStep,dt,nip,ntra,ntrb,nvlp1
    use fimnamelist,   only:ArchvTimeUnit,nvl
    use icosio_wrapper, only: maybe_write
    use module_initial_chem_namelists
    use module_wrf_control,only:num_chem,num_moist,nvl_gocart,nvl_gocart

    integer,intent(in)::its,nts
!sms$distribute (dh,2) begin
    real,intent(in)::&
      pr3d(nvlp1,nip),&
      tk3d(nvl,nip),tr(nvl,nip,ntra+ntrb)
    real:: qcct(nvl,nip),qict(nvl,nip),qrct(nvl,nip),qsct(nvl,nip),&
      rho_phys(nvl,nip),sea1(nvl,nip),sea2(nvl,nip),sea3(nvl,nip),sea4(nvl,nip),&
      trco(nvl,nip)
!sms$distribute end
!sms$distribute(dh,1) begin
    real,intent(in)::phys2dwrf(:,:) ! (nip,:)
!sms$distribute end
    integer::ichem_start,imoist_start,j,k
    real::dpsum

    if (mod(its,ArchvStep)==0.or.(its==nts.and.ArchvTimeUnit.eq.'ts')) then

      ichem_start=ntra+1
      imoist_start=ntra

      if ((.not.mp_physics==0).or.(.not.cu_physics==0)) then
        trco(:,:) = 0.
!       trco(:,:) = tr(:,:,5)
        trco(1,:) = phys2dwrf(:,1)
        trco(2,:) = phys2dwrf(:,2)
        trco(3,:) = phys2dwrf(:,3)
        trco(4,:) = phys2dwrf(:,4)
        trco(5,:) = phys2dwrf(:,5)
        trco(6,:) = phys2dwrf(:,6)
      endif
      if (mp_physics == 2) then
        qcct(:,:) = tr(:,:,imoist_start+p_qc)
        qrct(:,:) = tr(:,:,imoist_start+p_qr)
        qict(:,:) = tr(:,:,imoist_start+p_qi)
        qsct(:,:) = tr(:,:,imoist_start+p_qs)
      endif

      if (mp_physics.eq.2) then
        call maybe_write(its,'qcct',qcct,nvl)
        call maybe_write(its,'qrct',qrct,nvl)
        call maybe_write(its,'qict',qict,nvl)
        call maybe_write(its,'qsct',qsct,nvl)
      endif
      if ((mp_physics.ne.0).or.(cu_physics.ne.0)) then
        call maybe_write(its,'trco',trco,nvl)
      endif
    endif

  end subroutine wrf_output

end module module_wrf_output
