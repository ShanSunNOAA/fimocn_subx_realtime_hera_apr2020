module hycom_dissip
!*********************************************************************
!    del4prep, dissip 
!       Routines for biharmonic momentum dissipation
!       R. Bleck                                        Jul 2013
!*********************************************************************
contains

  subroutine del4prep (nstep, u_vel, v_vel)
    use module_control,   only: nip
    use module_constants, only: nprox, prox, cs, sn, perm
    use hycom_constants,  only: wet,slip
    use hycom_control,    only: kgap
    use fimnamelist,      only: kdm,diag_intvl,itest
    use stencilprint,     only: stencl

! --- compute a quantity proportional to the (negative) laplacian by
! --- subtracting from each -u,v- value the average of its 5 or 6 neighbors.

    implicit none
! Arguments
    integer, intent(in) :: nstep
!SMS$DISTRIBUTE (dh,2) BEGIN
    real, intent(inout) :: u_vel(kdm,nip)	! field to be diffused
    real, intent(inout) :: v_vel(kdm,nip)	! field to be diffused
! Local variables:
    real    :: worku(kdm,nip)
    real    :: workv(kdm,nip)
!SMS$DISTRIBUTE END
    integer :: k	        ! layer index
    integer :: i		! icos point index
    integer :: ix		! neighbor across shared edge
    integer :: edg		! icos edge index
    real    :: factor
    real    :: uedge,vedge,ungbor,vngbor,avgu(kdm),avgv(kdm)
    character :: text*32
    logical :: vrbos
    real,parameter :: bihar=1.	! 0 for laplacian, 1 for biharmonic mom.mixing

!SMS$EXCHANGE(u_vel,v_vel) 
!SMS$PARALLEL (dh,i) BEGIN
!SMS$HALO_COMP(<1,1>) BEGIN
!$OMP PARALLEL DO
    do i=1,nip
     worku(:,i)=u_vel(:,i)
     workv(:,i)=v_vel(:,i)
    end do
!$OMP END PARALLEL DO
!SMS$HALO_COMP END

!$OMP PARALLEL DO PRIVATE (vrbos,avgu,avgv,edg,ix,k,uedge,vedge,ungbor,vngbor,factor)
    do i=1,nip
     vrbos=i == itest .and. mod(nstep,diag_intvl) == 0
     if (wet(i) > 0) then
      avgu(:)=0.
      avgv(:)=0.
      do edg=1,nprox(i)				! loop through edges
       ix=prox(edg,i)				! neighbor across shared edge
       do k=1,kdm				! verical loop

! --- rotate neighboring vector into local cell-centered coordinates.
! --- the available -sn,cs- transformation coefficients force us to make this
! --- a 2-step process: first, transform halfway to edge, then to center

        if (wet(ix) > 0) then
         uedge= cs(2,edg,i)*worku(k,ix)+sn(2,edg,i)*workv(k,ix)
         vedge=-sn(2,edg,i)*worku(k,ix)+cs(2,edg,i)*workv(k,ix)
         ungbor=cs(1,edg,i)*uedge-sn(1,edg,i)*vedge
         vngbor=sn(1,edg,i)*uedge+cs(1,edg,i)*vedge
        else				! sidewall
         ungbor=slip*worku(k,i)
         vngbor=slip*workv(k,i)
        end if
        avgu(k)=avgu(k)+ungbor
        avgv(k)=avgv(k)+vngbor

        if (vrbos .and. mod(k,kgap).eq.1 .and. wet(ix) > 0) then
         print '(2i8,2i3,a,3p,2f7.2,2x,2f7.2)',nstep,perm(i),k,edg,	&
           ' (del4prep) u,v orig/rot:',worku(k,ix),workv(k,ix),		&
           ungbor,vngbor
        end if

       end do				! vertical loop
      end do				! loop through edges

      do k=1,kdm
       factor=bihar/nprox(i)		! bihar=1 for full biharmonic mixing

       if (vrbos .and. mod(k,kgap).eq.1) then
        print '(2i8,i3,a,3p,2f7.2,2x,2f7.2)',nstep,perm(i),k,		&
          ' (del4prep) u,v,del^2 u/v:',u_vel(k,i),v_vel(k,i),		&
          u_vel(k,i)-avgu(k)*factor,v_vel(k,i)-avgv(k)*factor
       end if

       u_vel(k,i)=u_vel(k,i)-avgu(k)*factor
       v_vel(k,i)=v_vel(k,i)-avgv(k)*factor
      end do
     end if
    end do				! horizontal loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END

    if (mod(nstep,diag_intvl).eq.0) then
     write (text,'(a,i8)') '(ocn_del4prep) step',nstep
     call stencl(u_vel,kdm,1.e3,trim(text)//'  -del^2 u (mm/sec)')
     call stencl(v_vel,kdm,1.e3,trim(text)//'  -del^2 v (mm/sec)')
    end if

    return
  end subroutine del4prep


  subroutine dissip (nstep, leap, worku, workv, dp)

! dissipation of -u,v- based on del^2(u), del^2(v) (thickness-weighted).
! to convert to biharmonic dissipation [i.e. switch to del^4(u), del^4(v)],
! precede this call by a call to del4prep which subtracts from each -u,v-
! value the average of the 5 or 6 surrounding values.

    use module_control,   only: nip,npp
    use module_constants, only: nprox, prox, rarea, sideln, cs, sn,	&
                                perm, nedge
    use hycom_constants,  only: wet,slip,land_spval
    use hycom_control,    only: kgap
    use fimnamelist,      only: kdm,diag_intvl,itest
    use stencilprint,     only: stencl
    use hycom_edgvar,     only: blokf

    implicit none
    integer,intent(in) :: nstep,leap
!SMS$DISTRIBUTE (dh,2) BEGIN
! Arguments
    real, intent(inout) :: worku(kdm,nip)	! on input: field to be diffused
    real, intent(inout) :: workv(kdm,nip)	! on input: field to be diffused
    real, intent(in)    :: dp   (kdm,nip,2)	! layer thickness
! Local variables
    real    :: pres(kdm+1,nip)
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE (dh,3) BEGIN
    real    :: ufxrot(kdm,npp,nip)	! u momentum flux across edges
    real    :: vfxrot(kdm,npp,nip)	! v momentum flux across edges
!SMS$DISTRIBUTE END
    integer :: k		! layer index
    integer :: i		! icos point index
    integer :: ix		! neighbor across shared edge
    integer :: edg		! icos edge index
    real    :: urot,vrot,ungbor,vngbor,uflux,vflux,factor
    real    :: blokd
    character :: string*4,text*32
    logical:: vrbos
    real, parameter:: thshld=1.e-11

! Statement function
    real      :: hfharm,a,b	! 0.5 * harmonic average
    hfharm(a,b)=a*b/(a+b)	! (see Appx.D, 1992 MICOM paper)

! do k=1,kdm,kgap
!  write (string,'(a,i2)') 'k=',k
!  call findmxmn2(worku,kdm,nip,k,'(dissip) u-in '//string,wet)
!  call findmxmn2(workv,kdm,nip,k,'(dissip) v-in '//string,wet)
! end do

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (k)
    do i=1,nip
     pres(:,i)=land_spval
     if (wet(i).gt.0) then
      pres(1,i)=0.
      do k=1,kdm
       pres(k+1,i)=pres(k,i)+dp(k,i,leap)
      end do
     end if
    end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END

!SMS$EXCHANGE(worku,workv,dp,pres)

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (vrbos,edg,ix,k,urot,vrot,ungbor,vngbor,blokd,factor)
    do i=1,nip
     vrbos=i == itest .and. mod(nstep,diag_intvl) == 0
     if (wet(i) > 0) then
      do edg=1,nedge(i)
       ix=prox(edg,i)
       do k=1,kdm

! --- Transform u,v at neighboring icos pt to coord.system centered on edge.
! --- cs and sn are coordinate transformation constants.

        urot= cs(1,edg,i)*worku(k,i)+sn(1,edg,i)*workv(k,i)
        vrot=-sn(1,edg,i)*worku(k,i)+cs(1,edg,i)*workv(k,i)

        if (wet(ix) > 0) then
         ungbor= cs(2,edg,i)*worku(k,ix)+sn(2,edg,i)*workv(k,ix)
         vngbor=-sn(2,edg,i)*worku(k,ix)+cs(2,edg,i)*workv(k,ix)

! --- modify neighboring u/v to account for possible presence of sidewall.
         blokd=blokf(dp(k,i,leap),pres(k+1,i),pres(kdm+1,ix),pres(k,ix))
         ungbor=(1.-blokd)*ungbor+blokd*urot*slip
         vngbor=(1.-blokd)*vngbor+blokd*vrot*slip

        else			! neighbor is land point
         ungbor=urot*slip
         vngbor=vrot*slip
        end if			! neighbor is ocean/land point

! --- momentum fluxes (pos.inward) in edge-centric coord.system:
        ufxrot(k,edg,i)=(ungbor-urot)*sideln(edg,i)		&
           *2.*hfharm(max(dp(k,i ,leap),thshld),		&
                      max(dp(k,ix,leap),thshld))		! N/sec
        vfxrot(k,edg,i)=(vngbor-vrot)*sideln(edg,i)		&
           *2.*hfharm(max(dp(k,i ,leap),thshld),		&
                      max(dp(k,ix,leap),thshld))		! N/sec

        if (vrbos .and. mod(k,kgap) == 1) then
         print 101,'orig u,v at',perm(i),perm(ix),k,edg, 	&
           worku(k,i),workv(k,i),worku(k,ix),workv(k,ix)
         print 101,' rot u,v at',perm(i),perm(ix),k,edg,   	&
           urot,vrot,ungbor,vngbor
101      format ('(ocn dissip) ',a,2i7,2i3,3p,3(f10.2,f8.2))
         factor=rarea(i)/max(dp(k,i,leap),thshld)
         print 102,' u/vflx in rotated system',perm(i),		&
           perm(ix),k,edg,ufxrot(k,edg,i)*factor,		&
           vfxrot(k,edg,i)*factor
102      format ('(ocn dissip) ',a,2i7,2i3,2es11.3)
        end if
       end do				! vertical loop
      end do				! loop through edges
     end if
    end do				! horizontal loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (edg,ix,k,factor,uflux,vflux)
    do i=1,nip
     if (wet(i).gt.0) then
      worku(:,i)=0.			! on output: dissip.tdcy (1/sec)
      workv(:,i)=0.			! on output: dissip.tdcy (1/sec)
      do edg=1,nprox(i)
       ix=prox(edg,i)			! neighbor across shared edge
       do k=1,kdm			! vertical loop
        factor=rarea(i)/max(dp(k,i,leap),thshld)
! --- rotate momentum fluxes back to lat/lon coord.system
        uflux=cs(1,edg,i)*ufxrot(k,edg,i)			&
             -sn(1,edg,i)*vfxrot(k,edg,i)
        vflux=sn(1,edg,i)*ufxrot(k,edg,i)			&
             +cs(1,edg,i)*vfxrot(k,edg,i)
        worku(k,i)=worku(k,i)+uflux*factor			! 1/sec
        workv(k,i)=workv(k,i)+vflux*factor			! 1/sec
       end do				! vertical loop
      end do				! loop through edges
     end if				! ocean point
    end do				! horizontal loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END

    if (mod(nstep,diag_intvl).eq.0) then
     write (text,'(a,i8)') '(ocn_dissip)  step',nstep
     call stencl(worku,kdm,1.e8,trim(text)//' -u- tdcy (1/sec) x 10^8')
     call stencl(workv,kdm,1.e8,trim(text)//' -v- tdcy (1/sec) x 10^8')
    end if

! do k=1,kdm,kgap
!  write (string,'(a,i2)') 'k=',k
!  call findmxmn2(worku,kdm,nip,k,'(dissip) u-out '//string,wet)
!  call findmxmn2(workv,kdm,nip,k,'(dissip) v-out '//string,wet)
! end do

    return
  end subroutine dissip
end module hycom_dissip
