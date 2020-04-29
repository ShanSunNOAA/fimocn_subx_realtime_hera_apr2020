!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutines in this module are computational only. They contain no SMS directives
! except for IGNORE. Thus they are good candidates for migration to MIC or GPU.
!
! momtum_compute: contains all the computational routines for momentum calculations.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module momtum_compute
  implicit none

  private
  public :: del4prep, dissip, momtum1

contains

  subroutine del4prep (u_vel, v_vel)
    use fimnamelist,      only: nvl, PrintIpnDiag
    use module_constants, only: nprox, prox, cs, sn, biharw
    use global_bounds,    only: ims, ime, ips, ipe, ihe
    use stencilprint,     only: stencl

! --- subtract from each -u,v- value the average of its 5 or 6 neighbors.
! --- blend laplacian and biharmonic mixing by multiplying neighboring
! --- values by a factor between 0 and 1 before subtracting.

! Arguments
    real, intent(inout) :: u_vel(nvl,ims:ime)	! field to be diffused
    real, intent(inout) :: v_vel(nvl,ims:ime)	! field to be diffused

! Local variables:
    integer :: k	        ! layer index
    integer :: ipn		! icos point index
    integer :: ipx		! neighbor across shared edge
    integer :: edg		! icos edge index
    real    :: worku(nvl,ims:ime)
    real    :: workv(nvl,ims:ime)
    real    :: factor
    real    :: uxy2
    real    :: vxy2
    real    :: avgu(nvl)
    real    :: avgv(nvl)

! Halo comp so horizontal loop indices are ips,ihe
!sms$ignore begin
!$OMP PARALLEL DO PRIVATE (k) SCHEDULE (static)
    do ipn=ips,ihe
      do k=1,nvl
        worku(k,ipn) = u_vel(k,ipn)
        workv(k,ipn) = v_vel(k,ipn)
      end do
    end do
!$OMP END PARALLEL DO
!sms$ignore end

!sms$ignore begin
!$OMP PARALLEL DO PRIVATE (k,edg,ipx,uxy2,vxy2,avgu,avgv,factor) SCHEDULE (runtime)
    do ipn=ips,ipe
      do k=1,nvl
        avgu(k) = 0.
        avgv(k) = 0.
      end do

      do edg=1,nprox(ipn)			! loop through edges
        ipx = prox(edg,ipn)			! neighbor across shared edge
        do k=1,nvl				! verical loop

! --- rotate neighboring vector into local cell-centered coordinates.
! --- the available -sn,cs- transformation coefficients force us to make this
! --- a 2-step process: first, transform halfway to edge, then to center

          uxy2 = cs(2,edg,ipn)*worku(k,ipx) + sn(2,edg,ipn)*workv(k,ipx)
          vxy2 =-sn(2,edg,ipn)*worku(k,ipx) + cs(2,edg,ipn)*workv(k,ipx)

          avgu(k) = avgu(k) + cs(1,edg,ipn)*uxy2 - sn(1,edg,ipn)*vxy2
          avgv(k) = avgv(k) + sn(1,edg,ipn)*uxy2 + cs(1,edg,ipn)*vxy2
        end do				! vertical loop
      end do				! loop through edges

      do k=1,nvl
        factor = biharw(k)/nprox(ipn)
        u_vel(k,ipn) = u_vel(k,ipn) - avgu(k)*factor
        v_vel(k,ipn) = v_vel(k,ipn) - avgv(k)*factor
      end do
    end do				! horizontal loop
!$OMP END PARALLEL DO
!sms$ignore end

    if (PrintIpnDiag > 0) then
      call stencl(u_vel,nvl,100.,'(atm del4prep) -u- output x 100')
      call stencl(v_vel,nvl,100.,'(atm del4prep) -v- output x 100')
    end if

    return

  end subroutine del4prep

! dissip: dissipation of -u,v- based on del^2 u, del^2 v (thickness-weighted).
! to substitute biharmonic for laplacian dissipation, precede this call
! by a call to del4prep which subtracts the average of the 5 or 6
! surrounding values from each -u,v- value.

  subroutine dissip (worku, workv, delp)
    use module_control,   only: nip, npp
    use fimnamelist,      only: nvl, PrintIpnDiag
    use module_constants, only: nprox, prox, rarea, sideln, cs, sn,	&
                                perm, nedge, dfflen
    use stencilprint,     only: stencl
    use global_bounds,    only: ims, ime, ips, ipe
    use findmaxmin2

! Arguments

    real, intent(inout) :: worku(nvl,ims:ime)	! field to be diffused
    real, intent(inout) :: workv(nvl,ims:ime)	! field to be diffused
    real, intent(in)    :: delp (nvl,ims:ime)	! layer thickness

! Local variables
    real    :: ufxrot(nvl,npp,ims:ime)	! u momentum flux across edges
    real    :: vfxrot(nvl,npp,ims:ime)	! v momentum flux across edges
    integer :: k		! layer index
    integer :: ipn		! icos point index
    integer :: ipx		! neighbor across shared edge
    integer :: edg		! icos edge index
    real    :: factor
    real    :: uxy1
    real    :: uxy2
    real    :: vxy1
    real    :: vxy2
    real    :: uflux
    real    :: vflux
    character :: string*4
    logical :: vrbos
    real, parameter :: thshld = 1.e-11
    real      :: hfharm,a,b	! 0.5 * harmonic average

! Statement funcion
    hfharm(a,b) = a*b/(a+b)	! (see Appx.D, 1992 MICOM paper)

!sms$ignore begin

#ifdef DEBUGPRINT
  do k=1,nvl,7
   write (string,'(a,i2)') 'k=',k
   call findmxmn2(worku,nvl,nip,k,'(dissip) u-in '//string)
   call findmxmn2(workv,nvl,nip,k,'(dissip) v-in '//string)
  end do
#endif

!$OMP PARALLEL DO PRIVATE (vrbos,edg,ipx,k,uxy1,vxy1,uxy2,vxy2,factor) SCHEDULE (runtime)
    do ipn=ips,ipe
      vrbos = ipn == PrintIpnDiag
      do edg=1,nedge(ipn)
        ipx = prox(edg,ipn)
        do k=1,nvl

! --- Transform u,v at neighboring icos pt to coord.system centered on edge.
! --- cs and sn are coordinate transformation constants.
! --- uxy,vxy are values of u and v rotated into local system.

          uxy1 = cs(1,edg,ipn)*worku(k,ipn) + sn(1,edg,ipn)*workv(k,ipn)
          vxy1 =-sn(1,edg,ipn)*worku(k,ipn) + cs(1,edg,ipn)*workv(k,ipn)
          uxy2 = cs(2,edg,ipn)*worku(k,ipx) + sn(2,edg,ipn)*workv(k,ipx)
          vxy2 =-sn(2,edg,ipn)*worku(k,ipx) + cs(2,edg,ipn)*workv(k,ipx)

! --- momentum fluxes (pos.inward) in local (rotated) coord.system:
          ufxrot(k,edg,ipn) = (uxy2-uxy1)*dfflen(k)			&
               *sideln(edg,ipn)*2.*hfharm(max(delp(k,ipn),thshld)	&
               ,max(delp(k,ipx),thshld))
          vfxrot(k,edg,ipn) = (vxy2-vxy1)*dfflen(k)			&
               *sideln(edg,ipn)*2.*hfharm(max(delp(k,ipn),thshld)	&
               ,max(delp(k,ipx),thshld))
#ifdef DEBUGPRINT
          if (vrbos .and. mod(k,7) == 1) then
            print 101,'orig u,v at',perm(ipn),perm(ipx),k,edg, 		&
                 worku(k,ipn),workv(k,ipn),worku(k,ipx),workv(k,ipx)
            print 101,' rot u,v at',perm(ipn),perm(ipx),k,edg,   	&
                 uxy1,vxy1,uxy2,vxy2
101         format ('(dissip) ',a,2i7,2i3,3(f10.2,f8.2))
            factor = rarea(ipn)/max(delp(k,ipn),thshld)
            print 102,' u/vflx in rotated system',perm(ipn),		&
                 perm(ipx),k,edg,ufxrot(k,edg,ipn)*factor,		&
                 vfxrot(k,edg,ipn)*factor
102         format ('(dissip) ',a,2i7,2i3,2f7.2)
          end if
#endif
        end do				! vertical loop
      end do				! loop through edges
    end do				! horizontal loop
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE (edg,ipx,k,factor,uflux,vflux) SCHEDULE (runtime)
    do ipn=ips,ipe
      worku(:,ipn) = 0.
      workv(:,ipn) = 0.
      do edg=1,nprox(ipn)
        ipx = prox(edg,ipn)		! neighbor across shared edge
        do k=1,nvl			! vertical loop
          factor = rarea(ipn)/max(delp(k,ipn),thshld)
! --- rotate momentum fluxes back to lat/lon coord.system
          uflux = cs(1,edg,ipn)*ufxrot(k,edg,ipn)			&
                - sn(1,edg,ipn)*vfxrot(k,edg,ipn)
          vflux = sn(1,edg,ipn)*ufxrot(k,edg,ipn)			&
                + cs(1,edg,ipn)*vfxrot(k,edg,ipn)
          worku(k,ipn) = worku(k,ipn) + uflux*factor
          workv(k,ipn) = workv(k,ipn) + vflux*factor
        end do				! vertical loop
      end do				! loop through edges
    end do				! horizontal loop
!$OMP END PARALLEL DO

    if (PrintIpnDiag > 0) then
      call stencl(worku,nvl,100.,'(atm dissip) -u- increment x 100')
      call stencl(workv,nvl,100.,'(atm dissip) -v- increment x 100')
    end if

#ifdef DEBUGPRINT
  do k=1,nvl,7
   write (string,'(a,i2)') 'k=',k
   call findmxmn2(worku,nvl,nip,k,'(dissip) u-out '//string)
   call findmxmn2(workv,nvl,nip,k,'(dissip) v-out '//string)
  end do
#endif


    return

!sms$ignore end

  end subroutine dissip

  subroutine momtum1 (its, u_velo, v_velo, exner, relvor,       	&
                      u_edg, v_edg, trcr_edg, geop_edg, mont_edg,	&
                      uvsq_edg, u_tndcy, v_tndcy, dp3d)

    use module_control,   only: nip, nvlp1, npp, nabl, ntra, dt
    use fimnamelist,      only: nvl,rleigh,rlthresh,rltaper,		&
                                PrintIpnDiag,janjic
    use module_constants, only: nprox, rarea, sidevec_c, sidevec_e,	&
                                corio, perm
    use module_variables, only: nf, of, vof, adbash1, adbash2, adbash3
    use global_bounds,    only: ims, ime, ips, ipe
    use stencilprint,     only: stencl
    use findmaxmin3

! Arguments
    integer, intent(in) :: its			! model time step
    real, intent(in)    :: u_edg   (nvl,npp,ims:ime)
    real, intent(in)    :: v_edg   (nvl,npp,ims:ime)
    real, intent(inout) :: u_velo  (nvl,ims:ime)
    real, intent(inout) :: v_velo  (nvl,ims:ime)
    real, intent(in)    :: exner   (nvlp1,ims:ime)
    real, intent(in)    :: trcr_edg(nvl,npp,ims:ime,ntra)
    real, intent(in)    :: geop_edg(nvl,npp,ims:ime)
    real, intent(in)    :: mont_edg(nvl,npp,ims:ime)
    real, intent(in)    :: uvsq_edg(nvl,npp,ims:ime)
    real, intent(inout) :: u_tndcy (nvl,ims:ime,nabl)
    real, intent(inout) :: v_tndcy (nvl,ims:ime,nabl)
    real, intent(out)   :: relvor  (nvl,ims:ime)
    real, intent(in)    :: dp3d    (nvl,ims:ime)     ! layer thickness

! Local variables
    real      :: pgfx(nvl,ims:ime)
    real      :: pgfy(nvl,ims:ime)
    real      :: pgfx_gz(nvl),pgfy_gz(nvl)	! PGF based on geopot. (Janjic)
    real      :: pgfx_mp(nvl),pgfy_mp(nvl)	! PGF based on montg.pot.
    integer   :: k		! layer index
    integer   :: ipn		! icos point index
    integer   :: edg		! icos edge index
    integer   :: ndamp		! number of damped layers
    character :: string*24
    logical   :: vrbos
    real      :: factor,wgt
    real      :: uzeta,vzeta
    real      :: uold,vold,qq

!sms$ignore begin

    ndamp = 0.25*nvl		! number of layers subjected to damping

!$OMP PARALLEL DO PRIVATE (vrbos,k,edg,uzeta,vzeta,uold,vold,factor,qq,	&
!$OMP          wgt,pgfx_gz,pgfy_gz,pgfx_mp,pgfy_mp) SCHEDULE (runtime)
    do ipn=ips,ipe
      vrbos = ipn == PrintIpnDiag

      ! Initialize line integrals:
      do k=1,nvl
        pgfx_gz(k) = 0.
        pgfy_gz(k) = 0.
        pgfx_mp(k) = 0.
        pgfy_mp(k) = 0.
        relvor(k,ipn) = 0.
      end do

      ! loop through edges and compute line integrals of bernoulli function,
      ! potential temperature, and pressure gradient

      do edg=1,nprox(ipn)
        do k=1,nvl

          !  Janjic single-term pressure gradient
          if (janjic.gt.0) then
            pgfx_gz(k) = pgfx_gz(k)					&
              -(geop_edg(k,edg,ipn)+uvsq_edg(k,edg,ipn))		&
              *sidevec_c(2,edg,ipn)
            pgfy_gz(k) = pgfy_gz(k)					&
              +(geop_edg(k,edg,ipn)+uvsq_edg(k,edg,ipn))		&
              *sidevec_c(1,edg,ipn)
          end if

          !  traditional 2-term pressure gradient
            pgfx_mp(k) = pgfx_mp(k)					&
!              -bnll_edg(k,edg,ipn)*sidevec_c(2,edg,ipn)		&
               -(mont_edg(k,edg,ipn)					&
                +uvsq_edg(k,edg,ipn))*sidevec_c(2,edg,ipn)		&
               +.5*(exner(k,ipn)+exner(k+1,ipn))			&
               *trcr_edg(k,edg,ipn,1)*sidevec_c(2,edg,ipn)
            pgfy_mp(k) = pgfy_mp(k)					&
!              +bnll_edg(k,edg,ipn)*sidevec_c(1,edg,ipn)		&
               +(mont_edg(k,edg,ipn)					&
                +uvsq_edg(k,edg,ipn))*sidevec_c(1,edg,ipn)		&
               -.5*(exner(k,ipn)+exner(k+1,ipn))			&
               *trcr_edg(k,edg,ipn,1)*sidevec_c(1,edg,ipn)

          !  Vorticity is calculated as a line integral of the tangential wind
          !  component given by the dot product of the wind vector with sidevec.
          !  Sidevec is a vectorial representation of the edge.
 
          relvor(k,ipn) = relvor(k,ipn)					&
               + ((sidevec_e(1,edg,ipn)*u_edg(k,edg,ipn)		&
               +   sidevec_e(2,edg,ipn)*v_edg(k,edg,ipn))) 
        end do
      end do

! --- blend Janjic and non-Janjic PGF formulae

      do k=1,nvl                ! vertical loop
        if (janjic.eq.0) then
          pgfx(k,ipn)=pgfx_mp(k)
          pgfy(k,ipn)=pgfy_mp(k)
        else if (janjic.eq.1) then		! Janjic in lower layers
          wgt=max(0.,min(1.,(.5*nvl-k)/(.1*nvl))) ! wgt=1 near srfc, 0 near top
          pgfx(k,ipn)=wgt*pgfx_gz(k)+(1.-wgt)*pgfx_mp(k)
          pgfy(k,ipn)=wgt*pgfy_gz(k)+(1.-wgt)*pgfy_mp(k)
        else					! Janjic in all layers
          pgfx(k,ipn)=pgfx_gz(k)
          pgfy(k,ipn)=pgfy_gz(k)
        end if
      end do

#ifdef DEBUGPRINT
      if (vrbos) then
        write (*,100) its,perm(ipn)
100     format ('m o m t u m   time step',i6,'  ipn =',i8/		&
             4x,'   uold   unew  utdcy  gradp  corio',			&
             4x,'   vold   vnew  vtdcy  gradp  corio')
      end if
#endif

      do k=1,nvl			! loop through layers
        pgfx  (k,ipn) = pgfx  (k,ipn)*rarea(ipn)
        pgfy  (k,ipn) = pgfy  (k,ipn)*rarea(ipn)
        relvor(k,ipn) = relvor(k,ipn)*rarea(ipn)

        uzeta = (corio(ipn) + relvor(k,ipn))*u_velo(k,ipn)
        vzeta = (corio(ipn) + relvor(k,ipn))*v_velo(k,ipn)

        !  u/v tendcy is the sum of bernoulli fct. gradient and coriolis term
        u_tndcy(k,ipn,nf) = pgfx(k,ipn) + vzeta
        v_tndcy(k,ipn,nf) = pgfy(k,ipn) - uzeta

        ! advance velocity field in time
        uold = u_velo(k,ipn)
        vold = v_velo(k,ipn)

        u_velo(k,ipn) = u_velo(k,ipn)				&
             +adbash1*u_tndcy(k,ipn, nf)			&
             +adbash2*u_tndcy(k,ipn, of)			&
             +adbash3*u_tndcy(k,ipn,vof)
        v_velo(k,ipn) = v_velo(k,ipn)				&
             +adbash1*v_tndcy(k,ipn, nf)			&
             +adbash2*v_tndcy(k,ipn, of)			&
             +adbash3*v_tndcy(k,ipn,vof)
#ifdef DEBUGPRINT
        if (vrbos .and. mod(k,7) == 1) then
          write (*,101) k, &
            uold,u_velo(k,ipn),u_tndcy(k,ipn,nf)*dt,pgfx(k,ipn)*dt, vzeta*dt, &
            vold,v_velo(k,ipn),v_tndcy(k,ipn,nf)*dt,pgfy(k,ipn)*dt,-uzeta*dt
101       format (i4,2f7.1,3f7.2,4x,2f7.1,3f7.2)
        end if
#endif
      end do
    
! --- Rayleigh damping of u,v near model top.
! --- Damping coefficient is proportional to top layer wind speed
! --- weighted by a factor dependent on the excess of wind speed over
! --- a threshold value. The factor is zero if speed < threshold.

#ifdef DEBUGPRINT      
      if (vrbos) then
        write (*,107) perm(ipn),'u,v bfore Rayleigh damping:',	&
             (k,u_velo(k,ipn),v_velo(k,ipn),k=nvl-ndamp,nvl)
107     format ('ipn=',i8,3x,a/(i15,2f9.2))
      end if
#endif

      qq=max(u_velo(nvl  ,ipn)**2 + v_velo(nvl  ,ipn)**2,	&
             u_velo(nvl-1,ipn)**2 + v_velo(nvl-1,ipn)**2,	&
             u_velo(nvl-2,ipn)**2 + v_velo(nvl-2,ipn)**2)/rlthresh**2
      if (qq > 1.) then 
        qq=sqrt((qq-1.)*qq)
        do k=nvl-ndamp,nvl
          factor=1. - (dt/86400.)*rleigh*qq * 2.**(rltaper*(k-nvl))
          u_velo(k,ipn) = u_velo(k,ipn)*factor
          v_velo(k,ipn) = v_velo(k,ipn)*factor
        end do
      end if
#ifdef DEBUGPRINT
      if (vrbos) then
        write (*,107) perm(ipn),'u,v after Rayleigh damping:',	&
             (k,u_velo(k,ipn),v_velo(k,ipn),k=nvl-ndamp,nvl)
      end if
#endif
    end do		                ! horizontal loop
!$OMP END PARALLEL DO

#ifdef DEBUGPRINT
    if (PrintIpnDiag > 0) then
      write (string,'(a,i6,2x)') 'step',its
      call stencl(u_velo,nvl,1.,string(1:12)//'(atm momtum) new u')
      call stencl(v_velo,nvl,1.,string(1:12)//'(atm momtum) new v')
!   call stencl(pgfx,nvl,1.e3,string(1:12)//'(atm momtum) pgfx x 1000')
!   call stencl(pgfy,nvl,1.e3,string(1:12)//'(atm momtum) pgfy x 1000')
    end if
#endif

#ifdef DEBUGPRINT
    do k=1,nvl,7
     write (string,'(a,i3,a)') 'k=',k,' u_tndcy'
     call findmxmn3(u_tndcy,nvl,nip,3,k,nf,string)
  
     write (string,'(a,i3,a)') 'k=',k,' v_tndcy'
     call findmxmn3(v_tndcy,nvl,nip,3,k,nf,string)
    end do
#endif

    return

!sms$ignore end

  end subroutine momtum1

end module momtum_compute
