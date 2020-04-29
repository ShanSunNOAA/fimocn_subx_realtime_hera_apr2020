module module_edgvar1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Nothing below needs to be processed by SMS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use stenedgprint
  use global_bounds,   only: ims, ime, ips, ipe, ihe
  use module_control  ,only: nvlp1, npp, ntra
  use fimnamelist     ,only: nvl, PrintIpnDiag, janjic, eqwgt, hiOrderFluxComp
  use module_constants,only: cs, sn, nedge, permedge, perm,		&
                             prox, nprox, idxtrl, idxldg, wgttrl, wgtldg

  implicit none

#include <gptl.inc>
  
contains

!*********************************************************************
!	edgvar1
!	Interpolates data to cell edges in the local stereographic grid
!	J. Lee                    Sep 2005
!	A. E. MacDonald		  Nov 2005  fim conversion
!	R. Bleck                  Apr 2008  cosmetic changes
!*********************************************************************

  subroutine edgvar1 (its, u_vel, v_vel, delp, pres, exner,              &
                     geop, montg, tracr, u_edg, v_edg, dp_edg,          &
                     lp_edg, geop_edg, mont_edg, uvsq_edg, trc_edg)
! Arguments
    integer, intent(in) :: its				! model time step
    real, intent(in) :: u_vel (nvl,ims:ime)		! west wind on s level
    real, intent(in) :: v_vel (nvl,ims:ime)		! south wind on s level
    real, intent(in) :: delp  (nvl,ims:ime)		! layer thickness
    real, intent(in) :: pres  (nvlp1,ims:ime)		! pres on interfaces
    real, intent(in) :: geop  (nvlp1,ims:ime)		! geopotential
    real, intent(in) :: exner (nvlp1,ims:ime)		! exner fct on interfc
    real, intent(in) :: montg (nvl,ims:ime)		! montgomery potential
    real, intent(in) :: tracr (nvl,ims:ime,ntra)	! mass field tracers
    real, intent(out)   :: u_edg   (nvl,npp,ims:ime)	! west wind on edges
    real, intent(out)   :: v_edg   (nvl,npp,ims:ime)	! south wind on edges
    real, intent(out)   :: dp_edg  (nvl,npp,ims:ime)	! layer thkns on edges
    real, intent(out)   :: lp_edg  (nvl,npp,ims:ime)	! midlayer prs on edges
    real, intent(out)   :: trc_edg (nvl,npp,ims:ime,ntra) ! tracers on edges
    real, intent(out)   :: geop_edg(nvl,npp,ims:ime)	! geopot. on edges
    real, intent(out)   :: mont_edg(nvl,npp,ims:ime)	! montg.pot. on edges
    real, intent(out)   :: uvsq_edg(nvl,npp,ims:ime)	! vel.-squared on edges

! Local workspace
    real    :: exmid(nvlp1)
    real    :: geopx(nvlp1,npp,ims:ime)		! geopot at neighboring points
!   real    :: work1d(ims:ime)
!   real    :: work2d(npp,ims:ime)
    integer :: edg                   ! icos edge index
    integer :: ipx                   ! neighbor across edge with index 'edg'
    integer :: em1,ep1               ! edges edg-1 and edg+1 (cyclic)
    integer :: im1,ip1               ! neighbors across edg-1 and edg+1
    integer :: i1trl,i2trl,i3trl     ! cell indices for interpolation to corners
    integer :: i1ldg,i2ldg,i3ldg     ! cell indices for interpolation to corners
    integer :: ipn                   ! icos index
    integer :: k,kx,ka,kb 	     ! layer Index
    integer :: nt                    ! tracer index: 1=theta, 2=qv, ....
    integer :: edgcount              ! count of icos edges
!    These are u and v at neighboring icos points, NOT on edge
    real    :: utrl,vtrl,uldg,vldg   ! u,v on leading/trailing endpts on edges
    real    :: u1trl,u2trl,u3trl,v1trl,v2trl,v3trl
    real    :: u1ldg,u2ldg,u3ldg,v1ldg,v2ldg,v3ldg
    real    :: m_1,m_2,m_1ldg,m_2ldg,m_3trl,m_3ldg
    real    :: var,exbot,extop,thbot,thtop,thav
    logical :: vrbos
    character(len=32) :: string
    logical,parameter :: simpson = .false.	! use simpson integral rule
    real, parameter :: divby18=1./18.
    real, parameter :: divby6=1./6.
    real, parameter :: divby3=1./3.
    integer :: ret
    real, parameter :: eps=0.001

! NOTE: Committed code has the following compile-time settings (as opposed to run-time settings)
! to allow the compiler to optimize:
! parameter :: janjic = 0
! parameter :: eqwgt = .true.
! parameter :: simpson = .false.

    ret = gptlstart ('edgvar1')
    if (janjic.gt.0) then

! --- integrate hydrostatic eqn at neighboring points to p levels in mid-column

! Halo comp means horizontal loop indices are ips,ihe
!$OMP PARALLEL DO PRIVATE (vrbos,edg,ipn,ipx,k,kx,ka,kb,var,		&
!$OMP             thbot,thtop,thav,exbot,extop,exmid) SCHEDULE (runtime)
    do ipn=ips,ipe
      vrbos = ipn == PrintIpnDiag

      exmid(:)=exner(:,ipn)
      do edg = 1,nprox(ipn)		! loop through edges
        ipx = prox(edg,ipn)
        do k=1,nvlp1
! --- avoid underground line integrations
          exmid(k)=min(exmid(k),exner(1,ipx))
        end do
      end do

      do edg=1,nprox(ipn)
        ipx = prox(edg,ipn)
        kx=1
        do k=1,nvlp1
          geopx(k,edg,ipn)=geop(k,ipn)
  2       if (exner(kx+1,ipx).le.exmid(k) .or. kx.eq.nvl) go to 1
          kx=kx+1
          go to 2
  1       if (exner(kx,ipx).ge.exmid(k)) then
! --- option 1: keep theta constant within each layer
!!          thav=tracr(kx,ipx,1)
! --- option 2: allow theta to vary linearly (use PLM limiters)
            ka=max(  1,kx-1)
            kb=min(nvl,kx+1)
            var=.5*min(tracr(kb,ipx,1)-tracr(kx,ipx,1),		&
                       tracr(kx,ipx,1)-tracr(ka,ipx,1))
            thbot=tracr(kx,ipx,1)-var
            thtop=tracr(kx,ipx,1)+var
            exbot=exner(kx  ,ipx)
            extop=exner(kx+1,ipx)
            thav=.5*(thbot+(thbot*(exmid(k)-extop)		&
                           +thtop*(exbot-exmid(k)))		&
                                 /(exbot   -extop))
            geopx(k,edg,ipn)=geop(kx,ipx)+thav*(exner(kx,ipx)-exmid(k))
          else 
            print '(2(i7,i3,a),f11.2)',ipn,k,'  no values in column ',	&
             perm(ipx),kx,'  bracketing exner=',exmid(k)
          end if
        end do			!  vertical loop
#ifdef DEBUGPRINT
        if (vrbos) then
          print '(2(a,i8),a/(i5,2(f15.2,f11.1),f15.1))',		&
          '(edgvar1)   neighb column',perm(ipx),'    central column',	&
          ipn,'    interpol geop',					&
          (k,exner(k,ipx),geop(k,ipx),exner(k,ipn),geop(k,ipn),		&
          geopx(k,edg,ipn),k=1,nvlp1)
        end if
#endif
      end do			! loop over edges
    end do			! horiz. loop
!$OMP END PARALLEL DO

    if (PrintIpnDiag > 0) then
      write (string,'(a)') '(atm edgvar1) geopotential'
      call stenedg(geop,geopx,nvlp1,trim(string)//', cell & edge')
    end if
  end if			! janjic > 0

!print*, 'in edgvar1, min and max of wgttrl', minval(wgttrl), maxval(wgttrl)
!print*, 'in edgvar1, wgttrl(10)',wgttrl(:,:,10) 
! Halo comp means horizontal loop indices are ips,ihe
!$OMP PARALLEL DO PRIVATE (vrbos,edgcount,edg,em1,ep1,ipx,im1,ip1,k,	&
!$OMP                      i1trl,i2trl,i3trl,i1ldg,i2ldg,i3ldg,		&
!$OMP                      u1trl,u2trl,u3trl,u1ldg,u2ldg,u3ldg,		&
!$OMP                      v1trl,v2trl,v3trl,v1ldg,v2ldg,v3ldg,		&
!$OMP                      m_1,m_2,m_1ldg,m_2ldg,m_3trl,m_3ldg,         &
!$OMP                      utrl,vtrl,uldg,vldg) SCHEDULE (runtime)
  do ipn=ips,ihe
    vrbos = ipn == PrintIpnDiag .and. ipn == perm(ipn)
    do edgcount=1,nedge(ipn)
  
! --- edge quantities are interpolated from 4 icos cells -- the cells on
! --- either side of the edge plus the 2 immediate neighbors of this pair.
  
      edg = permedge(edgcount,ipn)
      i1trl=idxtrl(edg,1,ipn)		! same as ipn
      i2trl=idxtrl(edg,2,ipn)		! cell across edge 'edg' from ipn
      i3trl=idxtrl(edg,3,ipn)		! cell across edge 'edg-1' from ipn
      i1ldg=idxldg(edg,1,ipn)		! same as ipn
      i2ldg=idxldg(edg,2,ipn)		! cell across edge 'edg' from ipn
      i3ldg=idxldg(edg,3,ipn)		! cell across edge 'edg+1' from ipn
#ifdef DEBUGPRINT
      if (vrbos) then
        print 100,its,perm(ipn),edg,					&
         'glb:', perm(i1trl),perm(i2trl),perm(i3trl),			&
                 perm(i1ldg),perm(i2ldg),perm(i3ldg),			&
         'loc:',        i1trl ,       i2trl ,       i3trl ,		&
                        i1ldg ,       i2ldg ,       i3ldg 
 100    format ('its=',i5,' ipn=',i8,' (edgvar1) edge',i2,		&
          '  interpolation based on cells...',2(/a,2(i14,2i8)))
        print 103,'wgt:',						&
         wgttrl(edg,1,ipn),wgttrl(edg,2,ipn),wgttrl(edg,3,ipn),		&
         wgtldg(edg,1,ipn),wgtldg(edg,2,ipn),wgtldg(edg,3,ipn)
 103    format (a,2(f14.4,2f8.4))
      end if
#endif
      do k=1,nvl

        !    Transform u,v at neighboring icos pt to local coord.system.
        !    cs and sn are coordinate transformation constants.
        !    u_xy,v_xy are values of u and v rotated into local system.
        !    (Unrolled to allow vectorization of the k loop)
 
        !   interpolate rotated wind components to edges
 
        if (.not. eqwgt) then
          u1trl =  cs(1,edg,ipn)*u_vel(k,i1trl) + sn(1,edg,ipn)*v_vel(k,i1trl)
          v1trl = -sn(1,edg,ipn)*u_vel(k,i1trl) + cs(1,edg,ipn)*v_vel(k,i1trl)
          u2trl =  cs(2,edg,ipn)*u_vel(k,i2trl) + sn(2,edg,ipn)*v_vel(k,i2trl)
          v2trl = -sn(2,edg,ipn)*u_vel(k,i2trl) + cs(2,edg,ipn)*v_vel(k,i2trl)
        endif
        u3trl =  cs(3,edg,ipn)*u_vel(k,i3trl) + sn(3,edg,ipn)*v_vel(k,i3trl)
        v3trl = -sn(3,edg,ipn)*u_vel(k,i3trl) + cs(3,edg,ipn)*v_vel(k,i3trl)

        if (.not. eqwgt) then
          utrl = u1trl*wgttrl(edg,1,ipn) + u2trl*wgttrl(edg,2,ipn)	&
               + u3trl*wgttrl(edg,3,ipn)
          vtrl = v1trl*wgttrl(edg,1,ipn) + v2trl*wgttrl(edg,2,ipn)	&
               + v3trl*wgttrl(edg,3,ipn)
        endif

        u1ldg =  cs(1,edg,ipn)*u_vel(k,i1ldg) + sn(1,edg,ipn)*v_vel(k,i1ldg)
        v1ldg = -sn(1,edg,ipn)*u_vel(k,i1ldg) + cs(1,edg,ipn)*v_vel(k,i1ldg)
        u2ldg =  cs(2,edg,ipn)*u_vel(k,i2ldg) + sn(2,edg,ipn)*v_vel(k,i2ldg)
        v2ldg = -sn(2,edg,ipn)*u_vel(k,i2ldg) + cs(2,edg,ipn)*v_vel(k,i2ldg)
        u3ldg =  cs(4,edg,ipn)*u_vel(k,i3ldg) + sn(4,edg,ipn)*v_vel(k,i3ldg)
        v3ldg = -sn(4,edg,ipn)*u_vel(k,i3ldg) + cs(4,edg,ipn)*v_vel(k,i3ldg)

        if (.not. eqwgt) then
          uldg = u1ldg*wgtldg(edg,1,ipn) + u2ldg*wgtldg(edg,2,ipn)	&
               + u3ldg*wgtldg(edg,3,ipn)
          vldg = v1ldg*wgtldg(edg,1,ipn) + v2ldg*wgtldg(edg,2,ipn)	&
               + v3ldg*wgtldg(edg,3,ipn)
        endif

! --- trapezoidal:
        if (eqwgt) then
          u_edg(k,edg,ipn) = (2.*(u1ldg+u2ldg)+u3trl+u3ldg)*divby6
          v_edg(k,edg,ipn) = (2.*(v1ldg+v2ldg)+v3trl+v3ldg)*divby6
          m_1ldg = sqrt(u1ldg*u1ldg+v1ldg*v1ldg)
          m_2ldg = sqrt(u2ldg*u2ldg+v2ldg*v2ldg)
          m_3trl = sqrt(u3trl*u3trl+v3trl*v3trl)
          m_3ldg = sqrt(u3ldg*u3ldg+v3ldg*v3ldg)
          m_1 = sqrt(u_edg(k,edg,ipn)*u_edg(k,edg,ipn)+v_edg(k,edg,ipn)*v_edg(k,edg,ipn))
          m_2 = (2.*(m_1ldg+m_2ldg)+m_3trl+m_3ldg)*divby6
          u_edg(k,edg,ipn) = u_edg(k,edg,ipn)*m_2/max(eps,m_1)
          v_edg(k,edg,ipn) = v_edg(k,edg,ipn)*m_2/max(eps,m_1)
        else
          u_edg(k,edg,ipn) = .5*(utrl+uldg)
          v_edg(k,edg,ipn) = .5*(vtrl+vldg)
        end if

! --- simpson (assuming -i1,i2- are cells nearest edge midpoint):
        if (simpson) then
          if (eqwgt) then
            u_edg(k,edg,ipn) = (8.*(u1ldg+u2ldg)+u3trl+u3ldg)*divby18
            v_edg(k,edg,ipn) = (8.*(v1ldg+v2ldg)+v3trl+v3ldg)*divby18
          else
            u_edg(k,edg,ipn) = (u_edg(k,edg,ipn)+u1trl+u2trl)*divby3
            v_edg(k,edg,ipn) = (v_edg(k,edg,ipn)+v1trl+v2trl)*divby3
          end if
        end if

#ifdef DEBUGPRINT
        if (.not. eqwgt) then
          if (u1ldg.ne.u1trl .or. u2ldg.ne.u2trl .or.			&
              v1ldg.ne.v1trl .or. v2ldg.ne.v2trl) then
            print 104,'(edgvar) error 1:',perm(ipn),k,edg,			&
              u1ldg,u1trl,u2ldg,u2trl,v1ldg,v1trl,v2ldg,v2trl
 104        format (a,i8,2i3/8es14.6)
            stop '(edgvar1 error 1)'
          end if
          if (vrbos .and. mod(k,9).eq.1) then
            print 102,perm(ipn),'  (edgvar1) k,edg=',k,edg,'-trl',		&
              'i1,i2,i3=',perm(i1trl),perm(i2trl),perm(i3trl),		&
              'u1trl,u2trl,u3trl=',u1trl,u2trl,u3trl,			&
              'v1trl,v2trl,u3trl=',v1trl,v2trl,v3trl,			&
              '        utrl,vtrl=',utrl,vtrl
            print 102,perm(ipn),'  (edgvar1) k,edg=',k,edg,'-ldg',		&
              'i1,i2,i3=',perm(i1ldg),perm(i2ldg),perm(i3ldg),		&
              'u1ldg,u2ldg,u3ldg=',u1ldg,u2ldg,u3ldg,			&
              'v1ldg,v2ldg,u3ldg=',v1ldg,v2ldg,v3ldg,			&
              '        uldg,vldg=',uldg,vldg
 102        format (i7,a,2i4,a,4x,a,3i8/(a,3f9.2))
          end if
        end if   ! .not. eqwgt
#endif

      !   interpolate layer thickness to edges
        if (.not. hiOrderFluxComp) then 
! --- trapezoidal:
        if (eqwgt) then
          dp_edg(k,edg,ipn)=					&
            (2.*(delp(k,i1ldg)+delp(k,i2ldg))			&
                +delp(k,i3trl)+delp(k,i3ldg))*divby6
        else
          dp_edg(k,edg,ipn) =					&
             ( delp(k,i1trl)*wgttrl(edg,1,ipn)			&
             + delp(k,i2trl)*wgttrl(edg,2,ipn)			&
             + delp(k,i3trl)*wgttrl(edg,3,ipn)			&
             + delp(k,i1ldg)*wgtldg(edg,1,ipn)			&
             + delp(k,i2ldg)*wgtldg(edg,2,ipn)			&
             + delp(k,i3ldg)*wgtldg(edg,3,ipn))*.5
        end if

! --- simpson (assuming -i1,i2- are cells nearest edge midpoint):
        if (simpson) then
          if (eqwgt) then
            dp_edg(k,edg,ipn)=					&
              (8.*(delp(k,i1ldg)+delp(k,i2ldg))			&
                  +delp(k,i3trl)+delp(k,i3ldg))*divby18
          else
            dp_edg(k,edg,ipn) = (dp_edg(k,edg,ipn)		&
              +delp(k,i1trl)+delp(k,i2trl))*divby3
          end if
        end if
        endif

      !   interpolate mid-layer pressure to edges
 
! --- trapezoidal:
        if (eqwgt) then
          lp_edg(k,edg,ipn)=					&
                (pres(k,i1ldg)+pres(k+1,i1ldg)			&
                +pres(k,i2ldg)+pres(k+1,i2ldg)			&
            +.5*(pres(k,i3trl)+pres(k+1,i3trl)			&
                +pres(k,i3ldg)+pres(k+1,i3ldg)))*divby6
        else
          lp_edg(k,edg,ipn) =						&
             ((pres(k,i1trl)+pres(k+1,i1trl))*wgttrl(edg,1,ipn)		&
             +(pres(k,i2trl)+pres(k+1,i2trl))*wgttrl(edg,2,ipn)		&
             +(pres(k,i3trl)+pres(k+1,i3trl))*wgttrl(edg,3,ipn)		&
             +(pres(k,i1ldg)+pres(k+1,i1ldg))*wgtldg(edg,1,ipn)		&
             +(pres(k,i2ldg)+pres(k+1,i2ldg))*wgtldg(edg,2,ipn)		&
             +(pres(k,i3ldg)+pres(k+1,i3ldg))*wgtldg(edg,3,ipn))*.25    
        end if

! --- simpson (assuming -i1,i2- are cells nearest edge midpoint):
        if (simpson) then
          if (eqwgt) then
            lp_edg(k,edg,ipn)=					&
              (4.*(pres(k  ,i1ldg)+pres(k  ,i2ldg)		&
                  +pres(k+1,i1ldg)+pres(k+1,i2ldg))		&
              +.5*(pres(k  ,i3trl)+pres(k  ,i3ldg)		&
                  +pres(k+1,i3trl)+pres(k+1,i3ldg)))*divby18
          else
            lp_edg(k,edg,ipn) = (lp_edg(k,edg,ipn)		&
              +.5*(pres(k,i1trl)+pres(k+1,i1trl)		&
                  +pres(k,i2trl)+pres(k+1,i2trl)))*divby3
          end if
        end if
 
      !   interpolate montgomery potential to edges

! --- trapezoidal:
        if (eqwgt) then
          mont_edg(k,edg,ipn)=					&
            (2.*(montg(k,i1ldg)+montg(k,i2ldg))			&
                +montg(k,i3trl)+montg(k,i3ldg))*divby6
        else
          mont_edg(k,edg,ipn) =					&
            ( montg(k,i1trl)*wgttrl(edg,1,ipn)			&
            + montg(k,i2trl)*wgttrl(edg,2,ipn)			&
            + montg(k,i3trl)*wgttrl(edg,3,ipn)			&
            + montg(k,i1ldg)*wgtldg(edg,1,ipn)			&
            + montg(k,i2ldg)*wgtldg(edg,2,ipn)			&
            + montg(k,i3ldg)*wgtldg(edg,3,ipn))*.5
        end if

! --- simpson (assuming -i1,i2- are cells nearest edge midpoint):
        if (simpson) then
          if (eqwgt) then
            mont_edg(k,edg,ipn)=				&
              (8.*(montg(k,i1ldg)+montg(k,i2ldg))		&
                  +montg(k,i3trl)+montg(k,i3ldg))*divby18
          else
            mont_edg(k,edg,ipn) = (mont_edg(k,edg,ipn)		&
              +montg(k,i1trl)+montg(k,i2trl))*divby3
          end if
        end if

#ifdef DEBUGPRINT
        if (delp(k,i1trl).ne.delp(k,i1ldg) .or.			&
            delp(k,i2trl).ne.delp(k,i2ldg)) then
         print 104,'(edgvar) error 2:',perm(ipn),k,edg,		&
          delp(k,i1trl),delp(k,i1ldg),delp(k,i2trl),delp(k,i2ldg)
         stop '(edgvar1 error 2)'
        end if
        if (pres(k,i1trl)+pres(k+1,i1trl) .ne.			&
            pres(k,i1ldg)+pres(k+1,i1ldg) .or.			&
            pres(k,i2trl)+pres(k+1,i2trl) .ne.			&
            pres(k,i2ldg)+pres(k+1,i2ldg)) then
         print 104,'(edgvar) error 3:',perm(ipn),k,edg,		&
          pres(k,i1trl)+pres(k+1,i1trl),			&
          pres(k,i1ldg)+pres(k+1,i1ldg),			&
          pres(k,i2trl)+pres(k+1,i2trl),			&
          pres(k,i2ldg)+pres(k+1,i2ldg)
          stop '(edgvar1 error 3)'
        end if
        if (montg(k,i1trl).ne.montg(k,i1ldg) .or.		&
            montg(k,i2trl).ne.montg(k,i2ldg)) then
         print 104,'(edgvar) error 4:',perm(ipn),k,edg,		&
          montg(k,i1trl),montg(k,i1ldg),montg(k,i2trl),montg(k,i2ldg)
         stop '(edgvar1 error 4)'
        end if
        if (vrbos .and. mod(k,9).eq.1) then
          print 101,perm(ipn),'  (edgvar1) k,edg=',		&
            k,edg,'-trl  ipn,montg,wgt 1/2/3=',			&
            perm(i1trl),montg(k,i1trl),wgttrl(edg,1,ipn),	&
            perm(i2trl),montg(k,i2trl),wgttrl(edg,2,ipn),	&
            perm(i3trl),montg(k,i3trl),wgttrl(edg,3,ipn)
            print 101,perm(ipn),'  (edgvar1) k,edg=',		&
            k,edg,'-ldg  ipn,montg,wgt 1/2/3=',			&
            perm(i1ldg),montg(k,i1ldg),wgtldg(edg,1,ipn),	&
            perm(i2ldg),montg(k,i2ldg),wgtldg(edg,2,ipn),	&
            perm(i3ldg),montg(k,i3ldg),wgtldg(edg,3,ipn)
        end if
#endif

        !   interpolate kinetic energy to edges
 
        uvsq_edg(k,edg,ipn) = .5*(u_edg(k,edg,ipn)**2		&
                                + v_edg(k,edg,ipn)**2)
      end do			! vertical loop
    end do			! loop over edges
  end do			! horizontal loop
!$OMP END PARALLEL DO

  !   interpolate tracer to edges
 
  do nt=1,ntra			! loop through tracers
    if (hiOrderFluxComp .and. nt > 1) exit 

! Halo comp means horizontal loop indices are ips,ihe
!$OMP PARALLEL DO PRIVATE (vrbos,k,edgcount,edg,em1,ep1,ipn,ipx,	&
!$OMP          im1,ip1,i1trl,i2trl,i3trl,i1ldg,i2ldg,i3ldg)		&
!$OMP          SCHEDULE (runtime)
    do ipn=ips,ihe
      vrbos=ipn.eq.PrintIpnDiag .and. ipn.eq.perm(ipn)
      do k=1,nvl
        trc_edg(k,npp,ipn,nt) = 0.
      end do

      do edgcount=1,nedge(ipn)	! loop through edges  
        edg = permedge(edgcount,ipn)
        i1trl=idxtrl(edg,1,ipn)		! same as ipn
        i2trl=idxtrl(edg,2,ipn)		! cell across edge 'edg' from ipn
        i3trl=idxtrl(edg,3,ipn)		! cell across edge 'edg-1' from ipn
        i1ldg=idxldg(edg,1,ipn)		! same as ipn
        i2ldg=idxldg(edg,2,ipn)		! cell across edge 'edg' from ipn
        i3ldg=idxldg(edg,3,ipn)		! cell across edge 'edg+1' from ipn

        do k=1,nvl
! --- trapezoidal:
          if (eqwgt) then
            trc_edg(k,edg,ipn,nt)=					&
              (2.*(tracr(k,i1ldg,nt)+tracr(k,i2ldg,nt))			&
                  +tracr(k,i3trl,nt)+tracr(k,i3ldg,nt))*divby6
          else
            trc_edg(k,edg,ipn,nt) =					&
                ( tracr(k,i1trl,nt)*wgttrl(edg,1,ipn)			&
                + tracr(k,i2trl,nt)*wgttrl(edg,2,ipn)			&
                + tracr(k,i3trl,nt)*wgttrl(edg,3,ipn)			&
                + tracr(k,i1ldg,nt)*wgtldg(edg,1,ipn)			&
                + tracr(k,i2ldg,nt)*wgtldg(edg,2,ipn)			&
                + tracr(k,i3ldg,nt)*wgtldg(edg,3,ipn))*.5
          end if

! --- simpson (assuming -i1,i2- are cells nearest edge midpoint):
          if (simpson) then
            if (eqwgt) then
              trc_edg(k,edg,ipn,nt)=					&
                (8.*(tracr(k,i1ldg,nt)+tracr(k,i2ldg,nt))		&
                    +tracr(k,i3trl,nt)+tracr(k,i3ldg,nt))*divby18
            else
              trc_edg(k,edg,ipn,nt) = (trc_edg(k,edg,ipn,nt)		&
                +tracr(k,i1trl,nt)+tracr(k,i2trl,nt))*divby3
            end if
          end if
#ifdef DEBUGPRINT   
          if (tracr(k,i1trl,nt).ne.tracr(k,i1ldg,nt) .or.		&
              tracr(k,i2trl,nt).ne.tracr(k,i2ldg,nt)) then
           print 104,'(edgvar) error 5:',perm(ipn),k,edg,		&
            tracr(k,i1trl,nt),tracr(k,i1ldg,nt),			&
            tracr(k,i2trl,nt),tracr(k,i2ldg,nt)
           stop '(edgvar1 error 5)'
          end if
          if (vrbos .and. mod(k,9).eq.1 .and. nt.eq.1) then             
           print 101,perm(ipn),'  (edgvar1) k,edg=',			&
             k,edg,'-trl  ipn,tracer 1,wgt 1/2/3=',			&
             perm(i1trl),tracr(k,i1trl,nt),wgttrl(edg,1,ipn),		&
             perm(i2trl),tracr(k,i2trl,nt),wgttrl(edg,2,ipn),		&
             perm(i3trl),tracr(k,i3trl,nt),wgttrl(edg,3,ipn)          
           print 101,perm(ipn),'  (edgvar1) k,edg=',			&
             k,edg,'-ldg  ipn,tracer 1,wgt 1/2/3=',			&
             perm(i1ldg),tracr(k,i1ldg,nt),wgtldg(edg,1,ipn),		&
             perm(i2ldg),tracr(k,i2ldg,nt),wgtldg(edg,2,ipn),		&
             perm(i3ldg),tracr(k,i3ldg,nt),wgtldg(edg,3,ipn)
          end if
 101      format (i7,a,2i3,a/3(i10,f9.1,f7.3))
#endif

        end do
      end do			!  loop over edges
    end do			!  horizontal loop
!$OMP END PARALLEL DO
 
!   if (nt.eq.1) then
!     write (string,'(a,i2)') '(atm edgvar1) tracer',nt
!     do k=1,nvl,7
!       work1d(:)=tracr(k,:,nt)
!       work2d(:,:)=trc_edg(k,:,:,nt)
!       call stenedg(work1d,work2d,1,trim(string)//', cell & edge')
!     end do
!   end if
 
  end do			!  loop through tracers

  if (janjic.gt.0) then

          !   interpolate geopotential to edges
  
!$OMP PARALLEL DO PRIVATE (vrbos,edg,em1,ep1,ipn,k,ipx,ip1,im1) SCHEDULE (runtime)
    do ipn=ips,ipe
      do edg = 1,nprox(ipn)		! loop through edges
        vrbos = ipn == PrintIpnDiag
        ipx = prox(edg,ipn)
        em1 = mod(edg-2+nprox(ipn),nprox(ipn))+1
        ep1 = mod(edg             ,nprox(ipn))+1
        im1 = prox(em1,ipn)
        ip1 = prox(ep1,ipn)

        i1trl=idxtrl(edg,1,ipn)		! same as ipn
        i2trl=idxtrl(edg,2,ipn)		! cell across edge 'edg' from ipn
        i3trl=idxtrl(edg,3,ipn)		! cell across edge 'edg-1' from ipn
        i1ldg=idxldg(edg,1,ipn)		! same as ipn
        i2ldg=idxldg(edg,2,ipn)		! cell across edge 'edg' from ipn
        i3ldg=idxldg(edg,3,ipn)		! cell across edge 'edg+1' from ipn

        do k=1,nvl
! --- trapezoidal:
          if (eqwgt) then
            geop_edg(k,edg,ipn)=					&
                 (geop (k,    ipn)+geop (k+1,    ipn)			&
                 +geopx(k,edg,ipn)+geopx(k+1,edg,ipn)			&
             +.5*(geopx(k,em1,ipn)+geopx(k+1,em1,ipn)			&
                 +geopx(k,ep1,ipn)+geopx(k+1,ep1,ipn)))*divby6
          else
            geop_edg(k,edg,ipn) =					&
	      ((geop (k,    ipn)+geop (k+1,    ipn))*wgttrl(edg,1,ipn)	&
	      +(geopx(k,edg,ipn)+geopx(k+1,edg,ipn))*wgttrl(edg,2,ipn)	&
	      +(geopx(k,em1,ipn)+geopx(k+1,em1,ipn))*wgttrl(edg,3,ipn)	&
	      +(geop (k,    ipn)+geop (k+1,    ipn))*wgtldg(edg,1,ipn)	&
	      +(geopx(k,edg,ipn)+geopx(k+1,edg,ipn))*wgtldg(edg,2,ipn)	&
	      +(geopx(k,ep1,ipn)+geopx(k+1,ep1,ipn))*wgtldg(edg,3,ipn))*.25
          end if
  
! --- simpson:
          if (simpson) then
            if (eqwgt) then
              geop_edg(k,edg,ipn)=					&
                (4.*(geop (k,    ipn)+geop (k+1,    ipn)		&
                    +geopx(k,edg,ipn)+geopx(k+1,edg,ipn))		&
                +.5*(geopx(k,em1,ipn)+geopx(k+1,em1,ipn)		&
                    +geopx(k,ep1,ipn)+geopx(k+1,ep1,ipn)))*divby18
            else
              geop_edg(k,edg,ipn) = (geop_edg(k,edg,ipn)		&
                +.5*(geop (k,    ipn)+geop (k+1,    ipn)		&
                    +geopx(k,edg,ipn)+geopx(k+1,edg,ipn)))*divby3
            end if
          end if
#ifdef DEBUGPRINT
          if (vrbos .and. mod(k,9).eq.1) then
            print 101,ipn,'  (edgvar1) k,edg=',				&
             k,edg,'-trl  ipn,geopot,wgt 1/2/3=',			&
             perm(ipn),geop (k,    ipn),wgttrl(edg,1,ipn),		&
             perm(ipx),geopx(k,edg,ipn),wgttrl(edg,2,ipn),		&
             perm(im1),geopx(k,em1,ipn),wgttrl(edg,3,ipn)
             print 101,ipn,'  (edgvar1) k,edg=',				&
             k,edg,'-ldg  ipn,geopot,wgt 1/2/3=',			&
             perm(ipn),geop (k,    ipn),wgtldg(edg,1,ipn),		&
             perm(ipx),geopx(k,edg,ipn),wgtldg(edg,2,ipn),		&
             perm(ip1),geopx(k,ep1,ipn),wgtldg(edg,3,ipn)
          end if
#endif
        end do			!  vertical loop
      end do			!  loop over edges
    end do			!  horizontal loop
!$OMP END PARALLEL DO

  end if			! janjic > 0
  ret = gptlstop ('edgvar1')
end subroutine edgvar1

end module module_edgvar1
