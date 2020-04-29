module module_hystat
  use findmaxmin1
  use module_control  ,only: nvlp1,ntra
  use fimnamelist     ,only: nvl,ptop,PrintIpnDiag
  use module_constants,only: p1000,cp,rd,qvmin,qwmin,perm
  use global_bounds,   only: ims, ime, ips, ipe

  implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End of SMS region. Nothing below needs to be processed by SMS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains
!*********************************************************************
!     hystat
!	Hydrostatic equation
!	Alexander E. MacDonald  11/14/2005
!	J.Lee                   01/04/2006
!	R.Bleck                 12/06/2007
!*********************************************************************

  subroutine hystat (its, ph3d, ex3d, mp3d, dp3d, &
                     tr3d, trdp, psrf, ptdcy)
! Arguments
    integer, intent (IN) :: its                    ! model time step
    real, intent (INOUT) :: ph3d(nvlp1,ims:ime)	   ! geopotential
    real, intent (IN)    :: ex3d(nvlp1,ims:ime)	   ! exner fct
    real, intent (OUT)   :: mp3d(nvl,ims:ime)	   ! montgomery potential
    real, intent (IN)    :: dp3d(nvl,ims:ime)	   ! layer thickness
    real, intent (INOUT) :: tr3d(nvl,ims:ime,ntra) ! mass field tracers
    real, intent (OUT)   :: trdp(nvl,ims:ime,ntra) ! tracer x thickness
    real, intent (INOUT) :: psrf(ims:ime)          ! surface pressure
!JR Cannot have intent(out) for ptdcy because this routine only sets half of the array, and Lahey
!JR assigns a "bad sequence of bits" to all intent(out) variables.
    real, intent (INOUT) :: ptdcy(ims:ime,2)	! srf.pres.tdcy, 2 time lvls

! local variables
    real :: work(ims:ime)
    integer :: ipn	! icos index
    integer :: k	! layer index
    integer :: ns	! tracer index
    logical :: vrbos	! switch for 'verbose' mode

!  Note that tracers, velocity, and montgomery potential, mp3d
!  are all constant through the layers.  Phi (ph3d) and pressure (dp3d)
!  vary through the layer.

!  Layer variables:  tr3d,mp3d
!  Level variables:  ex3d,ph3d

!sms$ignore begin
!$OMP PARALLEL DO PRIVATE (vrbos,k,ns) SCHEDULE (runtime)
    do ipn=ips,ipe
      vrbos = ipn == PrintIpnDiag

! --- srf.prs tendency is needed to evaluate model noise (diagnoise.F90)
      work(ipn) = psrf(ipn)
      psrf(ipn) = dp3d(nvl,ipn)
      do k=nvl-1,1,-1
        psrf(ipn) = psrf(ipn) + dp3d(k,ipn)
      end do
      ns = mod(its,2)+1
      if (its > 0) then
        ptdcy(ipn,ns) = psrf(ipn) - work(ipn)
      end if

!.........................................................
!   Determine bottom layer values
!.........................................................

      mp3d(1,ipn) = ex3d(1,ipn)*tr3d(1,ipn,1) + ph3d(1,ipn)	! montpot, lyr 1
      ph3d(2,ipn) = mp3d(1,ipn) - tr3d(1,ipn,1)*ex3d(2,ipn)	! geopot, lev 2

      do k=2,nvl
! Hydrostatic eqn:  d mp/d theta = exner fct, tr3d(.,.,1) = theta
        mp3d(k,ipn) = mp3d(k-1,ipn) + ex3d(k,ipn)*		&
          (tr3d(k,ipn,1) - tr3d(k-1,ipn,1))

! get geopotential from identity  montg pot = geopot + theta*exner
        ph3d(k+1,ipn) = mp3d(k,ipn) - tr3d(k,ipn,1)*ex3d(k+1,ipn)
      end do

! keep vapor mixing ratio and water content above prescribed lower limits
      do k=1,nvl
        tr3d(k,ipn,2) = max(qvmin,tr3d(k,ipn,2))
        tr3d(k,ipn,3) = max(qwmin,tr3d(k,ipn,3))
      end do

! compute tracer amount per unit area (tracer concentration x layer thickness)
      do ns=1,ntra
        do k=1,nvl
          trdp(k,ipn,ns) = tr3d(k,ipn,ns)*dp3d(k,ipn)
        end do
      end do

#ifdef DEBUGPRINT
      if (vrbos) then
        write (6,100) its,perm(ipn),(k,1.e3*(ex3d(k,ipn)/cp)**(cp/rd),	&
             ex3d(k,ipn),ph3d(k,ipn),mp3d(k,ipn),tr3d(k,ipn,1),		&
             k=1,nvl),nvlp1,1.e3*(ex3d(nvlp1,ipn)/cp)**(cp/rd),		&
             ex3d(nvlp1,ipn),ph3d(nvlp1,ipn)
100     format ('its,ipn=',i6,i8,					&
         '  HYSTAT   pres  exn.fct    geopot     montg    theta'/	&
         (i28,2f9.2,2f10.1,f9.2))
      end if
#endif

    end do		! horizontal loop 
!$OMP END PARALLEL DO
!sms$ignore end

    return

  end subroutine hystat

end module module_hystat
