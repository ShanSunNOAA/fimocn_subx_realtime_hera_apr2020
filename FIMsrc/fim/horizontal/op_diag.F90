module module_op_diag
  use module_constants,       only: p1000,rd,cp,grvity,sigak,sigbk,	&
                                    grav98,deg_lat,deg_lon,spval_p,perm,&
				    ratio_h20_dry, rovrm1_p
  use module_control,         only: nvlsig,nvlp1,nip,ntra,ntrb,		&
                                    nvlp,nvarp,nvar2d,nvarsig,dt,	&
                                    pres_pa
  use fimnamelist,            only: nvl,glvl,curve,yyyymmddhhmm,	&
                                    gz_smooth,ArchvTimeUnit,		&
				    PrintIpnDiag,pure_sig,enkfio_out,	&
 			 	    TimingBarriers,atmonly
  use module_wrf_control,     only: num_chem

  use machine,                only: kind_rad
  use funcphys,               only: fpvs, fpvsl 
  use module_outqv,           only: outqv
  use module_outqv_mn,        only: outqv_mn
  use module_outqv_mn_lat,    only: outqv_mn_lat
  use module_outqv_mn_lat_abs,only: outqv_mn_lat_abs
  use module_out4d_mn,        only: out4d_mn
  use dffusn,                 only: dffusn_lev
  use global_bounds,          only: ims, ime, ips, ipe
  use findmaxmin3

  implicit none

#include <gptl.inc>

contains
!*********************************************************************
!	op_diag
!	Calculate derived variables for output from FIM global model
!	S. Benjamin    Feb 2008
!       S. Benjamin  - Apr 2008
!                      -  mods for use of theta-v in th3d prog variable
!                         instead of previous non-virt theta in th3d
!       S. Benjamin  - May 2008
!                      - calculation of isobaric 25-mb fields 
!                         on native icosahedral horizontal grid
!                         added in array g3p
!                        Will allow improved isobaric fields to be
!                         output from FIMpost, with improved reduction
!                         from full multivariate calculations using
!                         all fields (which are available in FIM but not FIMpost)
!                    - April 2012
!                      - add output of variables on GFS sigma vertical levels
!        R. Bleck    - June 2012
!                      - changed interpolation to pressure levels
!*********************************************************************

  subroutine op_diag (its, time,                    &
                      us3d, vs3d, dp3d, sdot, pr3d, &
                      ex3d, mp3d, tr, ws3d,         &
                      ph3d, rn2d, rc2d, ts2d, us2d, &
                      hf2d, qf2d, rsds, rlds, st3d, &
                      qv3d, t2m2d,rn_xh,rn_6h,      &
! Below are output variables from op_diag
                      tk3d, rh3d, td3d, pw2d, pq2d, mslp,	&
                      g3p, g3p_chem, g3sig, g2d, t2m_dif,	&
                      sn_xh,sn_6h,wt_int_rev, k_int_rev)

! Arguments:
    integer, intent(in) :: its
    integer, intent(in) :: time

    real, intent(in) :: us3d(nvl  ,ims:ime)       ! west wind
    real, intent(in) :: vs3d(nvl  ,ims:ime)       ! south wind
    real, intent(in) :: dp3d(nvl  ,ims:ime)       ! delta-p
    real, intent(in) :: sdot(nvlp1,ims:ime)       ! vertical velocity
    real, intent(in) :: pr3d(nvlp1,ims:ime)       ! pressure
    real, intent(in) :: ex3d(nvlp1,ims:ime)       ! exner function
    real, intent(in) :: mp3d(nvl  ,ims:ime)       ! mont pot
    real, intent(in) :: ws3d(nvl  ,ims:ime)
    real, intent(in) :: ph3d(nvlp1,ims:ime)       ! geopotential
    real, intent(in) :: rn2d(ims:ime)             ! accumulated precip.
    real, intent(in) :: rn_xh(ims:ime)            ! accum precip. in last xh
    real, intent(in) :: rn_6h(ims:ime)            ! accum precip. in last 6h
    real, intent(in) :: rc2d(ims:ime)             ! rainfall
    real, intent(in) :: tr(nvl,ims:ime,ntra+ntrb) ! tracers

    real, intent(out) :: sn_xh(ims:ime)            ! accum snowfall in last xh
    real, intent(out) :: sn_6h(ims:ime)            ! accum snowfall in last xh

    real, intent(out) :: rh3d(nvl,ims:ime)
    real, intent(out) :: tk3d(nvl,ims:ime)

!JFM Let SMS distribute g3p because it is exchanged
!sms$distribute (dh,2) begin
    real, intent(out) :: g3p(nvlp,nip,nvarp)
!SB - g3p_chem for chem vars output on isobaric levels
    real, intent(out) :: g3p_chem(nvlp,nip,num_chem)
    real, intent(out) :: g3sig(nvlsig+1,nip,nvarsig)
!sms$distribute end

    integer, parameter :: nscalvar = 6   ! # non-chem, non-ht, non-theta vars to interp to isobaric
    real fldlyr(nvl,max(num_chem,nscalvar)),php(nvlp),thp(nvlp),fldp(nvlp,max(num_chem,nscalvar))
!       nvarsig indices
!       1  -  pressure (Pa) at interface levels
!       2  -  pressure thickness (Pa) in layers
!       3  -  temperature (K)
!       4/5-  u/v wind components (m/s)
!       6  -  specific humidity (g/g)
!       7  -  hydrometeor condensate (g/g) ("cloud water" in GFS-speak, but actually more)
!       8  -  ozone mixing ratio (g/g)
!       9+ -  future additional variables (aerosols, etc.)
    real, intent(out) :: g2d(ims:ime,nvar2d)
    real, intent(out) :: t2m_dif(ims:ime)         ! moved to here from output for OpenMP
    REAL, INTENT(inout) :: wt_int_rev(nvl,ims:ime)! weights for reverse interpolation from
!                                                 ! GFS sig levels back to hybrid levels
!                                           For FIM ensemble data assimilation
 INTEGER, INTENT(inout) :: k_int_rev(nvl,ims:ime)! k levs for reverse interpolation from

!       nvar2d indices
!       1  -  height cloud base (in meters above sea level (ASL))
!       2  -  pressure cloud base (Pa)
!       3  -  height cloud top (m ASL)
!       4  -  pressure cloud top (Pa)
!       5  -  column relative humidity with respect to saturated column precipitable water (0-1)
!       6  -  saturated column precipitable water
!       7  -  wind speed (m/s) at 80m above ground level

    real, intent(out) :: td3d(nvl,ims:ime)
    real, intent(out) :: pw2d(ims:ime)
    real, intent(out) :: pq2d(ims:ime)
    real, intent(out) :: mslp(ims:ime)

    real, intent(in) :: ts2d(ims:ime)
    real, intent(in) :: us2d(ims:ime)
    real, intent(in) :: hf2d(ims:ime)
    real, intent(in) :: qf2d(ims:ime)
    real, intent(in) :: rsds(ims:ime)
    real, intent(in) :: rlds(ims:ime)
    real, intent(in) :: st3d(4,ims:ime)
    real, intent(in) :: qv3d(nvl,ims:ime)
    real, intent(in) :: t2m2d(ims:ime)

! Local variables
! Variable names for outqv* printout
    character(len=4), parameter :: varname(7) = (/'Hgt ','Temp','RH  ','Uwnd','Vwnd', &
        'Vvel','qhyd'/)     		! qhyd = hydrometeor (condensate) mixing ratio (cloud water/ice/snow/rain/graupel combined)
    character(len=4), parameter :: varsig(nvarsig) = (/'Pres','Dprs','Temp','Uwnd','Vwnd','sphu','cldq','ozon'/)

    real, parameter :: cloud_def_p = 0.00001
    real, parameter :: gamd = 0.0100  ! 10  K/km - Dry adiabatic lapse rate
    real, parameter :: gams = 0.0065  ! 6.5 K/km - Standard lapse rate

    real :: dummy(nvl,ims:ime)
!  pkap = (p/p0)**(R/Cp)  i.e., p(in bars) to the kappa (Rd/Cp) power
!  pkap = Exner / Cp
    real :: pkap(nvl)
    real :: hgt(nvl)
    real :: pres(nvl)
    real :: tvirt(nvl)
    real :: thv(nvl)
    real :: thnv(nvl)
    real :: pklv(nvlp1)
    real :: pkapp(nvlp)
    real :: prsp(nvlp)
    real :: tsfc, tsfc7
    real :: pkapsig(nvlsig)
    real :: rhpw(ims:ime)
    real :: satpw(ims:ime)
    real :: cloudbase(ims:ime)
    real :: pres1(ims:ime)
    real :: hgt1 (ims:ime)
    real :: tsfc1(ims:ime)
    real :: tminustd(nvl,ims:ime)

    integer :: k     ! index over model levels
    integer :: ksig
    integer :: ipn   ! horizontal index
    integer :: nv    ! index over tracers
    integer :: kk    ! index over (interpolated) pressure levels
    integer :: ivar
    integer :: lb
    integer :: iter
    integer :: ret   ! return code from GPTL calls
    logical :: vrbos
    character :: text*6

    real :: drh
    real :: du
    real :: dv
    real :: dx
    real :: dtemp,dq,dqw,doz
    real :: rh
    real :: esw
    real :: qsw
    real :: es
    real :: esln
    real :: thet
    real :: dthet
    REAL :: dpkap, dpkapsig
    real :: dhgt
    real :: dvv
    real :: gam
    real :: exls, exlsinv
    real :: tt1
    real :: t6
    real :: thbar
    real :: zcldbase
    real :: pcldbase
    real :: zcldtop
    real :: pcldtop
    real :: watericemax
    real :: smoo_coef(nvlp)
    integer :: nsm_var0(6), nsm_var(6)
    integer :: k700,k850		! index on isobaric level list for 700mb and 850mb
    real :: wsp80
!   nvarp = 7 for Z, T, RH, u, v, vert.vel., cloud/hydrometeor condensate 
!                 --------------------       
    data nsm_var0/5 ,3 , 2 ,5 ,5 , 3      /  ! no smooth for condensate

    real (kind=kind_rad) :: tk

! for 80 m wind
    integer nl
    real u1,u2,v1,v2,fact,u80,v80

    real arg,ex2p,p2ex
    ex2p(arg) = p1000*(arg/cp)**(cp/rd)   !  convert Pi => p
    p2ex(arg) = cp*(arg/p1000)**(rd/cp)   !  convert p  => Pi

    ret = gptlstart ('op_diag')

! Calculate grid length, used later to calculate no. of smoothing passes
      dx = sqrt(5.e14/float(nip))		! 27,500m (27.5km) at g8
!     smoo_coef = 0.2*dx			! 5.5km at g8

! Calculate isobaric levels every 25 hPa
    k700=-999
    k850=-999
    do kk=1,nvlp
! prsp(kk) = p1000 -(float(kk-1) * 2500.)     ! in Pa
      prsp(kk) = float(pres_pa(kk))           ! in Pa
!   pkapp(kk) = (1.- float(kk-1)* 0.025)**(rd/cp)
      pkapp(kk) = (prsp(kk)/100000.)**(rd/cp)
!      print *, "in op_diag calc levels: prsp(", kk, "): ", prsp(kk), " pkapp(",kk,"): ", pkapp(kk)
      smoo_coef(kk) = 0.2*dx
      if (abs(prsp(kk)-70000.).le. .01) k700=kk
      if (abs(prsp(kk)-85000.).le. .01) k850=kk
    end do
    if (min(k700,k850).le.0) stop 'wrong value for k700/k850'
   
      print *, "in op_diag calc : dx=(",dx,"), smoo_coef=(", smoo_coef(1), ") "
      do ivar=1,6
        nsm_var(ivar) = max(1,nint(15000./dx*float(nsm_var0(ivar))) )
      end do

      print *,' No. of smoothing passes for Z/T/RH/u/v/vv =',(nsm_var(kk),kk=1,6)

    exls    = rd*gams/grvity
    exlsinv = 1./exls

! Threaded loop
!sms$ignore begin
!$OMP PARALLEL DO PRIVATE (vrbos,k,kk,nv,pkap,pres,hgt,tvirt,thv,thnv,tk,esw,qsw,es,esln,nl,pklv,&
!$OMP                      tt1,t6,gam,tsfc,dthet,dpkap,drh,du,dv,dvv,thet,rh,dhgt,	&
!$OMP                      thbar,watericemax,zcldbase,zcldtop,pcldbase,pcldtop,ksig,	&
!$OMP                      pkapsig,dtemp,dq,dqw,doz,fldlyr,php,thp,fldp,    &
!$OMP                      u1,u2,v1,v2,fact,u80,v80,wsp80,dpkapsig)	&
!$OMP             SCHEDULE (runtime)
    do ipn=ips,ipe
      vrbos=ipn.eq.PrintIpnDiag

! --- check consistency of pressure, exner fct, pot.temp., geopotential

#ifdef DEBUGPRINT
      if (vrbos) then
        do k=1,nvl
          print 101,its,'  ipn,k=',perm(ipn),k,				&
            '  theta, d(phi)/d(ex)',tr(k,ipn,1),			&
            (ph3d(k+1,ipn)-ph3d(k,ipn))/(ex3d(k,ipn)-ex3d(k+1,ipn))
        end do
 101    format ('(op_diag)',i6,a,i8,i3,a,2f11.3)
      end if
#endif
      pw2d(ipn)  = 0.
      pq2d(ipn)  = 0.
      satpw(ipn) = 0.
!     qv3d(:,ipn)=tr(:,ipn,2)
!     qw3d(:,ipn)=tr(:,ipn,3)
!     oz3d(:,ipn)=tr(:,ipn,4)

!  Initialize isobaric fields as -999999. = -spval_p
      g3p     (:,ipn,:) = -spval_p
      g3p_chem(:,ipn,:) = 0.
      g3p     (:,ipn,6) = 0.      
      g3p     (:,ipn,7) = 0.      
      g2d     (ipn,:)   = -spval_p
      fldp    (:,:)   = -spval_p
!  Initialize sigma fields as -999999. = -spval_p
      g3sig   (:,ipn,:) = -spval_p

      do k=1,nvl
!   PW calculation here uses water vapor only as of May 2014 (water vapor and condensate before)
!       pw2d(ipn) = pw2d(ipn) + dp3d(k,ipn)*(tr(k,ipn,2) + tr(k,ipn,3))/grvity
!   PQ calculation is vertical integration of hydrometeor condensate
        pq2d(ipn) = pq2d(ipn) + dp3d(k,ipn)*(tr(k,ipn,3))/grvity
        pw2d(ipn) = pw2d(ipn) + dp3d(k,ipn)*(tr(k,ipn,2))/grvity

!rb   pres(k) = 0.5*(pr3d(k,ipn)+pr3d(k+1,ipn))
!rb   hgt1(k)  = 0.5*(ph3d(k,ipn)+ph3d(k+1,ipn)) / grav98
!rb   pkap(k) = (pres(k)/p1000)**(rd/cp)

        if (pr3d(k,ipn) > pr3d(k+1,ipn) + 0.1) then  ! 0.1 - make sure pr3d decreases monotonically
          pkap(k) = (ex3d(k  ,ipn)*pr3d(k  ,ipn)                     &
                    -ex3d(k+1,ipn)*pr3d(k+1,ipn))/                   &
                    ((cp+rd)*(pr3d(k,ipn) - pr3d(k+1,ipn)))
        else
          pkap(k) = 0.5*(ex3d(k,ipn) + ex3d(k+1,ipn))/cp
        end if
        pres(k)  = p1000*pkap(k)**(cp/rd)
        hgt(k)   = (ph3d(k,ipn) + (ex3d(k,ipn) - cp*pkap(k))*tr(k,ipn,1)) / grav98
        if (k.eq.1) then
          pres1(ipn) = pres(k)
          hgt1 (ipn) = hgt (k)
        end if
        tvirt(k) = tr(k,ipn,1) * pkap(k)                    ! virtual temp
!   temperature
        tk3d(k,ipn) = tvirt(k)/(1. + rovrm1_p*tr(k,ipn,2))    ! temperature, rovrm1_p=0.6078
        thv(k)      = tr(k,ipn,1)                           ! virtual pot temp
        thnv(k)     = tr(k,ipn,1)/(1. + rovrm1_p*tr(k,ipn,2)) ! non-virt pot temp 
        tk = tk3d(k,ipn)
!  sat vapor pressure w.r.t. water/ice combination
!     esw=fpvs(tk)
!  sat vapor pressure w.r.t. water (liquid) only
        esw = fpvsl(tk)
        esw = min (real(fpvsl(tk)) , pres(k)) ! qsw <= 1, as it has to be
!   Calculation of mixing ratio from vapor pressure
        qsw = ratio_h20_dry*esw/(pres(k) - (1.-ratio_h20_dry)*esw)
        qsw = max(1.e-8,min(qsw,0.1))
        dummy(k,ipn) = qsw
        rh3d(k,ipn) = 100.*max(0., min(1., tr(k,ipn,2)/qsw))
! rh for 0-100, as done in GRIB for other models
!   PW in saturated column
        satpw(ipn) = satpw(ipn)+(dp3d(k,ipn)*qsw/grvity)
        es    = pres(k)*(tr(k,ipn,2) + 1.e-8)/(ratio_h20_dry+(tr(k,ipn,2) + 1.e-8))
        esln  = log(es)
        td3d(k,ipn) = (35.86*esln - 4947.2325)/(esln - 23.6837)   ! Teten's equation
!  http://storms.meas.ncsu.edu/users/mdparker/courses/mea712.2006/Teten_s_equation.pdf
        tminustd(k,ipn) = tk3d(k,ipn) - td3d(k,ipn)
      end do

      rhpw(ipn) = pw2d(ipn) / satpw(ipn)

! Compute 80 m wind
      do k=2,nvl
         if(hgt(k)-ph3d(1,ipn)/grav98 .ge. 80.)then
           nl=k-1 ! level right below 80 m
           go to 50
         else
           nl=0
         endif
      enddo
   50 continue
       if(nl .ge. 1) then
            u1 = us3d(nl,ipn)
            u2 = us3d(nl+1,ipn)
            v1 = vs3d(nl,ipn)
            v2 = vs3d(nl+1,ipn)
         fact=(80.+ph3d(1,ipn)/grav98-hgt(nl))/(hgt(nl+1)-hgt(nl))
       else
            u1 = 0.
            u2 = us3d(1,ipn)
            v1 = 0.
            v2 = vs3d(1,ipn)
         fact = 80./hgt(1)
       endif
!  if(fact.gt.1.) then
!         print*,'u1,u2,v1,v2,fact,nl,hgt(nl),hgt(nl+1),u80,v80,wsp80' &
!                ,u1,u2,v1,v2,fact,nl,hgt(nl),hgt(nl+1),u80,v80,wsp80,ipn
!  endif
       fact = min(1.,fact)
       u80 = u1 + (u2-u1)*fact
       v80 = v1 + (v2-v1)*fact
! 80 m wind speed
       wsp80      = sqrt(u80**2 + v80**2)
!  if(wsp80.gt.100.) then
!         print*,'u1,u2,v1,v2,fact,nl,hgt(nl),hgt(nl+1),u80,v80,wsp80' &
!                ,u1,u2,v1,v2,fact,nl,hgt(nl),hgt(nl+1),u80,v80,wsp80,ipn
!  endif

!   Exner-like variables
!   ===================
!  pkap = (p/p0)**(R/Cp)  i.e., p(in bars) to the kappa (Rd/Cp) power
!  pkap = Exner / Cp
!   ===================
!   pklv  - Exner fn / Cp (0-1) on nvl+1 LEVELS (not layer midpoints)
!   pkap  -  "      "        "     nvl   LAYER midpoints
!   pkapp -  "      "        "  on nvlp  ISOBARIC levels

!   temperature variables
!   ===================
!   th3d, thv - virtual  potential temperature on nvl LAYER midpoints
!   thnv      - NON-virt     "         "        " "   "      "
!   tvirt     - virtual           temperature on nvl LAYER midpoint
!   tk3d      - nonvirtual           "        "   "   "      "

      do k=1,nvlp1
        pklv(k) = ex3d(k,ipn)/cp
      end do

#ifdef DEBUGPRINT
      if (vrbos) then
        do k=1,nvl
          print 1002,its,			&
           ' pkap(k), pklv(k), ex3d(k,ipn)',	&
             pkap(k), pklv(k), ex3d(k,ipn)
        end do
1002    format ('(op_diag)',i6,a,3f11.3)
      end if
#endif

!............................................................
!   Calculate isobaric values
!............................................................

!    -- Step 1
!   Calculate lapse rate (in virtual temp!) near surface:
!      Level 6 should be about 60 hPa above surface.
!      To avoid excessive effect (warm or cold) from lowest level,
!        use theta(k=3) with pres(k=1) to obtain
!        estimated temp at level 1 (also used in RUC post - hb2p.f)

      tt1  = thnv(3)*pkap(1)
      t6   = thnv(6)*pkap(6)
      gam  = (tt1-t6)/(hgt(6)-hgt(1))
      gam  = min (gamd, max(gam,gams))
      tsfc = thnv(3)*(pr3d(1,ipn)/p1000)**(rd/cp)  ! tsfc only for reduction here
      tsfc1(ipn) = tsfc

!   tsfc      - virtual temp at lowest LEVEL
!   tt1       - virtual temp at lowest LAYER midpoint

! Step 2 - obtain isobaric variables IF terrain elevation allows

      do k=1,nvl
       fldlyr(k,1)=thnv(k)		! non-virtual pot.temperature
      end do

! --------------------------------------------------------
! Obtain isobaric height and virtual pot temperature
      call lyr2prs_th_gz(its,nvl, ex3d(:,ipn),ph3d(:,ipn),fldlyr,	&
                             nvlp,prsp,       php,        fldp,         &
                             vrbos,ipn)

      do k=1,nvlp
        if (g3p(k,ipn,1).ne.spval_p) then
          g3p(k,ipn,1) = php(k)/grav98			! geopotential meters
          g3p(k,ipn,2) = fldp(k,1)*pkapp(k)		! virt.temperature (K)
        end if     ! spval_p
      end do

! --------------------------------------------------------
! Obtain isobaric non-chem scalar fields
      do k=1,nvl
!      fldlyr(k,1)=rh3d(k,ipn)		! rel.humidity
       fldlyr(k,1)=qv3d(k,ipn)		! spe.humidity
       fldlyr(k,2)=us3d(k,ipn)		! u
       fldlyr(k,3)=vs3d(k,ipn)		! v
       fldlyr(k,4)=ws3d(k,ipn)		! omega
       fldlyr(k,5)=tr  (k,ipn,3)	! cloud/hydrometeor condensate (g/g)
       fldlyr(k,6)=thnv(k)		! non-virt pot temp
      end do

      call lyr2prs(its,nvl, ex3d(:,ipn),nscalvar,fldlyr,	&
                       nvlp,prsp,                fldp,		&
                       vrbos,ipn)

      do k=1,nvlp
        if (g3p(k,ipn,1).ne.spval_p) then
          g3p(k,ipn,2) = fldp(k,6)*pkapp(k)	! non-virtual temperature
          g3p(k,ipn,3) = fldp(k,1)		! spe.humidity
          g3p(k,ipn,4) = fldp(k,2)		! u
          g3p(k,ipn,5) = fldp(k,3)		! v
          g3p(k,ipn,6) = fldp(k,4)		! omega
          g3p(k,ipn,7) = fldp(k,5)		! cloud/hydrometeor condensate (g/g)
        end if     ! spval_p
      end do


! --------------------------------------------------------
! Obtain isobaric chem scalar fields
!SB - uncomment below if isobaric output of chem vars needed
     if(ntrb > 0 ) then
       do nv=1,ntrb
         do k=1,nvl
           fldlyr(k,nv)=tr(k,ipn,nv+ntra)
         enddo
       enddo

     call lyr2prs(its,nvl, ex3d(:,ipn),ntrb    ,fldlyr,	&
                      nvlp,prsp,                fldp,		&
                      vrbos,ipn)

       do nv=1,ntrb
         do k=1,nvlp
           if (g3p(k,ipn,1).ne.-spval_p) then
             g3p_chem(k,ipn,nv)=fldp(k,nv)
           endif
         enddo
       enddo
     endif




! Step 3 - REDUCE from atmosphere above to obtain heights/temps below terrain surface
!   Set other isobaric variables also below terrain surface
      do kk=nvlp-1,1,-1
        if (g3p(kk,ipn,6).eq.spval_p) g3p(kk,ipn,6) = 0.
        if (g3p(kk,ipn,7).eq.spval_p) g3p(kk,ipn,7) = 0.
        if (pkap(1) < pkapp(kk)) then
          g3p(kk,ipn,2) = tt1*(prsp(kk)/pres(1))**exls
          g3p(kk,ipn,1) = ph3d(1,ipn)/grav98 - (g3p(kk,ipn,2)-tsfc)/gams
          g3p(kk,ipn,3) = g3p(kk+1,ipn,3)
          g3p(kk,ipn,4) = g3p(kk+1,ipn,4)
          g3p(kk,ipn,5) = g3p(kk+1,ipn,5)
          g3p(kk,ipn,6) = 0.
          g3p(kk,ipn,7) = 0.
          if (ntrb > 0) then
            do nv=1,ntrb
              if (g3p_chem(kk,ipn,nv).eq.spval_p) g3p_chem(kk,ipn,nv) = 0.
              g3p_chem(kk,ipn,nv) = max(0.,g3p_chem(kk+1,ipn,nv))
            end do
          end if
        end if
      end do

!     do kk = nvlp-3, nvlp
! Step 5 - EXTRAPOLATE to obtain height/temp above top native level
!     Set other isobaric variables also above top native level
!       if (g3p(kk,ipn,2) < 10.) then
!          g3p(kk,ipn,2) = tk3d(nvl,ipn)     !  isothermal lapse rate
!          g3p(kk,ipn,1) = hgt(nvl) + rd*tk3d(nvl,ipn)/grvity &
!                  * alog(pres(nvl)/prsp(kk))
!          g3p(kk,ipn,3) = g3p(kk-1,ipn,3)
!          g3p(kk,ipn,4) = g3p(kk-1,ipn,4)
!          g3p(kk,ipn,5) = g3p(kk-1,ipn,5)
!         if(ntrb.gt.0) then
!          do nv=ntra+2,ntra+ntrb
!            g3p(kk,ipn,nv)    = g3p(kk-1,ipn,nv)
!          enddo
!         endif
!       end if
!     end do


! ================================================================
! ================================================================
!   Calculate values of FIM data on GFS sigma levels
!............................................................
!   If FIM is running in pure sigma coordinate mode, then simply
!     copy FIM prognostic fields into the g3sig array.
! -------------------------------------------
      if (pure_sig .and. nvl.eq.nvlsig) then
        do k=1,nvlp1
          g3sig(k,ipn,1) = pr3d(k,ipn)
        end do
        do k=1,nvl
          g3sig(k,ipn,2) = dp3d(k,ipn)
          g3sig(k,ipn,3) = tk3d(k,ipn)
          g3sig(k,ipn,4) = us3d(k,ipn)
          g3sig(k,ipn,5) = vs3d(k,ipn)
          g3sig(k,ipn,6) = tr(k,ipn,2)
          g3sig(k,ipn,7) = tr(k,ipn,3)
          g3sig(k,ipn,8) = tr(k,ipn,4)
        end do

! -------------------------------------------
      else   !   pure_sig = .false. (using hybrid isentropic-sigma coordinate)
! -------------------------------------------

! Step 1 - obtain pressure at all levels
        g3sig(  1,ipn,1) = pr3d(1,ipn)
        g3sig(nvlsig+1,ipn,1) = pr3d(nvlp1,ipn)
        
!     do k=1,nvl
!      pkap(k) = (pres(k)/p1000)**(rd/cp)
!     end do

        do k = 2,nvlsig
          g3sig(k,ipn,1) = (sigak(k)+sigbk(k)*pr3d(  1,ipn))  ! pressure
          g3sig(k,ipn,2) = g3sig(k-1,ipn,1) - g3sig(k,ipn,1)  ! pressure thicknesses
          pkapsig(k) = (0.5*(g3sig(k-1,ipn,1)+g3sig(k,ipn,1))/100000.)**(rd/cp)
        end do
! Ensure that pkapsig at lowest level is within range of pkap values in that column
        pkapsig(2) = min(pkap(1)-0.00001,pkapsig(2))

        if (vrbos) then
          write (6,*) 'pkap vs pkapsig'
          do k=1,nvl
            write (6,*) pkap(k),pkapsig(k)
          end do
        end if

! Set all lowest sigma values = lowest hybrid values -
        g3sig(1  ,ipn,3) = tk3d(1,ipn)
        g3sig(1  ,ipn,4) = us3d(1,ipn)
        g3sig(1  ,ipn,5) = vs3d(1,ipn)
        g3sig(1  ,ipn,6) = tr(1,ipn,2)
        g3sig(1  ,ipn,7) = tr(1,ipn,3)
        g3sig(1  ,ipn,8) = tr(1,ipn,4)
! Set all top sigma values = top hybrid values - not exactly right
        g3sig(nvlsig,ipn,3) = tk3d(nvl,ipn)
        g3sig(nvlsig,ipn,4) = us3d(nvl,ipn)
        g3sig(nvlsig,ipn,5) = vs3d(nvl,ipn)
        g3sig(nvlsig,ipn,6) = tr(nvl,ipn,2)
        g3sig(nvlsig,ipn,7) = tr(nvl,ipn,3)
        g3sig(nvlsig,ipn,8) = tr(nvl,ipn,4)

! Step 2 - obtain sigma variables (except height) IF terrain elevation allows
        do ksig = 2,nvlsig
          do k=2,nvl
            dtemp = tk3d(k,ipn) - tk3d(k-1,ipn)
            dpkap  = pkap(k)  - pkap(k-1)
            dq    = tr(k,ipn,2) - tr(k-1,ipn,2)
            dqw   = tr(k,ipn,3) - tr(k-1,ipn,3)
            doz   = tr(k,ipn,4) - tr(k-1,ipn,4)
            du    = us3d(k,ipn) - us3d(k-1,ipn)
            dv    = vs3d(k,ipn) - vs3d(k-1,ipn)
            if (pkap(k) .lt. pkapsig(ksig) .and.pkap(k-1).ge.pkapsig(ksig) ) then
              g3sig(ksig,ipn,3) = tk3d (k-1,ipn) &
                   + (pkapsig(ksig)-pkap(k-1)) * dtemp/dpkap
!         g3sig(k,ipn,3) = tvirt (k-1) &
!                  + (pkapsig(ksig)-pkap(k-1)) * dtvirt/dpkap
!         rh    = rh3d(k-1,ipn) &
!                  + (pkapp(ivlp)-pkap(k-1)) * drh  /dpkap
              g3sig(ksig,ipn,4)    = us3d(k-1,ipn) &
                   + (pkapsig(ksig)-pkap(k-1)) * du   /dpkap
              g3sig(ksig,ipn,5)    = vs3d(k-1,ipn) &
                   + (pkapsig(ksig)-pkap(k-1)) * dv   /dpkap
              g3sig(ksig,ipn,6)    = tr(k-1,ipn,2) &
                   + (pkapsig(ksig)-pkap(k-1)) * dq   /dpkap
              g3sig(ksig,ipn,7)    = tr(k-1,ipn,3) &
                   + (pkapsig(ksig)-pkap(k-1)) * dqw  /dpkap
              g3sig(ksig,ipn,8)    = tr(k-1,ipn,4) &
                   + (pkapsig(ksig)-pkap(k-1)) * doz  /dpkap
!         if(ntrb.gt.0)then
!           do nv=ntra+1,ntra+ntrb
!              dv    = tr(k,ipn,nv) - tr(k-1,ipn,nv)
!              g3sig(ksig,ipn,6+nv-ntra-1)    = tr(k-1,ipn,nv) &
!                    + (pkapsig(ksig)-pkap(k-1)) * dv   /dpkap
!           enddo
!         endif  ! ntrb
            end if
          end do  ! k

#ifdef DEBUGPRINT
          if (vrbos.and.ksig.eq.2) then
!!      do k=1,nvlp1
!!        print 101,its,'  ipn,k=',ipn,k,'  pr => ex',			&
!!          ex3d(k,ipn),p2ex(pr3d(k,ipn))
!!        print 101,its,'  ipn,k=',ipn,k,'  ex => pr',			&
!!          pr3d(k,ipn),ex2p(ex3d(k,ipn))
!!      end do
            print 1001,its,	&
                 'dpkap,pkapsig(ksig), pkap(k-1)', &
                 dpkap,pkapsig(ksig), pkap(k-1)
1001        format ('(op_diag)',i6,a,3f11.3)
          end if
#endif

        end do    ! ksig
      end if      ! pure_sig

!  Calculate reverse weights for interpolation from GFS sig back to nativeFIM levels

      IF (enkfio_out) THEN
        IF (.NOT. pure_sig) THEN
          wt_int_rev(1,ipn) = 0.
          k_int_rev (1,ipn) = 1
          
          wt_int_rev(nvl,ipn) = 1.
          k_int_rev(nvl,ipn) = nvlsig-1
          
          DO k=2,nvl-1
            ksig=2
            IF (pkapsig(ksig) < pkap(k)) THEN
              wt_int_rev(k,ipn) = 1.
              k_int_rev(k,ipn) = ksig-1
            ENDIF

            DO ksig = 3,nvlsig
              dpkapsig  = pkapsig(ksig)  - pkapsig(ksig-1)
              IF (pkapsig(ksig) < pkap(k) .AND. &
                  pkapsig(ksig-1) >= pkap(k) ) THEN
                wt_int_rev(k,ipn) = (pkap(k)-pkapsig(ksig-1))/dpkapsig
                k_int_rev(k,ipn) = ksig-1
              END IF
            END DO
          END DO
        ENDIF  ! if .not. pure_sig
      ENDIF    ! if enkfio_out
! ================================================================
!  End of calculating for output on GFS vertical levels
! ================================================================
! ================================================================

!===============================================================
!  Cloud base/top diagnosis
!===============================================================

! g2d(ipn,1) = -spval_p
      watericemax = 0.
      zcldbase    = -spval_p
      zcldtop     = -spval_p
      pcldbase    = 0.
      pcldtop     = 0.

      do k=1,nvl
        watericemax = max(watericemax, tr(k,ipn,3) )
      end do

      if (watericemax >= cloud_def_p) then
! At surface?
        if (tr(1,ipn,3) > cloud_def_p) then
          zcldbase = hgt (1)
          pcldbase = pres(1)
        else
! Cloud aloft?
          do k=3,nvl
            if (tr(k,ipn,3) > cloud_def_p) then
              zcldbase = hgt(k-1) + (tr(k,ipn,3)-cloud_def_p)     &
                         *(hgt(k  ) - hgt(k-1))                   &
                         / max(cloud_def_p,(tr (k, ipn,3)-tr(k-1,ipn,3)))
!             zcldbase = hgt(k)

              pcldbase = pres(k-1) + (tr(k,ipn,3)-cloud_def_p)  &
                         *(pres(k  ) - pres(k-1))               &
                         / max(cloud_def_p,(tr (k, ipn,3)-tr(k-1,ipn,3)))
!             pcldbase = pres(k)
              exit
            end if
          end do
        end if

        do k=nvl-2,2,-1
          if (tr(k,ipn,3) > cloud_def_p) then
            zcldtop = hgt(k) + (tr(k,ipn,3)-cloud_def_p)  &
                     *(hgt(k+1) - hgt(k))             &
                     / max(cloud_def_p,(tr(k, ipn,3)-tr(k+1,ipn,3)))
!           zcldtop = hgt(k)
            pcldtop = pres(k) + (tr(k,ipn,3)-cloud_def_p)  &
                     *(pres(k+1) - pres(k))             &
                     / max(cloud_def_p,(tr(k,ipn,3)-tr(k+1,ipn,3)))
!           pcldtop = pres(k)
            exit
          end if
        end do
      end if

      cloudbase(ipn) = zcldbase
      g2d(ipn,1) = zcldbase
      g2d(ipn,2) = pcldbase
      g2d(ipn,3) = zcldtop
      g2d(ipn,4) = pcldtop
      g2d(ipn,5) = rhpw(ipn)
      g2d(ipn,6) = satpw(ipn)
      g2d(ipn,7) = wsp80
      t2m_dif(ipn) = t2m2d(ipn) - tk3d(1,ipn)
!===============================================================
!  Snow accumulation diagnostic - purely based on 850 hPa temp and k=1 temp
!     Criterion:  (850 T < 0 deg C) .and. (k=1 T < 0.84 deg C)
!     NOTE:  Do not use 850 T check unless psfc > 870 hPa (small pressure
!         spacing from 850 hPa)
!==============================================================
      sn_xh(ipn)  = 0.
      sn_6h(ipn)  = 0.
      if (pr3d(1,ipn).gt.87000.) then
        if (g3p(k850,ipn,2) < 273.15 .and. tk3d(1,ipn) < 274.) then
          sn_xh(ipn) = rn_xh(ipn)
          sn_6h(ipn) = rn_6h(ipn)
        end if
      else
        if (tk3d(1,ipn) < 273.15 ) then
          sn_xh(ipn) = rn_xh(ipn)
          sn_6h(ipn) = rn_6h(ipn)
        end if
      end if

    end do     ! ipn loop - primary horizontal loop

!$OMP END PARALLEL DO
!sms$ignore end

#ifdef DEBUGPRINT
    do nv=1,4
      do k=1,nvlp
       write (text,'(i1,a,i2)') nv,' k=',k
       call findmxmn3(g3p,nvlp,nip,nvarp,k,nv,'(op_diag) g3p'//text)
      end do
      print *
    end do
#endif

! -----------------------------------------------------------------------
!   max/min/mean values for FIM data on GFS vertical levels
! -----------------------------------------------------------------------
  if (atmonly) then	! to reduce printout in the coupled run
!sms$ignore begin
    LB = LBOUND(g3sig,2)
!sms$ignore end
    do ivar=1,nvarsig
      call outqv_mn(g3sig(1,LB,ivar),nvlsig+1,                1.,varsig(ivar) // '-sigma')
      call outqv   (g3sig(1,LB,ivar),nvlsig+1,deg_lat,deg_lon,1.,varsig(ivar) // '-sigma  - max/min')
    end do

    call outqv_mn (tminustd, nvl,                   1., 'T-Td')
    call outqv    (tminustd, nvl, deg_lat, deg_lon, 1., 'T-Td  - max/min')
  end if

! -----------------------------------------------------------------------
! Smooth isobaric variables (Z/T/RH/u/v/vertical velocity) at all levels.
! -----------------------------------------------------------------------
!sms$ignore begin
    LB = LBOUND(g3p,2)
!sms$ignore end
    if (TimingBarriers) then
      ret = gptlstart ('op_diag_barrier')
!SMS$BARRIER
      ret = gptlstop ('op_diag_barrier')
    end if

    do ivar = 1,6   ! only covering vars 1-6 of nvarp

    do iter=1,nsm_var(ivar)
     ret = gptlstart ('op_diag_exchange')
!SMS$EXCHANGE(g3p(:,:,ivar))
     ret = gptlstop  ('op_diag_exchange')
     ret = gptlstart ('dffusn_lev')
     call dffusn_lev (g3p(:,:,ivar), smoo_coef ,nvlp, 1, nvlp)
     ret = gptlstop  ('dffusn_lev')
    end do
    end do

  if (atmonly) then	! to reduce printout in the coupled run
    write (6,'(a)') ' After height smoothing'
    ivar = 1
    call outqv_mn (g3p(1,LB,ivar), nvlp,                   1., varname(ivar)//'-isobaric')
    call outqv    (g3p(1,LB,ivar), nvlp, deg_lat, deg_lon, 1., varname(ivar)//'-isobaric  - max/min')
  end if

!sms$ignore begin
    do ipn=ips,ipe
! -----------------------------------------------------------------------
! Calculate MAPS sea-level pressure reduction  (Benjamin and Miller, 1992, Mon. Wea. Rev.)
!               g3p(13,ipn,2) = temperature at 700 hPa (13th level starting at 1000 hPa every 25 hPa)
! -----------------------------------------------------------------------

      tsfc7= g3p(k700,ipn,2)* (pres1(ipn)/70000.)**exls ! est. temp at sfc extrapolating down (up) from 700 hPa
      mslp(ipn) = pres1(ipn)*((tsfc7+gams*hgt1(ipn))/tsfc7)**exlsinv

    end do     ! ipn loop - primary horizontal loop
!sms$ignore end
      
  if (atmonly) then	! to reduce printout in the coupled run
! -----------------------------------------------------------------------
!  Print out max/min/mean values for different fields
! -----------------------------------------------------------------------
    write(6,*)'op_diag: Beginning outqv* calls for time step=', its
    LB = LBOUND(g2d,1)
    call outqv_mn (mslp,      1,                     0.01, 'MAPS SLP')
    call outqv    (mslp,      1,   deg_lat, deg_lon, 0.01, 'MAPS SLP  - max/min')
    call outqv_mn (pw2d,      1,                     1.,   'Precipitable water')
    call outqv    (pw2d,      1,   deg_lat, deg_lon, 1.,   'Precipitable water  - max/min')
    call outqv_mn (pq2d,      1,                     1.,   'Integrated condensate')
    call outqv    (pq2d,      1,   deg_lat, deg_lon, 1.,   'Integrated condensate  - max/min')
    call outqv_mn_lat (pw2d,  1,   deg_lat, deg_lon,30.,1.,'Prec water - lat/lon')
    call outqv_mn_lat (pq2d,  1,   deg_lat, deg_lon,30.,1.,'Integ condensate - lat/lon')
    call outqv_mn (g2d(LB,6), 1,                     1.,   'SatPW')
    call outqv    (g2d(LB,6), 1,   deg_lat, deg_lon, 1.,   'SatPW  - max/min')
    call outqv_mn (g2d(LB,5), 1,                     1.,   'RH w.r.t. PW')
    call outqv    (g2d(LB,5), 1,   deg_lat, deg_lon, 1.,   'RH w.r.t. PW  - max/min')
    call outqv_mn_lat (g2d(LB,5), 1,deg_lat,deg_lon,30.,1.,'RH w.r.t. PW')
    call outqv_mn (tk3d,      nvl,                   1.,   'Temperature')
    call outqv    (tk3d,      nvl, deg_lat, deg_lon, 1.,   'Temperature  - max/min')
    call outqv_mn (rh3d,      nvl,                   1.,   'RH')
    call outqv    (rh3d,      nvl, deg_lat, deg_lon, 1.,   'RH  - max/min')
    call outqv_mn (dummy,     nvl,                   1.,   'qsw')
    call outqv    (dummy,     nvl, deg_lat, deg_lon, 1.,   'qsw  - max/min')
    call outqv    (cloudbase, 1,   deg_lat, deg_lon, 1.,   'Cloudbase  - max/min')
    call outqv    (g2d(LB,1) ,1,   deg_lat, deg_lon, 1.,   'Cloud base height  - max/min')
    call outqv    (g2d(LB,2) ,1,   deg_lat, deg_lon, 1.,   'Cloud base pressure  - max/min')
    call outqv    (g2d(LB,3) ,1,   deg_lat, deg_lon, 1.,   'Cloud top height  - max/min')
    call outqv    (g2d(LB,4) ,1,   deg_lat, deg_lon, 1.,   'Cloud top pressure  - max/min')
!   call outqv_mn (sm3d,      4,                     1.,   'Soil moisture')
!   call outqv    (sm3d,      4,   deg_lat, deg_lon, 1.,   'Soil moisture  - max/min')

    call outqv_mn (hf2d,      1,                     1.,         'Sensible heat flux')
    call outqv_mn_lat (hf2d,  1,   deg_lat, deg_lon, 30., 1000., 'Sensible heat flux')

    call outqv_mn (qf2d,      1,                     1.,         ' Latent heat flux')
    call outqv_mn_lat (qf2d,  1,   deg_lat, deg_lon, 30., 1000., ' Latent heat flux')

! -----------------------------------------------------------------------
!   max/min/mean values for isobaric grids
! -----------------------------------------------------------------------
!sms$ignore begin
    LB = LBOUND(g3p,2)
!sms$ignore end
    do ivar=1,size(varname)
      call outqv_mn (g3p(1,LB,ivar), nvlp,                   1., varname(ivar) // '-isobaric')
      call outqv    (g3p(1,LB,ivar), nvlp, deg_lat, deg_lon, 1., varname(ivar) // '-isobaric  - max/min')
    end do

  end if

    ret = gptlstop ('op_diag')

    return
  end subroutine op_diag


  subroutine lyr2prs_th_gz(its,kkinp,exninp,gzinp,fldinp,	&
                               kkout,prsout,gzout,thout,	&
                           vrbos,ipn)

! --- interpolate layer data (single column) to prescribed pressure levels
!     -- For height and potential temperature only

! use module_constants,only: rd,cp,p1000,perm

  integer, intent(IN) :: its		! model time step
  integer, intent(IN) :: kkinp		! # of layers in input column
  integer, intent(IN) :: kkout		! # of output levels
  integer, intent(IN) :: ipn		! horiz.grid location
  logical, intent(IN) :: vrbos

  real   , intent(IN) :: prsout(kkout)
  real   , intent(IN) :: exninp(kkinp+1),gzinp(kkinp+1),	&
                         fldinp(kkinp)
  real   ,intent(OUT) :: gzout(kkout)
  real thinp(kkinp),thout(kkout)
  real pkinp(kkinp+1),pkout(kkout)		! (p/p0)**(R/cp)
  real var,thav,pk1,pk2
  integer ki,ko,ka,kb,n

   pkinp(:)=exninp(:)/cp
   pkout(:)=(prsout(:)/p1000)**(rd/cp)
! --- derive hydrostatically constrained pot.temperature
   do ki=1,kkinp
    thinp(ki)=(gzinp(ki+1)-gzinp(ki))/(exninp(ki)-exninp(ki+1))
   end do

   do ko=1,kkout	! loop over output levels
    gzout(ko)   =spval_p
    do ki=1,kkinp	! loop over input layers
     pk1=pkinp(ki  )
     pk2=pkinp(ki+1)
     if (pkout(ko).le.pk1 .and. pkout(ko).ge.pk2) then
      ka=max(    1,ki-1)
      kb=min(kkinp,ki+1)

! --- interpolate hydrostatically constrained pot.temperature
      var=0.5*min(thinp(kb)-thinp(ki),				&
                 thinp(ki)-thinp(ka))
! --- assume linear variation within layer (as in PLM)
      thout(ko)=((thinp(ki)-var)*(pkout(ko)-pk2)		&
                +(thinp(ki)+var)*(pk1-pkout(ko)))/(pk1-pk2)

! --- interpolate geopotential
      thav=0.5*(thinp(ki)-var+thout(ko))
!     thav=thinp(ki)				! traditional method
      gzout(ko)=gzinp(ki)+thav*cp*(pk1-pkout(ko))
     end if
    end do  !   ki
   end do   !   ko

  return
  end subroutine lyr2prs_th_gz

  subroutine lyr2prs(its,kkinp,exninp,nfields,fldinp,	&
                         kkout,prsout,        fldout,	&
                         vrbos,ipn)

! --- interpolate layer data (single column) to prescribed pressure levels

! use module_constants,only: rd,cp,p1000,perm

  integer, intent(IN) :: its		! model time step
  integer, intent(IN) :: kkinp		! # of layers in input column
  integer, intent(IN) :: kkout		! # of output levels
  integer, intent(IN) :: nfields 	! # of fields
  integer, intent(IN) :: ipn		! horiz.grid location
  logical, intent(IN) :: vrbos

  real   , intent(IN) :: prsout(kkout)
  real   , intent(IN) :: exninp(kkinp+1),	&
                         fldinp(kkinp,nfields)
  real   ,intent(OUT) :: fldout(kkout,nfields)

  real pkinp(kkinp+1),pkout(kkout)		! (p/p0)**(R/cp)
  real var,thav,pk1,pk2
  integer ki,ko,ka,kb,n

#ifdef DEBUGPRINT
   if (vrbos) then
!SMS$IGNORE BEGIN
    print 99,its,perm(ipn),deg_lat(ipn),deg_lon(ipn),		&
     '  lyr2prs   i n p u t   profile:',			&
     'intfc.exn     field1      field2      field3      field4      field5'
    print 100,(ki,exninp(ki),          				&
     (fldinp(ki,n),n=1,min(5,nfields)),ki=1,kkinp),kkinp+1,exninp(kkinp+1)
 99 format ('its,ipn=',i6,i8,'  lat/lon=',2f7.1,a/5x,a)
100 format (i4,f10.2,5es12.4)
!SMS$IGNORE END
   end if
#endif

   pkinp(:)=exninp(:)/cp
   pkout(:)=(prsout(:)/p1000)**(rd/cp)

   do ko=1,kkout	! loop over output levels
    fldout(ko,:)=spval_p
    do ki=1,kkinp	! loop over input layers
     pk1=pkinp(ki  )
     pk2=pkinp(ki+1)
     if (pkout(ko).le.pk1 .and. pkout(ko).ge.pk2) then
      ka=max(    1,ki-1)
      kb=min(kkinp,ki+1)

! --- interpolate fields
      do n=1,nfields
       if (fldinp(ka,n).le.fldinp(kb,n)) then
        var=0.5*min(max(0.,fldinp(kb,n)-fldinp(ki,n)),		&
                   max(0.,fldinp(ki,n)-fldinp(ka,n)))
       else
        var=0.5*max(min(0.,fldinp(kb,n)-fldinp(ki,n)),		&
                   min(0.,fldinp(ki,n)-fldinp(ka,n)))
       end if
! --- assume linear variation within layer (as in PLM)
       fldout(ko,n)=((fldinp(ki,n)-var)*(pkout(ko)-pk2)		&
                    +(fldinp(ki,n)+var)*(pk1-pkout(ko)))/(pk1-pk2)
      end do
! --- assume fldout_1 = theta. convert to T.
!     fldout(ko,1)=fldout(ko,1)*pkout(ko)
      exit
     end if
    end do		! loop over input layers
   end do		! loop over output levels

#ifdef DEBUGPRINT
   if (vrbos) then
!SMS$IGNORE BEGIN
    print 99,its,perm(ipn),deg_lat(ipn),deg_lon(ipn),		&
     '  lyr2prs   o u t p u t   profile:',			&
     'intfc.exn     field1      field2      field3      field4      field5'
    print 100,(ko,.01*prsout(ko),(fldout(ko,n),n=1,min(5,nfields)),ko=1,kkout)
!SMS$IGNORE END
   end if
#endif

  return
  end subroutine lyr2prs
end module module_op_diag
