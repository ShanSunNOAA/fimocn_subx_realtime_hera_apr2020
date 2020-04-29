! ----------------------------------------------------------------------
! --- k-profile vertical mixing model (HYCOM's KPP version, June 2008)
! ---   (ref: Large, McWilliams, Doney)
!
! --- iocnmx
! ---    0       no action aside from applying surface forcing
! ---    1       full-column kpp
! ---    2       partial column kpp
! ----------------------------------------------------------------------
module hycom_mxkpp

contains
   subroutine mxkppaij(nstep,temp,saln,u,v,dp,p,theta,			& !input
                       surflx,sswflx,salflx,ustar,ustarb,jerlov,	& !input
                       deg_lat,corio,akpar,dlt,i,vrbos,			& !input
                       zgrid,dift,difs,vcty,ghats,dpmixl,klist)		!output
! --- hycom version 1.0
! -------------------------------------------------------------
! --- kpp vertical diffusion, single j-row (part A)
! --- vertical coordinate is z negative below the ocean surface
!
! --- Large, W.C., J.C. McWilliams, and S.C. Doney, 1994: Oceanic
! --- vertical mixing: a review and a model with a nonlocal
! --- boundary layer paramterization. Rev. Geophys., 32, 363-403.
!
! --- quadratic interpolation and variable Cv from a presentation
! --- at the March 2003 CCSM Ocean Model Working Group Meeting
! --- on KPP Vertical Mixing by Gokhan Danabasoglu and Bill Large
! --- http://www.ccsm.ucar.edu/working_groups/Ocean/agendas/030320.html
! --- quadratic interpolation implemented here by 3-pt collocation,
! --- which is slightly different to the Danabasoglu/Large approach.
! -------------------------------------------------------------

   use module_constants,only: grvity
   use hycom_control   ,only: latdiw,iocnmx
   use fimnamelist     ,only: kdm
   use hycom_constants, only: onemm, onem, tencm, spcifh, epsil, thref
   use hycom_kpp_constants, only: difsiw, difmiw,			&
       jerlv0, locsig, shinst, dbdiff, cekman, epsilon,			&
       difm0, difs0, rrho0, dsfmax, bldmax, bldmin, vonk, hblflg,	&
       ricr, cv, c11, cg, nonloc, cmonob, vtc, dp0enh, bblkpp,		&
       betabl,betard,redfac
   use hycom_sigetc

  implicit none
  integer,intent(IN)    :: i,nstep,jerlov
  logical,intent(IN)    :: vrbos
  real*8 ,intent(IN)    :: theta(kdm), dp(kdm)
  real*8 ,intent(INOUT) :: temp(kdm), saln(kdm), u(kdm), v(kdm), p(kdm+1)
  real*8 ,intent(OUT)   :: zgrid(kdm+1),vcty(kdm+1),difs(kdm+1),	&
                           dift(kdm+1),ghats(kdm+1)
  real   ,intent(IN)    :: surflx,sswflx,salflx,ustar,ustarb,corio,	&
                           deg_lat,akpar,dlt
  real*8 ,intent(INOUT) :: dpmixl
  integer,intent(OUT)   :: klist

  real hmonob0,dpbl0,dpbbl0,mixflx0,buoflx0,bhtflx0     ! potential output fields
  real, parameter :: difmax = 9999.0e-4  !maximum diffusion/viscosity
  real, parameter :: dp0bbl =   20.0     !truncation dist. for bot. b.l.
  real, parameter :: ricrb  =    0.45    !critical bulk Ri for bot. b.l.
  real, parameter :: cv_max =    2.1     !maximum cv
  real, parameter :: cv_min =    1.7     !minimum cv
  real, parameter :: cv_bfq =  200.0     !cv scale factor

! local variables for kpp mixing
  real delta               ! fraction hbl lies beteen zgrid neighbors
  real zrefmn              ! nearsurface reference z, minimum
  real zref                ! nearsurface reference z
  real wref,qwref          ! nearsurface reference width,inverse
  real uref                ! nearsurface reference u
  real vref                ! nearsurface reference v
  real bref                ! nearsurface reference buoyancy
  real swfrac(kdm+1)       ! fractional surface shortwave radiation flux
  real shsq(kdm+1)         ! velocity shear squared
  real alfadt(kdm+1)       ! t contribution to density jump
  real betads(kdm+1)       ! s contribution to density jump
  real swfrml              ! fractional surface sw rad flux at ml base
  real ritop(kdm)          ! numerator of bulk richardson number
  real dbloc(kdm+1)        ! buoyancy jump across interface 
  real dvsq(kdm)           ! squared current shear for bulk richardson no.
  real zgridb(kdm+1)       ! zgrid for bottom boundary layer
  real hwide(kdm)          ! layer thicknesses in m (minimum 1mm)
  real dpmm(kdm)           !     max(onemm,dp(k,i,leapn))
  real qdpmm(kdm)          ! 1.0/max(onemm,dp(k,i,leapn))
  real case                ! 1 in case A; =0 in case B
  real hbl                 ! boundary layer depth
  real hbbl                ! bottom boundary layer depth
  real rib(3)              ! bulk richardson number
  real rrho                ! double diffusion parameter
  real diffdd              ! double diffusion diffusivity scale
  real prandtl             ! prandtl number
  real rigr                ! local richardson number
  real fri                 ! function of Rig for KPP shear instability
  real stable              ! = 1 in stable forcing; =0 in unstable
  real dkm1(3)             ! boundary layer diffusions at nbl-1 level
  real gat1(3)             ! shape functions at dnorm=1
  real dat1(3)             ! derivative of shape functions at dnorm=1
  real blmc(kdm+1,3)       ! boundary layer mixing coefficients
  real bblmc(kdm+1,3)      ! boundary layer mixing coefficients
  real wm                  ! momentum velocity scale
  real ws                  ! scalar velocity scale
  real dnorm               ! normalized depth
  real*8 tmn               ! time averaged SST
  real*8 smn               ! time averaged SSS
  real dsgdt               ! dsigdt(tmn,smn)
  real buoyfs              ! salinity  surface buoyancy (into atmos.)
  real buoyfl              ! total     surface buoyancy (into atmos.)
  real buoysw              ! shortwave surface buoyancy (into atmos.)
  real bfsfc               ! surface buoyancy forcing   (into atmos.)
  real bfbot               ! bottom buoyancy forcing
  real hekmanb             ! bottom ekman layer thickness
  real cormn4              ! = 4 x min. coriolis magnitude (at 4N, 4S)
  real dflsiw              ! lat.dep. internal wave diffusivity
  real dflmiw              ! lat.dep. internal wave viscosity
  real bfq                 ! buoyancy frequency
  real cvk                 ! ratio of buoyancy frequencies
  real ahbl,bhbl,chbl,dhbl ! coefficients for quadratic hbl calculation
  real hekman              ! Ekman depth
  real th3d(kdm)

  logical lhbl             ! safe to use quadratic hbl calculation

  integer nbl              ! layer containing boundary layer base
  integer nbbl             ! layer containing bottom boundary layer base
  integer kup2,kdn2,kup,kdn! bulk richardson number indices
  integer niter            ! iterations for semi-implicit soln. (2 recomended for KPP)

! --- local 1-d arrays for matrix solution
  real*8 :: u1do(kdm+1),u1dn(kdm+1),v1do(kdm+1),v1dn(kdm+1),t1do(kdm+1),    &
            t1dn(kdm+1),s1do(kdm+1),s1dn(kdm+1),diffm(kdm+1),difft(kdm+1),  &
            diffs(kdm+1),ghat(kdm+1),zm(kdm+1),hm(kdm),dzb(kdm)

! --- local 1-d arrays for iteration loops
      real*8 :: uold(kdm+1),vold(kdm+1)
      real*8 :: told(kdm+1),sold(kdm+1),thold(kdm+1)

! --- tridiagonal matrix solution arrays
      real*8 tri(kdm,0:1)       ! dt/dz/dz factors in trid. matrix
      real*8 tcu(kdm),        & ! upper coeff for (k-1) on k line of trid.matrix
             tcc(kdm),        & ! central ...     (k  ) ..
             tcl(kdm),        & ! lower .....     (k-1) ..
             rhs(kdm)           ! right-hand-side terms

      real*8 dtemp,dsaln,wq,wt,ratio,q,ghatflux,dvdzup,dvdzdn,viscp,difsp,   &
           diftp,f1,sigg,aa1,aa2,aa3,gm,gs,gt,dkmp2,dstar,hblmin,hblmax,     &
           sflux1,vtsq,vctyh,difsh,difth,zrefo,qspcifh,hbblmin,hbblmax,      &
           beta_b,beta_r,frac_b,frac_r,swfbqp,x0,x1,x2,y0,y1,y2

      integer k,ka,kb,nlayer,ksave,iter,jrlv,nblmin

      real*8,parameter :: qrinfy=1./0.7  ! (max.gradient richardson number)^-1

      cormn4 = 4.0e-5  !4 x min. coriolis magnitude (at 4N, 4S)

      th3d(:)=0.
      told(:)=0.
      sold(:)=0.
      thold(:)=0.
      uold(:)=0.
      vold(:)=0.
      zgrid(:)=0.
      vcty(:)=0.
      dift(:)=0.
      difs(:)=0.
      ghat(:)=0.

      diffm(:)=0.
      difft(:)=0.
      diffs(:)=0.
      hekman=ustar*(cekman*4.0)/max(cormn4,corio)		! Ekman depth
      if (latdiw) then

! ---   latitude dependent internal wave diffusion/viscosity
! ---   Gregg et. al. (2003): Reduced mixing from the breaking of
! ---   internal waves in equatorial waters, Nature 422 pp 513-515.
! ---   q is a quadratic fit to eqn 2 of Gregg, assuming N=N0 (lat=1-45).

        q      = max(1.0,min(45.0,abs(deg_lat)))
        q      = 0.0424 + 0.04*q - 26.e-05*q**2  !q=1 at 30degN
        dflsiw = q*difsiw  !difsiw is ref.value at 30degN
        dflmiw = q*difmiw  !difmiw is ref.value at 30degN
      else
! ---   constant internal wave diffusion/viscosity
        dflsiw =   difsiw
        dflmiw =   difmiw
      endif

! --- locate lowest substantial mass-containing layer.
      do k=1,kdm
        dpmm( k)  =max(onemm,dp(k))
        qdpmm(k)  =1.0/dpmm(k)
        p(k+1)=p(k)+dp(k)
      enddo
      do k=kdm,1,-1
        if (dpmm(k).gt.tencm) then
          exit
        endif
      enddo
      klist=max(k,2)  !always consider at least 2 layers

! --- forcing of t,s by surface fluxes. flux positive into ocean.
! --- shortwave flux penetration depends on kpar or jerlov water type.

  if (jerlv0.eq.0) then
    beta_r = 2.0/onem
    beta_b = akpar/onem
    beta_b = max( betabl(1), beta_b)  !time interp. beta_b can be -ve
    frac_b = max( 0.27, 0.695 - 5.7*onem*beta_b )
    frac_r = 1.0 - frac_b
  else
    jrlv   = jerlov
    beta_r = betard(jrlv)
    beta_b = betabl(jrlv)
    frac_r = redfac(jrlv)
    frac_b = 1.0 - frac_r
  endif
  qspcifh=1.0/spcifh

! --- evenly re-distribute the flux below the bottom
  k = klist
  if (-p(k+1)*beta_r.gt.-10.0) then
    swfbqp=frac_r*exp(-p(k+1)*beta_r)+frac_b*exp(-p(k+1)*beta_b)
  elseif (-p(k+1)*beta_b.gt.-10.0) then
    swfbqp=frac_b*exp(-p(k+1)*beta_b)
  else
    swfbqp=0.0
  endif
  swfbqp = swfbqp/p(k+1)

  if (vrbos) then
    write (*,'(a,4f10.4,i2)')  'frac[rb],beta[rb] =',	&
        frac_r,frac_b,onem*beta_r,onem*beta_b,jrlv
   endif

  do k=1,kdm
      if (-p(k+1)*beta_r.gt.-10.0) then
        swfrac(k+1)=frac_r*exp(-p(k+1)*beta_r)+		&
                    frac_b*exp(-p(k+1)*beta_b)
      elseif (-p(k+1)*beta_b.gt.-10.0) then
        swfrac(k+1)=frac_b*exp(-p(k+1)*beta_b)
      else
        swfrac(k+1)=0.0
      endif
      swfrac(k+1)=swfrac(k+1)-swfbqp*p(k+1)  !spread out bottom frac
      if (k.eq.1) then
        sflux1=surflx-sswflx
        dtemp=(sflux1+(1.-swfrac(k+1))*sswflx)*		&
               dlt*grvity*qspcifh*qdpmm(k)
        dsaln=salflx*dlt*grvity*qdpmm(k) 
        if (vrbos) then
          write (*,101) nstep,k,i,1.0,swfrac(k+1),dtemp,dsaln
         endif
       elseif (k.le.klist) then
         dtemp=(swfrac(k)-swfrac(k+1))*sswflx*dlt*grvity*qspcifh*qdpmm(k)
         dsaln=0.0
         if (vrbos) then
           write (*,101) nstep,k,i,swfrac(k),swfrac(k+1),dtemp
         endif
      else			! k.gt.klist
         dtemp=0.0
         dsaln=0.0
      endif 

! --- modify t and s; set old value arrays at p points for initial iteration
    if (k.le.klist) then
      temp(k)=temp(k)+dtemp
      saln(k)=saln(k)+dsaln
      th3d(k)=sigocn(temp(k),saln(k))
      told (k)=temp(k)
      sold (k)=saln(k)
      if (locsig) then
        if (k.eq.1) then
          thold(k)=th3d(k)
        else
          ka=k-1
          alfadt(k)=0.5*(dsiglocdt(told(ka),sold(ka),p(k))+	&
                         dsiglocdt(told(k ),sold(k ),p(k)))*	&
                         (told(ka)-told(k))
          betads(k)=0.5*(dsiglocds(told(ka),sold(ka),p(k))+	&
                         dsiglocds(told(k ),sold(k ),p(k)))*	&
                         (sold(ka)-sold(k))
          thold(k)=thold(ka)-alfadt(k)-betads(k)
        endif
      else
        thold(k)=th3d(k)
      endif
      uold (k)=u(k)
      vold (k)=v(k)
    endif
  enddo

  if (iocnmx.eq.0) return			! skip mixing

    k=klist
    ka=k+1
    kb=min(ka,kdm)
    told (ka)=temp(kb)
    sold (ka)=saln(kb)
    if (locsig) then
      alfadt(ka)=0.5*(dsiglocdt(told(k ),sold(k ),p(ka))+	&
                      dsiglocdt(told(ka),sold(ka),p(ka)))*	&
                       (told(k)-told(ka))
      betads(ka)=0.5*(dsiglocds(told(k ),sold(k ),p(ka))+	&
                      dsiglocds(told(ka),sold(ka),p(ka)))*	&
                       (sold(k)-sold(ka))
      thold(ka)=thold(k)-alfadt(ka)-betads(ka)
    else
      thold(ka)=th3d(kb)
    end if
    uold (ka)=u(k)
    vold (ka)=v(k)

! --- calculate z at vertical grid levels - this array is the z values in m
! --- at the mid-depth of each micom layer except for index klist+1, where it
! --- is the z value of the bottom

! --- calculate layer thicknesses in m
    do k=1,kdm
      if (k.eq.1) then
        hwide(k)=dpmm(k)/onem
        zgrid(k)=-.5*hwide(k)
      else if (k.lt.klist) then
        hwide(k)=dpmm(k)/onem
        zgrid(k)=zgrid(k-1)-.5*(hwide(k-1)+hwide(k))
      else if (k.eq.klist) then
        hwide(k)=dpmm(k)/onem
        zgrid(k)=zgrid(k-1)-.5*(hwide(k-1)+hwide(k))
        zgrid(k+1)=zgrid(k)-.5*hwide(k)
      else
        hwide(k)=0.
      endif
    enddo

! --- perform niter iterations to execute the semi-implicit solution of the
! --- diffusion equation. at least two iterations are recommended

    niter=2
    if (iocnmx.gt.4) niter=1
    do iter=1,niter

! --- calculate layer variables required to estimate bulk richardson number

! --- calculate nearsurface reference variables,
! --- averaged over -2*epsilon*zgrid, but no more than 8m.
        zrefmn = -4.0
        zrefo  =  1.0  ! impossible value
        do k=1,klist
          zref=max(epsilon*zgrid(k),zrefmn)  ! nearest to zero
          if (zref.ne.zrefo) then  ! new zref
            wref =-2.0*zref
            qwref=1.0/wref
            wq=min(hwide(1),wref)*qwref
            uref=uold(1)*wq
            vref=vold(1)*wq
            bref=-grvity*thref*thold(1)*wq
            wt=0.0
            do ka=2,k
              wt=wt+wq
              if (wt.ge.1.0) then
                exit
              endif
              wq=min(1.0-wt,hwide(ka)*qwref)
              uref=uref+uold(ka)*wq
              vref=vref+vold(ka)*wq
              bref=bref-grvity*thref*thold(ka)*wq
            enddo
          endif
          zrefo=zref

          ritop(k)=(zref-zgrid(k))*(bref+grvity*thref*thold(k))
          dvsq(k)=(uref-uold(k))**2+(vref-vold(k))**2

!         if (vrbos) then
!           if     (k.eq.1) then
!             write(*,'(3a)')
!    &          ' k        z  zref',
!    &          '      u   uref      v   vref',
!    &          '      b   bref    ritop   dvsq'
!           endif
!           write(*,'(i2,f9.2,f6.2,4f7.3,2f7.3,f9.4,f7.4)')
!    &         k,zgrid(k),zref,
!    &         uold(k),uref,vold(k),vref,
!    &         -grvity*thref*thold(k),bref,
!    &         ritop(k),dvsq(k)
!         endif

          if (zgrid(k)*onem*beta_r.gt.-10.0) then
            swfrac(k)=frac_r*exp(zgrid(k)*onem*beta_r)+		&
                      frac_b*exp(zgrid(k)*onem*beta_b)
          elseif (zgrid(k)*onem*beta_b.gt.-10.0) then
            swfrac(k)=frac_b*exp(zgrid(k)*onem*beta_b)
          else
            swfrac(k)=0.0
          endif
          swfrac(k)=swfrac(k)-swfbqp*zgrid(k)*onem  !spread out bottom frac
          if (vrbos) then
            write (*,'(i8,i4,i8,a,f8.2,f8.3)')  nstep,k,i,	&
                '  z,swfrac =',zgrid(k),swfrac(k)
          endif
        enddo  !k=1,klist

! --- calculate interface variables required to estimate interior diffusivities
        do k=1,klist
          kb=k+1
          ka=min(kb,kdm)
          shsq  (kb)=(uold(k)-uold(kb))**2+(vold(k)-vold(kb))**2
          if (.not.locsig) then
            alfadt(kb)=.5*(dsigdt(told(k ),sold(k ))+		&
                           dsigdt(told(kb),sold(kb)))*		&
                         (told(k)-told(kb))
            betads(kb)=.5*(dsigds(told(k ),sold(k ))+		&
                           dsigds(told(kb),sold(kb)))*		&
                         (sold(k)-sold(kb))
            dbloc(kb)=-grvity* thref*(thold(k)-thold(ka))
          else
            dbloc(kb)=-grvity*thref*(alfadt(kb)+betads(kb))
          endif
        enddo

! --- zero 1-d arrays for viscosity/diffusivity calculations

        do k=1,kdm+1
          vcty (k)  =0.0
          dift (k)  =0.0
          difs (k)  =0.0
          ghats(k)  =0.0
          blmc (k,1)=0.0
          blmc (k,2)=0.0
          blmc (k,3)=0.0
          bblmc(k,1)=0.0
          bblmc(k,2)=0.0
          bblmc(k,3)=0.0
        enddo

! --- determine interior diffusivity profiles throughout the water column

! --- shear instability plus background internal wave contributions
        do k=2,klist
          if (shinst) then
            q   =zgrid(k-1)-zgrid(k) !0.5*(hwide(k-1)+hwide(k))
            rigr=max(0.0,dbloc(k)*q/(shsq(k)+epsil))
            ratio=min(rigr*qrinfy,1.0)
            fri=(1.0-ratio*ratio)
            fri=fri*fri*fri
            vcty(k)=min(difm0*fri+dflmiw,difmax)
            difs(k)=min(difs0*fri+dflsiw,difmax)
          else
            vcty(k)=dflmiw
            difs(k)=dflsiw
          endif
          dift(k)=difs(k)
        enddo 

! --- double-diffusion (salt fingering and diffusive convection)
        if (dbdiff) then
          do k=2,klist

! --- salt fingering case
            if (-alfadt(k).gt.betads(k) .and. betads(k).gt.0.) then
              rrho= min(-alfadt(k)/betads(k),rrho0)
              diffdd=1.-((rrho-1.)/(rrho0-1.))**2
              diffdd=dsfmax*diffdd*diffdd*diffdd
              dift(k)=dift(k)+0.7*diffdd
              difs(k)=difs(k)+diffdd

! --- diffusive convection case
            else if ( alfadt(k).gt.0.0 .and. betads(k).lt.0.0		&
               .and. -alfadt(k).gt.betads(k)) then
              rrho=-alfadt(k)/betads(k)
              diffdd=1.5e-6*9.*.101*exp(4.6*exp(-.54*(1./rrho-1.)))
              if (rrho.gt.0.5) then
                prandtl=(1.85-.85/rrho)*rrho
              else
                prandtl=.15*rrho
              endif
              dift(k)=dift(k)+diffdd
              difs(k)=difs(k)+prandtl*diffdd
            endif
          enddo
        endif

        if (vrbos) then
           write (*,102) (nstep,iter,k,i,				&
         hwide(k),1.e4*vcty(k),1.e4*dift(k),1.e4*difs(k),		&
           k=1,kdm  ) 
        endif

        if (iocnmx.gt.4) then
          hbl=dpmixl/onem	! use Kraus-Turner ML depth

        else

! --- calculate boundary layer diffusivity profiles and match these to the
! --- previously-calculated interior diffusivity profiles

! --- diffusivities within the surface boundary layer are parameterized
! --- as a function of boundary layer thickness times a depth-dependent
! --- turbulent velocity scale (proportional to ustar) times a third-order
! --- polynomial shape function of depth. boundary layer diffusivities depend
! --- on surface forcing (the magnitude of this forcing and whether it is
! --- stabilizing or de-stabilizing) and the magnitude and gradient of interior
! --- mixing at the boundary layer base. boundary layer diffusivity profiles
! --- are smoothly matched to interior diffusivity profiles at the boundary
! --- layer base (the profiles and their first derivatives are continuous
! --- at z=-hbl). the turbulent boundary layer depth is diagnosed first, the
! --- boundary layer diffusivity profiles are calculated, then the boundary
! --- and interior diffusivity profiles are combined.

! --- minimum hbl is  top mid-layer + 1 cm or bldmin,
! --- maximum hbl is bottom mid-layer - 1 cm or bldmax.

          hblmin=max(hwide(1)+0.01,bldmin)
          hblmax=min(-zgrid(klist)-0.01,bldmax)

! --- buoyfl = total buoyancy flux (m**2/sec**3) into atmos.
! --- note: surface density increases (column is destabilized) if buoyfl > 0
! --- buoysw = shortwave radiation buoyancy flux (m**2/sec**3) into atmos.
! --- salflx, sswflx and surflx are positive into the ocean
          tmn=temp(1)
          smn=saln(1)
          dsgdt=          dsigdt(tmn,smn)
          buoyfs=grvity*thref*(dsigds(tmn,smn)*salflx*thref)
          buoyfl=buoyfs+						&
                 grvity*thref*(dsgdt          *surflx*thref/spcifh)
          buoysw=grvity*thref*(dsgdt          *sswflx*thref/spcifh)
! 
! --- diagnose the new boundary layer depth as the depth where a bulk
! --- richardson number exceeds ric

! --- initialize hbl and nbl to bottomed out values
          kup2=1
          kup =2
          kdn =3
          rib(kup2)=0.0
          rib(kup) =0.0
          nbl=klist
          hbl=hblmax

! --- find nbl (=nblmin) associated with minimum mixed layer depth
          do k=2,nbl
            nblmin=k
            if (zgrid(k).lt.-hblmin) exit
          end do

! --- diagnose hbl and nbl
          do k=2,nbl
            case=-zgrid(k)
            bfsfc=buoyfl-swfrac(k)*buoysw
            if     (bfsfc.le.0.0) then
              stable=1.0
              dnorm =1.0
            else
              stable=0.0
              dnorm =epsilon
            endif

! --- compute turbulent velocity scales at dnorm, for
! --- hbl = case = -zgrid(k)
            call wscale(i,case,dnorm,bfsfc,wm,ws,1,ustar,ustarb)

! --- compute the turbulent shear contribution to rib
            if     (max(dbloc(k),dbloc(k+1)).gt.0.0) then
              bfq=0.5*(dbloc(k  )/(zgrid(k-1)-zgrid(k  ))+		&
                       dbloc(k+1)/(zgrid(k  )-zgrid(k+1)) )
              if     (bfq.gt.0.0) then
                bfq=sqrt(bfq)
              else
                bfq=0.0  !neutral or unstable
              endif
            else
              bfq=0.0  !neutral or unstable
            endif
            if     (bfq.gt.0.0) then
              if     (cv.ne.0.0) then
                cvk=cv
              else !frequency dependent version
                cvk=max(cv_max-cv_bfq*bfq,cv_min) !between cv_min and cv_max
              endif
              vtsq=-zgrid(k)*ws*bfq*vtc*cvk
            else
              vtsq=0.0
            endif !bfq>0:else

! --- compute bulk richardson number at new level
            rib(kdn)=ritop(k)/(dvsq(k)+vtsq+epsil)
            if (nbl.eq.klist.and.rib(kdn).ge.ricr) then
! ---     interpolate to find hbl as the depth where rib = ricr
              if     (k.eq.2 .or. hblflg.eq.0) then  !nearest interface
                hbl = -zgrid(k-1)+0.5*hwide(k-1)
              elseif (k.lt.4 .or. hblflg.eq.1) then  !linear
                hbl = -zgrid(k-1)+ (zgrid(k-1)-zgrid(k))*		&
                         (ricr-rib(kup))/(rib(kdn)-rib(kup)+epsil)
              else !quadratic

! ---           Determine the coefficients A,B,C of the polynomial
! ---             Y(X) = A * (X-X2)**2 + B * (X-X2) + C
! ---           which goes through the data: (X[012],Y[012])

                x0 = zgrid(k-2)
                x1 = zgrid(k-1)
                x2 = zgrid(k)
                y0 = rib(kup2)
                y1 = rib(kup)
                y2 = rib(kdn)
                ahbl = ( (y0-y2)*(x1-x2) -		&
                         (y1-y2)*(x0-x2)  )/		&
                       ( (x0-x2)*(x1-x2)*(x0-x1) )
                bhbl = ( (y1-y2)*(x0-x2)**2 -		&
                         (y0-y2)*(x1-x2)**2  ) /	&
                       ( (x0-x2)*(x1-x2)*(x0-x1) )
                if     (abs(bhbl).gt.epsil) then
                  lhbl = abs(ahbl)/abs(bhbl).gt.epsil
                else
                  lhbl = .true.
                endif
                if     (lhbl) then !quadratic
! ---             find root of Y(X)-RICR nearest to X2
                  chbl = y2 - ricr
                  dhbl = bhbl**2 - 4.0*ahbl*chbl
                  if     (dhbl.lt.0.0) then !linear
                    hbl = -(x2 + (x1-x2)*(y2-ricr)/(y2-y1+epsil))
                  else
                    dhbl = sqrt(dhbl)
                    if     (abs(bhbl+dhbl).ge.abs(bhbl-dhbl)    ) then
                      hbl = -(x2 - 2.0*chbl/(bhbl+dhbl))
                    else
                      hbl = -(x2 - 2.0*chbl/(bhbl-dhbl))
                    endif !nearest root
                  endif !bhbl**2-4.0*ahbl*chbl.lt.0.0:else
                else !linear
                  hbl = -(x2 + (x1-x2)*(y2-ricr)/(y2-y1+epsil))
                endif !quadratic:linear
              endif !linear:quadratic
              nbl=k
              if (hbl.lt.hblmin) then
                hbl=hblmin
                nbl=nblmin
              endif
              if (hbl.gt.hblmax) then
                hbl=hblmax
                nbl=klist
              endif
              exit !k-loop
            endif

            ksave=kup2
            kup2=kup
            kup =kdn
            kdn =ksave
          enddo  !k=1,nbl

! --- calculate swfrml, the fraction of solar radiation left at depth hbl
          if     (-hbl*onem*beta_r.gt.-10.0) then
            swfrml=frac_r*exp(-hbl*onem*beta_r)+		&
                   frac_b*exp(-hbl*onem*beta_b)
          elseif (-hbl*onem*beta_b.gt.-10.0) then
            swfrml=frac_b*exp(-hbl*onem*beta_b)
          else
            swfrml=0.0
          endif
          swfrml=swfrml-swfbqp*hbl*onem  !spread out bottom frac
          if (vrbos)						&
            write (*,'(2i8,i5,a,4es10.2)') nstep,i,nbl,		&
                '  hbl,swfrml =',hbl,swfrml

! --- limit check on hbl for negative (stablizing) surface buoyancy forcing
          bfsfc=buoyfl-swfrml*buoysw
          if (bfsfc.le.0.0) then
            bfsfc=bfsfc-epsil  !insures bfsfc never=0
            hmonob0=min(-cmonob*ustar**3/(vonk*bfsfc), hblmax)
            hbl=max(hblmin,min(hbl,hekman,hmonob0))
          else
            hmonob0=hblmax
          endif
          dpmixl=hbl*onem
        end if		! KT or Ri-based ML depth

! --- find new nbl and re-calculate swfrml
        nbl=klist
        do k=2,klist
          if (-zgrid(k).gt.hbl) then
            nbl=k
            exit
          endif
        enddo
        if     (-hbl*onem*beta_r.gt.-10.0) then
          swfrml=frac_r*exp(-hbl*onem*beta_r)+			&
                 frac_b*exp(-hbl*onem*beta_b)
        elseif (-hbl*onem*beta_b.gt.-10.0) then
          swfrml=frac_b*exp(-hbl*onem*beta_b)
        else
          swfrml=0.0
        endif
        swfrml=swfrml-swfbqp*hbl*onem  !spread out bottom frac
        if (vrbos) then
          write (*,'(2i8,i5,a,4es10.2)') nstep,i,nbl,		&
              '  hbl,swfrml =',hbl,swfrml
        end if

! --- find forcing stability and buoyancy forcing for final hbl values
! --- determine case (for case=0., hbl lies between -zgrid(i,nbl)
! --- and the interface above. for case=1., hbl lies between 
! --- -zgrid(i,nbl-1) and the interface below)

! --- velocity scales at hbl
        bfsfc=buoyfl-swfrml*buoysw
        if     (bfsfc.le.0.0) then
          bfsfc=bfsfc-epsil  !insures bfsfc never=0
          stable=1.0
          dnorm =1.0
        else
          stable=0.0
          dnorm =epsilon
        endif
!JR sign function requires matching types so change .5 to .5_8
        case=.5+sign(.5_8,-zgrid(nbl)-.5*hwide(nbl)-hbl)

        buoflx0=bfsfc                          !mixed layer buoyancy
        bhtflx0=bfsfc-buoyfs                   !buoyancy from heat flux
        mixflx0=surflx-swfrml*sswflx !mixed layer heat flux

        call wscale(i,hbl,dnorm,bfsfc,wm,ws,1,ustar,ustarb)

! --- compute the boundary layer diffusivity profiles. first, find interior
! --- viscosities and their vertical derivatives at hbl
        ka=nint(case)*(nbl-1)+(1-nint(case))*nbl
        q=(hbl*onem-p(ka))*qdpmm(ka)
        vctyh=vcty(ka)+q*(vcty(ka+1)-vcty(ka))
        difsh=difs(ka)+q*(difs(ka+1)-difs(ka))
        difth=dift(ka)+q*(dift(ka+1)-dift(ka))

        q=(hbl+zgrid(nbl-1))/(zgrid(nbl-1)-zgrid(nbl))
        dvdzup=(vcty(nbl-1)-vcty(nbl  ))/hwide(nbl-1)
        dvdzdn=(vcty(nbl  )-vcty(nbl+1))/hwide(nbl  )
        viscp=.5*((1.-q)*(dvdzup+abs(dvdzup))+q*(dvdzdn+abs(dvdzdn)))
        dvdzup=(difs(nbl-1)-difs(nbl  ))/hwide(nbl-1)
        dvdzdn=(difs(nbl  )-difs(nbl+1))/hwide(nbl  )
        difsp=.5*((1.-q)*(dvdzup+abs(dvdzup))+q*(dvdzdn+abs(dvdzdn)))
        dvdzup=(dift(nbl-1)-dift(nbl  ))/hwide(nbl-1) 
        dvdzdn=(dift(nbl  )-dift(nbl+1))/hwide(nbl  )
        diftp=.5*((1.-q)*(dvdzup+abs(dvdzup))+q*(dvdzdn+abs(dvdzdn)))

        f1=-stable*c11*bfsfc/(ustar**4+epsil) 

        gat1(1)=vctyh/hbl/(wm+epsil)
        dat1(1)=min(0.,-viscp/(wm+epsil)+f1*vctyh)

        gat1(2)=difsh/hbl/(ws+epsil)
        dat1(2)=min(0.,-difsp/(ws+epsil)+f1*difsh) 

        gat1(3)=difth/hbl/(ws+epsil)
        dat1(3)=min(0.,-diftp/(ws+epsil)+f1*difth)

! --- compute turbulent velocity scales on the interfaces
        do k=2,kdm+1
          if (k.le.min(nbl,klist)) then
            sigg=p(k)/(hbl*onem)
            dnorm=stable*sigg+(1.-stable)*min(sigg,epsilon)

            call wscale(i,hbl,dnorm,bfsfc,wm,ws,1,ustar,ustarb)

! --- compute the dimensionless shape functions at the interfaces
            aa1=sigg-2.
            aa2=3.-2.*sigg
            aa3=sigg-1.

            gm=aa1+aa2*gat1(1)+aa3*dat1(1) 
            gs=aa1+aa2*gat1(2)+aa3*dat1(2)
            gt=aa1+aa2*gat1(3)+aa3*dat1(3)

! --- compute boundary layer diffusivities at the interfaces
            blmc(k,1)=hbl*wm*sigg*(1.+sigg*gm)
            blmc(k,2)=hbl*ws*sigg*(1.+sigg*gs)
            blmc(k,3)=hbl*ws*sigg*(1.+sigg*gt)

! --- compute nonlocal transport forcing term = ghats * <ws>o
            if (nonloc) then
              ghats(k)=(1.-stable)*cg/(ws*hbl+epsil)
            endif
          endif !k.le.min(nbl,klist)
        enddo !k

! --- enhance diffusivities on the interface closest to hbl

! --- first compute diffusivities at nbl-1 grid level 
        sigg=-zgrid(nbl-1)/hbl
        dnorm=stable*sigg+(1.-stable)*min(sigg,epsilon)

        call wscale(i,hbl,dnorm,bfsfc,wm,ws,1,ustar,ustarb)

        sigg=-zgrid(nbl-1)/hbl
        aa1=sigg-2.
        aa2=3.-2.*sigg
        aa3=sigg-1.
        gm=aa1+aa2*gat1(1)+aa3*dat1(1)
        gs=aa1+aa2*gat1(2)+aa3*dat1(2)
        gt=aa1+aa2*gat1(3)+aa3*dat1(3)
        dkm1(1)=hbl*wm*sigg*(1.+sigg*gm)
        dkm1(2)=hbl*ws*sigg*(1.+sigg*gs)
        dkm1(3)=hbl*ws*sigg*(1 +sigg*gt)

! --- now enhance diffusivity at interface nbl

! --- this procedure was altered for hycom to reduce diffusivity enhancement
! --- if the interface in question is located more than dp0enh below hbl.
! --- this prevents enhanced boundary layer mixing from penetrating too far
! --- below hbl when hbl is located in a very thick layer
        k=nbl-1
        ka=k+1
        delta=(hbl+zgrid(k))/(zgrid(k)-zgrid(ka))

        dkmp2=case*vcty(ka)+(1.-case)*blmc(ka,1)
        dstar=(1.-delta)**2*dkm1(1)+delta**2*dkmp2      
        blmc(ka,1)=(1.-delta)*vcty(ka)+delta*dstar

        dkmp2=case*difs(ka)+(1.-case)*blmc(ka,2)
        dstar=(1.-delta)**2*dkm1(2)+delta**2*dkmp2    
        blmc(ka,2)=(1.-delta)*difs(ka)+delta*dstar

        dkmp2=case*dift(ka)+(1.-case)*blmc(ka,3)
        dstar=(1.-delta)**2*dkm1(3)+delta**2*dkmp2     
        blmc(ka,3)=(1.-delta)*dift(ka)+delta*dstar

        if (case.eq.1.) then
          q=1.-case*max(0.,min(1.,(p(ka)-hbl*onem-dp0enh)/dp0enh))
          blmc(ka,1)=max(vcty(ka),q*blmc(ka,1))
          blmc(ka,2)=max(difs(ka),q*blmc(ka,2))
          blmc(ka,3)=max(dift(ka),q*blmc(ka,3))
        endif

        if (nonloc) then
          ghats(ka)=(1.-case)*ghats(ka)
        endif

! --- combine interior and boundary layer coefficients and nonlocal term
        if (.not.bblkpp) then
          do k=2,nbl
            vcty(k)=max(vcty(k),min(blmc(k,1),difmax))
            difs(k)=max(difs(k),min(blmc(k,2),difmax))
            dift(k)=max(dift(k),min(blmc(k,3),difmax))
          enddo
          do k=nbl+1,klist
            ghats(k)=0.0
          enddo
          do k=klist+1,kdm+1
            vcty(k)=dflmiw
            difs(k)=dflsiw
            dift(k)=dflsiw
            ghats(k)=0.0
          enddo
        endif !.not.bblkpp

        if (vrbos) then
          write (*,103) (nstep,iter,k,i,			&
          hwide(k),1.e4*vcty(k),1.e4*dift(k),1.e4*difs(k),	&
          ghats(k),k=1,kdm  )			! TNL
        endif

! --- save array dpbl=onem*hbl for ice, output and diagnosis
        dpbl0=onem*hbl

        if (bblkpp) then

! ------------------------------------------------
!
! --- begin bottom boundary layer parameterization
!
! ------------------------------------------------

! --- this bottom boundary algorithm follows the kpp algorithm included
! --- in the rutgers roms model. it is essentially an adaptation of the
! --- algorithm used for the surface boundary layer, involving diagnosis
! --- of the bottom boundary layer thickness hbbl using a bulk
! --- richardson number

! --- calculate zgridb
        do k=klist+1,1,-1
          zgridb(k)=zgrid(k)-zgrid(klist+1)
        enddo

! --- calculate bottom boundary layer diffusivity profiles and match these
! --- to the existing profiles

! --- minimum hbbl is 1 m, maximum is distance between bottom and one meter
! --- below the base of model layer 1

        hbblmin=1.0
        hbblmax=zgridb(2)
!     hbblmax=min(-hbl-zgrid(klist+1),2.0*thkbot)
!     
! --- buoyfl = buoyancy flux (m**2/sec**3) into bottom due to heating by
! ---          the penetrating shortwave radiation
! --- note: bottom density increases (column is destabilized) if buoyfl < 0
! --- buoysw = shortwave radiation buoyancy flux (m**2/sec**3) at the surface

! --- NOTE: the convention for the bottom bl in the roms model was to use the
! --- surface net (turbulent plus radiative) heat flux to represent buoyfl
! --- this convention is not used here - instead, bottom turbulent heat
! --- flux arises entirely due to heating of the bottom by the penetrating
! --- shortwave radiation. as a result, net heat flux (buoyfl) at the bottom
! --- is zero since upward turbulent heat flux due to bottom heating is
! --- opposed by the downward shortwave radiative heat flux (it is presently
! --- assumed that no heat is absorbed by the bottom). Moving upward from
! --- the bottom, penetrating shortwave radiation acts to stabilize the
! --- water column.

        if (locsig) then
          dsgdt=dsiglocdt(temp(k),saln(k),0.5*(p(klist)+p(klist+1)))
        else
          dsgdt=dsigdt(temp(k),saln(k))
        endif
        buoysw=-grvity*thref*dsgdt*sswflx*thref/spcifh
        buoyfl=-swfrac(klist+1)*buoysw
!     
! --- diagnose the new boundary layer depth as the depth where a bulk
! --- richardson number exceeds ric

! --- initialize hbbl and nbbl to extreme values
        kdn2=1
        kdn =2
        kup =3
        rib(kdn2)=0.0
        rib(kdn) =0.0
        nbbl=2
        hbbl=hbblmax

! --- nearbottom reference values of model variables are handled
! --- differently from the surface layer because the surface
! --- procedure does not work properly with the highly-uneven
! --- layer thicknesses often present near the bottom.
! --- reference values are chosen as the values present at the
! --- bottom = hence, uref, vref are zero and not used while
! --- bottom buoyancy is estimated assuming a linear vertical
! --- profile across the bottom layer

        bref=grvity*thref*(0.5*(3.0*thold(klist)-thold(klist-1))+35.)  ! 35 stands for thbase
!     
! --- diagnose hbbl and nbbl
        do k=klist,nbbl,-1
          ritop(k)=max(zgridb(k)*(bref-grvity*thref*thold(k)),epsil)
          dvsq(k)=uold(k)**2+vold(k)**2

          case=zgridb(k)
          bfbot=0.0
          stable=1.0
          dnorm =1.0

! --- compute turbulent velocity scales at dnorm, for
! --- hbbl = case = zgridb(k)
          call wscale(i,case,dnorm,bfbot,wm,ws,2,ustar,ustarb)

! --- compute the turbulent shear contribution to rib
          if     (max(dbloc(k),dbloc(k+1)).gt.0.0) then
            bfq=0.5*(dbloc(k  )/(zgrid(k-1)-zgrid(k ))+		&
                     dbloc(k+1)/(zgrid(k  )-zgrid(k+1)) )
            if     (bfq.gt.0.0) then
              bfq=sqrt(bfq)
            else
              bfq=0.0  !neutral or unstable
            endif
          else
            bfq=0.0  !neutral or unstable
          endif
          if     (bfq.gt.0.0) then
            if     (cv.ne.0.0) then
              cvk=cv
            else !frequency dependent version
              cvk=max(cv_max-cv_bfq*bfq,cv_min) !between cv_min and cv_max
            endif
            vtsq=zgridb(k)*ws*bfq*vtc*cvk
          else
            vtsq=0.0
          endif !bfq>0:else

! --- compute bulk richardson number at new level
! --- interpolate to find hbbl as the depth where rib = ricrb
! --- in stable or neutral conditions, hbbl can be no thicker than the
! --- bottom ekman layer
! --- ustarb is estimated in momtum.f

          rib(kup)=ritop(k)/(dvsq(k)+vtsq+epsil)
          if (rib(kup).ge.ricrb) then
            hekmanb=ustarb*cekman/max( cormn4,abs(corio))
            if     (hblflg.eq.0) then                         !nearest intf.
              hbbl = zgridb(k+1)-0.5*hwide(k+1)
            elseif (k.gt.klist-2 .or. hblflg.eq.1) then  !linear
              hbbl = zgridb(k+1)-				&
                       (zgridb(k+1)-zgridb(k))*			&
                       (ricrb-rib(kdn))/(rib(kup)-rib(kdn)+epsil)
            else                                              !quadratic

! ---         Determine the coefficients A,B,C of the polynomial
! ---           Y(X) = A * (X-X2)**2 + B * (X-X2) + C
! ---         which goes through the data: (X[012],Y[012])

              x0 = -zgridb(k+2)
              x1 = -zgridb(k+1)
              x2 = -zgridb(k)
              y0 = rib(kdn2)
              y1 = rib(kdn)
              y2 = rib(kup)
              ahbl = ( (y0-y2)*(x1-x2) -			&
                       (y1-y2)*(x0-x2)  )/			&
                     ( (x0-x2)*(x1-x2)*(x0-x1) )
              bhbl = ( (y1-y2)*(x0-x2)**2 -			&
                       (y0-y2)*(x1-x2)**2  ) /			&
                     ( (x0-x2)*(x1-x2)*(x0-x1) )
              if     (abs(bhbl).gt.epsil) then
                lhbl = abs(ahbl)/abs(bhbl).gt.epsil
              else
                lhbl = .true.
              endif
              if     (lhbl) then !quadratic
! ---           find root of Y(X)-RICR nearest to X2
                chbl = y2 - ricrb
                dhbl = bhbl**2 - 4.0*ahbl*chbl
                if     (dhbl.lt.0.0) then !linear
                  hbbl = -(x2 + (x1-x2)*(y2-ricrb)/(y2-y1+epsil))
                else
                  dhbl = sqrt(dhbl)
                  if     (abs(bhbl+dhbl) .ge. abs(bhbl-dhbl)) then
                    hbbl = -(x2 - 2.0*chbl/(bhbl+dhbl))
                  else
                    hbbl = -(x2 - 2.0*chbl/(bhbl-dhbl))
                  endif !nearest root
                endif !bhbl**2-4.0*ahbl*chbl.lt.0.0:else
              else !linear
                hbbl = -(x2 + (x1-x2)*(y2-ricrb)/(y2-y1+epsil))
              endif !quadratic:linear
            endif
            hbbl=max(hbblmin,min(hekmanb,hbblmax,hbbl))
            exit !k-loop
          endif

          ksave=kdn2
          kdn2=kdn
          kdn =kup
          kup =ksave
        enddo  !k=klist,nbbl,-1

!--- find new nbbl
        nbbl=2
        do k=klist,2,-1
          if (zgridb(k).gt.hbbl) then
            nbbl=k
            exit
          endif
        enddo  !k=klist,2,-1

! --- do not execute the remaining bottom boundary layer algorithm if
! --- vertical resolution is not available in the boundary layer

!     if (hbbl.lt.0.5*hwide(klist)) then
        if (hbbl.lt.0.5*    hwide(klist)+		&
                    0.5*min(hwide(klist  ),		&
                            hwide(klist-1) )) then
          go to 201  !one grid point in bbl and probably not a "plume"
        endif

! --- calculate swfrml, the fraction of solar radiation absorbed by depth hbbl
        q=(zgridb(nbbl-1)-hbbl)/(zgridb(nbbl-1)-zgridb(nbbl))
        swfrml=swfrac(nbbl-1)+q*(swfrac(nbbl)-swfrac(nbbl-1))

! --- find forcing stability and buoyancy forcing for final hbbl values
! --- determine case (for case=0., hbbl lies between -zgridb(nbbl)
! --- and the interface below. for case=1., hbbl lies between 
! --- -zgrid(nbbl+1) and the interface above)

! --- velocity scales at hbbl
        bfbot=buoyfl+swfrml*buoysw
        if     (bfbot.ge.0.0) then
          bfbot=bfbot+epsil  !insures bfbot never=0
          stable=1.0
          dnorm =1.0
        else
          stable=0.0
          dnorm =epsilon
        endif
        case=.5+sign(.5,zgridb(nbbl)-.5*hwide(nbbl)-hbbl)

        call wscale(i,hbbl,dnorm,bfbot,wm,ws,2,ustar,ustarb)

! --- compute the boundary layer diffusivity profiles. first, find interior
! --- viscosities and their vertical derivatives at hbbl
        ka=nint(case)*(nbbl+1)+(1-nint(case))*nbbl
        q=(p(klist+1)-p(ka)-hbbl*onem)*qdpmm(ka)
        vctyh=vcty(ka)+q*(vcty(ka+1)-vcty(ka))
        difsh=difs(ka)+q*(difs(ka+1)-difs(ka))
        difth=dift(ka)+q*(dift(ka+1)-dift(ka))

        q=(hbbl-zgridb(nbbl+1))/(zgridb(nbbl)-zgridb(nbbl+1))
        dvdzup=-(vcty(nbbl  )-vcty(nbbl+1))/hwide(nbbl  )
        dvdzdn=-(vcty(nbbl+1)-vcty(nbbl+2))/hwide(nbbl+1)
        viscp=.5*((1.-q)*(dvdzup+abs(dvdzup))+q*(dvdzdn+abs(dvdzdn)))
        dvdzup=-(difs(nbbl  )-difs(nbbl+1))/hwide(nbbl  )
        dvdzdn=-(difs(nbbl+1)-difs(nbbl+2))/hwide(nbbl+1)
        difsp=.5*((1.-q)*(dvdzup+abs(dvdzup))+q*(dvdzdn+abs(dvdzdn)))
        dvdzup=-(dift(nbbl  )-dift(nbbl+1))/hwide(nbbl) 
        dvdzdn=-(dift(nbbl+1)-dift(nbbl+2))/hwide(nbbl+1)
        diftp=.5*((1.-q)*(dvdzup+abs(dvdzup))+q*(dvdzdn+abs(dvdzdn)))

        f1=stable*c11*bfbot/(ustarb**4+epsil) 

        gat1(1)=vctyh/hbbl/(wm+epsil)
        dat1(1)=min(0.,-viscp/(wm+epsil)+f1*vctyh)

        gat1(2)=difsh/hbbl/(ws+epsil)
        dat1(2)=min(0.,-difsp/(ws+epsil)+f1*difsh) 

        gat1(3)=difth/hbbl/(ws+epsil)
        dat1(3)=min(0.,-diftp/(ws+epsil)+f1*difth)

! --- compute turbulent velocity scales on the interfaces
        do k=klist,nbbl+1,-1
          sigg=(p(klist+1)-p(k))/(hbbl*onem)
          dnorm=stable*sigg+(1.-stable)*min(sigg,epsilon)

          call wscale(i,hbbl,dnorm,bfbot,wm,ws,2,ustar,ustarb)

! --- compute the dimensionless shape functions at the interfaces
          aa1=sigg-2.
          aa2=3.-2.*sigg
          aa3=sigg-1.

          gm=aa1+aa2*gat1(1)+aa3*dat1(1) 
          gs=aa1+aa2*gat1(2)+aa3*dat1(2)
          gt=aa1+aa2*gat1(3)+aa3*dat1(3)

! --- compute boundary layer diffusivities at the interfaces
          bblmc(k,1)=hbbl*wm*sigg*(1.+sigg*gm)
          bblmc(k,2)=hbbl*ws*sigg*(1.+sigg*gs)
          bblmc(k,3)=hbbl*ws*sigg*(1.+sigg*gt)
        enddo !k=klist,nbbl+1,-1

! --- if the model interface nbbl+1 is located more than the distance
! --- dp0bbl above hbbl, reduce diffusivity to prevent bottom mixing
! --- from penetrating too far into the interior

        k=nbbl+1
        delta=(p(klist+1)-p(nbbl+1))/onem-hbbl
        if (delta.gt.dp0bbl) then
          dstar=max(0.0,1.0+(dp0bbl-delta)/dp0bbl)
          bblmc(k,1)=bblmc(k,1)*dstar
          bblmc(k,2)=bblmc(k,2)*dstar
          bblmc(k,3)=bblmc(k,3)*dstar
        endif

 201  continue  ! skip bbl algorithm due to poor vertical resolution

! --- save array dpbbl=onem*hbbl for output and diagnosis, and for momtum.f
        dpbbl0=onem*hbbl

! --- select maximum viscosity/diffusivity at all interfaces
        do k=2,klist
          if (k.le.klist) then
            vcty(k)=min(difmax,max(vcty(k),blmc(k,1),bblmc(k,1)))
            difs(k)=min(difmax,max(difs(k),blmc(k,2),bblmc(k,2)))
            dift(k)=min(difmax,max(dift(k),blmc(k,3),bblmc(k,3)))
          else
            vcty(k)=dflmiw
            difs(k)=dflsiw
            dift(k)=dflsiw
          endif
          if (k.ge.nbl+1) then
            ghats(k)=0.0
          endif
        enddo

        if (vrbos) then
          write (*,103) (nstep,iter,k,i,			&
        hwide(k),1.e4*vcty(k),1.e4*dift(k),1.e4*difs(k),	&
        ghats(k),k=kdm,1,-1)
        endif
        if (vrbos .and. mod(nstep,20).eq.0) then
          print *,'nbbl,hbbl',nbbl,hbbl
        endif

        endif    ! bblkpp

! <><><><><><><><>  p a r t i a l   c o l u m n   k p p  <><><><><><><><>
! --- to turn off mixing below mixed layer, reduce -klist- to mxlyr depth
        if (mod(iocnmx,4).eq.2) then
          q=0.
          do k=1,klist
            q=q+dpmm(k)
            if (q.ge.dpmixl) exit
          end do
!c        if (k.gt.klist) then
!c          print *,'warning: klist redefinition problem, k,i,klist ='
!c   .       ,k,i,klist
!c          print *,'q,dpmixl:',q,dpmixl
!c        end if
          klist=min(k+1,klist)
          if (vrbos) then
            write (*,'(2i8,a,i3)') nstep,i,			&
                     '  no diffusion beyond layer',klist
          end if
        end if                  ! iocnmx = 2
! <><><><><><><><>  p a r t i a l   c o l u m n   k p p  <><><><><><><><>

        if (iocnmx.le.2) then

! --- perform the vertical mixing at p points

          do k=1,klist
            difft(k+1)=dift(k+1)
            diffs(k+1)=difs(k+1)
            diffm(k+1)=vcty(k+1)
            ghat(k+1)=ghats(k+1)
            t1do(k)=temp(k)
            s1do(k)=saln(k)
            u1do(k)=   u(k)
            v1do(k)=   v(k)
            hm(k) =hwide(k)
            zm(k) =zgrid(k)
          enddo

          nlayer=klist
          k=nlayer+1
          ka=min(k,kdm)
          difft(k)=0.0
          diffs(k)=0.0
          diffm(k)=0.0
          ghat(k)=0.0
          t1do(k)=temp(ka)
          s1do(k)=saln(ka)
          u1do(k)=u1do(k-1)
          v1do(k)=v1do(k-1)
          zm(k)=zgrid(k)

! --- compute factors for coefficients of tridiagonal matrix elements.
!         tri(k=1:NZ,0) : dt/hwide(k)/ dzb(k-1)=z(k-1)-z(k)=dzabove)
!         tri(k=1:NZ,1) : dt/hwide(k)/(dzb(k  )=z(k)-z(k+1)=dzbelow)

          do k=1,nlayer
            dzb(k)=zm(k)-zm(k+1)
          enddo

          tri(1,1)=dlt/(hm(1)*dzb(1))
          tri(1,0)=0.
          do k=2,nlayer
            tri(k,1)=dlt/(hm(k)*dzb(k))
            tri(k,0)=dlt/(hm(k)*dzb(k-1))
          enddo

! --- solve the diffusion equation
! --- salflx, sswflx and surflx are positive into the ocean

! --- t solution
          ghatflux=-(surflx-sswflx)*thref/spcifh
          call tridcof(difft,tri,nlayer,tcu,tcc,tcl)
          call tridrhs(hm,t1do,difft,ghat,ghatflux,nlayer,rhs,dlt)
          call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,t1do,t1dn,difft)

! --- s solution
          ghatflux=-salflx*thref
          call tridcof(diffs,tri,nlayer,tcu,tcc,tcl)
          call tridrhs(hm,s1do,diffs,ghat,ghatflux,nlayer,rhs,dlt)
          call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,s1do,s1dn,diffs)

          if (vrbos) then
            write (*,104) (nstep,iter,k,i,hm(k),t1do(k),t1dn(k),	&
              s1do(k),s1dn(k),0.0,0.0,k=1,nlayer)
          endif

! --- u solution
          call tridcof(diffm,tri,nlayer,tcu,tcc,tcl)
          do k=1,nlayer
            rhs(k)=u1do(k)
          enddo
          call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,u1do,u1dn,diffm)

! --- v solution
          do k=1,nlayer
            rhs(k)=v1do(k)
          enddo
          call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,v1do,v1dn,diffm)

          if (vrbos) then
            write (*,105) (nstep,iter,k,i,hm(k),u1do(k),u1dn(k),	&
                            v1do(k),v1dn(k),k=1,nlayer)
          endif

! --- reset old variables in preparation for next iteration
          do k=1,nlayer+1
            told(k)=t1dn(k)
            sold(k)=s1dn(k)
            if (locsig) then
              if (k.eq.1) then
                thold(k)=sigocn(told(k),sold(k))
              else
                ka=k-1
                alfadt(k)=0.5*					&
                         (dsiglocdt(told(ka),sold(ka),p(k))+	&
                          dsiglocdt(told(k ),sold(k ),p(k)))*	&
                         (told(ka)-told(k))
                betads(k)=0.5*					&
                         (dsiglocds(told(ka),sold(ka),p(k))+	&
                          dsiglocds(told(k ),sold(k ),p(k)))*	&
                         (sold(ka)-sold(k))
                thold(k)=thold(ka)-alfadt(k)-betads(k)
              endif
            else
              thold(k)=sigocn(told(k),sold(k))
            endif
!rb         if (iter.lt.niter) then
              uold(k)=u1dn(k)
              vold(k)=v1dn(k)
!rb         endif
          enddo
        end if		! iocnmx = 1 or 2
      end do		! iteration loop

 101  format(i8,i3,i8,'  swfrac,dn,dtemp,dsaln ',2f8.3,2f11.6)
 102  format(26x,'   thick    viscty    t_diff    s_diff  '		&
           /(i8,2i3,i8,2x,4f10.2))
 103  format(26x,'   thick    viscty    t_diff    s_diff   nonlocal'	&
           /(i8,2i3,i8,2x,4f10.2,f11.6))
 104  format(26x,							&
           ' thick   t_old   t_new   s_old   s new trc_old trc_new'	&
           /(i8,2i3,i8,f10.2,4f8.3,2f8.4))
 105  format(27x,'  thick   u_old   u_new   v_old   v_new'		&
           /(i8,2i3,i8,2x,f10.2,4f8.3))

      return
      end subroutine mxkppaij

      subroutine mxkppbij (nstep,zgrid,dift,difs,ghats,dp,pbot,		&
                           surflx,sswflx,salflx,dlt,klist,i,		&
                           vrbos,dotrcr,deg_lat,temp,saln,passv_tr)

! --- hycom version 1.0
! -------------------------------------------------------------
! --- k-profile vertical diffusion, single j-row (part B)
! --- vertical coordinate is z negative below the ocean surface
! -------------------------------------------------------------

      use module_constants,only: grvity
      use hycom_constants ,only: thref, spcifh, onem, onemm
      use hycom_control   ,only: numtr, trcfrq, bclin_frq
      use fimnamelist     ,only: kdm, ocnmx_factor_s, ocnmx_factor_t	!factor to reduce difs/dift
      use hycom_sigetc

      implicit none
      integer,intent(IN) :: i,nstep,klist
      logical,intent(IN) :: vrbos,dotrcr
      real   ,intent(IN) :: dlt, deg_lat
      real*8 ,intent(IN) :: zgrid(kdm+1),dift(kdm+1),difs(kdm+1),ghats(kdm+1)
      real*8 ,intent(IN) :: dp(kdm),pbot  
      real, intent(IN) :: surflx, sswflx, salflx
      real*8, intent(INOUT) :: temp(kdm),saln(kdm)
      real, intent(INOUT),optional:: passv_tr(kdm,numtr)

! --- perform the final vertical mixing at p points

! --- local 1-d arrays for matrix solution
      real*8 t1do(kdm+1),t1dn(kdm+1),s1do(kdm+1),s1dn(kdm+1),	&
             tr1do(kdm+1,numtr),tr1dn(kdm+1,numtr),		&
             difft(kdm+1),diffs(kdm+1),difftr(kdm+1),		&
             ghat(kdm+1),zm(kdm+1),hm(kdm),dzb(kdm),		&
             totemo,totemn,tosalo,tosaln,tndcyt,tndcys,		&
             totrco(numtr),totrcn(numtr),trscal(numtr)

! --- tridiagonal matrix solution arrays
      real*8 tri(kdm,0:1)      ! dt/dz/dz factors in trid. matrix
      real*8 tcu(kdm),      &  ! upper coeff for (k-1) on k line of trid.matrix
             tcc(kdm),      &  ! central ...     (k  ) ..
             tcl(kdm),      &  ! lower .....     (k-1) ..
             rhs(kdm),      &  ! right-hand-side terms
             ghatflux,salin,sumthk,frac
      integer k,ka,ktr,nlayer

      real,   parameter :: acurcy=1.e-7
      real,   parameter :: difriv =   60.0e-4  !river diffusion
      real,   parameter :: tofset = 0.
      real,   parameter :: sofset = 0.
      integer,parameter :: tsofrq = 10
      real,   parameter :: fresh = 1.		! smallest allowed srf.salinity
      real              :: sqlat

      difft(:)=0. 
      diffs(:)=0. 
      difftr(:)=0. 
      ghat(:)=0. 
      nlayer=klist

      sqlat=(.1*deg_lat)**2
      do k=1,nlayer
        difft( k+1)=dift(k+1)*(ocnmx_factor_t+sqlat)/(1.+sqlat)
        diffs( k+1)=difs(k+1)*(ocnmx_factor_s+sqlat)/(1.+sqlat)
        difftr(k+1)=difs(k+1)
        ghat(k+1)= ghats(k+1)
        t1do(k)=temp(k)
        s1do(k)=saln(k)
        t1dn(k)=t1do(k)
        s1dn(k)=s1do(k)
        if (dotrcr .and. present(passv_tr)) then
          do ktr= 1,numtr
            tr1do(k,ktr)=passv_tr(k,ktr)
            tr1dn(k,ktr)=tr1do(k,ktr)
          enddo
        endif
        hm(k)=max(onemm,dp(k))/onem
        zm(k)=zgrid(k)
      end do !k

      do k=nlayer+1,kdm+1
        ka=min(k,kdm)
        difft( k)=0.0
        diffs( k)=0.0
        difftr(k)=0.0
        ghat(k)=0.0
        t1do(k)=temp(ka)
        s1do(k)=saln(ka)
        t1dn(k)=t1do(k)
        s1dn(k)=s1do(k)
        if (dotrcr .and. present(passv_tr)) then
          do ktr=1,numtr
            tr1do(k,ktr)=passv_tr(ka,ktr)
            tr1dn(k,ktr)=tr1do(k,ktr)
          end do
        end if
        zm(k)=zm(k-1)-0.001
      end do

!     if (vrbos) then
!         write (*,102) (nstep,k,i,hm(k),t1do(k),		&
!               t1dn(k),s1do(k),s1dn(k),k=1,nlayer)
!102    format(26x,'thick   t old   t ijo   s old   s ijo'	&
!          /(i8,i4,i8,2x,f9.2,4f8.3))
!     endif !test

! --- do rivers here because difs is also used for tracers.
!cc      if     (thkriv.gt.0.0 .and. rivers(i,1).ne.0.0) then
!cc        do k=1,nlayer-1
!cc          if     (-zm(k)+0.5*hm(k).lt.thkriv) then !interface<thkriv
!cc            diffs(k+1) = max(diffs(k+1),difriv)
!cc          endif
!cc        enddo !k
!cc      endif !river

! --- compute factors for coefficients of tridiagonal matrix elements.
!     tri(k=1:NZ,0) : dt/hwide(k)/ dzb(k-1)=z(k-1)-z(k)=dzabove)
!     tri(k=1:NZ,1) : dt/hwide(k)/(dzb(k  )=z(k)-z(k+1)=dzbelow)

      do k=1,nlayer
        dzb(k)=zm(k)-zm(k+1)
      enddo

      tri(1,1)=dlt/(hm(1)*dzb(1))
      tri(1,0)=0.
      do k=2,nlayer
        tri(k,1)=dlt/(hm(k)*dzb(k))
        tri(k,0)=dlt/(hm(k)*dzb(k-1))
      enddo

! --- solve the diffusion equation
! --- salflx, sswflx and surflx are positive into the ocean

! --- t solution

      ghatflux=-(surflx-sswflx)*thref/spcifh
      call tridcof(difft,tri,nlayer,tcu,tcc,tcl)
      call tridrhs(hm,t1do,difft,ghat,ghatflux,nlayer,rhs,dlt)
      call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,t1do,t1dn,difft)

! --- t-like tracer solution

      if (dotrcr .and. present(passv_tr)) then
        do ktr= 1,numtr
!cc        if     (bclin_frq(ktr).eq.2) then
          ghatflux=-(surflx-sswflx)*thref/spcifh
          call tridrhs(hm,tr1do(1,ktr),difft,ghat,ghatflux,	&
                       nlayer,rhs,dlt*trcfrq)
          call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,		&
                       tr1do(1,ktr),tr1dn(1,ktr),difft)
!cc        endif
        enddo
      end if

! --- s solution

      ghatflux=-salflx*thref
      call tridcof(diffs,tri,nlayer,tcu,tcc,tcl)
      call tridrhs(hm,s1do,diffs,ghat,ghatflux,nlayer,rhs,dlt)
      call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,s1do,s1dn,diffs)

! --- if surf.salinity gets too low, import some salt from lower layers
      if (s1dn(1).lt.fresh) then
        tosaln=s1dn(1)*hm(1)
        sumthk=        hm(1)
        print 105,'(mxkprfbij) warning: surf.salinity too low at',	&
          i,(k,hm(k),s1dn(k),k=1,5)
 105    format (a,i8/(i5,2f8.3))
        do k=2,kdm
          salin=tosaln/sumthk
          if (s1dn(k).gt.fresh .and. hm(k).gt.0. .and.			&
            tosaln+s1dn(k)*hm(k).gt.fresh*(sumthk+hm(k))) then
            frac=sumthk*(fresh-salin)/(hm(k)*(s1dn(k)-fresh))
!           if (frac.lt.0. .or. frac.gt.1.) then
!             print '(a,2i5,a,f9.4)','(mxkprfbij)',i,			&
!             '  illegal frac value:',frac
!             stop '(frac outside range)'
!           end if
            s1dn(k)=(1.-frac)*s1dn(k)+frac*fresh
            salin=fresh
            exit
          else          ! insufficient amount of salt in layer k
            tosaln=tosaln+s1dn(k)*hm(k)
            sumthk=sumthk+        hm(k)
            salin=tosaln/sumthk
          end if
        end do
        do ka=1,k-1
          s1dn(ka)=salin
        end do
        print 105,'(mxkprfbij) adjusted salinity profile at',		&
          i,(k,hm(k),s1dn(k),k=1,5)
      end if                            ! s1dn(1) < fresh

      if (vrbos) then
          write (*,103) (nstep,k,i,				&
          hm(max(1,k-1)),1.e4*difft(k),1.e4*diffs(k),		&
            ghat(k),k=1,nlayer+1)
          write (*,104) (nstep,k,i,				&
            hm(k),t1do(k),t1dn(k),s1do(k),s1dn(k),		&
            k=1,nlayer)
 103    format(26x,'thick    t_diff    s_diff   nonlocal'	&
           /(i8,i4,i8,1x,3f10.2,f11.6))
 104    format(26x,'thick   t_old   t_new   s_old   s_new'	&
           /(i8,i4,i8,2x,f9.2,4f8.3))
      endif

      if (dotrcr .and. present(passv_tr)) then
! --- tracer solution (using time step trcfrq*dlt)

        tri(1,1)=trcfrq*dlt/(hm(1)*dzb(1))
        tri(1,0)=0.
        do k=2,nlayer
          tri(k,1)=trcfrq*dlt/(hm(k)*dzb(k))
          tri(k,0)=trcfrq*dlt/(hm(k)*dzb(k-1))
        enddo

        call tridcof(difftr,tri,nlayer,tcu,tcc,tcl)

        do ktr= 1,numtr
         if (vrbos) write (*,'(a,i2/(8es10.3))')			&
           'old tracer',ktr,(tr1do(k,ktr),k=1,nlayer)
         ghatflux=0.
         call tridrhs(hm,tr1do(1,ktr),difftr,ghat,ghatflux,		&
                      nlayer,rhs,trcfrq*dlt)
         call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,			&
                      tr1do(1,ktr),tr1dn(1,ktr),difftr) 
         if (vrbos) write (*,'(a,i2/(8es10.3))')			&
           'new tracer',ktr,(tr1dn(k,ktr),k=1,nlayer)
        end do
      end if			! dotrcr

! --- check conservation of column integrals
      totemo=t1do(1)*dp(1)
      totemn=t1dn(1)*dp(1)
      tosalo=s1do(1)*dp(1)
      tosaln=s1dn(1)*dp(1)
      if (dotrcr .and. present(passv_tr)) then
        do ktr=1,numtr
          totrco(ktr)=tr1do(1,ktr)*dp(1)
          totrcn(ktr)=tr1dn(1,ktr)*dp(1)
          trscal(ktr)=max(abs(tr1do(1,ktr)),abs(tr1dn(1,ktr)))
        end do
      end if                  ! dotrcr

      do k=2,kdm
        totemo=totemo+t1do(k)*dp(k)
        totemn=totemn+t1dn(k)*dp(k)
        tosalo=tosalo+s1do(k)*dp(k)
        tosaln=tosaln+s1dn(k)*dp(k)
        if (dotrcr .and. present(passv_tr)) then
          do ktr=1,numtr
            totrco(ktr)=totrco(ktr)+tr1do(k,ktr)*dp(k)
            totrcn(ktr)=totrcn(ktr)+tr1dn(k,ktr)*dp(k)
            trscal(ktr)=max(trscal(ktr),abs(tr1do(k,ktr)),		&
                                        abs(tr1dn(k,ktr)))
          end do
        end if                        !  dotrcr
      end do
      tndcyt=totemn-totemo
      tndcys=tosaln-tosalo
      totemn=10.*pbot
      tosaln=35.*pbot
      if (abs(tndcyt).gt.acurcy*totemn) write (*,101) i,		&
        '  mxkpp - bad temp.intgl.',totemo,tndcyt,tndcyt/totemn
      if (abs(tndcys).gt.acurcy*tosaln) write (*,101) i,		&
        '  mxkpp - bad saln.intgl.',tosalo,tndcys,tndcys/tosaln

      if (dotrcr .and. present(passv_tr)) then
        do ktr=1,numtr
          tndcyt=totrcn(ktr)-totrco(ktr)
          if (abs(tndcyt).lt.1.e-33) tndcyt=0.
          totemn=trscal(ktr)*pbot
          if (abs(tndcyt).gt.acurcy*totemn) write (*,101) i,		&
          '  mxkpp - bad trcr.intgl.',totrco(ktr),tndcyt,tndcyt/totemn
        end do
      end if                         ! dotrcr
 101  format (i8,a,2es16.8,es9.1)

! --- adjust t, s, th, arrays

      if     ((tofset.eq.0.0 .and. sofset.eq.0.0 ) .or.			&
              (mod(nstep  ,tsofrq).ne.0 .and.				&
               mod(nstep+1,tsofrq).ne.0      )     ) then
        do k=1,klist
          temp(k)=t1dn(k)
          saln(k)=s1dn(k)
          if (dotrcr .and. present(passv_tr)) then
            do ktr= 1,numtr
              passv_tr(k,ktr)=tr1dn(k,ktr)
            end do !ktr
          end if
        end do !k
      else  !include [ts]ofset drift correction
        do k=1,klist
          temp(k)=t1dn(k) + dlt*max(2,tsofrq)*tofset
          saln(k)=s1dn(k) + dlt*max(2,tsofrq)*sofset
          if (dotrcr .and. present(passv_tr)) then
            do ktr= 1,numtr
              passv_tr(k,ktr)=tr1dn(k,ktr)
            end do !ktr
          end if
        enddo !k
      endif !without:with [ts]ofset

      return
      end subroutine mxkppbij

      subroutine mxkppcijuv(nstep,pbot,dpmixl,dp,vcty,u,v,dlt,i,vrbos)

! --- hycom version 1.0
! -------------------------------------------------------------------------
! --- k-profile vertical diffusion, single j-row, momentum at u grid points
! --- vertical coordinate is z negative below the ocean surface
! -------------------------------------------------------------------------

  use fimnamelist    ,only: kdm
  use hycom_control  ,only: iocnmx
  use hycom_constants,only: onem,onemm,tencm
  implicit none
  integer,intent(IN) :: i,nstep
  logical,intent(IN) :: vrbos
  real   ,intent(IN) :: dlt
  real*8,intent(IN)    :: vcty(kdm+1)
  real*8,intent(IN)    :: pbot,dpmixl
  real*8,intent(IN)    :: dp (kdm)          ! layer thickness
  real*8,intent(INOUT) :: u  (kdm)          ! u velocity
  real*8,intent(INOUT) :: v  (kdm)          ! v velocity

! --- local 1-d arrays for matrix solution
  real*8 u1do(kdm+1),u1dn(kdm+1),v1do(kdm+1),v1dn(kdm+1),diffm(kdm+1)	&
         ,zm(kdm+1),hm(kdm),dzb(kdm)

! --- tridiagonal matrix solution arrays
  real*8 tri(kdm,0:1)       ! dt/dz/dz factors in trid. matrix
  real*8 tcu(kdm),       &  ! upper coeff for (k-1) on k line of trid.matrix
         tcc(kdm),       &  ! central ...     (k  ) ..
         tcl(kdm),       &  ! lower .....     (k-1) ..
         rhs(kdm)           ! right-hand-side terms

  real*8 p,dpthmx
  integer k,nlayer

  nlayer=1
  p=0.
! --- dpthmx = depth range over which to mix momentum
  dpthmx=pbot-tencm

! <><><><><><><><>  p a r t i a l   c o l u m n   k p p  <><><><><><><><>
  if (mod(iocnmx,4).eq.2) dpthmx=min(dpthmx,dpmixl)
! <><><><><><><><>  p a r t i a l   c o l u m n   k p p  <><><><><><><><>

! do k=1,kdm+1    ! bug
  do k=1,kdm
    if (p.lt.dpthmx) then
      diffm(k+1)=vcty(k+1)
      u1do(k)=u(k)
      v1do(k)=v(k)
      hm(k)=max(onemm,dp(k))/onem
      if (k.eq.1) then
        zm(k)=-.5*hm(k)
      else
        zm(k)=zm(k-1)-.5*(hm(k-1)+hm(k))
      endif
      p=p+pbot
      nlayer=k
    else if (k.eq.nlayer+1) then
      diffm(k)=0.
      u1do(k)=u1do(k-1)
      v1do(k)=v1do(k-1)
      zm(k)=zm(k-1)-.5*hm(k-1)
      exit
    endif
  enddo

! --- compute factors for coefficients of tridiagonal matrix elements.
  do k=1,nlayer
    dzb(k)=zm(k)-zm(k+1)
  enddo

  tri(1,1)=dlt/(hm(1)*dzb(1))
  tri(1,0)=0.
  do k=2,nlayer
    tri(k,1)=dlt/(hm(k)*dzb(k))
    tri(k,0)=dlt/(hm(k)*dzb(k-1))
  enddo

! --- solve the diffusion equation
  call tridcof(diffm,tri,nlayer,tcu,tcc,tcl)
  do k=1,nlayer
    rhs(k)= u1do(k)
  enddo
  call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,u1do,u1dn,diffm)

  do k=1,nlayer
    rhs(k)= v1do(k)
  enddo
  call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,v1do,v1dn,diffm)
  do k=1,nlayer
    u(k)=u1dn(k)
    v(k)=v1dn(k)
  enddo

  if (vrbos) then
    write (*,106) (nstep,k,i,hm(k),u1do(k),u1dn(k),k=1,nlayer)
    write (*,107) (nstep,k,i,hm(k),v1do(k),v1dn(k),k=1,nlayer)
  endif
  return
 106  format(25x,'thick   u old   u new'/(i8,i3,i8,1x,f10.3,2f8.3))
 107  format(25x,'thick   v old   v new'/(i8,i3,i8,1x,f10.3,2f8.3))
  end subroutine mxkppcijuv

  subroutine wscale(i,zlevel,dnorm,bflux,wm,ws,isb,ustar,ustarb)
! -------------------------------------------------------------------------
! --- subroutine to compute turbulent velocity scales for kpp mixing scheme
! --- vertical coordinate is z negative below the ocean surface
! -------------------------------------------------------------------------
 
  use hycom_kpp_constants, only: nzehat, nustar, wmt, wst, c11,		&
      deltaz, deltau, umin,zmin,vonk,zmax,zmin,umin
  implicit none
  integer,intent(IN) :: i
 
  real,   intent(IN) :: zlevel,dnorm,bflux,ustar,ustarb
  integer,intent(IN) :: isb
  real,  intent(OUT) :: wm,ws

! --- isb determines whether the calculation is for the surface or
! --- bottom boundary layer

  real    zdiff,udiff,zfrac,ufrac,ust,wam,wbm,was,wbs,ucube,zehat
  integer iz,izp1,ju,jup1

! --- use lookup table for zehat < zmax  only;  otherwise use stable formulae

  if (isb.eq.1) then
    ust=ustar
    zehat=-vonk*dnorm*zlevel*bflux
  else
    ust=ustarb
    zehat= vonk*dnorm*zlevel*bflux
  endif
  if (zehat.le.zmax) then
    zdiff=zehat-zmin
    iz=int(zdiff/deltaz)
    iz=max(min(iz,nzehat),0)
    izp1=iz+1

    udiff=ust-umin
    ju=int(udiff/deltau)
    ju=max(min(ju,nustar),0)
    jup1=ju+1

    zfrac=zdiff/deltaz-iz
    ufrac=udiff/deltau-ju

    wam=(1.-zfrac)*wmt(iz,jup1)+zfrac*wmt(izp1,jup1)
    wbm=(1.-zfrac)*wmt(iz,ju  )+zfrac*wmt(izp1,ju  )
    wm =(1.-ufrac)*wbm         +ufrac*wam

    was=(1.-zfrac)*wst(iz,jup1)+zfrac*wst(izp1,jup1)
    wbs=(1.-zfrac)*wst(iz,ju  )+zfrac*wst(izp1,ju  )
    ws =(1.-ufrac)*wbs         +ufrac*was
  else
    ucube=ust**3
    wm=vonk*ust*ucube/(ucube+c11*zehat)
    ws=wm
  endif

  return
  end subroutine wscale


  subroutine inikpp
! -------------------------------------------------------------------
! --- initialize large, mc williams, doney kpp vertical mixing scheme
! -------------------------------------------------------------------

  use hycom_constants, only: epsil
  use hycom_kpp_constants, only: nzehat, nustar, wmt, wst,		&
      vonk,zmin,zmax,umin,umax,epsilon,vtc,cs0,cg,dp0enh,deltau,	&
      deltaz,c11,ricr,cstar,dp00,difs0,difm0,qdif0,qdifiw,		&
      difmiw,difsiw
  implicit none

  real zehat,zeta,usta,afourth,athird,ahalf
  parameter (afourth=.25,athird=1./3,ahalf=.5)
  integer i,j

  real, parameter :: am=1.257,cm=8.380,c22=16.0,zetam=-0.2,		&
	             as=-28.86, c33=16.0, zetas=-1.0

! --- 'vonk'      = von karman constant
! --- 'zmin,zmax' = zehat limits for velocity scale lookup table, m**3/s**3
! --- 'umin,umax' = ustar limits for velocity scale lookup table
! --- 'epsilon'   = vertical coordinate scale factor

  vonk   =  0.4
  zmin   = -0.4e-6
  zmax   =  0.0
  umin   =  0.0
  umax   =  0.16
  epsilon=  0.1

! --- construct the velocity-scale lookup tables

  deltaz = (zmax-zmin)/(nzehat+1)
  deltau = (umax-umin)/(nustar+1)

  do i=0,nzehat+1
  zehat=deltaz*i+zmin
  do j=0,nustar+1
  usta=deltau*j+umin
  zeta=zehat/(usta**3+epsil)
  if (zehat.ge.0.) then
    wmt(i,j)=vonk*usta/(1.+c11*zeta)
    wst(i,j)=wmt(i,j)
  else
    if (zeta.gt.zetam) then
      wmt(i,j)=vonk*usta*(1.-c22*zeta)**afourth
    else
      wmt(i,j)=vonk*(am*usta**3-cm*zehat)**athird
    endif
    if (zeta.gt.zetas) then
      wst(i,j)=vonk*usta*(1.-c33*zeta)**ahalf
    else
      wst(i,j)=vonk*(as*usta**3-cs0*zehat)**athird
    endif
  endif
  enddo
  enddo

! --- set derived constants
  vtc=sqrt(.2/cs0/epsilon)/vonk**2/ricr
  cg=cstar*vonk*(cs0*vonk*epsilon)**athird
  dp0enh=2.0*dp00

  qdif0 =difm0 /difs0
  qdifiw=difmiw/difsiw

!diag write(*,*) 'shown below: lookup table wst'
!diag call zebra(wst(0,0),nzehat+2,nzehat+2,nustar+2)
!diag write(*,*) 'shown below: lookup table wmt'
!diag call zebra(wmt(0,0),nzehat+2,nzehat+2,nustar+2)

  call initurb

  return
  end subroutine inikpp

  subroutine initurb
  use hycom_constants, only: onem,odepth,wet
  use hycom_kpp_constants, only: jerlv0, betard, betabl, redfac
  use module_constants,only: deg_lat
  use module_control,only: nip
  use hycom_variables, only: jerlov
  implicit none
  integer i

! --- map shallow depths ('brown' water) to high jerlov numbers
!SMS$PARALLEL(dh,i) BEGIN
!$OMP PARALLEL DO
  do i=1,nip
  if (wet(i) > 0 ) then
    jerlov(i)=6-max(1,min(5,int(odepth(i)/15.0)))
    jerlov(i)=max(jerlv0,jerlov(i))
! --- reduce turbidity of subtropical oceans
    if (abs(abs(deg_lat(i))-25.).lt.20.) jerlov(i)=1
  end if
  end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END

! --- red and blue light extinction coefficients (1/pressure units)
! --- for jerlov water types 1 to 5 - fraction of penetrating red light
  betard(1) = 1.0/( 0.35*onem)
  betard(2) = 1.0/( 0.6 *onem)
  betard(3) = 1.0/( 1.0 *onem)
  betard(4) = 1.0/( 1.5 *onem)
  betard(5) = 1.0/( 1.4 *onem)
  betabl(1) = 1.0/(23.0 *onem)
  betabl(2) = 1.0/(20.0 *onem)
  betabl(3) = 1.0/(17.0 *onem)
  betabl(4) = 1.0/(14.0 *onem)
  betabl(5) = 1.0/( 7.9 *onem)
  redfac(1) = 0.58
  redfac(2) = 0.62
  redfac(3) = 0.67
  redfac(4) = 0.77
  redfac(5) = 0.78
  print *,'turbidity parameters initialized'

  return
  end subroutine initurb


  subroutine tridcof(diff,tri,nlayer,tcu,tcc,tcl)
! ------------------------------------------------------------------
! --- matrix inversion subroutines for implicit solution of vertical
! --- diffusion equation - tri-diagonal matrix
! ------------------------------------------------------------------

  use fimnamelist, only : kdm
  implicit none

! --- compute coefficients for tridiagonal matrix (dimension=kdm).
! --- Note: tcu(1) = 0. and tcl(kdm+1) = 0. are necessary conditions.
!     
! --- input
  real*8 ,intent(IN) :: diff(kdm+1)	! diffusivity profile on interfaces
  real*8, intent(IN) :: tri(kdm,0:1)	! dt/dz/dz factors in trid. matrix
  integer,intent(IN) :: nlayer

! --- output
  real*8,intent(OUT) ::	&
      tcu(kdm),		&! upper coeff. for (k-1) on k line of trid.matrix
      tcc(kdm),		&! central ...      (k  ) ..
      tcl(kdm)		 ! lower .....      (k-1) ..

  integer :: k

! --- common tridiagonal factors

! --- in the surface layer
  tcu(1)=0.
  tcc(1)=1.+tri(1,1)*diff(2)		! 1.+ delt1/h(1)/dzb(1)*diff(2)
  tcl(1)=  -tri(1,1)*diff(2)		!   - delt1/h(1)/dzb(1)*diff(2)

! --- inside the domain
  do 10 k=2,nlayer
    tcu(k)=  -tri(k,0)*diff(k  )
    tcc(k)=1.+tri(k,1)*diff(k+1)+tri(k,0)*diff(k)
    tcl(k)=  -tri(k,1)*diff(k+1)
 10 continue

! --- in the bottom layer
  tcl(nlayer)= 0.
  return
  end subroutine tridcof


  subroutine tridrhs(h,yo,diff,ghat,ghatflux,nlayer,rhs,delt1)

  use fimnamelist, only : kdm
  implicit none

! --- compute right hand side of tridiagonal matrix for scalar fields:
! --- =  yo (old field) 
! ---  + flux-divergence of ghat
! ---  + flux-divergence of non-turbulant fluxes

! --- note: if surface and bottom fluxes are nonzero, the following must apply
! ---    sfc. lyr. needs +delt1/h(1)*surfaceflux
! ---    bot. lyr. needs +delt1/h(nlayer)*diff(nlayer+1)/
! ---                     dzb(nlayer)*yo(nlayer+1)

! --- input
  real*8,intent(IN) ::		&
          h(kdm),		& ! layer thickness
          yo(kdm+1),		& ! old profile
          diff(kdm+1),		& ! diffusivity profile on interfaces
          ghat(kdm+1),		& ! ghat turbulent flux   
          ghatflux		! surface flux for ghat: includes solar flux
  real   ,intent(IN) :: delt1	! time step
  integer,intent(IN) :: nlayer

! --- output
  real*8,intent(OUT) :: rhs(kdm)	! right hand side

  integer k 

! --- in the top layer
  rhs(1)=yo(1)+delt1/h(1)*(ghatflux*diff(2)*ghat(2))

! --- inside the domain 
  do 10 k=2,nlayer-1
  rhs(k)=yo(k)+delt1/h(k)*     &
      (ghatflux*(diff(k+1)*ghat(k+1)-diff(k)*ghat(k)))
 10   continue

! --- in the bottom layer     
  k=nlayer
  rhs(k)=yo(k)+delt1/h(k)*     &
       (ghatflux*(diff(k+1)*ghat(k+1)-diff(k)*ghat(k)))

  return
  end subroutine tridrhs


  subroutine tridmat(tcu,tcc,tcl,nlayer,h,rhs,yo,yn,diff)

! --- solve tridiagonal matrix for new vector yn, given right hand side
! --- vector rhs.

! --- note: if surface and bottom fluxes are nonzero, the following must apply
! ---    surface layer needs +delt1*surfaceflux/(h(1)*bet)
! ---    bottom  layer needs +tri(nlayer,1)*diff(nlayer+1)*yo(nlayer+1))/bet

  use fimnamelist, only : kdm,itest
  implicit none
  integer k

! --- input
  real*8 tcu (kdm),    & ! upper coeff. for (k-1) on k line of tridmatrix
         tcc (kdm),    & ! central ...      (k  ) ..
         tcl (kdm),    & ! lower .....      (k-1) ..
          h  (kdm),    & ! layer thickness
          rhs(kdm),    & ! right hand side
           yo(kdm+1),  & ! old field
         diff(kdm+1),  & ! diffusivity profile
          gam(kdm)  ,  & ! temporary array for tridiagonal solver
          bet            ! ...
  integer nlayer

! --- output
  real*8 yn(kdm+1)      ! new field

! --- solve tridiagonal matrix.
  bet=tcc(1)
  yn(1)=rhs(1)/bet                               ! surface
  do 21 k=2,nlayer
    gam(k)=tcl(k-1)/bet
    bet=tcc(k)-tcu(k)*gam(k)
    if (bet.eq.0.) then
      write(*,*) 
      write(*,*) '** algorithm for solving tridiagonal matrix fails'
      write(*,*) '** bet=',bet,itest
      write(*,*) '** k=',k,' tcc=',tcc(k),' tcu=',tcu(k),' gam=',gam(k)
      stop '(tridmat)'
!     bet=1.E-12
    endif
    yn(k) =      (rhs(k)  - tcu(k)  *yn(k-1)  )/bet
!     to avoid "Underflow" at single precision on the sun
!     yni   =      (rhs(k)  - tcu(k)  *yn(k-1)  )/bet 
!     if (yni.lt.0.) then
!       yn(k) = min( (rhs(k)  - tcu(k)  *yn(k-1)  )/bet ,-1.E-12 )
!     else
!       yn(k) = max( (rhs(k)  - tcu(k)  *yn(k-1)  )/bet , 1.E-12 )
!     endif
 21   continue

    do 22 k=nlayer-1,1,-1
      yn(k)=yn(k)-gam(k+1)*yn(k+1)
 22 continue

    yn(nlayer+1)=yo(nlayer+1)

    return
    end subroutine tridmat

  subroutine mxkpp(nstep,leapn,theta,dens,dp,pres,temp,saln,uvel,vvel,	&
                   surflx,sswflx,pmevap,ustar,ustarb,jerlov,deg_lat,	&
                   corio,dpmx)

  use module_control  ,only: nip
  use module_constants,only: grvity,perm
  use fimnamelist     ,only: kdm,itest,diag_intvl
  use hycom_control   ,only: numtr,iocnmx,bclin_frq
  use hycom_constants ,only: wet,thref,onem,onecm,onemm,batrop,spcifh
  use hycom_sigetc
  use findmaxmin1
  use findmaxmin2
  use findmaxmin3

  implicit none
  integer,intent(IN) :: nstep                   ! model time step
  integer,intent(IN) :: leapn                   ! leapfrog time slot (n='new')
  real   ,intent(IN) :: theta(kdm)              ! target densities
!SMS$DISTRIBUTE  (dh,1) BEGIN
  real,intent(IN)    :: surflx(nip)             ! rad.+sensib+latnt heat flux
  real,intent(IN)    :: sswflx(nip)             ! net shortwave heat flux
  real,intent(IN)    :: pmevap(nip)             ! downwd freshwater flux (p-e)
  real,intent(IN)    :: ustar (nip)             ! friction velocity
  real,intent(IN)    :: ustarb(nip)             ! bottom friction velocity
  integer,intent(IN) :: jerlov(nip)             ! jerlov water type 1-5
  real,intent(IN)    :: deg_lat(nip)            ! latitude in degree
  real,intent(IN)    :: corio(nip)              ! Coriolis acceleration
  real,intent(INOUT) :: dpmx(nip)               ! mxlyr depth (m)
!
  real*8 :: umix(nip),vmix(nip),thmix(nip)
  real*8 :: tmix(nip),smix(nip),dpmixl(nip)
  real*8 :: coltemo(nip),colsalo(nip),coltemn(nip),colsaln(nip)
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE  (dh,2) BEGIN
  real,intent(INOUT) :: dens    (kdm,nip  )     ! density
  real,intent(IN)    :: dp      (kdm,nip,2)     ! layer thickness
  real,intent(IN)    :: pres    (kdm+1,nip)     ! interface pressure
  real,intent(INOUT) :: temp    (kdm,nip)       ! temperature
  real,intent(INOUT) :: saln    (kdm,nip)       ! salinity
  real,intent(INOUT) :: uvel  (kdm,nip,2)       ! u velocity
  real,intent(INOUT) :: vvel  (kdm,nip,2)       ! v velocity
!SMS$DISTRIBUTE END

  integer :: k,i,klist
  character :: string*16
  logical :: dotrcr,vrbos

  real*8:: t1d(kdm),s1d(kdm),u1d(kdm),v1d(kdm),dp1d(kdm),p1d(kdm+1),theta8(kdm)
  real*8:: zgrid(kdm+1),vcty(kdm+1),dift(kdm+1),difs(kdm+1),ghats(kdm+1),pbot
  real  ::  salflx,delp,dlt,dtem,dsal

  dotrcr=.false.
  theta8(:)=theta(:)
  dpmixl(:)=0.
  t1d(:)=0.
  s1d(:)=0.
  tmix(:)=0.
  smix(:)=0.

  dlt=batrop*bclin_frq				! routine-specific time step
  if (mod(nstep,diag_intvl).eq.0) then
    call findmxmn1(temp(1,:),nip,'(mxkpp1) sst',wet)
    call findmxmn1(saln(1,:),nip,'(mxkpp1) sss',wet)
  end if

!SMS$PARALLEL (dh,i) BEGIN
!JR All args on 2nd OMP line are private due to intent(out) in mxkppaij and intent(in) in mxkppbij
!$OMP PARALLEL DO PRIVATE (vrbos,k,t1d,s1d,u1d,v1d,dp1d,p1d,pbot,salflx,dtem,dsal, &
!$OMP&                     zgrid,vcty,difs,dift,ghats,klist)
  do i=1,nip
    if (wet(i) > 0 ) then
    vrbos=i.eq.itest .and. mod(nstep,diag_intvl).eq.0

    coltemo(i)=0.
    colsalo(i)=0.
    do k=1,kdm
      t1d(k)=temp(k,i)
      s1d(k)=saln(k,i)
      u1d(k)=uvel(k,i,leapn)
      v1d(k)=vvel(k,i,leapn)
      dp1d(k)=dp(k,i,leapn)
      p1d(k)=pres(k,i)
      coltemo(i)=coltemo(i)+t1d(k)*dp1d(k)
      colsalo(i)=colsalo(i)+s1d(k)*dp1d(k)
    end do
    p1d(kdm+1)=pres(kdm+1,i)
    pbot=p1d(kdm+1)
    salflx = -s1d(1)*pmevap(i)/thref

    if (vrbos) then
     write(*,'(i8,i7,3x,a)') nstep,i,					&
      'mxkpp IN:   thkns    temp    saln     u (cm/s) v'
     write(*,'(24x,i3,f8.1,4f8.3)')(k,dp1d(k)/onem,t1d(k),s1d(k),	&
      u1d(k)*100.,v1d(k)*100.,k=1,kdm)
     write (*,'(i8,i7,a,3es11.3)') nstep,i,' energy & salt flux, P-E:',	&
      surflx(i),salflx,pmevap(i)
    end if

    call mxkppaij (nstep,t1d,s1d,u1d,v1d,dp1d,p1d,theta8,		& !input
                   surflx(i),sswflx(i),salflx,ustar(i),ustarb(i),	& !input
                   jerlov(i),deg_lat(i),corio(i),0.,dlt,perm(i),vrbos,	& !input
                   zgrid,dift,difs,vcty,ghats,dpmixl(i),klist)		!output
    call mxkppbij (nstep,zgrid,dift,difs,ghats,dp1d,pbot,		& !input
                   surflx(i),sswflx(i),salflx,dlt,klist,perm(i),vrbos,	& !input
                   dotrcr,deg_lat(i),					& !input
                   t1d,s1d)						!output
    call mxkppcijuv(nstep,pbot,dpmixl(i),dp1d,vcty,u1d,v1d,dlt,perm(i),	&
                    vrbos)
    coltemn(i)=0.
    colsaln(i)=0.
    do k=1,kdm
      temp(k,i)=t1d(k)
      saln(k,i)=s1d(k)
      dens(k,i)=sigocn(t1d(k),s1d(k))               ! update density here
      uvel(k,i,leapn)=u1d(k)
      vvel(k,i,leapn)=v1d(k)
      coltemn(i)=coltemn(i)+t1d(k)*dp1d(k)
      colsaln(i)=colsaln(i)+s1d(k)*dp1d(k)
    end do

    if (vrbos) then
     print '(i8,a,2es13.6,2es9.2/8x,a,2es13.6,2es9.2)',i,		&
      ' (mxkpp) Told/new,diff,flux=',coltemo(i),coltemn(i),		&
     coltemn(i)-coltemo(i),surflx(i)*grvity/spcifh*dlt,			&
     '         Sold/new,diff,flux=',colsalo(i),colsaln(i),		&
     colsaln(i)-colsalo(i),salflx*grvity*dlt

     write(*,'(2i8,a,3f7.1/a)') nstep,i,'  dpmixl,heatfx,shwflx=',	&
     dpmixl(i)/onem,surflx(i),sswflx(i),				&
     'mxkpp OUT: thkns    temp    saln    dens     u (cm/s) v   dift   difs vcty*10^4'
     write(*,'(i8,f8.1,5f8.3,3f7.2)') (k,dp1d(k)/onem,t1d(k),s1d(k),	&
     dens(k,i),u1d(k)*100.,v1d(k)*100.,					&
     dift(k+1)*1.e4,difs(k+1)*1.e4,vcty(k+1)*1.e4,k=1,kdm)
    end if    ! if vrbos

! --- to overcome numerical precision problems, restore column integrals
    dtem=(coltemo(i)-coltemn(i)+surflx(i)*grvity/spcifh*dlt)/pres(kdm+1,i)
    dsal=(colsalo(i)-colsaln(i)+salflx   *grvity       *dlt)/pres(kdm+1,i)
    do k=1,kdm
      temp(k,i)=temp(k,i)+dtem
      saln(k,i)=saln(k,i)+dsal
    end do

    end if    ! if wet
  end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END

  if (mod(nstep,diag_intvl).eq.0) then
    call findmxmn1(temp(1,:),nip,'(mxkpp2) sst',wet)
    call findmxmn1(saln(1,:),nip,'(mxkpp2) sss',wet)

! --- calculate bulk mixed layer t, s, theta, u & v

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (p1d,delp)
    do i=1,nip
     if (wet(i) > 0 ) then
      if (iocnmx.le.3) then
        p1d(kdm+1)=pres(kdm+1,i)
        dpmixl(i)=min(dpmixl(i),p1d(kdm+1))
      end if
      p1d(2)=pres(2,i)
      delp=min(p1d(2),dpmixl(i))
      tmix(i)=delp*temp(1,i)
      smix(i)=delp*saln(1,i)
      umix(i)=delp*uvel(1,i,2)
      vmix(i)=delp*vvel(1,i,2)
      do k=2,kdm
      p1d(k+1)=pres(k+1,i)
      delp=min(p1d(k+1),dpmixl(i))-min(p1d(k),dpmixl(i))
      tmix(i)=tmix(i)+delp*temp(k,i)
      smix(i)=smix(i)+delp*saln(k,i)
      umix(i)=umix(i)+delp*uvel(k,i,2)
      vmix(i)=vmix(i)+delp*vvel(k,i,2)
      enddo
      tmix(i)=tmix(i)/dpmixl(i)
      smix(i)=smix(i)/dpmixl(i)
      umix(i)=umix(i)/dpmixl(i)
      vmix(i)=vmix(i)/dpmixl(i)
      thmix(i)=sigocn(tmix(i),smix(i))

      if (vrbos) then
       write (*,'(i8,2i5,a,f9.2)') nstep,k,i,'   dpmixl =',dpmixl(i)/onem
      end if				! if vrbos
     end if				! if wet
    end do				! do i
!$OMP END PARALLEL DO
!SMS$PARALLEL END
  end if				! if diagno

  dpmx(:)=dpmixl(:)/onem

!!!SMS$PARALLEL (dh,i) BEGIN
!!   do i=1,nip
!!   if (wet(i) > 0 ) then
!!     dpmxav(i)=dpmxav(i)+dpmixl(i)
!!   end if 
!!   end do 
!!!SMS$PARALLEL END

  return
  end subroutine mxkpp
end module hycom_mxkpp
