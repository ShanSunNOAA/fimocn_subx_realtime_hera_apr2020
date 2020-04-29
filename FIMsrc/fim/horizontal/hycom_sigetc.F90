module hycom_sigetc
contains

  real*8 function sigocn(t,s)

! --- density anomaly (i.e. actual density minus reference density)

  use hycom_constants ,only: pref
  implicit none
  real*8, intent(in)  :: t,s
! --- based on Jackett et al. 2006, J. of atm & oceanic technology.
! --- rho(t=25, s=35, 2000) = 31. 650 560 565 76
! --- rho(t=20, s=20, 1000) = 17. 728 868 019 64
! --- rho(t=12, s=40, 8000) = 62. 952 798 206 31
        
  real*8, parameter ::							&
        k01=9.9984085444849347e2,    a1=7.3471625860981584,		&
        a2=-5.3211231792841769e-2,   a3=3.6492439109814549e-4,		&
        b1=2.588057102399139,        b2=-6.7168282786692355e-3,		&
        b3=1.9203202055760151e-3,    cc1=1.1798263740430364e-2,		&
        cc2=9.8920219266399117e-8,   cc3=4.699664277175473e-6,		&
        cc4=-2.5862187075154352e-8,  cc5=-3.2921414007960662e-12

  real*8, parameter:: k02=1.0,						&
        d1=7.2815210113327091e-3,    d2=-4.4787265461983921e-5,		&
        d3=3.3851002965802430e-7,    d4=1.3651202389758572e-10,		&
        e1=1.7632126669040377e-3,    e2=-8.8066583251206474e-6,		&
        e3=-1.8832689434804897e-10,  e4=5.7463776745432097e-6,		&
        e5=1.4716275472242334e-9,    f1=6.7103246285651894e-6,		&
        f2=-2.4461698007024582e-17,  f3=-9.1534417604289062e-18

  real*8 pn,pd,p

  p=pref*1.e-4					! p in dbar
  pn=k01+t*(a1+t*(a2+a3*t))+s*(b1+b2*t+b3*s)			&
     +p*(cc1+cc2*t*t+cc3*s+p*(cc4+cc5*t*t))

  pd=k02+t*(d1+d2*t+t*t*(d3+d4*t))				&
        +s*(e1+e2*t+e3*t*t*t+sqrt(max(0.,s))*(e4+e5*t*t))	&
        +p*(f1+p*t*(f2*t*t+f3*p))

  sigocn = pn/pd - 1000.			! subtract reference density
  return
  end function sigocn
!
  real*8 function dsigdt(t,s)
  use hycom_constants ,only: c1,c2,c3,c4,c5,c6,c7,c8,c9
  implicit none
  real*8 t,s
  dsigdt=(c2+s*(c5+c9*s)+2.*t*(c4+c7*s+1.5*c6*t))           !  SI c9
  return
  end function dsigdt
!
  real*8 function dsigds(t,s)
  use hycom_constants ,only: c1,c2,c3,c4,c5,c6,c7,c8,c9
  implicit none
  real*8 t,s
  dsigds=c3+2*s*(c8+c9*t)+t*(c5+t*c7)  ! SI c9
  return
  end function dsigds
!
  real function tofsig(sigm,salin)
!
! --- temp (deg c) as a function of sigma and salinity (this routine mimics
! --- the statement function of the same name in state_eqn.h.)
!
  use hycom_constants ,only: c1,c2,c3,c4,c5,c6,c7,c8,c9
  implicit none
  real:: sigm,salin,sqq,a0,a1,a2,cubr,cubq,cuban,cubrl,cubim
  real, parameter:: athird=1./3.
!
  a0=(c1+salin*(c3+c8*salin))/c6
  a1=(c2+salin*(c5+c9*salin))/c6
  a2=(c4+c7*salin)/c6
  cubq=athird*a1-(athird*a2)**2
  cubr=athird*(.5*a1*a2-1.5*(a0-sigm/c6))-(athird*a2)**3
!
! --- if q**3+r**2>0, water is too dense to yield real root at given
! --- salinitiy. setting q**3+r**2=0 in that case is equivalent to
! --- lowering sigma until a double real root is obtained.
!
  cuban=athird*atan2(sqrt(max(0.,-(cubq**3+cubr**2))),cubr)
  sqq=sqrt(-cubq)
  cubrl=sqq*cos(cuban)
  cubim=sqq*sin(cuban)
  tofsig=-cubrl+sqrt(3.)*cubim-athird*a2
!!!      if (abs(sig(tofsig,salin)-sigm).gt.1.e-6) write (*,100)
!!!     .   tofsig,salin,sigm,sig(tofsig,salin)
 100  format ('tofsig,sal,old/new sig =',2f9.3,3p,2f9.3)
  return
  end function tofsig
!
  real*8 function sofsig(sigma,tem)
  use hycom_constants ,only: c1,c2,c3,c4,c5,c6,c7,c8,c9
  implicit none
  real*8 :: sigma,tem,aa,bb,cc

  aa=c8+c9*tem
  bb=c3+c5*tem+c7*tem*tem
  cc=c1+c2*tem+c4*tem*tem+c6*tem*tem*tem-sigma
  sofsig=(-bb+sqrt(bb*bb-4.*aa*cc))/(2.*aa)
  return
  end function sofsig
!
  real function qsatur(t)
! --- saturation specific humidity (lowe, j.appl.met., 16, 100-103, 1976)
  implicit none
  real, intent(in) :: t
  qsatur=.622e-3*(6.107799961e+00+t*(4.436518521e-01      &
        +t*(1.428945805e-02+t*(2.650648471e-04            &
        +t*(3.031240396e-06+t*(2.034080948e-08            &
        +t* 6.136820929e-11))))))
  return
  end function qsatur
!
  real*8 function sigma_pref_star(t,s,p)
  use hycom_constants ,only: c1,c2,c3,c4,c5,c6,c7,c8,c9,pref
  implicit none
  real*8, intent(in) :: s,t,p       ! pres in unit of dbar
  real*8 :: kappa,sig9,sdif,tdif
!
! --- coefficients for kappa^(theta) fit towards JM06
! --- new values (w.r.t. t-toff,s-soff) from Shan Sun Oct.2012
  real*8, parameter ::   &  ! 8 coeffi. for kappa_theta, fit range: [-2:30],[30:38],[0:5500m]
    qttt=-3.282380E-09,  qtt= 4.459121E-07,   qt=-2.793139E-05,		&
    qs=-7.587727E-06,    qst= 7.516168E-08,   qpt= 9.971622E-10,	&
    qpst= 9.377664E-12,  qptt=-1.338658E-11,				&
    reftem= 3.0d0, refsal= 3.45d1

    if (t<-5. .or. t>35. .or. s<0. .or. s>42. .or. p*1.e-4 >10000.) then
      print *,' t(deg),s(ppt),p(dbar)=',t,s,p*1.e-4
      stop 'wrong t,s,p in sig_jm06'
    end if

    sig9=c1+s*(c3+c8*s)+t*(c2+c5*s+t*(c4+c7*s+c6*t))+c9*t*s*s
    tdif=t-reftem
    sdif=s-refsal
    kappa=(tdif*(qt+tdif*(qtt+tdif*qttt))+sdif*(qs+tdif*qst)		&
      +.5*tdif*(qpt+sdif*qpst+tdif*qptt)*(p+pref)*1.e-4)*(p-pref)*1.e-4		! 1.e-8: convert pres to dbar
!-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>
!   sigma_pref_star=sig9 + kappa
!-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>
    sigma_pref_star=sig9				! no thermobaricity
!-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>

  return
  end function sigma_pref_star
!
  real*8 function sigloc(t,s,prs)
! --- locally referenced sigma, a fit towards Jackett & McDougall (1995)
! --- t: potential temperature; s: psu; prs: pressure unit, converted to dbar
!
  use hycom_constants, only: alphap, betap, gammap
  implicit none
  real*8 t,s,prs,c1p,c2p,c3p,c4p,c5p,c6p,c7p,c8p,c9p
  c1p=alphap(1)+1.e-4*prs*(betap(1)+1.e-4*prs*gammap(1))
  c2p=alphap(2)+1.e-4*prs*(betap(2)+1.e-4*prs*gammap(2))
  c3p=alphap(3)+1.e-4*prs*(betap(3)+1.e-4*prs*gammap(3))
  c4p=alphap(4)+1.e-4*prs*(betap(4)+1.e-4*prs*gammap(4))
  c5p=alphap(5)+1.e-4*prs*(betap(5)+1.e-4*prs*gammap(5))
  c6p=alphap(6)+1.e-4*prs*(betap(6)+1.e-4*prs*gammap(6))
  c7p=alphap(7)+1.e-4*prs*(betap(7)+1.e-4*prs*gammap(7))
  c8p=alphap(8)+1.e-4*prs*(betap(8)+1.e-4*prs*gammap(8))
  c9p=alphap(9)+1.e-4*prs*(betap(9)+1.e-4*prs*gammap(9))
  sigloc=c1p+s*(c3p+c8p*s)+t*(c2p+c5p*s+t*(c4p+c7p*s+c6p*t))+c9p*t*s*s
  return
  end function sigloc
!
!
  real*8 function dsiglocdt(t,s,prs)
! --- locally referenced sigma, a fit towards Jackett & McDougall (1995)
! --- t: potential temperature; s: psu; prs: pressure
!
  use hycom_constants, only: alphap, betap, gammap
  implicit none
  real*8 :: t,s,prs,c2p,c4p,c5p,c6p,c7p,c9p
! c1p=alphap(1)+1.e-4*prs*(betap(1)+1.e-4*prs*gammap(1))
  c2p=alphap(2)+1.e-4*prs*(betap(2)+1.e-4*prs*gammap(2))
! c3p=alphap(3)+1.e-4*prs*(betap(3)+1.e-4*prs*gammap(3))
  c4p=alphap(4)+1.e-4*prs*(betap(4)+1.e-4*prs*gammap(4))
  c5p=alphap(5)+1.e-4*prs*(betap(5)+1.e-4*prs*gammap(5))
  c6p=alphap(6)+1.e-4*prs*(betap(6)+1.e-4*prs*gammap(6))
  c7p=alphap(7)+1.e-4*prs*(betap(7)+1.e-4*prs*gammap(7))
! c8p=alphap(8)+1.e-4*prs*(betap(8)+1.e-4*prs*gammap(8))
  c9p=alphap(9)+1.e-4*prs*(betap(9)+1.e-4*prs*gammap(9))
!
  dsiglocdt=c2p+s*(c5p+c9p*s)+2.*t*(c4p+c7p*s+1.5*c6p*t)    !  SI c9
  return
  end function dsiglocdt
!
  real*8 function dsiglocds(t,s,prs)
! --- locally referenced sigma, a fit towards Jackett & McDougall (1995)
! --- t: potential temperature; s: psu; prs: pressure
!
  use hycom_constants, only: alphap, betap, gammap
  implicit none
  real*8 :: t,s,prs,c3p,c5p,c7p,c8p,c9p
! c1p=alphap(1)+1.e-4*prs*(betap(1)+1.e-4*prs*gammap(1))
! c2p=alphap(2)+1.e-4*prs*(betap(2)+1.e-4*prs*gammap(2))
  c3p=alphap(3)+1.e-4*prs*(betap(3)+1.e-4*prs*gammap(3))
! c4p=alphap(4)+1.e-4*prs*(betap(4)+1.e-4*prs*gammap(4))
  c5p=alphap(5)+1.e-4*prs*(betap(5)+1.e-4*prs*gammap(5))
! c6p=alphap(6)+1.e-4*prs*(betap(6)+1.e-4*prs*gammap(6))
  c7p=alphap(7)+1.e-4*prs*(betap(7)+1.e-4*prs*gammap(7))
  c8p=alphap(8)+1.e-4*prs*(betap(8)+1.e-4*prs*gammap(8))
  c9p=alphap(9)+1.e-4*prs*(betap(9)+1.e-4*prs*gammap(9))
!
  dsiglocds=c3p+2*s*(c8p+c9p*t)+t*(c5p+t*c7p)                       
  return
  end function dsiglocds
!
!     subroutine totals(dp1,field1,dp2,field2,text)
!     implicit none
!
! --- compute volume integral of 2 fields (field1,field2), each associated
! --- with its own layer thickness field (dp1,dp2)
!
!     integer i,j,k,l
!
!     real dp1(idm,jdm,kdm),dp2(idm,jdm,kdm),field1(idm,jdm,kdm),
!    .     field2(idm,jdm,kdm),sum1j(jdm),sum2j(jdm),sum1,sum2
!     character text*(*)
!
!!$OMP PARALLEL DO SCHEDULE(STATIC,jchunk)
!     do 1 j=1,jj
!     sum1j(j)=0.
!     sum2j(j)=0.
!     do 1 k=1,kk
!     do 1 l=1,isp(j)
!     do 1 i=ifp(j,l),ilp(j,l)
!     sum1j(j)=sum1j(j)+dp1(i,j,k)*field1(i,j,k)*scp2(i,j)
!1    sum2j(j)=sum2j(j)+dp2(i,j,k)*field2(i,j,k)*scp2(i,j)
!!$OMP END PARALLEL DO
!
!     sum1=0.
!     sum2=0.
!!$OMP PARALLEL DO REDUCTION(+:sum1,sum2) SCHEDULE(STATIC,jchunk)
!     do 2 j=1,jj
!     sum1=sum1+sum1j(j)
!2    sum2=sum2+sum2j(j)
!
!     write (*,'(a,1p,2e19.9)') text,sum1,sum2
!     return
!     end
end module hycom_sigetc
