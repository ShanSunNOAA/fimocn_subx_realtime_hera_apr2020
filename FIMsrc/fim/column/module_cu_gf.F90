!WRF:MODEL_LAYER:PHYSICS
!

MODULE module_cu_gf
use physcons, g => con_g, cp => con_cp, xlv => con_hvap   &
     &,             r_v => con_rv
    real, parameter:: c1=.0002,c1_shal=0. ! .0005
    real, parameter:: c1_max=.002,c0_shal=.002
    real, parameter:: betajb=1.0

CONTAINS

!          call GFDRV(tottracer,ntrac,cactiv,ipn,lonf,IX,levs,DTP,forcet       &
!     &      ,forceq,phil,RAIN1
!     &      ,GQ0,GT0,cld1d,gu0,gv0,tGRS,VVEL,QGRS,prsl,pgr
!     &      ,kbot,ktop,SLMSK,hflx,evap,clw
!     &      ,hpbl,ud_mf,dd_mf,dt_mf,sdiaga,sdiagb,cnvw,cnvc,ishal_cnv)
!
! t2di is temp after advection, but before physics
! t = current temp (t2d+physics up to now)
!-------------------------------------------------------------
   SUBROUTINE GFDRV(tottracer,ntrac,cactiv,ipn,lonf,IX,KM,DT         &
              ,forcet,forceq,phil,RAINCV,q,t,cld1d                   &
              ,US,VS,t2di,W,q2di,p2di,psuri                          &
              ,hbot,htop,kcnv,XLAND,hfx,qfx,clw                           &
              ,pbl,ud_mf,dd_mf,dt_mf,gdc,gdc2,cnvw,cnvc,ishal_cnv)
!-------------------------------------------------------------
   IMPLICIT NONE
! autoconv, 1=old c0, 2=berry c0
      integer, parameter :: autoconv=1
!aeroevap, 1=old,2=?, 3=average
      integer, parameter :: aeroevap=1
      integer, parameter :: training=0
      integer, parameter :: trainingm=0
      integer, parameter :: use_excess=0
      integer, parameter :: use_excess_sh=0
      integer, parameter :: use_excess_m=0
      integer, parameter :: maxiens=1
      integer, parameter :: maxens=1
      integer, parameter :: maxens2=1
      integer, parameter :: maxens3=16
      integer, parameter :: ensdim=16
      integer            :: ishallow_g3	! 0 or 1 depending on ishal_cnv
      integer, parameter :: imid_gf=1
      integer, parameter :: ideep=1
      integer, parameter :: ichoice=0	! 0 2 5 8 11 GG
      integer, parameter :: ichoicem=5	! 0 2 5 8 11 GG
      integer, parameter :: ichoice_s=0	! 0 2 5 8 11 GG
      real, parameter :: ccnclean=250.
      real, parameter :: aodccn=0.1
      real :: beta,betam,dts,fpi,fp
      integer, parameter :: dicycle=2 !- diurnal cycle flag
!-------------------------------------------------------------
   INTEGER      :: its,ite, jts,jte, kts,kte 
   INTEGER,      INTENT(IN   ) :: ipn,IX,KM,lonf,ntrac,tottracer


   REAL,  DIMENSION( ix , km )         ,    &
          INTENT(IN   ) ::                                      &
                 forcet,forceq,w,phil
   REAL,  DIMENSION( ix , km )         ,    &
          INTENT(INOUT   ) ::                                      &
                                                          t,    &
                                                 us,vs
   REAL,  DIMENSION( ix , km, 11 )         ,    &
          INTENT(INOUT   ) ::                                   &
               GDC,GDC2
   REAL,  DIMENSION( ix , km )         ,    &
          INTENT(INOUT   ) ::                                   &
               cnvw,cnvc
   REAL,  DIMENSION( ix , km,tottracer+2 )         ,    &
          INTENT(INOUT   ) ::            clw

   REAL, DIMENSION( ix ),INTENT(IN) :: pbl,XLAND
!   integer, DIMENSION( 1 ),INTENT(INOUT) :: cactiv
   integer, INTENT(INOUT) :: cactiv(ix)
   INTEGER, INTENT(INOUT), DIMENSION( ix ) :: hbot,htop,kcnv
   INTEGER, DIMENSION( ix ) :: kpbl,k22_shallow, &
                 kbcon_shallow,ktop_shallow
!
   REAL, INTENT(IN   ) :: DT
   integer, intent(in) :: ishal_cnv
   REAL :: DX
!

   REAL, DIMENSION( ix ),                        &
         INTENT(INOUT) ::           RAINCV,cld1d
!+lxz

  REAL, DIMENSION( ix ) :: aoddd,ht

! LOCAL VARS
!    real,    dimension ( ix , 1:ensdim) ::      &
!       massflni,xfi_ens,pri_ens
     real,    dimension (ix,km) ::                    &
        OUTT,OUTQ,OUTQC,phh,subm,cupclw,cupclws,dhdt,         &
        zu,zus,zd,outts,outqs,phf,outqcs,outu,outv,outus,outvs,         &
        ud_mf,dd_mf,dt_mf,zum,zdm,outum,outvm
     real,    dimension (ix,km) ::                    &
        OUTTm,OUTQm,OUTQCm,submm,cupclwm,cnvwt,cnvwts,cnvwtm
     real,    dimension (ix,km) :: hco,hcdo,zdo,hcom,hcdom,zdom
     real,    dimension (km) ::                    &
        zh
     real,    dimension (ix)         ::                    &
        tau_ecmwf,edt,edtm,zws,ztexec,zqexec,pret, prets,pretm,ter11, aa0, xlandi
     real,    dimension (ix)         :: hexec
     real,    dimension (ix,10)         :: forcing
!+lxz
     integer, dimension (ix) ::                            &
        kbcon, ktop,ierr,ierrs,ierrm,kpbli,k22s,kbcons,ktops, &
        k22,jmin,jminm,kbconm,ktopm,k22m
!.lxz
     integer :: iens,ibeg,iend,jbeg,jend,n
     integer :: ibegh,iendh,jbegh,jendh
     integer :: ibegc,iendc,jbegc,jendc,kstop
   real rho_dryar,temp
   real :: PTEN,PQEN,PAPH,ZRHO,PAHFS,PQHFL,ZKHVFL,PGEOH


!
! basic environmental input includes moisture convergence (mconv)
! omega (omeg), windspeed (us,vs), and a flag (aaeq) to turn off
! convection for this call only and at that particular gridpoint
!
     real,    dimension (ix,km) ::                    &
        qcheck,zo,t2di,T2d,q2d,PO,P2d,p2di,rhoi,tn,qo,tshall,qshall
     real,    dimension (ix,km,ntrac) ::q2di,q
     real,    dimension (ix,km) ::                    &
        dz8w,omeg
     real, dimension (ix)            ::                    &
        ccn,Z1,psuri,PSUR,AAEQ,aaeqm,aaeqs,cuten,cutens,cutenm,umean, &
        vmean,pmean,xmbs,xmbs2,xmb,xmbm,xmb_dumm
     real, dimension (ix)     ::                    &
        mconv

   INTEGER :: i,j,k,ICLDCK,ipr,jpr,jpr_deep,ipr_deep
   REAL    :: excess,tcrit,tscl_KF,dp,dq,sub_spread,subcenter
   INTEGER :: itf,jtf,ktf,iss,jss,nbegin,nend
   INTEGER :: high_resolution
   REAL    :: dsubclw,dsubclws,dsubclwm,ztm,ztq,hfm,qfm,aodtest,rkbcon,rktop        !-lxz
   REAL, dimension (ix)    :: flux_tun,tun_rad_mid,tun_rad_shall,tun_rad_deep
   character*50 :: ierrc(ix),ierrcm(ix)
   character*50 :: ierrcs(ix)
! ruc variable
   REAL, dimension( ix  ) :: hfx,qfx
   real tem,tem1,tf, tcr, tcrf
   parameter (tf=233.16, tcr=263.16, tcrf=1.0/(tcr-tf))

   its=1
   ite=ix
   jts=1
   jte=1
   kts=1
   kte=km
   ktf=kte-1
!
! tuning constants for radiation coupling
!
   tun_rad_shall(:)=.0
   tun_rad_mid(:)=.0
   tun_rad_deep(:)=.0
   flux_tun(:)=1.5
  ! dx for scale awareness
   dx=40075000./float(lonf)
   tscl_kf=dx/25000.
   ccn(its:ite)=150.
  !
   ishallow_g3 = 0
   if (ishal_cnv == 2 .or. ishal_cnv == 3) ishallow_g3 = 1
   high_resolution=0
   subcenter=0.
   iens=1
!
! these can be set for debugging
!
   ipr=0
   jpr=0
   ipr_deep=0
   jpr_deep= 0 !53322 ! 528196 !0 ! 1136 !0 !421755 !3536
!
!
   ibeg=its
   iend=ite
   tcrit=258.

   itf=ite
   ktf=kte-1
   jtf=jte
   ztm=0.
   ztq=0.
   hfm=0.
   qfm=0.
   ud_mf =0.
   dd_mf =0.
   dt_mf =0.
   tau_ecmwf(:)=0.
!                                                                      
       j=1
       ht(:)=phil(:,1)/g
       do i=its,ite
        cld1d(i)=0.
        zo(i,:)=phil(i,:)/g
        dz8w(i,1)=zo(i,2)-zo(i,1)
        zh(1)=0.
        kpbli(i)=2
        do k=kts+1,ktf
          dz8w(i,k)=zo(i,k+1)-zo(i,k)
        enddo
        do k=kts+1,ktf
          zh(k)=zh(k-1)+dz8w(i,k-1)
          if(zh(k).gt.pbl(i))then
           kpbli(i)=max(2,k)
           exit
          endif
        enddo
       enddo
     DO I= its,itf
        forcing(i,:)=0.
        ccn(i)=100.
        HBOT(I)  =KTE
        HTOP(I)  =KTS
        raincv(i)=0.
        xlandi(i)=xland(i)
!       if(abs(xlandi(i)-1.).le.1.e-3) tun_rad_shall(i)=.15     
!       if(abs(xlandi(i)-1.).le.1.e-3) flux_tun(i)=1.5     
     ENDDO
!    DO I= its,itf
!       aodtest=min(aoddd(i),1.)
!       ccn(i) = max( 50., ( 370.37*aoddd(i) )**1.555 )
!    enddo
     DO I= its,itf
        mconv(i)=0.
     ENDDO
     do k=kts,ktf
     DO I= its,itf
         omeg(i,k)=0.
         zu(i,k)=0.
         zum(i,k)=0.
         zus(i,k)=0.
         zd(i,k)=0.
         zdm(i,k)=0.
     ENDDO
     ENDDO

     psur(:)=0.01*psuri(:)
     DO I=ITS,ITF
         TER11(I)=max(0.,HT(i))
!        DO K=kts,ktf
!        zo(i,k)=phil(i,k)/g
!        enddo
     ENDDO
!   if(ipn.eq.jpr_deep )write(12,*)'j,psu,ter11,kpbl,zo(1) = ',j,psur(1),ter11(1),kpbli(1),zo(1,kts)
     DO K=kts,kte
     DO I=ITS,ITe
         cnvw(i,k)=0.
         cnvc(i,k)=0.
         gdc(i,k,1)=0.
         gdc(i,k,2)=0.
         gdc(i,k,3)=0.
         gdc(i,k,4)=0.
         gdc(i,k,7)=0.
         gdc(i,k,8)=0.
         gdc(i,k,9)=0.
         gdc(i,k,10)=0.
         gdc2(i,k,1)=0.
     ENDDO
     ENDDO
     ierr(:)=0
     ierrm(:)=0
     ierrs(:)=0
     cuten(:)=0.
     cutenm(:)=0.
     cutens(:)=1.
     if(ishallow_g3.eq.0)cutens(:)=0.
     ierrc(:)=" "

     kbcon(:)=0
     kbcons(:)=0
     kbconm(:)=0

     ktop(:)=0
     ktops(:)=0
     ktopm(:)=0

     xmb(:)=0.
     xmb_dumm(:)=0.
     xmbm(:)=0.
     xmbs(:)=0.
     xmbs2(:)=0.

     k22s(:)=0
     k22m(:)=0
     k22(:)=0

     jmin(:)=0
     jminm(:)=0

     ZTEXEC(:) = 0.
     ZQEXEC(:) = 0.

     aaeq(:)=0.
     aaeqm(:)=0.
     aaeqs(:)=0.

     pret(:)=0.
     prets(:)=0.
     pretm(:)=0.

     umean(:)=0.
     vmean(:)=0.
     pmean(:)=0.

     cupclw(:,:)=0.
     cupclwm(:,:)=0.
     cupclws(:,:)=0.

     cnvwt(:,:)=0.
     cnvwts(:,:)=0.

     hco(:,:)=0.
     hcom(:,:)=0.
     hcdo(:,:)=0.
     hcdom(:,:)=0.

     OUTT(:,:)=0.
     OUTTS(:,:)=0.
     OUTTm(:,:)=0.

     OUTU(:,:)=0.
     OUTUS(:,:)=0.
     OUTUM(:,:)=0.

     OUTV(:,:)=0.
     OUTVS(:,:)=0.
     OUTVM(:,:)=0.

     OUTQ(:,:)=0.
     OUTQS(:,:)=0.
     outqm(:,:)=0.

     OUTQC(:,:)=0.
     OUTQCS(:,:)=0.
     OUTQCm(:,:)=0.

     subm(:,:)=0.
     dhdt(:,:)=0.
     DO K=kts,ktf
     DO I=ITS,ITF
         p2d(i,k)=0.01*p2di(i,k)
         po(i,k)=p2d(i,k) !*.01
         rhoi(i,k) = 100.*p2d(i,k)/(287.04*(t2di(i,k)*(1.+0.608*q2di(i,k,1))))
         qcheck(I,K)=q(i,k,1)
         Tshall(I,K)=t(i,k)
         qshall(I,K)=q(i,k,1)
!        tn(i,k)=t(i,k)+forcet(i,k)*dt
!        qo(i,k)=q(i,k)+forceq(i,k)*dt
         tn(i,k)=t(i,k)!+forcet(i,k)*dt
         qo(i,k)=max(1.e-16,q(i,k,1))!+forceq(i,k)*dt
         t2d(i,k)=t2di(i,k)-forcet(i,k)*dt
         q2d(i,k)=max(1.e-16,q2di(i,k,1)-forceq(i,k)*dt)
         IF(Qo(I,K).LT.1.E-16)Qo(I,K)=1.E-16
!        if(ipn.eq.jpr_deep)then
!         write(12,123)k,p2d(i,k),t2d(i,k),tn(i,k),q2d(i,k),QO(i,k),rhoi(i,k)
!        endif
     ENDDO
     ENDDO
123  format(1x,i2,1x,f8.0,1x,2(1x,f8.3),3(1x,e13.5))
!    if(use_excess.ge.1 .or. use_excess_sh.ge.1 .or. use_excess_m.ge.1 )then
     DO I=ITS,ITF
       PTEN = t2d(i,1)
       PQEN = q2d(I,1)
       PAPH = 100.*psur(i)
       ZRHO  = PAPH/(287.04*(t2d(i,1)*(1.+0.608*q2d(i,1))))
       !- LE and H fluxes 
       ! must be < 0 for daytime over the land and W/m2
       PAHFS=-hfx(i) ! hfx in w/m2
       ! must be < 0 for daytime over the land and kg/m2/s
       PQHFL=-qfx(i)/xlv ! qfx in w/m2
       !- buoyancy flux (H+LE)
       ZKHVFL= (PAHFS/1004.64+0.608*PTEN*PQHFL)/ZRHO
       !- height of the 1st level
       PGEOH = dz8w(i,1)*g ! zo(i,1)-ht(i) !zt(2)*rtgt(i)
       !-convective-scale velocity w*
       ZWS(i) = max(0.,0.001-1.5*0.41*ZKHVFL*PGEOH/PTEN)
!      if(j.eq.jpr_deep)write(0,*)'excess1=' ,j,hfx(i),qfx(i),pgeoh,zws(i)
!      if(j.eq.jpr_deep)write(0,*)'excess1=' ,j,hfx(i),qfx(i),pgeoh,zws(i)

       if(ZWS(i) > TINY(PGEOH)) then
         !-convective-scale velocity w*
         ZWS(i) = 1.2*ZWS(i)**.3333
!        !- temperature excess 
!        ZTEXEC(i)     = MAX(0.,-flux_tun(i)*PAHFS/(ZRHO*ZWS(i)*1004.64))
         ZTEXEC(i)     = -flux_tun(i)*PAHFS/(ZRHO*ZWS(i)*1004.64)
         !- moisture  excess
!        ZQEXEC(i)     = MAX(0.,-flux_tun(i)*PQHFL/(ZRHO*ZWS(i)))
         ZQEXEC(i)     = -flux_tun(i)*PQHFL/(ZRHO*ZWS(i))
         hexec(i)=ztexec(i)*cp +zqexec(i)*2.5e6
!      if(j.eq.jpr_deep)write(0,*)'excess2=' ,i,ZTEXeC(i)  ,ZQEXeC(i)
       !srf  ZWS(I)=zws(i)*zrho
        
       endif   ! zws > 0
       !- ZWS for shallow convection closure (Grant 2001)
       !- height of the pbl
       PGEOH = pbl(i)*g
       !-convective-scale velocity w*
!       ZWS(i) = max(0.,0.001-flux_tun(i)*0.41*ZKHVFL*PGEOH/PTEN)
       ZWS(i) = max(0.,0.001-1.5*0.41*ZKHVFL*PGEOH/PTEN)
       ZWS(i) = 1.2*ZWS(i)**.3333
       ZWS(i) = ZWS(i)*zrho !check if zrho is correct

!      if(ZTEXeC(i).gt.ztm)ztm=ZTEXeC(i)
!      if(ZqEXeC(i).gt.ztq)ztq=ZqEXeC(i)
!      if(-hfx(i).gt.hfm)hfm=-hfx(i)
!      if(-qfx(i).gt.qfm)qfm=-qfx(i)


!srf     if(xlandi(i).lt.0.5 .or. xlandi(i).gt.1.5)then
!    if(use_excess.ge.1 .or. use_excess_sh.ge.1 .or. use_excess_m.ge.1 )then
!        ztexec(i)=0.
!        zqexec(i)=0.
!srf     endif
     ENDDO
!    endif  ! use_excess > 0
     DO I=ITS,ITF
     DO K=kts,kpbli(i)
         tn(i,k)=t(i,k)+float(use_excess)*ztexec(i) !+forcet(i,k)*dt
         qo(i,k)=max(1.e-16,q(i,k,1)+float(use_excess)*zqexec(i)) !+forceq(i,k)*dt
     ENDDO
     ENDDO
     nbegin=0
     nend=0
         DO I=ITS,ITF
         DO K=kts,kpbli(i)
         DHDT(I,K)=cp*(forcet(i,k)+(t(i,k)-t2di(i,k))/dt+float(use_excess_sh)*ztexec(i)/dt)+ &
                   XLV*(forceq(i,k)+(q(i,k,1)-q2di(i,k,1))/dt+float(use_excess_sh)*zqexec(i)/dt)
         TSHALL(I,K)=t(i,k)+forcet(i,k)*dt+float(use_excess_sh)*ztexec(i)
         QSHALL(I,K)=q(i,k,1)+forceq(i,k)*dt+float(use_excess_sh)*zqexec(i)
        enddo
        enddo
      do k=  kts+1,ktf-1
      DO I = its,itf
         if((p2d(i,1)-p2d(i,k)).gt.150.and.p2d(i,k).gt.300)then
            dp=-.5*(p2d(i,k+1)-p2d(i,k-1))
            umean(i)=umean(i)+us(i,k)*dp
            vmean(i)=vmean(i)+vs(i,k)*dp
            pmean(i)=pmean(i)+dp
         endif
      enddo
      enddo
      DO K=kts,ktf-1
      DO I = its,itf
        omeg(I,K)= w(i,k) ! -g*rhoi(i,k)*w(i,k)
      ENDDO
      ENDDO
!
!---- CALL CUMULUS PARAMETERIZATION
!
       if(ishallow_g3.eq.1)then
!
! this turns off shallow convection when deep convection is active
!
       do i=its,ite
        aaeq(i)=0.
        aaeqm(i)=0.
        aaeqs(i)=0.
        ierrs(i)=0
        ierrm(i)=0
!       if(pret(i).gt.0.)then
!           ierrs(i)=1
!           aaeq(i)=-100.
!       endif  ! pret > 0
       enddo
!
! call shallow
!
   call CUP_gf_sh(cnvwts,zus,zws,xmbs,xmbs2,zo,OUTQCs,ipn,AAEQ,T2D,Q2D,TER11, &
              prets,Tshall,Qshall,P2d,P2d,OUTTS,OUTQS,DT,PSUR,US,VS,          &
              TCRIT,ztexec,zqexec,ccn,ccnclean,rhoi,dx,dhdt,                  &
              cupclws,kpbli,kbcons,ktops,k22s,xlandi,                         &
              tscl_kf,ichoice_s,ipr,jpr,ierrs,ierrcs,                         &
              autoconv,itf,jtf,ktf,                                           &
              use_excess_sh,its,ite, jts,jte, kts,kte                         &
                                                                              )
       do i=its,itf
        if(xmbs(i).le.0.)cutens(i)=0.
       enddo
       CALL neg_check('shallow',ipn,dt,qcheck,outqs,outts,outus,outvs,   &
                                 outqcs,prets,its,ite,kts,kte,itf,ktf)
       endif

!!sms$ignore begin
!      if(ipn.eq.jpr_deep)then
!        i=1
!        write(12,129)ipn,ter11(i),dx,psur(i),kpbli(i),dt,ztexec(i),zqexec(i),zws(i)
!      do k=kts,ktf
!        write(12,128)k,t2d(i,k),tshall(i,k),q2d(i,k),qshall(i,k),p2d(i,k),zo(i,k),us(i,k),vs(i,k),omeg(i,k)
!      end do
!      endif
!!sms$ignore end
!      jpr_deep=0
!      ipr_deep=0
!        if(ipn.eq.jpr_deep)then
!!sms$ignore begin
!         i=1
!         j=1
!         write(13,129)ipn,ter11(i),dx,psur(i),kpbli(i),dt,ztexec(i),zqexec(i),mconv(i)
!       do k=kts,ktf
!         write(13,127)k,t2d(i,k),tn(i,k),us(i,k),vs(i,k),q2d(i,k),qo(i,k),p2d(i,k),zo(i,k),rhoi(i,k),omeg(i,k)
!       end do
!        jpr=1
!        ipr=1
!!sms$ignore end
!        endif
128   format(1x,i2,2f7.2,4e13.6,1x,f8.1,1x,f8.1,1x,e11.3)
127   format(1x,i2,4f7.2,2e13.6,1x,2(1x,f8.1),1x,2e13.6)
129   format(1x,i8,1x,e13.6,1x,f7.0,1x,f8.1,1x,i2,1x,f6.0,1x,f6.2,1x,2e13.6)
! memory to GF!!
      if(imid_gf.eq.1)then
!
! - change tshall, qshall, and dhdt with just calculated shallow tendencies
!
!            DO I=its,itf
!            DO K=kts,ktf
!              qshall(i,k)=max(1.e-10,qshall(i,k) +outqs(i,k)*dt)
!              tshall(i,k)=tshall(i,k) +outts(i,k)*dt
!            enddo
!            enddo
!            DO I=ITS,ITF
!            DO K=kts,kpbli(i)
!               DHDT(I,K)=dhdt(i,k)+cp*outts(i,k)+XLV*outqs(i,k)
!            enddo
!            enddo
      CALL CUP_gf(cactiv,cnvwtm,zum,zdm,edtm,kpbli,zws,dhdt, &
           imid_gf,xmbm,xmb_dumm,xmbs,zo,forcing,outqcm,  &
           ipn,AAEQm,T2d,Q2d,TER11,subm,Tshall,Qshall,PO,PRETm,&
           outum,outvm,P2d,OUTTm,outqm,DT,PSUR,US,VS,tcrit,iens, &
           ztexec,zqexec,ccn,ccnclean,rhoi,dx,mconv,omeg,  &
           maxiens,maxens,maxens2,maxens3,ensdim,                 &
           kbconm,ktopm,cupclwm, XLANDi,        &
           ichoicem,0,jpr_deep,ierrm,ierrcm,    &
           betam,autoconv,aeroevap,itf,jtf,ktf,training, &
           use_excess_m,its,ite, jts,jte, kts,kte , &
           tau_ecmwf,2,hcom,hcdom,k22m,jminm  )
!       if(ipn.eq.jpr_deep)then
!!sms$ignore begin
!      do i=its,itf
!      write(13,*)'k22,jmin,kbcon,ktop = ',k22m(i),jminm(i),kbconm(i),ktopm(i)
!       write(13,*)'xmb,ztex = ',xmbm(i),edtm(i) !xfi_ens(i,1:ensdim)
!       if(ktopm(i).gt.2)then
!       do k=1,ktopm(i)+1
!         write(13,333)k,zum(i,k),zdm(i,k),hcom(i,k),hcdom(i,k), &
!                 (outtm(i,k))*86400., &
!                 xlv/cp*86400.*(outqm(i,k)), &
!                   outqcm(i,k)*xlv/cp*86400.
!       enddo
!      endif
!      enddo
 333    format(1x,i4,1x,2(1x,F5.2),2e13.4,3(1x,f9.2))
!!sms$ignore end
!      endif
            DO I=its,itf
            DO K=kts,ktf
              qcheck(i,k)=q(i,k,1) +outqs(i,k)*dt
            enddo
            enddo
      CALL neg_check('mid',ipn,dt,qcheck,outqm,outtm,outum,outvm,   &
                     outqcm,pretm,its,ite,kts,kte,itf,ktf)
!!sms$ignore begin
!           do i=its,itf
!!         if(pretm(i).gt.0. .and. ipn.lt.100000)write(14,*)ipn,pretm(i),xmbm(i)
!         write(14,*)ipn,pretm(i),xmbm(i)
!         enddo
!        call flush(14)
!!sms$ignore end
      endif
      if(ideep.eq.1)then
      CALL CUP_gf(cactiv,cnvwt,zu,zd,edt,kpbli,zws,dhdt,0,xmb,xmbm,xmbs,zo,   &
           forcing,outqc,ipn,AAEQ,T2d,Q2d,TER11,subm,TN,QO,PO,PRET,&
           outu,outv,P2d,OUTT,OUTQ,DT,PSUR,US,VS,tcrit,iens, &
           ztexec,zqexec,ccn,ccnclean,rhoi,dx,mconv,omeg,  &
           maxiens,maxens,maxens2,maxens3,ensdim,                 &
           kbcon,ktop,cupclw,XLANDi,        &
           ichoice,ipr,jpr,ierr,ierrc,    &
           beta,autoconv,aeroevap,itf,jtf,ktf,training, &
           use_excess,its,ite, jts,jte, kts,kte , &
           tau_ecmwf,dicycle,hco,hcdo,k22,jmin  )
        jpr=0
        ipr=0
!       if(ipn.eq.jpr_deep)then
!!sms$ignore begin
!      do i=its,itf
!      write(13,*)'k22,jmin,kbcon,ktop = ',k22(i),jmin(i),kbcon(i),ktop(i)
!       write(13,*)'xmb,ztex = ',xmb(i),hexec(i) !xfi_ens(i,1:ensdim)
!       do k=1,ktop(i)+1
!         write(13,333)k,zu(i,k),zd(i,k),hco(i,k),hcdo(i,k), &
!                 (outt(i,k))*86400., &
!                 xlv/cp*86400.*(outq(i,k)), &
!                   outqc(i,k)*xlv/cp*86400.
!       enddo
!      enddo
! 333    format(1x,i4,1x,2(1x,F5.2),2e13.4,3(1x,f9.2))
!!sms$ignore end
!      endif

            DO I=its,itf
            DO K=kts,ktf
              qcheck(i,k)=q(i,k,1) +(outqs(i,k)+outqm(i,k))*dt
            enddo
            enddo
      CALL neg_check('deep',ipn,dt,qcheck,outq,outt,outu,outv,   &
                                         outqc,pret,its,ite,kts,kte,itf,ktf)
!
      endif
! since kinetic energy is being dissipated, add heating accordingly (from ECMWF)
!
      do i=its,itf
          if(pret(i).gt.0) then
             dts=0.
             fpi=0.
!!sms$ignore begin
!             write(13,*)'om output for i,j,ktop,pret = ',i,j,ktop(i),pret(i)
             do k=kts+1,ktop(i)
                dp=-.5*(p2d(i,k+1)-p2d(i,k-1))*100.
!total KE dissiptaion estimate
                dts= dts +(outu(i,k)*us(i,k)+outv(i,k)*vs(i,k))*dp/g
! fpi needed for calcualtion of conversion to pot. energyintegrated 
                fpi = fpi  +sqrt(outu(i,k)*outu(i,k) + outv(i,k)*outv(i,k))*dp
             enddo
             if(fpi.gt.0.)then
             do k=kts+1,ktop(i)
                fp= sqrt((outu(i,k)*outu(i,k)+outv(i,k)*outv(i,k)))/fpi
                outt(i,k)=outt(i,k)+fp*dts*g/cp
!                write(13,*)'k,fpi,dts,outt,fp =',k,fpi,dts,outt(i,k),fp
             enddo
!!sms$ignore end
             endif
          endif
      enddo




!!sms$ignore begin
!!               if(cactiv(1).gt.5)write(12,*)'ipn,pret(1),cactiv(1)=',ipn,pret(1),cactiv(1)
!!sms$ignore end
!           if(j.lt.jbegc.or.j.gt.jendc)go to 100
            DO I=its,itf
              if(pret(i).gt.0.)then
                 cuten(i)=1.
                 cld1d(i)=aaeq(i)
              else 
                 kbcon(i)=0
                 ktop(i)=0
                 cuten(i)=0.
              endif   ! pret > 0
              if(pretm(i).gt.0.)then
                 cutenm(i)=1.
!                cld1d(i)=aaeq(i)
              else 
                 kbconm(i)=0
                 ktopm(i)=0
                 cutenm(i)=0.
              endif   ! pret > 0
            ENDDO
!!sms$ignore begin
            kstop=kts
            DO I=its,itf
            if(ktopm(i).gt.kts .or. ktop(i).gt.kts)kstop=max(ktopm(i),ktop(i))
            if(ktops(i).gt.kts)kstop=max(kstop,ktops(i))
            if(kstop.gt.2)then
            htop(i)=kstop
            kcnv(i)=1
            if(kbcon(i).gt.2)hbot(i)=jmin(i) !kbcon(i)
!            if(kbconm(i).gt.2)hbot(i)=min(kbconm(i),hbot(i))
            DO K=kts,kstop
               cnvc(i,k) = 0.04 * log(1. + 675. * zu(i,k) * xmb(i)) +   &
                           0.04 * log(1. + 675. * zum(i,k) * xmbm(i)) + &
                           0.04 * log(1. + 675. * zus(i,k) * xmbs(i))
               cnvc(i,k) = min(cnvc(i,k), 0.6)
               cnvc(i,k) = max(cnvc(i,k), 0.0)
!               cnvw(i,k)=cuten(i)*cnvwt(i,k)*xmb(i)*dt+cutens(i)*cnvwts(i,k)*xmbs(i)*dt
               cnvw(i,k)=cnvwt(i,k)*xmb(i)*dt+cnvwts(i,k)*xmbs(i)*dt+cnvwtm(i,k)*xmbm(i)*dt
!               cnvw(i,k) = max(cnvw(i,k), 0.0)
               ud_mf(i,k)=cuten(i)*zu(i,k)*xmb(i)*dt
               dd_mf(i,k)=cuten(i)*zd(i,k)*edt(i)*xmb(i)*dt
               T(I,K)=t(i,k)+dt*(cutens(i)*outts(i,k)+cutenm(i)*outtm(i,k)+outt(i,k)*cuten(i))
               Q(I,K,1)=max(1.e-16,q(i,k,1)+dt*(cutens(i)*outqs(i,k)+cutenm(i)*outqm(i,k)+outq(i,k)*cuten(i)))
               gdc(i,k,7)=sqrt(us(i,k)**2 +vs(i,k)**2)
               us(i,k)=us(i,k)+outu(i,k)*cuten(i)*dt +outum(i,k)*cutenm(i)*dt
               vs(i,k)=vs(i,k)+outv(i,k)*cuten(i)*dt +outvm(i,k)*cutenm(i)*dt

!               T(I,K)=t(i,k)+dt*(cutens(i)*outts(i,k)+(outt(i,k))*cuten(i))
!               Q(I,K,1)=q(i,k,1)+dt*(cutens(i)*outqs(i,k)+(outq(i,k))*cuten(i))
!               us(i,k)=us(i,k)+outu(i,k)*cuten(i)*dt
!               vs(i,k)=vs(i,k)+outv(i,k)*cuten(i)*dt
               gdc(i,k,1)= max(0.,tun_rad_shall(i)*cupclws(i,k)*cutens(i))	! my mod
               gdc2(i,k,1)=max(0.,tun_rad_deep(i)*(cupclwm(i,k)*cutenm(i)+cupclw(i,k)*cuten(i)))
                gdc(i,k,2)=(outt(i,k))*86400.
                gdc(i,k,3)=(outtm(i,k))*86400. 
                gdc(i,k,4)=(outts(i,k))*86400.
               gdc(i,k,7)=-(gdc(i,k,7)-sqrt(us(i,k)**2 +vs(i,k)**2))/dt
                gdc(i,k,8)=(outq(i,k))*86400.*xlv/cp
                gdc(i,k,10)=(outqs(i,k))*86400.*xlv/cp

!
! subsidence effect on clw
!
               dsubclw=0.
               dsubclwm=0.
               dsubclws=0.
               dp=100.*(p2d(i,k)-p2d(i,k+1))
               if (clw(i,k,2) .gt. -999.0 .and. clw(i,k+1,2) .gt. -999.0 )then
                  dsubclw=((-edt(i)*zd(i,k+1)+zu(i,k+1))*clw(i,k+1,1)   &
                       -(-edt(i)*zd(i,k)  +zu(i,k))  *clw(i,k,1)  )*g/dp
                  dsubclwm=((-edtm(i)*zdm(i,k+1)+zum(i,k+1))*clw(i,k+1,1)   &
                       -(-edtm(i)*zdm(i,k)  +zum(i,k))  *clw(i,k,1)  )*g/dp
                  dsubclws=(zus(i,k+1)*clw(i,k+1,1)-zus(i,k)*clw(i,k,1))*g/dp
                  dsubclw=dsubclw+(zu(i,k+1)*clw(i,k+1,2)-zu(i,k)*clw(i,k,2))*g/dp 
                  dsubclwm=dsubclwm+(zum(i,k+1)*clw(i,k+1,2)-zum(i,k)*clw(i,k,2))*g/dp 
                  dsubclws=dsubclws+(zus(i,k+1)*clw(i,k+1,2)-zus(i,k)*clw(i,k,2))*g/dp 
               endif
!               if(k.eq.ktop(i))dsubclw=0
!               if(k.eq.ktopm(i))dsubclwm=0
!               if(k.eq.ktops(i))dsubclws=0
!
!
!               tem  = outqcs(i,k)*cutens(i)*dt+outqc(i,k)*dt*cuten(i) &
!                      +dsubclw*xmb(i)*dt+dsubclws*xmbs(i)*dt
               tem  = dt*(outqcs(i,k)*cutens(i)+outqc(i,k)*cuten(i)       &
                      +outqcm(i,k)*cutenm(i)                           &
!                      +dsubclw*xmb(i)+dsubclws*xmbs(i)+dsubclwm*xmbm(i) &
                      )
                gdc(i,k,9)=dsubclw*xmb(i)*86400.+dsubclws*xmbs(i)*86400. &
                           +dsubclwm*xmbm(i)*86400.
               tem1 = max(0.0, min(1.0, (tcr-t(i,k))*tcrf))
               if (clw(i,k,2) .gt. -999.0) then
                clw(i,k,1) = max(0.,clw(i,k,1) + tem * tem1)            ! ice
                clw(i,k,2) = max(0.,clw(i,k,2) + tem *(1.0-tem1))       ! water
              else
                clw(i,k,1) = max(0.,clw(i,k,1) + tem)
              endif

!               if(i.eq.ipr_deep.and.j.eq.jpr_deep)then
!                 write(0,*)'output',j,zo(i,k),86400.*(outt(i,k)),              &
!                  86400.*xlv/cp*(outq(i,k))
!               endif
            ENDDO
!        call flush(14)
!              gdc(i,ktf-2)=forcing(i,9)
!              gdc(i,ktf-3)=forcing(i,3)
!              gdc(i,ktf-4)=forcing(i,4)
!              gdc(i,ktf-5)=forcing(i,5)
!               gdc2(i,ktf)=zws(i)
!               gdc2(i,ktf-1)=xmbm(i)
!               gdc2(i,ktf-2)=xmb(i)
!               gdc2(i,ktf-3)=pretm(i)
!               gdc2(i,ktf-4)=hfx(i)
!               gdc2(i,ktf-5)=ztexec(i)
            if(ktop(i).gt.2 .and.pret(i).gt.0.)dt_mf(i,ktop(i)-1)=ud_mf(i,ktop(i))
            endif
            ENDDO
!!sms$ignore end
            DO I=its,itf
              if(pret(i).gt.0.)then
                 cactiv(i)=1
                 raincv(i)=.001*(cutenm(i)*pretm(i)+cutens(i)*prets(i)+cuten(i)*pret(i))*dt
!                 if(raincv(i).gt.100.)then
!                   write(0,*)'rc,ipn = ',raincv(i),ipn,cuten(i),cutens(i),pret(i),xmb(i)
!                 endif
!                 rkbcon = kte+kts - kbcon(i)
!                 rktop  = kte+kts -  ktop(i)
!                 if (ktop(i)  > HTOP(i)) HTOP(i) = ktop(i)
!                 if (kbcon(i) < HBOT(i)) HBOT(i) = kbcon(i)
              else
                 cactiv(i)=0
                 if(pretm(i).gt.0)raincv(i)=.001*cutenm(i)*pretm(i)*dt
              endif   ! pret > 0
            ENDDO
 100    continue

!       write(0,*)'ztm,ztq = ',ztm,ztq,hfm,qfm

   END SUBROUTINE GFDRV


   SUBROUTINE CUP_gf(csum,cnvwt,zuo,zdo,edto,kpbl,zws,dhdt,imid,      &
              xmb_out,xmbm_in,xmbs_in,zo,forcing,OUTQC,J,AAEQ,T,Q,Z1,sub_mas,   &
              TN,QO,PO,PRE,outu,outv,P,OUTT,OUTQ,DTIME,PSUR,US,VS,    &
              TCRIT,iens,                                        &
              ztexec,zqexec,ccn,ccnclean,rho,dx,mconv,                               &
              omeg,maxiens,                          &
              maxens,maxens2,maxens3,ensdim,                           &
              kbcon,ktop,cupclw,         &   !-lxz
              xland,               &
              ichoice,ipr,jpr,ierr,ierrc,         &
              beta,autoconv,aeroevap,itf,jtf,ktf,training,               &
#if ( WRF_DFI_RADAR == 1 )
              do_capsuppress,cap_suppress_j,                         &             
#endif
              use_excess,its,ite, jts,jte, kts,kte                                &
              ,tau_ecmwf_out,dicycle,hco,hcdo,k22,jmin                                   )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        autoconv,aeroevap,itf,jtf,ktf,training,use_excess,        &
        its,ite, jts,jte, kts,kte,ipr,jpr,imid
     integer, intent (in   )              ::                           &
        j,ensdim,maxiens,maxens,maxens2,maxens3,ichoice,iens
  !
  ! 
  !
      real,    dimension (its:ite,1:ensdim) :: xf_ens,pr_ens
    real, dimension( its:ite  )                               &
          :: weight_GR,weight_W,weight_MC,weight_ST,weight_AS
  ! outtem = output temp tendency (per s)
  ! outq   = output q tendency (per s)
  ! outqc  = output qc tendency (per s)
  ! pre    = output precip
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout  )                   ::                           &
        cnvwt,outu,outv,OUTT,OUTQ,OUTQC,sub_mas,cupclw
     real,    dimension (its:ite)                                      &
        ,intent (inout  )                   ::                           &
        pre,xmb_out
     real,    dimension (its:ite)                                      &
        ,intent (in  )                   ::                           &
        xmbm_in,xmbs_in
     integer,    dimension (its:ite)                                   &
        ,intent (inout  )                   ::                           &
        kbcon,ktop
     integer,    dimension (its:ite)                                   &
        ,intent (in  )                   ::                           &
        kpbl
  !
  ! basic environmental input includes moisture convergence (mconv)
  ! omega (omeg), windspeed (us,vs), and a flag (aaeq) to turn off
  ! convection for this call only and at that particular gridpoint
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        dhdt,rho,T,PO,P,US,VS,tn
     real,    dimension (its:ite,kts:kte)                       &
        ,intent (inout   )                   ::                           &
        omeg
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout)                   ::                           &
         Q,QO
     real, dimension (its:ite)                                         &
        ,intent (in   )                   ::                           &
        ccn,Z1,PSUR,xland
     real, dimension (its:ite)                                         &
        ,intent (inout   )                   ::                           &
        zws,ztexec,zqexec,aaeq
     real, dimension (its:ite)                                         &
        ,intent (inout   )                   ::                           &
        mconv

       
       real                                                            &
        ,intent (in   )                   ::                           &
        dx,ccnclean,dtime,tcrit
       real                                                            &
        ,intent (inout   )                   ::                           &
        beta


!
!  local ensemble dependent variables in this routine
!
     real,    dimension (its:ite,1:maxens)  ::                         &
        xaa0_ens
     real,    dimension (1:maxens)  ::                                 &
        mbdt_ens
     real,    dimension (1:maxens2) ::                                 &
        edt_ens
     real,    dimension (its:ite,1:maxens2) ::                         &
        edtc
     real,    dimension (its:ite,kts:kte,1:maxens2) ::                 &
        dellat_ens,dellaqc_ens,dellaq_ens,pwo_ens
!
!
!
!***************** the following are your basic environmental
!                  variables. They carry a "_cup" if they are
!                  on model cloud levels (staggered). They carry
!                  an "o"-ending (z becomes zo), if they are the forced
!                  variables. They are preceded by x (z becomes xz)
!                  to indicate modification by some typ of cloud
!
  ! z           = heights of model levels
  ! q           = environmental mixing ratio
  ! qes         = environmental saturation mixing ratio
  ! t           = environmental temp
  ! p           = environmental pressure
  ! he          = environmental moist static energy
  ! hes         = environmental saturation moist static energy
  ! z_cup       = heights of model cloud levels
  ! q_cup       = environmental q on model cloud levels
  ! qes_cup     = saturation q on model cloud levels
  ! t_cup       = temperature (Kelvin) on model cloud levels
  ! p_cup       = environmental pressure
  ! he_cup = moist static energy on model cloud levels
  ! hes_cup = saturation moist static energy on model cloud levels
  ! gamma_cup = gamma on model cloud levels
!
!
  ! hcd = moist static energy in downdraft
  ! zd normalized downdraft mass flux
  ! dby = buoancy term
  ! entr = entrainment rate
  ! zd   = downdraft normalized mass flux
  ! entr= entrainment rate
  ! hcd = h in model cloud
  ! bu = buoancy term
  ! zd = normalized downdraft mass flux
  ! gamma_cup = gamma on model cloud levels
  ! qcd = cloud q (including liquid water) after entrainment
  ! qrch = saturation q in cloud
  ! pwd = evaporate at that level
  ! pwev = total normalized integrated evaoprate (I2)
  ! entr= entrainment rate
  ! z1 = terrain elevation
  ! entr = downdraft entrainment rate
  ! jmin = downdraft originating level
  ! kdet = level above ground where downdraft start detraining
  ! psur        = surface pressure
  ! z1          = terrain elevation
  ! pr_ens = precipitation ensemble
  ! xf_ens = mass flux ensembles
  ! massfln = downdraft mass flux ensembles used in next timestep
  ! omeg = omega from large scale model
  ! mconv = moisture convergence from large scale model
  ! zd      = downdraft normalized mass flux
  ! zu      = updraft normalized mass flux
  ! dir     = "storm motion"
  ! mbdt    = arbitrary numerical parameter
  ! dtime   = dt over which forcing is applied
  ! iact_gr_old = flag to tell where convection was active
  ! kbcon       = LFC of parcel from k22
  ! k22         = updraft originating level
  ! icoic       = flag if only want one closure (usually set to zero!)
  ! dby = buoancy term
  ! ktop = cloud top (output)
  ! xmb    = total base mass flux
  ! hc = cloud moist static energy
  ! hkb = moist static energy at originating level

     real,    dimension (its:ite,kts:kte) ::                           &
        entr_rate_2d,mentrd_rate_2d,he,hes,qes,z,                      &
        heo,heso,qeso,zo,                                              &
        xhe,xhes,xqes,xz,xt,xq,                                        &

        qes_cup,q_cup,he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup,      &
        qeso_cup,qo_cup,heo_cup,heso_cup,zo_cup,po_cup,gammao_cup,     &
        tn_cup,                                                        &
        xqes_cup,xq_cup,xhe_cup,xhes_cup,xz_cup,xp_cup,xgamma_cup,     &
        xt_cup,hcot,                                             &

        dby,qc,qrcd,pwd,pw,hcd,qcd,dbyd,hc,qrc,zu,zd,clw_all,   &
        dbyo,qco,qrcdo,pwdo,pwo,hcdo,qcdo,dbydo,hco,qrco,zuo,zdo,      &
        xdby,xqc,xqrcd,xpwd,xpw,xhcd,xqcd,xhc,xqrc,xzu,xzd,            &

  ! cd  = detrainment function for updraft
  ! cdd = detrainment function for downdraft
  ! dellat = change of temperature per unit mass flux of cloud ensemble
  ! dellaq = change of q per unit mass flux of cloud ensemble
  ! dellaqc = change of qc per unit mass flux of cloud ensemble

        cd,cdd,DELLAH,DELLAQ,DELLAT,DELLAQC,          &
        u_cup,v_cup,uc,vc,ucd,vcd,dellu,dellv

  ! aa0 cloud work function for downdraft
  ! edt = epsilon
  ! aa0     = cloud work function without forcing effects
  ! aa1     = cloud work function with forcing effects
  ! xaa0    = cloud work function with cloud effects (ensemble dependent)
  ! edt     = epsilon

     real,    dimension (its:ite) ::                                   &
       edt,edto,edtx,AA1,AA0,XAA0,HKB,                          &
       HKBO,XHKB,QKB,QKBO,                                    &
       XMB,XPWAV,XPWEV,PWAV,PWEV,PWAVO,                                &
       PWEVO,BU,BUD,BUO,cap_max,                                   &
       cap_max_increment,closure_n,psum,psumh,sig,sigd,zuhe
     real,    dimension (its:ite) ::                                   &
        axx,edtmax,edtmin,entr_rate
     integer,    dimension (its:ite) ::                                &
       kzdown,KDET,K22,KB,JMIN,kstabi,kstabm,K22x,xland1,        &   !-lxz
       ktopdby,KBCONx,KBx,KTOPx,ierr,ierr2,ierr3,KBMAX

     integer,    dimension (its:ite), intent(in) :: csum
     integer                              ::                           &
       iloop,nall,iedt,nens,nens3,ki,I,K,KK,iresult
     real                                 ::                           &
      day,dz,dzo,mbdt,radius, &
      zcutdown,depth_min,zkbmax,z_detr,zktop,      &
      massfld,dh,cap_maxs,trash,frh,radiusd,frhd
      real detdo1,detdo2,entdo,dp,subin,detdo,entup,                &
      detup,subdown,entdoj,entupk,detupk,totmas
      real :: zusafe,xxx,xx1,xx2,zutop,power_entr,dzm1,dzp1


     integer :: jprnt,k1,k2,kbegzu,kdefi,kfinalzu,kstart,jmini
     logical :: keep_going,flg(its:ite)
     real tot_time_hr,xff_shal(9),blqe,xkshal
     character*50 :: ierrc(its:ite)
     real,    dimension (its:ite,kts:kte) ::                           &
       up_massentr,up_massdetr,dd_massentr,dd_massdetr,c1d             &
      ,up_massentro,up_massdetro,dd_massentro,dd_massdetro
     real,    dimension (its:ite,kts:kte) ::                           &
       up_massentru,up_massdetru
     real dts,fp,fpi,lambau,pgcon,pgc,pgcd,up_massent,up_massdet
     real,    dimension (kts:kte) :: smth
     real :: xff_mid(its:ite,2)
! rainevap from sas
      integer irainevap
      real zuh2(40)
      real, dimension (its:ite) :: deltv,rntot,delqev,delq2,delt,delq,qevap,rn,qcond
      real :: rain,t1,q1,elocp,evef,el2orc,evfact,evfactl,g_rain,e_dn,c_up
!srf- begin - DICYCLE
      integer :: iversion=1,ifim
      real :: denom,h_entr,umean,t_star
      real  :: outvars(its:ite,kts:kte,10)
      integer, intent(IN) :: DICYCLE
      real,    dimension (its:ite) :: aa1_bl,hkbo_bl,tau_bl,tau_ecmwf,wmean,tau_ecmwf_out
      real, dimension (its:ite,kts:kte) :: tn_bl, qo_bl, qeso_bl, heo_bl, heso_bl &
                                          ,qeso_cup_bl,qo_cup_bl, heo_cup_bl,heso_cup_bl&
                                          ,gammao_cup_bl,tn_cup_bl,hco_bl,DBYo_bl
      real, dimension(its:ite) :: xf_dicycle
      real, dimension(its:ite,10) :: forcing
      integer :: masscon
!srf- end

      xff_mid(:,:)=0.	! initialization is needed
      irainevap=1
      elocp=2.5e6/cp
      el2orc=2.5e6*2.5e6/(r_v*cp)
      evfact=.3
      evfactl=.3
      ifim=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!proportionality constant to estimate pressure gradient of updraft (Zhang and Wu, 2003, JAS
!
! ECMWF
       pgcon=0.
       lambau=2.
! SAS
!     lambau=0.
!     pgcon=-.55
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      pgcd=1.
      power_entr=2. ! 1.2
      day=86400.
!     cap_maxs=225.
!     if(imid.eq.1)cap_maxs=150.
      cap_maxs=150.
!     if(imid.eq.1)cap_maxs=100.
      do i=its,itf
        edto(i)=0.
        closure_n(i)=16.
        xmb_out(i)=0.
        cap_max(i)=cap_maxs
        cap_max_increment(i)=20.
        if(imid.eq.1)cap_max_increment(i)=10.
!
! for water or ice
!
        xland1(i)=int(xland(i)+.0001) ! 1.
        if(xland(i).gt.1.5 .or. xland(i).lt.0.5)then
            xland1(i)=0
            if(imid.eq.0)cap_max(i)=cap_maxs-25.
            if(imid.eq.1)cap_max(i)=cap_maxs-50.
            cap_max_increment(i)=20.
        else
            if(ztexec(i).gt.0.)cap_max(i)=cap_max(i)+25.
            if(ztexec(i).lt.0.)cap_max(i)=cap_max(i)-25.
        endif
        ierrc(i)=" "
!       cap_max_increment(i)=1.
      enddo
!
!--- initial entrainment rate (these may be changed later on in the
!--- program
!
      do i=its,ite
         c1d(i,:)= c1 ! 0. ! c1 ! max(.003,c1+float(csum(i))*.0001)
         entr_rate(i)=7.e-5   - min(20.,float(csum(i))) * 3.e-6
         if(xland1(i) == 0)entr_rate(i)=7.e-5
         if(imid.eq.1)entr_rate(i)=3.e-4
         if(imid.eq.1)c1d(i,:)=c1
         radius=.2/entr_rate(i)
         frh=3.14*radius*radius/dx/dx
         if(frh .gt. 0.7)then
            frh=.7
            radius=sqrt(frh*dx*dx/3.14)
            entr_rate(i)=.2/radius
         endif
         sig(i)=(1.-frh)**2
      enddo

      
!
!--- entrainment of mass
!
!
!--- initial detrainmentrates
!
      do k=kts,ktf
      do i=its,itf
        cnvwt(i,k)=0.
        zuo(i,k)=0.
        zdo(i,k)=0.
        hcot(i,k)=0.
        z(i,k)=zo(i,k)
        xz(i,k)=zo(i,k)
        cupclw(i,k)=0.
        cd(i,k)=1.e-9 ! 1.*entr_rate
        cdd(i,k)=1.e-9
        hcdo(i,k)=0.
        qrcdo(i,k)=0.
        dellaqc(i,k)=0.
      enddo
      enddo
!
!--- max/min allowed value for epsilon (ratio downdraft base mass flux/updraft
!    base mass flux
!
      edtmax(:)=1.
      if(imid.eq.1)edtmax(:)=.15
      edtmin(:)=.1
      if(imid.eq.1)edtmin(:)=.05
!
!--- minimum depth (m), clouds must have
!
      depth_min=1000.
      if(imid.eq.1)depth_min=500.
!
!--- maximum depth (mb) of capping 
!--- inversion (larger cap = no convection)
!
      DO i=its,itf
        kbmax(i)=1
        aa0(i)=0.
        aa1(i)=0.
        edt(i)=0.
        kstabm(i)=ktf-1
        IERR2(i)=0
        IERR3(i)=0
 enddo
!     do i=its,itf
!         cap_max(i)=cap_maxs
!         cap_max3(i)=25.
        iresult=0

!     enddo
!      write(11,*)'ipr,ips = ',ipr,its,cap_max(its)
!
!--- max height(m) above ground where updraft air can originate
!
      zkbmax=4000.
      if(imid.eq.1)zkbmax=2000.
!
!--- height(m) above which no downdrafts are allowed to originate
!
      zcutdown=3000.
!
!--- depth(m) over which downdraft detrains all its mass
!
      z_detr=1000.
!     if(imid.eq.1)z_detr=800.
!
      do nens=1,maxens
         mbdt_ens(nens)=(float(nens)-3.)*dtime*1.e-3+dtime*5.E-03
      enddo
         mbdt_ens(1)=10. !dtime*1.E-02
      do nens=1,maxens2
         edt_ens(nens)=.95-float(nens)*.01
      enddo
!
!--- environmental conditions, FIRST HEIGHTS
!
      do i=its,itf
         if(ierr(i).ne.20)then
            do k=1,maxens*maxens2*maxens3
               xf_ens(i,(iens-1)*maxens*maxens2*maxens3+k)=0.
               pr_ens(i,(iens-1)*maxens*maxens2*maxens3+k)=0.
            enddo
         endif
      enddo
!
!--- calculate moist static energy, heights, qes
!
      call cup_env(z,qes,he,hes,t,q,p,z1, &
           psur,ierr,tcrit,-1,   &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      call cup_env(zo,qeso,heo,heso,tn,qo,po,z1, &
           psur,ierr,tcrit,-1,   &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)

!
!--- environmental values on cloud levels
!
      call cup_env_clev(t,qes,q,he,hes,z,p,qes_cup,q_cup,he_cup, &
           hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur, &
           ierr,z1,          &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      call cup_env_clev(tn,qeso,qo,heo,heso,zo,po,qeso_cup,qo_cup, &
           heo_cup,heso_cup,zo_cup,po_cup,gammao_cup,tn_cup,psur,  &
           ierr,z1,          &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      if(ifim.eq.1.and.ierr(1).ne.0)return
      do i=its,itf
        if(ierr(i).eq.0)then
          if(kpbl(i).gt.5 .and. imid.eq.1)cap_max(i)=po_cup(i,kpbl(i))
          u_cup(i,kts)=us(i,kts)
          v_cup(i,kts)=vs(i,kts)
          do k=kts+1,ktf
           u_cup(i,k)=.5*(us(i,k-1)+us(i,k))
           v_cup(i,k)=.5*(vs(i,k-1)+vs(i,k))
          enddo
        endif
      enddo
      do i=its,itf
        if(ierr(i).eq.0)then
        if(aaeq(i).lt.-0.1)then
           ierr(i)=20
        endif
!     if(ierr(i).eq.0)then
!
      do k=kts,ktf
        if(zo_cup(i,k).gt.zkbmax+z1(i))then
          kbmax(i)=k
          go to 25
        endif
      enddo
 25   continue
!
!--- level where detrainment for downdraft starts
!
      do k=kts,ktf
        if(zo_cup(i,k).gt.z_detr+z1(i))then
          kdet(i)=k
          go to 26
        endif
      enddo
 26   continue
!
      endif
      enddo
!       if(j.eq.jpr)then
!       i=ipr
!!sms$ignore begin
!         write(12,*)'k,z(i,k),he(i,k),hes(i,k)'
!       do k=kts,ktf
!         write(12,*)k,z(i,k),he(i,k),hes(i,k)
!       enddo
!         write(12,*)'k,zo(i,k),heo(i,k),heso(i,k)'
!      do k=kts,ktf
!        write(12,*)k,zo(i,k),heo(i,k),heso(i,k)
!      enddo
!      call flush(12)
!!sms$ignore end
!       endif

!
!
!
!------- DETERMINE LEVEL WITH HIGHEST MOIST STATIC ENERGY CONTENT - K22
!
      k22(:)=1
      ktopdby(:)=0
      CALL cup_MAXIMI(HEO_CUP,1,KBMAX,K22,ierr, &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
       DO 36 i=its,itf
         IF(ierr(I).eq.0.)THEN
!         k22(i)=max(2,k22(i))
         IF(K22(I).GE.KBMAX(i))then
           ierr(i)=2
           ierrc(i)="could not find k22"
           ktop(i)=0
           k22(i)=0
           kbcon(i)=0
         endif
         endif
 36   CONTINUE
      if(ifim.eq.1.and.ierr(1).ne.0)return
!
!--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
!

      do i=its,itf
       IF(ierr(I).eq.0.)THEN
         hkb(i)=he_cup(i,k22(i))
         hkbo(i)=heo_cup(i,k22(i)) ! +float(use_excess)*(xl*zqexec(i)+cp*ztexec(i))
       endif ! ierr
      enddo
!     if(j.eq.jpr)write(0,*)'k22,kbmax,cap_inc,cap =',k22(ipr),kbmax(ipr),cap_max_increment,cap_max(ipr)
      jprnt=0
!      if(j.eq.jpr)jprnt=1
      iloop=1
!     if(imid.eq.1)iloop=5
      call cup_kbcon(ierrc,cap_max_increment,iloop,k22,kbcon,heo_cup,heso_cup, &
           hkbo,ierr,kbmax,po_cup,cap_max, &
           ztexec,zqexec,use_excess,       &
           jprnt,itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte, &
           z_cup,entr_rate,heo,imid)
!
!--- increase detrainment in stable layers
!
      CALL cup_minimi(HEso_cup,Kbcon,kstabm,kstabi,ierr,  &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      DO i=its,itf
         if(kbcon(i).lt.2)then
            ierr(i)=43
         endif
         if(kstabi(i).lt.kbcon(i))then
           kbcon(i)=1
           ierr(i)=42
         endif
         do k=kts,ktf
            entr_rate_2d(i,k)=entr_rate(i)
         enddo
         IF(ierr(I).eq.0.)THEN
            hkb(i)=he_cup(i,k22(i))
            hkbo(i)=heo_cup(i,k22(i)) ! +float(use_excess)*(xl*zqexec(i)+cp*ztexec(i))
            kdefi=0
            do k=kts,ktf
   ! Saulo's approach
               frh = min(qo_cup(i,k)/qeso_cup(i,k),1.)
               entr_rate_2d(i,k)=entr_rate(i) *(1.3-frh)
            enddo
          endif
       enddo
!
! get entrainment and detrainmentrates for updraft
!
      i=0
      if(j.eq.jpr)i=1
      call rates_up_pdf('deep',ktop,ierr,po_cup,entr_rate_2d,hkbo,heo,heso_cup,zo_cup, &
           xland1,i,kstabi,k22,kbcon,its,ite,itf,kts,kte,ktf,zuo,kpbl,ktopdby,csum)
! for mid level clouds we do not allow clouds taller than where stability
! changes
      if(imid.eq.1)call rates_up_pdf('mid',ktop,ierr,po_cup,entr_rate_2d,hkbo,heo,heso_cup,zo_cup, &
           xland1,i,kstabi,k22,kbcon,its,ite,itf,kts,kte,ktf,zuo,kpbl,ktopdby,csum)
!       if(j.eq.jpr)then
!!sms$ignore begin
!      do i=its,itf
!        write(12,*)"k,ktop(i),kstabi(i),kbcon(i),zuo(i,k)"
!      do k=kts,ktf
!        write(12,*)k,ktop(i),kstabi(i),kbcon(i),zuo(i,k)
!      enddo
!      enddo
!      call flush(12)
!!sms$ignore end
!       endif


!
! calculate mass entrainment and detrainment
!
      do k=kts,ktf
      do i=its,itf
         uc(i,k)=0.
         vc(i,k)=0.
         hc(i,k)=0.
         DBY(I,K)=0.
         hco(i,k)=0.
         DBYo(I,K)=0.
         up_massentro(i,k)=0.
         up_massdetro(i,k)=0.
         up_massentr(i,k)=0.
         up_massdetr(i,k)=0.
      enddo
      enddo
      do i=its,itf
       IF(ierr(I).eq.0.)THEN
         do k=1,kbcon(i)
            uc(i,k)=u_cup(i,k)
            vc(i,k)=v_cup(i,k)
         enddo
         do k=1,max(1,k22(i) -1)
            hc(i,k)=he_cup(i,k)
            hco(i,k)=heo_cup(i,k)
!           uc(i,k)=u_cup(i,k22(i))
!           vc(i,k)=v_cup(i,k22(i))
         enddo
       endif ! ierr
      enddo
!
!
      do i=its,itf
         if(ierr(i).eq.0)then
            hc(i,k22(i))=hkb(i)
            hco(i,k22(i))=hkbo(i)
!-- formulation 2
         if(k22(i).gt.1)then
            do k=1,k22(i) -1
              zuo(i,k)=0.
              zu(i,k)=0.
              xzu(i,k)=0.
            enddo
         endif
         do k=k22(i),ktop(i)
          xzu(i,k)= zuo(i,k)
          zu (i,k)= zuo(i,k)
          enddo
         do k=ktop(i)+1,ktf
           zuo(i,k)=0.
           zu(i,k)=0.
           xzu(i,k)=0.
          enddo

         !- mass entrainment and detrinament is defined on model levels
        
         do k=2,maxloc(zuo(i,:),1)
         !=> below maximum value zu -> change entrainment
           dz=zo_cup(i,k)-zo_cup(i,k-1)
        
           up_massdetro(i,k-1)=cd(i,k-1)*dz*zuo(i,k-1)
           up_massentro(i,k-1)=zuo(i,k)-zuo(i,k-1)+up_massdetro(i,k-1)
           if(up_massentro(i,k-1).lt.0.)then
              up_massentro(i,k-1)=0.
              up_massdetro(i,k-1)=zuo(i,k-1)-zuo(i,k)
              if(zuo(i,k-1).gt.0.)cd(i,k-1)=up_massdetro(i,k-1)/(dz*zuo(i,k-1))
           endif
           if(zuo(i,k-1).gt.0.)entr_rate_2d(i,k-1)=(up_massentro(i,k-1))/(dz*zuo(i,k-1))
         enddo
         do k=maxloc(zuo(i,:),1)+1,ktop(i)
         !=> above maximum value zu -> change detrainment
           dz=zo_cup(i,k)-zo_cup(i,k-1)
           up_massentro(i,k-1)=entr_rate_2d(i,k-1)*dz*zuo(i,k-1)
           up_massdetro(i,k-1)=zuo(i,k-1)+up_massentro(i,k-1)-zuo(i,k)
           if(up_massdetro(i,k-1).lt.0.)then
              up_massdetro(i,k-1)=0.
              up_massentro(i,k-1)=zuo(i,k)-zuo(i,k-1)
              if(zuo(i,k-1).gt.0.)entr_rate_2d(i,k-1)=(up_massentro(i,k-1))/(dz*zuo(i,k-1))
           endif
        
           if(zuo(i,k-1).gt.0.)cd(i,k-1)=up_massdetro(i,k-1)/(dz*zuo(i,k-1))
         enddo
        
         do k=2,ktf-1
          up_massentr(i,k-1)=up_massentro(i,k-1)
          up_massdetr(i,k-1)=up_massdetro(i,k-1)
          up_massentru(i,k-1)=up_massentro(i,k-1)+lambau*up_massdetro(i,k-1)
          up_massdetru(i,k-1)=up_massdetro(i,k-1)+lambau*up_massdetro(i,k-1)
         enddo
!
!   NOTE: Ktop here already includes overshooting, ktopdby is without
!   overshooting
!
        do k=k22(i)  +1,ktop(i)  !mass cons option
          denom=zuo(i,k-1)-.5*up_massdetro(i,k-1)+up_massentro(i,k-1)
          if(denom.lt.1.e-8)then
!          write(14,*)j,k,denom,zuo(i,k-1),up_massdetro(i,k-1),up_massentro(i,k-1)
           ierr(i)=51
           exit
          endif

          hc(i,k)=(hc(i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k-1)*hc(i,k-1)+ &
                         up_massentr(i,k-1)*he(i,k-1))   /            &
                         (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
          uc(i,k)=(uc(i,k-1)*zu(i,k-1)-.5*up_massdetru(i,k-1)*uc(i,k-1)+ &
                          up_massentru(i,k-1)*us(i,k-1) &
                         -pgcon*.5*(zu(i,k)+zu(i,k-1))*(u_cup(i,k)-u_cup(i,k-1))) /     &
                         (zu(i,k-1)-.5*up_massdetru(i,k-1)+up_massentru(i,k-1))
          vc(i,k)=(vc(i,k-1)*zu(i,k-1)-.5*up_massdetru(i,k-1)*vc(i,k-1)+                &
                          up_massentru(i,k-1)*vs(i,k-1)      &
                         -pgcon*.5*(zu(i,k)+zu(i,k-1))*(v_cup(i,k)-v_cup(i,k-1))) /    &
                         (zu(i,k-1)-.5*up_massdetru(i,k-1)+up_massentru(i,k-1))
          dby(i,k)=hc(i,k)-hes_cup(i,k)
          hco(i,k)=(hco(i,k-1)*zuo(i,k-1)-.5*up_massdetro(i,k-1)*hco(i,k-1)+ &
                         up_massentro(i,k-1)*heo(i,k-1))   /            &
                         (zuo(i,k-1)-.5*up_massdetro(i,k-1)+up_massentro(i,k-1))
          dbyo(i,k)=hco(i,k)-heso_cup(i,k)
         enddo
41       continue
         if(ktop(i).lt.kbcon(i)+2)then
            ierr(i)=5
            ierrc(i)='ktop too small'
            ktop(i)=0
         endif
         do k=ktop(i)+1,ktf
           HC(i,K)=hes_cup(i,k)
           UC(i,K)=u_cup(i,k)
           VC(i,K)=v_cup(i,k)
           HCo(i,K)=heso_cup(i,k)
           DBY(I,K)=0.
           DBYo(I,K)=0.
           zu(i,k)=0.
           zuo(i,k)=0.
           cd(i,k)=0.
           entr_rate_2d(i,k)=0.
           up_massentr(i,k)=0.
           up_massdetr(i,k)=0.
           up_massentro(i,k)=0.
           up_massdetro(i,k)=0.
         enddo
      endif
      enddo
!
      if(ifim.eq.1.and.ierr(1).ne.0)return
      DO 37 i=its,itf
         kzdown(i)=0
         if(ierr(i).eq.0)then
            zktop=(zo_cup(i,ktop(i))-z1(i))*.6
            if(imid.eq.1)zktop=(zo_cup(i,ktop(i))-z1(i))*.4
            zktop=min(zktop+z1(i),zcutdown+z1(i))
            do k=kts,ktf
              if(zo_cup(i,k).gt.zktop)then
                 kzdown(i)=k
                 kzdown(i)=min(kzdown(i),kstabi(i)-1)  !
                 go to 37
              endif
              enddo
         endif
 37   CONTINUE
!
!--- DOWNDRAFT ORIGINATING LEVEL - JMIN
!
      call cup_minimi(HEso_cup,K22,kzdown,JMIN,ierr, &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      DO 100 i=its,itf
         IF(ierr(I).eq.0.)THEN
!
!--- check whether it would have buoyancy, if there where
!--- no entrainment/detrainment
!
         jmini = jmin(i)
         keep_going = .TRUE.
         do while ( keep_going )
           keep_going = .FALSE.
           if ( jmini - 1 .lt. kdet(i)   ) kdet(i) = jmini-1
           if ( jmini     .ge. ktop(i)-1 ) jmini = ktop(i) - 2
           ki = jmini
           hcdo(i,ki)=heso_cup(i,ki)
           DZ=Zo_cup(i,Ki+1)-Zo_cup(i,Ki)
           dh=0.
           do k=ki-1,1,-1
             hcdo(i,k)=heso_cup(i,jmini)
             DZ=Zo_cup(i,K+1)-Zo_cup(i,K)
             dh=dh+dz*(HCDo(i,K)-heso_cup(i,k))
             if(dh.gt.0.)then
               jmini=jmini-1
               if ( jmini .gt. 5 ) then
                 keep_going = .TRUE.
               else
                 ierr(i) = 9
                 ierrc(i) = "could not find jmini9"
                 exit
               endif
             endif
           enddo
         enddo
         jmin(i) = jmini 
         if ( jmini .le. 5 ) then
           ierr(i)=4
           ierrc(i) = "could not find jmini4"
         endif
       ENDIF
100   continue
!
      if(ifim.eq.1.and.ierr(1).ne.0)return
! - Must have at least depth_min m between cloud convective base
!     and cloud top.
!
      do i=its,itf
         IF(ierr(I).eq.0.)THEN
            if ( jmin(i) - 1 .lt. kdet(i)   ) kdet(i) = jmin(i)-1
            IF(-zo_cup(I,KBCON(I))+zo_cup(I,KTOP(I)).LT.depth_min)then
               ierr(i)=6
               ierrc(i)="cloud depth very shallow"
            endif
         endif
      enddo
      if(ifim.eq.1.and.ierr(1).ne.0)return

!
!--- normalized downdraft mass flux profile,also work on bottom detrainment
!--- in this routine
!
      do k=kts,ktf
      do i=its,itf
       zd(i,k)=0.
       zdo(i,k)=0.
       cdd(i,k)=0.
       dd_massentr(i,k)=0.
       dd_massdetr(i,k)=0.
       dd_massentro(i,k)=0.
       dd_massdetro(i,k)=0.
       hcdo(i,k)=heso_cup(i,k)
       ucd(i,k)=u_cup(i,k)
       vcd(i,k)=v_cup(i,k)
       dbydo(i,k)=0.
      enddo
      enddo
      do i=its,itf
        beta=max(.025,.055-float(csum(i))*.0015)  !.02
        mentrd_rate_2d(i,:)=entr_rate(i)
        if(imid.eq.0 .and. xland1(i) == 0)then
!             beta=.01
              edtmax(i)=max(0.1,.4-float(csum(i))*.015) !.3)
        endif
        if(imid.eq.1)beta=.025
        bud(i)=0.
        IF(ierr(I).eq.0)then
        cdd(i,1:jmin(i))=1.e-9
        cdd(i,jmin(i))=0.
        dd_massdetro(i,:)=0.
        dd_massentro(i,:)=0.
           call get_zu_zd_pdf_fim(xland1(i),zuh2,"DOWN",ierr(i),kdet(i),jmin(i),zdo(i,:),kts,kte,ktf,beta,kpbl(i),csum(i))
        if(zdo(i,jmin(i)) .lt.1.e-8)then
          zdo(i,jmin(i))=0.
          jmin(i)=jmin(i)-1
          if(zdo(i,jmin(i)) .lt.1.e-8)then
             ierr(i)=876
             exit
          endif
        endif
        xzd(i,:)= zdo(i,:)
        zd (i,:)= zdo(i,:)
        do ki=jmin(i)  ,maxloc(zdo(i,:),1),-1
!-srf mass cons
          !=> from jmin to maximum value zd -> change entrainment
          dzo=zo_cup(i,ki+1)-zo_cup(i,ki)
          dd_massdetro(i,ki)=cdd(i,ki)*dzo*zdo(i,ki+1)
          dd_massentro(i,ki)=zdo(i,ki)-zdo(i,ki+1)+dd_massdetro(i,ki)
          if(dd_massentro(i,ki).lt.0.)then
             dd_massentro(i,ki)=0.
             dd_massdetro(i,ki)=zdo(i,ki+1)-zdo(i,ki)
             if(zdo(i,ki+1).gt.0.)cdd(i,ki)=dd_massdetro(i,ki)/(dzo*zdo(i,ki+1))
          endif
          if(zdo(i,ki+1).gt.0.)mentrd_rate_2d(i,ki)=dd_massentro(i,ki)/(dzo*zdo(i,ki+1))
        enddo
        mentrd_rate_2d(i,1)=0.
        do ki=maxloc(zdo(i,:),1)-1,1,-1
          !=> from maximum value zd to surface -> change detrainment
          dzo=zo_cup(i,ki+1)-zo_cup(i,ki)
          dd_massentro(i,ki)=mentrd_rate_2d(i,ki)*dzo*zdo(i,ki+1)
          dd_massdetro(i,ki) = zdo(i,ki+1)+dd_massentro(i,ki)-zdo(i,ki)
          if(dd_massdetro(i,ki).lt.0.)then
            dd_massdetro(i,ki)=0.
            dd_massentro(i,ki)=zdo(i,ki)-zdo(i,ki+1)
            if(zdo(i,ki+1).gt.0.)mentrd_rate_2d(i,ki)=dd_massentro(i,ki)/(dzo*zdo(i,ki+1))
          endif
          if(zdo(i,ki+1).gt.0.)cdd(i,ki)= dd_massdetro(i,ki)/(dzo*zdo(i,ki+1))
        enddo


!         write(92,*)'k,zdo(i,k),dd_massentro(i,k),dd_massdetro(i,k)'
        do k=jmin(i),1,-1
          xzd(i,k)= zdo(i,k)
          zd (i,k)= zdo(i,k)
          dd_massentr(i,k)=dd_massentro(i,k)
          dd_massdetr(i,k)=dd_massdetro(i,k)

!         write(92,*)k,zdo(i,k),dd_massentro(i,k),dd_massdetro(i,k)
        enddo
!        do k=1,jmin(i)
!         c1d(i,k)=0.
!        enddo
!        do k=jmin(i)-1,ktop(i)
!          dz=zo_cup(i,ktop(i))-zo_cup(i,jmin(i))
!          c1d(i,k)=c1+c1*40.*(zo_cup(i,k+1)-zo_cup(i,k))/dz
!          if(imid.eq.1)c1d(i,k)=c1+c1*10.*(zo_cup(i,k+1)-zo_cup(i,k))/dz
!          c1d(i,k)=max(0.,c1d(i,k))
!          c1d(i,k)=min(.002,c1d(i,k))
!        enddo
        do k=kts,ktop(i)+1
          if(p_cup(i,k).gt.600.)then
            c1d(i,k)=0.
          elseif(p_cup(i,k).gt.400. .and. p_cup(i,k).le.600.)then
            c1d(i,k)=c1_max*(600.-p_cup(i,k))/200.
            c1d(i,k)=max(0.,c1d(i,k))
          elseif(p_cup(i,k).le.400.)then
            c1d(i,k)=c1_max
          endif
        enddo


! downdraft moist static energy + moisture budget
            dbydo(i,jmin(i))=hcdo(i,jmin(i))-heso_cup(i,jmin(i))
            bud(i)=dbydo(i,jmin(i))*(zo_cup(i,jmin(i)+1)-zo_cup(i,jmin(i)))
            do ki=jmin(i)  ,1,-1
!-srf mass cons
             dzo=zo_cup(i,ki+1)-zo_cup(i,ki)
!             h_entr=heo(i,ki)
             h_entr=.5*(heo(i,ki)+.5*(hco(i,ki)+hco(i,ki+1)))
             ucd(i,ki)=(ucd(i,ki+1)*zdo(i,ki+1)                       &
                         -.5*dd_massdetro(i,ki)*ucd(i,ki+1)+ &
                        dd_massentro(i,ki)*us(i,ki)          &
                        -pgcon*zdo(i,ki+1)*(us(i,ki+1)-us(i,ki)))   /     &
                        (zdo(i,ki+1)-.5*dd_massdetro(i,ki)+dd_massentro(i,ki))
             vcd(i,ki)=(vcd(i,ki+1)*zdo(i,ki+1)                       &
                         -.5*dd_massdetro(i,ki)*vcd(i,ki+1)+ &
                        dd_massentro(i,ki)*vs(i,ki)            &
                        -pgcon*zdo(i,ki+1)*(vs(i,ki+1)-vs(i,ki)))   /  &
                        (zdo(i,ki+1)-.5*dd_massdetro(i,ki)+dd_massentro(i,ki))
             hcdo(i,ki)=(hcdo(i,ki+1)*zdo(i,ki+1)                       &
                         -.5*dd_massdetro(i,ki)*hcdo(i,ki+1)+ &
                        dd_massentro(i,ki)*h_entr)   /            &
                        (zdo(i,ki+1)-.5*dd_massdetro(i,ki)+dd_massentro(i,ki))
             dbydo(i,ki)=hcdo(i,ki)-heso_cup(i,ki)
!            if(i.eq.ipr.and.j.eq.jpr)write(0,*)'ki,bud = ',ki,bud(i),hcdo(i,ki)
             bud(i)=bud(i)+dbydo(i,ki)*dzo
            enddo
          endif

        if(bud(i).gt.0)then
          ierr(i)=7
          ierrc(i)='downdraft is not negatively buoyant '
        endif
      enddo
      if(ifim.eq.1.and.ierr(1).ne.0)return
!          if(j.eq.jpr)then
!!sms$ignore begin
!          do i=its,itf
!            write(12,*)'k22,jmin,kbcon,ktop,ierr =',k22(i),jmin(i),kbcon(i),ktop(i),ierr(i)
!            write(12,*)'k,p_cup(i,k),zuo(i,k),hco(i,k)'
!            do k=kts,ktop(i)
!              write(12,*)k,p_cup(i,k),zuo(i,k),hco(i,k)
!            enddo
!            write(12,*)'k,p_cup(i,k),heo_cup(i,k),heso_cup(i,k)'
!            do k=kts,ktop(i)
!              write(12,*)k,p_cup(i,k),heo_cup(i,k),heso_cup(i,k)
!            enddo
!            write(12,*)'k,p_cup(i,k),zdo(i,k),hcdo(i,k)'
!            do k=jmin(i),1
!              write(12,*)k,p_cup(i,k),zdo(i,k),hcdo(i,k)
!            enddo
!         enddo
!      call flush(12)
!!sms$ignore end
!         endif
!
!--- calculate moisture properties of downdraft
!
      call cup_dd_moisture_new(ierrc,zdo,hcdo,heso_cup,qcdo,qeso_cup, &
           pwdo,qo_cup,zo_cup,dd_massentro,dd_massdetro,jmin,ierr,gammao_cup, &
           pwevo,bu,qrcdo,qo,heo,tn_cup,1, &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      if(ifim.eq.1.and.ierr(1).ne.0)return
!
!--- calculate moisture properties of updraft
!
      call cup_up_moisture('deep',ierr,zo_cup,qco,qrco,pwo,pwavo, &
           ccnclean,p_cup,kbcon,ktop,cd,dbyo,clw_all, &
           t_cup,qo,GAMMAo_cup,zuo,qeso_cup,k22,qo_cup,        &
           ZQEXEC,use_excess,ccn,rho,c1d,up_massentr,up_massdetr,psum,psumh,&
           autoconv,aeroevap,1,itf,jtf,ktf,j,ipr,jpr, &
           its,ite, jts,jte, kts,kte)
      if(ifim.eq.1.and.ierr(1).ne.0)return
      do i=its,itf
      if(ierr(i).eq.0)then
      do k=kts+1,ktop(i)
         dp=100.*(po_cup(i,1)-po_cup(i,2))
         cupclw(i,k)=qrco(i,k)	! my mod
         cnvwt(i,k)=zuo(i,k)*cupclw(i,k)*g/dp
      enddo
      endif
      enddo
!
!--- calculate workfunctions for updrafts
!
      call cup_up_aa0(aa0,z,zu,dby,GAMMA_CUP,t_cup, &
           kbcon,ktop,ierr,           &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      call cup_up_aa0(aa1,zo,zuo,dbyo,GAMMAo_CUP,tn_cup, &
           kbcon,ktop,ierr,           &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      do i=its,itf
         if(ierr(i).eq.0)then
           if(aa1(i).eq.0.)then
               ierr(i)=17
               ierrc(i)="cloud work function zero"
           endif
           aaeq(i)=aa1(i)
         endif
      enddo
      if(ifim.eq.1.and.ierr(1).ne.0)return
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=================================================================
!-srf- begin
      !--- AA1 from boundary layer (bl) processes only
      aa1_bl    (:) = 0.0
      xf_dicycle(:) = 0.0
      tau_ecmwf(:)  = 0.
      tau_ecmwf_out(:)  = 0.
      !- way to calculate the fraction of cape consumed by shallow convection
      iversion=1 ! ecmwf  
      !iversion=0 ! orig    
      !
      ! Betchold et al 2008 time-scale of cape removal
!
! wmean is of no meaning over land....
! still working on replacing it over water
!
      DO i=its,itf
            if(ierr(i).eq.0)then
                !- mean vertical velocity 
                wmean(i) = 7.0 ! m/s ! in the future change for Wmean == integral( W dz) / cloud_depth
                !- time-scale cape removal from  Betchold et al. 2008
                tau_ecmwf(i)=( zo_cup(i,ktopdby(i))- zo_cup(i,kbcon(i)) ) / wmean(i) 
                tau_ecmwf(i)= tau_ecmwf(i) * (1.0061 + 1.23E-2 * (dx/1000.))! dx must be in meters 
                tau_ecmwf_out(i)= tau_ecmwf(i)
            endif
      enddo
      tau_bl(:)     = 0.
      !
      IF(dicycle == 1) then
        DO i=its,itf
            
            if(ierr(i).eq.0)then
                if(xland1(i) ==  0 ) then
                  !- over water
                  umean= 2.0+sqrt(2.0*(US(i,1)**2+VS(i,1)**2+US(i,kbcon(i))**2+VS(i,kbcon(i))**2))
                  tau_bl(i) = (zo_cup(i,kbcon(i))- z1(i)) /umean        
                else
                  !- over land
                  tau_bl(i) =( zo_cup(i,ktopdby(i))- zo_cup(i,kbcon(i)) ) / wmean(i)
                endif

            endif
        ENDDO

        if(iversion == 1) then 
	!-- version ecmwf
        t_star=4.  !original =1

           !-- calculate pcape from BL forcing only
            call cup_up_aa1bl(aa1_bl,t,tn,q,qo,dtime, &
	  		    zo_cup,zuo,dbyo_bl,GAMMAo_CUP_bl,tn_cup_bl, &
          		    kbcon,ktop,ierr,	       &
          		    itf,jtf,ktf,its,ite, jts,jte, kts,kte)

            DO i=its,itf

            if(ierr(i).eq.0)then

!- only for convection rooting in the PBL
               if(zo_cup(i,kbcon(i))-z1(i) > zo(i,kpbl(i)-1)) then 
                  aa1_bl(i) = 0.0
               else
!	!- multiply aa1_bl the " time-scale" - tau_bl
                  aa1_bl(i) = max(0.,aa1_bl(i)/t_star* tau_bl(i))
               endif 
            endif
            ENDDO
            
	else
	
	  !- version for real cloud-work function
	  
          !-get the profiles modified only by bl tendencies
          DO i=its,itf
           tn_bl(i,:)=0.;qo_bl(i,:)=0.
           if ( ierr(i) == 0 )then
            !below kbcon -> modify profiles
            tn_bl(i,1:kbcon(i)) = tn(i,1:kbcon(i))
            qo_bl(i,1:kbcon(i)) = qo(i,1:kbcon(i))
     	    !above kbcon -> keep environment profiles
            tn_bl(i,kbcon(i)+1:ktf) = t(i,kbcon(i)+1:ktf)
            qo_bl(i,kbcon(i)+1:ktf) = q(i,kbcon(i)+1:ktf)
           endif 
          ENDDO
          !--- calculate moist static energy, heights, qes, ... only by bl tendencies
          call cup_env(zo,qeso_bl,heo_bl,heso_bl,tn_bl,qo_bl,po,z1, &
                     psur,ierr,tcrit,-1,   &
                     itf,jtf,ktf, its,ite, jts,jte, kts,kte)
          !--- environmental values on cloud levels only by bl tendencies
          call cup_env_clev(tn_bl,qeso_bl,qo_bl,heo_bl,heso_bl,zo,po,qeso_cup_bl,qo_cup_bl, &
          		    heo_cup_bl,heso_cup_bl,zo_cup,po_cup,gammao_cup_bl,tn_cup_bl,psur,  &
          		    ierr,z1,	       &
          		    itf,jtf,ktf,its,ite, jts,jte, kts,kte)
          DO i=its,itf
            IF(ierr(I).eq.0.)THEN
               hkbo_bl(i)=heo_cup_bl(i,k22(i)) ! +float(use_excess)*(xl*zqexec(i)+cp*ztexec(i))
            endif ! ierr
          ENDDO
          DO k=kts,ktf
           do i=its,itf
             hco_bl (i,k)=0.
             DBYo_bl(i,k)=0.
           enddo
          ENDDO
          DO i=its,itf
            IF(ierr(I).eq.0.)THEN
             do k=1,kbcon(i)-1
              hco_bl(i,k)=hkbo_bl(i)
             enddo
             k=kbcon(i)
             hco_bl (i,k)=hkbo_bl(i)
             DBYo_bl(i,k)=Hkbo_bl(i) - HESo_cup_bl(i,k)
            ENDIF
          ENDDO
!	  
!	  
          DO i=its,itf
            if(ierr(i).eq.0)then
               do k=kbcon(i)+1,ktop(i)
          	  hco_bl(i,k)=(hco_bl(i,k-1)*zuo(i,k-1)-.5*up_massdetro(i,k-1)*hco_bl(i,k-1)+ &
          		     up_massentro(i,k-1)*heo_bl(i,k-1))   /	       &
          		     (zuo(i,k-1)-.5*up_massdetro(i,k-1)+up_massentro(i,k-1))
          	  dbyo_bl(i,k)=hco_bl(i,k)-heso_cup_bl(i,k)
               enddo
               do k=ktop(i)+1,ktf
                  hco_bl (i,k)=heso_cup_bl(i,k)
                  dbyo_bl(i,k)=0.0
               enddo
            endif
          ENDDO
        
          !--- calculate workfunctions for updrafts
          call cup_up_aa0(aa1_bl,zo,zuo,dbyo_bl,GAMMAo_CUP_bl,tn_cup_bl, &
                        kbcon,ktop,ierr,           &
                        itf,jtf,ktf,its,ite, jts,jte, kts,kte)

          DO i=its,itf
            
            if(ierr(i).eq.0)then
                !- get the increment on AA0 due the BL processes
                aa1_bl(i) = aa1_bl(i) - aa0(i)
                !- only for convection rooting in the PBL
                !if(zo_cup(i,kbcon(i))-z1(i) > 500.0) then !- instead 500 -> zo_cup(kpbl(i))
                !   aa1_bl(i) = 0.0
                !else
                !   !- multiply aa1_bl the "normalized time-scale" - tau_bl/ model_timestep
                   aa1_bl(i) = aa1_bl(i)* tau_bl(i)/ dtime
                !endif 
                print*,'aa0,aa1bl=',aa0(i),aa1_bl(i),aa0(i)-aa1_bl(i),tau_bl(i)!,dtime,xland(i)     
            endif
           ENDDO
	ENDIF
     ENDIF  ! version of implementation
	!i=maxloc(aa1_bl,1)
        !if(aa1_bl(i)>0. .and. ierr(i) == 0) print*,'aa0,aa1bl=',j,i,aa0(i),aa1_bl(i), xland(i) &
	!                    ,tau_bl(i) ; call flush(6)!,tau_bl(i),dtime,xland(i)     
	
        !aa1_bl(:)=0.0 !tmp
!-srf -end
!=================================================================


       axx(:)=aa1(:)

!
!--- DETERMINE DOWNDRAFT STRENGTH IN TERMS OF WINDSHEAR
!
      call cup_dd_edt(ierr,us,vs,zo,ktop,kbcon,edt,po,pwavo, &
           pwo,ccn,pwevo,edtmax,edtmin,maxens2,edtc,psum,psumh, &
           ccnclean,rho,aeroevap,itf,jtf,ktf,j,ipr,jpr, &
           its,ite, jts,jte, kts,kte)
      if(ifim.eq.1.and.ierr(1).ne.0)return
      do 250 iedt=1,maxens2
        do i=its,itf
         if(ierr(i).eq.0)then
!        edt(i)=edtc(i,iedt)
         edto(i)=edtc(i,iedt)
         edtx(i)=edtc(i,iedt)
         if(maxens2.eq.3)then
!           edt(i)=edtc(i,3)
            edto(i)=edtc(i,3)
            edtx(i)=edtc(i,3)
         endif
         endif
        enddo
        do k=kts,ktf
        do i=its,itf
           dellat_ens(i,k,iedt)=0.
           dellaq_ens(i,k,iedt)=0.
           dellaqc_ens(i,k,iedt)=0.
           pwo_ens(i,k,iedt)=0.
        enddo
        enddo
!
!--- change per unit mass that a model cloud would modify the environment
!
!--- 1. in bottom layer
!
      do k=kts,ktf
      do i=its,itf
        dellu(i,k)=0.
        dellv(i,k)=0.
        dellah(i,k)=0.
        dellat(i,k)=0.
        dellaq(i,k)=0.
        dellaqc(i,k)=0.
      enddo
      enddo
!
!----------------------------------------------  cloud level ktop
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level ktop-1
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!
!----------------------------------------------  cloud level k+2
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level k+1
!
!----------------------------------------------  cloud level k+1
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level k
!
!----------------------------------------------  cloud level k
!
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!
!----------------------------------------------  cloud level 3
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level 2
!
!----------------------------------------------  cloud level 2
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level 1

      do i=its,itf
        if(ierr(i).eq.0)then
         dp=100.*(po_cup(i,1)-po_cup(i,2))
         dellu(i,1)=pgcd*(edto(i)*zdo(i,2)*ucd(i,2)   &
                     -edto(i)*zdo(i,2)*u_cup(i,2))*g/dp
         dellv(i,1)=pgcd*(edto(i)*zdo(i,2)*vcd(i,2)   &
                     -edto(i)*zdo(i,2)*v_cup(i,2))*g/dp

         do k=kts+1,ktop(i)
            ! these three are only used at or near mass detrainment and/or entrainment levels
            pgc=pgcon
            entupk=0.
            detupk=0.
            entdoj=0.
            ! detrainment and entrainment for fowndrafts
            detdo=edto(i)*dd_massdetro(i,k)
            entdo=edto(i)*dd_massentro(i,k)
            ! entrainment/detrainment for updraft
            entup=up_massentro(i,k)
            detup=up_massdetro(i,k)
            ! subsidence by downdrafts only
            subin=-zdo(i,k+1)*edto(i)
            subdown=-zdo(i,k)*edto(i)
!
            !         SPECIAL LEVELS
            ! downdraft originating level, similiar to k22-1 for updraft
!-srf- 21jul2015 mass con            
!            if(k.eq.jmin(i))then
!              entdoj=edto(i)*zdo(i,k)
!           endif
!-srf- mass con
            if(k.eq.ktop(i))then
               detupk=zuo(i,ktop(i))
               subin=0.
               subdown=0.
               detdo=0.
               entdo=0.
               entup=0.
               detup=0.
            endif
            totmas=subin-subdown+detup-entup-entdo+ &
                   detdo-entupk-entdoj+detupk+zuo(i,k+1)-zuo(i,k)
            if(abs(totmas).gt.1.e-6)then
               write(0,*)'*********************',i,j,k,totmas
123     formAT(1X,i2,8E12.4)
            endif
            dp=100.*(po_cup(i,k)-po_cup(i,k+1))
             pgc=pgcon
            if(k.ge.ktop(i))pgc=0.

             dellu(i,k) =-(zuo(i,k+1)*(uc (i,k+1)-u_cup(i,k+1) ) - &
                            zuo(i,k  )*(uc (i,k  )-u_cup(i,k  ) ) )*g/dp &
                         +(zdo(i,k+1)*(ucd(i,k+1)-u_cup(i,k+1) ) - &
                            zdo(i,k  )*(ucd(i,k  )-u_cup(i,k  ) ) )*g/dp*edto(i)*pgcd &
                     +(pgc*.5*(zuo(i,k)+zuo(i,k+1))*(u_cup(i,k+1)-u_cup(i,k)) &
                     + pgcd*pgc*.5*edto(i)*(zdo(i,k)+zdo(i,k+1))*(u_cup(i,k+1)-u_cup(i,k)))*g/dp
             dellv(i,k) =-(zuo(i,k+1)*(vc (i,k+1)-v_cup(i,k+1) ) - &
                            zuo(i,k  )*(vc (i,k  )-v_cup(i,k  ) ) )*g/dp &
                         +(zdo(i,k+1)*(vcd(i,k+1)-v_cup(i,k+1) ) - &
                            zdo(i,k  )*(vcd(i,k  )-v_cup(i,k  ) ) )*g/dp*edto(i)*pgcd &
                     +(pgc*.5*(zuo(i,k)+zuo(i,k+1))*(v_cup(i,k+1)-v_cup(i,k)) &
                     + pgcd*pgc*.5*edto(i)*(zdo(i,k)+zdo(i,k+1))*(v_cup(i,k+1)-v_cup(i,k)))*g/dp
 
       enddo   ! k

      endif
    enddo


    do i=its,itf
        !trash  = 0.0
        !trash2 = 0.0
        if(ierr(i).eq.0)then

	 dp=100.*(po_cup(i,1)-po_cup(i,2))

	 dellah(i,1)=(edto(i)*zdo(i,2)*hcdo(i,2)   &
                     -edto(i)*zdo(i,2)*heo_cup(i,2))*g/dp

	 dellaq (i,1)=(edto(i)*zdo(i,2)*qcdo(i,2)   &
                      -edto(i)*zdo(i,2)*qo_cup(i,2))*g/dp
	
	 G_rain=  0.5*(pwo (i,1)+pwo (i,2))*g/dp
	 E_dn  = -0.5*(pwdo(i,1)+pwdo(i,2))*g/dp*edto(i)  ! pwdo < 0 and E_dn must > 0
         dellaq(i,1) = dellaq(i,1)+ E_dn-G_rain
         
	 !--- conservation check
	 !- water mass balance
	 !trash = trash  + (dellaq(i,1)+dellaqc(i,1)+G_rain-E_dn)*dp/g 	 
	 !- H  budget
	 !trash2 = trash2+ (dellah(i,1))*dp/g
	 

!           if(J.eq.312772)write(13,*)'k,dellah(i,k),dellaq(i,k),edto(i,k),dellaqc(i,k)'
	 do k=kts+1,ktop(i)
            dp=100.*(po_cup(i,k)-po_cup(i,k+1))
            ! these three are only used at or near mass detrainment and/or entrainment levels

            dellah(i,k) =-(zuo(i,k+1)*(hco (i,k+1)-heo_cup(i,k+1) ) -                 &
                           zuo(i,k  )*(hco (i,k  )-heo_cup(i,k  ) ) )*g/dp            &
			 +(zdo(i,k+1)*(hcdo(i,k+1)-heo_cup(i,k+1) ) -                 &
                           zdo(i,k  )*(hcdo(i,k  )-heo_cup(i,k  ) ) )*g/dp*edto(i)
			

	    !- check H conservation 
	    ! trash2 = trash2+ (dellah(i,k))*dp/g
	
	
            !-- take out cloud liquid water for detrainment
!            detup=up_massdetro(i,k)
            dz=zo_cup(i,k+1)-zo_cup(i,k)
!           if(imid.eq.1)then
	    dellaqc(i,k) = zuo(i,k)*c1d(i,k)*qrco(i,k)*dz/dp*g !detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp
!           else
!            dellaqc(i,k)=  detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp
!           endif
            if(k.eq.ktop(i))dellaqc(i,k)= detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp
	    !---
	    G_rain=  0.5*(pwo (i,k)+pwo (i,k+1))*g/dp
	    E_dn  = -0.5*(pwdo(i,k)+pwdo(i,k+1))*g/dp*edto(i) ! pwdo < 0 and E_dn must > 0
            !-- condensation source term = detrained + flux divergence of
            !-- cloud liquid water (qrco) + converted to rain
	
            C_up = dellaqc(i,k)+(zuo(i,k+1)* qrco(i,k+1) -       &
                                 zuo(i,k  )* qrco(i,k  )  )*g/dp + G_rain
            !-- water vapor budget
	    !-- = flux divergence z*(Q_c - Q_env)_up_and_down &
	    !--   - condensation term + evaporation
            dellaq(i,k) =-(zuo(i,k+1)*(qco (i,k+1)-qo_cup(i,k+1) ) -                 &
                           zuo(i,k  )*(qco (i,k  )-qo_cup(i,k  ) ) )*g/dp            &
			 +(zdo(i,k+1)*(qcdo(i,k+1)-qo_cup(i,k+1) ) -                 &
                           zdo(i,k  )*(qcdo(i,k  )-qo_cup(i,k  ) ) )*g/dp*edto(i)    &
                         - C_up + E_dn
	    !- check water conservation liq+condensed (including rainfall)
            ! trash= trash+ (dellaq(i,k)+dellaqc(i,k)+ G_rain-E_dn)*dp/g
!            if(J.eq.312772)write(13,444)k,dellaq(i,k),dellaqc(i,k),c_up,g_rain,e_dn,qrco(i,k)

         enddo   ! k
        endif

      enddo
444   format(1x,i2,1x,6e12.4) !,1x,f7.2,2x,e13.5)
!
!--- using dellas, calculate changed environmental profiles
!
!     do 200 nens=1,maxens
      mbdt=mbdt_ens(1)
      do i=its,itf
      xaa0_ens(i,:)=0.
      enddo

      do i=its,itf
         if(ierr(i).eq.0)then
           do k=kts,ktf
            XHE(I,K)=DELLAH(I,K)*MBDT+HEO(I,K)
            XQ(I,K)=max(1.e-16,DELLAQ(I,K)*MBDT+QO(I,K))
            DELLAT(I,K)=(1./cp)*(DELLAH(I,K)-xlv*DELLAQ(I,K))
            XT(I,K)= DELLAT(I,K)*MBDT+TN(I,K)
            xt(i,k)=max(190.,xt(i,k))
           enddo
         ENDIF
      enddo
      do i=its,itf
      if(ierr(i).eq.0)then
      XHKB(I)=DELLAH(I,K22(i))*MBDT+HKBO(I)
      XHE(I,ktf)=HEO(I,ktf)
      XQ(I,ktf)=QO(I,ktf)
      XT(I,ktf)=TN(I,ktf)
      endif
      enddo
!
!--- calculate moist static energy, heights, qes
!
      call cup_env(xz,xqes,xhe,xhes,xt,xq,po,z1, &
           psur,ierr,tcrit,-1,   &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
!
!--- environmental values on cloud levels
!
      call cup_env_clev(xt,xqes,xq,xhe,xhes,xz,po,xqes_cup,xq_cup, &
           xhe_cup,xhes_cup,xz_cup,po_cup,gamma_cup,xt_cup,psur,   &
           ierr,z1,          &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
!
!
!**************************** static control
!
!--- moist static energy inside cloud
!
      do k=kts,ktf
      do i=its,itf
         xhc(i,k)=0.
         xDBY(I,K)=0.
      enddo
      enddo
      do i=its,itf
        if(ierr(i).eq.0)then
         do k=1,kbcon(i)-1
            xhc(i,k)=xhkb(i)
         enddo
          k=kbcon(i)
          xhc(i,k)=xhkb(i)
          xDBY(I,Kbcon(i))=xHkb(I)-xHES_cup(I,K)
        endif !ierr
      enddo
!
!
      do i=its,itf
      if(ierr(i).eq.0)then
     do k=k22(i)  +1,ktop(i)
       xhc(i,k)=(xhc(i,k-1)*xzu(i,k-1)-.5*up_massdetro(i,k-1)*xhc(i,k-1)+ &
                         up_massentro(i,k-1)*xhe(i,k-1))   /            &
                         (xzu(i,k-1)-.5*up_massdetro(i,k-1)+up_massentro(i,k-1))
       xdby(i,k)=xhc(i,k)-xhes_cup(i,k)
      enddo
      do k=ktop(i)+1,ktf
           xHC(i,K)=xhes_cup(i,k)
           xDBY(I,K)=0.
      enddo
      endif
      enddo

!
!--- workfunctions for updraft
!
      call cup_up_aa0(xaa0,xz,xzu,xdby,GAMMA_CUP,xt_cup, &
           kbcon,ktop,ierr,           &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      do 200 nens=1,maxens
      do i=its,itf 
         if(ierr(i).eq.0)then
           xaa0_ens(i,nens)=xaa0(i)
           nall=(iens-1)*maxens3*maxens*maxens2 &
                +(iedt-1)*maxens*maxens3 &
                +(nens-1)*maxens3
           do k=kts,ktf
              if(k.le.ktop(i))then
                 do nens3=1,maxens3
                 if(nens3.eq.7)then
!--- b=0
                 pr_ens(i,nall+nens3)=pr_ens(i,nall+nens3)  &
                                    +pwo(i,k) 
!--- b=beta
                 else if(nens3.eq.8)then
                 pr_ens(i,nall+nens3)=pr_ens(i,nall+nens3)+ &
                                    pwo(i,k)
!--- b=beta/2
                 else if(nens3.eq.9)then
                 pr_ens(i,nall+nens3)=pr_ens(i,nall+nens3)  &
                                 +  pwo(i,k)
                 else
                 pr_ens(i,nall+nens3)=pr_ens(i,nall+nens3)+ &
                                    pwo(i,k) 
                 endif
                 enddo
              endif
           enddo
         if(pr_ens(i,nall+7).lt.1.e-6)then
            ierr(i)=18
            ierrc(i)="total normalized condensate too small"
            do nens3=1,maxens3
               pr_ens(i,nall+nens3)=0.
            enddo
         endif
         do nens3=1,maxens3
           if(pr_ens(i,nall+nens3).lt.1.e-4)then
            pr_ens(i,nall+nens3)=0.
           endif
         enddo
         endif
      enddo
 200  continue
      if(ifim.eq.1.and.ierr(1).ne.0)return
!
!--- LARGE SCALE FORCING
!
!
!------- CHECK wether aa0 should have been zero, assuming this 
!        ensemble is chosen
!
!
      do i=its,itf
         ierr2(i)=ierr(i)
         ierr3(i)=ierr(i)
         k22x(i)=k22(i)
      enddo
       if(maxens.gt.0)then
      CALL cup_MAXIMI(HEO_CUP,3,KBMAX,K22x,ierr, &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      call cup_kbcon(ierrc,cap_max_increment,2,k22x,kbconx,heo_cup, &
           heso_cup,hkbo,ierr2,kbmax,po_cup,cap_max, &
           ztexec,zqexec,use_excess,       &
           0,itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte, &
           z_cup,entr_rate,heo,imid)
      call cup_kbcon(ierrc,cap_max_increment,3,k22x,kbconx,heo_cup, &
           heso_cup,hkbo,ierr3,kbmax,po_cup,cap_max, &
           ztexec,zqexec,use_excess,       &
           0,itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte, &
           z_cup,entr_rate,heo,imid)
      endif
!
!--- calculate cloud base mass flux
!
      do i=its,itf
         mconv(i)=0.
         if(ierr(i).eq.0)then
            do k=k22(i),ktop(i)
                mconv(i)=mconv(i)-omeg(i,k)*(qo(i,k+1)-qo(i,k))
            enddo
         endif
      enddo
      call cup_forcing_ens_3d(closure_n,xland1,aa0,aa1,xaa0_ens,mbdt_ens,dtime,   &
           ierr,ierr2,ierr3,xf_ens,j,'deeps',axx,forcing,                 &
           maxens,iens,iedt,maxens2,maxens3,mconv,            &
           po_cup,ktop,omeg,zdo,k22,zuo,pr_ens,edto,kbcon,    &
           ensdim,ichoice,     &
           imid,ipr,jpr,itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte, &
	   dicycle,tau_ecmwf,aa1_bl,xf_dicycle)
!
      if(ifim.eq.1.and.ierr(1).ne.0)return
      do k=kts,ktf
      do i=its,itf
        if(ierr(i).eq.0)then
           dellat_ens(i,k,iedt)=dellat(i,k)
           dellaq_ens(i,k,iedt)=dellaq(i,k)
           dellaqc_ens(i,k,iedt)=dellaqc(i,k)
           pwo_ens(i,k,iedt)=pwo(i,k)+edto(i)*pwdo(i,k)
        else 
           dellat_ens(i,k,iedt)=0.
           dellaq_ens(i,k,iedt)=0.
           dellaqc_ens(i,k,iedt)=0.
           pwo_ens(i,k,iedt)=0.
        endif
      enddo
      enddo
 250  continue
!
!--- FEEDBACK
!
       if(imid.eq.1 .and. ichoice .le.2)then
         do i=its,itf
! boundary layer QE (from Saulo Freitas)
          xff_mid(i,1)=0.
          xff_mid(i,2)=0.
          if(ierr(i).eq.0)then
          blqe=0.
          trash=0.
          if(k22(i).lt.kpbl(i)+1)then
             do k=1,kpbl(i)
                blqe=blqe+100.*dhdt(i,k)*(po_cup(i,k)-po_cup(i,k+1))/g
             enddo
             trash=max((hco(i,kbcon(i))-heo_cup(i,kbcon(i))),1.e1)
             xff_mid(i,1)=max(0.,blqe/trash)
             xff_mid(i,1)=min(0.1,xff_mid(i,1))
           endif
           xff_mid(i,2)=.03*zws(i)
         endif
         enddo
       endif
       call cup_output_ens_3d(xff_mid,xf_ens,ierr,dellat_ens,dellaq_ens, &
            dellaqc_ens,outt,     &
            outq,outqc,zuo,sub_mas,pre,pwo_ens,xmb,ktop,      &
            j,'deep',maxens2,maxens,iens,ierr2,ierr3,         &
            pr_ens,maxens3,ensdim,                    &
            sig,closure_n,xland1,xmbm_in,xmbs_in,   &
            weight_GR,weight_W,weight_MC,weight_ST,weight_AS,training, &
            ichoice,imid,ipr,jpr,itf,jtf,ktf,                        &
            its,ite, jts,jte, kts,kte, &
	    dicycle,xf_dicycle )
      k=1
      do i=its,itf
          if(ierr(i).eq.0 .and.pre(i).gt.0.) then
             PRE(I)=MAX(PRE(I),0.)
             xmb_out(i)=xmb(i)
             do k=kts,ktop(i)
               outu(i,k)=dellu(i,k)*xmb(i)
               outv(i,k)=dellv(i,k)*xmb(i)
               sub_mas(i,k)=zuo(i,k)*xmb(i)
             enddo
          elseif(ierr(i).ne.0 .or. pre(i).eq.0.)then
             ktop(i)=0
             do k=kts,kte
               outt(i,k)=0.
               outq(i,k)=0.
               outqc(i,k)=0.
               outu(i,k)=0.
               outv(i,k)=0.
             enddo
          endif
      enddo
! rain evaporation as in SAS
!
      if(irainevap.eq.1)then
      do i = its,itf
       rntot(i) = 0.
       delqev(i) = 0.
       delq2(i) = 0.
       rn(i)    = 0.
       rntot(i)    = 0.
       rain=0.
       if(ierr(i).eq.0)then
         do k = ktop(i), 1, -1
              rain =  pwo(i,k) + edto(i) * pwdo(i,k)
              rntot(i) = rntot(i) + rain * xmb(i)* .001 * dtime
         enddo
       endif
      enddo
      do i = its,itf
         deltv(i) = 0.
         delq(i) = 0.
         qevap(i) = 0.
         flg(i) = .true.
         if(ierr(i).eq.0)then
         evef = edt(i) * evfact
         if(xland(i).gt.0.5 .and. xland(i).lt.1.5) evef=edt(i) * evfactl
!         write(32,*)'evef,edt = ',evef,edt(i)
!            write(32,*)'k,q1*1000.,t1,rn(i),qcond(i),qevap(i)'
         do k = ktop(i), 1, -1
              rain =  pwo(i,k) + edto(i) * pwdo(i,k)
              rn(i) = rn(i) + rain * xmb(i) * .001 * dtime
              if(flg(i))then
              q1=qo(i,k)+(outq(i,k))*dtime
              t1=tn(i,k)+(outt(i,k))*dtime
              qcond(i) = evef * (q1 - qeso(i,k))            &
     &                 / (1. + el2orc * qeso(i,k) / t1**2)
              dp = -100.*(p_cup(i,k+1)-p_cup(i,k))
              if(rn(i).gt.0. .and. qcond(i).lt.0.) then
                qevap(i) = -qcond(i) * (1.-exp(-.32*sqrt(dtime*rn(i))))
                qevap(i) = min(qevap(i), rn(i)*1000.*g/dp)
                delq2(i) = delqev(i) + .001 * qevap(i) * dp / g
              endif
              if(rn(i).gt.0..and.qcond(i).lt.0..and. &
     &           delq2(i).gt.rntot(i)) then
                qevap(i) = 1000.* g * (rntot(i) - delqev(i)) / dp
                flg(i) = .false.
              endif
              if(rn(i).gt.0..and.qevap(i).gt.0.) then
                outq(i,k) = outq(i,k) + qevap(i)/dtime
                outt(i,k) = outt(i,k) - elocp * qevap(i)/dtime
                rn(i) = max(0.,rn(i) - .001 * qevap(i) * dp / g)
                pre(i) = pre(i) - qevap(i) * dp /g/dtime
                PRE(I)=MAX(PRE(I),0.)
                deltv(i) = - elocp*qevap(i)/dtime
                delqev(i) = delqev(i) + .001*dp*qevap(i)/g
              endif
!            write(32,321)k,q1*1000.,t1,rn(i),qcond(i),qevap(i)
          endif
        enddo
!       pre(i)=1000.*rn(i)/dtime
      endif
      enddo
      endif
!321   format(1x,i2,2(1x,f6.1),3(1x,e13.4))
!      write(32,*)'pre new= ',pre(1)
!      i=1
!      do k=kts,ktop(i)
!      write(32,*)k,(outt(i,k))*86400.,(outq(i,k))*86400.*2.5e6/1004.
!      enddo


!
!---------------------------done------------------------------
!

   END SUBROUTINE CUP_gf


   SUBROUTINE cup_dd_edt(ierr,us,vs,z,ktop,kbcon,edt,p,pwav, &
              pw,ccn,pwev,edtmax,edtmin,maxens2,edtc,psum2,psumh, &
              ccnclean,rho,aeroevap,itf,jtf,ktf,j,ipr,jpr,          &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        j,ipr,jpr,aeroevap,itf,jtf,ktf,           &
        its,ite, jts,jte, kts,kte
     integer, intent (in   )              ::                           &
        maxens2
  !
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        rho,us,vs,z,p,pw
     real,    dimension (its:ite,1:maxens2)                            &
        ,intent (out  )                   ::                           &
        edtc
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        edt
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        pwav,pwev,ccn,psum2,psumh,edtmax,edtmin
     real                                                              &
        ,intent (in   )                   ::                           &
        ccnclean
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        ktop,kbcon
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
!
!  local variables in this routine
!

     integer i,k,kk
     real    einc,pef,pefb,prezk,zkbc
     real,    dimension (its:ite)         ::                           &
      vshear,sdp,vws
     real :: prop_c,pefc,aeroadd,alpha3,beta3,rhoc
     prop_c=8. !10.386
     alpha3 = 1.9
     beta3  = -1.13
     pefc=0.

!
!--- DETERMINE DOWNDRAFT STRENGTH IN TERMS OF WINDSHEAR
!
! */ calculate an average wind shear over the depth of the cloud
!
       do i=its,itf
        edt(i)=0.
        vws(i)=0.
        sdp(i)=0.
        vshear(i)=0.
       enddo
       do k=1,maxens2
       do i=its,itf
        edtc(i,k)=0.
       enddo
       enddo
       do kk = kts,ktf-1
         do 62 i=its,itf
          IF(ierr(i).ne.0)GO TO 62
          if (kk .le. min0(ktop(i),ktf) .and. kk .ge. kbcon(i)) then
             vws(i) = vws(i)+ &
              (abs((us(i,kk+1)-us(i,kk))/(z(i,kk+1)-z(i,kk))) &
          +   abs((vs(i,kk+1)-vs(i,kk))/(z(i,kk+1)-z(i,kk)))) * &
              (p(i,kk) - p(i,kk+1))
            sdp(i) = sdp(i) + p(i,kk) - p(i,kk+1)
          endif
          if (kk .eq. ktf-1)vshear(i) = 1.e3 * vws(i) / sdp(i)
   62   continue
       end do
      do i=its,itf
         IF(ierr(i).eq.0)then
            pef=(1.591-.639*VSHEAR(I)+.0953*(VSHEAR(I)**2) &
               -.00496*(VSHEAR(I)**3))
            if(pef.gt.0.9)pef=0.9
            if(pef.lt.0.1)pef=0.1
!
!--- cloud base precip efficiency
!
            zkbc=z(i,kbcon(i))*3.281e-3
            prezk=.02
            if(zkbc.gt.3.)then
               prezk=.96729352+zkbc*(-.70034167+zkbc*(.162179896+zkbc &
               *(- 1.2569798E-2+zkbc*(4.2772E-4-zkbc*5.44E-6))))
            endif
            if(zkbc.gt.25)then
               prezk=2.4
            endif
            pefb=1./(1.+prezk)
!           write(11,*)pefb,prezk,zkbc
            if(pefb.gt.0.9)pefb=0.9
            if(pefb.lt.0.1)pefb=0.1
            EDT(I)=1.-.5*(pefb+pef)
            if(aeroevap.gt.1)then
               aeroadd=(ccnclean**beta3)*((psumh(i))**(alpha3-1)) !*1.e6
!              if(i.eq.ipr.and.j.eq.jpr)write(0,*)'edt',ccnclean,psumh(i),aeroadd
!              prop_c=.9/aeroadd
               prop_c=.5*(pefb+pef)/aeroadd
               aeroadd=(ccn(i)**beta3)*((psum2(i))**(alpha3-1)) !*1.e6
!              if(i.eq.ipr.and.j.eq.jpr)write(0,*)'edt',ccn(i),psum2(i),aeroadd,prop_c
               aeroadd=prop_c*aeroadd
               pefc=aeroadd
               if(pefc.gt.0.9)pefc=0.9
               if(pefc.lt.0.1)pefc=0.1
               EDT(I)=1.-pefc
               if(aeroevap.eq.2)EDT(I)=1.-.25*(pefb+pef+2.*pefc)
            endif


!--- edt here is 1-precipeff!
            einc=.2*edt(i)
            do k=1,maxens2
                edtc(i,k)=edt(i)+float(k-2)*einc
            enddo
         endif
      enddo
      do i=its,itf
         IF(ierr(i).eq.0)then
            do k=1,maxens2
               EDTC(I,K)=-EDTC(I,K)*pwav(I)/PWEV(I)
               IF(EDTC(I,K).GT.edtmax(i))EDTC(I,K)=edtmax(i)
               IF(EDTC(I,K).LT.edtmin(i))EDTC(I,K)=edtmin(i)
            enddo
         endif
      enddo

   END SUBROUTINE cup_dd_edt


   SUBROUTINE cup_dd_moisture_new(ierrc,zd,hcd,hes_cup,qcd,qes_cup,    &
              pwd,q_cup,z_cup,dd_massentr,dd_massdetr,jmin,ierr,            &
              gamma_cup,pwev,bu,qrcd,                        &
              q,he,t_cup,iloop,           &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  itf,jtf,ktf,           &
                                  its,ite, jts,jte, kts,kte
  ! cdd= detrainment function 
  ! q = environmental q on model levels
  ! q_cup = environmental q on model cloud levels
  ! qes_cup = saturation q on model cloud levels
  ! hes_cup = saturation h on model cloud levels
  ! hcd = h in model cloud
  ! bu = buoancy term
  ! zd = normalized downdraft mass flux
  ! gamma_cup = gamma on model cloud levels
  ! mentr_rate = entrainment rate
  ! qcd = cloud q (including liquid water) after entrainment
  ! qrch = saturation q in cloud
  ! pwd = evaporate at that level
  ! pwev = total normalized integrated evaoprate (I2)
  ! entr= entrainment rate 
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        zd,t_cup,hes_cup,hcd,qes_cup,q_cup,z_cup,                      &
        dd_massentr,dd_massdetr,gamma_cup,q,he 
     integer                                                           &
        ,intent (in   )                   ::                           &
        iloop
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        jmin
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        qcd,qrcd,pwd
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        pwev,bu
     character*50 :: ierrc(its:ite)
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k,ki
     real                                 ::                           &
        dh,dz,dqeva

      do i=its,itf
         bu(i)=0.
         pwev(i)=0.
      enddo
      do k=kts,ktf
      do i=its,itf
         qcd(i,k)=0.
         qrcd(i,k)=0.
         pwd(i,k)=0.
      enddo
      enddo
!
!
!
      do 100 i=its,itf
      IF(ierr(I).eq.0)then
      k=jmin(i)
      DZ=Z_cup(i,K+1)-Z_cup(i,K)
      qcd(i,k)=q_cup(i,k)
      DH=HCD(I,k)-HES_cup(I,K)
      if(dh.lt.0)then
        QRCD(I,K)=(qes_cup(i,k)+(1./XLV)*(GAMMA_cup(i,k) &
                  /(1.+GAMMA_cup(i,k)))*DH)
        else
          qrcd(i,k)=qes_cup(i,k)
        endif
      pwd(i,jmin(i))=zd(i,jmin(i))*min(0.,qcd(i,k)-qrcd(i,k))
      qcd(i,k)=qrcd(i,k)
      pwev(i)=pwev(i)+pwd(i,jmin(i)) ! *dz
!
      bu(i)=dz*dh
      do ki=jmin(i)-1,1,-1
         DZ=Z_cup(i,Ki+1)-Z_cup(i,Ki)
!        QCD(i,Ki)=(qCD(i,Ki+1)*(1.-.5*CDD(i,Ki+1)*DZ) &
!                 +entr*DZ*q(i,Ki) &
!                )/(1.+entr*DZ-.5*CDD(i,Ki+1)*DZ)
!        dz=qcd(i,ki)
         qcd(i,ki)=(qcd(i,ki+1)*zd(i,ki+1)                          &
                  -.5*dd_massdetr(i,ki)*qcd(i,ki+1)+ &
                  dd_massentr(i,ki)*q(i,ki))   /            &
                  (zd(i,ki+1)-.5*dd_massdetr(i,ki)+dd_massentr(i,ki))
!        write(0,*)'qcd in dd_moi = ',qcd(i,ki)

!
!--- to be negatively buoyant, hcd should be smaller than hes!
!--- ideally, dh should be negative till dd hits ground, but that is not always
!--- the case
!
         DH=HCD(I,ki)-HES_cup(I,Ki)
         bu(i)=bu(i)+dz*dh
         QRCD(I,Ki)=qes_cup(i,ki)+(1./XLV)*(GAMMA_cup(i,ki) &
                  /(1.+GAMMA_cup(i,ki)))*DH
         dqeva=qcd(i,ki)-qrcd(i,ki)
         if(dqeva.gt.0.)then
          dqeva=0.
          qrcd(i,ki)=qcd(i,ki)
         endif
         pwd(i,ki)=zd(i,ki)*dqeva
         qcd(i,ki)=qrcd(i,ki)
         pwev(i)=pwev(i)+pwd(i,ki) ! *dz
!        if(iloop.eq.1.and.i.eq.102.and.j.eq.62)then
!         print *,'in cup_dd_moi ', hcd(i,ki),HES_cup(I,Ki),dh,dqeva
!        endif
      enddo
!
!--- end loop over i
       if(pwev(I).eq.0.and.iloop.eq.1)then
!        print *,'problem with buoy in cup_dd_moisture',i
         ierr(i)=7
         ierrc(i)="problem with buoy in cup_dd_moisture"
       endif
       if(BU(I).GE.0.and.iloop.eq.1)then
!        print *,'problem with buoy in cup_dd_moisture',i
         ierr(i)=7
         ierrc(i)="problem2 with buoy in cup_dd_moisture"
       endif
      endif
100    continue

   END SUBROUTINE cup_dd_moisture_new

   SUBROUTINE cup_env(z,qes,he,hes,t,q,p,z1,                 &
              psur,ierr,tcrit,itest,                   &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,jtf,ktf,           &
        its,ite, jts,jte, kts,kte
  !
  ! ierr error value, maybe modified in this routine
  ! q           = environmental mixing ratio
  ! qes         = environmental saturation mixing ratio
  ! t           = environmental temp
  ! tv          = environmental virtual temp
  ! p           = environmental pressure
  ! z           = environmental heights
  ! he          = environmental moist static energy
  ! hes         = environmental saturation moist static energy
  ! psur        = surface pressure
  ! z1          = terrain elevation
  ! 
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        p,t,q
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        he,hes,qes
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout)                   ::                           &
        z
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        psur,z1
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     integer                                                           &
        ,intent (in   )                   ::                           &
        itest
!
!  local variables in this routine
!

     integer                              ::                           &
       i,k,iph
      real, dimension (1:2) :: AE,BE,HT
      real, dimension (its:ite,kts:kte) :: tv
      real :: tcrit,e,tvbar
!      real, external :: satvap
!      real :: satvap


      HT(1)=XLV/CP
      HT(2)=2.834E6/CP
      BE(1)=.622*HT(1)/.286
      AE(1)=BE(1)/273.+ALOG(610.71)
      BE(2)=.622*HT(2)/.286
      AE(2)=BE(2)/273.+ALOG(610.71)
!      print *, 'TCRIT = ', tcrit,its,ite
      DO k=kts,ktf
      do i=its,itf
        if(ierr(i).eq.0)then
!Csgb - IPH is for phase, dependent on TCRIT (water or ice)
        IPH=1
        IF(T(I,K).LE.TCRIT)IPH=2
!       print *, 'AE(IPH),BE(IPH) = ',AE(IPH),BE(IPH),AE(IPH)-BE(IPH),T(i,k),i,k
!       E=EXP(AE(IPH)-BE(IPH)/T(I,K))
!       print *, 'P, E = ', P(I,K), E
!       QES(I,K)=.622*E/(100.*P(I,K)-E)
        e=satvap(t(i,k))
        qes(i,k)=0.622*e/max(1.e-8,(p(i,k)-e))
        IF(QES(I,K).LE.1.E-16)QES(I,K)=1.E-16
        IF(QES(I,K).LT.Q(I,K))QES(I,K)=Q(I,K)
!       IF(Q(I,K).GT.QES(I,K))Q(I,K)=QES(I,K)
        TV(I,K)=T(I,K)+.608*Q(I,K)*T(I,K)
        endif
      enddo
      enddo
!
!--- z's are calculated with changed h's and q's and t's
!--- if itest=2
!
      if(itest.eq.1 .or. itest.eq.0)then
         do i=its,itf
           if(ierr(i).eq.0)then
             Z(I,1)=max(0.,Z1(I))-(ALOG(P(I,1))- &
                 ALOG(PSUR(I)))*287.*TV(I,1)/9.81
           endif
         enddo

! --- calculate heights
         DO K=kts+1,ktf
         do i=its,itf
           if(ierr(i).eq.0)then
              TVBAR=.5*TV(I,K)+.5*TV(I,K-1)
              Z(I,K)=Z(I,K-1)-(ALOG(P(I,K))- &
               ALOG(P(I,K-1)))*287.*TVBAR/9.81
           endif
         enddo
         enddo
      else if(itest.eq.2)then
         do k=kts,ktf
         do i=its,itf
           if(ierr(i).eq.0)then
             z(i,k)=(he(i,k)-1004.*t(i,k)-2.5e6*q(i,k))/9.81
             z(i,k)=max(1.e-3,z(i,k))
           endif
         enddo
         enddo
      else if(itest.eq.-1)then
      endif
!
!--- calculate moist static energy - HE
!    saturated moist static energy - HES
!
       DO k=kts,ktf
       do i=its,itf
         if(ierr(i).eq.0)then
         if(itest.le.0)HE(I,K)=9.81*Z(I,K)+1004.*T(I,K)+2.5E06*Q(I,K)
         HES(I,K)=9.81*Z(I,K)+1004.*T(I,K)+2.5E06*QES(I,K)
         IF(HE(I,K).GE.HES(I,K))HE(I,K)=HES(I,K)
         endif
      enddo
      enddo

   END SUBROUTINE cup_env


   SUBROUTINE cup_env_clev(t,qes,q,he,hes,z,p,qes_cup,q_cup,   &
              he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur, &
              ierr,z1,                                &
              itf,jtf,ktf,                       &
              its,ite, jts,jte, kts,kte                       )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,jtf,ktf,           &
        its,ite, jts,jte, kts,kte
  !
  ! ierr error value, maybe modified in this routine
  ! q           = environmental mixing ratio
  ! q_cup       = environmental mixing ratio on cloud levels
  ! qes         = environmental saturation mixing ratio
  ! qes_cup     = environmental saturation mixing ratio on cloud levels
  ! t           = environmental temp
  ! t_cup       = environmental temp on cloud levels
  ! p           = environmental pressure
  ! p_cup       = environmental pressure on cloud levels
  ! z           = environmental heights
  ! z_cup       = environmental heights on cloud levels
  ! he          = environmental moist static energy
  ! he_cup      = environmental moist static energy on cloud levels
  ! hes         = environmental saturation moist static energy
  ! hes_cup     = environmental saturation moist static energy on cloud levels
  ! gamma_cup   = gamma on cloud levels
  ! psur        = surface pressure
  ! z1          = terrain elevation
  ! 
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        qes,q,he,hes,z,p,t
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        qes_cup,q_cup,he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        psur,z1
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
!
!  local variables in this routine
!

     integer                              ::                           &
       i,k


      do k=kts,ktf
      do i=its,itf
        qes_cup(i,k)=0.
        q_cup(i,k)=0.
        hes_cup(i,k)=0.
        he_cup(i,k)=0.
        z_cup(i,k)=0.
        p_cup(i,k)=0.
        t_cup(i,k)=0.
        gamma_cup(i,k)=0.
      enddo
      enddo
      do k=kts+1,ktf
      do i=its,itf
        if(ierr(i).eq.0)then
        qes_cup(i,k)=.5*(qes(i,k-1)+qes(i,k))
        q_cup(i,k)=.5*(q(i,k-1)+q(i,k))
        hes_cup(i,k)=.5*(hes(i,k-1)+hes(i,k))
        he_cup(i,k)=.5*(he(i,k-1)+he(i,k))
        if(he_cup(i,k).gt.hes_cup(i,k))he_cup(i,k)=hes_cup(i,k)
        z_cup(i,k)=.5*(z(i,k-1)+z(i,k))
        p_cup(i,k)=.5*(p(i,k-1)+p(i,k))
        t_cup(i,k)=.5*(t(i,k-1)+t(i,k))
        gamma_cup(i,k)=(xlv/cp)*(xlv/(r_v*t_cup(i,k) &
                       *t_cup(i,k)))*qes_cup(i,k)
        endif
      enddo
      enddo
      do i=its,itf
        if(ierr(i).eq.0)then
        qes_cup(i,1)=qes(i,1)
        q_cup(i,1)=q(i,1)
!       hes_cup(i,1)=hes(i,1)
!       he_cup(i,1)=he(i,1)
        hes_cup(i,1)=9.81*z1(i)+1004.*t(i,1)+2.5e6*qes(i,1)
        he_cup(i,1)=9.81*z1(i)+1004.*t(i,1)+2.5e6*q(i,1)
        z_cup(i,1)=.5*(z(i,1)+z1(i))
        p_cup(i,1)=.5*(p(i,1)+psur(i))
        z_cup(i,1)=z1(i)
        p_cup(i,1)=psur(i)
        t_cup(i,1)=t(i,1)
        gamma_cup(i,1)=xlv/cp*(xlv/(r_v*t_cup(i,1) &
                       *t_cup(i,1)))*qes_cup(i,1)
        endif
      enddo

   END SUBROUTINE cup_env_clev

   SUBROUTINE cup_forcing_ens_3d(closure_n,xland,aa0,aa1,xaa0,mbdt,dtime,ierr,ierr2,ierr3,&
              xf_ens,j,name,axx,forcing,maxens,iens,iedt,maxens2,maxens3,mconv,    &
              p_cup,ktop,omeg,zd,k22,zu,pr_ens,edt,kbcon,      &
              ensdim,icoic,            &
              imid,ipr,jpr,itf,jtf,ktf,               &
              its,ite, jts,jte, kts,kte, &
	      dicycle,tau_ecmwf,aa1_bl,xf_dicycle  )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        imid,ipr,jpr,itf,jtf,ktf,           &
        its,ite, jts,jte, kts,kte
     integer, intent (in   )              ::                           &
        j,ensdim,maxens,iens,iedt,maxens2,maxens3
  !
  ! ierr error value, maybe modified in this routine
  ! pr_ens = precipitation ensemble
  ! xf_ens = mass flux ensembles
  ! massfln = downdraft mass flux ensembles used in next timestep
  ! omeg = omega from large scale model
  ! mconv = moisture convergence from large scale model
  ! zd      = downdraft normalized mass flux
  ! zu      = updraft normalized mass flux
  ! aa0     = cloud work function without forcing effects
  ! aa1     = cloud work function with forcing effects
  ! xaa0    = cloud work function with cloud effects (ensemble dependent)
  ! edt     = epsilon
  ! dir     = "storm motion"
  ! mbdt    = arbitrary numerical parameter
  ! dtime   = dt over which forcing is applied
  ! iact_gr_old = flag to tell where convection was active
  ! kbcon       = LFC of parcel from k22
  ! k22         = updraft originating level
  ! icoic       = flag if only want one closure (usually set to zero!)
  ! name        = deep or shallow convection flag
  !
     real,    dimension (its:ite,1:ensdim)                     &
        ,intent (inout)                   ::                           &
        pr_ens
     real,    dimension (its:ite,1:ensdim)                     &
        ,intent (inout  )                   ::                           &
        xf_ens
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        zd,zu,p_cup
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        omeg
     real,    dimension (its:ite,1:maxens)                             &
        ,intent (in   )                   ::                           &
        xaa0
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        aa1,edt
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        mconv,axx
     real,    dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        aa0,closure_n
     real,    dimension (1:maxens)                                     &
        ,intent (in   )                   ::                           &
        mbdt
     real                                                              &
        ,intent (in   )                   ::                           &
        dtime
     integer, dimension (its:ite)                                      &
        ,intent (inout   )                   ::                           &
        k22,kbcon,ktop
     integer, dimension (its:ite)                                      &
        ,intent (in      )                   ::                           &
        xland
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr,ierr2,ierr3
     integer                                                           &
        ,intent (in   )                   ::                           &
        icoic
      character *(*), intent (in)         ::                           &
       name
!-srf begin
      integer, intent(IN) :: DICYCLE
      real,    intent(IN)   , dimension (its:ite) :: aa1_bl,tau_ecmwf
      real,    intent(INOUT), dimension (its:ite) :: xf_dicycle
      !- local var
      real  :: xff_dicycle
!-srf end
!
!  local variables in this routine
!

     real,    dimension (1:maxens3)       ::                           &
       xff_ens3
     real,    dimension (1:maxens)        ::                           &
       xk
     integer                              ::                           &
       i,k,nall,n,ne,nens,nens3,iresult,iresultd,iresulte,mkxcrt,kclim
     parameter (mkxcrt=15)
     real                                 ::                           &
       a1,massfld,a_ave,xff0,xff00,xxx,xomg,aclim1,aclim2,aclim3,aclim4
     real,    dimension(1:mkxcrt)         ::                           &
       pcrit,acrit,acritt

     integer :: nall2,ixxx,irandom
     integer,  dimension (8) :: seed
     real, dimension (its:ite) :: ens_adj
     real, dimension (its:ite,10) :: forcing


      DATA PCRIT/850.,800.,750.,700.,650.,600.,550.,500.,450.,400.,    &
                 350.,300.,250.,200.,150./
      DATA ACRIT/.0633,.0445,.0553,.0664,.075,.1082,.1521,.2216,       &
                 .3151,.3677,.41,.5255,.7663,1.1686,1.6851/
!  GDAS DERIVED ACRIT
      DATA ACRITT/.203,.515,.521,.566,.625,.665,.659,.688,             &
                  .743,.813,.886,.947,1.138,1.377,1.896/

!
       ens_adj=1.
       seed=0
       do i=its,itf
        if(ierr(i).eq.0)then
          seed(1)=int(aa0(i))
          seed(2)=int(aa1(i))
          exit
        endif
       enddo

       nens=0
       irandom=0

!--- LARGE SCALE FORCING
!
       DO 100 i=its,itf
          if(name.eq.'deeps'.and.ierr(i).gt.995)then
           aa0(i)=0.
           ierr(i)=0
          endif
          IF(ierr(i).eq.0)then
          ens_adj(i)=1.
           if(ierr2(i).gt.0.and.ierr3(i).eq.0)ens_adj(i)=0.666 ! 2./3.
           if(ierr2(i).gt.0.and.ierr3(i).gt.0)ens_adj(i)=0.333
!
!---
!
             if(name.eq.'deeps')then
!
                a_ave=0.
                  a_ave=axx(i)
!               if(i.eq.ipr.and.j.eq.jpr)write(0,*)'in forcing, a_ave,axx(i,ne)
!               = ',a_ave,axx(i,ne),maxens,xland(i)
                a_ave=max(0.,a_ave)
                a_ave=min(a_ave,aa1(i))
                a_ave=max(0.,a_ave)
                do ne=1,16
                  xff_ens3(ne)=0.
                enddo
                xff0= (AA1(I)-AA0(I))/DTIME
                xff_ens3(1)=max(0.,(AA1(I)-AA0(I))/dtime)
                xff_ens3(2)=max(0.,(AA1(I)-AA0(I))/dtime)
                forcing(i,1)=xff_ens3(2)
!               xff_ens3(2)=max(0.,(a_ave-AA0(I))/dtime)

!               if(i.eq.ipr.and.j.eq.jpr)write(0,*)AA1(I),AA0(I),xff_ens3(1),xff_ens3(2),dtime
                   xff_ens3(3)=max(0.,(AA1(I)-AA0(I))/dtime)
                   xff_ens3(16)=max(0.,(AA1(I)-AA0(I))/dtime)
!   
!--- omeg is in bar/s, mconv done with omeg in Pa/s
!     more like Brown (1979), or Frank-Cohen (199?)
!
! average aaround kbcon
!
                   xomg=0.
		if(zu(i,kbcon(i)-1).gt.0.) then
                   do k=kbcon(i)-1,kbcon(i)+1
                     xomg=xomg-omeg(i,k)/9.81/zu(i,k)
                   enddo
                end if
                   xff_ens3(4)=xomg/3.

!
! max below kbcon
                   xff_ens3(6)=-omeg(i,k22(i))/9.81
                   do k=k22(i),kbcon(i)
                     xomg=-omeg(i,k)/9.81
                     if(xomg.gt.xff_ens3(6))xff_ens3(6)=xomg
                   enddo
!
                   xff_ens3(6)=betajb*xff_ens3(6)/zu(i,kbcon(i))
                   xff_ens3(4)=betajb*xff_ens3(4)
                   xff_ens3(5)=betajb*xff_ens3(4)
                   if(xff_ens3(4).lt.0.)xff_ens3(4)=0.
                   if(xff_ens3(5).lt.0.)xff_ens3(5)=0.
                   if(xff_ens3(6).lt.0.)xff_ens3(6)=0.
                   xff_ens3(14)=betajb*xff_ens3(4)
                forcing(i,2)=xff_ens3(6)
!               if(i.eq.ipr.and.j.eq.jpr)write(0,*)xff_ens3(4),xff_ens3(5)
!
!--- more like Krishnamurti et al.; pick max and average values
!
                xff_ens3(7)=mconv(i)
                xff_ens3(8)=mconv(i)
                xff_ens3(9)=mconv(i)
                forcing(i,3)=xff_ens3(8)
!               if(i.eq.ipr.and.j.eq.jpr)write(0,*)xff_ens3(7),xff_ens3(8)
!
                xff_ens3(15)=mconv(i)
!
!
!- more like Bechtold et al. (JAS 2014)
                xff_ens3(10)=AA1(i)/tau_ecmwf(i)
                xff_ens3(11)=AA1(I)/tau_ecmwf(i)
                xff_ens3(12)=AA1(I)/tau_ecmwf(i)
                xff_ens3(13)=(AA1(i))/tau_ecmwf(i) !(60.*15.) !tau_ecmwf(i)
                forcing(i,4)=xff_ens3(10)
                xff_dicycle = max(0.,AA1_BL(i)/tau_ecmwf(i)) !(60.*30.) !tau_ecmwf(i)
!-srf end  
!
!gtest
                if(icoic.eq.0)then
                if(xff0.lt.0.)then
                     xff_ens3(1)=0.
                     xff_ens3(2)=0.
                     xff_ens3(3)=0.
                     xff_ens3(10)=0.
                     xff_ens3(11)=0.
                     xff_ens3(12)=0.
		     xff_ens3(13)= 0.
		     xff_ens3(16)= 0.
		     xff_dicycle = 0.
                endif
                  if(xff0.lt.0 .and. xland(i).lt.0.1 .and. imid.eq.0)then
                     xff_ens3(:)=0.
                  endif
                endif

!                  if(i.eq.ipr.and.j.eq.jpr)write(0,*)'xff_ens =
!                  ',i,ipr,jpr,xff_ens3


                do nens=1,maxens
                   XK(nens)=(XAA0(I,nens)-AA1(I))/MBDT(1)
                   forcing(i,5)=xk(nens)
!                  if(i.eq.ipr.and.j.eq.jpr)write(0,*)'xks =
!                  ',xk(nens),XAA0(I,nens),AA1(I),mbdt
                   if(xk(nens).le.0.and.xk(nens).gt.-.01*mbdt(1)) &
                           xk(nens)=-.01*mbdt(1)
                   if(xk(nens).gt.0.and.xk(nens).lt.1.e-2) &
                           xk(nens)=1.e-2
!                forcing(i,5)=xk(nens)
                enddo
!
!--- add up all ensembles
!
                do 350 ne=1,maxens
!
!--- for every xk, we have maxens3 xffs
!--- iens is from outermost ensemble (most expensive!
!
!--- iedt (maxens2 belongs to it)
!--- is from second, next outermost, not so expensive
!
!--- so, for every outermost loop, we have maxens*maxens2*3
!--- ensembles!!! nall would be 0, if everything is on first
!--- loop index, then ne would start counting, then iedt, then iens....
!
                   iresult=0
                   iresultd=0
                   iresulte=0
                   nall=(iens-1)*maxens3*maxens*maxens2 &
                        +(iedt-1)*maxens*maxens3 &
                        +(ne-1)*maxens3
!                 if(i.eq.ipr.and.j.eq.jpr)write(0,*)'maxens',ne,nall,iens,maxens3,maxens,maxens2,iedt
!
! over water, enfor!e small cap for some of the closures
!
                if(maxens.gt.0 .and. xland(i).lt.0.1)then
                 if(ierr2(i).gt.0.or.ierr3(i).gt.0)then
                      xff_ens3(1) =ens_adj(i)*xff_ens3(1)
                      xff_ens3(2) =ens_adj(i)*xff_ens3(2)
                      xff_ens3(3) =ens_adj(i)*xff_ens3(3)
                      xff_ens3(4) =ens_adj(i)*xff_ens3(4)
                      xff_ens3(5) =ens_adj(i)*xff_ens3(5)
                      xff_ens3(6) =ens_adj(i)*xff_ens3(6)
                      xff_ens3(7) =ens_adj(i)*xff_ens3(7)
                      xff_ens3(8) =ens_adj(i)*xff_ens3(8)
                      xff_ens3(9) =ens_adj(i)*xff_ens3(9)
                      xff_ens3(10) =ens_adj(i)*xff_ens3(10)
                      xff_ens3(11) =ens_adj(i)*xff_ens3(11)
                      xff_ens3(12) =ens_adj(i)*xff_ens3(12)
                      xff_ens3(13) =ens_adj(i)*xff_ens3(13)
                      xff_ens3(14) =ens_adj(i)*xff_ens3(14)
                      xff_ens3(15) =ens_adj(i)*xff_ens3(15)
                      xff_ens3(16) =ens_adj(i)*xff_ens3(16)
		      !srf
         	      xff_dicycle = ens_adj(i)*xff_dicycle
                      !srf end
!                      xff_ens3(7) =0.
!                      xff_ens3(8) =0.
!                      xff_ens3(9) =0.
                 endif
                endif
!
! end water treatment
!
!
!--- check for upwind convection
!                  iresult=0
                   massfld=0.

                   IF(XK(ne).lt.0.and.xff0.gt.0.)iresultd=1
                   iresulte=max(iresult,iresultd)
                   iresulte=1
                   if(iresulte.eq.1)then
!
!--- special treatment for stability closures
!
!                      if(i.eq.ipr.and.j.eq.jpr)write(0,*)'xffs =
!                      ',xff_ens3(1:16)

!                     if(xff0.gt.0.)then
                          if(XK(ne).lt.0.)then
                         if(xff_ens3(1).gt.0)xf_ens(i,nall+1)=max(0.,-xff_ens3(1)/xk(ne))
                         if(xff_ens3(2).gt.0)xf_ens(i,nall+2)=max(0.,-xff_ens3(2)/xk(ne))
                         if(xff_ens3(3).gt.0)xf_ens(i,nall+3)=max(0.,-xff_ens3(3)/xk(ne))
                         if(xff_ens3(16).gt.0)xf_ens(i,nall+16)=max(0.,-xff_ens3(16)/xk(ne))
                         forcing(i,6)=xf_ens(i,nall+2)
                      else
                         xff_ens3(1)=0
                         xff_ens3(2)=0
                         xff_ens3(3)=0
                         xff_ens3(16)=0
                      endif
!
!--- if iresult.eq.1, following independent of xff0
!
                         xf_ens(i,nall+4)=max(0.,xff_ens3(4))
                         xf_ens(i,nall+5)=max(0.,xff_ens3(5))
                         xf_ens(i,nall+6)=max(0.,xff_ens3(6))
                         xf_ens(i,nall+14)=max(0.,xff_ens3(14))
                         a1=max(1.e-3,pr_ens(i,nall+7))
                         xf_ens(i,nall+7)=max(0.,xff_ens3(7)/a1)
!                      if(i.eq.ipr.and.j.eq.jpr)write(0,*)'a1 =
!                      ',xff_ens3(7),a1,xf_ens(i,nall+7)
                         a1=max(1.e-3,pr_ens(i,nall+8))
                         xf_ens(i,nall+8)=max(0.,xff_ens3(8)/a1)
                         forcing(i,7)=xf_ens(i,nall+8)
                         a1=max(1.e-3,pr_ens(i,nall+9))
                         xf_ens(i,nall+9)=max(0.,xff_ens3(9)/a1)
                         a1=max(1.e-3,pr_ens(i,nall+15))
                         xf_ens(i,nall+15)=max(0.,xff_ens3(15)/a1)
                         if(XK(ne).lt.0.)then
                            xf_ens(i,nall+10)=max(0.,-xff_ens3(10)/xk(ne))
                            xf_ens(i,nall+11)=max(0.,-xff_ens3(11)/xk(ne))
                            xf_ens(i,nall+12)=max(0.,-xff_ens3(12)/xk(ne))
                            xf_ens(i,nall+13)=max(0.,-xff_ens3(13)/xk(ne))
                            forcing(i,8)=xf_ens(i,nall+11)
                         endif
!srf-begin
                          if(XK(ne).lt.0.)then
			    xf_dicycle(i)      =  max(0.,-xff_dicycle /xk(ne))
                         forcing(i,9)=xf_dicycle(i)
                          else
                            xf_dicycle(i)      = 0.
                          endif
!srf-end
                      if(icoic.ge.1)then
                      closure_n(i)=0.
                      xf_ens(i,nall+1)=xf_ens(i,nall+icoic)
                      xf_ens(i,nall+2)=xf_ens(i,nall+icoic)
                      xf_ens(i,nall+3)=xf_ens(i,nall+icoic)
                      xf_ens(i,nall+4)=xf_ens(i,nall+icoic)
                      xf_ens(i,nall+5)=xf_ens(i,nall+icoic)
                      xf_ens(i,nall+6)=xf_ens(i,nall+icoic)
                      xf_ens(i,nall+7)=xf_ens(i,nall+icoic)
                      xf_ens(i,nall+8)=xf_ens(i,nall+icoic)
                      xf_ens(i,nall+9)=xf_ens(i,nall+icoic)
                      xf_ens(i,nall+10)=xf_ens(i,nall+icoic)
                      xf_ens(i,nall+11)=xf_ens(i,nall+icoic)
                      xf_ens(i,nall+12)=xf_ens(i,nall+icoic)
                      xf_ens(i,nall+13)=xf_ens(i,nall+icoic)
                      xf_ens(i,nall+14)=xf_ens(i,nall+icoic)
                      xf_ens(i,nall+15)=xf_ens(i,nall+icoic)
                      xf_ens(i,nall+16)=xf_ens(i,nall+icoic)
                      endif
!
! 16 is a randon pick from the oher 15
!
                if(irandom.eq.1)then
                   call random_number (xxx)
                   ixxx=min(15,max(1,int(15.*xxx+1.e-8)))
                   xf_ens(i,nall+16)=xf_ens(i,nall+ixxx)
!               else
!                  xf_ens(i,nall+16)=xf_ens(i,nall+1)
                endif
!
!
!--- do some more on the caps!!! ne=1 for 175, ne=2 for 100,....
!
!     do not care for caps here for closure groups 1 and 5,
!     they are fine, do not turn them off here
!
!!!!    NOT USED FOR "NORMAL" APPLICATION (maxens=1)
!
                if(maxens.gt.1)then
                if(ne.eq.2.and.ierr2(i).gt.0)then
                      xf_ens(i,nall+1) =0.
                      xf_ens(i,nall+2) =0.
                      xf_ens(i,nall+3) =0.
                      xf_ens(i,nall+4) =0.
                      xf_ens(i,nall+5) =0.
                      xf_ens(i,nall+6) =0.
                      xf_ens(i,nall+7) =0.
                      xf_ens(i,nall+8) =0.
                      xf_ens(i,nall+9) =0.
                      xf_ens(i,nall+10)=0.
                      xf_ens(i,nall+11)=0.
                      xf_ens(i,nall+12)=0.
                      xf_ens(i,nall+13)=0.
                      xf_ens(i,nall+14)=0.
                      xf_ens(i,nall+15)=0.
                      xf_ens(i,nall+16)=0.
                endif
                if(ne.eq.3.and.ierr3(i).gt.0)then
                      xf_ens(i,nall+1) =0.
                      xf_ens(i,nall+2) =0.
                      xf_ens(i,nall+3) =0.
                      xf_ens(i,nall+4) =0.
                      xf_ens(i,nall+5) =0.
                      xf_ens(i,nall+6) =0.
                      xf_ens(i,nall+7) =0.
                      xf_ens(i,nall+8) =0.
                      xf_ens(i,nall+9) =0.
                      xf_ens(i,nall+10)=0.
                      xf_ens(i,nall+11)=0.
                      xf_ens(i,nall+12)=0.
                      xf_ens(i,nall+13)=0.
                      xf_ens(i,nall+14)=0.
                      xf_ens(i,nall+15)=0.
                      xf_ens(i,nall+16)=0.
                endif
                endif

                   endif
 350            continue
                if(maxens.gt.1)then
! ne=1, cap=175
!
                   nall=(iens-1)*maxens3*maxens*maxens2 &
                        +(iedt-1)*maxens*maxens3
! ne=2, cap=100
!
                   nall2=(iens-1)*maxens3*maxens*maxens2 &
                        +(iedt-1)*maxens*maxens3 &
                        +(2-1)*maxens3
                      xf_ens(i,nall+4) = xf_ens(i,nall2+4)
                      xf_ens(i,nall+5) =xf_ens(i,nall2+5)
                      xf_ens(i,nall+6) =xf_ens(i,nall2+6)
                      xf_ens(i,nall+14) =xf_ens(i,nall2+14)
                      xf_ens(i,nall+7) =xf_ens(i,nall2+7)
                      xf_ens(i,nall+8) =xf_ens(i,nall2+8)
                      xf_ens(i,nall+9) =xf_ens(i,nall2+9)
                      xf_ens(i,nall+15) =xf_ens(i,nall2+15)
                      xf_ens(i,nall+10)=xf_ens(i,nall2+10)
                      xf_ens(i,nall+11)=xf_ens(i,nall2+11)
                      xf_ens(i,nall+12)=xf_ens(i,nall2+12)
!                     if(i.eq.ipr.and.j.eq.jpr)write(0,*)'should not be here'
                   endif
!         if(i.eq.ipr.and.j.eq.jpr)write(0,*)'xff_ens3=',xff_ens3
                go to 100
             endif
          elseif(ierr(i).ne.20.and.ierr(i).ne.0)then
             do n=1,ensdim
               xf_ens(i,n)=0.
	       !srf begin
	       xf_dicycle(i) = 0.
	       !srf end	      
             enddo
          endif
!         if(i.eq.ipr.and.j.eq.jpr)write(0,*)'xff_ens3=',xff_ens3
 100   continue

   END SUBROUTINE cup_forcing_ens_3d

   SUBROUTINE cup_kbcon(ierrc,cap_inc,iloop_in,k22,kbcon,he_cup,hes_cup, &
              hkb,ierr,kbmax,p_cup,cap_max,                         &
              ztexec,zqexec,use_excess,       &
              jprnt,itf,jtf,ktf,                        &
              its,ite, jts,jte, kts,kte, &
              z_cup,entr_rate,heo,imid                        )

   IMPLICIT NONE
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        use_excess,jprnt,itf,jtf,ktf,imid,           &
        its,ite, jts,jte, kts,kte
  ! 
  ! 
  ! 
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        he_cup,hes_cup,p_cup
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        entr_rate,ztexec,zqexec,cap_inc,cap_max
     real,    dimension (its:ite)                                      &
        ,intent (inout   )                   ::                           &
        hkb !,cap_max
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbmax
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        kbcon,k22,ierr
     integer                                                           &
        ,intent (in   )                   ::                           &
        iloop_in
     character*50 :: ierrc(its:ite)
     real, dimension (its:ite,kts:kte),intent (in) :: z_cup,heo
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k,k1,k2,iloop
     real                                 ::                           &
        pbcdif,plus,hetest,dz
     real, dimension (its:ite,kts:kte) ::hcot
!
!--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
!
      iloop=iloop_in
       DO 27 i=its,itf
      kbcon(i)=1
!
! reset iloop for mid level convection
      if(cap_max(i).gt.200 .and. imid.eq.1)iloop=5
!
      IF(ierr(I).ne.0)GO TO 27
      KBCON(I)=K22(I)+1
      if(iloop.eq.5)KBCON(I)=K22(I)
       !== including entrainment for hetest
        hcot(i,1:k22(i)) = HKB(I)
        do k=k22(i)+1,KBMAX(i)+3
           dz=z_cup(i,k)-z_cup(i,k-1)
           hcot(i,k)= ( (1.-0.5*entr_rate(i)*dz)*hcot(i,k-1)   &
                         + entr_rate(i)*dz*heo(i,k-1)       )/ &
                      (1.+0.5*entr_rate(i)*dz)
        enddo
       !==

      GO TO 32
 31   CONTINUE
      KBCON(I)=KBCON(I)+1
      IF(KBCON(I).GT.KBMAX(i)+2)THEN
         if(iloop.ne.4)then
                ierr(i)=3
                ierrc(i)="could not find reasonable kbcon in cup_kbcon"
         endif
        GO TO 27
      ENDIF
 32   CONTINUE
      hetest=hcot(i,kbcon(i)) !hkb(i) ! HE_cup(I,K22(I))
!      if(iloop.eq.5)then
!       hetest=HKB(I)
!      endif
      IF(HETEST.LT.HES_cup(I,KBCON(I)))then
!       if(jprnt.eq.1)write(0,*)'htest',k22(i),kbcon(i),HETEST,-P_cup(I,KBCON(I))+P_cup(I,K22(I))
        GO TO 31
      endif

!     cloud base pressure and max moist static energy pressure
!     i.e., the depth (in mb) of the layer of negative buoyancy
      if(KBCON(I)-K22(I).eq.1)go to 27
      if(iloop.eq.5 .and. (KBCON(I)-K22(I)).le.2)go to 27
      PBCDIF=-P_cup(I,KBCON(I))+P_cup(I,K22(I))
      plus=max(25.,cap_max(i)-float(iloop-1)*cap_inc(i))
      if(iloop.eq.4)plus=cap_max(i)
!
! for shallow convection, if cap_max is greater than 25, it is the pressure at pbltop
      if(iloop.eq.5)plus=150.
      if(iloop.eq.5.and.cap_max(i).gt.200)pbcdif=-P_cup(I,KBCON(I))+cap_max(i)
      IF(PBCDIF.GT.plus)THEN
!      if(jprnt.eq.1)write(0,*)'htest',k22(i),kbcon(i),plus,-P_cup(I,KBCON(I))+P_cup(I,K22(I))
        K22(I)=K22(I)+1
        KBCON(I)=K22(I)+1
        hkb(i)=he_cup(i,k22(i))
       !== including entrainment for hetest
        hcot(i,1:k22(i)) = HKB(I)
        do k=k22(i)+1,KBMAX(i)+3
           dz=z_cup(i,k)-z_cup(i,k-1)

           hcot(i,k)= ( (1.-0.5*entr_rate(i)*dz)*hcot(i,k-1)   &
                         + entr_rate(i)*dz*heo(i,k-1)       )/ &
                      (1.+0.5*entr_rate(i)*dz)
        enddo
       !==

        if(iloop.eq.5)KBCON(I)=K22(I)
        IF(KBCON(I).GT.KBMAX(i)+2)THEN
            if(iloop.ne.4)then
                ierr(i)=3
                ierrc(i)="could not find reasonable kbcon in cup_kbcon"
            endif
            GO TO 27
        ENDIF
        GO TO 32
      ENDIF
 27   CONTINUE

   END SUBROUTINE cup_kbcon


   SUBROUTINE cup_ktop(ierrc,ilo,dby,kbcon,ktop,ierr,              &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,jtf,ktf,           &
        its,ite, jts,jte, kts,kte
  ! dby = buoancy term
  ! ktop = cloud top (output)
  ! ilo  = flag
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout)                   ::                           &
        dby
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon
     integer                                                           &
        ,intent (in   )                   ::                           &
        ilo
     integer, dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        ktop
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     character*50 :: ierrc(its:ite)
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k
!
        DO 42 i=its,itf
        ktop(i)=1
         IF(ierr(I).EQ.0)then
          DO 40 K=KBCON(I)+1,ktf-1
            IF(DBY(I,K).LE.0.)THEN
                KTOP(I)=K-1
                GO TO 41
             ENDIF
  40      CONTINUE
          if(ilo.eq.1)ierr(i)=5
          if(ilo.eq.1)ierrc(i)="problem with defining ktop"
!         if(ilo.eq.2)ierr(i)=998
          GO TO 42
  41     CONTINUE
         do k=ktop(i)+1,ktf
           dby(i,k)=0.
         enddo
         if(kbcon(i).eq.ktop(i))then
            ierr(i)=55
            ierrc(i)="kbcon == ktop "
         endif
         endif
  42     CONTINUE

   END SUBROUTINE cup_ktop


   SUBROUTINE cup_MAXIMI(ARRAY,KS,KE,MAXX,ierr,              &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
         itf,jtf,ktf,                                    &
         its,ite, jts,jte, kts,kte
  ! array input array
  ! x output array with return values
  ! kt output array of levels
  ! ks,kend  check-range
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
         array
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
         ierr,ke
     integer                                                           &
        ,intent (in   )                   ::                           &
         ks
     integer, dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
         maxx
     real,    dimension (its:ite)         ::                           &
         x
     real                                 ::                           &
         xar
     integer                              ::                           &
         i,k

       DO 200 i=its,itf
       MAXX(I)=KS
       if(ierr(i).eq.0)then
      X(I)=ARRAY(I,KS)
!
       DO 100 K=KS,KE(i)
         XAR=ARRAY(I,K)
         IF(XAR.GE.X(I)) THEN
            X(I)=XAR
            MAXX(I)=K
         ENDIF
 100  CONTINUE
      endif
 200  CONTINUE

   END SUBROUTINE cup_MAXIMI


   SUBROUTINE cup_minimi(ARRAY,KS,KEND,KT,ierr,              &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
         itf,jtf,ktf,                                    &
         its,ite, jts,jte, kts,kte
  ! array input array
  ! x output array with return values
  ! kt output array of levels
  ! ks,kend  check-range
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
         array
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
         ierr,ks,kend
     integer, dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
         kt
     real,    dimension (its:ite)         ::                           &
         x
     integer                              ::                           &
         i,k,kstop

       DO 200 i=its,itf
      KT(I)=KS(I)
      if(ierr(i).eq.0)then
      X(I)=ARRAY(I,KS(I))
       KSTOP=MAX(KS(I)+1,KEND(I))
!
       DO 100 K=KS(I)+1,KSTOP
         IF(ARRAY(I,K).LT.X(I)) THEN
              X(I)=ARRAY(I,K)
              KT(I)=K
         ENDIF
 100  CONTINUE
      endif
 200  CONTINUE

   END SUBROUTINE cup_MINIMI


   SUBROUTINE cup_up_aa0(aa0,z,zu,dby,GAMMA_CUP,t_cup,       &
              kbcon,ktop,ierr,                               &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,jtf,ktf,                                     &
        its,ite, jts,jte, kts,kte
  ! aa0 cloud work function
  ! gamma_cup = gamma on model cloud levels
  ! t_cup = temperature (Kelvin) on model cloud levels
  ! dby = buoancy term
  ! zu= normalized updraft mass flux
  ! z = heights of model levels 
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        z,zu,gamma_cup,t_cup,dby
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon,ktop
!
! input and output
!


     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        aa0
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k
     real                                 ::                           &
        dz,da
!
        do i=its,itf
         aa0(i)=0.
        enddo
        DO 100 k=kts+1,ktf
        DO 100 i=its,itf
         IF(ierr(i).ne.0)GO TO 100
         IF(K.LT.KBCON(I))GO TO 100
         IF(K.Gt.KTOP(I))GO TO 100
         DZ=Z(I,K)-Z(I,K-1)
         da=zu(i,k)*DZ*(9.81/(1004.*( &
                (T_cup(I,K)))))*DBY(I,K-1)/ &
             (1.+GAMMA_CUP(I,K))
!         IF(K.eq.KTOP(I).and.da.le.0.)go to 100
         AA0(I)=AA0(I)+max(0.,da)
         if(aa0(i).lt.0.)aa0(i)=0.
100     continue

   END SUBROUTINE cup_up_aa0

!====================================================================
   SUBROUTINE neg_check(name,j,dt,q,outq,outt,outu,outv,   &
                                    outqc,pret,its,ite,kts,kte,itf,ktf)

   INTEGER,      INTENT(IN   ) ::            j,its,ite,kts,kte,itf,ktf

     real, dimension (its:ite,kts:kte  )                    ,                 &
      intent(inout   ) ::                                                     &
       outq,outt,outqc,outu,outv
     real, dimension (its:ite,kts:kte  )                    ,                 &
      intent(inout   ) ::                                                     &
       q
     real, dimension (its:ite  )                            ,                 &
      intent(inout   ) ::                                                     &
       pret
      character *(*), intent (in)         ::                           &
       name
     real                                                                     &
        ,intent (in  )                   ::                                   &
        dt
     real :: names,scalef,thresh,qmem,qmemf,qmem2,qtest,qmem1
     integer :: icheck
!
! first do check on vertical heating rate
!
      icheck=0
      thresh=300.01
      thresh=200.01	!ss
      thresh=250.01
      names=1.
      if(name == 'shallow')then
        thresh=48.01
        names=2.
      endif
      scalef=86400.
      do i=its,itf
      qmemf=1.
      qmem=0.
      do k=kts,ktf
         qmem=(outt(i,k))*86400.
         if(qmem.gt.thresh)then
           qmem2=thresh/qmem
           qmemf=min(qmemf,qmem2)
      icheck=1
!
!
!          print *,'1',' adjusted massflux by factor ',i,j,k,qmem,qmem2,qmemf,dt
         endif
         if(qmem.lt.-.5*thresh*names)then
           qmem2=-.5*names*thresh/qmem
           qmemf=min(qmemf,qmem2)
      icheck=2
!
!
!!sms$ignore begin
!          write(0,*)name, icheck,' adjusted massflux by factor ',qmem2, 'at j,k = ',j,k
!!sms$ignore end
         endif
      enddo
!    if(name == 'shallow')then
 !    if(j.eq.586981)then
!!sms$ignore begin
!      if(qmemf.lt.0.5)then
!           write(12,*)'1',' adjusted massflux by factor ',j,qmemf,icheck,name
!        do k=kts,ktf-5
!         write(12,*)k,scalef*(outq(i,k)),scalef*qmemf*(outt(i,k)+)
!        enddo
!      endif
!!sms$ignore end
!!     endif
!      endif
      do k=kts,ktf
         outq(i,k)=outq(i,k)*qmemf
         outt(i,k)=outt(i,k)*qmemf
         outu(i,k)=outu(i,k)*qmemf
         outv(i,k)=outv(i,k)*qmemf
         outqc(i,k)=outqc(i,k)*qmemf
      enddo
      pret(i)=pret(i)*qmemf 
      enddo
!
! check whether routine produces negative q's. This can happen, since 
! tendencies are calculated based on forced q's. This should have no
! influence on conservation properties, it scales linear through all
! tendencies
!
      thresh=1.e-10
      do i=its,itf
      qmemf=1.
      do k=kts,ktf-1
         qmem=outq(i,k)
         if(abs(qmem).gt.0.)then
         qtest=q(i,k)+(outq(i,k))*dt
         if(qtest.lt.thresh)then
!
! qmem2 would be the maximum allowable tendency
!
           qmem1=outq(i,k)
           qmem2=(thresh-q(i,k))/dt
           qmemf=min(qmemf,qmem2/qmem1)
           qmemf=max(0.,qmemf)
!          write(0,*)'4 adjusted tendencies ',i,k,qmem,qmem2,qmemf
!          write(0,*)'4 adjusted tendencies ',i,j,k,q(i,k),qmem1,qmemf
         endif
         endif
      enddo
!     if(qmemf.lt.1.)write(0,*)'4 adjusted tendencies ',i,j,qmemf
      do k=kts,ktf
         outq(i,k)=outq(i,k)*qmemf
         outt(i,k)=outt(i,k)*qmemf
         outu(i,k)=outu(i,k)*qmemf
         outv(i,k)=outv(i,k)*qmemf
         outqc(i,k)=outqc(i,k)*qmemf
      enddo
      pret(i)=pret(i)*qmemf 
      enddo

   END SUBROUTINE neg_check


   SUBROUTINE cup_output_ens_3d(xff_mid,xf_ens,ierr,dellat,dellaq,dellaqc,  &
              outtem,outq,outqc,     &
              zu,sub_mas,pre,pw,xmb,ktop,                 &
              j,name,nx,nx2,iens,ierr2,ierr3,pr_ens,             &
              maxens3,ensdim,                            &
              sig,closure_n,xland1,xmbm_in,xmbs_in,    &
              weight_GR,weight_W,weight_MC,weight_ST,weight_AS,training,  &
	      ichoice,imid,ipr,jpr,itf,jtf,ktf, &
              its,ite, jts,jte, kts,kte, &
	      dicycle,xf_dicycle )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        ichoice,imid,ipr,jpr,itf,jtf,ktf,     &
        its,ite, jts,jte, kts,kte
     integer, intent (in   )              ::                           &
        j,ensdim,nx,nx2,iens,maxens3,training
  ! xf_ens = ensemble mass fluxes
  ! pr_ens = precipitation ensembles
  ! dellat = change of temperature per unit mass flux of cloud ensemble
  ! dellaq = change of q per unit mass flux of cloud ensemble
  ! dellaqc = change of qc per unit mass flux of cloud ensemble
  ! outtem = output temp tendency (per s)
  ! outq   = output q tendency (per s)
  ! outqc  = output qc tendency (per s)
  ! pre    = output precip
  ! xmb    = total base mass flux
  ! xfac1  = correction factor
  ! pw = pw -epsilon*pd (ensemble dependent)
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,1:ensdim)                     &
        ,intent (inout)                   ::                           &
       xf_ens,pr_ens
!srf ------
     real, dimension( its:ite )                      &
         ,intent(in) :: weight_gr,weight_w,weight_mc,weight_st,weight_as
!-srf---
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout  )                   ::                           &
        outtem,outq,outqc,sub_mas
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in  )                   ::                           &
        zu
     real,   dimension (its:ite)                                      &
         ,intent (in  )                   ::                           &
        sig,xmbm_in,xmbs_in
     real,   dimension (its:ite,2)                                      &
         ,intent (in  )                   ::                           &
        xff_mid
     real,    dimension (its:ite)                                      &
        ,intent (inout  )                   ::                           &
        pre,xmb
     real,    dimension (its:ite)                                      &
        ,intent (inout  )                   ::                           &
        closure_n
     real,    dimension (its:ite,kts:kte,1:nx)                     &
        ,intent (in   )                   ::                           &
       dellat,dellaqc,dellaq,pw
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        ktop,xland1
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr,ierr2,ierr3
!-srf begin
      integer, intent(IN) :: DICYCLE
      real,    intent(IN), dimension (its:ite) :: xf_dicycle
!-srf end
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k,n,ncount
     real                                 ::                           &
        outtes,ddtes,dtt,dtq,dtqc,dtpw,prerate,clos_wei,xmbhelp
     real                                 ::                           &
        dtts,dtqs
     real,    dimension (its:ite)         ::                           &
       xfac1,xfac2
     real,    dimension (its:ite)::                           &
       xmb_ske,xmb_ave,xmb_std,xmb_cur,xmbweight
     real,    dimension (its:ite)::                           &
       pr_ske,pr_ave,pr_std,pr_cur
     real,    dimension (its:ite)::                           &
               pr_gr,pr_w,pr_mc,pr_st,pr_as,pr_capma,     &
               pr_capme,pr_capmi
     real, dimension (5) :: weight,wm,wm1,wm2,wm3
     real, dimension (its:ite,5) :: xmb_w

!
      character *(*), intent (in)        ::                           &
       name

!
     weight(1) = -999.  !this will turn off weights
     wm(1)=-999.

!
!
      DO k=kts,kte
      do i=its,ite
        outtem(i,k)=0.
        outq(i,k)=0.
        outqc(i,k)=0.
        sub_mas(i,k)=0.
      enddo
      enddo
      do i=its,itf
        pre(i)=0.
        xmb(i)=0.
         xfac1(i)=0.
         xfac2(i)=0.
        xmbweight(i)=1.
      enddo
      do i=its,itf
        IF(ierr(i).eq.0)then
        do n=(iens-1)*nx*nx2*maxens3+1,iens*nx*nx2*maxens3
           if(pr_ens(i,n).le.0.)then
!            if(i.eq.ipr.and.j.eq.jpr)write(0,*)'pr_ens',n,pr_ens(i,n),xf_ens(i,n)
             xf_ens(i,n)=0.
           endif
        enddo
        endif
      enddo
!
!--- calculate ensemble average mass fluxes
!
       
       xmb_w=0.
!
!-- now do feedback
!
      ddtes=100.
      if(imid.eq.0)then
      do i=its,itf
        if(ierr(i).eq.0)then
         k=0
         xmb_ave(i)=0.
         do n=(iens-1)*nx*nx2*maxens3+1,iens*nx*nx2*maxens3
          k=k+1
          xmb_ave(i)=xmb_ave(i)+xf_ens(i,n)
!!sms$ignore begin
!          if(i.eq.ipr.and.j.eq.jpr)write(13,*)'xmb_ave = ',k,xmb_ave(i),xf_ens(i,n)
!!sms$ignore end
         enddo
         xmb_ave(i)=xmb_ave(i)/float(k)
	 !srf begin
	 !print*,"XMB=",xmb_ave(i),xf_dicycle(i),xmb_ave(i) - xf_dicycle(i);call  flush(6)
         if(dicycle == 2 )then
            xmb_ave(i)=max(0.,xmb_ave(i)-xmbm_in(i)-xmbs_in(i))
         else if (dicycle == 1) then
            xmb_ave(i)=min(xmb_ave(i),xmb_ave(i) - xf_dicycle(i))
            xmb_ave(i)=max(0.,xmb_ave(i))
         endif
         xmb_ave(i)=min(xmb_ave(i),100.)
         xmb(i)=sig(i)*xmb_ave(i)
! --- Now use proper count of how many closures were actually
!       used in cup_forcing_ens (including screening of some
!       closures over water) to properly normalize xmb
           clos_wei=16./max(1.,closure_n(i))

           if(xmb(i).eq.0.)then
              ierr(i)=19
           endif
           xfac1(i)=xmb(i)
           xfac2(i)=xmb(i)

        endif
      ENDDO
      else  ! imid == 1
         do i=its,itf
         xmb_ave(i)=0.
         IF(ierr(i).eq.0)then
! ! first get xmb_ves, depend on icoice
!
           xmb_ave(i)=.5*sig(i)*(xff_mid(i,1)+xff_mid(i,2))
           if(ichoice.eq.1 .or. ichoice.eq.2)xmb_ave(i)=sig(i)*xff_mid(i,ichoice)
           if(ichoice.gt.2)then
              k=0
              do n=(iens-1)*nx*nx2*maxens3+1,iens*nx*nx2*maxens3
                    k=k+1
                    xmb_ave(i)=xmb_ave(i)+xf_ens(i,n)
              enddo
              xmb_ave(i)=xmb_ave(i)/float(k)
           endif   ! ichoice gt.2
! which dicycle method
           if(dicycle == 2 )then
              xmb(i)=max(0.,xmb_ave(i)-xmbs_in(i))
           else if (dicycle == 1) then
              xmb(i)=min(xmb_ave(i),xmb_ave(i) - xf_dicycle(i))
              xmb(i)=max(0.,xmb_ave(i))
           else if (dicycle == 0) then
              xmb(i)=max(0.,xmb_ave(i))
           endif   ! dicycle=1,2
         endif     ! ierr >0
         enddo     ! i
      endif        ! imid=1

      do i=its,itf
        IF(ierr(i).eq.0)then
            DO k=kts,ktop(i)
            dtt =0.
            dtts=0.
            dtq =0.
            dtqs=0.
            dtqc=0.
            dtpw=0.
            do n=1,nx
              dtt =dtt  + dellat  (i,k,n)
              dtq =dtq  + dellaq  (i,k,n)
              dtqc=dtqc + dellaqc (i,k,n)
              dtpw=dtpw + pw      (i,k,n)
           enddo
           OUTTEM(I,K)= XMB(I)* dtt /float(nx)
           OUTQ  (I,K)= XMB(I)* dtq /float(nx)
           OUTQC (I,K)= XMB(I)* dtqc/float(nx)
	   PRE(I)=PRE(I)+XMB(I)*dtpw/float(nx)
           xf_ens(i,:)=sig(i)*xf_ens(i,:) !*dtpw/float(nx)
           sub_mas(i,k)=zu(i,k)*xmb(i)
          enddo
        endif
      enddo

!     do i=its,itf
!       if(ierr(i).eq.0)then
!       do k=(iens-1)*nx*nx2*maxens3+1,iens*nx*nx2*maxens3
!         xf_ens(i,k)=xf_ens(i,k)*xfac1(i)
!       enddo
!       endif
!     ENDDO

   END SUBROUTINE cup_output_ens_3d
!-------------------------------------------------------
   SUBROUTINE cup_up_moisture(name,ierr,z_cup,qc,qrc,pw,pwav,     &
              ccnclean,p_cup,kbcon,ktop,cd,dby,clw_all,&
              t_cup,q,GAMMA_cup,zu,qes_cup,k22,qe_cup,         &
              ZQEXEC,use_excess,ccn,rho,c1d, &
              up_massentr,up_massdetr,psum,psumh,                 &
              autoconv,aeroevap,itest,itf,jtf,ktf,j,ipr,jpr,                &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
  real, parameter :: BDISPM = 0.366       !Berry--size dispersion (martime)
  REAL, PARAMETER :: BDISPC = 0.146       !Berry--size dispersion (continental)
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  use_excess,itest,autoconv,aeroevap,itf,jtf,ktf,           &
                                  its,ite, jts,jte,j,ipr,jpr, kts,kte
  ! cd= detrainment function 
  ! q = environmental q on model levels
  ! qe_cup = environmental q on model cloud levels
  ! qes_cup = saturation q on model cloud levels
  ! dby = buoancy term
  ! cd= detrainment function 
  ! zu = normalized updraft mass flux
  ! gamma_cup = gamma on model cloud levels
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        t_cup,p_cup,rho,q,zu,gamma_cup,qe_cup,                         &
        up_massentr,up_massdetr,dby,qes_cup,z_cup,cd
     real,    dimension (its:ite)                              &
        ,intent (in   )                   ::                           &
        zqexec
  ! entr= entrainment rate 
     real                                                              &
        ,intent (in   )                   ::                           &
        ccnclean
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon,ktop,k22
!
! input and output
!

   ! ierr error value, maybe modified in this routine

     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
      character *(*), intent (in)        ::                           &
       name
   ! qc = cloud q (including liquid water) after entrainment
   ! qrch = saturation q in cloud
   ! qrc = liquid water content in cloud after rainout
   ! pw = condensate that will fall out at that level
   ! pwav = totan normalized integrated condensate (I1)
   ! c0 = conversion rate (cloud to rain)

     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        qc,qrc,pw,clw_all
     real,    dimension (its:ite,kts:kte) ::                           &
        qch,qrcb,pwh,clw_allh,c1d
     real,    dimension (its:ite)         ::                           &
        pwavh
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        pwav,psum,psumh
     real,    dimension (its:ite)                                      &
        ,intent (in  )                   ::                           &
        ccn
!
!  local variables in this routine
!

     integer                              ::                           &
        iounit,iprop,iall,i,k,k1,k2
     real                                 ::                           &
        prop_ave,qrcb_h,bdsp,dp,g,rhoc,dh,qrch,c0,dz,radius,berryc0,q1,berryc
     real                                 ::                           &
        denom
     real,    dimension (kts:kte)         ::                           &
        prop_b
!
        prop_b(kts:kte)=0
        iall=0
        c0=.002
        g=9.81
        bdsp=BDISPM
!
!--- no precip for small clouds
!
        if(name.eq.'shallow')then
            c0=0.002
!           iall=1
        endif
        do i=its,itf
          pwav(i)=0.
          pwavh(i)=0.
          psum(i)=0.
          psumh(i)=0.
        enddo
        do k=kts,ktf
        do i=its,itf
          pw(i,k)=0.
          pwh(i,k)=0.
          qc(i,k)=0.
          if(ierr(i).eq.0)qc(i,k)=qe_cup(i,k)
          if(ierr(i).eq.0)qch(i,k)=qe_cup(i,k)
          clw_all(i,k)=0.
          clw_allh(i,k)=0.
          qrc(i,k)=0.
          qrcb(i,k)=0.
        enddo
        enddo
!     if(use_excess < 2 ) then
      do i=its,itf
      if(ierr(i).eq.0.)then
!
!  initialize below originating air
!
      do k=2,k22(i)
        DZ=Z_cup(i,K)-Z_cup(i,K-1)
        qc(i,k)=qe_cup(i,k)!+float(use_excess)*zqexec(i)
        qch(i,k)=qe_cup(i,k)!+float(use_excess)*zqexec(i)
      enddo
      endif
      enddo

       DO 100 i=its,itf
         IF(ierr(i).eq.0)then

! below LFC, but maybe above LCL
!
            DO k=k22(i)+1,kbcon(i)-1
              qc(i,k)=   (qc(i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k-1)* qc(i,k-1)+ &
                         up_massentr(i,k-1)*q(i,k-1))   /            &
                         (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
!              QRCH=QES_cup(I,K)
               QRCH=QES_cup(I,K)+(1./XLV)*(GAMMA_cup(i,k) &
                 /(1.+GAMMA_cup(i,k)))*DBY(I,K)
              if(qc(i,k).gt.qrch)then
                DZ=Z_cup(i,K)-Z_cup(i,K-1)
                QRC(I,K)=(QC(I,K)-QRCH)/(1.+(c0+c1d(i,k))*DZ)
                PW(i,k)=c0*dz*QRC(I,K)*zu(i,k)
                qc(i,k)=qrch+qrc(i,k)
              endif
            enddo
!
!now do the rest
!
            DO k=kbcon(i),ktop(i)
               denom=zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1)
               if(denom.lt.1.e-8)then
!          write(14,*)j,k,denom,zuo(i,k-1),up_massdetro(i,k-1),up_massentro(i,k-1)
                     ierr(i)=51
                exit
               endif

   
               rhoc=.5*(rho(i,k)+rho(i,k-1))
               DZ=Z_cup(i,K)-Z_cup(i,K-1)
               DP=p_cup(i,K)-p_cup(i,K-1)
!
!--- saturation  in cloud, this is what is allowed to be in it
!
               QRCH=QES_cup(I,K)+(1./XLV)*(GAMMA_cup(i,k) &
                 /(1.+GAMMA_cup(i,k)))*DBY(I,K)
!
!------    1. steady state plume equation, for what could
!------       be in cloud without condensation
!
!
               qc(i,k)=   (qc(i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k-1)* qc(i,k-1)+ &
                         up_massentr(i,k-1)*q(i,k-1))   /            &
                         (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
               qch(i,k)= (qch(i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k-1)*qch(i,k-1)+ &
                         up_massentr(i,k-1)*q(i,k-1))   /            &
                         (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))

               if(qc(i,k).le.qrch)then
                 qc(i,k)=qrch
               endif
               if(qch(i,k).le.qrch)then
                 qch(i,k)=qrch
               endif
!
!------- Total condensed water before rainout
!
               clw_all(i,k)=max(0.,QC(I,K)-QRCH)
               QRC(I,K)=max(0.,(QC(I,K)-QRCH)) ! /(1.+C0*DZ*zu(i,k))
               clw_allh(i,k)=max(0.,QCH(I,K)-QRCH)
               QRCB(I,K)=max(0.,(QCH(I,K)-QRCH)) ! /(1.+C0*DZ*zu(i,k))
               IF(autoconv.eq.2) then


! 
! normalized berry
!
! first calculate for average conditions, used in cup_dd_edt!
! this will also determine proportionality constant prop_b, which, if applied,
! would give the same results as c0 under these conditions
!
                 q1=1.e3*rhoc*qrcb(i,k)  ! g/m^3 ! g[h2o]/cm^3
                 berryc0=q1*q1/(60.0*(5.0 + 0.0366*CCNclean/ &
                    ( q1 * BDSP)  ) ) !/(
!     if(i.eq.ipr.and.j.eq.jpr)write(0,*)'cupm',k,rhoc,rho(i,k)
!        qrcb_h=qrcb(i,k)/(1.+c0*dz)
                 qrcb_h=((QCH(I,K)-QRCH)*zu(i,k)-qrcb(i,k-1)*(.5*up_massdetr(i,k-1)))/ &
                   (zu(i,k)+.5*up_massdetr(i,k-1)+c0*dz*zu(i,k))
                 prop_b(k)=c0*qrcb_h*zu(i,k)/(1.e-3*berryc0)
!     if(i.eq.ipr.and.j.eq.jpr)write(0,*)'cupm',berryc0,prop_b(k),qrcb_h
                 pwh(i,k)=zu(i,k)*1.e-3*berryc0*dz*prop_b(k) ! 2.
                 berryc=qrcb(i,k)
                 qrcb(i,k)=((QCh(I,K)-QRCH)*zu(i,k)-pwh(i,k)-qrcb(i,k-1)*(.5*up_massdetr(i,k-1)))/ &
                       (zu(i,k)+.5*up_massdetr(i,k-1))
                 if(qrcb(i,k).lt.0.)then
                   berryc0=(qrcb(i,k-1)*(.5*up_massdetr(i,k-1))-(QCh(I,K)-QRCH)*zu(i,k))/zu(i,k)*1.e-3*dz*prop_b(k)
                   pwh(i,k)=zu(i,k)*1.e-3*berryc0*dz*prop_b(k)
                   qrcb(i,k)=0.
                 endif
!     if(i.eq.ipr.and.j.eq.jpr)write(0,*)'cupm',zu(i,k),pwh(i,k),dz,qrch,qrcb(i,k),clw_allh(i,k)
                 QCh(I,K)=QRCb(I,K)+qrch
                 PWAVH(I)=PWAVH(I)+pwh(I,K)
                 Psumh(I)=Psumh(I)+clw_allh(I,K)*zu(i,k) *dz
        !
! then the real berry
!
                 q1=1.e3*rhoc*qrc(i,k)  ! g/m^3 ! g[h2o]/cm^3
                 berryc0=q1*q1/(60.0*(5.0 + 0.0366*CCN(i)/ &
                    ( q1 * BDSP)  ) ) !/(
                 berryc0=1.e-3*berryc0*dz*prop_b(k) ! 2.
                 berryc=qrc(i,k)
                 qrc(i,k)=((QC(I,K)-QRCH)*zu(i,k)-zu(i,k)*berryc0-qrc(i,k-1)*(.5*up_massdetr(i,k-1)))/ &
                       (zu(i,k)+.5*up_massdetr(i,k-1))
                 if(qrc(i,k).lt.0.)then
                    berryc0=((QC(I,K)-QRCH)*zu(i,k)-qrc(i,k-1)*(.5*up_massdetr(i,k-1)))/zu(i,k)
                    qrc(i,k)=0.
                 endif
                 pw(i,k)=berryc0*zu(i,k)
                 QC(I,K)=QRC(I,K)+qrch
!
!  if not running with berry at all, do the following
!
               ELSE       !c0=.002
                 if(iall.eq.1)then
                   qrc(i,k)=0.
                   pw(i,k)=(QC(I,K)-QRCH)*zu(i,k)
                   if(pw(i,k).lt.0.)pw(i,k)=0.
                 else
                   QRC(I,K)=(QC(I,K)-QRCH)/(1.+(c1d(i,k)+C0)*DZ)
                   PW(i,k)=c0*dz*QRC(I,K)*zu(i,k)
          !        qrc(i,k)=((QC(I,K)-QRCH)*zu(i,k)-qrc(i,k-1)*(.5*up_massdetr(i,k-1)))/ &
!                  (zu(i,k)+.5*up_massdetr(i,k-1)+c0*dz*zu(i,k))
!        PW(i,k)=c0*dz*qrc(i,k)*zu(i,k)
                   if(qrc(i,k).lt.0)then
                     qrc(i,k)=0.
                     pw(i,k)=0.
                   endif
                 endif
                 QC(I,K)=QRC(I,K)+qrch
               endif !autoconv
               PWAV(I)=PWAV(I)+PW(I,K)
               Psum(I)=Psum(I)+clw_all(I,K)*zu(i,k) *dz
            enddo ! k=kbcon,ktop
      endif ! ierr
!
!--- integrated normalized ondensate
!
 100     CONTINUE
       prop_ave=0.
       iprop=0
       do k=kts,kte
        prop_ave=prop_ave+prop_b(k)
        if(prop_b(k).gt.0)iprop=iprop+1
       enddo
       iprop=max(iprop,1)
!      write(11,*)'prop_ave = ',prop_ave/float(iprop)
!      print *,'pwav = ',pwav(1)

   END SUBROUTINE cup_up_moisture
!====================================================================

!--------------------------------------------------------------------

      real function satvap(temp2)
      implicit none
      real :: temp2, temp, toot, toto, eilog, tsot,  &
     &        ewlog, ewlog2, ewlog3, ewlog4
      temp = temp2-273.155
      if (temp.lt.-20.) then   !!!! ice saturation
        toot = 273.16 / temp2
        toto = 1 / toot
        eilog = -9.09718 * (toot - 1) - 3.56654 * (log(toot) / &
     &    log(10.)) + .876793 * (1 - toto) + (log(6.1071) / log(10.))
        satvap = 10 ** eilog
      else
        tsot = 373.16 / temp2
        ewlog = -7.90298 * (tsot - 1) + 5.02808 * &
     &             (log(tsot) / log(10.))
        ewlog2 = ewlog - 1.3816e-07 * &
     &             (10 ** (11.344 * (1 - (1 / tsot))) - 1)
        ewlog3 = ewlog2 + .0081328 * &
     &             (10 ** (-3.49149 * (tsot - 1)) - 1)
        ewlog4 = ewlog3 + (log(1013.246) / log(10.))
        satvap = 10 ** ewlog4
      end if
      return
      end function
   SUBROUTINE CUP_gf_sh(cnvwt,zuo,zws,xmb_out,xmbs2,zo,OUTQC,J,AAEQ,T,Q,Z1,       &
              pre,TN,QO,PO,P,OUTT,OUTQ,DTIME,PSUR,US,VS,                  &
              TCRIT,ztexec,zqexec,ccn,ccnclean,rho,dx,dhdt,           &
              cupclw,kpbl,kbcon,ktop,k22,xland,         &   !-lxz
              tscl_kf,ichoice,ipr,jpr,ierr,ierrc,         &
              autoconv,itf,jtf,ktf,               &
              use_excess,its,ite, jts,jte, kts,kte &
                                                )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        autoconv,itf,jtf,ktf,use_excess,        &
        its,ite, jts,jte, kts,kte,ipr,jpr
     integer, intent (in   )              ::                           &
        j,ichoice
  !
  ! 
  !
  ! outtem = output temp tendency (per s)
  ! outq   = output q tendency (per s)
  ! outqc  = output qc tendency (per s)
  ! pre    = output precip
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout  )                   ::                           &
        cnvwt,OUTT,OUTQ,OUTQC,cupclw
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        xmb_out,xmbs2
     integer,    dimension (its:ite)                                   &
        ,intent (out  )                   ::                           &
        kbcon,ktop,k22
     integer,    dimension (its:ite)                                   &
        ,intent (in  )                   ::                           &
        kpbl
  !
  ! basic environmental input includes moisture convergence (mconv)
  ! omega (omeg), windspeed (us,vs), and a flag (aaeq) to turn off
  ! convection for this call only and at that particular gridpoint
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        rho,T,PO,P,US,VS,tn,dhdt
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout)                   ::                           &
         Q,QO
     real, dimension (its:ite)                                         &
        ,intent (in   )                   ::                           &
        xland,zws,ztexec,zqexec,ccn,Z1,PSUR,AAEQ
       
       real                                                            &
        ,intent (in   )                   ::                           &
        tscl_kf,dx,ccnclean,dtime,tcrit


!
!
!***************** the following are your basic environmental
!                  variables. They carry a "_cup" if they are
!                  on model cloud levels (staggered). They carry
!                  an "o"-ending (z becomes zo), if they are the forced
!                  variables. They are preceded by x (z becomes xz)
!                  to indicate modification by some typ of cloud
!
  ! z           = heights of model levels
  ! q           = environmental mixing ratio
  ! qes         = environmental saturation mixing ratio
  ! t           = environmental temp
  ! p           = environmental pressure
  ! he          = environmental moist static energy
  ! hes         = environmental saturation moist static energy
  ! z_cup       = heights of model cloud levels
  ! q_cup       = environmental q on model cloud levels
  ! qes_cup     = saturation q on model cloud levels
  ! t_cup       = temperature (Kelvin) on model cloud levels
  ! p_cup       = environmental pressure
  ! he_cup = moist static energy on model cloud levels
  ! hes_cup = saturation moist static energy on model cloud levels
  ! gamma_cup = gamma on model cloud levels
!
!
  ! hcd = moist static energy in downdraft
  ! zd normalized downdraft mass flux
  ! dby = buoancy term
  ! entr = entrainment rate
  ! zd   = downdraft normalized mass flux
  ! entr= entrainment rate
  ! hcd = h in model cloud
  ! bu = buoancy term
  ! zd = normalized downdraft mass flux
  ! gamma_cup = gamma on model cloud levels
  ! qcd = cloud q (including liquid water) after entrainment
  ! qrch = saturation q in cloud
  ! pwd = evaporate at that level
  ! pwev = total normalized integrated evaoprate (I2)
  ! entr= entrainment rate
  ! z1 = terrain elevation
  ! entr = downdraft entrainment rate
  ! jmin = downdraft originating level
  ! kdet = level above ground where downdraft start detraining
  ! psur        = surface pressure
  ! z1          = terrain elevation
  ! pr_ens = precipitation ensemble
  ! xf_ens = mass flux ensembles
  ! massfln = downdraft mass flux ensembles used in next timestep
  ! omeg = omega from large scale model
  ! mconv = moisture convergence from large scale model
  ! zd      = downdraft normalized mass flux
  ! zu      = updraft normalized mass flux
  ! dir     = "storm motion"
  ! mbdt    = arbitrary numerical parameter
  ! dtime   = dt over which forcing is applied
  ! iact_gr_old = flag to tell where convection was active
  ! kbcon       = LFC of parcel from k22
  ! k22         = updraft originating level
  ! icoic       = flag if only want one closure (usually set to zero!)
  ! dby = buoancy term
  ! ktop = cloud top (output)
  ! xmb    = total base mass flux
  ! hc = cloud moist static energy
  ! hkb = moist static energy at originating level

     real,    dimension (its:ite,kts:kte) ::                           &
        entr_rate_2d,mentrd_rate_2d,he,hes,qes,z,                      &
        heo,heso,qeso,zo,                                              &
        xhe,xhes,xqes,xz,xt,xq,                                        &

        qes_cup,q_cup,he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup,      &
        qeso_cup,qo_cup,heo_cup,heso_cup,zo_cup,po_cup,gammao_cup,     &
        tn_cup,                                                        &
        xqes_cup,xq_cup,xhe_cup,xhes_cup,xz_cup,xp_cup,xgamma_cup,     &
        xt_cup,                                                        &

        dby,qc,qrcd,pwd,pw,hcd,qcd,dbyd,hc,qrc,zu,zd,clw_all,   &
        dbyo,qco,qrcdo,pwdo,pwo,hcdo,qcdo,dbydo,hco,qrco,zuo,zdo,      &
        xdby,xqc,xqrcd,xpwd,xpw,xhcd,xqcd,xhc,xqrc,xzu,xzd,            &

  ! cd  = detrainment function for updraft
  ! cdd = detrainment function for downdraft
  ! dellat = change of temperature per unit mass flux of cloud ensemble
  ! dellaq = change of q per unit mass flux of cloud ensemble
  ! dellaqc = change of qc per unit mass flux of cloud ensemble

        cd,cdd,DELLAH,DELLAQ,DELLAT,DELLAQC

  ! aa0 cloud work function for downdraft
  ! edt = epsilon
  ! aa0     = cloud work function without forcing effects
  ! aa1     = cloud work function with forcing effects
  ! xaa0    = cloud work function with cloud effects (ensemble dependent)
  ! edt     = epsilon

     real,    dimension (its:ite) ::                                   &
       pre,edt,edto,edtx,AA1,AA0,XAA0,HKB,                          &
       HKBO,XHKB,QKB,QKBO,                                    &
       xmbmax,XMB,XPWAV,XPWEV,PWAV,PWEV,PWAVO,                                &
       PWEVO,BU,BUD,BUO,cap_max,entr_rate,                                    &
       cap_max_increment,closure_n,psum,psumh,sig,zuhe
     integer,    dimension (its:ite) ::                                &
       kzdown,KDET,KB,JMIN,kstabi,kstabm,K22x,        &   !-lxz
       xland1,KBCONx,KBx,KTOPx,ierr,ierr2,ierr3,KBMAX

     integer                              ::                           &
       nall,iedt,nens,nens3,ki,I,K,KK,iresult
     real                                 ::                           &
      day,dz,dzo,mbdt,radius,  &
      zcutdown,edtmax,edtmin,depth_min,zkbmax,z_detr,zktop,      &
      massfld,dh,cap_maxs,trash,trash2,trashq,frh,fsum
      
      real detdo1,detdo2,entdo,dp,subin,detdo,entup,                &
      detup,subdown,entdoj,entupk,detupk,totmas
      real :: zufinals,power_entr,dzm1,dzp1


     integer :: icount,tun_lim,jprnt,k1,k2,kbegzu,kfinalzu,kstart,jmini
     logical :: keep_going
     real xff_shal(9),blqe,xkshal
     character*50 :: ierrc(its:ite)
     real,    dimension (its:ite,kts:kte) ::                           &
       up_massentr,up_massdetr,dd_massentr,dd_massdetr                 &
      ,up_massentro,up_massdetro,dd_massentro,dd_massdetro
     real,    dimension (kts:kte) :: smth
      real zuh2(40)
     real :: C_up
      power_entr=2.
      icount=0
      day=86400.
      do i=its,itf
        xland1(i)=int(xland(i)+.001) ! 1.
        if(xland(i).gt.1.5 .or. xland(i).lt.0.5)then
            xland1(i)=0
        endif
        pre(i)=0.
        xmb_out(i)=0.
        xmbs2(i)=0.
        cap_max_increment(i)=25.
        ierrc(i)=" "
        entr_rate(i) =1.75e-3 ! 1.2e-3 ! .2/50.
      enddo
!
!--- initial entrainment rate (these may be changed later on in the
!--- program
!
      tun_lim=10
      
!
!--- initial detrainmentrates
!
      do k=kts,ktf
      do i=its,itf
        up_massentro(i,k)=0.
        up_massdetro(i,k)=0.
        z(i,k)=zo(i,k)
        xz(i,k)=zo(i,k)
        qrco(i,k)=0.
        pwo(i,k)=0.
        cd(i,k)=1.*entr_rate(i)
        dellaqc(i,k)=0.
        cupclw(i,k)=0.
      enddo
      enddo
!
!--- max/min allowed value for epsilon (ratio downdraft base mass flux/updraft
!
!--- minimum depth (m), clouds must have
!
      depth_min=50.
!
!--- maximum depth (mb) of capping 
!--- inversion (larger cap = no convection)
!
      cap_maxs=125.
      DO i=its,itf
        kbmax(i)=1
        aa0(i)=0.
        aa1(i)=0.
      enddo
      do i=its,itf
          cap_max(i)=cap_maxs
        iresult=0
      enddo
!
!--- max height(m) above ground where updraft air can originate
!
      zkbmax=3000.
!
!--- calculate moist static energy, heights, qes
!
      call cup_env(z,qes,he,hes,t,q,p,z1, &
           psur,ierr,tcrit,-1,   &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      call cup_env(zo,qeso,heo,heso,tn,qo,po,z1, &
           psur,ierr,tcrit,-1,   &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)

!
!--- environmental values on cloud levels
!
      call cup_env_clev(t,qes,q,he,hes,z,p,qes_cup,q_cup,he_cup, &
           hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur, &
           ierr,z1,          &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      call cup_env_clev(tn,qeso,qo,heo,heso,zo,po,qeso_cup,qo_cup, &
           heo_cup,heso_cup,zo_cup,po_cup,gammao_cup,tn_cup,psur,  &
           ierr,z1,          &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      do i=its,itf
        if(ierr(i).eq.0)then
        if(aaeq(i).lt.-0.1)then
           ierr(i)=20
        endif
!
      do k=kts,ktf
        if(zo_cup(i,k).gt.zkbmax+z1(i))then
          kbmax(i)=k
          go to 25
        endif
      enddo
 25   continue
!
      kbmax(i)=min(kbmax(i),ktf-4)
      endif
      enddo

!
!
!
!------- DETERMINE LEVEL WITH HIGHEST MOIST STATIC ENERGY CONTENT - K22
!
      k22(:)=2
      ktopx(:)=0
      CALL cup_MAXIMI(HEO_CUP,2,KBMAX,K22,ierr, &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
       DO 36 i=its,itf
         if(kpbl(i).gt.3)cap_max(i)=po_cup(i,kpbl(i))
         IF(ierr(I).eq.0.)THEN
         IF(K22(I).GT.KBMAX(i))then
           ierr(i)=2
           ierrc(i)="could not find k22"
           ktop(i)=0
           k22(i)=0
           kbcon(i)=0
         endif
         endif
 36   CONTINUE
!
!--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
!

      do i=its,itf
       IF(ierr(I).eq.0.)THEN
         if(use_excess == 2) then
             k1=max(1,k22(i)-1)
             k2=k22(i)+1
             hkb(i) =sum(he_cup(i,k1:k2))/float(k2-k1+1)
             hkbo(i)=sum(heo_cup(i,k1:k2))/float(k2-k1+1) !+xl*zqexec(i)+cp*ztexec(i)
             qkbo(i)=sum(qo_cup(i,k1:k2))/float(k2-k1+1)
!            write(0,*)sum(heo_cup(i,k1:k2))/float(k2-k1+1),heo_cup(i,k1),heo(i,k1:k2)
        else if(use_excess <= 1) then
             hkb(i)=he_cup(i,k22(i))
             hkbo(i)=heo_cup(i,k22(i)) !+float(use_excess)*(xl*zqexec(i)+cp*ztexec(i))
             if(k22(i).gt.kpbl(i))hkbo(i)=hkbo(i) !+xlv*zqexec(i)+cp*ztexec(i)
             qkbo(i)=qo_cup(i,k22(i))
        endif  ! excess
       endif ! ierr
      enddo
!JOE-Georg and Saulo's new idea:
      do i=its,itf
      do k=kts,ktf
          dbyo(i,k)=hkbo(i)-heso_cup(i,k)
      enddo
      enddo

      call cup_kbcon(ierrc,cap_max_increment,5,k22,kbcon,heo_cup,heso_cup, &
           hkbo,ierr,kbmax,po_cup,cap_max, &
           ztexec,zqexec,use_excess, &
           0,itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte, &
           z_cup,entr_rate,heo,0)
!
!
      DO i=its,itf
         entr_rate_2d(i,:)=entr_rate(i)
         IF(ierr(I).eq.0.)THEN
            if(kbcon(i).gt.ktf-4)then
                ierr(i)=231
            endif
            do k=kts,ktf
               frh = 2.*min(qo_cup(i,k)/qeso_cup(i,k),1.)
               entr_rate_2d(i,k)=entr_rate(i)*(2.3-frh)
               cd(i,k)=entr_rate_2d(i,k)
            enddo
         endif
      enddo
      call rates_up_pdf('shallow',ktop,ierr,po_cup,entr_rate_2d,hkbo,heo,heso_cup,zo_cup, &
           xland1,0,kstabi,k22,kbcon,its,ite,itf,kts,kte,ktf,zuo,kpbl,ktopx,ktopx)
      do i=its,itf
         if(ierr(i).eq.0)then
         do k=maxloc(zuo(i,:),1),1,-1 ! ktop(i)-1,1,-1
             if(zuo(i,k).lt.1.e-6)then
               k22(i)=k+1
               exit
             endif
         enddo
         if(k22(i).gt.1)then
         do k=1,k22(i)-1
           zuo(i,k)=0.
           zu(i,k)=0.
           xzu(i,k)=0.
         enddo
         endif
         do k=maxloc(zuo(i,:),1),ktop(i)
             if(zuo(i,k).lt.1.e-6)then
               ktop(i)=k-1
               exit
             endif
         enddo
         do k=k22(i),ktop(i)
          xzu(i,k)= zuo(i,k)
          zu (i,k)= zuo(i,k)
         enddo
         do k=ktop(i)+1,ktf
           zuo(i,k)=0.
           zu(i,k)=0.
           xzu(i,k)=0.
         enddo
         k22(i)=max(2,k22(i))
         endif
      enddo
!
! calculate mass entrainment and detrainment
!
      do k=kts,ktf
      do i=its,itf
         hc(i,k)=0.
         qco(i,k)=0.
         qrco(i,k)=0.
         DBY(I,K)=0.
         hco(i,k)=0.
         DBYo(I,K)=0.
         up_massentro(i,k)=0.
         up_massdetro(i,k)=0.
         up_massentr(i,k)=0.
         up_massdetr(i,k)=0.
      enddo
      enddo
      do i=its,itf
       IF(ierr(I).eq.0.)THEN
         do k=1,k22(i)
            hc(i,k)=he_cup(i,k)
            hco(i,k)=heo_cup(i,k)
            qco(i,k)=qo_cup(i,k)
         enddo

!        k=k22(i)
!        hc(i,k)=hkb(i)
!        qco(i,k)=qkbo(i)
!        hco(i,k)=hkbo(i)
       endif ! ierr
      enddo
!
!
      do 42 i=its,itf
         if(ierr(i).eq.0)then
         !- mass entrainment and detrinament is defined on model levels
         do k=2,maxloc(zuo(i,:),1)
         !=> below maximum value zu -> change entrainment
           dz=zo_cup(i,k)-zo_cup(i,k-1)
           up_massdetro(i,k-1)=cd(i,k-1)*dz*zuo(i,k-1)
           up_massentro(i,k-1)=zuo(i,k)-zuo(i,k-1)+up_massdetro(i,k-1)
           if(zuo(i,k-1).gt.0)entr_rate_2d(i,k-1)=(up_massentro(i,k-1))/(dz*zuo(i,k-1))
         enddo
         do k=maxloc(zuo(i,:),1)+1,ktop(i)
         !=> above maximum value zu -> change detrainment
           dz=zo_cup(i,k)-zo_cup(i,k-1)
           up_massentro(i,k-1)=entr_rate_2d(i,k-1)*dz*zuo(i,k-1)
           !-special treatment for ktop
           up_massdetro(i,k-1)=zuo(i,k-1)+up_massentro(i,k-1)-zuo(i,k)
!          write(0,*)k,dz,zuo(i,k-1),i
           if(zuo(i,k-1).gt.0.)cd(i,k-1)=up_massdetro(i,k-1)/(dz*zuo(i,k-1))
         enddo
         up_massdetro(i,ktop(i))=zuo(i,ktop(i))
         up_massentro(i,ktop(i))=0.
         do k=2,ktf-1
           up_massentr(i,k-1)=up_massentro(i,k-1)
          up_massdetr(i,k-1)=up_massdetro(i,k-1)
         enddo
         do k=kbcon(i)+1,ktf-1
           dbyo(i,k)=-1.
         enddo
         do k=k22(i)+1,ktop(i)
          hc(i,k)=(hc(i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k-1)*hc(i,k-1)+ &
                         up_massentr(i,k-1)*he(i,k-1))   /            &
                         (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
          dby(i,k)=hc(i,k)-hes_cup(i,k)
          hco(i,k)=(hco(i,k-1)*zuo(i,k-1)-.5*up_massdetro(i,k-1)*hco(i,k-1)+ &
                         up_massentro(i,k-1)*heo(i,k-1))   /            &
                         (zuo(i,k-1)-.5*up_massdetro(i,k-1)+up_massentro(i,k-1))
          dbyo(i,k)=hco(i,k)-heso_cup(i,k)
         enddo
         if(ktop(i).lt.kbcon(i)+1)then
            ierr(i)=5
            ierrc(i)='ktop is less than kbcon+1'
             go to 42
         endif
         if(ktop(i).gt.ktf-2)then
             ierr(i)=5
             ierrc(i)="ktop is larger than ktf-2"
             go to 42
         endif
         do k=k22(i)+1,ktop(i)
          trash=QESo_cup(I,K)+(1./XLV)*(GAMMAo_cup(i,k) &
              /(1.+GAMMAo_cup(i,k)))*DBYo(I,K)
          !- total water liq+vapour
          trash2  = qco(i,k-1) ! +qrco(i,k-1)
          qco (i,k)=   (trash2* ( zuo(i,k-1)-0.5*up_massdetr(i,k-1)) + &
                       up_massentr(i,k-1)*qo(i,k-1))   /            &
                       (zuo(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))

          if(qco(i,k)>=trash) then
              DZ=Z_cup(i,K)-Z_cup(i,K-1)
              ! cloud liquid water
              qrco(i,k)= (qco(i,k)-trash)/(1.+(c0_shal+c1_shal)*dz)
              pwo(i,k)=c0_shal*dz*qrco(i,k)*zuo(i,k)
              ! cloud water vapor 
              qco (i,k)= trash+qrco(i,k)
        
          else
              qrco(i,k)= 0.0
          endif 
          cupclw(i,k)=qrco(i,k)
         enddo
      do k=kts+1,ktop(i)
         dp=100.*(po_cup(i,k)-po_cup(i,k+1))
         cnvwt(i,k)=zuo(i,k)*cupclw(i,k)*g/dp
      enddo
 
         do k=ktop(i)+1,ktf-1
           HC(i,K)=hes_cup(i,k)
           HCo(i,K)=heso_cup(i,k)
           qco (i,k)=QESo_cup(I,K)
           qrco(i,k)=0.0
           DBY(I,K)=0.
           DBYo(I,K)=0.
           zu(i,k)=0.
           xzu(i,k)=0.
           zuo(i,k)=0.
           cd(i,k)=0.
           entr_rate_2d(i,k)=0.
           up_massentr(i,k)=0.
           up_massdetr(i,k)=0.
           up_massentro(i,k)=0.
           up_massdetro(i,k)=0.
         enddo
!        if(i.eq.ipr.and.j.eq.jpr)then
!           write(12,*)'hcnew = '
!           do k=1,ktop(i)
!             write(12,*)k,hco(i,k),dbyo(i,k)
!           enddo
!        endif
      endif
42    continue
!
!--- calculate workfunctions for updrafts
!
      call cup_up_aa0(aa0,z,zu,dby,GAMMA_CUP,t_cup, &
           kbcon,ktop,ierr,           &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      call cup_up_aa0(aa1,zo,zuo,dbyo,GAMMAo_CUP,tn_cup, &
           kbcon,ktop,ierr,           &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      do i=its,itf
         if(ierr(i).eq.0)then
           if(aa1(i).eq.0.)then
               ierr(i)=17
               ierrc(i)="cloud work function zero"
           endif
         endif
      enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!--- change per unit mass that a model cloud would modify the environment
!
!--- 1. in bottom layer
!
      do k=kts,ktf
      do i=its,itf
        dellah(i,k)=0.
        dellaq(i,k)=0.
      enddo
      enddo
!
!----------------------------------------------  cloud level ktop
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level ktop-1
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!
!----------------------------------------------  cloud level k+2
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level k+1
!
!----------------------------------------------  cloud level k+1
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level k
!
!----------------------------------------------  cloud level k
!
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!
!----------------------------------------------  cloud level 3
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level 2
!
!----------------------------------------------  cloud level 2
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level 1
      trashq=0.
      trash2=0.
      do i=its,itf
        if(ierr(i).eq.0)then
         do k=kts+1,ktop(i)
! entrainment/detrainment for updraft
            entup=up_massentro(i,k)
            detup=up_massdetro(i,k)
            totmas=detup-entup+zuo(i,k+1)-zuo(i,k)
            if(abs(totmas).gt.1.e-6)then
               write(0,*)'*********************',i,j,k,totmas
               write(0,*)k22(i),kbcon(i),ktop(i)
!              write(0,123)k,detup,entup,detupk,zuo(i,k+1),zuo(i,k)
!              123     format(1X,i2,5E12.4)
            !        call wrf_error_fatal ( 'totmas .gt.1.e-6' )
            endif
            dp=100.*(po_cup(i,k)-po_cup(i,k+1))
            dellah(i,k) =-(zuo(i,k+1)*(hco(i,k+1)-heo_cup(i,k+1) )-     &
                           zuo(i,k  )*(hco(i,k  )-heo_cup(i,k  ) ))*g/dp

            !-- take out cloud liquid water for detrainment
            dz=zo_cup(i,k+1)-zo_cup(i,k)
            if(k.lt.ktop(i))then
            dellaqc(i,k)= zuo(i,k)*c1_shal*qrco(i,k)*dz/dp*g !  detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp
            else
            dellaqc(i,k)=   detup*qrco(i,k) *g/dp
            endif

            !-- condensation source term = detrained + flux divergence of 
            !-- cloud liquid water (qrco)
            C_up = dellaqc(i,k)+(zuo(i,k+1)* qrco(i,k+1) -       &
                                 zuo(i,k  )* qrco(i,k  )  )*g/dp
            !-- water vapor budget (flux divergence of Q_up-Q_env - condensation
            !term)
            dellaq(i,k) =-(zuo(i,k+1)*(qco(i,k+1)-qo_cup(i,k+1) ) -      &
                           zuo(i,k  )*(qco(i,k  )-qo_cup(i,k  ) ) )*g/dp &
                           - C_up - 0.5*(pwo (i,k)+pwo (i,k+1))*g/dp
            !- check water conservation liq+vapor

!
          enddo   ! k

        endif
      enddo

!
!--- using dellas, calculate changed environmental profiles
!
      mbdt=3.e-4

      do k=kts,ktf
      do i=its,itf
         dellat(i,k)=0.
         if(ierr(i).eq.0)then
            XHE(I,K)=DELLAH(I,K)*MBDT+HEO(I,K)
!           XQ(I,K)=(DELLAQ(I,K)+dellaqc(i,k))*MBDT+QO(I,K)
!           DELLAT(I,K)=(1./cp)*(DELLAH(I,K)-xlv*(DELLAQ(I,K)+dellaqc(i,k)))
!           dellaq(i,k)=DELLAQ(I,K)+dellaqc(i,k)
!           dellaqc(i,k)=0.
            XQ(I,K)=max(1.e-16,(DELLAQ(I,K)+dellaqc(i,k))*MBDT+QO(I,K))
            DELLAT(I,K)=(1./cp)*(DELLAH(I,K)-xlv*(DELLAQ(I,K)))
            XT(I,K)= (-dellaqc(i,k)*xlv/cp+DELLAT(I,K))*MBDT+TN(I,K)
            xt(i,k)=max(190.,xt(i,k))
         ENDIF
      enddo
      enddo
      do i=its,itf
      if(ierr(i).eq.0)then
      xhkb(i)=hkbo(i)+(DELLAH(I,K22(i)))*MBDT
      XHE(I,ktf)=HEO(I,ktf)
      XQ(I,ktf)=QO(I,ktf)
      XT(I,ktf)=TN(I,ktf)
      endif
      enddo
!
!--- calculate moist static energy, heights, qes
!
      call cup_env(xz,xqes,xhe,xhes,xt,xq,po,z1, &
           psur,ierr,tcrit,-1,   &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
!
!--- environmental values on cloud levels
!
      call cup_env_clev(xt,xqes,xq,xhe,xhes,xz,po,xqes_cup,xq_cup, &
           xhe_cup,xhes_cup,xz_cup,po_cup,gamma_cup,xt_cup,psur,   &
           ierr,z1,          &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
!
!
!**************************** static control
      do k=kts,ktf
      do i=its,itf
         xhc(i,k)=0.
         xDBY(I,K)=0.
      enddo
      enddo
      do i=its,itf
        if(ierr(i).eq.0)then
         do k=1,kbcon(i)-1
            xhc(i,k)=xhkb(i)
         enddo
          k=kbcon(i)
          xhc(i,k)=xhkb(i)
          xDBY(I,Kbcon(i))=xHkb(I)-xHES_cup(I,K)
        endif !ierr
      enddo
!
!
      do i=its,itf
      if(ierr(i).eq.0)then
      xzu(i,1:ktf)=zuo(i,1:ktf)	!ss 2/19/14
      do k=kbcon(i)+1,ktop(i)
       xhc(i,k)=(xhc(i,k-1)*xzu(i,k-1)-.5*up_massdetro(i,k-1)*xhc(i,k-1)+ &
                         up_massentro(i,k-1)*xhe(i,k-1))   /            &
                         (xzu(i,k-1)-.5*up_massdetro(i,k-1)+up_massentro(i,k-1))
       xdby(i,k)=xhc(i,k)-xhes_cup(i,k)
      enddo
      do k=ktop(i)+1,ktf
           xHC(i,K)=xhes_cup(i,k)
           xDBY(I,K)=0.
           xzu(i,k)=0.
      enddo
      endif
      enddo

!
!--- workfunctions for updraft
!
      call cup_up_aa0(xaa0,xz,xzu,xdby,GAMMA_CUP,xt_cup, &
           kbcon,ktop,ierr,           &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
!
! now for shallow forcing
!
       do i=its,itf
        xmb(i)=0.
        xff_shal(1:9)=0.
        if(ierr(i).eq.0)then
          xmbmax(i)=1.0  
          xkshal=(xaa0(i)-aa1(i))/mbdt
          if(xkshal.ge.0.)xkshal=+1.e6
          if(xkshal.gt.-1.e-4 .and. xkshal.lt.0.)xkshal=-1.e-4
!         xff_shal(1)=max(0.,-(aa1(i)-aa0(i))/(xkshal*dtime))
!         xff_shal(1)=min(xmbmax(i),xff_shal(1))
!         xff_shal(2)=max(0.,-(aa1(i)-aa0(i))/(xkshal*dtime))
!         xff_shal(2)=min(xmbmax(i),xff_shal(2))
!         xff_shal(3)=max(0.,-(aa1(i)-aa0(i))/(xkshal*dtime))
!         xff_shal(3)=min(xmbmax(i),xff_shal(3))
!         if(aa1(i).le.0)then
!          xff_shal(1)=0.
!          xff_shal(2)=0.
!          xff_shal(3)=0.
!         endif
!         if(aa1(i)-aa0(i).le.0.)then
!          xff_shal(1)=0.
!          xff_shal(2)=0.
!          xff_shal(3)=0.
!         endif
           xff_shal(1)=.03*zws(i)
           xff_shal(2)=.03*zws(i)
           xff_shal(3)=.03*zws(i)
! boundary layer QE (from Saulo Freitas)
          blqe=0.
          trash=0.
!         if(k22(i).lt.kpbl(i)+1)then
             do k=1,kbcon(i)-1
                blqe=blqe+100.*dhdt(i,k)*(po_cup(i,k)-po_cup(i,k+1))/g
             enddo
             trash=max((hc(i,kbcon(i))-he_cup(i,kbcon(i))),1.e1)
             xff_shal(4)=max(0.,blqe/trash)
             xff_shal(4)=min(xmbmax(i),xff_shal(4))
             xff_shal(5)=xff_shal(4)
             xff_shal(6)=xff_shal(4)
             xmbs2(i)=xff_shal(4)
!         else
!         else
!            xff_shal(7)=0.
!         endif
!         if(xkshal.lt.-1.1e-04)then ! .and.  &
!            ((aa1(i)-aa0(i).gt.0.) .or. (xff_shal(7).gt.0)))then
!         xff_shal(4)=max(0.,-aa0(i)/(xkshal*tscl_KF))
!         xff_shal(4)=min(xmbmax(i),xff_shal(4))
!         xff_shal(5)=xff_shal(4)
!          xff_shal(4)=0.
!          xff_shal(5)=0.
!          xff_shal(6)=0.
!         endif
!         write(0,888)'i0=',i,j,kpbl(i),blqe,xff_shal(7)
!888       format(a3,3(1x,i3),2e12.4)
!         xff_shal(8)= xff_shal(7)
!         xff_shal(9)= xff_shal(7)
!         xff_shal(4)= xff_shal(7)
!         xff_shal(5)= xff_shal(1)
!         xff_shal(6)= (xff_shal(1)+xff_shal(7))*.5
!         fsum=0.
!         do k=1,9
!          xmb(i)=xmb(i)+xff_shal(k)
!          fsum=fsum+1.
!         enddo
          xmb(i)=max(xff_shal(2),xff_shal(5))
          xmb(i)=min(xmbmax(i),xmb(i))
          if(ichoice.gt.0)xmb(i)=min(xmbmax(i),xff_shal(ichoice))
          if(i.eq.ipr.and.j.eq.jpr)write(0,*)',ierr,xffs',ierr(i),xff_shal(1:9),xmb(i),xmbmax(i)
          if(xmb(i).eq.0.)ierr(i)=22
          if(xmb(i).eq.0.)ierrc(i)="22"
          if(xmb(i).lt.0.)then
             ierr(i)=21
             ierrc(i)="21"
!            write(0,*)'neg xmb,i,j,xmb for shallow = ',i,j,k22(i),ierr(i)
          endif
        endif
        if(ierr(i).ne.0)then
           k22(i)=0
           kbcon(i)=0
           ktop(i)=0
           xmb(i)=0
           do k=kts,ktf
              outt(i,k)=0.
              outq(i,k)=0.
              outqc(i,k)=0.
           enddo
        else if(ierr(i).eq.0)then
!
! got the mass flux, sanity check, first for heating rates
!
          trash=0.
!         kmaxx=0
!         do k=2,ktop(i)
!          trash=max(trash,86400.*(dellat(i,k))*xmb(i))
!         enddo
!         if(trash.gt.100.)then
!            xmb(i)=xmb(i)*100./trash
!         endif
!         trash=0.
!         do k=2,ktop(i)
!          trash=min(trash,86400.*(dellat(i,k))*xmb(i))
!         enddo
!         if(trash.lt.-100.)then
!             xmb(i)=-xmb(i)*100./trash
!         endif
!
! sanity check on moisture tendencies: do not allow anything that may allow neg
! tendencies
!
!         do k=2,ktop(i)
!          trash=q(i,k)+(dellaq(i,k))*xmb(i)*dtime
!         if(trash.lt.1.e-12)then
! max allowable tendency over tendency that would lead to too small mix ratios
!
!           trash=(1.e-12 -q(i,k))/((dellaq(i,k))*dtime)
!           xmb(i)=(1.e-12 -q(i,k))/((dellaq(i,k))*dtime)
!         endif
!         enddo
          xmb_out(i)=xmb(i)
! 
! final tendencies
!
          do k=2,ktop(i)
           outt(i,k)=(dellat(i,k))*xmb(i)
           outq(i,k)=(dellaq(i,k))*xmb(i)
!          outqc(i,k)= 0. !dellaqc(i,k)*xmb(i), dellaqc is in dellat
           outqc(i,k)= dellaqc(i,k)*xmb(i)
           pre(i)=pre(i)+xmb(i)*pwo(i,k)
          enddo
        endif
       enddo
!      
! done shallow
!--------------------------done------------------------------
!

   END SUBROUTINE CUP_gf_sh
 SUBROUTINE rates_up_pdf(name,ktop,ierr,p_cup,entr_rate_2d,hkbo,heo,heso_cup,z_cup, &
               xland,jpr,kstabi,k22,kbcon,its,ite,itf,kts,kte,ktf,zuo,kpbl,ktopdby,csum)
     implicit none
     integer, intent(in) :: jpr,its,ite,itf,kts,kte,ktf
     real, dimension (its:ite,kts:kte),intent (inout) :: entr_rate_2d,zuo
     real, dimension (its:ite,kts:kte),intent (in) ::p_cup, heo,heso_cup,z_cup
     real, dimension (its:ite),intent (in) :: hkbo
     integer, dimension (its:ite),intent (in) :: kstabi,k22,kbcon,kpbl,csum,xland
     integer, dimension (its:ite),intent (inout) :: ierr,ktop,ktopdby
     real, dimension (its:ite,kts:kte) :: hcot
      character *(*), intent (in)         ::                           &
       name
     real :: dz,dh, &
             dbythresh,dzh2
     real :: dby(kts:kte)
      real zuh2(40),zh2(40)
     integer :: kk,kbegin,i,k,ipr,kdefi,kstart,kbegzu,kfinalzu
     !integer :: lev_adj,power_entr
     dbythresh=0.6
     dbythresh=0.85
     dbythresh=0.95
     dby(:)=0.
!     dbythresh=1.
     if(name == 'shallow')dbythresh=1.
     if(name == 'mid')dbythresh=1.
     DO i=its,itf
      zuo(i,:)=0.
      dby(:)=0.
      ktop(i)=0
      if(ierr(i).eq.0)then
        hcot(i,kbcon(i))=hkbo(i)
        dz=z_cup(i,kbcon(i))-z_cup(i,kbcon(i)-1)
        dby(kbcon(i))=(hcot(i,kbcon(i))-heso_cup(i,kbcon(i)))*dz
        do k=kbcon(i)+1,ktf-2
           dz=z_cup(i,k)-z_cup(i,k-1)

           hcot(i,k)=( (1.-0.5*entr_rate_2d(i,k-1)*dz)*hcot(i,k-1) &
                      + entr_rate_2d(i,k-1)*dz*heo(i,k-1))/ &
                      (1.+0.5*entr_rate_2d(i,k-1)*dz)
           dby(k)=dby(k-1)+(hcot(i,k)-heso_cup(i,k))*dz
       enddo
!      do k=kbcon(i),ktf-2
       ktopdby(i)=maxloc(dby(:),1)
       do k=maxloc(dby(:),1)+1,ktf-2
!          if((hcot(i,k)-heso_cup(i,k)).lt.dbythresh)then
          if(dby(k).lt.dbythresh*maxval(dby))then
              kfinalzu=k  - 1
              ktop(i)=kfinalzu
              go to 412
          endif
       enddo
       kfinalzu=ktf-2
       ktop(i)=kfinalzu
412    continue
!
! at least overshoot by one level
!
       kfinalzu=min(max(kfinalzu,ktopdby(i)+1),ktopdby(i)+2)
       ktop(i)=kfinalzu
       if( name == 'deep' ) then
          if(kfinalzu.le.kbcon(i)+2)then
              ierr(i)=41
              ktop(i)= 0
          else
           call get_zu_zd_pdf_fim(xland(i),zuh2,"UP",ierr(i),k22(i),kfinalzu,zuo(i,kts:kte),kts,kte,ktf,0.,kpbl(i),csum(i))
          endif
       elseif ( name == 'mid' ) then
!         kfinalzu=min(kfinalzu,kstabi(i)+1)
!         ktopdby(i)=kfinalzu
!         ktop(i)=kfinalzu+1
          call get_zu_zd_pdf_fim(xland(i),zuh2,"MID",ierr(i),k22(i),kfinalzu,zuo(i,kts:kte),kts,kte,ktf,0.,kpbl(i),csum(i))
           kbegin=0
           dzh2=(z_cup(i,ktop(i))-z_cup(i,k22(i)))/40.
           zh2(1)=z_cup(i,k22(i))
           if(zuh2(1).gt.0.1 .and. kbegin.eq.0)kbegin=1
           do k=2,40
             zh2(k)=zh2(k-1)+dzh2
             if(zuh2(k).gt.0.1 .and. kbegin.eq.0)kbegin=k
           enddo
           zuo(i,k22(i))=zuh2(kbegin)
           do k=k22(i)+1,kfinalzu+1
             do kk=kbegin,39
             if(z_cup(i,k).gt.zh2(kk) .and.  z_cup(i,k).le.zh2(kk+1)) then
                 zuo(i,k)=zuh2(kk)
                 exit
              endif
             enddo
           enddo
           if(zuo(i,ktop(i)).lt.1.e-4)ktop(i)=ktop(i)-1
!          if(jpr.eq.1)then
!
!!sms$ignore begin
!           do k=kts,kfinalzu
!             write(12,*)k,zuo(i,k),zuh2(k)
!           enddo
!           call flush(12)
!!sms$ignore end
!          endif
       elseif ( name == 'shallow' ) then
          if(kfinalzu.le.kbcon(i)+1)then
              ierr(i)=41
              ktop(i)= 0
          else
           call get_zu_zd_pdf_fim(xland(i),zuh2,"SH3",ierr(i),k22(i),kfinalzu,zuo(i,kts:kte),kts,kte,ktf,0.,kpbl(i),csum(i))
           kbegin=0
           dzh2=(z_cup(i,kfinalzu)-z_cup(i,k22(i)))/40.
           zh2(1)=z_cup(i,k22(i))
           if(zuh2(1).gt.0.1 .and. kbegin.eq.0)kbegin=1
           do k=2,40
             zh2(k)=zh2(k-1)+dzh2
             if(zuh2(k).gt.0.1 .and. kbegin.eq.0)kbegin=k
           enddo
           zuo(i,k22(i))=zuh2(kbegin)
           do k=k22(i)+1,kfinalzu
             do kk=kbegin,39
             if(z_cup(i,k).gt.zh2(kk) .and.  z_cup(i,k).le.zh2(kk+1)) then
                 zuo(i,k)=zuh2(kk)
                 exit
              endif
             enddo
           enddo
         endif
      ENDIF
      ENDIF
     ENDDO

  END SUBROUTINE rates_up_pdf
!-------------------------------------------------------------------------
subroutine get_zu_zd_pdf_fim(xland,zuh2,draft,ierr,kb,kt,zu,kts,kte,ktf,max_mass,kpbli,csum)

implicit none
integer, intent(in) ::xland,kb,kt,kts,kte,ktf,kpbli,csum
real, intent(in) ::max_mass
real, intent(inout) :: zu(kts:kte)
real  :: zuh(kts:kte),zuh2(1:40)
integer, intent(inout) :: ierr
character*(*), intent(in) ::draft

!- local var
integer :: kk,add,i,nrec=0,k,kb_adj,kpbli_adj
real ::krmax,zumax,ztop_adj
real ::beta, alpha,kratio,tunning,FZU
!- kb cannot be at 1st level

   !-- fill zu with zeros
   zu=0.0
   zuh=0.0

IF(draft == "UP") then
   kb_adj=max(kb,2)
  !beta=4.  !=> must larger than 1
            !=> higher makes the profile sharper
            !=> around the maximum zu 
  !- 2nd approach for beta and alpha parameters
  !- the tunning parameter must be between 0.5 (low  level max zu)
  !-                                   and 1.5 (high level max zu)
!   tunning = 0.9 ! mar 1, first scores
!  tunning = 0.5
!  tunning = 1.
!  tunning = 0.6+min(0.8,float(csum)*.04)
  if(xland == 1 ) then
     tunning = 1.2-min(0.9,float(csum)*.01)
     beta    = 2. ! 2.0/tunning
  else
     tunning = 0.6
     beta    = 2. ! 2.0/tunning
  endif
  if(draft == "MID")then
   tunning = .6
   beta = 2.
  endif
  alpha   = tunning*beta
  fzu=1.
  do k=kts,min(kte,kt+1)
      kratio= float(k)/float(kt+1)
      zu(k) = FZU*kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
   enddo
   if(maxval(zu(kts:min(ktf,kt+1))).gt.0.)  &
      zu(kts:min(ktf,kt+1))= zu(kts:min(ktf,kt+1))/maxval(zu(kts:min(ktf,kt+1)))
     do k=maxloc(zu(:),1),1,-1
       if(zu(k).lt.1.e-6)then
         kb_adj=k+1
         exit
       endif
     enddo
     kb_adj=max(2,kb_adj)
     do k=kts,kb_adj-1
       zu(k)=0.
     enddo

ELSEIF(draft == "SH") then

  tunning = 0.6
  beta    =2.2/tunning
  if(draft == "MID")then
   tunning = .6
   beta = 2.
  endif
  alpha   = tunning*beta
  fzu=1.
  do k=kts,min(kt+1,ktf)
      kratio= float(k)/float(kt+1)
      zuh(k) = FZU*kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
   enddo
   do k=maxloc(zuh(:),1),1,-1
      kk=kb+k-maxloc(zuh(:),1)
      if(kk.gt.1)zu(kk)=zuh(k)
   enddo
   do k=maxloc(zuh(:),1)+1,kt
      kk=kb+k-maxloc(zuh(:),1)
      if(kk.le.kt)zu(kk)=zuh(k)
   enddo
   if(maxval(zu(kts:min(ktf,kt+1))).gt.0.)  &
      zu(kts:min(ktf,kt+1))= zu(kts:min(ktf,kt+1))/maxval(zu(kts:min(ktf,kt+1)))
ELSEIF(draft == "SH3" .or. draft == "MID") then
  tunning = 0.6
  beta    =2.2/tunning
  if(draft == "MID")then
   tunning = .6
   beta = 2.
  endif
  alpha   = tunning*beta
  fzu=1.
  do k=1,40
      kratio= float(k)/float(40)
      zuh2(k) = FZU*kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
   enddo
   if(maxval(zuh2(1:40)).gt.0.)  &
      zuh2(:)= zuh2(:)/ maxval(zuh2(1:40))
ELSEIF(draft == "SH2") then
   kb_adj=2 ! level where mass flux starts
   kpbli_adj=kpbli
   if(kpbli_adj < kb_adj) then
      print*,"kpbli < kb", kpbli_adj, kb_adj
      kpbli_adj = kb_adj + 1
   endif
   !- location of the maximum Zu
   krmax=float(kpbli_adj-kb_adj+1)/float(kt+1)
   !
   beta= 6.
   !- this alpha imposes the maximum zu at kpbli
   alpha=1.+krmax*(beta-1.)/(1.-krmax)

   !- to check if dZu/dk = 0 at k=kpbli_adj
   !kratio=krmax
   !dzudk=(alpha-1.)*(kratio**(alpha-2.)) * (1.-kratio)**(beta-1.) - &
   !      (kratio**(alpha-1.))*((1.-kratio)**(beta-2.))*(beta-1.)

   !- Beta PDF
   do k=kts+kb_adj-1,min(kte,kt+1)
      kratio=float(k+1-kb_adj)/float(kt+1)  !-kb_adj+1)
      zu(k)=kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
   enddo
     !do k=kts,min(kt+1,ktf)
     !  print*,"shallow2",k,zu(k),kt
     !- normalize ZU
     zu(kts:min(kte,kt+1))= zu(kts:min(kte,kt+1))/ maxval(zu(kts:min(kte,kt+1)))
     !enddo
  !   Zu(kt)=0.1
ELSEIF(draft == "DOWN" .or. draft == "DOWNM") then

! tunning = 0.8
! beta    = 3.0/tunning
  tunning = 0.8
  beta    =2.0/tunning
  alpha   = tunning*beta
  fzu=1.
  zuh(:)=0.
  zu(:)=0.
  do k=kts,min(kt+1,ktf)
      kratio= float(k)/float(kt+1)
      zuh(k) = FZU*kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
!      write(0,*)k,zuh(k)
   enddo
   if(maxloc(zuh(:),1).ge.kpbli)then
      do k=maxloc(zuh(:),1),1,-1
         kk=kpbli+k-maxloc(zuh(:),1)
         if(kk.gt.1)zu(kk)=zuh(k)
!         if(kk.gt.1)write(0,*)kk,zu(kk)
      enddo
      do k=maxloc(zuh(:),1)+1,kt
         kk=kpbli+k-maxloc(zuh(:),1)
         if(kk.le.kt)zu(kk)=zuh(k)
!         if(kk.ge.1)write(0,*)kk,zu(kk)
      enddo
   else
      do k=2,kt ! maxloc(zuh(:),1)
        zu(k)=zuh(k-1)
!         write(0,*)k,zu(k)
      enddo
   endif
   fzu=maxval(zu(kts:min(ktf,kt+1)))
   if(fzu.gt.0.)  &
      zu(kts:min(ktf,kt+1))= zu(kts:min(ktf,kt+1))/fzu
    if(zu(2).gt.max_mass)fzu=max_mass/zu(2) ! max(0.,zu(2)-max_mass)
     do k=2,kt+1
       zu(k)=fzu*zu(k)
     enddo
     zu(1)=0.


ENDIF
   !- normalize ZU
!  if(maxval(zu(kts:min(ktf,kt+1))).gt.0.)  &
!     zu(kts:min(ktf,kt+1))= zu(kts:min(ktf,kt+1))/ maxval(zu(kts:min(ktf,kt+1)))



return


end subroutine get_zu_zd_pdf_fim



!-------------------------------------------------------------------------
   SUBROUTINE cup_up_aa1bl(aa0,t,tn,q,qo,dtime, &
              z,zu,dby,GAMMA_CUP,t_cup,         &
              kbcon,ktop,ierr,                  &
              itf,jtf,ktf,                      &
              its,ite, jts,jte, kts,kte         )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,jtf,ktf,                                     &
        its,ite, jts,jte, kts,kte
  ! aa0 cloud work function
  ! gamma_cup = gamma on model cloud levels
  ! t_cup = temperature (Kelvin) on model cloud levels
  ! dby = buoancy term
  ! zu= normalized updraft mass flux
  ! z = heights of model levels 
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        z,zu,gamma_cup,t_cup,dby,t,tn,q,qo
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon,ktop
     real, intent(in) :: dtime
!
! input and output
!


     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        aa0
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k
     real                                 ::                           &
        dz,dA
!
        DO i=its,itf
         AA0(I)=0.
        ENDDO
        DO 100 k=kts+1,ktf
        DO 100 i=its,itf
         IF(ierr(i).ne.0 )GO TO 100
         IF(k.gt.KBCON(i))GO TO 100

         DZ=Z(I,K)-Z(I,K-1)
	 !print*,"dz=",i,k,z(i,k),Z(I,K-1),dz         
         !da=zu(i,k)*DZ*(9.81/(1004.*( &
         !        (T_cup(I,K)))))*DBY(I,K-1)/ &
         !     (1.+GAMMA_CUP(I,K))
         ! IF(K.eq.KTOP(I).and.da.le.0.)go to 100

	 dA=  DZ*9.81*( tn(i,k)-t(i,k) + 0.608*(qo(i,k)-q(i,k)))/dtime
         AA0(I)=AA0(I)+dA
100     CONTINUE

   END SUBROUTINE cup_up_aa1bl
!---------------------------------------------------------------------- 

END MODULE module_cu_gf
