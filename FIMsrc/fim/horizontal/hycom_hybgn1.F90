module hycom_hybgn1
   implicit none
contains

! --- ---------------------------
! --- hybrid ocean grid generator (based on Bleck 2002 paper)
! --- ---------------------------

   subroutine hybgn1(nstep,leap,targt,dens,dp,pres,temp,saln,uvel,vvel,	&
                     latij,passv_tr)

! --- hycom version 0.9.12
   use module_control  ,only: nip
   use module_constants,only: perm
   use fimnamelist    ,only: kdm,itest,diag_intvl
   use hycom_constants ,only: wet,r_onem
   use hycom_control   ,only: numtr

   integer,intent(IN) :: nstep			! model time step
   integer,intent(IN) :: leap			! leapfrog time slot
   real   ,intent(IN) :: targt(kdm)		! target densities

!SMS$DISTRIBUTE (dh,2) BEGIN
   real,intent(INOUT) :: uvel (kdm,nip,2)	! u velocity
   real,intent(INOUT) :: vvel (kdm,nip,2)	! v velocity
   real,intent(INOUT) :: dens (kdm,nip  )	! density
   real,intent(INOUT) :: dp   (kdm,nip,2)	! layer thickness
   real,intent(INOUT) :: pres (kdm+1,nip)	! interface pressure
   real,intent(INOUT) :: temp (kdm,nip)		! temperature
   real,intent(INOUT) :: saln (kdm,nip)		! salinity
   real,intent(INOUT),optional :: passv_tr(kdm,nip,numtr)	! passive tracer
!SMS$DISTRIBUTE END

!SMS$DISTRIBUTE (dh,1) BEGIN
   real,intent(IN)    :: latij(nip)		! latitude in degrees
   real dispmx(nip)
!SMS$DISTRIBUTE END

   real*8 dns1d(kdm),tem1d(kdm),sal1d(kdm),u1d(kdm),v1d(kdm),		&
          pr1d(kdm+1),trc1d(kdm,numtr)
   real   glbdisp
   logical :: vrbos,dotrcr
   logical,parameter :: useppm=.true.
   integer i,k,k1,nt,nwrk2,nwrk3,ntot2,ntot3

   dotrcr=.false.
   if (present(passv_tr)) dotrcr=.true.

!  print *,'entering hybgn1, step',nstep,' itest=',itest,' dotrcr=',dotrcr

   if (itest > 0) then
!SMS$PARALLEL (dh,i) BEGIN
    do i=1,nip
     vrbos=i==itest .and. wet(i)>0
     if (vrbos) then
      write (*,103) nstep,perm(i),					&
      ' ocn_hybgn1 IN:  temp   saln   dens  thkns   dpth     u      v',	&
      (k,temp(k,i),saln(k,i),dens(k,i),					&
      dp(k,i,leap)*r_onem,pres(k+1,i)*r_onem,uvel(k,i,leap)*100.,	&
      vvel(k,i,leap)*100.,k=1,kdm)
 103  format (2i7,a/(i28,3f7.3,f7.2,f7.1,2f7.2))

      write (*,99) nstep,perm(i),'   o l d   profile'
      do k=1,kdm,10
       write (*,100) (pres (k1,i)*r_onem,k1=k,min(kdm+1,k+10))
       write (*,101) (dens (k1,i),k1=k,min(kdm,k+9))
       write (*,101) (targt(k1),  k1=k,min(kdm,k+9))
       write (*,102) (temp (k1,i),k1=k,min(kdm,k+9))
       write (*,102) (saln (k1,i),k1=k,min(kdm,k+9))
      end do
 99   format (i9,i8,a,'  (5-line groups: dpth,dens,targt,T,S)')
 100  format (11f7.1)
 101  format (4x,10f7.2)
 102  format (4x,10f7.2)
     end if
    end do
!SMS$PARALLEL END
   end if			! itest > 0

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE(dns1d,tem1d,sal1d,pr1d,u1d,v1d,trc1d)		&
!$OMP SCHEDULE(runtime)							&
!$OMP REDUCTION(max:glbdisp) REDUCTION(+:nwrk2,nwrk3,ntot2,ntot3)
   do 2 i=1,nip
   if (wet(i) > 0) then
    vrbos=i==itest

! --- extract t,s,rho column from 3-d grid 
    pr1d(1)=pres(1,i)
    do k=1,kdm
     dns1d(k)=dens(k,i)
     tem1d(k)=temp(k,i)
     sal1d(k)=saln(k,i)
     u1d(k)=uvel(k,i,leap)
     v1d(k)=vvel(k,i,leap)
     if (dotrcr) then
      do nt=1,numtr
       trc1d(k,nt)=passv_tr(k,i,nt)
      end do
     end if
     pr1d(k+1)=pr1d(k)+dp(k,i,leap)
    end do

    call hybgn1_1d(nstep,kdm,targt,pr1d,tem1d,sal1d,dns1d,u1d,v1d,	&
                   trc1d,numtr,dotrcr,ntot2,ntot3,nwrk2,nwrk3,		&
                   dispmx(i),latij(i),vrbos,perm(i))

! --- put 1-d column back into 3-d grid
    do k=1,kdm
     dens(k,i)=dns1d(k)
     temp(k,i)=tem1d(k)
     saln(k,i)=sal1d(k)
     uvel(k,i,leap)=u1d(k)
     vvel(k,i,leap)=v1d(k)
     dp(k,i,leap)=pr1d(k+1)-pr1d(k)
     if (dotrcr) then
      do nt=1,numtr
       passv_tr(k,i,nt)=trc1d(k,nt)
      end do
     end if
     pres(k+1,i)=pr1d(k+1)
    end do

    glbdisp=dispmx(i)
   end if				! ocean point
 2 continue
!$OMP END PARALLEL DO
!SMS$PARALLEL END

   if (itest > 0) then
!SMS$PARALLEL (dh,i) BEGIN
    do i=1,nip
     vrbos=i==itest .and. wet(i)>0
     if (vrbos) then
      write (*,99) nstep,perm(i),'   n e w   profile'
      do k=1,kdm,10
       write (*,100) (pres (k1,i)*r_onem,k1=k,min(kdm+1,k+10))
       write (*,101) (dens (k1,i),k1=k,min(kdm,k+9))
       write (*,101) (targt(k1),  k1=k,min(kdm,k+9))
       write (*,102) (temp (k1,i),k1=k,min(kdm,k+9))
       write (*,102) (saln (k1,i),k1=k,min(kdm,k+9))
      end do

      write (*,103) nstep,perm(i),					&
      ' ocn_hybgn1 OUT: temp   saln   dens  thkns   dpth     u      v',	&
      (k,temp(k,i),saln(k,i),dens(k,i),					&
      dp(k,i,leap)*r_onem,pres(k+1,i)*r_onem,uvel(k,i,leap)*100.,	&
      vvel(k,i,leap)*100.,k=1,kdm)
     end if
    end do
!SMS$PARALLEL END
   end if			! itest > 0

   if (mod(nstep,diag_intvl).eq.0) then

!SMS$REDUCE (glbdisp,MAX)
    if (glbdisp.gt.0.) then
!SMS$PARALLEL (dh,i) BEGIN
     do i=1,nip
      if (wet(i) > 0) then
       if (dispmx(i).eq.glbdisp) then
        print '(i7,a,f7.1,a,i8,i3)',nstep,' (hybgn1) max.displ',	&
          glbdisp*r_onem,'m  at i=',perm(i)
       end if
      end if
     end do
!SMS$PARALLEL END
    end if

!SMS$REDUCE (nwrk2,nwrk3,ntot2,ntot3,SUM)
    write (*,'(a,f6.1,a,i9,a)') 'hybgen - grid restoration at',		&
      100.*float(nwrk3)/float(ntot3),' per cent of',ntot3,' points'
    write (*,'(a,f6.1,a,i9,a)') 'hybgen - new bottom layer at',		&
      100.*float(nwrk2)/float(ntot2),' per cent of',ntot2,' points'
   end if
!  print *,'...exiting hybgn1'
   return
   end subroutine hybgn1


   subroutine hybgn1_1d(nstep,kdm,targt,pr1d,tem1d,sal1d,dns1d,		&
                        u1d,v1d,trc1d,numtr,dotrcr,ntot2,ntot3,		&
                        nwrk2,nwrk3,dispmx,lati,vrbos,i)

! --- ------------------------------------------
! --- 1-D version of hybrid ocean grid generator (based on Bleck 2002 paper)
! --- ------------------------------------------

   use hycom_control   ,only: bclin_frq
   use hycom_constants ,only: wet,epsil,onem,tencm,onecm,onemm,r_onem,	&
                              salmin,dpmin,huuge,batrop
   use hycom_sigetc

   integer,intent(IN)    :: nstep,kdm,numtr,i
   integer,intent(OUT)   :: ntot2,ntot3,nwrk2,nwrk3
   real,   intent(OUT)   :: dispmx
   real*8, intent(INOUT) :: dns1d(kdm),tem1d(kdm),sal1d(kdm),		&
                            u1d(kdm),v1d(kdm),trc1d(kdm,numtr),		&
                            pr1d(kdm+1)
   real,   intent(IN)    :: lati,targt(kdm)		! latitude, target dens
   logical,intent(IN)    :: vrbos,dotrcr

   real*8 delp,dp0,dp0abv,dpsum,phi,plo,dsgdt,dsgds,slak,q1,q2,		&
     rho_lo,tem_lo,sal_lo,rho_up,tem_up,sal_up,tem,sal,p_hat,q,		&
     coltemo,coltemn,colsalo,colsaln,coltrco,coltrcn,			&
     scale,displ(kdm+1),sumrho,sumtem,sumsal,try,dlt,glbdisp
   real*8 uold(kdm),vold(kdm),trold(kdm,numtr),pold(kdm+1),		&
          told(kdm),sold(kdm),unew(kdm),vnew(kdm),trnew(kdm,numtr)
   logical :: event
   logical,parameter :: useppm=.true.
!  real cushn
!  external cushn
   integer k,k1,k2,kp,nt
   character info*16
   real,parameter :: tfreez=-2.2			! abyssal freezing temp
   real,parameter :: sigjmp=.01, acurcy=1.e-10
   real,parameter :: scalt=-30.,scals=10.		! oceanic t/s range
   integer,parameter :: miglim=-1	! intfc migration limit (<0: no limit)

!  parameter (slak=.5/86400.)	! intfc nudging time scale: 2 days
   parameter (slak=1./86400.)	! intfc nudging time scale: 1 day
!  parameter (slak=2./86400.)	! intfc nudging time scale: 12 hrs
!  parameter (slak=4./86400.)	! intfc nudging time scale: 6 hrs

! --- linear taper functions (latitude and depth-dependent) for slak
   real*8 tapr,wgtf
   real slakf,z
   tapr(q)=1.+9.*max(0.,1.-.01*r_onem*q)                ! q = pressure (Pa)
!  wgtf(z)=max(0.,min(1.,(abs(z)-50.)*.1))              ! 0->1 for z=50->60
!  slakf(z)=min(0.7,max( tapr(p_hat)*slak*dlt,          ! z = latitude (deg)
! .    0.7*wgtf(z)+tapr(p_hat)*slak*dlt*(1.-wgtf(z))))
!  slakf(z)=0.7*wgtf(z)+tapr(p_hat)*slak*dlt*(1.-wgtf(z))  ! z = lat.(deg)
   slakf(z)=tapr(p_hat)*slak*dlt                        ! no lat.dependence
!  slakf(z)=slak*dlt                                    ! no lat/dpth depend'ce

! --- slak function for 'no slak':
!  slakf(z)=min(0.7,max( tapr(p_hat),			& ! z = latitude (deg)
!     0.7*wgtf(z)+tapr(p_hat)*(1.-wgtf(z))))

   dlt=batrop*bclin_frq
   do k=1,kdm
    told(k)=tem1d(k)
    sold(k)=sal1d(k)
    uold(k)=u1d(k)
    vold(k)=v1d(k)
    if (dotrcr) then
     do nt=1,numtr
      trold(k,nt)=trc1d(k,nt)
     end do
    end if
   end do

   ntot2=0
   ntot3=0
   nwrk2=0
   nwrk3=0

   coltemo=0.
   colsalo=0.
   do k=1,kdm
    coltemo=coltemo+tem1d(k)*(pr1d(k+1)-pr1d(k))
    colsalo=colsalo+sal1d(k)*(pr1d(k+1)-pr1d(k))
   end do

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- do convective adjustment to eliminate static instabilities

!cc 30 event=.false.
!cc   do k=kdm/5,kdm-1
!cc    dpsum=pr1d(k+1)-pr1d(k)
!cc    if (dpsum.gt.0.) then
!cc     sumrho=dns1d(k)*dpsum
!cc      sumtem=tem1d(k)*dpsum
!cc     sumsal=sal1d(k)*dpsum
!cc     do k1=k+1,kdm+1
!cc      if (k1.gt.kdm. or. (pr1d(k1+1).gt.pr1d(k1) .and.		&
!cc        sumrho.lt.(dns1d(k1)+.001*sigjmp)*dpsum)) exit
!cc      if (vrbos) write (*,'(i8,2(a,i3))') i,			&
!cc        '  entrain layer',k1,'  into layer',k
!cc      delp=pr1d(k1+1)-pr1d(k1)
!cc      if (delp.gt.0.) then
!cc       event=.true.
!cc       dpsum=dpsum+delp
!cc       sumrho=sumrho+dns1d(k1)*delp
!cc       sumtem=sumtem+tem1d(k1)*delp
!cc       sumsal=sumsal+sal1d(k1)*delp
!cc      end if
!cc     end do
!cc     if (event) then
!cc      if (vrbos) write (*,'(i8,3x,2(a,i3),a,f9.3)') i,		&
!cc        'replacing rho in layers',k,' -',k1-1,'  by',sumrho/dpsum
!cc      sal1d(k)=sumsal/dpsum
!cc      tem1d(k)=sumtem/dpsum
!cc      dns1d(k)=sigocn(tem1d(k),sal1d(k))
!cc      do k2=k+1,k1-1
!cc       dns1d(k2)=dns1d(k)
!cc       tem1d(k2)=tem1d(k)
!cc       sal1d(k2)=sal1d(k)
!cc      end do
!cc      go to 30
!cc     end if
!cc    end if
!cc   end do
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! --- absorb near-massless layers on sea floor in layer above

   kp=1
   do 4 k=2,kdm
!  if (pr1d(k).lt.pr1d(kdm+1)-onecm) then
   if (pr1d(k).lt.pr1d(kdm+1)-tencm) then
    kp=k
   else
    q1=max(epsil,pr1d(k  )-pr1d(kp))
    q2=max(   0.,pr1d(k+1)-pr1d(k ))
    q=q1/(q1+q2)
    if (q.lt.0. .or. q.gt.1.)						&
      write (*,*) 'error hybgn1 - i,q1,q2,q=',i,q1,q2,q
    tem1d(kp)=tem1d(kp)*q+tem1d(k)*(1.-q)
    sal1d(kp)=sal1d(kp)*q+sal1d(k)*(1.-q)
    tem1d(k)=tem1d(kp)
    sal1d(k)=max(sal1d(kp),salmin(k))
    if (dotrcr) then
     do nt=1,numtr
      trold(kp,nt)=trold(kp,nt)*q+trold(k,nt)*(1.-q)
      trold(k,nt)=trold(kp,nt)
     end do
    end if

    if (vrbos) then
     write (*,'(i9,i8,a,i3,a,i3,5x,a,3f7.3)') nstep,i,			&
      '  absorb lyr',k,' in',kp,'new t,s,rho:',tem1d(kp),sal1d(kp),	&
     sigocn(tem1d(kp),sal1d(kp))
    end if

   end if
 4 continue

   dns1d(kp)=sigocn(tem1d(kp),sal1d(kp))
   do k=kp+1,kdm
    dns1d(k)=dns1d(kp)
    pr1d(k)=pr1d(kdm+1)
   end do

   coltemn=0.
   colsaln=0.
   do k=1,kdm
    coltemn=coltemn+tem1d(k)*(pr1d(k+1)-pr1d(k))
    colsaln=colsaln+sal1d(k)*(pr1d(k+1)-pr1d(k))
   end do

   if (abs(coltemn-coltemo).gt.acurcy*10.*pr1d(kdm+1))			&
     write (*,104) i,'  ocn_hybgn1 - bad temp.intgl',coltemo,		&
     coltemn,(coltemn-coltemo)/(10.*pr1d(kdm+1))
   if (abs(colsaln-colsalo).gt.acurcy*35.*pr1d(kdm+1))			&
     write (*,104) i,'  ocn_hybgn1 - bad saln.intgl',colsalo,		&
     colsaln,(colsaln-colsalo)/(35.*pr1d(kdm+1))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- check whether density values exceed range of 'targt' values
! --- homogenize layers to eliminate this condition

!c   do k=kdm,2,-1
!c    if (dns1d(k).le.targt(kdm)+.01*sigjmp) then
!c     exit				!  no action required
!c    else
!c c --- average of layers k...kdm denser than densest target
!c     if (vrbos)							&
!c      write (*,108) i,'lyrs',k-1,kdm,' bfore retargeting:',		&
!c       (pr1d(k2)*r_onem,dns1d(k2),k2=k-1,kdm),pr1d(kdm+1)*r_onem
!c c --- compute density obtained by homogenizing layers k-1,...,kdm
!c     q1=max(   0.,pr1d(k  )-pr1d(k-1))
!c     q2=max(epsil,pr1d(k+1)-pr1d(k  ))
!c     q=1./(q1+q2)
!c     q1=q1*q
!c     q2=q2*q
!c     try=dns1d(k-1)*q1+dns1d(k)*q2
!c     if (try.gt.targt(kdm)-.01*sigjmp) then
!c c --- average of layers k-1,...,kdm still too dense. continue adding layers
!c      pr1d(k)=pr1d(kdm+1)
!c      tem1d(k-1)=tem1d(k-1)*q1+tem1d(k)*q2
!c      sal1d(k-1)=sal1d(k-1)*q1+sal1d(k)*q2
!c      dns1d(k-1)=sigocn(tem1d(k-1),sal1d(k-1))
!c      if (dotrcr) then
!c       trold(k-1,:)=trold(k-1,:)*q1+trold(k,:)*q2
!c      end if
!c      do k1=k,kdm
!c       dns1d(k1)=dns1d(k-1)
!c       tem1d(k1)=tem1d(k-1)
!c       sal1d(k1)=sal1d(k-1)
!c       if (dotrcr) then
!c        trold(k1,:)=trold(k-1,:)
!c       end if
!c      end do
!c      if (vrbos)							&
!c       write (*,108) i,'lyrs',k-1,kdm,' after retargeting:',		&
!c        (pr1d(k2)*r_onem,dns1d(k2),k2=k-1,kdm),pr1d(kdm+1)*r_onem
!c     else
!c c --- adding all of layer k-1 is overkill. entrain only part of layer k-1
!c      p_hat=max(pr1d(k-1),min(pr1d(k+1),				&
!c       (pr1d(k+1)*(targt(kdm)-dns1d(k  ))				&
!c       +pr1d(k  )*(  dns1d(k)-dns1d(k-1)))				&
!c                /(targt(kdm)-dns1d(k-1))))
!c      q1=max(   0.,pr1d(k  )-p_hat)
!c      q2=max(epsil,pr1d(k+1)-pr1d(k))
!c      q=1./(q1+q2)
!c      q1=q1*q
!c      q2=q2*q
!c      pr1d(k)=p_hat
!c      tem1d(k)=tem1d(k-1)*q1+tem1d(k)*q2
!c      sal1d(k)=sal1d(k-1)*q1+sal1d(k)*q2
!c      dns1d(k)=sigocn(tem1d(k-1),sal1d(k-1))
!c      if (dotrcr) then
!c       trold(k,:)=trold(k-1,:)*q1+trold(k,:)*q2
!c      end if
!c      do k1=k+1,kdm
!c       dns1d(k1)=dns1d(k)
!c       tem1d(k1)=tem1d(k)
!c       sal1d(k1)=sal1d(k)
!c       if (dotrcr) then
!c        trold(k1,:)=trold(k,:)
!c       end if
!c      end do
!c      if (vrbos)							&
!c       write (*,108) i,'lyrs',k-1,kdm,' after retargeting:',
!c       (pr1d(k2)*r_onem,dns1d(k2),k2=k-1,kdm),pr1d(kdm+1)*r_onem
!c      exit
!c     end if
!c    end if
!c   end do
!c 108 format (i8,2x,a,i3,'-',i2,a,25(f7.1,f8.4))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- is layer touching sea floor to light?
   if (kp.eq.1) go to 10
   k=kp
   if (dns1d(k).le.max(dns1d(k-1),targt(k-1)) .or.			&
       dns1d(k).ge.targt(k)) go to 10

! --- split layer -kp- into 2 sublayers such that density in lower sublayer
! --- matches target density of lyr -kp-. combine upper sublayer with lyr kp-1.

   tem=tem1d(k)
   sal=sal1d(k)
   dsgdt=dsigdt(tem,sal)*scalt
   dsgds=dsigds(tem,sal)*scals
   q=1./(dsgdt+dsgds)

! --- set properties in lower sublayer:
   sal_lo=min(sal+(targt(k)-dns1d(k))*q*scals,max(sal1d(k-1),		&
              sal+.01*scals))
   tem_lo=max(tem+(targt(k)-dns1d(k))*q*scalt,tfreez)
   rho_lo=sigocn(tem_lo,sal_lo)

! --- set properties in upper sublayer:
   rho_up=max(dns1d(k-1),min(targt(k-1),dns1d(k)))
   q=(rho_lo-dns1d(k))/max(epsil,rho_lo-rho_up)
   if (q.lt.0. .or. q.gt.1.) go to 10
   p_hat=pr1d(k)*(1.-q)+pr1d(k+1)*q
   q=(pr1d(k+1)-p_hat)/max(epsil,p_hat-pr1d(k))
   sal_up=sal*(1.+q)-sal_lo*q
   tem_up=tem*(1.+q)-tem_lo*q

   if (vrbos) then
    write (*,'(i9,i8,i3,a,3f7.3,f8.2)') nstep,i,k,			&
      '  t,s,th,dp in upper sblyr:',tem_up,sal_up,			&
      sigocn(tem_up,sal_up),(p_hat-pr1d(k))*r_onem
    write (*,'(20x,a,3f7.3,f8.2)') 					&
      '  t,s,th,dp in lower sblyr:',tem_lo,sal_lo,			&
      rho_lo,(pr1d(k+1)-p_hat)*r_onem
    write (*,'(20x,a,2es11.3)') '  scalt,scals =',scalt,scals
   end if

! --- combine upper sublayer with layer k-1
   q=(p_hat-pr1d(k))/max(p_hat-pr1d(k-1),epsil)
   if (q.lt.-0.00000001 .or. q.gt.1.00000001) 				&
    write (*,*) 'q out of range - i,k,p_hat,pr1d(k),q=',		&
      i,k,p_hat,pr1d(k),q
   sal1d(k-1)=sal_up*q+sal1d(k-1)*(1.-q)
   tem1d(k-1)=tem_up*q+tem1d(k-1)*(1.-q)
   if (dotrcr) then
    do nt=1,numtr
     trold(k-1,nt)=trold(k,nt)*q+trold(k-1,nt)*(1.-q)
    end do
   end if

   if (vrbos) then
    write (*,'(22x,a,2f7.3)')						&
     '  old/new th(k-1):',dns1d(k-1),sigocn(tem1d(k-1),sal1d(k-1)),	&
     '  old/new th(k  ):',dns1d(k  ),rho_lo
   end if

   dns1d(k-1)=sigocn(tem1d(k-1),sal1d(k-1))
   tem1d(k)=tem_lo
   sal1d(k)=sal_lo
   dns1d(k)=rho_lo
   pr1d(k)=p_hat
   nwrk2=nwrk2+1
10 continue
   ntot2=ntot2+1

   do k=1,kdm+1
    pold(k)=pr1d(k)
   end do

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- try to restore isopycnic conditions by moving layer interfaces
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   dpsum=0.
   dp0=huuge
   do 8 k=1,kdm
   ntot3=ntot3+1

! --- set lower limits for layer thknss (dp0) and depth of lower intfc (dpsum)

   dp0abv=dp0
   dp0=dpmin(k)*onem

! --- optional: reduce spacing of z layers near equator, but hide transition
! --- in a subtropical latitude band where z layers are least likely to exist
   if (k.gt.1) dp0=dp0*max(.6,min(1.,(abs(lati)+5.)*.04))

! --- reduce layer thickness in shallow spots, creating sigma coord. effect
   if (4*k.lt.kdm) dp0=dp0*min(1.,pr1d(kdm+1)/(200.*onem)+.4)
   dpsum=dpsum+dp0

! --- maintain constant thickness in layer 1
   if (k.eq.1) then
    p_hat=dp0
    if (p_hat.gt.pr1d(2)) then
! --- layer 1 is too thin. entrain water from layers below
     p_hat=min(p_hat,pr1d(2)+						&
      max(onecm,slakf(lati)*(p_hat-pr1d(2))))
     info='layer too thin  '
     go to 5
    else if (p_hat.lt.pr1d(2)) then
! --- layer 1 is too thick. expell layer 1 water into layer 2
     p_hat=max(p_hat,pr1d(2)+						&
       min(-onecm,slakf(lati)*(p_hat-pr1d(2))))
     info='layer too thick '

     if (vrbos) then
      write (*,105) i,k,info,						&
       'lower intfc',pr1d(k+1)*r_onem,'=>',p_hat*r_onem
105  format (i8,i3,2x,a,'  try moving',2(1x,a,f9.3))
     end if

     q=(pr1d(2)-p_hat)/max(pr1d(3)-p_hat,epsil)
     if (q.lt.0. .or. q.gt.1.) then
      write (*,*) 'i,k,pr1d(2),p_hat,q=',i,k,pr1d(2),p_hat,q
     end if
     sal1d(2)=sal1d(2)*(1.-q)+sal1d(1)*q
     tem1d(2)=tem1d(2)*(1.-q)+tem1d(1)*q
     dns1d(2)=sigocn(tem1d(2),sal1d(2))

     if (vrbos) then
      write (*,107) i,k,'detrain into layer below: lower intfc',	&
       pr1d(k+1)*r_onem,'=>',p_hat*r_onem
107  format (i8,i3,1x,2(1x,a,f9.3))
     end if

     pr1d(2)=p_hat
     nwrk3=nwrk3+1
    end if
    go to 8
   end if				!  k = 1

! --- k > 1  hereafter

! --- are we dealing with a near-massless layer on the sea floor?
   if (pr1d(k).eq.pr1d(kdm+1)) dns1d(k)=max(targt(k),dns1d(k))
   if (pr1d(k).gt.pr1d(kdm+1)-onecm) go to 8

! --- is lower intfc too close to the surface?
   p_hat=dpsum
   if (k.lt.kdm .and. p_hat.gt.pr1d(k+1)) then
    p_hat=min(p_hat,pr1d(k+1)+						&
      max(onecm,slakf(lati)*(p_hat-pr1d(k+1))))
    info='too close to srf'
    go to 5
   end if

! --- is density noticeably different from target value?
   if (abs(dns1d(k)-targt(k)).lt..1*sigjmp) then
!c   if (vrbos) print '(i8,i3,a,2f8.3)',i,k,				&
!c       '  actual density matches target:',dns1d(k),targt(k)
    go to 8
   end if

!c if (vrbos) print '(i8,i3,a,2f8.3)',i,k,				&
!c   '  actual versus target density: ',dns1d(k),targt(k)
   if (dns1d(k).le.targt(k)) go to 7		!  layer too light

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- water in layer k is too  d e n s e . dilute with water from layer k-1
!                              ^^^^^^^^^
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   if (k.eq.2) go to 6			!  don't touch layer 1
   q=(targt(k)-dns1d(k))/max(targt(k)-dns1d(k-1),sigjmp*10.)
   p_hat=pr1d(k)*(1.-q)+pr1d(k+1)*q

! --- maintain minimum thickness of layer k-1
!! p_hat=pr1d(k-1)+cushn(p_hat-pr1d(k-1),dp0abv)
   p_hat=max(p_hat,pr1d(k-1)+dp0abv)				! no cushn
   p_hat=min(p_hat,.5*(pr1d(k-1)+pr1d(k+1)))
   info='layer too dense '

   if (vrbos) then
    write (*,105) i,k,info,'upper intfc',			&
     pr1d(k)*r_onem,'=>',p_hat*r_onem
   end if

   if (p_hat.lt.pr1d(k)) then

! --- upper intfc moves up. entrain layer k-1 water into layer k

    p_hat=max(p_hat,pr1d(k-1),pr1d(k)+					&
      min(-onecm,slakf(lati)*(p_hat-pr1d(k))))
    if (useppm .and. 							&
      abs(dns1d(k-1)-targt(k-1)).gt..1*sigjmp) then		!  use ppm
     displ(1)=0.
     displ(2)=0.
     displ(3)=p_hat-pr1d(k)
     displ(4)=0.
     call ppmad3(pr1d(k-2),displ,sal1d(k-2),i)
     call ppmad3(pr1d(k-2),displ,tem1d(k-2),i)
     dns1d(k-1)=sigocn(tem1d(k-1),sal1d(k-1))
     dns1d(k  )=sigocn(tem1d(k  ),sal1d(k  ))
    else							!  use pcm
     q=(pr1d(k)-p_hat)/max(pr1d(k+1)-p_hat,epsil)
     if (q.lt.0. .or. q.gt.1.) 						&
      write (*,*) 'i,k,p_hat,pr1d(k),q=',i,k,p_hat,pr1d(k),q
     sal1d(k)=sal1d(k)*(1.-q)+sal1d(k-1)*q
     tem1d(k)=tem1d(k)*(1.-q)+tem1d(k-1)*q
     dns1d(k)=sigocn(tem1d(k),sal1d(k))
    end if

    if (vrbos) then
     write (*,107) i,k,'entrain from layer above: upper intfc',		&
      pr1d(k)*r_onem,'=>',p_hat*r_onem
    end if

    pr1d(k)=p_hat
    nwrk3=nwrk3+1

   else if (p_hat.gt.pr1d(k)) then		!  p_hat > pr1d(k)

    if (vrbos) then
     print '(i8,i3,a)',i,k,'  lyr can''t expand upward'
    end if

! --- layer k-1 is too thin for allowing upper intfc to move up.  instead,
! --- move upper interface down and entrain layer k water into layer k-1

    p_hat=min(p_hat,pr1d(k+1),pr1d(k)+					&
      max(onecm,slakf(lati)*(p_hat-pr1d(k))))
    if (useppm .and. k.lt.kdm) then				!  use ppm
     displ(1)=0.
     displ(2)=p_hat-pr1d(k)
     displ(3)=0.
     displ(4)=0.
     call ppmad3(pr1d(k-1),displ,sal1d(k-1),i)
     call ppmad3(pr1d(k-1),displ,tem1d(k-1),i)
     dns1d(k-1)=sigocn(tem1d(k-1),sal1d(k-1))
     dns1d(k  )=sigocn(tem1d(k  ),sal1d(k  ))
    else							!  use pcm
     q=(p_hat-pr1d(k))/max(p_hat-pr1d(k-1),epsil)
     if (q.lt.0. .or. q.gt.1.) then
      write (*,*) 'i,k,p_hat,pr1d(k),q=',i,k,p_hat,pr1d(k),q
     end if
     sal1d(k-1)=sal1d(k-1)*(1.-q)+sal1d(k)*q
     tem1d(k-1)=tem1d(k-1)*(1.-q)+tem1d(k)*q
     dns1d(k-1)=sigocn(tem1d(k-1),sal1d(k-1))
    end if

    if (vrbos) then
     write (*,107) i,k,'detrain into layer above: upper intfc',		&
      pr1d(k)*r_onem,'=>',p_hat*r_onem
    end if

    pr1d(k)=p_hat
    nwrk3=nwrk3+1
   end if

! --- do we need to inflate layer k by lowering  l o w e r  interface?

 6 p_hat=pr1d(k)+dpmin(2)*onem
   if (k.lt.kdm .and. pr1d(k+1).lt.pr1d(kdm+1)-onemm .and.		&
    pr1d(k+1).lt.p_hat) then
    p_hat=min(p_hat,pr1d(k+1)+						&
      max(onecm,slakf(lati)*(p_hat-pr1d(k+1))))
    info='layer too thin  '
    go to 5
   else

    if (vrbos) then
     print '(i8,i3,a)',i,k,'  lyr can''t expand dnward'
    end if

   end if
   go to 8

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- water in layer k is too  l i g h t . dilute with water from layer k+1
!                              ^^^^^^^^^
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 7 if (k.ge.kdm .or. pr1d(k+1).gt.pr1d(kdm+1)-onemm) then
    if (vrbos) print '(i8,i3,a)',i,k,'  lyr can''t expand dnward'
    go to 8
   end if

   q=(dns1d(k)-targt(k))/max(dns1d(k+1)-targt(k),sigjmp*10.)
   p_hat=pr1d(k+1)*(1.-q)+pr1d(k)*q

! --- curtail downward growth of layers (esp. lowest hybrid layer)
   p_hat=min(p_hat,max(.5*(pr1d(k)+pr1d(k+2)),2.*pr1d(k+1)-pr1d(k)))

   p_hat=min(p_hat,pr1d(k+1)+						&
      max(onecm,slakf(lati)*(p_hat-pr1d(k+1))))
   info='layer too light '

 5 p_hat=min(p_hat,pr1d(k+2))
   if (p_hat.gt.pr1d(k+1)+onemm) then

    if (vrbos) then
     write (*,105) i,k,info,'lower intfc',				&
      pr1d(k+1)*r_onem,'=>',p_hat*r_onem
    end if

    if (useppm .and. k.lt.kdm-1 .and.					&
      abs(dns1d(k+1)-targt(k+1)).gt..1*sigjmp) then		!  use ppm
     displ(1)=0.
     displ(2)=p_hat-pr1d(k+1)
     displ(3)=0.
     displ(4)=0.
     call ppmad3(pr1d(k),displ,sal1d(k),i)
     call ppmad3(pr1d(k),displ,tem1d(k),i)
     dns1d(k  )=sigocn(tem1d(k  ),sal1d(k  ))
     dns1d(k+1)=sigocn(tem1d(k+1),sal1d(k+1))
    else							!  use pcm
     q=(p_hat-pr1d(k+1))/max(p_hat-pr1d(k),epsil)
     if (q.lt.0. .or. q.gt.1.) then
      write (*,*) 'i,k,p_hat,pr1d(k+1),q=',i,k,p_hat,pr1d(k+1),q
     end if
     sal1d(k)=sal1d(k)*(1.-q)+sal1d(k+1)*q
     tem1d(k)=tem1d(k)*(1.-q)+tem1d(k+1)*q
     dns1d(k)=sigocn(tem1d(k),sal1d(k))
    end if

    if (vrbos) then
     write (*,107) i,k,'entrain from layer below: lower intfc',		&
     pr1d(k+1)*r_onem,'=>',p_hat*r_onem
    end if

    pr1d(k+1)=p_hat
    nwrk3=nwrk3+1
   end if
 8 continue				! main k loop

! --- suppress irregularities in interface spacing 
   call equilb(nstep,pr1d,dns1d,tem1d,sal1d,kdm,vrbos,i)

   dispmx=0.
   do k=1,kdm
    dispmx=max(dispmx,abs(pr1d(k+1)-pold(k+1)))
   end do

   do k=1,kdm-2
    if (pr1d(k+1).lt.pr1d(k) .or. pr1d(k+1).gt.pr1d(k+2)) then
     write (*,'(a,3i5)') 'trcr monotonicity problems at',i,k
     write (*,'(a/(8f9.0))') 'pold:',pold
     write (*,'(a/(8f9.0))') 'pnew:',pr1d
!    stop '(monotonicity problem)'
     pr1d(k+1)=max(pr1d(k),min(pr1d(k+1),pold(k+2)))
    end if
   end do

! --- evaluate effect of regridding on -u,v-

   call pcwise(1,kdm,pold,uold,pr1d,unew,miglim,vrbos,i)
   call pcwise(1,kdm,pold,vold,pr1d,vnew,miglim,vrbos,i)

   if (dotrcr) then

! --- evaluate effect of regridding on tracer field(s)

    do nt=1,numtr
     scale=1.e-33
     coltrco=0.
     do k=1,kdm
      coltrco=coltrco+trold(k,nt)*(pold(k+1)-pold(k))
      scale=max(scale,abs(trold(k,nt)))
     end do

     call pcwise(3,kdm,pold,trold(1,nt),pr1d,trnew(1,nt),miglim,	&
                 vrbos,i)

     coltrcn=0.
     do k=1,kdm
      coltrcn=coltrcn+trnew(k,nt)*(pr1d(k+1)-pr1d(k))
     end do

!    if (vrbos) then
      if (abs(coltrcn-coltrco).gt.acurcy*scale*pr1d(kdm+1))		&
       write (*,104) i,'  ocn_hybgn1 - bad trcr.intgl',coltrco,		&
       coltrcn,(coltrcn-coltrco)/(scale*pr1d(kdm+1))
!    end if

    end do				!  numtr
   end if				!  dotrcr

   coltemn=0.
   colsaln=0.
   do k=1,kdm
    coltemn=coltemn+tem1d(k)*(pr1d(k+1)-pr1d(k))
    colsaln=colsaln+sal1d(k)*(pr1d(k+1)-pr1d(k))
   end do

!  if (vrbos) then
    if (abs(coltemn-coltemo).gt.acurcy*10.*pr1d(kdm+1))		&
     write (*,104) i,'  ocn_hybgn1 - bad temp.intgl',coltemo,	&
     coltemn,(coltemn-coltemo)/(10.*pr1d(kdm+1))
    if (abs(colsaln-colsalo).gt.acurcy*35.*pr1d(kdm+1))		&
     write (*,104) i,'  ocn_hybgn1 - bad saln.intgl',colsalo,	&
     colsaln,(colsaln-colsalo)/(35.*pr1d(kdm+1))
 104  format (i8,a,2es15.7,es9.1)
!  end if

! --- reduce salt content wherever density > targt(kdm)

!  do k=1,kdm
!   if (pr1d(k).lt.pr1d(kdm+1)) then
!    q=sofsig(targt(kdm),tem1d(k))
!    if (q.lt.sal1d(k)) then
!    if (vrbos .or. sal1d(k)-q.gt..2)					&
!      write (*,'(i9,i8,i3,4(a,f7.3))') nstep,i,k,			&
!     ' reduce s=',sal1d(k),' ->',q,' to reduce rho=',			&
!      dns1d(k),' ->',targt(kdm)
!    sal1d(k)=q
!    dns1d(k)=targt(kdm)
!    end if
!   end if
!  end do

   do k=1,kdm
    u1d(k)=unew(k)
    v1d(k)=vnew(k)
    if (dotrcr) then
     do nt=1,numtr
      trc1d(k,nt)=trnew(k,nt)
     end do
    end if
   end do

   return
   end subroutine hybgn1_1d


   subroutine ppmad3(x,dx,y,i)

! --- advection by piecewise parabolic method
! --- this is a special recursive version taylored to nmax = 3
!                       ^^^^^^^^^
! --- input variables:
! --- y(nmax)    - function values at cell midpoints
! --- x(nmax+1)  - cell boundaries
! --- dx(nmax+1) - displacement of cell boundaries during 1 time step
! ---              (= neg.flow velocity x time step)

! --- output variables:
! --- y(nmax)    - function values after advection (overwriting input)

   integer,parameter     :: nmax=3,ndim=3
   real*8 ,intent(IN)    :: x(nmax+1),dx(nmax+1)
   real*8 ,intent(INOUT) :: y(nmax)
   integer,intent(IN)    :: i
   real*8 :: total,tndcy,scale,dxlft,dxmid,dxrgt,wdth,slab,		&
          yold(ndim),ytmp(ndim),dxnew,a,b,c,yl,yr
   integer :: n
   logical :: diagno=.false.
   real,parameter :: acurcy=1.e-11 ,onemu=.098
   real,parameter :: athird=1./3.

   total=0.
   scale=1.e-33
   do 3 n=1,nmax
   if (x(n).gt.x(n+1)+onemu) then
!SMS$IGNORE BEGIN
    write (*,'(i8,a,4f9.1)') i,' ppmad3 error: x not monotonic',x
!SMS$IGNORE END
!   stop '(ppmad3)'
   end if

   yold(n)=y(n)
   ytmp(n)=y(n)*(x(n+1)-x(n))
   total=total+ytmp(n)
 3 scale=scale+abs(ytmp(n))
   scale=scale/float(nmax)

! --- construct parabola whose integral over [-.5,+.5] equals y(2) and
! --- which passes though points yl,yr at [-.5,+.5] resp.

   dxlft=max(onemu,x(2)-x(1))
   dxmid=max(onemu,x(3)-x(2))
   dxrgt=max(onemu,x(4)-x(3))
   yl=(dxlft*yold(2)+dxmid*yold(1))/(dxlft+dxmid)
   yr=(dxrgt*yold(2)+dxmid*yold(3))/(dxrgt+dxmid)

   a=1.5*yold(2)-.25*(yl+yr)
   b=yr-yl
   c=6.*(.5*(yl+yr)-yold(2))
   if (abs(yr-yl) .lt. 6.*abs(.5*(yl+yr)-yold(2))) then

! --- apex of parabola lies inside interval [-.5,+.5].
! --- => need to change curve to prevent over- and undershoots

    if (abs(yr-yl) .gt. 2.*abs(.5*(yl+yr)-yold(2))) then
! --- put apex of parabola on edge of interval [-.5,+.5]
     if ((yr-yl)*(.5*(yl+yr)-yold(2)) .gt. 0.) then
! --- apex at x=-.5
      a=.25*(3.*yold(2)+yl)
      c=3.*(yold(2)-yl)
      b=c
     else
! --- apex at x=+.5
      a=.25*(3.*yold(2)+yr)
      c=3.*(yold(2)-yr)
      b=-c
     end if
    else			!  -1/6 < x < +1/6
! --- can't put apex on edge of interval. only option is to flatten curve
     a=yold(2)
     b=0.
     c=0.
    end if
   end if

! --- now transport -y- across cell boundaries

   if (dx(2).gt.0.) then
! --- velocity at left edge is negative. export slab to cell on left
    wdth= dx(2)/(x(3)-x(2))
    slab=(x(3)-x(2))*							&
      wdth*(a+b*.5*(wdth-1.)+c*(.25-wdth*(.5-wdth*athird)))
    ytmp(2)=ytmp(2)-slab
    ytmp(1)=ytmp(1)+slab
   end if
   if (dx(3).lt.0.) then
! --- velocity at right edge is positive. export slab to cell on right
    wdth=-dx(3)/(x(3)-x(2))
    slab=(x(3)-x(2))*							&
      wdth*(a+b*.5*(1.-wdth)+c*(.25-wdth*(.5-wdth*athird)))
    ytmp(2)=ytmp(2)-slab
    ytmp(3)=ytmp(3)+slab
   end if

   tndcy=-total
   do 5 n=1,nmax
   dxnew=x(n+1)-x(n)+dx(n+1)-dx(n)
   if (dxnew.gt.0.) y(n)=ytmp(n)/dxnew
 5 tndcy=tndcy+y(n)*dxnew
   if (abs(tndcy).gt.acurcy*max(abs(total),scale*x(nmax+1))) then
!SMS$IGNORE BEGIN
    write (*,'(i8,a,2es11.3)') i,' ppmad3 - bad intgl',		&
     total,tndcy
!SMS$IGNORE END
   end if

   if (diagno) then
!SMS$IGNORE BEGIN
    write (*,100) 'ppmad3:',					&
     (n,x(n),yold(n),dx(n),y(n),n=1,nmax)
100 format (a/(i3,f14.1,es17.7,f14.1,es17.7))
!SMS$IGNORE END
   end if

   return
   end subroutine ppmad3


   subroutine pcwise(remap_optn,kk,xold,fldold,xnew,fldnew,miglim,vrbos,ipn)

!-- 1-dim PCM/PLM/PPM transport routine, extracted from HYCOM's fct3d.f.
!-- if miglim.ne.1, |xold-xnew| can exceed cell width (i.e. no CFL constraints)

!-- there are 3 advection/remapping choices:
!-- remap_optn=1			! PCM (piecewise constant)
!-- remap_optn=2			! PLM (piecewise linear)
!-- remap_optn=3			! PPM (piecewise parabolic)

   integer  ,intent(IN)  :: remap_optn
   integer  ,intent(IN)  :: kk		! vert.dimension
   integer  ,intent(IN)  :: miglim	! intfc migration limit (<0: no limit)
   logical  ,intent(IN)  :: vrbos	!  switch for 'verbose' mode
   integer  ,intent(IN)  :: ipn	! grid cell number

!-- xold/new	- old/new cell boundaries
!-- fldold/new	- mixing ratio of dep.variable before and after transport

   real*8 ,intent(IN)  :: xold(kk+1),xnew(kk+1),fldold(kk)
   real*8 ,intent(OUT) :: fldnew(kk)

   real*8 zold(kk+1),znew(kk+1),delx(kk+1),delz(kk+1),fco(kk),fcn(kk),	&
        vertfx(kk+1),vertdv(kk)
   real*8 a(kk),b(kk),c(kk),dx,fcdx,yl,yr
   real*8 amount,bfore,after,dpth,scale,slab,dslab
   integer k,lyr
   real, parameter :: athird=1./3.
!  real, parameter :: small=1.e-9
   real, parameter :: small=1.e-11
   real, parameter :: acurcy=1.e-11

   delx=xnew(:)-xold(:)
!-- make sure -x- increases with -k-. change sign of -xold/new- if necessary
   if (xold(1) < xold(kk+1)) then
     zold(:)=xold(:)
     znew(:)=xnew(:)
   else
     zold(:)=-xold(:)
     znew(:)=-xnew(:)
   end if
   delz=znew(:)-zold(:)

   if (vrbos) then
!SMS$IGNORE BEGIN
     write (*,100) ipn,							&
     'entering pcwise:   old_indep    d(indep)   depnd.var',		&
     (k,xold(k),delx(k),fldold(k),k=1,kk),kk+1,xold(kk+1),delx(kk+1)
100 format (i8,3x,a/(i27,2f12.2,es12.3))
!SMS$IGNORE END
    end if

!-- deduce old and new cell width from -zold,znew-
   do 15 k=1,kk
   fco(k)=max(0.,zold(k+1)-zold(k))
15 fcn(k)=max(0.,znew(k+1)-znew(k))

   bfore=0.
   scale=1.e-33
   dpth=0.
   do k=1,kk
     bfore=bfore+fldold(k)*fco(k)
     dpth=dpth+fco(k)
     scale=scale+abs(fldold(k))
   end do
   fldnew=fldold

!-- start by filling zero-width cells with data from neighboring cells

   do 17 k=kk-1,1,-1
17 fldnew(k)=(fldnew(k)*fco(k)+fldnew(k+1)*small)			&
            /(          fco(k)+            small)
   do 18 k=2,kk
18 fldnew(k)=(fldnew(k)*fco(k)+fldnew(k-1)*small)			&
            /(          fco(k)+            small)

!-- fit 0th, 1st, or 2nd deg. polynomial to -fldnew- in each cell
   a(1 )=fldnew(1 )
   b(1 )=0.
   c(1 )=0.
   a(kk)=fldnew(kk)
   b(kk)=0.
   c(kk)=0.

   do 16 k=2,kk-1

   if (remap_optn.eq.1) then
!-- piecewise constant method:
     a(k)=fldnew(k)
     b(k)=0.
     c(k)=0.

   else if (remap_optn.eq.2) then
!-- piecewise linear method:
!-- fit linear function a+bx to tracer in each cell (-.5 < x < +.5)
     a(k)=fldnew(k)
     b(k)=0.
     if (fldnew(k) <= min(fldnew(k-1),fldnew(k+1)) .or.			&
         fldnew(k) >= max(fldnew(k-1),fldnew(k+1))) then
       b(k)=0.
     else if ((fldnew(k+1)-fldnew(k-1))*(fldnew(k-1)+fldnew(k+1)	&
       -2.*fldnew(k)) > 0.) then
       b(k)=fldnew(k)-fldnew(k-1)
     else
       b(k)=fldnew(k+1)-fldnew(k)
     end if
     c(k)=0.

   else if (remap_optn.eq.3) then
!-- piecewise parabolic method:
!-- construct parabola a+bx+cx^2  whose integral over [-.5,+.5] equals
!-- fldnew(k) and which passes though points yl,yr at [-.5,+.5] resp.
!!   yl=.5*(fldnew(k-1)+fldnew(k))
!!   yr=.5*(fldnew(k+1)+fldnew(k))
     yl=(max(small,fco(k-1))*fldnew(k)+max(small,fco(k))*fldnew(k-1))/	&
        (max(small,fco(k-1))          +max(small,fco(k)))
     yr=(max(small,fco(k+1))*fldnew(k)+max(small,fco(k))*fldnew(k+1))/	&
        (max(small,fco(k+1))          +max(small,fco(k)))
     a(k)=1.5*fldnew(k)-.25*(yl+yr)
     b(k)=yr-yl
     c(k)=6.*(.5*(yl+yr)-fldnew(k))
     if (abs(yr-yl) < 6.*abs(.5*(yl+yr)-fldnew(k))) then
!-- apex of parabola lies inside interval [-.5,+.5], implying an over-
!-- or undershoot situation. change curve to prevent over/undershoots.
       if (abs(yr-yl) > 2.*abs(.5*(yl+yr)-fldnew(k))) then
!-- put apex of parabola on edge of interval [-.5,+.5]
         if ((yr-yl)*(.5*(yl+yr)-fldnew(k)) > 0.) then
!-- apex at x=-.5
           a(k)=.25*(3.*fldnew(k)+yl)
           c(k)=3.*(fldnew(k)-yl)
           b(k)=c(k)
         else
!-- apex at x=+.5
           a(k)=.25*(3.*fldnew(k)+yr)
           c(k)=3.*(fldnew(k)-yr)
           b(k)=-c(k)
         end if
       else			!  -1/6 < x < +1/6
!-- moving apex won't help. replace parabola by constant.
         a(k)=fldnew(k)
         b(k)=0.
         c(k)=0.
       end if
     end if
   else
!SMS$IGNORE BEGIN
     print *,'subr.pcwise -- illegal remap_optn=',remap_optn
!SMS$IGNORE END
     stop '(pcwise)'
   end if
16 continue

!-- get flux by summing -fldnew- over upstream slab of thickness -delz-

   if (miglim.eq.1) then
     do 27 k=2,kk
     slab=0.
     amount=0.
     vertfx(k)=0.
     if (delz(k) > 0.) then			! interface moves in +k dir.
       lyr=k
       if (fco(lyr) > 0.) then
         dslab=min(slab+fco(lyr), delz(k))	&
              -min(slab         , delz(k))
         dx=dslab/fco(lyr)
         fcdx=a(lyr)				&
             +b(lyr)*.5*(dx-1.)			& !  not needed in pcm
             +c(lyr)*(.25-dx*(.5-dx*athird))	  !  not needed in pcm,plm
         amount=amount+fcdx*dslab
         slab=slab+dslab
       end if
     else if (delz(k) < 0.) then		! interface moves in -k dir.
       lyr=k-1
       if (fco(lyr) > 0.) then
         dslab=min(slab+fco(lyr),-delz(k))	&
              -min(slab         ,-delz(k))
         dx=dslab/fco(lyr)
         fcdx=a(lyr)				&
             +b(lyr)*.5*(1.-dx)			& !  not needed in pcm
             +c(lyr)*(.25-dx*(.5-dx*athird))	  !  not needed in pcm,plm
         amount=amount+fcdx*dslab
         slab=slab+dslab
       end if
     end if
     if (slab.ne.0.) vertfx(k)=-delz(k)*amount/slab
27   continue

   else				! miglim .ne. 1
     do 22 k=2,kk
     slab=0.
     amount=0.
     vertfx(k)=0.
     if (delz(k) > 0.) then			! interface moves in +k dir.
       lyr=k-1
24     lyr=lyr+1
       if (slab >= delz(k)) goto 23
       if (fco(lyr) > 0.) then
         dslab=min(slab+fco(lyr), delz(k))	&
              -min(slab         , delz(k))
         dx=dslab/fco(lyr)
         fcdx=a(lyr)				&
             +b(lyr)*.5*(dx-1.)			& !  not needed in pcm
             +c(lyr)*(.25-dx*(.5-dx*athird))	  !  not needed in pcm,plm
         amount=amount+fcdx*dslab
         slab=slab+dslab
       end if
       if (lyr < kk) go to 24
     else if (delz(k) < 0.) then		! interface moves in -k dir.
       lyr=k
25     lyr=lyr-1
       if (slab >= -delz(k)) goto 23
       if (fco(lyr) > 0.) then
         dslab=min(slab+fco(lyr),-delz(k))	&
              -min(slab         ,-delz(k))
         dx=dslab/fco(lyr)
         fcdx=a(lyr)				&
             +b(lyr)*.5*(1.-dx)			& !  not needed in pcm
             +c(lyr)*(.25-dx*(.5-dx*athird))	  !  not needed in pcm,plm
         amount=amount+fcdx*dslab
         slab=slab+dslab
       end if
       if (lyr > 2) go to 25
     end if
23   if (slab.ne.0.) vertfx(k)=-delz(k)*amount/slab
22   continue
   end if

   vertfx(   1)=0.			!  don't allow flux through lower bdry
   vertfx(kk+1)=0.			!  don't allow flux through upper bdry
   do 26 k=1,kk
26 vertdv(k)=vertfx(k+1)-vertfx(k)

   if (vrbos) then
!SMS$IGNORE BEGIN
    write (*,'(a/(i3,4es12.3))')					&
     'pcwise:   flux  flx.div/thk    old_thk     new_thk',		&
     (k,vertfx(k),vertdv(k)/max(small,fcn(k)),fco(k),fcn(k),k=1,kk),	&
     kk+1,vertfx(kk+1)
!SMS$IGNORE END
   end if

   do 4 k=1,kk
   amount=fldnew(k)*fco(k)-vertdv(k)
4  fldnew(k)=(fldnew(k)*small+amount)/(small+fcn(k))

   after=0.
   do k=1,kk
     after=after+fldnew(k)*fcn(k)
   end do

   if (abs(bfore-after)*kk > acurcy*scale*dpth) 			&
     write (*,104) ipn,'pcwise - bad column intgl.:',bfore,after
104  format (i8,3x,a,2es15.7)

   if (vrbos) then
!SMS$IGNORE BEGIN
    write (*,100) ipn,							&
     'exiting pcwise:     d(indep)   new_indep   depnd.var',		&
     (k,delx(k),xnew(k),fldnew(k),k=1,kk),kk+1,delx(kk+1),xnew(kk+1)
!SMS$IGNORE END
   end if

   return
   end subroutine pcwise


   real function cushn(delp,dp0)

! --- c u s h i o n   function (from Bleck & Benjamin, 1993):

!     /  1        for x < 2-x1
!     |
!     |   (x + x1 - 2)**2
! --- cushn = <   1 + ---------------  for 2-x1 < x < x1
!     |  4 ( x1 - 1)
!     |
!     \  x        for x > x1

! --- if x = delp/dp0 >>  0, cushn*dp0 returns -delp-
! --- if x = delp/dp0 <<  0, cushn*dp0 returns -dp0-

   real delp,dp0,qq,x1,factor
!cc   parameter (x1=4.)     !  used in Bleck&Benjamin 1993
!cc   parameter (x1=6.)     !  used in Bleck 2002
!cc   parameter (x1=8.)
   parameter (x1=6.828427125)  !  x1=4+sqrt(8) yields cushn(0)=2

   parameter (factor=.25/(x1-1.))

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!cc   qq=max(-1.,min(2.,delp/(2.*dp0)))
!cc   cushn=dp0*(1.+athird*(qq+1.)**2)
!cc  .   *max(1.,delp/(2.*dp0))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!cc   qq=max(-.75,min(1.25,delp/(4.*dp0)))
!cc   cushn=dp0*(1.+(qq+.75)**2)
!cc  .   *max(1.,delp/(5.*dp0))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!cc   qq=max(-1.,min(1.5,delp/(4.*dp0)))
!cc   cushn=dp0*(1+.8*(qq+1.)**2)
!cc  .   *max(1.,delp/(6.*dp0))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   qq=max(2.-x1,min(x1,delp/dp0))
   cushn=dp0*(1.+factor*(qq+x1-2.)**2) * max(1.,delp/(x1*dp0))
   return
   end function cushn


   subroutine equilb(nstep,pres,thet,temp,saln,kdm,vrbos,i)

! --- expand thin layers at the expense of thick layers above and below.
! --- do this without changing targt

   use hycom_constants ,only: tenm,onem,tencm,onemm,onemu,epsil

   integer,intent(IN) :: nstep		! model time step
   integer,intent(IN) :: kdm,i		! no. of layers, test point location
   logical,intent(IN) :: vrbos		! if true, print results at test point
   real*8,intent(INOUT) :: pres(kdm+1)	! Exner fcn on interfaces
   real*8,intent(INOUT) :: thet(kdm)	! pot.density in layers
   real*8,intent(INOUT) :: temp(kdm)	! temperature in layers
   real*8,intent(INOUT) :: saln(kdm)	! salinity in layers
   integer k,k1,k2,ncount,iter
   real*8 dp1,dp2,dp3,dp4,dp5,th1,th2,th3,th4,dis2,dis3,dis4,		&
     ratio,goal,pold(kdm+1),pnew(kdm+1),pwrk(kdm+1),heatfx(kdm),	&
     plimit1,plimit2,cnv,told(kdm),sold(kdm),tdz(kdm),sdz(kdm),		&
     twrk(kdm),swrk(kdm),dtdz,dsdz,tbfor,taftr,sbfor,saftr,delp
   logical event
   real,parameter :: acurcy=1.e-11
   real,parameter :: slak=.25	! retardation coefficient
   real,parameter :: dffudt=.1	! therm.diffu.coeff x time step [m^2]

   ncount=0
   pold(:)=pres(:)
   pnew(:)=pres(:)
   pwrk(:)=pres(:)

   tbfor=0.
   sbfor=0.
   do k=1,kdm
    told(k)=temp(k)
    sold(k)=saln(k)
    twrk(k)=temp(k)
    swrk(k)=saln(k)
    tbfor=tbfor+(pres(k+1)-pres(k))*temp(k)
    sbfor=sbfor+(pres(k+1)-pres(k))*saln(k)
    tdz(k)=max(epsil,pres(k+1)-pres(k))*temp(k)
    sdz(k)=max(epsil,pres(k+1)-pres(k))*saln(k)
   end do

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- invoke heat diffusion (McDougall-Dewar) to inflate thin layers

   heatfx(kdm)=0.
   heatfx(  1)=0.
   do iter=1,3				! apply the scheme repeatedly

    do k=2,kdm-1
     heatfx(k)=0.
     if (pwrk(k  ).gt.pwrk(    1)+onem .and.			&
         pwrk(k+1).lt.pwrk(kdm+1)-onem) then
      cnv=onem				! pressure units per meter
!     heatfx(k)=dffudt*cnv*cnv*.5*(thet(k+1)-thet(k-1))			&
      heatfx(k)=dffudt*cnv*cnv*min(thet(k+1)-thet(k),thet(k)-thet(k-1))	&
       /max(.01*cnv,pwrk(k+1)-pwrk(k))
      plimit1=max(pwrk(k  )-tenm,.5*(pwrk(k  )+pwrk(k-1)))
      plimit2=min(pwrk(k+1)+tenm,.5*(pwrk(k+1)+pwrk(k+2)))
      heatfx(k)=min(max(0.,heatfx(k)),				&
       max(0.,(plimit2-pwrk(k+1))*(thet(k+1)-thet(k  ))),	&
       max(0.,(pwrk(k  )-plimit1)*(thet(k  )-thet(k-1))))
     end if
    end do

    do k=1,kdm-1
     if (thet(k+1).gt.thet(k)) then
      pnew(k+1)=pwrk(k+1)+(heatfx(k)-heatfx(k+1))/(thet(k+1)-thet(k))
      delp=pnew(k+1)-pwrk(k+1)
      if (delp.gt.0.) then
       dtdz=twrk(k+1)*delp
       dsdz=swrk(k+1)*delp
      else		! delp < 0
       dtdz=twrk(k  )*delp
       dsdz=swrk(k  )*delp
      end if
      tdz(k  )=tdz(k  )+dtdz
      tdz(k+1)=tdz(k+1)-dtdz
      sdz(k  )=sdz(k  )+dsdz
      sdz(k+1)=sdz(k+1)-dsdz
     end if
    end do

    event=.false.
    do k=2,kdm+1
     if (pnew(k).lt.pnew(k-1)-onemu) then
      event=.true.
!SMS$IGNORE BEGIN
      print '(i8,i3,a)',i,k,' (equilb) dp<0 due to heat diffusion'
!SMS$IGNORE END
     end if
    end do

    do k=1,kdm
     if (pnew(k+1).gt.pnew(k)) then
      twrk(k)=tdz(k)/(pnew(k+1)-pnew(k))
      swrk(k)=sdz(k)/(pnew(k+1)-pnew(k))
     end if
    end do

    if (vrbos .or. event) then
!SMS$IGNORE BEGIN
     print 99,i,' heat diffusion',				&
      '  (6-line groups: old/new p,pchg,rho,T,S): iter',iter
     do k2=1,kdm,9
      write (6,100) (pwrk(k1),k1=k2,min(kdm+1,k2+9) )
      write (6,100) (pnew(k1),k1=k2,min(kdm+1,k2+9) )
      write (6,102) (int(pnew(k1)-pwrk(k1)),k1=k2,min(kdm+1,k2+9) )
      write (6,101) (thet(k1),k1=k2,min(kdm  ,k2+8) )
      write (6,101) (twrk(k1),k1=k2,min(kdm  ,k2+8) )
      write (6,101) (swrk(k1),k1=k2,min(kdm  ,k2+8) )
      write (6,100)
     end do
!SMS$IGNORE END
    end if

    pwrk(:)=pnew(:)
   end do                           ! iter
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- scenario 1A: sequence of 5 thin-thick-thin-thick-thin layers

!  if (mod(nstep,2).eq.0) then		! alternate between scenarios 1A,1B

!   do 1 k=4,kdm-3
!   if (pnew(k-2).lt.pnew(    1)+tencm .or.			&
!       pnew(k+3).gt.pnew(kdm+1)-tencm) go to 1
!   dp1=pnew(k-1)-pnew(k-2)
!   dp2=pnew(k  )-pnew(k-1)
!   dp3=pnew(k+1)-pnew(k  )
!   dp4=pnew(k+2)-pnew(k+1)
!   dp5=pnew(k+3)-pnew(k+2)
!   th2=thet(k-1)
!   th3=thet(k  )
!   th4=thet(k+1)
! --- look for small dp1,dp3,dp5 in combination with large dp2,dp4
!   if (dp2.gt.dp1 .and. dp4.gt.dp5) then
!    goal=.5*(dp3+min(dp2,dp4))		!  desired thknss of lyr 3
!    if (dp2.gt.goal .and. dp4.gt.goal) then
! --- thin-thick-thin-thick-thin combination found -> inflate lyr 3
!     dis3=min(dp2-goal,goal-dp3) * slak
!     dis4=min(dp4-goal,goal-dp3) * slak
!     if (th3.gt.th2 .and. th4.gt.th3) then
! --- targt conservation requires   dis3*(th3-th2)+dis4*(th4-th3)=0
!      if (th4-th3.gt.th3-th2) then
!       dis4=dis3*(th3-th2)/(th4-th3)
!      else
!       dis3=dis4*(th4-th3)/(th3-th2)
!      end if
!     end if
! --- ready to expand middle layer
!     pnew(k  )=pnew(k  )-dis3
!     pnew(k+1)=pnew(k+1)+dis4

!     tdz(k-1)=tdz(k-1)-dis3*twrk(k-1)
!     tdz(k  )=tdz(k  )+dis3*twrk(k-1)+dis4*twrk(k+1)
!     tdz(k+1)=tdz(k+1)               -dis4*twrk(k+1)

!     sdz(k-1)=sdz(k-1)-dis3*swrk(k-1)
!     sdz(k  )=sdz(k  )+dis3*swrk(k-1)+dis4*swrk(k+1)
!     sdz(k+1)=sdz(k+1)               -dis4*swrk(k+1)

!     if (vrbos) then
!!SMS$IGNORE BEGIN
!      write (6,'(a,i3,a,-4p,5f9.3/5x,a,5f9.3)')		&
!       'k=',k,' thknss quintuplet',dp1,dp2,dp3,dp4,dp5,	&
!              '           becomes',dp1,pnew(k)-pnew(k-1),	&
!       pnew(k+1)-pnew(k),pnew(k+2)-pnew(k+1),dp5
!!SMS$IGNORE END
!     end if

!     ncount=ncount+1
!    end if
!   end if
!1  continue
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- scenario 1B: sequence of 5 thick-thin-thick-thin-thick layers

!  else				! mod(nstep,2)=1

!   do 5 k=4,kdm-3
!   if (pnew(k-2).lt.pnew(    1)+tencm .or.			&
!       pnew(k+3).gt.pnew(kdm+1)-tencm) go to 5
!   dp1=pnew(k-1)-pnew(k-2)
!   dp2=pnew(k  )-pnew(k-1)
!   dp3=pnew(k+1)-pnew(k  )
!   dp4=pnew(k+2)-pnew(k+1)
!   dp5=pnew(k+3)-pnew(k+2)
!   th2=thet(k-1)
!   th3=thet(k  )
!   th4=thet(k+1)
! --- look for large dp1,dp3,dp5 in combination with small dp2,dp4
!   if (dp2.lt.dp1 .and. dp4.lt.dp5) then
!    goal=.5*(dp3+max(dp2,dp4))		!  desired thknss of lyr 3
!    if (dp2.lt.goal .and. dp4.lt.goal) then
! --- thick-thin-thick-thin-thick combination found -> deflate lyr 3
!     dis3=min(dp3-goal,goal-dp2) * slak
!     dis4=min(dp3-goal,goal-dp4) * slak
!     if (th3.gt.th2 .and. th4.gt.th3) then
! --- targt conservation requires   dis3*(th3-th2)+dis4*(th4-th3)=0
!      if (th4-th3.gt.th3-th2) then
!       dis4=dis3*(th3-th2)/(th4-th3)
!      else
!       dis3=dis4*(th4-th3)/(th3-th2)
!      end if
!     end if
! --- ready to shrink middle layer
!     pnew(k  )=pnew(k  )+dis3
!     pnew(k+1)=pnew(k+1)-dis4

!     tdz(k-1)=tdz(k-1)+dis3*twrk(k-1)
!     tdz(k  )=tdz(k  )-dis3*twrk(k-1)-dis4*twrk(k+1)
!     tdz(k+1)=tdz(k+1)               +dis4*twrk(k+1)

!     sdz(k-1)=sdz(k-1)+dis3*swrk(k-1)
!     sdz(k  )=sdz(k  )-dis3*swrk(k-1)-dis4*swrk(k+1)
!     sdz(k+1)=sdz(k+1)               +dis4*swrk(k+1)

!     if (vrbos) then
!!SMS$IGNORE BEGIN
!      write (6,'(a,i3,a,-4p,5f9.3/5x,a,5f9.3)')		&
!       'k=',k,' thknss quintuplet',dp1,dp2,dp3,dp4,dp5,	&
!              '           becomes',dp1,pnew(k)-pnew(k-1),	&
!       pnew(k+1)-pnew(k),pnew(k+2)-pnew(k+1),dp5
!!SMS$IGNORE END
!     end if

!     ncount=ncount+1
!    end if
!   end if
!5  continue

!  end if			! nstep odd or even
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- scenario 2: sequence of 3 thin-thick-thin layers
!
!   do 2 k=3,kdm-2
!   if (pnew(k-1).lt.pnew(    1)+tencm .or.			&
!       pnew(k+2).gt.pnew(kdm+1)-tencm) go to 2
!   dp2=pnew(k  )-pnew(k-1)
!   dp3=pnew(k+1)-pnew(k  )
!   dp4=pnew(k+2)-pnew(k+1)
!   th2=thet(k-1)
!   th3=thet(k  )
!   th4=thet(k+1)
! --- look for small dp2,dp4 in combination with large dp3
!   if (dp3.gt.dp2 .and. dp3.gt.dp4) then
!    goal=.5*(dp3+max(dp2,dp4))		!  desired thknss of lyr 3
!    if (dp2.lt.goal .and. dp4.lt.goal) then
! --- thin-thick-thin combination found -> deflate lyr 3
!     dis3=.5*(goal-dp2) * slak
!     dis4=.5*(goal-dp4) * slak
!     if (th3.gt.th2 .and. th4.gt.th3) then
! --- targt conservation requires  dis3*(th3-th2)+dis4*(th4-th3)=0
!      ratio=(th4-th3)/(th3-th2)
!      if (ratio.gt.1.) then
!       dis4=dis3/ratio
!      else
!       dis3=dis4*ratio
!      end if
!     end if
! --- ready to shrink middle layer
!     pnew(k  )=pnew(k  )+dis3
!     pnew(k+1)=pnew(k+1)-dis4
!
!     tdz(k-1)=tdz(k-1)+dis3*twrk(k)
!     tdz(k  )=tdz(k  )-dis3*twrk(k)-dis4*twrk(k)
!     tdz(k+1)=tdz(k+1)             +dis4*twrk(k)
!
!     sdz(k-1)=sdz(k-1)+dis3*swrk(k)
!     sdz(k  )=sdz(k  )-dis3*swrk(k)-dis4*swrk(k)
!     sdz(k+1)=sdz(k+1)             +dis4*swrk(k)
!
!     if (vrbos) then
!!SMS$IGNORE BEGIN
!      write (6,'(a,i3,a,-4p,3f9.3/5x,a,3f9.3)')		&
!       'k=',k,' thknss triplet',dp2,dp3,dp4,			&
!              '        becomes',pnew(k)-pnew(k-1),		&
!       pnew(k+1)-pnew(k),pnew(k+2)-pnew(k+1)
!!SMS$IGNORE END
!     end if
!
!     ncount=ncount+1
!    end if
!   end if
! 2 continue
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- scenario 3: sequence of 3 thick-thin-thick layers
!
!   do 3 k=3,kdm-2
!   if (pnew(k-1).lt.pnew(    1)+tencm .or.			&
!       pnew(k+2).gt.pnew(kdm+1)-tencm) go to 3
!   dp2=pnew(k  )-pnew(k-1)
!   dp3=pnew(k+1)-pnew(k  )
!   dp4=pnew(k+2)-pnew(k+1)
!   th2=thet(k-1)
!   th3=thet(k  )
!   th4=thet(k+1)
! --- look for large dp2,dp4 in combination with small dp3
!   if (dp3.lt.dp2 .and. dp3.lt.dp4) then
! --- thick-thin-thick combination found -> inflate lyr 3
!    goal=.5*(dp3+min(dp2,dp4))		!  desired thknss of lyr 3
!    dis3=.5*(goal-dp3) * slak
!    dis4=dis3
!    if (th3.gt.th2 .and. th4.gt.th3) then
! --- targt conservation requires  dis3*(th3-th2)+dis4*(th4-th3)=0
!     ratio=(th4-th3)/(th3-th2)
!     if (ratio.gt.1.) then
!      dis4=dis3/ratio
!     else
!      dis3=dis4*ratio
!     end if
!    end if
! --- ready to inflate middle layer
!    pnew(k  )=pnew(k  )-dis3
!    pnew(k+1)=pnew(k+1)+dis4
!
!    tdz(k-1)=tdz(k-1)-dis3*twrk(k-1)
!    tdz(k  )=tdz(k  )+dis3*twrk(k-1)+dis4*twrk(k+1)
!    tdz(k+1)=tdz(k+1)               -dis4*twrk(k+1)
!
!    sdz(k-1)=sdz(k-1)-dis3*swrk(k-1)
!    sdz(k  )=sdz(k  )+dis3*swrk(k-1)+dis4*swrk(k+1)
!    sdz(k+1)=sdz(k+1)               -dis4*swrk(k+1)
!
!    if (vrbos) then
!!SMS$IGNORE BEGIN
!     write (6,'(a,i3,a,-4p,3f9.3/5x,a,3f9.3)')			&
!      'k=',k,' thknss triplet',dp2,dp3,dp4,			&
!             '        becomes',pnew(k)-pnew(k-1),		&
!      pnew(k+1)-pnew(k),pnew(k+2)-pnew(k+1)
!!SMS$IGNORE END
!    end if
!
!    ncount=ncount+1
!   end if
! 3 continue
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- scenario 4: sequence of thick-(thin-thin)-thick layers

   do 4 k=4,kdm-2
   if (pnew(k-2).lt.pnew(    1)+tencm .or.			&
       pnew(k+2).gt.pnew(kdm+1)-tencm) go to 4
   dp1=pnew(k-1)-pnew(k-2)
   dp2=pnew(k  )-pnew(k-1)
   dp3=pnew(k+1)-pnew(k  )
   dp4=pnew(k+2)-pnew(k+1)
   th1=thet(k-2)
   th2=thet(k-1)
   th3=thet(k  )
   th4=thet(k+1)
! --- look for small dp2,dp3 in combination with large dp1,dp4
   if (dp1.gt.dp2+dp3 .and. dp4.gt.dp2+dp3) then
! --- thick-thin-thin-thick combination found -> inflate lyrs 2,3
    goal=.5*(dp2+dp3+min(dp1,dp4))	!  desired thknss of lyr 2+3
    dis2=(goal-dp2-dp3) * slak
    dis4=dis2
    if (th2.gt.th1 .and. th3.gt.th2 .and. th4.gt.th3) then
! --- targt conservation requires  dis2*(th2-th1)+dis4*(th3-th4)=0
     ratio=(th4-th3)/(th2-th1)
     if (ratio.gt.1.) then
      dis4=dis2/ratio
     else
      dis2=dis4*ratio
     end if
    end if
! --- ready to expand middle layers
    pnew(k-1)=pnew(k-1)-dis2
    pnew(k+1)=pnew(k+1)+dis4

    tdz(k-2)=tdz(k-2)-dis2*twrk(k-2)
    tdz(k-1)=tdz(k-1)+dis2*twrk(k-2)
    tdz(k  )=tdz(k  )+dis4*twrk(k+1)
    tdz(k+1)=tdz(k+1)-dis4*twrk(k+1)

    sdz(k-2)=sdz(k-2)-dis2*swrk(k-2)
    sdz(k-1)=sdz(k-1)+dis2*swrk(k-2)
    sdz(k  )=sdz(k  )+dis4*swrk(k+1)
    sdz(k+1)=sdz(k+1)-dis4*swrk(k+1)

    if (vrbos) then
!SMS$IGNORE BEGIN
     write (6,'(a,i3,a,-4p,4f9.3/5x,a,4f9.3)')			&
      'k=',k,' thknss quadruplet',dp1,dp2,dp3,dp4,		&
             '           becomes',pnew(k-1)-pnew(k-2),		&
      pnew(k)-pnew(k-1),pnew(k+1)-pnew(k),pnew(k+2)-pnew(k+1)
!SMS$IGNORE END
    end if

    ncount=ncount+1
   end if
 4 continue
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   do k=1,kdm
    if (pnew(k+1)-pnew(k).gt.onemu) then
     temp(k)=tdz(k)/(pnew(k+1)-pnew(k))
     saln(k)=sdz(k)/(pnew(k+1)-pnew(k))
    end if
    pres(k)=pnew(k)
   end do

   event = ncount>20			!  find interesting cases

   do k=1,kdm
    if (pnew(k+1).lt.pnew(k)-onemm) then
     event=.true.
!SMS$IGNORE BEGIN
     write (*,'(i8,i5,a)') i,k,						&
       ' error: nonmonotonic pressure on return from equilb'
!SMS$IGNORE END
    end if
   end do

   if (event .or. vrbos) then
!SMS$IGNORE BEGIN
    write (6,99) i,'  equilb input profile',				&
     '  (4-line groups: p,rho,T,S):'
    do k2=1,kdm,10
     write (6,100) (pwrk(k1),k1=k2,min(kdm+1,k2+10) )
     write (6,101) (thet(k1),k1=k2,min(kdm  ,k2+9 ) )
     write (6,101) (twrk(k1),k1=k2,min(kdm  ,k2+9 ) )
     write (6,101) (swrk(k1),k1=k2,min(kdm  ,k2+9 ) )
     write (6,100)
    end do
    write (6,99) i,							&
     '  equilb output profile','  (5-line groups: p,dp,rho,T,S):',	&
      ncount,' infl.'
    do k2=1,kdm,10
     write (6,100) (pnew(k1),k1=k2,min(kdm+1,k2+10) )
     write (6,102) (int(.01*(pnew(k1)-pres(k1))),k1=k2,min(kdm+1,k2+10) )
     write (6,101) (thet(k1),k1=k2,min(kdm  ,k2+9 ) )
     write (6,101) (temp(k1),k1=k2,min(kdm  ,k2+9 ) )
     write (6,101) (saln(k1),k1=k2,min(kdm  ,k2+9 ) )
     write (6,100)
    end do
 99 format ('i=',i8,2a,i4,a)
100 format (-4p,11f7.1)
102 format (11i7)
101 format (5x,10f7.2)
!SMS$IGNORE END
   end if

   taftr=0.
   saftr=0.
   do k=1,kdm
    taftr=taftr+(pres(k+1)-pres(k))*temp(k)
    saftr=saftr+(pres(k+1)-pres(k))*saln(k)
   end do
   if (abs(taftr-tbfor).gt.acurcy*abs(tbfor)) then
!SMS$IGNORE BEGIN
    print 103,i,' (equilb) bad temp column intgl',			&
     tbfor,taftr-tbfor,(taftr-tbfor)/tbfor
!SMS$IGNORE END
   end if
   if (abs(saftr-sbfor).gt.acurcy*abs(sbfor)) then
!SMS$IGNORE BEGIN
    print 103,i,' (equilb) bad saln column intgl',			&
     tbfor,saftr-sbfor,(saftr-sbfor)/sbfor
103 format (i8,a,es14.6,2es11.3)
!SMS$IGNORE END
   end if

   return
   end subroutine equilb
end module hycom_hybgn1
