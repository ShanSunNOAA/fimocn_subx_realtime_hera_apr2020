module hycom_sverbal
use findmaxmin1
use stencilprint
use stenedgprint
contains
!*********************************************************************
!  ocn_flxsum, sverbal, getcurl
!
!  a collection of routines for comparing boundary current transport
!  to the transport implied by Sverdrup balance
!
!  R.Bleck				Dec. 2011
!*********************************************************************
   subroutine sverbal(nstep,curl_avg)

! --- compute sverdrup transport (zonally integrated wind stress curl)

   use module_control,  only: nip,npp
   use fimnamelist,     only: ArchvTimeUnit,itest,diag_intvl
   use module_constants,only: prox,nprox,area,deg_lat,deg_lon,		&
                              nedge,permedge,rarea,corio,omegx2,raddeg, pi
   use hycom_control   ,only: bdcurrsecs,nsecs
   use hycom_constants ,only: wet,thref
   use hycom_variables ,only: lgthset,cellset,ctrlat,ctrlon

   implicit none
   integer,intent(IN)  :: nstep
!SMS$DISTRIBUTE (dh,1) BEGIN
   real   ,intent(IN)  :: curl_avg(nip)
!SMS$DISTRIBUTE END
   integer n,nm1,np1,ns,lgthsub,i,len,maxlen
   integer,allocatable :: cellsub(:)
   real   ,allocatable :: latsub(:),lonsub(:)
   real                :: trnsp,beta,q,dlon,reflat,numbr
   character           :: text*16
   logical,parameter   :: vrbos=.false.
   real,   external    :: its2time

#if ( defined NEED_SINDCOSD )
!JR Define required statement functions for situations where
!JR compiler doesn't support them (e.g. gfortran)
   real :: val, cosd
   cosd(val) = cos(val*pi/180.)
#endif

   print *,'entering sverbal...'
   write (text,'(i9,a)') nstep,' curlav'
   call findmxmn1(curl_avg,nip,text,wet)

! --- extract list of icos cells from class 2 transects

   maxlen=maxval(lgthset)
   allocate (cellsub(maxlen),latsub(maxlen),lonsub(maxlen))
   do n=1,nsecs
    if (bdcurrsecs(n)(1:1).ne.' ') then
     cellsub(:)=cellset(:,n,2)
     latsub(:)=ctrlat(:,n,2)
     lonsub(:)=ctrlon(:,n,2)
     lgthsub=lgthset(n,2)
     nm1=1
     do ns=2,lgthsub
      if (cellsub(ns).ne.cellsub(nm1)) then
       nm1=nm1+1
       cellsub(nm1)=cellsub(ns)
       latsub(nm1)=latsub(ns)
       lonsub(nm1)=lonsub(ns)
      end if
     end do
     lgthsub=nm1

     if (vrbos) then
      print '(2a)','icos cells only (no edges) for transect ',		&
       trim(bdcurrsecs(n))
      print 100,(cellsub(len),latsub(len),lonsub(len),len=1,lgthsub)
100  format (3(i10,2f8.3))
     end if

! --- evaluate zonal integral of time-averaged wind stress curl

     trnsp=0.
     reflat=0.
     numbr=0.
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE(np1,nm1,dlon,q) REDUCTION(+:trnsp,reflat,numbr)
     do i=1,nip
      if (wet(i) > 0 ) then
       do len=1,lgthsub
        if (cellsub(len).eq.i) then	! add contribution of cell -i- to intgl
         numbr=numbr+1.
         np1=min(lgthsub,n+1)
         nm1=max(   1,n-1)
         dlon=(mod(lonsub(np1)-lonsub(nm1)+540.,360.)-180.)/(np1-nm1)
         q=curl_avg(i)*dlon
         trnsp=trnsp+q
         reflat=reflat+latsub(n)
        end if
       end do				! loop over cell set
      end if				! ocean point
     end do				! loop over longitude
!$OMP END PARALLEL DO
!SMS$PARALLEL END
!SMS$REDUCE (trnsp,reflat,SUM)
!SMS$REDUCE (numbr,SUM)
     reflat=reflat/numbr
     beta=omegx2*cosd(reflat)/(raddeg*111.2e3)
     print '(a,f6.1,a,es11.2)','beta at',reflat,' deg:',beta
     trnsp=trnsp*thref*111.2e3*cosd(reflat)/beta
     print '(a,f8.1,2(2x,a),f9.2)',ArchvTimeUnit,its2time(nstep),	&
       trim(bdcurrsecs(n)),'Sverdrup trnsp (Sv):',trnsp*1.e-6
    end if
   end do				! loop over transects
   deallocate (cellsub,latsub,lonsub)

   print *,'...exiting sverbal'
   return
   end subroutine sverbal


   subroutine getcurl(nstep,taux,tauy,curl,curl_avg)

! --- accumulate time integral of wind stress curl

   use module_control,  only: nip,npp
   use fimnamelist,     only: ArchvTimeUnit,itest,diag_intvl
   use module_constants,only: prox,nprox,sidevec_e,cs,sn,rarea,		&
                              deg_lat,deg_lon,nedge,permedge,perm
   use hycom_constants ,only: wet

   implicit none

!SMS$DISTRIBUTE (dh,1) BEGIN
   real,intent(IN)    :: taux(nip),tauy(nip)
   real,intent(OUT)   :: curl(nip)
   real,intent(INOUT) :: curl_avg(nip)
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE (dh,2) BEGIN
   real               :: taux_edg(npp,nip),tauy_edg(npp,nip)
!SMS$DISTRIBUTE END
   integer,intent(IN) :: nstep
   integer i,ix,im1,ip1,edg,edgcount
   logical vrbos
   real,parameter :: divby18=1./18.
!  The following are u and v at neighboring icos points, NOT on edge
   real :: u_xy1,u_xy2,u_xy3,u_xy4 ! taux on the xy local grid (m/s), at prox pt
   real :: v_xy1,v_xy2,v_xy3,v_xy4 ! tauy on the xy local grid (m/s), at prox pt
   real,   external    :: its2time
   character text*16

! --- in preparation for wind stress curl calculation, find tangential
! --- stresses around each cell

!SMS$PARALLEL (dh,i) BEGIN
!SMS$EXCHANGE(taux,tauy)
!$OMP PARALLEL DO PRIVATE (vrbos,edg,ix,im1,ip1,u_xy1,v_xy1,u_xy2,v_xy2,u_xy3,v_xy3,u_xy4,v_xy4)
   do i=1,nip				! horizontal loop
    if (wet(i) > 0) then
     vrbos=i.eq.itest .and. mod(nstep,diag_intvl).eq.0
     taux_edg(:,i)=0.
     tauy_edg(:,i)=0.
     do edg=1,nprox(i)
    
! --- edge quantities are interpolated from 4 icos cells -- the cells on
! --- either side of the edge plus the 2 joint neighbors of this pair.
    
      ix=prox(edg,i)
      im1=mod(edg-2+nprox(i),nprox(i))+1
      ip1=mod(edg           ,nprox(i))+1
      im1=prox(im1,i)
      ip1=prox(ip1,i)
    
!  Transform taux,tauy at neighboring icos pt to local coord.system.
!  cs and sn are coordinate transformation constants.
!  u_xy,v_xy are values of taux and tauy rotated into local system.
    
      u_xy1= cs(1,edg,i)*taux(i)+sn(1,edg,i)*tauy(i)
      v_xy1=-sn(1,edg,i)*taux(i)+cs(1,edg,i)*tauy(i)
      if (wet(ix) > 0 ) then
       u_xy2= cs(2,edg,i)*taux(ix)+sn(2,edg,i)*tauy(ix)
       v_xy2=-sn(2,edg,i)*taux(ix)+cs(2,edg,i)*tauy(ix)
       if (wet(i) == 0 ) then		! ix is ocean point but i is not
        u_xy1=u_xy2
        v_xy1=v_xy2
       end if
      else				! ix is land point
       u_xy2=u_xy1
       v_xy2=v_xy1
      end if
      if (wet(im1) > 0 ) then
       u_xy3= cs(3,edg,i)*taux(im1)+sn(3,edg,i)*tauy(im1)
       v_xy3=-sn(3,edg,i)*taux(im1)+cs(3,edg,i)*tauy(im1)
      else				! im1 is land point
       u_xy3=.5*(u_xy1+u_xy2)
       v_xy3=.5*(v_xy1+v_xy2)
      end if
      if (wet(ip1) > 0 ) then
       u_xy4= cs(4,edg,i)*taux(ip1)+sn(4,edg,i)*tauy(ip1)
       v_xy4=-sn(4,edg,i)*taux(ip1)+cs(4,edg,i)*tauy(ip1)
      else				! ip1 is land point
       u_xy4=.5*(u_xy1+u_xy2)
       v_xy4=.5*(v_xy1+v_xy2)
      end if
  
      if (vrbos) then
       print 100,perm(i),edg,perm(i),wet(i),perm(ix),wet(ix),		&
         perm(im1),wet(im1),perm(ip1),wet(ip1)
 100  format (i8,' edg=',i1,' intpol based on icos cells',		&
        2(i9,i2),' (lrg wgt),',/40x,2(i9,i2),' (sml wgt)')
       print '(a,8es9.2)','details:',					&
         cs(1,edg,i),taux(i ),sn(1,edg,i),tauy(i ),			&
         cs(2,edg,i),taux(ix),sn(2,edg,i),tauy(ix)
      end if
    
! --- interpolate rotated stress components to edges
  
      taux_edg(edg,i)=(8.*(u_xy1+u_xy2)+u_xy3+u_xy4)*divby18
      tauy_edg(edg,i)=(8.*(v_xy1+v_xy2)+v_xy3+v_xy4)*divby18
      if (vrbos) then
       print '(i8,a,i1,a,5es11.2/14x,a,5es11.2)',perm(i),' edg=',edg,	&
        '  taux_edg:',taux_edg(edg,i),u_xy1,u_xy2,u_xy3,u_xy4,		&
        '  tauy_edg:',tauy_edg(edg,i),v_xy1,v_xy2,v_xy3,v_xy4
      end if
     end do				! loop over edges
    end if				! ocean point
   end do				! horizontal loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END

! --- knowing tangential stresses, we can compute wind stress curl

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (edg)
   do i=1,nip
    curl(i)=0.
    if (wet(i) > 0 ) then
     do edg=1,nprox(i)			! loop through edges
      curl(i)=curl(i)+(sidevec_e(1,edg,i)*taux_edg(edg,i)		&
                      +sidevec_e(2,edg,i)*tauy_edg(edg,i))*rarea(i)
     end do
     curl_avg(i)=curl_avg(i)+curl(i)
    end if				! ocean point
   end do				! horizontal loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END

!  write (text,'(a,i9)') 'step',nstep
!  call stencl(taux,1,1.e4,trim(text)//'  taux x 1.e4')
!  call stencl(tauy,1,1.e4,trim(text)//'  tauy x 1.e4')
!  call stencl(curl,1,1.e9,trim(text)//'  wind stress curl x 10^9')
   if (mod(nstep,diag_intvl).eq.0) then
    write (text,'(a,f8.1,a)') ArchvTimeUnit,its2time(nstep),' curl'
    call findmxmn1(curl,nip,text,wet)
   end if
   return
   end subroutine getcurl


   subroutine ocn_flxsum(nstep,massflx)

! --- sum up mass fluxes normal to a specified line of icos edges

   use module_control,  only: nip,npp
   use fimnamelist,     only: ArchvTimeUnit,kdm,itest,diag_intvl
   use module_constants,only: prox,nprox,perm
   use hycom_control   ,only: bdcurrsecs,bylayrsecs,nsecs,thruflsecs
   use hycom_constants ,only: wet,r_onem,odepth
   use hycom_variables ,only: lgthset,cellset,edgeset,sense,	&
                              ctrlat,ctrlon,transp

   implicit none
   integer,intent(IN) :: nstep
!SMS$DISTRIBUTE (dh,3) BEGIN
   real   ,intent(IN) :: massflx(kdm,npp,nip)
!SMS$DISTRIBUTE END
   integer            :: n,i,ix,len,k
   real               :: sum
   logical,parameter  :: vrbos=.false.
   real,  external    :: its2time

!  print *,'entering ocn_flxsum...'
   print *,'-------ocn.flxsum--------ocn.flxsum--------ocn.flxsum--------ocn.flxsum-------'
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- flow through passages (class 1 transects):
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   do n=1,nsecs
    if (thruflsecs(n)(1:1).ne.' ') then
     print '(2a)','starting thruflow calculation for ',			&
       trim(thruflsecs(n))
     if (vrbos) then
      print '(2a)',trim(thruflsecs(n)),' -- cell/edge info:'
      print 100,(cellset(len,n,1),edgeset(len,n,1),sense(len,n,1),	&
       ctrlat(len,n,1),ctrlon(len,n,1),len=1,lgthset(n,1))
     end if

     transp(:,:,n,1)=0.
     do k=1,kdm
      do len=1,lgthset(n,1)
       sum=0.
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (ix) REDUCTION(+:sum)
       do i=1,nip
        if (cellset(len,n,1).eq.i) then
         if (wet(i) > 0 ) then
          ix=prox(edgeset(len,n,1),i)
          if (wet(ix) > 0 ) then
! --- are we circling the cell in clock- or counterclockwise direction?
           if (sense(len,n,1).eq.'clockw') then
            sum= massflx(k,edgeset(len,n,1),cellset(len,n,1))		&
                *r_onem*1.e-6					! -> Sv
           else if (sense(len,n,1).eq.'countr') then
            sum=-massflx(k,edgeset(len,n,1),cellset(len,n,1))		&
                *r_onem*1.e-6					! -> Sv
           else
            stop '(wrong characters in sense array)'
           end if
           if (vrbos) then
            if (sum.ge.0.) then
             print  101,k,i,'  exports', sum,' Sv across edge',		&
              edgeset(len,n,1),'    to',perm(ix)
            else
             print  101,k,i,'  imports',-sum,' Sv across edge',		&
              edgeset(len,n,1),'  from',perm(ix)
            end if
 101        format('(flxsum)  k=',i2,'  cell',i8,a,f9.2,a,i2,a,i8)
           end if
          end if		! neighbor is ocean cell
         end if			! i is ocean cell
        end if			! cell i belongs to transect
       end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END
!SMS$REDUCE (sum,SUM)
       if (len.eq.1) then
        transp(k,len,n,1)=sum
       else
        transp(k,len,n,1)=transp(k,len-1,n,1)+sum	! sum up transport
       end if
      end do			! loop over cells in transect

!     if (vrbos) print 103,'(flxsum)  k=',k,				&
!       'transport contribution (cumul.) from indiv.cell edges:'	&
!       ,(cellset(len,n,1),edgeset(len,n,1),transp(k,len,n,1),		&
!       len=1,lgthset(n,1))
103   format (a,i2,2x,a/(3(i11,i3,f9.2)))
      print 102,ArchvTimeUnit,its2time(nstep),'  k=',k,			&
        trim(thruflsecs(n)),'thruflow (Sv):',transp(k,lgthset(n,1),n,1)
102   format (a,f8.1,a,i2,2(2x,a),f9.2)
     end do			! vertical loop

     do k=1,kdm
      if (k.gt.1) then
       do len=1,lgthset(n,1)
        transp(k,len,n,1)=transp(k,len,n,1)+transp(k-1,len,n,1)
       end do
      end if
!     print 102,ArchvTimeUnit,its2time(nstep),'  k=1-',k,		&
!       trim(thruflsecs(n)),'thruflow (Sv):',transp(k,lgthset(n,1),n,1)
     end do			! vertical loop
     print 107,ArchvTimeUnit,its2time(nstep),trim(thruflsecs(n)),	&
       ' thruflow (cumul.)    ',(transp(k,lgthset(n,1),n,1),k=1,kdm)
     print 107,ArchvTimeUnit,its2time(nstep),trim(thruflsecs(n)),	&
       ' thruflow min,max,totl',minval(transp(:,lgthset(n,1),n,1)),	&
                                maxval(transp(:,lgthset(n,1),n,1)),	&
                                transp(kdm,lgthset(n,1),n,1)
107  format (a2,f8.1,1x,a18,a22,4f7.2/(2x,11f7.2))
    end if
   end do			! loop over transects

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- western boundary current transport (class 2 transects):
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   do n=1,nsecs
    if (bdcurrsecs(n)(1:1).ne.' ') then
     print '(2a)','starting boundary current calculation for ',		&
       trim(bdcurrsecs(n))
     if (vrbos) then
      print '(2a)',trim(bdcurrsecs(n)),' -- cell/edge info:'
      print 100,(cellset(len,n,2),edgeset(len,n,2),sense(len,n,2),	&
       ctrlat(len,n,2),ctrlon(len,n,2),len=1,lgthset(n,2))
100 format (2(i10,i2,1x,a6,2f8.2))
     end if

     transp(1,:,n,2)=0.
     do len=1,lgthset(n,2)
      sum=0.
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (ix) REDUCTION(+:sum)
      do i=1,nip
       if (cellset(len,n,2).eq.i) then
        if (wet(i) > 0 ) then
         ix=prox(edgeset(len,n,2),i)
         if (wet(ix) > 0 ) then
          if (sense(len,n,2).eq.'clockw') then
           do k=1,kdm
            sum=sum+massflx(k,edgeset(len,n,2),cellset(len,n,2))	&
                *r_onem*1.e-6					! -> Sv
           end do
          else if (sense(len,n,2).eq.'countr') then
           do k=1,kdm
            sum=sum-massflx(k,edgeset(len,n,2),cellset(len,n,2))	&
                *r_onem*1.e-6					! -> Sv
           end do
          else
           stop '(wrong characters in sense array)'
          end if
          if (vrbos) then
           if (sum.ge.0.) then
            print  106,i,'  exports', sum,' Sv across edge',		&
             edgeset(len,n,2),'    to',perm(ix)
           else
            print  106,i,'  imports',-sum,' Sv across edge',		&
             edgeset(len,n,2),'  from',perm(ix)
          end if
106        format('(flxsum)  cell',i8,a,f9.2,a,i2,a,i8)
          end if
         end if			! neighbor is ocean cell
        end if			! i is ocean cell
       end if			! cell i belongs to transect
      end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END
!SMS$REDUCE (sum,SUM)
      if (len.eq.1) then
       transp(1,len,n,2)=sum
      else
       transp(1,len,n,2)=transp(1,len-1,n,2)+sum
      end if
     end do			! loop over cells in transect
     print 105,								&
      '(flxsum) transport contribution (cumul.) from indiv.cell edges:'	&
       ,(cellset(len,n,2),edgeset(len,n,2),transp(1,len,n,2),		&
        len=1,lgthset(n,2))
105  format (a/(3(i11,i3,f9.2)))
     print 104,ArchvTimeUnit,its2time(nstep),trim(bdcurrsecs(n)),	&
       'bdry current transport (Sv):',maxval(transp(1,:,n,2))
104   format (a,f8.1,2(2x,a),f9.2)
    end if
   end do			! loop over transects

   print *,'-------ocn.flxsum--------ocn.flxsum--------ocn.flxsum--------ocn.flxsum-------'
!  print *,'...exiting ocn_flxsum'
   return
   end subroutine ocn_flxsum
end module hycom_sverbal
