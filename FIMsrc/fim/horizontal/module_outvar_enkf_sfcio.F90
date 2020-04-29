MODULE module_outvar_enkf_sfcio

  IMPLICIT NONE

  PRIVATE
  
  PUBLIC :: outvar_enkf_sfcio
  
CONTAINS
  
  SUBROUTINE outvar_enkf_sfcio(time,ph2d)
    
    USE slint, ONLY: slint_init,nn_interp
    USE sfcio_module
    USE module_constants       ,ONLY: lat,lon,grvity,pi
    USE module_control     ,ONLY: nip
    USE module_sfc_variables
    USE specmod, ONLY: gaulats
    USE kinds, ONLY: i_kind,r_double,r_kind,r_single
    USE units,          ONLY: getunit, returnunit
    USE module_outvar_enkf_spectral, ONLY : cycle_hour,&
         jcap,nlonsgauss,nlatsgauss,jcap,idate,fhour,timestr


! External variable declarations:
!sms$distribute (dh,1) begin
    REAL   ,INTENT(IN)    :: ph2d(nip)
!sms$distribute end

    INTEGER, INTENT(IN)    :: time

!local vars

    TYPE(sfcio_head) :: sfchead
    TYPE(sfcio_data) :: sfcdata

    REAL, DIMENSION(nip) :: topo
    REAL, DIMENSION(nip,2)      :: latlonicos
    REAL, DIMENSION(nlonsgauss*nlatsgauss,2) :: latlongauss
    REAL, ALLOCATABLE, DIMENSION(:) :: ug

    CHARACTER(len=250) :: filename,command
    INTEGER :: npts,lonb,latb,imo,jmo,idrt,lsoil
    INTEGER :: i,j,k,indx,big_endian_lun,ierr

    WRITE(6,*)'Writing out surface enkfio_out file at time=',timestr
    CALL FLUSH(6)
    
    filename='fim_out_sfcio_'//TRIM(timestr)

    imo=nlonsgauss
    jmo=nlatsgauss
    idrt=4
    lsoil=4

    npts = nlonsgauss*nlatsgauss

!sms$serial (<lat,lon,in>:default=ignore) begin

    CALL sfcio_alhead(sfchead,ierr,jmo,lsoil)

    sfchead%fhour = fhour
    sfchead%idate = idate
    sfchead%latb = jmo
    sfchead%lonb = imo
    sfchead%ivs = 200509
    sfchead%lsoil = lsoil
    sfchead%irealf = 1
    sfchead%lpl = imo
    sfchead%zsoil=(/-0.1,-0.4,-1.0,-2.0/)

    CALL sfcio_aldata(sfchead,sfcdata,ierr)

    latlonicos(:,1)=lat
    latlonicos(:,2)=lon

    ALLOCATE(ug(npts))
    
    DO i = 1, nlonsgauss
       DO j = 1, nlatsgauss
          indx=(j - 1) * nlonsgauss + i
          latlongauss(indx,1)=ASIN(gaulats(j))
          latlongauss(indx,2)=2.*pi*float(i-1)/nlonsgauss
       END DO
    END DO

    CALL slint_init(latlonicos,nip,latlongauss,npts)

!sms$serial end

!3d soil data

!sms$serial (<st3d,in>:default=ignore) begin

    DO k=1,lsoil

       CALL nn_interp(st3d(k,:),ug)
       DO j=1,nlatsgauss
          indx=(j - 1) * nlonsgauss
          sfcdata%stc(:,j,k)=ug(indx+1:indx+nlonsgauss)
       ENDDO

    ENDDO

!sms$serial end

!sms$serial (<sm3d,in>:default=ignore) begin

    DO k=1,lsoil

       CALL nn_interp(sm3d(k,:),ug)
       DO j=1,nlatsgauss
          indx=(j - 1) * nlonsgauss
          sfcdata%smc(:,j,k)=ug(indx+1:indx+nlonsgauss)
       ENDDO

    ENDDO

!sms$serial end

!sms$serial (<slc3d,in>:default=ignore) begin
 DO k=1,lsoil

       CALL nn_interp(slc3d(k,:),ug)
       DO j=1,nlatsgauss
          indx=(j - 1) * nlonsgauss
          sfcdata%slc(:,j,k)=ug(indx+1:indx+nlonsgauss)
       ENDDO

    ENDDO

!sms$serial end

!2d surface data

!do serials on 2d arrays in chunks of 5 because of memory issues

!sms$serial (<ts2d,sheleg2d,tg32d,zorl2d,cv2d,in>:default=ignore) begin

    CALL nn_interp(ts2d,ug)
    DO j=1,nlatsgauss
       indx=(j - 1) * nlonsgauss
       sfcdata%tsea(:,j)=ug(indx+1:indx+nlonsgauss)
    ENDDO
    
    sfcdata%tisfc=sfcdata%tsea

    CALL nn_interp(sheleg2d,ug)
    DO j=1,nlatsgauss
       indx=(j - 1) * nlonsgauss
       sfcdata%sheleg(:,j)=ug(indx+1:indx+nlonsgauss)
    ENDDO

    CALL nn_interp(tg32d,ug)
    DO j=1,nlatsgauss
       indx=(j - 1) * nlonsgauss
       sfcdata%tg3(:,j)=ug(indx+1:indx+nlonsgauss)
    ENDDO

    CALL nn_interp(zorl2d,ug)
    DO j=1,nlatsgauss
       indx=(j - 1) * nlonsgauss
       sfcdata%zorl(:,j)=ug(indx+1:indx+nlonsgauss)
    ENDDO

    CALL nn_interp(cv2d,ug)
    DO j=1,nlatsgauss
       indx=(j - 1) * nlonsgauss
       sfcdata%cv(:,j)=ug(indx+1:indx+nlonsgauss)
    ENDDO

!sms$serial end

!sms$serial (<cvb2d,cvt2d,alvsf2d,alvwf2d,alnsf2d,in>:default=ignore) begin

    CALL nn_interp(cvb2d,ug)
    DO j=1,nlatsgauss
       indx=(j - 1) * nlonsgauss
       sfcdata%cvb(:,j)=ug(indx+1:indx+nlonsgauss)
    ENDDO

    CALL nn_interp(cvt2d,ug)
    DO j=1,nlatsgauss
       indx=(j - 1) * nlonsgauss
       sfcdata%cvt(:,j)=ug(indx+1:indx+nlonsgauss)
    ENDDO

    CALL nn_interp(alvsf2d,ug)
    DO j=1,nlatsgauss
       indx=(j - 1) * nlonsgauss
       sfcdata%alvsf(:,j)=ug(indx+1:indx+nlonsgauss)
    ENDDO

    CALL nn_interp(alvwf2d,ug)
    DO j=1,nlatsgauss
       indx=(j - 1) * nlonsgauss
       sfcdata%alvwf(:,j)=ug(indx+1:indx+nlonsgauss)
    ENDDO

    CALL nn_interp(alnsf2d,ug)
    DO j=1,nlatsgauss
       indx=(j - 1) * nlonsgauss
       sfcdata%alnsf(:,j)=ug(indx+1:indx+nlonsgauss)
    ENDDO

!sms$serial end

!sms$serial (<alnwf2d,slmsk2d,vfrac2d,canopy2d,f10m2d,in>:default=ignore) begin    

    CALL nn_interp(alnwf2d,ug)
    DO j=1,nlatsgauss
       indx=(j - 1) * nlonsgauss
       sfcdata%alnwf(:,j)=ug(indx+1:indx+nlonsgauss)
    ENDDO

    CALL nn_interp(slmsk2d,ug)
    DO j=1,nlatsgauss
       indx=(j - 1) * nlonsgauss
       sfcdata%slmsk(:,j)=ug(indx+1:indx+nlonsgauss)
    ENDDO

    CALL nn_interp(vfrac2d,ug)
    DO j=1,nlatsgauss
       indx=(j - 1) * nlonsgauss
       sfcdata%vfrac(:,j)=ug(indx+1:indx+nlonsgauss)
    ENDDO

    CALL nn_interp(canopy2d,ug)
    DO j=1,nlatsgauss
       indx=(j - 1) * nlonsgauss
       sfcdata%canopy(:,j)=ug(indx+1:indx+nlonsgauss)
    ENDDO

    CALL nn_interp(f10m2d,ug)
    DO j=1,nlatsgauss
       indx=(j - 1) * nlonsgauss
       sfcdata%f10m(:,j)=ug(indx+1:indx+nlonsgauss)
    ENDDO

!sms$serial end

!sms$serial (<t2m2d,q2m2d,vtype2d,stype2d,facsf2d,in>:default=ignore) begin    

    CALL nn_interp(t2m2d,ug)
    DO j=1,nlatsgauss
       indx=(j - 1) * nlonsgauss
       sfcdata%t2m(:,j)=ug(indx+1:indx+nlonsgauss)
    ENDDO

    CALL nn_interp(q2m2d,ug)
    DO j=1,nlatsgauss
       indx=(j - 1) * nlonsgauss
       sfcdata%q2m(:,j)=ug(indx+1:indx+nlonsgauss)
    ENDDO

    CALL nn_interp(vtype2d,ug)
    DO j=1,nlatsgauss
       indx=(j - 1) * nlonsgauss
       sfcdata%vtype(:,j)=ug(indx+1:indx+nlonsgauss)
    ENDDO

    CALL nn_interp(stype2d,ug)
    DO j=1,nlatsgauss
       indx=(j - 1) * nlonsgauss
       sfcdata%stype(:,j)=ug(indx+1:indx+nlonsgauss)
    ENDDO

    CALL nn_interp(facsf2d,ug)
    DO j=1,nlatsgauss
       indx=(j - 1) * nlonsgauss
       sfcdata%facsf(:,j)=ug(indx+1:indx+nlonsgauss)
    ENDDO

!sms$serial end

!sms$serial (<facwf2d,uustar2d,ffmm2d,ffhh2d,hice2d,in>:default=ignore) begin

    CALL nn_interp(facwf2d,ug)
    DO j=1,nlatsgauss
       indx=(j - 1) * nlonsgauss
       sfcdata%facwf(:,j)=ug(indx+1:indx+nlonsgauss)
    ENDDO
 
    CALL nn_interp(uustar2d,ug)
    DO j=1,nlatsgauss
       indx=(j - 1) * nlonsgauss
       sfcdata%uustar(:,j)=ug(indx+1:indx+nlonsgauss)
    ENDDO

    CALL nn_interp(ffmm2d,ug)
    DO j=1,nlatsgauss
       indx=(j - 1) * nlonsgauss
       sfcdata%ffmm(:,j)=ug(indx+1:indx+nlonsgauss)
    ENDDO
 
    CALL nn_interp(ffhh2d,ug)
    DO j=1,nlatsgauss
       indx=(j - 1) * nlonsgauss
       sfcdata%ffhh(:,j)=ug(indx+1:indx+nlonsgauss)
    ENDDO

    CALL nn_interp(hice2d,ug)
    DO j=1,nlatsgauss
       indx=(j - 1) * nlonsgauss
       sfcdata%hice(:,j)=ug(indx+1:indx+nlonsgauss)
    ENDDO
 
!sms$serial end

!sms$serial (<fice2d,tprcp2d,srflag2d,snwdph2d,shdmin2d,shdmax2d,in>:default=ignore) begin

    CALL nn_interp(fice2d,ug)
    DO j=1,nlatsgauss
       indx=(j - 1) * nlonsgauss
       sfcdata%fice(:,j)=ug(indx+1:indx+nlonsgauss)
    ENDDO
 
    CALL nn_interp(tprcp2d,ug)
    DO j=1,nlatsgauss
       indx=(j - 1) * nlonsgauss
       sfcdata%tprcp(:,j)=ug(indx+1:indx+nlonsgauss)
    ENDDO
    
    CALL nn_interp(srflag2d,ug)
    DO j=1,nlatsgauss
       indx=(j - 1) * nlonsgauss
       sfcdata%srflag(:,j)=ug(indx+1:indx+nlonsgauss)
    ENDDO
 
    CALL nn_interp(snwdph2d,ug)
    DO j=1,nlatsgauss
       indx=(j - 1) * nlonsgauss 
       sfcdata%snwdph(:,j)=ug(indx+1:indx+nlonsgauss)
    ENDDO

    CALL nn_interp(shdmin2d,ug)
    DO j=1,nlatsgauss
       indx=(j - 1) * nlonsgauss 
       sfcdata%shdmin(:,j)=ug(indx+1:indx+nlonsgauss)
    ENDDO

    CALL nn_interp(shdmax2d,ug)
    DO j=1,nlatsgauss
       indx=(j - 1) * nlonsgauss 
       sfcdata%shdmax(:,j)=ug(indx+1:indx+nlonsgauss)
    ENDDO

!sms$serial end

    big_endian_lun=getunit(82)
    if (big_endian_lun.lt.0) then
      big_endian_lun=getunit(76)
      if (big_endian_lun.lt.0) then
        write (*,'(a)') &
          'outvar_enkf_sfcio: Could not obtain big-endian unit number.'
        stop
      endif
    endif

!sms$serial (<slope2d,snoalb2d,ph2d,in>:default=ignore) begin

    CALL nn_interp(slope2d,ug)
    DO j=1,nlatsgauss
       indx=(j - 1) * nlonsgauss 
       sfcdata%slope(:,j)=ug(indx+1:indx+nlonsgauss)
    ENDDO
    
    CALL nn_interp(snoalb2d,ug)
    DO j=1,nlatsgauss
       indx=(j - 1) * nlonsgauss
       sfcdata%snoalb(:,j)=ug(indx+1:indx+nlonsgauss)
    ENDDO

    topo=ph2d/grvity

    CALL nn_interp(topo,ug)
    DO j=1,nlatsgauss
       indx=(j - 1) * nlonsgauss 
       sfcdata%orog(:,j)=ug(indx+1:indx+nlonsgauss)
    ENDDO

    DEALLOCATE(ug)

    CALL sfcio_swohdc(big_endian_lun,filename,sfchead,sfcdata,ierr)
    CALL sfcio_axdata(sfcdata,ierr)

!sms$serial end

    CALL returnunit(big_endian_lun)

  END SUBROUTINE outvar_enkf_sfcio

END MODULE module_outvar_enkf_sfcio
