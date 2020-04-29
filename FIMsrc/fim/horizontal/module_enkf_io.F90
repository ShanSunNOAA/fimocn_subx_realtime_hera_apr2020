!to write out and read in state vars for enkf + soil
!modified from restart.F90

MODULE module_enkf_io
  
  USE module_control,   ONLY: nip,nvlsig,ntra,ntrb
  USE fimnamelist,      ONLY: nvl,curve,pure_sig,bckg_fname
  USE module_constants, ONLY: thetac, sigak, sigbk
  USE module_variables, ONLY: us3d,vs3d,pr3d,ph3d,tr3d,wt_int_rev,k_int_rev

  USE module_sfc_variables, ONLY : st3d,sm3d,slc3d

  USE units,          ONLY: getunit, returnunit

  IMPLICIT NONE
  
#include <gptl.inc>

  PRIVATE
  PUBLIC :: write_enkfio, read_enkfio,unitno_enkf

  INTEGER :: unitno_enkf
  
CONTAINS
  
! write_enkfio: open the output enkf file and write state vars + soil
  
  SUBROUTINE write_enkfio ()

    INTEGER           :: ret         ! function return

    INTEGER :: i,nsoil,status

!sms$ignore begin
    nsoil=SIZE(st3d,1)
!sms$ignore end

    ret = gptlstart ('write_enkfio')

    unitno_enkf = getunit ()
    IF (unitno_enkf < 0) THEN
       WRITE(6,*)'write_enkfio: Bad return from getunit'
       STOP
    END IF

    OPEN (unitno_enkf,file=TRIM(bckg_fname),form='unformatted',action='write',status='new',iostat=status)
    if (status.ne.0) then
      print *,'write_enkfio: Error opening unit ',unitno_enkf,'. Stopping...'
      stop
    endif

    ret = gptlprint_memusage ('before write_enkfio')

    WRITE (unitno_enkf,iostat=status) curve,nip,nvl,nvlsig,ntra,ntrb,nsoil,pure_sig
    if (status.ne.0) then
      print *,'write_enkfio: Write error 1. Stopping...'
      stop
    endif

    WRITE (unitno_enkf,iostat=status) sigak
    if (status.ne.0) then
      print *,'write_enkfio: Write error 2. Stopping...'
      stop
    endif

    WRITE (unitno_enkf,iostat=status) sigbk
    if (status.ne.0) then
      print *,'write_enkfio: Write error 3. Stopping...'
      stop
    endif

    CALL writearr32_r (us3d, nvl, unitno_enkf)
    CALL writearr32_r (vs3d, nvl, unitno_enkf)

    DO i=1,ntra+ntrb
       CALL writearr32_r (tr3d(:,:,i), nvl, unitno_enkf)
    ENDDO

    CALL writearr32_r (ph3d(1,:),1, unitno_enkf)

    IF (pure_sig ) THEN
       CALL writearr32_r (pr3d(1,:),1, unitno_enkf)
    ELSE
       CALL writearr32_r (pr3d,nvl+1, unitno_enkf)
       CALL writearr32_r (wt_int_rev, nvl, unitno_enkf)
       CALL writearr32_i (k_int_rev, nvl, unitno_enkf)
    ENDIF

    CALL writearr32_r (st3d,nsoil,unitno_enkf)
    CALL writearr32_r (sm3d,nsoil,unitno_enkf)
    CALL writearr32_r (slc3d,nsoil,unitno_enkf)
    
    ret = gptlprint_memusage ('after write_enkfio')

    CLOSE (unitno_enkf)

    CALL returnunit (unitno_enkf)

    ret = gptlstop ('write_enkfio')

    RETURN
    
  END SUBROUTINE write_enkfio
  
! read_enkfio: open the input enkf_io file and read state vars and soil
  
  SUBROUTINE read_enkfio ()
    INTEGER :: ret         ! function return
    
    INTEGER :: i,curvetmp,niptmp,nvltmp,nvlsigtmp,&
      ntratmp,ntrbtmp,nsoiltmp,status

    LOGICAL :: pure_sigtmp

    ret = gptlstart ('read_enkfio')

    ALLOCATE(sigak(nvlsig+1),sigbk(nvlsig+1))

    unitno_enkf = getunit ()
    
    IF (unitno_enkf < 0) THEN
       WRITE(6,*)'read_enkfio: Bad return from getunit'
       STOP
    END IF

    WRITE(6,*)'read_enkfio: trying to open file: ', &
         TRIM(bckg_fname),' on unit ', unitno_enkf

    OPEN (unitno_enkf,file=TRIM(bckg_fname),form='unformatted',action='read',status='old',iostat=status)
    if (status.ne.0) then
      print *,'read_enkfio: Error opening unit ',unitno_enkf,'. Stopping...'
      stop
    endif

    WRITE (6,*)'read_enkfio: opened read_enkfio file: ', bckg_fname

    ret = gptlprint_memusage ('before read_enkfio')
    
    READ (unitno_enkf,iostat=status) curvetmp,niptmp,nvltmp,nvlsigtmp,ntratmp,ntrbtmp,nsoiltmp,pure_sigtmp
    if (status.ne.0) then
      print *,'read_enkfio: Read error 1. Stopping...'
      stop
    endif

    IF (curvetmp /= curve .OR. niptmp /= nip .OR. nvltmp /= nvl .OR.&
         nvlsigtmp /= nvlsig .OR. &
         ntratmp /= ntra .OR. ntrbtmp /= ntrb .OR. &
         pure_sigtmp .neqv. pure_sig) THEN

       WRITE(6,*)'Wrong parameters on input in read_enkfio - Stopping'
       STOP
    ENDIF

    READ (unitno_enkf,iostat=status) sigak
    if (status.ne.0) then
      print *,'read_enkfio: Read error 2. Stopping...'
      stop
    endif

    READ (unitno_enkf,iostat=status) sigbk
    if (status.ne.0) then
      print *,'read_enkfio: Read error 3. Stopping...'
      stop
    endif

    CALL readarr32_r (us3d, nvl, unitno_enkf)
    CALL readarr32_r (vs3d, nvl, unitno_enkf)

    DO i=1,ntra+ntrb
       CALL readarr32_r (tr3d(:,:,i), nvl, unitno_enkf)
    ENDDO

    CALL readarr32_r (ph3d(1,:),1, unitno_enkf)

    IF ( pure_sig ) THEN
       CALL readarr32_r (pr3d(1,:),1, unitno_enkf)
    ELSE
       CALL readarr32_r (pr3d,nvl+1, unitno_enkf)
       CALL readarr32_r (wt_int_rev, nvl, unitno_enkf)
       CALL readarr32_i (k_int_rev, nvl, unitno_enkf)
    ENDIF

    ret = gptlprint_memusage ('after reading dyn vars in read_enkfio')
    
    ret = gptlstop ('read_enkfio')

    RETURN

  END SUBROUTINE read_enkfio
  
END MODULE module_enkf_io
