!WRF:MEDIATION_LAYER:PHYSICS
!
MODULE module_wrfradiation_driver
   USE MODULE_CONSTANTS,only: g=>grvity,cp,r_d=>rd,r_v=>rv
   REAL    , PARAMETER :: DPD          = 360./365.
   REAL    , PARAMETER :: DEGRAD       = 3.1415926/180.
CONTAINS
!BOP
! !IROUTINE: radiation_driver - interface to radiation physics options

! !INTERFACE:
   SUBROUTINE radiation_driver (                                          &
               itimestep,dt ,lw_physics,sw_physics                        &
              ,RTHRATENLW ,RTHRATENSW ,RTHRATEN                           &
              ,ACSWUPT,ACSWUPTC,ACSWDNT,ACSWDNTC                          & ! Optional
              ,ACSWUPB,ACSWUPBC,ACSWDNB,ACSWDNBC                          & ! Optional
              ,ACLWUPT,ACLWUPTC,ACLWDNT,ACLWDNTC                          & ! Optional
              ,ACLWUPB,ACLWUPBC,ACLWDNB,ACLWDNBC                          & ! Optional
              ,  SWUPT,  SWUPTC,  SWDNT,  SWDNTC                          & ! Optional
              ,  SWUPB,  SWUPBC,  SWDNB,  SWDNBC                          & ! Optional
              ,  LWUPT,  LWUPTC,  LWDNT,  LWDNTC                          & ! Optional
              ,  LWUPB,  LWUPBC,  LWDNB,  LWDNBC                          & ! Optional
              ,LWCF,SWCF,OLR                                              & ! Optional
              ,SWUPFLX,SWUPFLXC,SWDNFLX,SWDNFLXC                          & ! Optional
              ,LWUPFLX,LWUPFLXC,LWDNFLX,LWDNFLXC                          & ! Optional
!jdfcz        ,progn,prescribe                                            & ! Optional
              ,progn                                                      & ! Optional
              ,GLW, GSW, alb_uvb,alb_uvd,alb_irb,alb_ird,XLAT, XLONG      &
              ,EMISS, rho, p8w, p , pi , dz8w ,t, t8w, GMT                &
              ,XLAND, XICE, TSK, HTOP,HBOT,VEGFRA, SNOW     &
              ,JULDAY, xtime, RADT, ICLOUD, warm_rain     &
              ,ids, ide, jds, jde, kds, kde                               &
              ,ims, ime, jms, jme, kms, kme                               &
              ,its, ite, jts, jte, kts, kte                               &
              ,CURR_SECS                       &
              ,qv,qc,qr,qi,qs,qg,qndrop                                   &
              ,f_qv,f_qc,f_qr,f_qi,f_qs,f_qg,f_qndrop                     &
              ,CLDFRA                                                  &
               !aerosol shortwave
              ,ext_cof,sscal,asymp                                     &
               !aerosol longwave
              ,tauaerlw,qc_adjust ,qi_adjust                              & ! jm
              ,cu_rad_feedback, aer_ra_feedback                           & ! jm
                                                                          )


!-------------------------------------------------------------------------

! !USES:
!  USE module_state_description, ONLY : RRTMSCHEME, GFDLLWSCHEME        &
!                                      ,RRTMG_LWSCHEME, RRTMG_SWSCHEME  &
!                                      ,SWRADSCHEME, GSFCSWSCHEME       &
!                                      ,GFDLSWSCHEME, CAMLWSCHEME, CAMSWSCHEME &
!                                      ,HELDSUAREZ                      &
!                                      ,goddardlwscheme                 & 
!                                      ,goddardswscheme                   


!  USE module_model_constants
!  USE module_wrf_error , ONLY : wrf_err_message

! *** add new modules of schemes here

!  USE module_ra_sw         , ONLY : swrad
!  USE module_ra_gsfcsw     , ONLY : gsfcswrad
!  USE module_ra_rrtm       , ONLY : rrtmlwrad
   USE module_ra_rrtmg_lw   , ONLY : rrtmg_lwrad
   USE module_ra_rrtmg_sw   , ONLY : rrtmg_swrad
!  USE module_ra_cam        , ONLY : camrad
!  USE module_ra_gfdleta    , ONLY : etara
!#ifdef HWRF
!   USE module_ra_hwrf
!#endif
!   USE module_ra_hs         , ONLY : hsrad

!   USE module_ra_goddard    , ONLY : goddardrad  !Goddard Radiation Scheme by Toshihisa Matsui and Jainn Shi



   !  This driver calls subroutines for the radiation parameterizations.
   !
   !  short wave radiation choices:
   !  1. swrad (19??)
   !  4. rrtmg_sw - Added November 2008, MJIacono, AER, Inc.
   !
   !  long wave radiation choices:
   !  1. rrtmlwrad
   !  4. rrtmg_lw - Added November 2008, MJIacono, AER, Inc.
   !
!----------------------------------------------------------------------
   IMPLICIT NONE
!<DESCRIPTION>
!
! Radiation_driver is the WRF mediation layer routine that provides the interface to
! to radiation physics packages in the WRF model layer. The radiation
! physics packages to call are chosen by setting the namelist variable
! (Rconfig entry in Registry) to the integer value assigned to the 
! particular package (package entry in Registry). For example, if the
! namelist variable ra_lw_physics is set to 1, this corresponds to the
! Registry Package entry for swradscheme.  Note that the Package
! names in the Registry are defined constants (frame/module_state_description.F)
! in the CASE statements in this routine.
!
! Among the arguments is moist, a four-dimensional scalar array storing
! a variable number of moisture tracers, depending on the physics 
! configuration for the WRF run, as determined in the namelist.  The
! highest numbered index of active moisture tracers the integer argument
! n_moist (note: the number of tracers at run time is the quantity
! <tt>n_moist - PARAM_FIRST_SCALAR + 1</tt> , not n_moist. Individual tracers
! may be indexed from moist by the Registry name of the tracer prepended
! with P_; for example P_QC is the index of cloud water. An index 
! represents a valid, active field only if the index is greater than
! or equal to PARAM_FIRST_SCALAR.  PARAM_FIRST_SCALAR and the individual
! indices for each tracer is defined in module_state_description and
! set in <a href=set_scalar_indices_from_config.html>set_scalar_indices_from_config</a> defined in frame/module_configure.F.
!
! Physics drivers in WRF 2.0 and higher, originally model-layer 
! routines, have been promoted to mediation layer routines and they
! contain OpenMP threaded loops over tiles.  Thus, physics drivers
! are called from single-threaded regions in the solver. The physics
! routines that are called from the physics drivers are model-layer
! routines and fully tile-callable and thread-safe.
!</DESCRIPTION>
! 
!======================================================================
! Grid structure in physics part of WRF
!----------------------------------------------------------------------
! The horizontal velocities used in the physics are unstaggered
! relative to temperature/moisture variables. All predicted
! variables are carried at half levels except w, which is at full
! levels. Some arrays with names (*8w) are at w (full) levels.
!
!----------------------------------------------------------------------
! In WRF, kms (smallest number) is the bottom level and kme (largest
! number) is the top level.  In your scheme, if 1 is at the top level,
! then you have to reverse the order in the k direction.
!
!         kme      -   half level (no data at this level)
!         kme    ----- full level
!         kme-1    -   half level
!         kme-1  ----- full level
!         .
!         .
!         .
!         kms+2    -   half level
!         kms+2  ----- full level
!         kms+1    -   half level
!         kms+1  ----- full level
!         kms      -   half level
!         kms    ----- full level
!
!======================================================================
! Grid structure in physics part of WRF
! 
!-------------------------------------
! The horizontal velocities used in the physics are unstaggered 
! relative to temperature/moisture variables. All predicted 
! variables are carried at half levels except w, which is at full 
! levels. Some arrays with names (*8w) are at w (full) levels.
!
!==================================================================
! Definitions
!-----------
! Theta      potential temperature (K)
! Qv         water vapor mixing ratio (kg/kg)
! Qc         cloud water mixing ratio (kg/kg)
! Qr         rain water mixing ratio (kg/kg)
! Qi         cloud ice mixing ratio (kg/kg)
! Qs         snow mixing ratio (kg/kg)
!-----------------------------------------------------------------
!-- PM2_5_DRY     Dry PM2.5 aerosol mass for all species (ug m^-3)
!-- PM2_5_WATER   PM2.5 water mass (ug m^-3)
!-- PM2_5_DRY_EC  Dry PM2.5 elemental carbon aersol mass (ug m^-3)
!-- RTHRATEN      Theta tendency 
!                 due to radiation (K/s)
!-- RTHRATENLW    Theta tendency 
!                 due to long wave radiation (K/s)
!-- RTHRATENSW    Theta temperature tendency 
!                 due to short wave radiation (K/s)
!-- dt            time step (s)
!-- itimestep     number of time steps
!-- GLW           downward long wave flux at ground surface (W/m^2)
!-- GSW           net short wave flux at ground surface (W/m^2)
!-- SWDOWN        downward short wave flux at ground surface (W/m^2)
!-- XLAT          latitude, south is negative (degree)
!-- XLONG         longitude, west is negative (degree)
!-- CLDFRA        cloud fraction (between 0 and 1)
!-- EMISS         surface emissivity (between 0 and 1)
!-- rho_phy       density (kg/m^3)
!-- rr            dry air density (kg/m^3)
!-- moist         moisture array (4D - last index is species) (kg/kg)
!-- n_moist       number of moisture species
!-- qndrop        Cloud droplet number (#/kg)
!-- p8w           pressure at full levels (Pa)
!-- p_phy         pressure (Pa)
!-- pi_phy        exner function (dimensionless)
!-- dz8w          dz between full levels (m)
!-- t_phy         temperature (K)
!-- t8w           temperature at full levels (K)
!-- GMT           Greenwich Mean Time Hour of model start (hour)
!-- JULDAY        the initial day (Julian day)
!-- RADT          time for calling radiation (min)
!-- ra_call_offset -1 (old) means usually just before output, 0 after
!-- DEGRAD        conversion factor for 
!                 degrees to radians (pi/180.) (rad/deg)
!-- DPD           degrees per day for earth's 
!                 orbital position (deg/day)
!-- R_d           gas constant for dry air (J/kg/K)
!-- CP            heat capacity at constant pressure for dry air (J/kg/K)
!-- G             acceleration due to gravity (m/s^2)
!-- rvovrd        R_v divided by R_d (dimensionless)
!-- XTIME         time since simulation start (min)
!-- DECLIN        solar declination angle (rad)
!-- SOLCON        solar constant (W/m^2)
!-- ids           start index for i in domain
!-- ide           end index for i in domain
!-- jds           start index for j in domain
!-- jde           end index for j in domain
!-- kds           start index for k in domain
!-- kde           end index for k in domain
!-- ims           start index for i in memory
!-- ime           end index for i in memory
!-- jms           start index for j in memory
!-- jme           end index for j in memory
!-- kms           start index for k in memory
!-- kme           end index for k in memory
!-- i_start       start indices for i in tile
!-- i_end         end indices for i in tile
!-- j_start       start indices for j in tile
!-- j_end         end indices for j in tile
!-- kts           start index for k in tile
!-- kte           end index for k in tile
!-- num_tiles     number of tiles
!
!==================================================================
!
   INTEGER,      INTENT(IN   )    ::   ids,ide, jds,jde, kds,kde, &
                                       ims,ime, jms,jme, kms,kme, &
                                       its,ite, jts,jte, kts,kte

   INTEGER, INTENT(IN)            :: lw_physics, sw_physics

   INTEGER,      INTENT(IN   )    ::   ICLOUD

   LOGICAL,      INTENT(IN   )    ::   warm_rain

   REAL,      INTENT(IN   )       ::   RADT

   REAL, DIMENSION( ims:ime, jms:jme ),                           &
         INTENT(IN   )  ::                                 XLAND, &
                                                            XICE, &
                                                             TSK, &
                                                          VEGFRA, &
                                                            SNOW 


   REAL, DIMENSION( ims:ime, jms:jme ),                           &
         INTENT(INOUT)  ::                                  HTOP, &
                                                            HBOT

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                  &
         INTENT(IN ) ::                                     dz8w, &
                                                             p8w, &
                                                               p, &
                                                              pi, &
                                                               t, &
                                                             t8w, &
                                                             rho
!
   REAL, DIMENSION( kts:kte, jms:jme,14 ), OPTIONAL ,       &
         INTENT(IN ) ::  ext_cof,sscal,asymp
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), OPTIONAL ,       &
         INTENT(IN ) ::  qc_adjust, qi_adjust

   REAL, DIMENSION( kms:kme, jms:jme, 16 ), OPTIONAL ,       &
         INTENT(IN ) ::  tauaerlw

   LOGICAL, OPTIONAL :: cu_rad_feedback

   INTEGER, INTENT(IN   ), OPTIONAL  ::   aer_ra_feedback

!jdfcz   INTEGER, OPTIONAL, INTENT(IN   )    :: progn,prescribe
   INTEGER, OPTIONAL, INTENT(IN   )    :: progn
!
! variables for aerosols (only if running with chemistry)
!
!
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                  &
         INTENT(INOUT)  ::                              RTHRATEN, &
                                                      RTHRATENLW, &
                                                      RTHRATENSW

!  REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), OPTIONAL ,       &
!        INTENT(INOUT)  ::                                  SWUP, &
!                                                           SWDN, &
!                                                      SWUPCLEAR, &
!                                                      SWDNCLEAR, &
!                                                           LWUP, &
!                                                           LWDN, &
!                                                      LWUPCLEAR, &
!                                                      LWDNCLEAR

   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT) ::&
                      ACSWUPT,ACSWUPTC,ACSWDNT,ACSWDNTC,          &
                      ACSWUPB,ACSWUPBC,ACSWDNB,ACSWDNBC,          &
                      ACLWUPT,ACLWUPTC,ACLWDNT,ACLWDNTC,          &
                      ACLWUPB,ACLWUPBC,ACLWDNB,ACLWDNBC

! TOA and surface, upward and downward, total and clear fluxes
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT) ::&
                        SWUPT,  SWUPTC,  SWDNT,  SWDNTC,          &
                        SWUPB,  SWUPBC,  SWDNB,  SWDNBC,          &
                        LWUPT,  LWUPTC,  LWDNT,  LWDNTC,          &
                        LWUPB,  LWUPBC,  LWDNB,  LWDNBC

! Upward and downward, total and clear sky layer fluxes (W m-2)
   REAL, DIMENSION( ims:ime, kms:kme+2, jms:jme ),                &
         OPTIONAL, INTENT(INOUT) ::                               &
                               SWUPFLX,SWUPFLXC,SWDNFLX,SWDNFLXC, &
                               LWUPFLX,LWUPFLXC,LWDNFLX,LWDNFLXC

   REAL, DIMENSION( ims:ime, jms:jme ),          OPTIONAL ,       &
         INTENT(INOUT)  ::                                  SWCF, &
                                                            LWCF, &
                                                             OLR


!
   REAL, DIMENSION( ims:ime, jms:jme ),                           &
         INTENT(IN   )  ::                                  XLAT, &
                    alb_uvb,alb_uvd,alb_irb,alb_ird,XLONG,EMISS
!
   REAL, DIMENSION( ims:ime, jms:jme ),                           &
         INTENT(INOUT)  ::                                   GSW, &
                                                             GLW

!  REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT)  ::   SWDOWN
!
   REAL, INTENT(IN  )   ::                                GMT,dt, &
                                                   xtime
   REAL   ::     julian
!
   INTEGER, INTENT(IN  ) ::                    JULDAY, itimestep
   REAL, INTENT(IN ),OPTIONAL     ::                    CURR_SECS



!
! Optional 
!
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                  &
         OPTIONAL,                                                &
         INTENT(INOUT) ::                                 CLDFRA

!
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                  &
         OPTIONAL,                                                &
         INTENT(INOUT ) ::                                        &
                                        qv,qc,qr,qi,qs,qg,qndrop

   LOGICAL, OPTIONAL ::     f_qv,f_qc,f_qr,f_qi,f_qs,f_qg,f_qndrop
!
! LOCAL  VAR

   REAL, DIMENSION( ims:ime, jms:jme ) ::             GLAT,GLON
   REAL, DIMENSION( ims:ime, jms:jme ) ::             coszr

   REAL    ::    DECLIN,SOLCON,XXLAT,TLOCTM,XT24
   INTEGER ::    i,j,k,ij
   INTEGER ::    STEPABS
   LOGICAL ::    gfdl_lw,gfdl_sw
   LOGICAL ::    doabsems
!  LOGICAL, EXTERNAL :: wrf_dm_on_monitor

   REAL    ::    OBECL,SINOB,SXLONG,ARG,DECDEG,                  &
                 DJUL,RJUL,ECCFAC
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ) :: qi_temp,qc_temp
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ) :: qi_save,qc_save

   REAL    ::    next_rad_time
   LOGICAL ::    run_param
!------------------------------------------------------------------
! solar related variables are added to declaration
!-------------------------------------------------
!  REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme), INTENT(OUT) :: COSZEN
!  REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme), INTENT(OUT) :: HRANG
   julian=float(julday)
!------------------------------------------------------------------

   if (lw_physics .eq. 0 .and. sw_physics .eq. 0)         return

! ra_call_offset = -1 gives old method where radiation may be called just before output
! ra_call_offset =  0 gives new method where radiation may be called just after output
!                     and is also consistent with removal of offset in new XTIME

!
! Calculate whether or not to run the radiation step.
! If CURR_SECS is passed in, we will calculate based upon that.  If it
!   is not passed in, we'll do the old method of using the time STEP number
!   

! moved up and out of OMP loop because it only needs to be computed once
! and because it is not entirely thread-safe (XT24, TOLOCTM and XXLAT need
! their thread-privacy)  JM 20100217
     CALL radconst(XTIME,DECLIN,SOLCON,JULIAN,               &
                   DEGRAD,DPD                                )


!    IF(PRESENT(coszen).AND.PRESENT(hrang))THEN
! state arrays of hrang and coszen used in surface driver
!      XT24=MOD(XTIME+RADT*0.5,1440.)
!      DO j=jts,jte
!      DO i=its,ite
!         TLOCTM=GMT+XT24/60.+XLONG(I,J)/15.
!         HRANG(I,J)=15.*(TLOCTM-12.)*DEGRAD
!         XXLAT=XLAT(I,J)*DEGRAD
!         COSZEN(I,J)=SIN(XXLAT)*SIN(DECLIN)+COS(XXLAT)*COS(DECLIN)*COS(HRANG(I,J))
!      ENDDO
!      ENDDO
!    ENDIF

!---------------
! initialize data

     DO j=jts,jte
     DO i=its,ite
        GSW(I,J)=0.
        GLW(I,J)=0.
        GLAT(I,J)=XLAT(I,J)*DEGRAD
        GLON(I,J)=XLONG(I,J)*DEGRAD
     ENDDO
     ENDDO

     DO j=jts,jte
     DO k=kts,kte+1
     DO i=its,ite
        RTHRATEN(I,K,J)=0.
        RTHRATENLW(I,K,J)=0.
        RTHRATENSW(I,K,J)=0.
!        SWUP(I,K,J) = 0.0
!        SWDN(I,K,J) = 0.0
!        SWUPCLEAR(I,K,J) = 0.0
!        SWDNCLEAR(I,K,J) = 0.0
!        LWUP(I,K,J) = 0.0
!        LWDN(I,K,J) = 0.0
!        LWUPCLEAR(I,K,J) = 0.0
!        LWDNCLEAR(I,K,J) = 0.0
     ENDDO
     ENDDO
     ENDDO

     IF ( PRESENT( SWUPFLX ) ) THEN
        DO j=jts,jte
        DO k=kts,kte+2
        DO i=its,ite
           SWUPFLX(I,K,J) = 0.0
           SWDNFLX(I,K,J) = 0.0
           SWUPFLXC(I,K,J) = 0.0
           SWDNFLXC(I,K,J) = 0.0
           LWUPFLX(I,K,J) = 0.0
           LWDNFLX(I,K,J) = 0.0
           LWUPFLXC(I,K,J) = 0.0
           LWDNFLXC(I,K,J) = 0.0
        ENDDO
        ENDDO
        ENDDO
     ENDIF

! temporarily modify hydrometeors (currently only done for GD scheme and WRF-Chem)
!
     IF ( PRESENT( cu_rad_feedback ) ) THEN
       IF ( PRESENT( qc ) .AND. PRESENT( qc_adjust ) .AND. cu_rad_feedback ) THEN
          DO j=jts,jte
          DO k=kts,kte
          DO i=its,ite
            qc_save(i,k,j) = qc(i,k,j)
            qc(i,k,j) = qc(i,k,j) + qc_adjust(i,k,j)
          ENDDO
          ENDDO
          ENDDO
       ENDIF
       IF ( PRESENT( qi ) .AND. PRESENT( qi_adjust ) .AND. cu_rad_feedback ) THEN
          DO j=jts,jte
          DO k=kts,kte
          DO i=its,ite
            qi_save(i,k,j) = qi(i,k,j)
            qi(i,k,j) = qi(i,k,j) + qi_adjust(i,k,j)
          ENDDO
          ENDDO
          ENDDO
       ENDIF
     ENDIF


! Fill temporary water variable depending on micro package (tgs 25 Apr 2006)
     if(PRESENT(qc) .and. PRESENT(F_QC)) then
        DO j=jts,jte
        DO k=kts,kte
        DO i=its,ite
           qc_temp(I,K,J)=qc(I,K,J)
        ENDDO
        ENDDO
        ENDDO
     else
        DO j=jts,jte
        DO k=kts,kte
        DO i=its,ite
           qc_temp(I,K,J)=0.
        ENDDO
        ENDDO
        ENDDO
     endif
     if(PRESENT(qr) .and. PRESENT(F_QR)) then
        DO j=jts,jte
        DO k=kts,kte
        DO i=its,ite
           qc_temp(I,K,J) = qc_temp(I,K,J) + qr(I,K,J)
        ENDDO
        ENDDO
        ENDDO
     endif

!---------------

!    lwrad_cldfra_select: SELECT CASE(lw_physics)

!       CASE (RRTMG_LWSCHEME)

     IF ( PRESENT ( CLDFRA ) .AND.                           &
          PRESENT(F_QC) .AND. PRESENT ( F_QI ) ) THEN
        write(0,*) 'do wrf rrtmg '
! Call to cloud fraction routine based on Randall 1994 (Hong Pan 1998)

        CALL cal_cldfra2(CLDFRA,qv,qc,qi,qs,               &
                   F_QV,F_QC,F_QI,F_QS,t,p,                &
                   ids,ide, jds,jde, kds,kde,              &
                   ims,ime, jms,jme, kms,kme,              &
                   its,ite, jts,jte, kts,kte               )

     ENDIF
 
!       CASE DEFAULT

!    IF ( PRESENT ( CLDFRA ) .AND.                           &
!         PRESENT(F_QC) .AND. PRESENT ( F_QI ) ) THEN
!      CALL cal_cldfra(CLDFRA,qc,qi,F_QC,F_QI,               &
!                      ids,ide, jds,jde, kds,kde,            &
!                      ims,ime, jms,jme, kms,kme,            &
!                      its,ite, jts,jte, kts,kte             )
!    ENDIF

!    END SELECT lwrad_cldfra_select    

!    lwrad_select: SELECT CASE(lw_physics)



!       CASE (RRTMG_LWSCHEME)
!            CALL wrf_debug (100, 'CALL rrtmg_lw')

             CALL RRTMG_LWRAD(                                      &
                  RTHRATENLW=RTHRATEN,                              &
                  LWUPT=LWUPT,LWUPTC=LWUPTC,                        &
                  LWDNT=LWDNT,LWDNTC=LWDNTC,                        &
                  LWUPB=LWUPB,LWUPBC=LWUPBC,                        &
                  LWDNB=LWDNB,LWDNBC=LWDNBC,                        &
                  GLW=GLW,OLR=OLR,LWCF=LWCF,                     &
                  EMISS=EMISS,                                      &
                  P8W=p8w,P3D=p,PI3D=pi,DZ8W=dz8w,TSK=tsk,T3D=t,    &
                  T8W=t8w,RHO3D=rho,R=R_d,G=G,                      &
                  ICLOUD=icloud,WARM_RAIN=warm_rain,CLDFRA3D=CLDFRA,&
                  XLAND=XLAND,XICE=XICE,SNOW=SNOW,                  &
                  QV3D=QV,QC3D=QC,QR3D=QR,                          &
                  QI3D=QI,QS3D=QS,QG3D=QG,                          &
                  F_QV=F_QV,F_QC=F_QC,F_QR=F_QR,                    &
                  F_QI=F_QI,F_QS=F_QS,F_QG=F_QG,                    &
                  TAUAERLW=tauaerlw,                                &
                  aer_ra_feedback=aer_ra_feedback,                  &
!jdfcz            progn=progn,prescribe=prescribe,                   &
                  progn=progn,                                       &
!#endif
                  QNDROP3D=qndrop,F_QNDROP=f_qndrop,                &
                  IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde,&
                  IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme,&
                  ITS=its,ITE=ite, JTS=jts,JTE=jte, KTS=kts,KTE=kte,&
                  LWUPFLX=LWUPFLX,LWUPFLXC=LWUPFLXC,                &
                  LWDNFLX=LWDNFLX,LWDNFLXC=LWDNFLXC                 &
                                                                    )

!       CASE DEFAULT
  
!            WRITE( wrf_err_message , * ) 'The longwave option does not exist: lw_physics = ', lw_physics
!            CALL wrf_error_fatal ( wrf_err_message )
           
!    END SELECT lwrad_select    

!    IF (lw_physics .gt. 0 .and. .not.gfdl_lw) THEN
        DO j=jts,jte
        DO k=kts,kte
        DO i=its,ite
           RTHRATENLW(I,K,J)=RTHRATEN(I,K,J)
        ENDDO
        ENDDO
        ENDDO
!    ENDIF


!
!    swrad_cldfra_select: SELECT CASE(sw_physics)

 
!       CASE (RRTMG_SWSCHEME)

     IF ( PRESENT ( CLDFRA ) .AND.                           &
          PRESENT(F_QC) .AND. PRESENT ( F_QI ) ) THEN
! Call to cloud fraction routine based on Randall 1994 (Hong Pan 1998)

        CALL cal_cldfra2(CLDFRA,qv,qc,qi,qs,               &
                   F_QV,F_QC,F_QI,F_QS,t,p,                &
                   ids,ide, jds,jde, kds,kde,              &
                   ims,ime, jms,jme, kms,kme,              &
                   its,ite, jts,jte, kts,kte               )
     ENDIF

!       CASE DEFAULT

!    END SELECT swrad_cldfra_select    

!    swrad_select: SELECT CASE(sw_physics)


!       CASE (RRTMG_SWSCHEME)
!            CALL wrf_debug(100, 'CALL rrtmg_sw')
             CALL RRTMG_SWRAD(                                         &
                     RTHRATENSW=RTHRATENSW,                            &
                     SWUPT=SWUPT,SWUPTC=SWUPTC,                        &
                     SWDNT=SWDNT,SWDNTC=SWDNTC,                        &
                     SWUPB=SWUPB,SWUPBC=SWUPBC,                        &
                     SWDNB=SWDNB,SWDNBC=SWDNBC,                        &
                     SWCF=SWCF,GSW=GSW,                                &
                     XTIME=XTIME,GMT=GMT,XLAT=XLAT,XLONG=XLONG,        &
                     RADT=RADT,DEGRAD=DEGRAD,DECLIN=DECLIN,            &
                     COSZR=COSZR,JULDAY=JULDAY,SOLCON=SOLCON,          &
                     t3d=t,t8w=t8w,TSK=TSK,              &
                     p3d=p,p8w=p8w,pi3d=pi,rho3d=rho,                  &
                     dz8w=dz8w,CLDFRA3D=CLDFRA,R=R_D,G=G,              &
                     ICLOUD=icloud,WARM_RAIN=warm_rain,                &
                     XLAND=XLAND,XICE=XICE,SNOW=SNOW,                  &
                     alb_uvb=alb_uvb,alb_uvd=alb_uvd,                  &
                     alb_irb=alb_irb,alb_ird=alb_ird,    &
                     QV3D=qv,QC3D=qc,QR3D=qr,                          &
                     QI3D=qi,QS3D=qs,QG3D=qg,                          &
                     F_QV=f_qv,F_QC=f_qc,F_QR=f_qr,                    &
                     F_QI=f_qi,F_QS=f_qs,F_QG=f_qg,                    &
                     ext_cof=ext_cof,sscal=sscal,asymp=asymp,          &
                     aer_ra_feedback=aer_ra_feedback,                  &
                     progn=progn,                                      &
                     QNDROP3D=qndrop,F_QNDROP=f_qndrop,                &
                     IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde,&
                     IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme,&
                     ITS=its,ITE=ite, JTS=jts,JTE=jte, KTS=kts,KTE=kte,&
                     SWUPFLX=SWUPFLX,SWUPFLXC=SWUPFLXC,                &
                     SWDNFLX=SWDNFLX,SWDNFLXC=SWDNFLXC                 &
                                                                       )
             DO j=jts,jte
             DO k=kts,kte
             DO i=its,ite
                RTHRATEN(I,K,J)=RTHRATEN(I,K,J)+RTHRATENSW(I,K,J)
             ENDDO
             ENDDO
             ENDDO
            write(0,*)'maxgsw = ',maxval(gsw)


!       CASE DEFAULT

!            WRITE( wrf_err_message , * ) 'The shortwave option does not exist: sw_physics = ', sw_physics
!            CALL wrf_error_fatal ( wrf_err_message )

!    END SELECT swrad_select    


!    IF (sw_physics .gt. 0 .and. .not.gfdl_sw) THEN
        DO j=jts,jte
        DO k=kts,kte
        DO i=its,ite
           RTHRATENSW(I,K,J)=RTHRATEN(I,K,J)-RTHRATENLW(I,K,J)
        ENDDO
        ENDDO
        ENDDO

!       DO j=jts,jte
!       DO i=its,ite
!          SWDOWN(I,J)=GSW(I,J)/(1.-ALBEDO(I,J))
!       ENDDO
!       ENDDO

!    ENDIF

     IF ( PRESENT( cu_rad_feedback ) ) THEN
       IF ( PRESENT( qc  ) .AND. PRESENT( qc_adjust ) .AND. cu_rad_feedback ) THEN
           DO j=jts,jte
           DO k=kts,kte
           DO i=its,ite
             qc(i,k,j) = qc_save(i,k,j)
           ENDDO
           ENDDO
           ENDDO
        ENDIF
        IF ( PRESENT( qi  ) .AND. PRESENT( qi_adjust ) .AND. cu_rad_feedback ) THEN
           DO j=jts,jte
           DO k=kts,kte
           DO i=its,ite
             qi(i,k,j) = qi_save(i,k,j)
           ENDDO
           ENDDO
           ENDDO
        ENDIF
      ENDIF

!  ENDDO
!  !$OMP END PARALLEL DO

!   ENDIF Radiation_step

!    accumulate_lw_select: SELECT CASE(lw_physics)

!    CASE (CAMLWSCHEME,RRTMG_LWSCHEME)
   IF(PRESENT(LWUPTC))THEN
!  !$OMP PARALLEL DO   &
!  !$OMP PRIVATE ( ij ,i,j,k,its,ite,jts,jte)

!  DO ij = 1 , num_tiles
!    its = i_start(ij)
!    ite = i_end(ij)
!    jts = j_start(ij)
!    jte = j_end(ij)

        DO j=jts,jte
        DO i=its,ite
           ACLWUPT(I,J) = ACLWUPT(I,J) + LWUPT(I,J)*DT
           ACLWUPTC(I,J) = ACLWUPTC(I,J) + LWUPTC(I,J)*DT
           ACLWDNT(I,J) = ACLWDNT(I,J) + LWDNT(I,J)*DT
           ACLWDNTC(I,J) = ACLWDNTC(I,J) + LWDNTC(I,J)*DT
           ACLWUPB(I,J) = ACLWUPB(I,J) + LWUPB(I,J)*DT
           ACLWUPBC(I,J) = ACLWUPBC(I,J) + LWUPBC(I,J)*DT
           ACLWDNB(I,J) = ACLWDNB(I,J) + LWDNB(I,J)*DT
           ACLWDNBC(I,J) = ACLWDNBC(I,J) + LWDNBC(I,J)*DT
        ENDDO
        ENDDO
!  ENDDO
!  !$OMP END PARALLEL DO
   ENDIF
!    CASE DEFAULT
!    END SELECT accumulate_lw_select

!    accumulate_sw_select: SELECT CASE(sw_physics)

!    CASE (CAMSWSCHEME,RRTMG_SWSCHEME)
   IF(PRESENT(SWUPTC))THEN
!   !$OMP PARALLEL DO   &
!   !$OMP PRIVATE ( ij ,i,j,k,its,ite,jts,jte)

!   DO ij = 1 , num_tiles
!     its = i_start(ij)
!     ite = i_end(ij)
!     jts = j_start(ij)
!     jte = j_end(ij)

        DO j=jts,jte
        DO i=its,ite
           ACSWUPT(I,J) = ACSWUPT(I,J) + SWUPT(I,J)*DT
           ACSWUPTC(I,J) = ACSWUPTC(I,J) + SWUPTC(I,J)*DT
           ACSWDNT(I,J) = ACSWDNT(I,J) + SWDNT(I,J)*DT
           ACSWDNTC(I,J) = ACSWDNTC(I,J) + SWDNTC(I,J)*DT
           ACSWUPB(I,J) = ACSWUPB(I,J) + SWUPB(I,J)*DT
           ACSWUPBC(I,J) = ACSWUPBC(I,J) + SWUPBC(I,J)*DT
           ACSWDNB(I,J) = ACSWDNB(I,J) + SWDNB(I,J)*DT
           ACSWDNBC(I,J) = ACSWDNBC(I,J) + SWDNBC(I,J)*DT
        ENDDO
        ENDDO
!   ENDDO
!   !$OMP END PARALLEL DO
   ENDIF

!     CASE DEFAULT
!     END SELECT accumulate_sw_select

   END SUBROUTINE radiation_driver


!---------------------------------------------------------------------
!BOP
! !IROUTINE: radconst - compute radiation terms
! !INTERFAC:
   SUBROUTINE radconst(XTIME,DECLIN,SOLCON,JULIAN,                   &
                       DEGRAD,DPD                                    )
!---------------------------------------------------------------------
   IMPLICIT NONE
!---------------------------------------------------------------------

! !ARGUMENTS:
   REAL, INTENT(IN   )      ::       DEGRAD,DPD,XTIME,JULIAN
   REAL, INTENT(OUT  )      ::       DECLIN,SOLCON
   REAL                     ::       OBECL,SINOB,SXLONG,ARG,  &
                                     DECDEG,DJUL,RJUL,ECCFAC
!
! !DESCRIPTION:
! Compute terms used in radiation physics 
!EOP

! for short wave radiation

   DECLIN=0.
   SOLCON=0.

!-----OBECL : OBLIQUITY = 23.5 DEGREE.
        
   OBECL=23.5*DEGRAD
   SINOB=SIN(OBECL)
        
!-----CALCULATE LONGITUDE OF THE SUN FROM VERNAL EQUINOX:
        
   IF(JULIAN.GE.80.)SXLONG=DPD*(JULIAN-80.)
   IF(JULIAN.LT.80.)SXLONG=DPD*(JULIAN+285.)
   SXLONG=SXLONG*DEGRAD
   ARG=SINOB*SIN(SXLONG)
   DECLIN=ASIN(ARG)
   DECDEG=DECLIN/DEGRAD
!----SOLAR CONSTANT ECCENTRICITY FACTOR (PALTRIDGE AND PLATT 1976)
   DJUL=JULIAN*360./365.
   RJUL=DJUL*DEGRAD
   ECCFAC=1.000110+0.034221*COS(RJUL)+0.001280*SIN(RJUL)+0.000719*  &
          COS(2*RJUL)+0.000077*SIN(2*RJUL)
   SOLCON=1370.*ECCFAC
   
   END SUBROUTINE radconst

!---------------------------------------------------------------------
!BOP
! !IROUTINE: cal_cldfra - Compute cloud fraction
! !INTERFACE:
   SUBROUTINE cal_cldfra(CLDFRA,QC,QI,F_QC,F_QI,                     &
          ids,ide, jds,jde, kds,kde,                                 &
          ims,ime, jms,jme, kms,kme,                                 &
          its,ite, jts,jte, kts,kte                                  )
!---------------------------------------------------------------------
   IMPLICIT NONE
!---------------------------------------------------------------------
   INTEGER,  INTENT(IN   )   ::           ids,ide, jds,jde, kds,kde, &
                                          ims,ime, jms,jme, kms,kme, &
                                          its,ite, jts,jte, kts,kte

!
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(OUT  ) ::    &
                                                             CLDFRA

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::    &
                                                                 QI, &
                                                                 QC

   LOGICAL,INTENT(IN) :: F_QC,F_QI

   REAL thresh
   INTEGER:: i,j,k
! !DESCRIPTION:
! Compute cloud fraction from input ice and cloud water fields
! if provided.
!
! Whether QI or QC is active or not is determined from the indices of
! the fields into the 4D scalar arrays in WRF. These indices are
! P_QI and P_QC, respectively, and they are passed in to the routine
! to enable testing to see if QI and QC represent active fields in
! the moisture 4D scalar array carried by WRF.
!
! If a field is active its index will have a value greater than or
! equal to PARAM_FIRST_SCALAR, which is also an input argument to
! this routine.
!EOP
!---------------------------------------------------------------------
     thresh=1.0e-6

     IF ( f_qi .AND. f_qc ) THEN
        DO j = jts,jte
        DO k = kts,kte
        DO i = its,ite
           IF ( QC(i,k,j)+QI(I,k,j) .gt. thresh) THEN
              CLDFRA(i,k,j)=1.
           ELSE
              CLDFRA(i,k,j)=0.
           ENDIF
        ENDDO
        ENDDO
        ENDDO
     ELSE IF ( f_qc ) THEN
        DO j = jts,jte
        DO k = kts,kte
        DO i = its,ite
           IF ( QC(i,k,j) .gt. thresh) THEN
              CLDFRA(i,k,j)=1.
           ELSE
              CLDFRA(i,k,j)=0.
           ENDIF
        ENDDO
        ENDDO
        ENDDO
     ELSE
        DO j = jts,jte
        DO k = kts,kte
        DO i = its,ite
           CLDFRA(i,k,j)=0.
        ENDDO
        ENDDO
        ENDDO
     ENDIF

   END SUBROUTINE cal_cldfra

!BOP
! !IROUTINE: cal_cldfra2 - Compute cloud fraction
! !INTERFACE:
! cal_cldfra_xr - Compute cloud fraction.
! Code adapted from that in module_ra_gfdleta.F in WRF_v2.0.3 by James Done
!!
!!---  Cloud fraction parameterization follows Randall, 1994
!!     (see Hong et al., 1998)
!!     (modified by Ferrier, Feb '02)
!
   SUBROUTINE cal_cldfra2(CLDFRA, QV, QC, QI, QS,                     &
                         F_QV, F_QC, F_QI, F_QS, t_phy, p_phy,       &
          ids,ide, jds,jde, kds,kde,                                 &
          ims,ime, jms,jme, kms,kme,                                 &
          its,ite, jts,jte, kts,kte                                  )
!---------------------------------------------------------------------
   IMPLICIT NONE
!---------------------------------------------------------------------
   INTEGER,  INTENT(IN   )   ::           ids,ide, jds,jde, kds,kde, &
                                          ims,ime, jms,jme, kms,kme, &
                                          its,ite, jts,jte, kts,kte

!
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(OUT  ) ::    &
                                                             CLDFRA

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::    &
                                                                 QV, &
                                                                 QI, &
                                                                 QC, &
                                                                 QS, &
                                                              t_phy, &
                                                              p_phy


   LOGICAL,OPTIONAL,INTENT(IN) :: F_QC,F_QI,F_QV,F_QS

!  REAL thresh
   INTEGER:: i,j,k
   REAL    :: RHUM, tc, esw, esi, weight, qvsw, qvsi, qvs_weight, QIMID, QWMID, QCLD, DENOM, ARG, SUBSAT

   REAL    ,PARAMETER :: ALPHA0=100., GAMMA=0.49, QCLDMIN=1.E-12,    &
                                        PEXP=0.25, RHGRID=1.0
   REAL    , PARAMETER ::  SVP1=0.61078
   REAL    , PARAMETER ::  SVP2=17.2693882
   REAL    , PARAMETER ::  SVPI2=21.8745584
   REAL    , PARAMETER ::  SVP3=35.86
   REAL    , PARAMETER ::  SVPI3=7.66
   REAL    , PARAMETER ::  SVPT0=273.15
   REAL    , PARAMETER ::  r_d = 287.
   REAL    , PARAMETER ::  r_v = 461.6
   REAL    , PARAMETER ::  ep_2=r_d/r_v
! !DESCRIPTION:
! Compute cloud fraction from input ice and cloud water fields
! if provided.
!
! Whether QI or QC is active or not is determined from the indices of
! the fields into the 4D scalar arrays in WRF. These indices are 
! P_QI and P_QC, respectively, and they are passed in to the routine
! to enable testing to see if QI and QC represent active fields in
! the moisture 4D scalar array carried by WRF.
! 
! If a field is active its index will have a value greater than or
! equal to PARAM_FIRST_SCALAR, which is also an input argument to 
! this routine.
!EOP


!-----------------------------------------------------------------------
!---  COMPUTE GRID-SCALE CLOUD COVER FOR RADIATION
!     (modified by Ferrier, Feb '02)
!
!---  Cloud fraction parameterization follows Randall, 1994
!     (see Hong et al., 1998)
!-----------------------------------------------------------------------
! Note: ep_2=287./461.6 Rd/Rv
! Note: R_D=287.

! Alternative calculation for critical RH for grid saturation
!     RHGRID=0.90+.08*((100.-DX)/95.)**.5

! Calculate saturation mixing ratio weighted according to the fractions of
! water and ice.
! Following:
! Murray, F.W. 1966. ``On the computation of Saturation Vapor Pressure''  J. Appl. Meteor.  6 p.204
!    es (in mb) = 6.1078 . exp[ a . (T-273.16)/ (T-b) ]
!
!       over ice        over water
! a =   21.8745584      17.2693882
! b =   7.66            35.86

!---------------------------------------------------------------------

    DO j = jts,jte
    DO k = kts,kte
    DO i = its,ite
      tc         = t_phy(i,k,j) - SVPT0
      esw     = 1000.0 * SVP1 * EXP( SVP2  * tc / ( t_phy(i,k,j) - SVP3  ) )
      esi     = 1000.0 * SVP1 * EXP( SVPI2 * tc / ( t_phy(i,k,j) - SVPI3 ) )
      QVSW = EP_2 * esw / ( p_phy(i,k,j) - esw )
      QVSI = EP_2 * esi / ( p_phy(i,k,j) - esi )

      IF ( PRESENT(F_QI) .and. PRESENT(F_QC) .and. PRESENT(F_QS) ) THEN

! mji - For MP options 2, 4, 6, 7, 8, etc. (qc = liquid, qi = ice, qs = snow)
         IF ( F_QI .and. F_QC .and. F_QS) THEN
            QCLD = QI(i,k,j)+QC(i,k,j)+QS(i,k,j)
            IF (QCLD .LT. QCLDMIN) THEN
               weight = 0.
            ELSE
               weight = (QI(i,k,j)+QS(i,k,j)) / QCLD
            ENDIF
         ENDIF

! mji - For MP options 1 and 3, (qc only)
!  For MP=1, qc = liquid, for MP=3, qc = liquid or ice depending on temperature
         IF ( F_QC .and. .not. F_QI .and. .not. F_QS ) THEN
            QCLD = QC(i,k,j)
            IF (QCLD .LT. QCLDMIN) THEN
               weight = 0.
            ELSE
               if (t_phy(i,k,j) .gt. 273.15) weight = 0.
               if (t_phy(i,k,j) .le. 273.15) weight = 1.
            ENDIF
         ENDIF

! mji - For MP option 5; (qc = liquid, qs = ice)
         IF ( F_QC .and. .not. F_QI .and. F_QS ) THEN

! Mixing ratios of cloud water & total ice (cloud ice + snow).
! Mixing ratios of rain are not considered in this scheme.

           QIMID = QS(i,k,j)
           QWMID = QC(i,k,j)
!
!--- Total "cloud" mixing ratio, QCLD.  Rain is not part of cloud,
!    only cloud water + cloud ice + snow
!
           QCLD=QWMID+QIMID
           IF (QCLD .LT. QCLDMIN) THEN
              weight = 0.
           ELSE
              weight = 1.
           ENDIF
         ENDIF

      ELSE
         CLDFRA(i,k,j)=0.

      ENDIF !  IF ( F_QI .and. F_QC .and. F_QS)


      QVS_WEIGHT = (1-weight)*QVSW + weight*QVSI
      RHUM=QV(i,k,j)/QVS_WEIGHT   !--- Relative humidity
!
!--- Determine cloud fraction (modified from original algorithm)
!
      IF (QCLD .LT. QCLDMIN) THEN
!
!--- Assume zero cloud fraction if there is no cloud mixing ratio
!
        CLDFRA(i,k,j)=0.
      ELSEIF(RHUM.GE.RHGRID)THEN
!
!--- Assume cloud fraction of unity if near saturation and the cloud
!    mixing ratio is at or above the minimum threshold
!
        CLDFRA(i,k,j)=1.
      ELSE
!
!--- Adaptation of original algorithm (Randall, 1994; Zhao, 1995)
!    modified based on assumed grid-scale saturation at RH=RHgrid.
!
        SUBSAT=MAX(1.E-10,RHGRID*QVS_WEIGHT-QV(i,k,j))
        DENOM=(SUBSAT)**GAMMA
        ARG=MAX(-6.9, -ALPHA0*QCLD/DENOM)    ! <-- EXP(-6.9)=.001
! prevent negative values  (new)
        RHUM=MAX(1.E-10, RHUM)
        CLDFRA(i,k,j)=(RHUM/RHGRID)**PEXP*(1.-EXP(ARG))
!!              ARG=-1000*QCLD/(RHUM-RHGRID)
!!              ARG=MAX(ARG, ARGMIN)
!!              CLDFRA(i,k,j)=(RHUM/RHGRID)*(1.-EXP(ARG))
        IF (CLDFRA(i,k,j) .LT. .01) CLDFRA(i,k,j)=0.
      ENDIF          !--- End IF (QCLD .LT. QCLDMIN) ...
    ENDDO          !--- End DO i
    ENDDO          !--- End DO k
    ENDDO          !--- End DO j

   END SUBROUTINE cal_cldfra2



END MODULE module_wrfradiation_driver
