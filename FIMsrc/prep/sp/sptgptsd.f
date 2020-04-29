C-----------------------------------------------------------------------
      SUBROUTINE SPTGPTSD(IROMB,MAXWV,KMAX,NMAX,
     &                    KWSKIP,KGSKIP,NRSKIP,NGSKIP,
     &                    RLAT,RLON,WAVE,GP,XP,YP)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:  SPTGPTSD   TRANSFORM SPECTRAL SCALAR TO STATION POINTS
C   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-02-29
C
C ABSTRACT: THIS SUBPROGRAM PERFORMS A SPHERICAL TRANSFORM
C           FROM SPECTRAL COEFFICIENTS OF SCALAR QUANTITIES
C           TO SPECIFIED SETS OF STATION POINT VALUES
C           AND THEIR GRADIENTS ON THE GLOBE.
C           THE WAVE-SPACE CAN BE EITHER TRIANGULAR OR RHOMBOIDAL.
C           THE WAVE AND POINT FIELDS MAY HAVE GENERAL INDEXING,
C           BUT EACH WAVE FIELD IS IN SEQUENTIAL 'IBM ORDER',
C           I.E. WITH ZONAL WAVENUMBER AS THE SLOWER INDEX.
C           THE TRANSFORMS ARE ALL MULTIPROCESSED OVER STATIONS.
C           TRANSFORM SEVERAL FIELDS AT A TIME TO IMPROVE VECTORIZATION.
C           SUBPROGRAM CAN BE CALLED FROM A MULTIPROCESSING ENVIRONMENT.
C
C PROGRAM HISTORY LOG:
C   96-02-29  IREDELL
C 1998-12-15  IREDELL  OPENMP DIRECTIVES INSERTED
C 1999-08-18  IREDELL  OPENMP DIRECTIVE TYPO FIXED 
C
C USAGE:    CALL SPTGPTSD(IROMB,MAXWV,KMAX,NMAX,
C    &                    KWSKIP,KGSKIP,NRSKIP,NGSKIP,
C    &                    RLAT,RLON,WAVE,GP,XP,YP)
C   INPUT ARGUMENTS:
C     IROMB    - INTEGER SPECTRAL DOMAIN SHAPE
C                (0 FOR TRIANGULAR, 1 FOR RHOMBOIDAL)
C     MAXWV    - INTEGER SPECTRAL TRUNCATION
C     KMAX     - INTEGER NUMBER OF FIELDS TO TRANSFORM.
C     NMAX     - INTEGER NUMBER OF STATION POINTS TO RETURN
C     KWSKIP   - INTEGER SKIP NUMBER BETWEEN WAVE FIELDS
C                (DEFAULTS TO (MAXWV+1)*((IROMB+1)*MAXWV+2) IF KWSKIP=0)
C     KGSKIP   - INTEGER SKIP NUMBER BETWEEN STATION POINT SETS
C                (DEFAULTS TO NMAX IF KGSKIP=0)
C     NRSKIP   - INTEGER SKIP NUMBER BETWEEN STATION LATS AND LONS
C                (DEFAULTS TO 1 IF NRSKIP=0)
C     NGSKIP   - INTEGER SKIP NUMBER BETWEEN STATION POINTS
C                (DEFAULTS TO 1 IF NGSKIP=0)
C     RLAT     - REAL (*) STATION LATITUDES IN DEGREES
C     RLON     - REAL (*) STATION LONGITUDES IN DEGREES
C     WAVE     - REAL (*) WAVE FIELDS
C   OUTPUT ARGUMENTS:
C     GP       - REAL (*) STATION POINT SETS
C     XP       - REAL (*) STATION POINT X-GRADIENT SETS
C     YP       - REAL (*) STATION POINT Y-GRADIENT SETS
C
C SUBPROGRAMS CALLED:
C   SPWGET       GET WAVE-SPACE CONSTANTS
C   SPLEGEND     COMPUTE LEGENDRE POLYNOMIALS
C   SPSYNTH      SYNTHESIZE FOURIER FROM SPECTRAL
C   SPGRADY      COMPUTE Y-GRADIENT IN SPECTRAL SPACE
C   SPGRADX      COMPUTE X-GRADIENT IN FOURIER SPACE
C   SPFFTPT      COMPUTE FOURIER TRANSFORM TO GRIDPOINTS
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 77
C
C$$$
      REAL RLAT(*),RLON(*),WAVE(*)
      REAL GP(*),XP(*),YP(*)
      REAL EPS((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EPSTOP(MAXWV+1)
      REAL ENN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL ELONN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL EON((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EONTOP(MAXWV+1)
      INTEGER MP(2*KMAX)
      REAL W((MAXWV+1)*((IROMB+1)*MAXWV+2)/2*2,2*KMAX)
      REAL WTOP(2*(MAXWV+1),2*KMAX)
      REAL PLN((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),PLNTOP(MAXWV+1)
      REAL F(2*MAXWV+2,2,3*KMAX),G(3*KMAX)
      PARAMETER(PI=3.14159265358979)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  CALCULATE PRELIMINARY CONSTANTS
      CALL SPWGET(IROMB,MAXWV,EPS,EPSTOP,ENN1,ELONN1,EON,EONTOP)
      MX=(MAXWV+1)*((IROMB+1)*MAXWV+2)/2
      MXTOP=MAXWV+1
      MDIM=2*MX
      IDIM=2*MAXWV+2
      KW=KWSKIP
      KG=KGSKIP
      NR=NRSKIP
      NG=NGSKIP
      IF(KW.EQ.0) KW=2*MX
      IF(KG.EQ.0) KG=NMAX
      IF(NR.EQ.0) NR=1
      IF(NG.EQ.0) NG=1
      MP(1:KMAX)=10
      MP(KMAX+1:2*KMAX)=1
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  CALCULATE SPECTRAL WINDS
!JR Mods to get omp to work
C$OMP PARALLEL DO PRIVATE(KWS,KS,KY,i)
      DO K=1,KMAX
        KWS=(K-1)*KW
        KS=0*KMAX+K
        KY=1*KMAX+K
        DO I=1,2*MX
          W(I,KS)=WAVE(KWS+I)
        ENDDO
        DO I=1,2*MXTOP
          WTOP(I,KS)=0
        ENDDO
        CALL SPGRADY(IROMB,MAXWV,ENN1,EON,EONTOP,
     &               WAVE(KWS+1),W(1,KY),WTOP(1,KY))
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  CALCULATE STATION FIELDS
C$OMP PARALLEL DO PRIVATE(KS,KY,KX,SLAT1,CLAT1)
C$OMP&            PRIVATE(PLN,PLNTOP,F,G,NK)
      DO N=1,NMAX
        IF(ABS(RLAT((N-1)*NR+1)).GE.89.9995) THEN
          SLAT1=SIGN(1.,RLAT((N-1)*NR+1))
          CLAT1=0.
        ELSE
          SLAT1=SIN(PI/180*RLAT((N-1)*NR+1))
          CLAT1=COS(PI/180*RLAT((N-1)*NR+1))
        ENDIF
        CALL SPLEGEND(IROMB,MAXWV,SLAT1,CLAT1,EPS,EPSTOP,
     &                PLN,PLNTOP)
        CALL SPSYNTH(IROMB,MAXWV,2*MAXWV,IDIM,MDIM,2*MXTOP,2*KMAX,
     &               CLAT1,PLN,PLNTOP,MP,W,WTOP,F)
        CALL SPGRADX(MAXWV,IDIM,KMAX,MP,CLAT1,F(1,1,1),F(1,1,2*KMAX+1))
        CALL SPFFTPT(MAXWV,1,IDIM,1,3*KMAX,RLON((N-1)*NR+1),F,G)
        DO K=1,KMAX
          KS=0*KMAX+K
          KY=1*KMAX+K
          KX=2*KMAX+K
          NK=(N-1)*NG+(K-1)*KG+1
          GP(NK)=G(KS)
          XP(NK)=G(KX)
          YP(NK)=G(KY)
        ENDDO
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
