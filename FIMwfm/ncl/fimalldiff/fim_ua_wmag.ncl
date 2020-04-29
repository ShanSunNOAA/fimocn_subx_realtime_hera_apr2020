;
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_LIB/ncarg/nclscripts/csm/contributed.ncl"
modl = getenv("MODL")
;
load "$NCL_HOME/obj/cmap_diff.ncl"
;
begin
;
fname = asciiread("arw_file.txt",-1,"string")
grb_file = addfile(fname,"r")
fname2 = asciiread("arw_file2.txt",-1,"string")
grb_file2 = addfile(fname2,"r")
;
domainname = asciiread("domain.txt",-1,"string")
;
if (domainname .eq. "fim") then
;
  ny = filevardimsizes(grb_file,"lat_4")
  nx = filevardimsizes(grb_file,"lon_4")
;
  intm = grb_file->T_4_SFC@initial_time
  fctm = grb_file->T_4_SFC@forecast_time
  lat = grb_file->lat_4(:)
  lon = grb_file->lon_4(:)
;
end if
;
if (domainname .eq. "236") then
;
  fimsizes = filevardimsizes(grb_file,"gridlat_236")
  ny = fimsizes(0)
  nx = fimsizes(1)
;
  intm = grb_file->T_236_SFC@initial_time
  fctm = grb_file->T_236_SFC@forecast_time
  lat = grb_file->gridlat_236(:,:)
  lon = grb_file->gridlon_236(:,:)
  latin1 = grb_file->gridlon_236@Latin1
  latin2 = grb_file->gridlon_236@Latin2
  lov = grb_file->gridlon_236@Lov
;
end if
;
if (domainname .eq. "africa" .or. domainname .eq. "europe" .or. \
    domainname .eq. "w_pacific" .or. domainname .eq. "floating" .or. \
    domainname .eq. "e_pacific") then
;
  nx = filevardimsizes(grb_file,"g0_lon_2")
  ny = filevardimsizes(grb_file,"g0_lat_1")
;
  intm = grb_file->T_GDS0_SFC@initial_time
  fctm = grb_file->T_GDS0_SFC@forecast_time
  lat = grb_file->g0_lat_1(:)
  lon = grb_file->g0_lon_2(:)
;
end if
;
if (domainname .eq. "219") then
;
  fimsizes = filevardimsizes(grb_file,"gridlat_219")
  ny = fimsizes(0)
  nx = fimsizes(1)
;
  intm = grb_file->T_219_SFC@initial_time
  fctm = grb_file->T_219_SFC@forecast_time
  lat = grb_file->gridlat_219(:,:)
  lon = grb_file->gridlon_219(:,:)
;
end if
;
if (domainname .eq. "28") then
;
  fimsizes = filevardimsizes(grb_file,"gridlat_28")
  ny = fimsizes(0)
  nx = fimsizes(1)
;
  intm = grb_file->T_28_SFC@initial_time
  fctm = grb_file->T_28_SFC@forecast_time
  lat = grb_file->gridlat_28(:,:)
  lon = grb_file->gridlon_28(:,:)
;
end if
;
if (domainname .eq. "244") then
;
  ny = filevardimsizes(grb_file,"lat_244")
  nx = filevardimsizes(grb_file,"lon_244")
;
  intm = grb_file->T_244_SFC@initial_time
  fctm = grb_file->T_244_SFC@forecast_time
  lat = grb_file->lat_244(:)
  lon = grb_file->lon_244(:)
;
end if
;
year = stringtointeger(getSubString(intm,6,9))
mon  = stringtointeger(getSubString(intm,0,1))
day  = stringtointeger(getSubString(intm,3,4))
hour = stringtointeger(getSubString(intm,12,13))
time = (/0.0/)
time@units = "seconds since 1970-01-01 00:00:0.0"
date = ut_inv_calendar(year,mon,day,hour,0,0,time@units, 0)
fcst = (/date + stringtointeger(fctm)*3600/)
fcst@units = "seconds since 1970-01-01 00:00:0.0"
utc_date = ut_calendar(fcst, 0)
year   = floattointeger(utc_date(:,0))    ; Convert to integer for
month  = floattointeger(utc_date(:,1))    ; use in the sprinti function.
day    = floattointeger(utc_date(:,2))
hour   = floattointeger(utc_date(:,3))
minute = floattointeger(utc_date(:,4))
second = utc_date(:,5)
date_str = "Valid " + sprinti("%0.2i/", month) + sprinti("%0.2i/",day) + \
           sprinti("%0.4i", year) + " " + sprinti("%0.2i:00 ", hour)  + \
           "UTC"
;
xdim = nx - 1
ydim = ny - 1
rnx = 1.*xdim
rny = 1.*ydim
if (rny .ge. rnx) then
  vpheight = .88
  vpy = 1.-((1.-vpheight)/2.)
  vpwidth = (vpheight/rny)*rnx
  vpx = (1.-((vpheight/rny)*rnx))/2.
else
  vpwidth = .88
  vpx = (1.-vpwidth)/2.
  vpheight = (vpwidth/rnx)*rny
  vpy = 1.-((1.-((vpwidth/rnx)*rny))/2.)
end if
;
; Create an application object.
;
appid = create "mdss" appClass defaultapp
    "appUsrDir" : "./"
    "appDefaultParent" : True
end create
;
; Set only one of the following three params to 1, depending on the type
; of output desired (NCGM=metacode file, X11=screen, PS=PostScript file)
;
    NCGM=1
    X11=0
    PS=0
;
    if (NCGM .eq. 1) then
;
; Create an ncgmWorkstation object.
;
        wid = create "mdssWork" ncgmWorkstationClass defaultapp
           "wkMetaName" : "./ua_wmag.ncgm"
        end create
    else
      if (X11 .eq. 1) then
;
; Create an XWorkstation object.
;
        wid = create "mdssWork" xWorkstationClass defaultapp
            "wkPause" : True
        end create
      else
        if (PS .eq. 1) then
;
; Create a PSWorkstation object.
;
          wid = create "mdssWork" psWorkstationClass defaultapp
           "wkPSFileName" : "./ua_wmag.ps"
          end create
        end if
      end if
    end if
;
; Assign the colormap to the workstation.
setvalues wid
   "wkColorMap" : gen_colormap()
end setvalues
;
setvalues NhlGetWorkspaceObjectId()
   "wsMaximumSize": 247108864
end setvalues
;
load "$NCL_HOME/obj/fim_maps_blk.ncl"
;
levnam = grb_file->lv_ISBL3
lev = levnam  ; to initialize
do ii = 0, dimsizes(levnam)-1
  lev(ii) = ii
end do
;
do i = 0, 5
;
  if (i .eq. 0) then  ; 850 mb level
    do ii = 0, dimsizes(levnam)-1
      if (levnam(ii) .eq. 850) then
        levelname = levnam(ii)
        level = lev(ii)
      end if
    end do
  end if
  if (i .eq. 1) then  ; 250 mb level
    do ii = 0, dimsizes(levnam)-1
      if (levnam(ii) .eq. 250) then
        levelname = levnam(ii)
        level = lev(ii)
      end if
    end do
  end if
  if (i .eq. 2) then  ; 25 mb level
    do ii = 0, dimsizes(levnam)-1
      if (levnam(ii) .eq. 25) then
        levelname = levnam(ii)
        level = lev(ii)
      end if
    end do
  end if
  if (i .eq. 3) then  ; 20 mb level
    do ii = 0, dimsizes(levnam)-1
      if (levnam(ii) .eq. 20) then
        levelname = levnam(ii)
        level = lev(ii)
      end if
    end do
  end if
  if (i .eq. 4) then  ; 10 mb level
    do ii = 0, dimsizes(levnam)-1
      if (levnam(ii) .eq. 10) then
        levelname = levnam(ii)
        level = lev(ii)
      end if
    end do
  end if
  if (i .eq. 5) then  ; 5 mb level
    do ii = 0, dimsizes(levnam)-1
      if (levnam(ii) .eq. 5) then
        levelname = levnam(ii)
        level = lev(ii)
      end if
    end do
  end if
;
  if (domainname .eq. "fim") then
;
    uvar = grb_file->UW_4_ISBL(level,:,:) * 1.9429 ; convert to knots
    vvar = grb_file->VW_4_ISBL(level,:,:) * 1.9429 ; convert to knots
    uvar_2 = grb_file2->UW_4_ISBL(level,:,:) * 1.9429 ; convert to knots
    vvar_2 = grb_file2->VW_4_ISBL(level,:,:) * 1.9429 ; convert to knots
;
  end if
;
  if (domainname .eq. "236") then
;
    uvar = grb_file->UW_236_ISBL(level,:,:) * 1.9429 ; convert to knots
    vvar = grb_file->VW_236_ISBL(level,:,:) * 1.9429 ; convert to knots
    uvar_2 = grb_file2->UW_236_ISBL(level,:,:) * 1.9429 ; convert to knots
    vvar_2 = grb_file2->VW_236_ISBL(level,:,:) * 1.9429 ; convert to knots
;
  end if
;
  if (domainname .eq. "africa" .or. domainname .eq. "europe" .or. \
      domainname .eq. "w_pacific" .or. domainname .eq. "floating" .or. \
      domainname .eq. "e_pacific") then
;
    uvar = grb_file->UW_GDS0_ISBL(level,:,:) * 1.9429 ; convert to knots
    vvar = grb_file->VW_GDS0_ISBL(level,:,:) * 1.9429 ; convert to knots
    uvar_2 = grb_file2->UW_GDS0_ISBL(level,:,:) * 1.9429 ; convert to knots
    vvar_2 = grb_file2->VW_GDS0_ISBL(level,:,:) * 1.9429 ; convert to knots
;
  end if
;
  if (domainname .eq. "219") then
;
    uvar = grb_file->UW_219_ISBL(level,:,:) * 1.9429 ; convert to knots
    vvar = grb_file->VW_219_ISBL(level,:,:) * 1.9429 ; convert to knots
    uvar_2 = grb_file2->UW_219_ISBL(level,:,:) * 1.9429 ; convert to knots
    vvar_2 = grb_file2->VW_219_ISBL(level,:,:) * 1.9429 ; convert to knots
;
  end if
;
  if (domainname .eq. "28") then
;
    uvar = grb_file->UW_28_ISBL(level,:,:) * 1.9429 ; convert to knots
    vvar = grb_file->VW_28_ISBL(level,:,:) * 1.9429 ; convert to knots
    uvar = -uvar
    vvar = -vvar
    uvar_2 = grb_file2->UW_28_ISBL(level,:,:) * 1.9429 ; convert to knots
    vvar_2 = grb_file2->VW_28_ISBL(level,:,:) * 1.9429 ; convert to knots
    uvar_2 = -uvar_2
    vvar_2 = -vvar_2
;
  end if
;
  if (domainname .eq. "244") then
;
    uvar = grb_file->UW_244_ISBL(level,:,:) * 1.9429 ; convert to knots
    vvar = grb_file->VW_244_ISBL(level,:,:) * 1.9429 ; convert to knots
    uvar_2 = grb_file2->UW_244_ISBL(level,:,:) * 1.9429 ; convert to knots
    vvar_2 = grb_file2->VW_244_ISBL(level,:,:) * 1.9429 ; convert to knots
;
  end if
;
  uv2 = uvar * uvar
  vv2 = vvar * vvar
  wind = sqrt(uv2 + vv2)
  uv2_2 = uvar_2 * uvar_2
  vv2_2 = vvar_2 * vvar_2
  wind_2 = sqrt(uv2_2 + vv2_2)
;
  tempvar = wind
  wind = wind - wind_2
;
  print(levelname+"mb")
  print("wind mag diff")
  print(min(wind))
  print(max(wind))
;
  if (min(wind) .eq. 0 .and. max(wind) .eq. 0) then
    print("constant field of zero - exiting")
    exit
  end if
;
  wind(0,0) = 5.1  ; set this corner pt so that black plots are avoided
                   ; this can happen if all differences are < 5 kts
;
; Create ScalarField and VectorField objects.
;
if (domainname .eq. "fim" .or. domainname .eq. "244") then
;
    tfield = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : wind
        "sfXCStartV" : 0
        "sfYCStartV" : ydim
        "sfXCEndV" : xdim
        "sfYCEndV" : 0
        "sfMissingValueV" : wind@_FillValue
    end create
;
    tempfield = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : tempvar
        "sfXCStartV" : 0
        "sfYCStartV" : ydim
        "sfXCEndV" : xdim
        "sfYCEndV" : 0
        "sfMissingValueV" : tempvar@_FillValue
    end create
;
  end if
;
  if (domainname .eq. "236" .or. domainname .eq. "africa" .or. \
      domainname .eq. "219" .or. domainname .eq. "w_pacific" .or. \
      domainname .eq. "europe" .or. domainname .eq. "floating" .or. \
      domainname .eq. "e_pacific") then
;
    tfield = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : wind
        "sfXCStartV" : 0
        "sfYCStartV" : 0
        "sfXCEndV" : xdim
        "sfYCEndV" : ydim
        "sfMissingValueV" : wind@_FillValue
    end create
;
    tempfield = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : tempvar
        "sfXCStartV" : 0
        "sfYCStartV" : 0
        "sfXCEndV" : xdim
        "sfYCEndV" : ydim
        "sfMissingValueV" : tempvar@_FillValue
    end create
;
  end if
;
  if (domainname .eq. "28") then
;
    tfield = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : wind
        "sfXCStartV" : xdim
        "sfYCStartV" : ydim
        "sfXCEndV" : 0
        "sfYCEndV" : 0
        "sfMissingValueV" : wind@_FillValue
    end create
;
    tempfield = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : tempvar
        "sfXCStartV" : xdim
        "sfYCStartV" : ydim
        "sfXCEndV" : 0
        "sfYCEndV" : 0
        "sfMissingValueV" : tempvar@_FillValue
    end create
;
  end if
;
; Create ContourPlot objects.
;
; Wind 
;
  load "$NCL_HOME/obj/diff_-40-40by5.ncl"
  load "$NCL_HOME/obj/wind_20-180by20_contours.ncl"
;
  if (domainname .eq. "236" .or. domainname .eq. "africa" .or. \
      domainname .eq. "244" .or. domainname .eq. "w_pacific" .or. \
      domainname .eq. "europe" .or. domainname .eq. "28" .or. \
      domainname .eq. "floating" .or. domainname .eq. "e_pacific") then
;
    setvalues cnid
;      "cnRasterModeOn" : False
      "cnRasterSmoothingOn" : True
    end setvalues
;
  end if
;
  titlestring = levelname + "mb Wind (kt) (magnitude diff)"
  load "$NCL_HOME/obj/add_titles.ncl"
;
; Draw map
;
  draw(cnid)
  draw(mpid)
  draw(cnid4)
  frame(wid)
;
end do
;
end
;
