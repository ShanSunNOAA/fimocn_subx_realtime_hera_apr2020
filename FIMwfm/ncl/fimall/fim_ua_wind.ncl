;
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_LIB/ncarg/nclscripts/csm/contributed.ncl"
modl = getenv("MODL")
;
load "$NCL_HOME/obj/cmap1b.ncl"
;
begin
;
airport_id = asciiread("$NCL_HOME/Airports_id.txt",-1,"string")
airport_loc = asciiread("$NCL_HOME/Airports_locs.txt",(/dimsizes(airport_id),2/),"float")
;
fname = asciiread("arw_file.txt",-1,"string")
grb_file = addfile(fname,"r")
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
  wstrid = 10
;
end if
;
if (domainname .eq. "83") then
;
  ny = filevardimsizes(grb_file,"gridlat_83")
  nx = filevardimsizes(grb_file,"gridlon_83")
;
  intm = grb_file->T_83_SFC@initial_time
  fctm = grb_file->T_83_SFC@forecast_time
  lat = grb_file->gridlat_83(:,:)
  lon = grb_file->gridlon_83(:,:)
  wstrid = 10
;
end if
;
if (domainname .eq. "130") then
;
  fimsizes = filevardimsizes(grb_file,"gridlat_130")
  ny = fimsizes(0)
  nx = fimsizes(1)
;
  intm = grb_file->T_130_SFC@initial_time
  fctm = grb_file->T_130_SFC@forecast_time
  lat = grb_file->gridlat_130(:,:)
  lon = grb_file->gridlon_130(:,:)
  latin1 = grb_file->gridlon_130@Latin1
  latin2 = grb_file->gridlon_130@Latin2
  lov = grb_file->gridlon_130@Lov
  wstrid = 6
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
  wstrid = 3
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
  wstrid = 5
  if (domainname .eq. "floating") then
    wstrid = 2
  end if
  if (domainname .eq. "w_pacific" .or. domainname .eq. "e_pacific") then
    wstrid = 7
  end if
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
  wstrid = 10
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
  wstrid = 2
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
  wstrid = 7
;
end if
;
load "$NCL_HOME/obj/get_date_str_FIM.ncl"
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
           "wkMetaName" : "./ua_wind.ncgm"
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
           "wkPSFileName" : "./ua_wind.ps"
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
    hvar = grb_file->GH_4_ISBL(level,:,:) / 10. ; use decameters
    uvar = grb_file->UW_4_ISBL(level,:,:) * 1.9429 ; convert to knots
    vvar = grb_file->VW_4_ISBL(level,:,:) * 1.9429 ; convert to knots
    uv2 = uvar * uvar
    vv2 = vvar * vvar
    tvar = sqrt(uv2 + vv2)
;
  end if
;
  if (domainname .eq. "83") then
;
    hvar = grb_file->GH_83_ISBL(level,:,:) / 10. ; use decameters
    uvar = grb_file->UW_83_ISBL(level,:,:) * 1.9429 ; convert to knots
    vvar = grb_file->VW_83_ISBL(level,:,:) * 1.9429 ; convert to knots
    uv2 = uvar * uvar
    vv2 = vvar * vvar
    tvar = sqrt(uv2 + vv2)
;
  end if
;
  if (domainname .eq. "130") then
;
    hvar = grb_file->GH_130_ISBL(level,:,:) / 10. ; use decameters
    uvar = grb_file->UW_130_ISBL(level,:,:) * 1.9429 ; convert to knots
    vvar = grb_file->VW_130_ISBL(level,:,:) * 1.9429 ; convert to knots
    uv2 = uvar * uvar
    vv2 = vvar * vvar
    tvar = sqrt(uv2 + vv2)
;
  end if
;
  if (domainname .eq. "236") then
;
    hvar = grb_file->GH_236_ISBL(level,:,:) / 10. ; use decameters
    uvar = grb_file->UW_236_ISBL(level,:,:) * 1.9429 ; convert to knots
    vvar = grb_file->VW_236_ISBL(level,:,:) * 1.9429 ; convert to knots
    uv2 = uvar * uvar
    vv2 = vvar * vvar
    tvar = sqrt(uv2 + vv2)
;
  end if
;
  if (domainname .eq. "africa" .or. domainname .eq. "europe" .or. \
      domainname .eq. "w_pacific" .or. domainname .eq. "floating" .or. \
      domainname .eq. "e_pacific") then
;
    hvar = grb_file->GH_GDS0_ISBL(level,:,:) / 10. ; use decameters
    uvar = grb_file->UW_GDS0_ISBL(level,:,:) * 1.9429 ; convert to knots
    vvar = grb_file->VW_GDS0_ISBL(level,:,:) * 1.9429 ; convert to knots
    uv2 = uvar * uvar
    vv2 = vvar * vvar
    tvar = sqrt(uv2 + vv2)
;
  end if
;
  if (domainname .eq. "219") then
;
    hvar = grb_file->GH_219_ISBL(level,:,:) / 10. ; use decameters
    uvartmp = grb_file->UW_219_ISBL(level,:,:)
    vvartmp = grb_file->VW_219_ISBL(level,:,:)
;
    if (modl .eq. "EXPER_FIMZEUS-G8_C" .or. modl .eq. "EXPER_FIMZEUS-9_C" .or. \
        modl .eq. "EXPER_FIM-95_C" .or. modl .eq. "EXPER_FIM-8_C" .or. \
        modl .eq. "EXPER_FIM-9_C" .or. modl .eq. "FIMX_C") then  ; rotate winds
;
      cdf_file = addfile("gridrot_219.nc","r")    
      rot = cdf_file->gridrot_219(:,:)
;
      ugrid = cos(rot)*uvartmp - sin(rot)*vvartmp
      vgrid = sin(rot)*uvartmp + cos(rot)*vvartmp
;
      uvar = ugrid * 1.9429 ; convert to knots
      vvar = vgrid * 1.9429 ; convert to knots
;
    else
;
      uvar = uvartmp * 1.9429 ; convert to knots
      vvar = vvartmp * 1.9429 ; convert to knots
;
    end if
;
    uv2 = uvar * uvar
    vv2 = vvar * vvar
    tvar = sqrt(uv2 + vv2)
;
  end if
;
  if (domainname .eq. "28") then
;
    hvar = grb_file->GH_28_ISBL(level,:,:) / 10. ; use decameters
    uvar = grb_file->UW_28_ISBL(level,:,:) * 1.9429 ; convert to knots
    vvar = grb_file->VW_28_ISBL(level,:,:) * 1.9429 ; convert to knots
    uvar = -uvar
    vvar = -vvar
    uv2 = uvar * uvar
    vv2 = vvar * vvar
    tvar = sqrt(uv2 + vv2)
;
  end if
;
  if (domainname .eq. "244") then
;
    hvar = grb_file->GH_244_ISBL(level,:,:) / 10. ; use decameters
    uvar = grb_file->UW_244_ISBL(level,:,:) * 1.9429 ; convert to knots
    vvar = grb_file->VW_244_ISBL(level,:,:) * 1.9429 ; convert to knots
    uv2 = uvar * uvar
    vv2 = vvar * vvar
    tvar = sqrt(uv2 + vv2)
;
  end if
;
  print(levelname+"mb")
  print("gh")
  print(min(hvar))
  print(max(hvar))
  print("wind")
  print(min(tvar))
  print(max(tvar))
;
; Create ScalarField and VectorField objects.
;
  if (domainname .eq. "130") then
;
  load "$NCL_HOME/obj/ruc_tiles_loop.ncl"
;
    if (iq .eq. 0) then
      load "$NCL_HOME/obj/lambert_conformal_map_blk.ncl"
      setvalues mpid
        "mpUSStateLineThicknessF" : 2.0
        "mpGeophysicalLineThicknessF" : 2.0
        "mpNationalLineThicknessF" : 2.0
      end setvalues
    else
      load "$NCL_HOME/obj/lc_map_counties.ncl"
      setvalues mpid
        "mpNationalLineColor" : "foreground"
      end setvalues
    end if
;
    tfield = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : tvar
        "sfXCStartSubsetV" : llx
        "sfYCStartSubsetV" : lly
        "sfXCEndSubsetV" : urx
        "sfYCEndSubsetV" : ury
        "sfMissingValueV" : tvar@_FillValue
    end create
;
    hfield = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : hvar
        "sfXCStartSubsetV" : llx
        "sfYCStartSubsetV" : lly
        "sfXCEndSubsetV" : urx
        "sfYCEndSubsetV" : ury
        "sfMissingValueV" : hvar@_FillValue
    end create
;
    if (iq .eq. 0) then
      vfield = create "VectorField" vectorFieldClass appid
        "vfUDataArray" : uvar
        "vfVDataArray" : vvar
        "vfXCStartSubsetV" : llx
        "vfYCStartSubsetV" : lly
        "vfXCEndSubsetV" : urx
        "vfYCEndSubsetV" : ury
        "vfMissingVValueV" : vvar@_FillValue
        "vfMissingUValueV" : uvar@_FillValue
        "vfXCStride" : 7          ; i.e., plot a vector at every (n+1)th point.
        "vfYCStride" : 7
      end create
    else
      vfield = create "VectorField" vectorFieldClass appid
        "vfUDataArray" : uvar
        "vfVDataArray" : vvar
        "vfXCStartSubsetV" : llx
        "vfYCStartSubsetV" : lly
        "vfXCEndSubsetV" : urx
        "vfYCEndSubsetV" : ury
        "vfMissingVValueV" : vvar@_FillValue
        "vfMissingUValueV" : uvar@_FillValue
        "vfXCStride" : 5          ; i.e., plot a vector at every (n+1)th point.
        "vfYCStride" : 5
      end create
    end if
;
    vcid = create "vectorplot" vectorPlotClass wid
        "vcVectorFieldData":  vfield
        "vcGlyphStyle" : "WindBarb"
        "vcWindBarbLineThicknessF" : 1.2
        "vcWindBarbColor" : 0
        "vcRefAnnoOn" : False
        "vcWindBarbTickSpacingF" : 0.1875
        "vcMinMagnitudeF" : 0
        "pmTickMarkDisplayMode" : "NoCreate"
        "vpXF"      : vpx        ; set viewport
        "vpYF"      : vpy
        "vpWidthF"  : vpwidth
        "vpHeightF" : vpheight
    end create
;
    if (levelname .eq. 25 .or. levelname .eq. 20 .or. \
        levelname .eq. 10 .or. levelname .eq. 5) then
;
      load "$NCL_HOME/obj/wind_20-300by20.ncl"
      load "$NCL_HOME/obj/height_120m.ncl"
;
    end if
;
    if (levelname .eq. 250) then
;
      load "$NCL_HOME/obj/wind_20-180by20.ncl"
      load "$NCL_HOME/obj/height_120m.ncl"
;
    end if
;
    if (levelname .eq. 850) then
;
      load "$NCL_HOME/obj/wind_5-90by5.ncl"
      load "$NCL_HOME/obj/height_60m.ncl"
;
    end if
;
    setvalues cnid
      "cnRasterModeOn" : False
    end setvalues
;
    if (iq .eq. 0) then
      setvalues cnid4
        "cnLineThicknessF" : 2.0
      end setvalues
    else
      setvalues cnid4
        "cnLineThicknessF" : 3.0
      end setvalues
    end if
;
    titlestring = levelname + "mb Wind (kt)"
    load "$NCL_HOME/obj/add_titles.ncl"
    load "$NCL_HOME/obj/add_explabel_vcid.ncl"
;
    draw(cnid)
    draw(mpid)
    draw(cnid4)
;
    gsres = True
    if (iq .eq. 0) then
      do ii = 0,dimsizes(airport_id)-1
        gsres@gsMarkerIndex = 16
        gsres@gsMarkerColor = 20
        gsres@gsMarkerSizeF = 0.004
        gsn_polymarker(wid,mpid,airport_loc(ii,1),airport_loc(ii,0),gsres)
        gsres@gsMarkerIndex = 16
        gsres@gsMarkerColor = 1
        gsres@gsMarkerSizeF = 0.002
        gsn_polymarker(wid,mpid,airport_loc(ii,1),airport_loc(ii,0),gsres)
      end do
    else
      do ii = 0,dimsizes(airport_id)-1
        gsres@gsMarkerIndex = 16
        gsres@gsMarkerColor = 20
        gsres@gsMarkerSizeF = 0.01
        gsn_polymarker(wid,mpid,airport_loc(ii,1),airport_loc(ii,0),gsres)
        gsres@gsMarkerIndex = 16
        gsres@gsMarkerColor = 1
        gsres@gsMarkerSizeF = 0.005
        gsn_polymarker(wid,mpid,airport_loc(ii,1),airport_loc(ii,0),gsres)
      end do
    end if
;
    draw(vcid)
    frame(wid)
;
    end do  ; loop begins in ruc_tiles_loop.ncl
;
  else  ; domain is not 130
;
    load "$NCL_HOME/obj/fim_maps_blk.ncl"
;
    if (domainname .eq. "fim" .or. domainname .eq. "244") then
;
      tfield = create "ScalarField" scalarFieldClass appid
          "sfDataArray" : tvar
          "sfXCStartV" : 0
          "sfYCStartV" : ydim
          "sfXCEndV" : xdim
          "sfYCEndV" : 0
          "sfMissingValueV" : tvar@_FillValue
      end create
;
      hfield = create "ScalarField" scalarFieldClass appid
          "sfDataArray" : hvar
          "sfXCStartV" : 0
          "sfYCStartV" : ydim
          "sfXCEndV" : xdim
          "sfYCEndV" : 0
          "sfMissingValueV" : hvar@_FillValue
      end create
;
      vfield = create "VectorField" vectorFieldClass appid
          "vfUDataArray" : uvar
          "vfVDataArray" : vvar
          "vfXCStartV" : 0
          "vfYCStartV" : ydim
          "vfXCEndV" : xdim
          "vfYCEndV" : 0
          "vfMissingVValueV" : vvar@_FillValue
          "vfMissingUValueV" : uvar@_FillValue
          "vfXCStride" : wstrid      ; i.e., plot a vector at every (n+1)th point.
          "vfYCStride" : wstrid 
      end create
;
    end if
;
    if (domainname .eq. "236" .or. domainname .eq. "africa" .or. \
        domainname .eq. "219" .or. domainname .eq. "w_pacific" .or. \
        domainname .eq. "europe" .or. domainname .eq. "floating" .or. \
        domainname .eq. "e_pacific" .or. domainname .eq. "83") then
;
      tfield = create "ScalarField" scalarFieldClass appid
          "sfDataArray" : tvar
          "sfXCStartV" : 0
          "sfYCStartV" : 0
          "sfXCEndV" : xdim
          "sfYCEndV" : ydim
          "sfMissingValueV" : tvar@_FillValue
      end create
;
      hfield = create "ScalarField" scalarFieldClass appid
          "sfDataArray" : hvar
          "sfXCStartV" : 0
          "sfYCStartV" : 0
          "sfXCEndV" : xdim
          "sfYCEndV" : ydim
          "sfMissingValueV" : hvar@_FillValue
      end create
;
      vfield = create "VectorField" vectorFieldClass appid
          "vfUDataArray" : uvar
          "vfVDataArray" : vvar
          "vfXCStartV" : 0
          "vfYCStartV" : 0
          "vfXCEndV" : xdim
          "vfYCEndV" : ydim
          "vfMissingVValueV" : vvar@_FillValue
          "vfMissingUValueV" : uvar@_FillValue
          "vfXCStride" : wstrid      ; i.e., plot a vector at every (n+1)th point.
          "vfYCStride" : wstrid 
      end create
;
    end if
;
    if (domainname .eq. "28") then
;
      tfield = create "ScalarField" scalarFieldClass appid
          "sfDataArray" : tvar
          "sfXCStartV" : xdim
          "sfYCStartV" : ydim
          "sfXCEndV" : 0
          "sfYCEndV" : 0
          "sfMissingValueV" : tvar@_FillValue
      end create
;
      hfield = create "ScalarField" scalarFieldClass appid
          "sfDataArray" : hvar
          "sfXCStartV" : xdim
          "sfYCStartV" : ydim
          "sfXCEndV" : 0
          "sfYCEndV" : 0
          "sfMissingValueV" : hvar@_FillValue
      end create
;
      vfield = create "VectorField" vectorFieldClass appid
          "vfUDataArray" : uvar
          "vfVDataArray" : vvar
          "vfXCStartV" : xdim
          "vfYCStartV" : ydim
          "vfXCEndV" : 0
          "vfYCEndV" : 0
          "vfMissingVValueV" : vvar@_FillValue
          "vfMissingUValueV" : uvar@_FillValue
          "vfXCStride" : wstrid       ; i.e., plot a vector at every (n+1)th point.
          "vfYCStride" : wstrid
      end create
;
    end if
;
    vcid = create "vectorplot" vectorPlotClass wid
        "vcVectorFieldData":  vfield
        "vcGlyphStyle" : "WindBarb"
        "vcWindBarbLineThicknessF" : 1.2
        "vcWindBarbColor" : 0
        "vcRefAnnoOn" : False
        "vcMinMagnitudeF" : 0.0
        "vcWindBarbTickSpacingF" : 0.1875
        "pmTickMarkDisplayMode" : "NoCreate"
        "vpXF"      : vpx        ; set viewport
        "vpYF"      : vpy       
        "vpWidthF"  : vpwidth
        "vpHeightF" : vpheight
    end create
;
; Create ContourPlot objects.
;
    if (levelname .eq. 25 .or. levelname .eq. 20 .or. \
        levelname .eq. 10 .or. levelname .eq. 5) then
;
      load "$NCL_HOME/obj/wind_20-300by20.ncl"
      load "$NCL_HOME/obj/height_120m.ncl"
;
    end if
;
    if (levelname .eq. 250) then
;
      load "$NCL_HOME/obj/wind_20-180by20.ncl"
      load "$NCL_HOME/obj/height_120m.ncl"
;
    end if
;
    if (levelname .eq. 850) then
;
      load "$NCL_HOME/obj/wind_5-90by5.ncl"
      load "$NCL_HOME/obj/height_60m.ncl"
;
    end if
;
    if (domainname .eq. "236" .or. domainname .eq. "africa" .or. \
        domainname .eq. "244" .or. domainname .eq. "w_pacific" .or. \
        domainname .eq. "europe" .or. domainname .eq. "28" .or. \
        domainname .eq. "floating" .or. domainname .eq. "e_pacific" .or. \
        domainname .eq. "83") then
;
      setvalues cnid
        "cnRasterModeOn" : False
      end setvalues
;
    end if
;
    setvalues cnid4
      "cnLineThicknessF" : 2.0
    end setvalues
;
    if (domainname .eq. "fim") then
      setvalues cnid4
        "cnLineThicknessF" : 1
        "cnLineLabelDensityF" : 0.001
      end setvalues
    end if
;
    if (domainname .eq. "219") then
      setvalues cnid4
        "cnLineLabelDensityF" : 0.001
      end setvalues
    end if
;
    setvalues vcid
      "vcMinMagnitudeF" : 0
    end setvalues
;
    titlestring = levelname + "mb Wind (kt)"
    load "$NCL_HOME/obj/add_titles.ncl"
    load "$NCL_HOME/obj/add_explabel_vcid.ncl"
;
; Draw map
;
    draw(cnid)
    draw(mpid)
    draw(cnid4)
;
    gsres = True
    if (domainname .eq. "236") then
      do ii = 0,dimsizes(airport_id)-1
        gsres@gsMarkerIndex = 16
        gsres@gsMarkerColor = 20
        gsres@gsMarkerSizeF = 0.004
        gsn_polymarker(wid,mpid,airport_loc(ii,1),airport_loc(ii,0),gsres)
        gsres@gsMarkerIndex = 16
        gsres@gsMarkerColor = 1
        gsres@gsMarkerSizeF = 0.002
        gsn_polymarker(wid,mpid,airport_loc(ii,1),airport_loc(ii,0),gsres)
      end do
    end if
;
    draw(vcid)
    frame(wid)
;
  end if  ; domainname .eq. "130"
;
end do ; i = 0,5 (isobaric levels)
;
end
;
