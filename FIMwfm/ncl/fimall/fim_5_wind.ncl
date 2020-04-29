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
; Get the 5mb level index
;
do ll = 0,5
  lnam = "lv_ISBL"+ll
  if (isfilevar(grb_file,lnam)) then
    break
  end if
end do
;
levnam = grb_file->$lnam$
levunits = grb_file->$lnam$@units
if (levunits .eq. "Pa") then
  levnam = levnam / 100
end if
do ii = 0, dimsizes(levnam)-1
  if (levnam(ii) .eq. 5) then
    level = ii
    break
  end if
end do
;
domainname = asciiread("domain.txt",-1,"string")
;
load "$NCL_HOME/obj/fim_grib_names.ncl"
;
intm = grb_file->$TSFC$@initial_time
fctm = grb_file->$TSFC$@forecast_time
fctmunits = grb_file->$TSFC$@forecast_time_units
if (fctmunits .ne. "hours") then
  fctm = fctm * tointeger(fctmunits)
end if
;
if (domainname .eq. "fim") then
;
  ny = filevardimsizes(grb_file,LAT)
  nx = filevardimsizes(grb_file,LON)
;
  lat = grb_file->$LAT$(:)
  lon = grb_file->$LON$(:)
  hvar = grb_file->$GHISBL$(level,:,:) / 10. ; use decameters
  uvar = grb_file->$UISBL$(level,:,:) * 1.9429 ; convert to knots
  vvar = grb_file->$VISBL$(level,:,:) * 1.9429 ; convert to knots
  uv2 = uvar * uvar
  vv2 = vvar * vvar
  tvar = sqrt(uv2 + vv2)
  wstrid = 10
;
end if
;
if (domainname .eq. "129") then
;
  ny = filevardimsizes(grb_file,LAT)
  nx = filevardimsizes(grb_file,LON)
;
  lat = grb_file->$LAT$(:)
  lon = grb_file->$LON$(:)
  hvar = grb_file->$GHISBL$(level,:,:) / 10. ; use decameters
  uvar = grb_file->$UISBL$(level,:,:) * 1.9429 ; convert to knots
  vvar = grb_file->$VISBL$(level,:,:) * 1.9429 ; convert to knots
  uv2 = uvar * uvar
  vv2 = vvar * vvar
  tvar = sqrt(uv2 + vv2)
  wstrid = 20
;
end if
;
if (domainname .eq. "83") then
;
  fimsizes = filevardimsizes(grb_file,LAT)
  ny = fimsizes(0)
  nx = fimsizes(1)
;
  lat = grb_file->$LAT$(:,:)
  lon = grb_file->$LON$(:,:)
  hvar = grb_file->$GHISBL$(level,:,:) / 10. ; use decameters
  uvar = grb_file->$UISBL$(level,:,:) * 1.9429 ; convert to knots
  vvar = grb_file->$VISBL$(level,:,:) * 1.9429 ; convert to knots
  uv2 = uvar * uvar
  vv2 = vvar * vvar
  tvar = sqrt(uv2 + vv2)
  wstrid = 10
;
end if
;
if (domainname .eq. "130") then
;
  fimsizes = filevardimsizes(grb_file,LAT)
  ny = fimsizes(0)
  nx = fimsizes(1)
;
  lat = grb_file->$LAT$(:,:)
  lon = grb_file->$LON$(:,:)
  latin1 = grb_file->$LON$@Latin1
  latin2 = grb_file->$LON$@Latin2
  lov = grb_file->$LON$@Lov
  hvar = grb_file->$GHISBL$(level,:,:) / 10. ; use decameters
  uvar = grb_file->$UISBL$(level,:,:) * 1.9429 ; convert to knots
  vvar = grb_file->$VISBL$(level,:,:) * 1.9429 ; convert to knots
  uv2 = uvar * uvar
  vv2 = vvar * vvar
  tvar = sqrt(uv2 + vv2)
  wstrid = 6
;
end if
;
if (domainname .eq. "174") then
;
  ny = filevardimsizes(grb_file,LAT)
  nx = filevardimsizes(grb_file,LON)
;
  lat = grb_file->$LAT$(::-1)
  lon = grb_file->$LON$(:)
  lonorig = lon
  fliplon = lonFlip(lon)
  hvar = grb_file->$GHISBL$(level,::-1,:) / 10. ; use decameters
  hvarorig = hvar
  fliphvar = lonFlip(hvar)
  uvar = grb_file->$UISBL$(level,::-1,:) * 1.9429 ; convert to knots
  uvarorig = uvar
  flipuvar = lonFlip(uvar)
  vvar = grb_file->$VISBL$(level,::-1,:) * 1.9429 ; convert to knots
  vvarorig = vvar
  flipvvar = lonFlip(vvar)
  uv2 = uvar * uvar
  vv2 = vvar * vvar
  tvar = sqrt(uv2 + vv2)
  tvarorig = tvar
  fliptvar = lonFlip(tvar)
;
end if
;
if (domainname .eq. "236") then
;
  fimsizes = filevardimsizes(grb_file,LAT)
  ny = fimsizes(0)
  nx = fimsizes(1)
;
  lat = grb_file->$LAT$(:,:)
  lon = grb_file->$LON$(:,:)
  latin1 = grb_file->$LON$@Latin1
  latin2 = grb_file->$LON$@Latin2
  lov = grb_file->$LON$@Lov
  hvar = grb_file->$GHISBL$(level,:,:) / 10. ; use decameters
  uvar = grb_file->$UISBL$(level,:,:) * 1.9429 ; convert to knots
  vvar = grb_file->$VISBL$(level,:,:) * 1.9429 ; convert to knots
  uv2 = uvar * uvar
  vv2 = vvar * vvar
  tvar = sqrt(uv2 + vv2)
  wstrid = 3
;
end if
;
if (domainname .eq. "201" .or. domainname .eq. "224" .or. domainname .eq. "242") then
;
  fimsizes = filevardimsizes(grb_file,LAT)
  ny = fimsizes(0)
  nx = fimsizes(1)
;
  lat = grb_file->$LAT$(:,:)
  lon = grb_file->$LON$(:,:)
  hvar = grb_file->$GHISBL$(level,:,:) / 10. ; use decameters
  uvartmp = grb_file->$UISBL$(level,:,:)
  vvartmp = grb_file->$VISBL$(level,:,:)
  rot = grb_file->$ROT$(:,:)
  if (domainname .eq. "224") then
    rot = -rot
    uvartmp = -uvartmp
    vvartmp = -vvartmp
  end if
;
  ugrid = cos(rot)*uvartmp - sin(rot)*vvartmp
  vgrid = sin(rot)*uvartmp + cos(rot)*vvartmp
;
  uvar = ugrid * 1.9429 ; convert to knots
  vvar = vgrid * 1.9429 ; convert to knots
;
  uv2 = uvar * uvar
  vv2 = vvar * vvar
  tvar = sqrt(uv2 + vv2)
  wstrid = 4
;
end if
;
if (domainname .eq. "244") then
;
  ny = filevardimsizes(grb_file,LAT)
  nx = filevardimsizes(grb_file,LON)
;
  lat = grb_file->$LAT$(:)
  lon = grb_file->$LON$(:)
  hvar = grb_file->$GHISBL$(level,:,:) / 10. ; use decameters
  uvar = grb_file->$UISBL$(level,:,:) * 1.9429 ; convert to knots
  vvar = grb_file->$VISBL$(level,:,:) * 1.9429 ; convert to knots
  uv2 = uvar * uvar
  vv2 = vvar * vvar
  tvar = sqrt(uv2 + vv2)
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
           "wkMetaName" : "./5_wind.ncgm"
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
           "wkPSFileName" : "./5_wind.ps"
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
print(levnam(level)+"mb")
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
    wstrid = 7
  else
    wstrid = 5
  end if
;
  vfield = create "VectorField" vectorFieldClass appid
    "vfUDataArray" : uvar
    "vfVDataArray" : vvar
    "vfXCStartSubsetV" : llx
    "vfYCStartSubsetV" : lly
    "vfXCEndSubsetV" : urx
    "vfYCEndSubsetV" : ury
    "vfMissingVValueV" : vvar@_FillValue
    "vfMissingUValueV" : uvar@_FillValue
    "vfXCStride" : wstrid       ; i.e., plot a vector at every (n+1)th point.
    "vfYCStride" : wstrid
  end create
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
  if (levnam(level) .eq. 25 .or. levnam(level) .eq. 20 .or. \
      levnam(level) .eq. 10 .or. levnam(level) .eq. 5) then
;
    load "$NCL_HOME/obj/wind_20-300by20.ncl"
    load "$NCL_HOME/obj/height_120m.ncl"
;
  end if
;
  if (levnam(level) .eq. 250) then
;
    load "$NCL_HOME/obj/wind_20-180by20.ncl"
    load "$NCL_HOME/obj/height_120m.ncl"
;
  end if
;
  if (levnam(level) .eq. 850) then
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
  titlestring = levnam(level) + "mb Wind (kt)"
  load "$NCL_HOME/obj/add_titles.ncl"
  load "$NCL_HOME/obj/add_explabel_vcid.ncl"
;
  draw(cnid)
  draw(mpid)
;  draw(cnid4)
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
end if  ; domainname .eq. "130"
;
if (domainname .eq. "174") then
;
  load "$NCL_HOME/obj/fim_174subdomain_loop.ncl"
;
  if (idom .eq. 0 .or. idom .eq. 2) then
    load "$NCL_HOME/obj/fim_maps_cent0_blk.ncl"
  else
    load "$NCL_HOME/obj/fim_maps_cent180_blk.ncl"
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
  vfield = create "VectorField" vectorFieldClass appid
      "vfUDataArray" : uvar
      "vfVDataArray" : vvar
      "vfXCStartSubsetV" : llx
      "vfYCStartSubsetV" : lly
      "vfXCEndSubsetV" : urx
      "vfYCEndSubsetV" : ury
      "vfMissingVValueV" : vvar@_FillValue
      "vfMissingUValueV" : uvar@_FillValue
      "vfXCStride" : wstrid       ; i.e., plot a vector at every (n+1)th point.
      "vfYCStride" : wstrid
  end create
;
  vcid = create "vectorplot" vectorPlotClass wid
      "vcVectorFieldData":  vfield
      "vcGlyphStyle" : "WindBarb"
      "vcWindBarbLineThicknessF" : 1.2
      "vcWindBarbColor" : 0
      "vcRefAnnoOn" : False
      "vcWindBarbTickSpacingF" : 0.1875
      "pmTickMarkDisplayMode" : "NoCreate"
      "vpXF"      : vpx        ; set viewport
      "vpYF"      : vpy
      "vpWidthF"  : vpwidth
      "vpHeightF" : vpheight
  end create
;
  if (levnam(level) .eq. 25 .or. levnam(level) .eq. 20 .or. \
      levnam(level) .eq. 10 .or. levnam(level) .eq. 5) then
;
    load "$NCL_HOME/obj/wind_20-300by20.ncl"
    load "$NCL_HOME/obj/height_120m.ncl"
;
  end if
;
  if (levnam(level) .eq. 250) then
;
    load "$NCL_HOME/obj/wind_20-180by20.ncl"
    load "$NCL_HOME/obj/height_120m.ncl"
;
  end if
;
  if (levnam(level) .eq. 850) then
;
    load "$NCL_HOME/obj/wind_5-90by5.ncl"
    load "$NCL_HOME/obj/height_60m.ncl"
;
  end if
;
  titlestring = levnam(level) + "mb Wind (kt)"
  load "$NCL_HOME/obj/add_titles.ncl"
  load "$NCL_HOME/obj/add_explabel_vcid.ncl"
;
  draw(cnid)
  draw(mpid)
;  draw(cnid4)
  draw(vcid)
;
  frame(wid)
;
  end do  ; loop begins in fim_174subdomain_loop.ncl
;
end if  ; domainname .eq. "174"
;
if (domainname .eq. "242") then
;
  load "$NCL_HOME/obj/242_tiles_loop.ncl"
;
  centerlat = 90
  centerlon = grb_file->$LON$@Lov
  load "$NCL_HOME/obj/stereographic_map_blk.ncl"
;
  setvalues mpid
    "mpUSStateLineThicknessF" : 2
    "mpGeophysicalLineThicknessF" : 2
    "mpNationalLineThicknessF" : 2
  end setvalues
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
  vfield = create "VectorField" vectorFieldClass appid
      "vfUDataArray" : uvar
      "vfVDataArray" : vvar
      "vfXCStartSubsetV" : llx
      "vfYCStartSubsetV" : lly
      "vfXCEndSubsetV" : urx
      "vfYCEndSubsetV" : ury
      "vfMissingVValueV" : vvar@_FillValue
      "vfMissingUValueV" : uvar@_FillValue
      "vfXCStride" : wstride       ; i.e., plot a vector at every (n+1)th point.
      "vfYCStride" : wstride
  end create
;
  vcid = create "vectorplot" vectorPlotClass wid
      "vcVectorFieldData":  vfield
      "vcGlyphStyle" : "WindBarb"
      "vcWindBarbLineThicknessF" : 1.2
      "vcWindBarbColor" : 0
      "vcRefAnnoOn" : False
      "vcWindBarbTickSpacingF" : 0.1875
      "pmTickMarkDisplayMode" : "NoCreate"
      "vpXF"      : vpx        ; set viewport
      "vpYF"      : vpy
      "vpWidthF"  : vpwidth
      "vpHeightF" : vpheight
  end create
;
  if (levnam(level) .eq. 25 .or. levnam(level) .eq. 20 .or. \
      levnam(level) .eq. 10 .or. levnam(level) .eq. 5) then
;
    load "$NCL_HOME/obj/wind_20-300by20.ncl"
    load "$NCL_HOME/obj/height_120m.ncl"
;
  end if
;
  if (levnam(level) .eq. 250) then
;
    load "$NCL_HOME/obj/wind_20-180by20.ncl"
    load "$NCL_HOME/obj/height_120m.ncl"
;
  end if
;
  if (levnam(level) .eq. 850) then
;
    load "$NCL_HOME/obj/wind_5-90by5.ncl"
    load "$NCL_HOME/obj/height_60m.ncl"
;
  end if
;
  setvalues cnid
    "cnFillMode" : "AreaFill"
  end setvalues
;
  titlestring = levnam(level) + "mb Wind (kt)"
  load "$NCL_HOME/obj/add_titles.ncl"
  load "$NCL_HOME/obj/add_explabel_vcid.ncl"
;
  draw(cnid)
  draw(mpid)
;  draw(cnid4)
  draw(vcid)
;
  frame(wid)
;
  end do  ; loop begins in 242_tiles_loop.ncl
;
end if  ; domainname .eq. "242"
;
if (domainname .ne. "130" .and. domainname .ne. "174" .and. domainname .ne. "242") then
;
  load "$NCL_HOME/obj/fim_maps_blk.ncl"
;
  if (domainname .eq. "fim" .or. domainname .eq. "244" .or. domainname .eq. "129") then
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
  if (domainname .eq. "236" .or. domainname .eq. "201" .or. domainname .eq. "224" .or. \
      domainname .eq. "83") then
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
  if (levnam(level) .eq. 25 .or. levnam(level) .eq. 20 .or. \
      levnam(level) .eq. 10 .or. levnam(level) .eq. 5) then
;
    load "$NCL_HOME/obj/wind_20-300by20.ncl"
    load "$NCL_HOME/obj/height_120m.ncl"
;
  end if
;
  if (levnam(level) .eq. 250) then
;
    load "$NCL_HOME/obj/wind_20-180by20.ncl"
    load "$NCL_HOME/obj/height_120m.ncl"
;
  end if
;
  if (levnam(level) .eq. 850) then
;
    load "$NCL_HOME/obj/wind_5-90by5.ncl"
    load "$NCL_HOME/obj/height_60m.ncl"
;
  end if
;
  if (domainname .eq. "236" .or. domainname .eq. "244" .or. domainname .eq. "83" .or. \
      domainname .eq. "201" .or. domainname .eq. "224") then
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
  if (domainname .eq. "fim" .or. domainname .eq. "129") then
    setvalues cnid4
      "cnLineThicknessF" : 1
      "cnLineLabelDensityF" : 0.5
    end setvalues
  end if
;
  if (domainname .eq. "201") then
    setvalues cnid4
      "cnLineLabelDensityF" : 0.001
    end setvalues
  end if
;
  setvalues vcid
    "vcMinMagnitudeF" : 0
  end setvalues
;
  titlestring = levnam(level) + "mb Wind (kt)"
  load "$NCL_HOME/obj/add_titles.ncl"
  load "$NCL_HOME/obj/add_explabel_vcid.ncl"
;
; Draw map
;
  draw(cnid)
  draw(mpid)
;  draw(cnid4)
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
end if  ; domainname .ne. "130" .and. domainname .ne. "174" .and. domainname .ne. "242"
;
end
;
