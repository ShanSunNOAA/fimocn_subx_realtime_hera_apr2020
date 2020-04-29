;
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_LIB/ncarg/nclscripts/csm/contributed.ncl"
modl = getenv("MODL")
;
load "$NCL_HOME/obj/cmap_diff.ncl"
;
begin
;
airport_id = asciiread("$NCL_HOME/Airports_id.txt",-1,"string")
airport_loc = asciiread("$NCL_HOME/Airports_locs.txt",(/dimsizes(airport_id),2/),"float")
;
fname = asciiread("arw_file.txt",-1,"string")
grb_file = addfile(fname,"r")
fname2 = asciiread("arw_file2.txt",-1,"string")
grb_file2 = addfile(fname2,"r")
;
domainname = asciiread("domain.txt",-1,"string")
;
load "$NCL_HOME/obj/fim_grib_names.ncl"
;
; Get the 25mb level index
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
  if (levnam(ii) .eq. 25) then
    level = ii
  end if
end do
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
  uvar = grb_file->$UISBL$(level,:,:) * 1.9429 ; convert to knots
  vvar = grb_file->$VISBL$(level,:,:) * 1.9429 ; convert to knots
  uvar2 = grb_file2->$UISBL$(level,:,:) * 1.9429 ; convert to knots
  vvar2 = grb_file2->$VISBL$(level,:,:) * 1.9429 ; convert to knots
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
  uvar = grb_file->$UISBL$(level,:,:) * 1.9429 ; convert to knots
  vvar = grb_file->$VISBL$(level,:,:) * 1.9429 ; convert to knots
  uvar2 = grb_file2->$UISBL$(level,:,:) * 1.9429 ; convert to knots
  vvar2 = grb_file2->$VISBL$(level,:,:) * 1.9429 ; convert to knots
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
  uvar = grb_file->$UISBL$(level,:,:) * 1.9429 ; convert to knots
  vvar = grb_file->$VISBL$(level,:,:) * 1.9429 ; convert to knots
  uvar2 = grb_file2->$UISBL$(level,:,:) * 1.9429 ; convert to knots
  vvar2 = grb_file2->$VISBL$(level,:,:) * 1.9429 ; convert to knots
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
  uvar = grb_file->$UISBL$(level,:,:) * 1.9429 ; convert to knots
  vvar = grb_file->$VISBL$(level,:,:) * 1.9429 ; convert to knots
  uvar2 = grb_file2->$UISBL$(level,:,:) * 1.9429 ; convert to knots
  vvar2 = grb_file2->$VISBL$(level,:,:) * 1.9429 ; convert to knots
  wstrid = 6
;
end if
;
if (domainname .eq. "174") then
;
  ny = filevardimsizes(grb_file,LAT)
  nx = filevardimsizes(grb_file,LON)
;
  lat = grb_file->$LAT$(:)
  lon = grb_file->$LON$(:)
  lonorig = lon
  fliplon = lonFlip(lon)
  uvar = grb_file->$UISBL$(level,::-1,:) * 1.9429 ; convert to knots
  vvar = grb_file->$VISBL$(level,::-1,:) * 1.9429 ; convert to knots
  uvar2 = grb_file2->$UISBL$(level,::-1,:) * 1.9429 ; convert to knots
  vvar2 = grb_file2->$VISBL$(level,::-1,:) * 1.9429 ; convert to knots
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
  uvar = grb_file->$UISBL$(level,:,:) * 1.9429 ; convert to knots
  vvar = grb_file->$VISBL$(level,:,:) * 1.9429 ; convert to knots
  uvar2 = grb_file2->$UISBL$(level,:,:) * 1.9429 ; convert to knots
  vvar2 = grb_file2->$VISBL$(level,:,:) * 1.9429 ; convert to knots
  wstrid = 3
;
end if
;
if (domainname .eq. "201" .or. domainname .eq. "224") then
;
  fimsizes = filevardimsizes(grb_file,LAT)
  ny = fimsizes(0)
  nx = fimsizes(1)
;
  lat = grb_file->$LAT$(:,:)
  lon = grb_file->$LON$(:,:)
  uvartmp = grb_file->$UISBL$(level,:,:) * 1.9429 ; convert to knots
  vvartmp = grb_file->$VISBL$(level,:,:) * 1.9429 ; convert to knots
  rot = grb_file->$ROT$(:,:)
;
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
  uvartmp2 = grb_file2->$UISBL$(level,:,:) * 1.9429 ; convert to knots
  vvartmp2 = grb_file2->$VISBL$(level,:,:) * 1.9429 ; convert to knots
  rot2 = grb_file2->$ROT$(:,:)
;
  ugrid2 = cos(rot2)*uvartmp2 - sin(rot2)*vvartmp2
  vgrid2 = sin(rot2)*uvartmp2 + cos(rot2)*vvartmp2
;
  uvar2 = ugrid2 * 1.9429 ; convert to knots
  vvar2 = vgrid2 * 1.9429 ; convert to knots
;
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
  uvar = grb_file->$UISBL$(level,:,:) * 1.9429 ; convert to knots
  vvar = grb_file->$VISBL$(level,:,:) * 1.9429 ; convert to knots
  uvar2 = grb_file2->$UISBL$(level,:,:) * 1.9429 ; convert to knots
  vvar2 = grb_file2->$VISBL$(level,:,:) * 1.9429 ; convert to knots
  wstrid = 7
;
end if
;
load "$NCL_HOME/obj/get_date_str_FIM.ncl"
;
if (domainname .eq. "174") then
  uvar = uvar - uvar2
  uvarorig = uvar
  uvarorig2 = uvar2
  flipuvar = lonFlip(uvar)
  flipuvar2 = lonFlip(uvar2)
  vvar = vvar - vvar2
  vvarorig = vvar
  vvarorig2 = vvar2
  flipvvar = lonFlip(vvar)
  flipvvar2 = lonFlip(vvar2)
  uv2 = uvar * uvar
  vv2 = vvar * vvar
  tvar = sqrt(uv2 + vv2)
  tvarorig = tvar
  fliptvar = lonFlip(tvar)
else
  tempuvar = uvar
  tempuvar2 = uvar2
  tempvvar = vvar
  tempvvar2 = vvar2
  uvar = uvar - uvar2
  vvar = vvar - vvar2
  uv2 = uvar * uvar
  vv2 = vvar * vvar
  tvar = sqrt(uv2 + vv2)
  tvarorig = tvar
end if
;
print(levnam(level)+"mb vector wind diff")
print(min(tvar))
print(max(tvar))
;
if (min(tvar) .eq. 0 .and. max(tvar) .eq. 0) then
  print("constant field of zero - exiting")
  exit
end if
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
           "wkMetaName" : "./25_wind.ncgm"
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
           "wkPSFileName" : "./25_wind.ps"
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
if (domainname .eq. "130") then
;
  load "$NCL_HOME/obj/ruc_tiles_loop.ncl"
;
  tvar = tvarorig  ; this will undo any corner points set to avoid black plots
;
  if (iq .eq. 0) then
    load "$NCL_HOME/obj/lambert_conformal_map_blk.ncl"
    setvalues mpid
      "mpGeophysicalLineColor" : "background"
      "mpGeophysicalLineThicknessF" : 2
      "mpNationalLineColor" : "background"
      "mpNationalLineThicknessF" : 2
      "mpUSStateLineColor" : "background"
      "mpUSStateLineThicknessF" : 2
    end setvalues
  else
    load "$NCL_HOME/obj/lc_map_counties.ncl"
    setvalues mpid
      "mpGeophysicalLineColor" : "background"
      "mpGeophysicalLineThicknessF" : 2
      "mpNationalLineColor" : "background"
      "mpNationalLineThicknessF" : 2
      "mpUSStateLineColor" : "background"
      "mpUSStateLineThicknessF" : 2
      "mpCountyLineColor" : "background"
      "mpCountyLineThicknessF" : 1
    end setvalues
  end if
;
  if (max(tvar(lly:ury,llx:urx)) .le. 5.0 ) then
    tvar(lly+1,llx+1) = 5.2  ; set a corner pt to avoid black plot
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
  if (iq .eq. 0) then
    barbthk = 1
  else
    barbthk = 2
  end if
;
  vcid = create "vectorplot" vectorPlotClass wid
      "vcVectorFieldData":  vfield
      "vcGlyphStyle" : "WindBarb"
      "vcWindBarbLineThicknessF" : barbthk
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
  load "$NCL_HOME/obj/diff_-40-40by5.ncl"
;
  setvalues cnid
    "cnRasterModeOn" : False
  end setvalues
;
  titlestring = "25mb Wind (kt)"
  load "$NCL_HOME/obj/add_titles.ncl"
  load "$NCL_HOME/obj/add_explabel_vcid.ncl"
;
; Draw map
;
  draw(cnid)
  draw(mpid)
  draw(vcid)
;
  gsres = True
  if (iq .eq. 0) then
    do i = 0,dimsizes(airport_id)-1
      gsres@gsMarkerIndex = 16
      gsres@gsMarkerColor = 19
      gsres@gsMarkerSizeF = 0.004
      gsn_polymarker(wid,mpid,airport_loc(i,1),airport_loc(i,0),gsres)
      gsres@gsMarkerIndex = 16
      gsres@gsMarkerColor = 1
      gsres@gsMarkerSizeF = 0.002
      gsn_polymarker(wid,mpid,airport_loc(i,1),airport_loc(i,0),gsres)
    end do
  else
    do i = 0,dimsizes(airport_id)-1
      gsres@gsMarkerIndex = 16
      gsres@gsMarkerColor = 19
      gsres@gsMarkerSizeF = 0.01
      gsn_polymarker(wid,mpid,airport_loc(i,1),airport_loc(i,0),gsres)
      gsres@gsMarkerIndex = 16
      gsres@gsMarkerColor = 1
      gsres@gsMarkerSizeF = 0.005
      gsn_polymarker(wid,mpid,airport_loc(i,1),airport_loc(i,0),gsres)
    end do
  end if
;
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
  load "$NCL_HOME/obj/diff_-40-40by5.ncl"
;
  titlestring = "25mb Wind (kt)"
  load "$NCL_HOME/obj/add_titles.ncl"
  load "$NCL_HOME/obj/add_explabel_vcid.ncl"
;
  draw(cnid)
  draw(mpid)
  draw(vcid)
;
  frame(wid)
;
  end do  ; loop begins in fim_174subdomain_loop.ncl
;
end if  ; domainname .eq. "174"
;
if (domainname .ne. "130" .and. domainname .ne. "174") then
;
  load "$NCL_HOME/obj/fim_maps_blk.ncl"
  setvalues mpid
    "mpGeophysicalLineColor" : "background"
    "mpGeophysicalLineThicknessF" : 2
    "mpNationalLineColor" : "background"
    "mpNationalLineThicknessF" : 2
    "mpUSStateLineColor" : "background"
    "mpUSStateLineThicknessF" : 2
  end setvalues
;
; Create a ScalarField object.
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
  if (domainname .eq. "236" .or. domainname .eq. "83" .or. \
      domainname .eq. "201" .or. domainname .eq. "224") then
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
; Create a VectorPlot object.
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
; Create a ContourPlot object.
;
  load "$NCL_HOME/obj/diff_-40-40by5.ncl"
;
  if (domainname .eq. "236" .or. domainname .eq. "244" .or. domainname .eq. "83" .or. \
      domainname .eq. "201" .or. domainname .eq. "224") then
;
    setvalues cnid
;      "cnRasterModeOn" : False
      "cnRasterSmoothingOn" : True
    end setvalues
;
  end if
;
  titlestring = "25mb Wind (kt)"
  load "$NCL_HOME/obj/add_titles.ncl"
  load "$NCL_HOME/obj/add_explabel_vcid.ncl"
;
; Draw map
;
  draw(cnid)
  draw(mpid)
  draw(vcid)
  frame(wid)
;
end if  ; domainname .ne. "130" .and. domainname .ne. "174"
;
end
;
