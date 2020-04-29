;
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_LIB/ncarg/nclscripts/csm/contributed.ncl"
modl = getenv("MODL")
interval = getenv("FCST_INTERVAL")
interval = 6
;
load "$NCL_HOME/obj/cmap1_track2.ncl"
;
begin
;
airport_id = asciiread("$NCL_HOME/Airports_id.txt",-1,"string")
airport_loc = asciiread("$NCL_HOME/Airports_locs.txt",(/dimsizes(airport_id),2/),"float")
;
fname = asciiread("arw_file.txt",-1,"string") 
grb_file = addfile(fname,"r")
;
; Get the 500mb and 1000mb level indexes
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
  if (levnam(ii) .eq. 500) then
    level500 = ii
  end if
  if (levnam(ii) .eq. 1000) then
    level1000 = ii
  end if
end do
;
domainname = asciiread("domain.txt",-1,"string")
;
load "$NCL_HOME/obj/fim_grib_names.ncl"
;
if (modl .ne. "MPAS") then
  intm = grb_file->$TSFC$@initial_time
  fctm = grb_file->$TSFC$@forecast_time
  fctmunits = grb_file->$TSFC$@forecast_time_units
  if (fctmunits .ne. "hours") then
    fctm = fctm * tointeger(fctmunits)
  end if
end if
;
if (domainname .eq. "fim") then
;
  if (modl .eq. "MPAS") then
;
    ny = filevardimsizes(grb_file,"g0_lat_0")
    nx = filevardimsizes(grb_file,"g0_lon_1")
;
    intm = grb_file->PRES_GDS0_SFC@initial_time
    fctm = grb_file->PRES_GDS0_SFC@forecast_time
    lat = grb_file->g0_lat_0(:)
    lon = grb_file->g0_lon_1(:)
;
    if ( fctm .gt. 0 ) then
      tvar = grb_file->ACPCP_GDS0_SFC(::-1,:)
      tvar = tvar * 0.03937 ; convert to inches
      tvar = lonFlip(tvar)
    else
      exit
    end if
;
    mvar = grb_file->PRMSL_GDS0_MSL(::-1,:) / 100. ; convert from Pa to mb
    mvar = lonFlip(mvar)
    mint = 8
    mdens = 0.8
    mlinethick = 1.0
;
    h5 = grb_file->HGT_GDS0_ISBL(level500,::-1,:) / 10. ; use decameters
    h10 = grb_file->HGT_GDS0_ISBL(level1000,::-1,:) / 10. ; use decameters
    thick = h5 - h10
    thick = lonFlip(thick)
    thint = 6
    thdens = 0.1
    thlinethick = 1.0
;
  else
;
    ny = filevardimsizes(grb_file,LAT)
    nx = filevardimsizes(grb_file,LON)
;
    lat = grb_file->$LAT$(:)
    lon = grb_file->$LON$(:)
    if ( fctm .gt. 0 ) then
      if (isfilevar(grb_file,CPACC)) then
        zz = CPACC
      else
        zz = CPACC + interval + "h"
      end if
      tvar = grb_file->$zz$(:,:) * 0.03937 ; convert to inches
    else
      exit
    end if
;
    t32 = grb_file->$T2M$(:,:) - 273.15
    t32 = t32 * 9/5 + 32 ; Convert to F
;
    mvar = grb_file->$MSLP$(:,:) / 100. ; convert from Pa to mb
    mint = 8
    mdens = 0.8
    mlinethick = 1.0
;
    h5 = grb_file->$GHISBL$(level500,:,:) / 10. ; use decameters
    h10 = grb_file->$GHISBL$(level1000,:,:) / 10. ; use decameters
    thick = h5 - h10
    thint = 6
    thdens = 0.1
    thlinethick = 1.0
;
  end if
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
  if ( fctm .gt. 0 ) then
    if (isfilevar(grb_file,CPACC)) then
      zz = CPACC
    else
      zz = CPACC + interval + "h"
    end if
    tvar = grb_file->$zz$(:,:) * 0.03937 ; convert to inches
  else
    exit
  end if
;
  t32 = grb_file->$T2M$(:,:) - 273.15
  t32 = t32 * 9/5 + 32 ; Convert to F
;
  mvar = grb_file->$MSLP$(:,:) / 100. ; convert from Pa to mb
  mint = 8
  mdens = 0.8
  mlinethick = 1.0
;
  h5 = grb_file->$GHISBL$(level500,:,:) / 10. ; use decameters
  h10 = grb_file->$GHISBL$(level1000,:,:) / 10. ; use decameters
  thick = h5 - h10
  thint = 6
  thdens = 0.1
  thlinethick = 1.0
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
  if ( fctm .gt. 0 ) then
    if (isfilevar(grb_file,CPACC)) then
      zz = CPACC
    else
      zz = CPACC + interval + "h"
    end if
    tvar = grb_file->$zz$(:,:) * 0.03937 ; convert to inches
  else
    exit
  end if
;
  t32 = grb_file->$T2M$(:,:) - 273.15
  t32 = t32 * 9/5 + 32 ; Convert to F
;
  mvar = grb_file->$MSLP$(:,:) / 100. ; convert from Pa to mb
  mint = 8
  mdens = 0.8
  mlinethick = 1.0
;
  h5 = grb_file->$GHISBL$(level500,:,:) / 10. ; use decameters
  h10 = grb_file->$GHISBL$(level1000,:,:) / 10. ; use decameters
  thick = h5 - h10
  thint = 6
  thdens = 0.1
  thlinethick = 1.0
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

  if ( fctm .gt. 0 ) then
    if (isfilevar(grb_file,CPACC)) then
      zz = CPACC
    else
      zz = CPACC + interval + "h"
    end if
    tvar = grb_file->$zz$(:,:) * 0.03937 ; convert to inches
  else
    exit
  end if
;
  t32 = grb_file->$T2M$(:,:) - 273.15
  t32 = t32 * 9/5 + 32 ; Convert to F
;
  mvar = grb_file->$MSLP$(:,:) / 100. ; convert from Pa to mb
  mint = 2
  mdens = 0.5
  mlinethick = 1.0
;
  h5 = grb_file->$GHISBL$(level500,:,:) / 10. ; use decameters
  h10 = grb_file->$GHISBL$(level1000,:,:) / 10. ; use decameters
  thick = h5 - h10
  thint = 6
  thdens = 0.5
  thlinethick = 1.0
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
  if ( fctm .gt. 0 ) then
    if (isfilevar(grb_file,CPACC)) then
      zz = CPACC
    else
      zz = CPACC + interval + "h"
    end if
    tvar = grb_file->$zz$(::-1,:) * 0.03937 ; convert to inches
  else
    exit
  end if
;
  tvarorig = tvar
  fliptvar = lonFlip(tvar)
;
  t32 = grb_file->$T2M$(::-1,:) - 273.15
  t32 = t32 * 9/5 + 32 ; Convert to F
  t32orig = t32
  flipt32 = lonFlip(t32)
;
  mvar = grb_file->$MSLP$(::-1,:) / 100. ; convert from Pa to mb
  mvarorig = mvar
  flipmvar = lonFlip(mvar)
  mint = 8
  mdens = 1.0
  mlinethick = 1.0
;
  h5 = grb_file->$GHISBL$(level500,::-1,:) / 10. ; use decameters
  h10 = grb_file->$GHISBL$(level1000,::-1,:) / 10. ; use decameters
  thick = h5 - h10
  thickorig = thick
  flipthick = lonFlip(thick)
  thint = 6
  thdens = 1.0
  thlinethick = 1.0
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

  if ( fctm .gt. 0 ) then
    if (isfilevar(grb_file,CPACC)) then
      zz = CPACC
    else
      zz = CPACC + interval + "h"
    end if
    tvar = grb_file->$zz$(:,:) * 0.03937 ; convert to inches
  else
    exit
  end if
;
  t32 = grb_file->$T2M$(:,:) - 273.15
  t32 = t32 * 9/5 + 32 ; Convert to F
;
  mvar = grb_file->$MSLP$(:,:) / 100. ; convert from Pa to mb
  mint = 2
  mdens = 0.5
  mlinethick = 1.0
;
  h5 = grb_file->$GHISBL$(level500,:,:) / 10. ; use decameters
  h10 = grb_file->$GHISBL$(level1000,:,:) / 10. ; use decameters
  thick = h5 - h10
  thint = 6
  thdens = 0.5
  thlinethick = 1.0
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
  if ( fctm .gt. 0 ) then
    if (isfilevar(grb_file,CPACC)) then
      zz = CPACC
    else
      zz = CPACC + interval + "h"
    end if
    tvar = grb_file->$zz$(:,:) * 0.03937 ; convert to inches
  else
    exit
  end if
;
  t32 = grb_file->$T2M$(:,:) - 273.15
  t32 = t32 * 9/5 + 32 ; Convert to F
;
  mvar = grb_file->$MSLP$(:,:) / 100. ; convert from Pa to mb
  mint = 4
  mdens = 0.5
  mlinethick = 1.0
;
  h5 = grb_file->$GHISBL$(level500,:,:) / 10. ; use decameters
  h10 = grb_file->$GHISBL$(level1000,:,:) / 10. ; use decameters
  thick = h5 - h10
  thint = 6
  thdens = 0.5
  thlinethick = 1.0
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
  if ( fctm .gt. 0 ) then
    if (isfilevar(grb_file,CPACC)) then
      zz = CPACC
    else
      zz = CPACC + interval + "h"
    end if
    tvar = grb_file->$zz$(:,:) * 0.03937 ; convert to inches
  else
    exit
  end if
;
  t32 = grb_file->$T2M$(:,:) - 273.15
  t32 = t32 * 9/5 + 32 ; Convert to F
;
  mvar = grb_file->$MSLP$(:,:) / 100. ; convert from Pa to mb
  mint = 4
  mdens = 0.5
  mlinethick = 1.0
;
  h5 = grb_file->$GHISBL$(level500,:,:) / 10. ; use decameters
  h10 = grb_file->$GHISBL$(level1000,:,:) / 10. ; use decameters
  thick = h5 - h10
  thint = 6
  thdens = 0.5
  thlinethick = 1.0
;
end if
;
load "$NCL_HOME/obj/get_date_str_FIM.ncl"
;
print(date_str)
print("acpcp")
print(min(tvar))
print(max(tvar))
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
           "wkMetaName" : "./sfc_acpcp.ncgm"
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
           "wkPSFileName" : "./sfc_acpcp.ps"
          end create
        end if
      end if
    end if
;
; Assign the colormap to the workstation.
;
setvalues wid
   "wkColorMap" : gen_colormap()
end setvalues
;
setvalues NhlGetWorkspaceObjectId()
   "wsMaximumSize": 247108864
end setvalues
;
; Create a ScalarField object.
;
if (domainname .eq. "130") then
;
  load "$NCL_HOME/obj/ruc_tiles_loop.ncl"
;
  if (iq .eq. 0) then
    load "$NCL_HOME/obj/lambert_conformal_map_blk.ncl"
  else
    load "$NCL_HOME/obj/lc_map_counties.ncl"
    setvalues mpid
      "mpNationalLineColor" : "background"
      "mpCountyLineColor" : "background"
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
  t32field = create "ScalarField" scalarFieldClass appid
      "sfDataArray" : t32
      "sfXCStartSubsetV" : llx
      "sfYCStartSubsetV" : lly
      "sfXCEndSubsetV" : urx
      "sfYCEndSubsetV" : ury
      "sfMissingValueV" : t32@_FillValue
  end create
;
  mfield = create "ScalarField" scalarFieldClass appid
      "sfDataArray" : mvar
      "sfXCStartSubsetV" : llx
      "sfYCStartSubsetV" : lly
      "sfXCEndSubsetV" : urx
      "sfYCEndSubsetV" : ury
      "sfMissingValueV" : mvar@_FillValue
  end create
;
  zfield = create "ScalarField" scalarFieldClass appid
      "sfDataArray" : thick
      "sfXCStartSubsetV" : llx
      "sfYCStartSubsetV" : lly
      "sfXCEndSubsetV" : urx
      "sfYCEndSubsetV" : ury
      "sfMissingValueV" : thick@_FillValue
  end create
;
; 3-h convective precip
;
  load "$NCL_HOME/obj/precip_0-2_whtbkg.ncl"
;
  setvalues cnid
    "cnRasterModeOn" : False
  end setvalues
;
; 32F isotherm
;
  load "$NCL_HOME/obj/32F_isotherm.ncl"
;
  setvalues cnid5
      "cnLineColors" : (/24/)
  end setvalues
;
  if (iq .ne. 0) then
    mlinethick = 2.0
    thlinethick = 2.0
    setvalues cnid5
      "cnLineThicknessF" : 3.0
    end setvalues
  end if
;
  cnid4 = create "contourplot" contourPlotClass wid
      "cnScalarFieldData":  mfield
      "cnFillOn" : False
      "cnLevelSelectionMode" : "ManualLevels"
      "cnMinLevelValF" : 900
      "cnMaxLevelValF" : 1100
      "cnLevelSpacingF" : mint
      "cnLevelFlags" : (/"LineAndLabel"/)
      "cnInfoLabelOn" : False
      "cnLinesOn" : True
      "cnLineColor" : 23
      "cnLineLabelsOn" : True
      "cnLineLabelInterval" : 1
      "cnLineLabelFont" : 21
      "cnLineLabelFontHeightF" : 0.01
      "cnLineLabelFontColor" : 23 ; "Background"
      "cnLineLabelBackgroundColor" : 1
      "cnLineLabelPlacementMode" : "Computed"
      "cnLineLabelDensityF" : mdens
      "cnLineThicknessF" : mlinethick
      "pmTickMarkDisplayMode" : "NoCreate"
      "vpXF"      : vpx        ; set viewport
      "vpYF"      : vpy
      "vpWidthF"  : vpwidth
      "vpHeightF" : vpheight
  end create
;
  cnid6 = create "contourplot" contourPlotClass wid
      "cnScalarFieldData":  zfield
      "cnFillOn" : False
      "cnLevelSelectionMode" : "ManualLevels"
      "cnMinLevelValF" : 402
      "cnMaxLevelValF" : 600
      "cnLevelSpacingF" : thint
      "cnLevelFlags" : (/"LineAndLabel"/)
      "cnInfoLabelOn" : False
      "cnLinesOn" : True
      "cnLineColor" : 2
;      "cnLineDashPattern" : 5
      "cnLineLabelsOn" : True
      "cnLineLabelInterval" : 1
      "cnLineLabelFont" : 21
      "cnLineLabelFontHeightF" : 0.01
      "cnLineLabelFontColor" : 2 ; "Foreground"
      "cnLineLabelBackgroundColor" : 1
      "cnLineLabelPlacementMode" : "Computed"
      "cnLineLabelDensityF" : thdens
      "cnLineThicknessF" : thlinethick
      "pmTickMarkDisplayMode" : "NoCreate"
      "vpXF"      : vpx        ; set viewport
      "vpYF"      : vpy
      "vpWidthF"  : vpwidth
      "vpHeightF" : vpheight
  end create
;
  titlestring = interval + "h Conv Precip (in), MSLP (mb), 1000-500 Thickness (dm)"
  load "$NCL_HOME/obj/add_titles.ncl"
  load "$NCL_HOME/obj/add_explabel.ncl"
;
; Draw map
;
  draw(cnid)
  draw(cnid5)
  draw(cnid6)
  draw(cnid4)
  draw(mpid)
;
  gsres = True
  if (iq .eq. 0) then
    do i = 0,dimsizes(airport_id)-1
      gsres@gsMarkerIndex = 16
      gsres@gsMarkerColor = 20
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
      gsres@gsMarkerColor = 20
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
  t32field = create "ScalarField" scalarFieldClass appid
      "sfDataArray" : t32
      "sfXCStartSubsetV" : llx
      "sfYCStartSubsetV" : lly
      "sfXCEndSubsetV" : urx
      "sfYCEndSubsetV" : ury
      "sfMissingValueV" : t32@_FillValue
  end create
;
  mfield = create "ScalarField" scalarFieldClass appid
      "sfDataArray" : mvar
      "sfXCStartSubsetV" : llx
      "sfYCStartSubsetV" : lly
      "sfXCEndSubsetV" : urx
      "sfYCEndSubsetV" : ury
      "sfMissingValueV" : mvar@_FillValue
  end create
;
  zfield = create "ScalarField" scalarFieldClass appid
      "sfDataArray" : thick
      "sfXCStartSubsetV" : llx
      "sfYCStartSubsetV" : lly
      "sfXCEndSubsetV" : urx
      "sfYCEndSubsetV" : ury
      "sfMissingValueV" : thick@_FillValue
  end create
;
; 3-h convective precip
;
  load "$NCL_HOME/obj/precip_0-2_whtbkg.ncl"
;
  setvalues cnid
    "cnRasterModeOn" : False
  end setvalues
;
; 32F isotherm
;
  load "$NCL_HOME/obj/32F_isotherm.ncl"
;
  setvalues cnid5
      "cnLineColors" : (/24/)
  end setvalues
;
  if (idom .eq. 0 .or. idom .eq. 1 .or. idom .eq. 4) then
    mdens = 3.0
  else
    mdens = 1.0
  end if
;
  cnid4 = create "contourplot" contourPlotClass wid
      "cnScalarFieldData":  mfield
      "cnFillOn" : False
      "cnLevelSelectionMode" : "ManualLevels"
      "cnMinLevelValF" : 900
      "cnMaxLevelValF" : 1100
      "cnLevelSpacingF" : mint
      "cnLevelFlags" : (/"LineAndLabel"/)
      "cnInfoLabelOn" : False
      "cnLinesOn" : True
      "cnLineColor" : 23
      "cnLineLabelsOn" : True
      "cnLineLabelInterval" : 1
      "cnLineLabelFont" : 21
      "cnLineLabelFontHeightF" : 0.01
      "cnLineLabelFontColor" : 23 ; "Background"
      "cnLineLabelBackgroundColor" : 1
      "cnLineLabelPlacementMode" : "Computed"
      "cnLineLabelDensityF" : mdens
      "cnLineThicknessF" : mlinethick
      "pmTickMarkDisplayMode" : "NoCreate"
      "vpXF"      : vpx        ; set viewport
      "vpYF"      : vpy
      "vpWidthF"  : vpwidth
      "vpHeightF" : vpheight
  end create
;
  cnid6 = create "contourplot" contourPlotClass wid
      "cnScalarFieldData":  zfield
      "cnFillOn" : False
      "cnLevelSelectionMode" : "ManualLevels"
      "cnMinLevelValF" : 402
      "cnMaxLevelValF" : 600
      "cnLevelSpacingF" : thint
      "cnLevelFlags" : (/"LineAndLabel"/)
      "cnInfoLabelOn" : False
      "cnLinesOn" : True
      "cnLineColor" : 2
;      "cnLineDashPattern" : 5
      "cnLineLabelsOn" : True
      "cnLineLabelInterval" : 1
      "cnLineLabelFont" : 21
      "cnLineLabelFontHeightF" : 0.01
      "cnLineLabelFontColor" : 2 ; "Foreground"
      "cnLineLabelBackgroundColor" : 1
      "cnLineLabelPlacementMode" : "Computed"
      "cnLineLabelDensityF" : thdens
      "cnLineThicknessF" : thlinethick
      "pmTickMarkDisplayMode" : "NoCreate"
      "vpXF"      : vpx        ; set viewport
      "vpYF"      : vpy
      "vpWidthF"  : vpwidth
      "vpHeightF" : vpheight
  end create
;
  titlestring = interval + "h Conv Precip (in), MSLP (mb), 1000-500 Thickness (dm)"
  load "$NCL_HOME/obj/add_titles.ncl"
  load "$NCL_HOME/obj/add_explabel.ncl"
;
; Draw map
;
  draw(cnid)
  draw(cnid5)
  draw(cnid6)
  draw(cnid4)
  draw(mpid)
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
  t32field = create "ScalarField" scalarFieldClass appid
      "sfDataArray" : t32
      "sfXCStartSubsetV" : llx
      "sfYCStartSubsetV" : lly
      "sfXCEndSubsetV" : urx
      "sfYCEndSubsetV" : ury
      "sfMissingValueV" : t32@_FillValue
  end create
;
  mfield = create "ScalarField" scalarFieldClass appid
      "sfDataArray" : mvar
      "sfXCStartSubsetV" : llx
      "sfYCStartSubsetV" : lly
      "sfXCEndSubsetV" : urx
      "sfYCEndSubsetV" : ury
      "sfMissingValueV" : mvar@_FillValue
  end create
;
  zfield = create "ScalarField" scalarFieldClass appid
      "sfDataArray" : thick
      "sfXCStartSubsetV" : llx
      "sfYCStartSubsetV" : lly
      "sfXCEndSubsetV" : urx
      "sfYCEndSubsetV" : ury
      "sfMissingValueV" : thick@_FillValue
  end create
;
; 3-h convective precip
;
  load "$NCL_HOME/obj/precip_0-2_whtbkg.ncl"
;
  setvalues cnid
    "cnFillMode" : "AreaFill"
  end setvalues
;
  load "$NCL_HOME/obj/32F_isotherm.ncl"
;
  setvalues cnid5
      "cnLineColors" : (/24/)
  end setvalues
;
  mdens = 1.0
;
  cnid4 = create "contourplot" contourPlotClass wid
      "cnScalarFieldData":  mfield
      "cnFillOn" : False
      "cnLevelSelectionMode" : "ManualLevels"
      "cnMinLevelValF" : 900
      "cnMaxLevelValF" : 1100
      "cnLevelSpacingF" : mint
      "cnLevelFlags" : (/"LineAndLabel"/)
      "cnInfoLabelOn" : False
      "cnLinesOn" : True
      "cnLineColor" : 23
      "cnLineLabelsOn" : True
      "cnLineLabelInterval" : 1
      "cnLineLabelFont" : 21
      "cnLineLabelFontHeightF" : 0.01
      "cnLineLabelFontColor" : 23 ; "Background"
      "cnLineLabelBackgroundColor" : 1
      "cnLineLabelPlacementMode" : "Computed"
      "cnLineLabelDensityF" : mdens
      "cnLineThicknessF" : mlinethick
      "pmTickMarkDisplayMode" : "NoCreate"
      "vpXF"      : vpx        ; set viewport
      "vpYF"      : vpy
      "vpWidthF"  : vpwidth
      "vpHeightF" : vpheight
  end create
;
  cnid6 = create "contourplot" contourPlotClass wid
      "cnScalarFieldData":  zfield
      "cnFillOn" : False
      "cnLevelSelectionMode" : "ManualLevels"
      "cnMinLevelValF" : 402
      "cnMaxLevelValF" : 600
      "cnLevelSpacingF" : thint
      "cnLevelFlags" : (/"LineAndLabel"/)
      "cnInfoLabelOn" : False
      "cnLinesOn" : True
      "cnLineColor" : 2
;      "cnLineDashPattern" : 5
      "cnLineLabelsOn" : True
      "cnLineLabelInterval" : 1
      "cnLineLabelFont" : 21
      "cnLineLabelFontHeightF" : 0.01
      "cnLineLabelFontColor" : 2 ; "Foreground"
      "cnLineLabelBackgroundColor" : 1
      "cnLineLabelPlacementMode" : "Computed"
      "cnLineLabelDensityF" : thdens
      "cnLineThicknessF" : thlinethick
      "pmTickMarkDisplayMode" : "NoCreate"
      "vpXF"      : vpx        ; set viewport
      "vpYF"      : vpy
      "vpWidthF"  : vpwidth
      "vpHeightF" : vpheight
  end create
;
  titlestring = interval + "h Conv Precip (in), MSLP (mb), 1000-500 Thickness (dm)"
  load "$NCL_HOME/obj/add_titles.ncl"
  load "$NCL_HOME/obj/add_explabel.ncl"
;
; Draw map
;
  draw(cnid)
  draw(cnid5)
  draw(cnid6)
  draw(cnid4)
  draw(mpid)
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
; Create ScalarField objects.
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
    if (modl .ne. "MPAS") then
      t32field = create "ScalarField" scalarFieldClass appid
          "sfDataArray" : t32
          "sfXCStartV" : 0
          "sfYCStartV" : ydim
          "sfXCEndV" : xdim
          "sfYCEndV" : 0
          "sfMissingValueV" : t32@_FillValue
      end create
    end if
;
    mfield = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : mvar
        "sfXCStartV" : 0
        "sfYCStartV" : ydim
        "sfXCEndV" : xdim
        "sfYCEndV" : 0
        "sfMissingValueV" : mvar@_FillValue
    end create
;
    zfield = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : thick
        "sfXCStartV" : 0
        "sfYCStartV" : ydim
        "sfXCEndV" : xdim
        "sfYCEndV" : 0
        "sfMissingValueV" : thick@_FillValue
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
    t32field = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : t32
        "sfXCStartV" : 0
        "sfYCStartV" : 0
        "sfXCEndV" : xdim
        "sfYCEndV" : ydim
        "sfMissingValueV" : t32@_FillValue
    end create
;
    mfield = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : mvar
        "sfXCStartV" : 0
        "sfYCStartV" : 0
        "sfXCEndV" : xdim
        "sfYCEndV" : ydim
        "sfMissingValueV" : mvar@_FillValue
    end create
;
    zfield = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : thick
        "sfXCStartV" : 0
        "sfYCStartV" : 0
        "sfXCEndV" : xdim
        "sfYCEndV" : ydim
        "sfMissingValueV" : thick@_FillValue
    end create
;
  end if
;
  load "$NCL_HOME/obj/precip_0-2_whtbkg.ncl"
;
; 32F isotherm
;
  if (modl .ne. "MPAS") then
    load "$NCL_HOME/obj/32F_isotherm.ncl"
;
    setvalues cnid5
        "cnLineColors" : (/24/)
    end setvalues
  end if
;
  cnid4 = create "contourplot" contourPlotClass wid
      "cnScalarFieldData":  mfield
      "cnFillOn" : False
      "cnLevelSelectionMode" : "ManualLevels"
      "cnMinLevelValF" : 900
      "cnMaxLevelValF" : 1100
      "cnLevelSpacingF" : mint
      "cnLevelFlags" : (/"LineAndLabel"/)
      "cnInfoLabelOn" : False
      "cnLinesOn" : True
      "cnLineColor" : 23
      "cnLineLabelsOn" : True
      "cnLineLabelInterval" : 1
      "cnLineLabelFont" : 21
      "cnLineLabelFontHeightF" : 0.01
      "cnLineLabelFontColor" : 23 ; "Background"
      "cnLineLabelBackgroundColor" : 1
      "cnLineLabelPlacementMode" : "Computed"
      "cnLineLabelDensityF" : mdens
      "cnLineThicknessF" : mlinethick
      "pmTickMarkDisplayMode" : "NoCreate"
      "vpXF"      : vpx        ; set viewport
      "vpYF"      : vpy
      "vpWidthF"  : vpwidth
      "vpHeightF" : vpheight
  end create
;
  cnid6 = create "contourplot" contourPlotClass wid
      "cnScalarFieldData":  zfield
      "cnFillOn" : False
      "cnLevelSelectionMode" : "ManualLevels"
      "cnMinLevelValF" : 402
      "cnMaxLevelValF" : 600
      "cnLevelSpacingF" : thint
      "cnLevelFlags" : (/"LineAndLabel"/)
      "cnInfoLabelOn" : False
      "cnLinesOn" : True
      "cnLineColor" : 2
;      "cnLineDashPattern" : 5
      "cnLineLabelsOn" : True
      "cnLineLabelInterval" : 1
      "cnLineLabelFont" : 21
      "cnLineLabelFontHeightF" : 0.01
      "cnLineLabelFontColor" : 2 ; "Foreground"
      "cnLineLabelBackgroundColor" : 1 
      "cnLineLabelPlacementMode" : "Computed"
      "cnLineLabelDensityF" : thdens
      "cnLineThicknessF" : thlinethick
      "pmTickMarkDisplayMode" : "NoCreate"
      "vpXF"      : vpx        ; set viewport
      "vpYF"      : vpy
      "vpWidthF"  : vpwidth
      "vpHeightF" : vpheight
  end create
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
  if (domainname .eq. "fim" .or. domainname .eq. "28" .or. domainname .eq. "129") then
    titlestring = interval + "h Conv Precip (in), MSLP (mb)"
  else
    titlestring = interval + "h Conv Precip (in), MSLP (mb), 1000-500 Thickness (dm)"
  end if
;
  load "$NCL_HOME/obj/add_titles.ncl"
  load "$NCL_HOME/obj/add_explabel.ncl"
;
; Draw map of 3-h Acc Convective Precip, MSLP
;
  draw(cnid)
  if (domainname .eq. "fim" .or. domainname .eq. "28" .or. domainname .eq. "129") then
    if (modl .ne. "MPAS") then
      draw(cnid5)
    end if
  else
    if (modl .ne. "MPAS") then
      draw(cnid5)
    end if
    draw(cnid6)
  end if
  draw(cnid4)
  draw(mpid)
;
  gsres = True
  if (domainname .eq. "236") then
    do i = 0,dimsizes(airport_id)-1
      gsres@gsMarkerIndex = 16
      gsres@gsMarkerColor = 20
      gsres@gsMarkerSizeF = 0.004
      gsn_polymarker(wid,mpid,airport_loc(i,1),airport_loc(i,0),gsres)
      gsres@gsMarkerIndex = 16
      gsres@gsMarkerColor = 1
      gsres@gsMarkerSizeF = 0.002
      gsn_polymarker(wid,mpid,airport_loc(i,1),airport_loc(i,0),gsres)
    end do
  end if
;
  frame(wid)
;
end if  ; domainname .ne. "130" .and. domainname .ne. "174" .and. domainname .ne. "242"
;
end
;
