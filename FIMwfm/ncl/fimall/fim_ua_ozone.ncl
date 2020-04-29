;
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_LIB/ncarg/nclscripts/csm/contributed.ncl"
theos = systemfunc("uname")
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
level = 50  ; hybrid level
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
  tvar = grb_file->O3MR_4_HYBL(level,:,:) * 1000.  ; already scaled by 1000, scale up to 1e6
  hgt = grb_file->GH_4_HYBL(level,:,:) / 10. ; use decameters
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
  tvar = grb_file->O3MR_236_HYBL(level,:,:) * 1000.  ; already scaled by 1000, scale up to 1e6
  hgt = grb_file->GH_236_HYBL(level,:,:) / 10. ; use decameters
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
  tvar = grb_file->O3MR_GDS0_HYBL(level,:,:) * 1000.  ; already scaled by 1000, scale up to 1e6
  hgt = grb_file->GH_GDS0_HYBL(level,:,:) / 10. ; use decameters
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
  tvar = grb_file->O3MR_219_HYBL(level,:,:) * 1000.  ; already scaled by 1000, scale up to 1e6
  hgt = grb_file->GH_219_HYBL(level,:,:) / 10. ; use decameters
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
  tvar = grb_file->O3MR_28_HYBL(level,:,:) * 1000.  ; already scaled by 1000, scale up to 1e6
  hgt = grb_file->GH_28_HYBL(level,:,:) / 10. ; use decameters
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
  tvar = grb_file->O3MR_244_HYBL(level,:,:) * 1000.  ; already scaled by 1000, scale up to 1e6
  hgt = grb_file->GH_244_HYBL(level,:,:) / 10. ; use decameters
;
end if
;
load "$NCL_HOME/obj/get_date_str_FIM.ncl"
;
print(date_str) 
print("ozone mixing ratio")
print(min(tvar))
print(max(tvar))
print("hybl 50 hgt")
print(min(hgt))
print(max(hgt))
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
           "wkMetaName" : "./ua_ozone.ncgm"
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
           "wkPSFileName" : "./ua_ozone.ps"
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
; Create a ScalarField object.
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
      "sfDataArray" : hgt
      "sfXCStartV" : 0
      "sfYCStartV" : ydim
      "sfXCEndV" : xdim
      "sfYCEndV" : 0
      "sfMissingValueV" : hgt@_FillValue
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
      "sfDataArray" : tvar
      "sfXCStartV" : 0
      "sfYCStartV" : 0
      "sfXCEndV" : xdim
      "sfYCEndV" : ydim
      "sfMissingValueV" : tvar@_FillValue
  end create
;
  hfield = create "ScalarField" scalarFieldClass appid
      "sfDataArray" : hgt
      "sfXCStartV" : 0
      "sfYCStartV" : 0
      "sfXCEndV" : xdim
      "sfYCEndV" : ydim
      "sfMissingValueV" : hgt@_FillValue
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
      "sfDataArray" : hgt
      "sfXCStartV" : xdim
      "sfYCStartV" : ydim
      "sfXCEndV" : 0
      "sfYCEndV" : 0
      "sfMissingValueV" : hgt@_FillValue
  end create
;
end if
;
; Create a ContourPlot object.
;
; Precipitable Water
;
load "$NCL_HOME/obj/ozone.ncl"
;
if (domainname .eq. "236" .or. domainname .eq. "africa" .or. \
    domainname .eq. "244" .or. domainname .eq. "w_pacific" .or. \
    domainname .eq. "europe" .or. domainname .eq. "28" .or. \
    domainname .eq. "floating" .or. domainname .eq. "e_pacific") then
;
  setvalues cnid
    "cnRasterModeOn" : False
  end setvalues
;
end if
;
cnid2 = create "contourplot" contourPlotClass wid
    "cnScalarFieldData":  hfield
    "cnLevelSelectionMode" : "ManualLevels"
    "cnMinLevelValF" : 2000
    "cnMaxLevelValF" : 3000
    "cnLevelSpacingF" : 24
    "cnLinesOn" : True
    "cnMonoLineColor" : True
    "cnLineColor" : 1
    "cnLineThicknessF" : 1.25
    "cnLineLabelsOn" : True
    "cnLineLabelFontHeightF" : 0.01
    "cnLineLabelPlacementMode" : "Computed"
    "cnLineLabelDensityF" : 0.2
    "cnInfoLabelOn" : False
    "cnLineLabelInterval" : 1
    "cnLineLabelFont" : 21
    "cnLineLabelFontColor" : "Foreground"
    "cnLineLabelBackgroundColor" : "Background"
    "pmTickMarkDisplayMode" : "NoCreate"
    "cnFillOn" : False
    "cnGridBoundPerimOn" : False
    "vpXF"      : vpx        ; set viewport
    "vpYF"      : vpy
    "vpWidthF"  : vpwidth
    "vpHeightF" : vpheight
end create
;
if (domainname .eq. "219") then
;
  setvalues cnid2
    "cnLineLabelDensityF" : 1.0
  end setvalues
;
end if
;
titlestring = "Ozone Mixing Ratio (kg/kg)*1e6"
load "$NCL_HOME/obj/add_titles.ncl"
load "$NCL_HOME/obj/add_explabel.ncl"
;
; Draw map
;
draw(cnid)
draw(mpid)
;draw(cnid2)
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
end

