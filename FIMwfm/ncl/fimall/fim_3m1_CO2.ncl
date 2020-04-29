;
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_LIB/ncarg/nclscripts/csm/contributed.ncl"
theos = systemfunc("uname")
modl = getenv("MODL")
;
load "$NCL_HOME/obj/cmap_diff3.ncl"
;
begin
;
airport_id = asciiread("$NCL_HOME/Airports_id.txt",-1,"string")
airport_loc = asciiread("$NCL_HOME/Airports_locs.txt",(/dimsizes(airport_id),2/),"float")
;
fname = asciiread("arw_file.txt",-1,"string")
grb_file = addfile(fname,"r")
;
levnam = grb_file->lv_ISBL3
do ii = 0, dimsizes(levnam)-1
  if (levnam(ii) .eq. 500) then
    level = ii
  end if
end do
;
domainname = asciiread("domain.txt",-1,"string")
;
load "$NCL_HOME/obj/fim_grib_names.ncl"
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
  tvar = grb_file->CFNLF_4_HYBL(0,:,:)
  tvar2 = grb_file->CFNLF_4_HYBL(2,:,:)
  tvar = tvar2 - tvar
  hgt = grb_file->GH_4_ISBL(level,:,:) / 10. ; use decameters
  hybhgt = grb_file->GH_4_HYBL(:,:,:) * 3.28084 ; convert to feet
  wstrid = 10
;
  uvar = grb_file->UW_4_HYBL(0,:,:) * 1.9429 ; convert to knots
  vvar = grb_file->VW_4_HYBL(0,:,:) * 1.9429 ; convert to knots
;
end if
;
if (domainname .eq. "129") then
;
  ny = filevardimsizes(grb_file,"lat_129")
  nx = filevardimsizes(grb_file,"lon_129")
;
  intm = grb_file->T_129_SFC@initial_time
  fctm = grb_file->T_129_SFC@forecast_time
  lat = grb_file->lat_129(:)
  lon = grb_file->lon_129(:)
  tvar = grb_file->CFNLF_129_HYBL(0,:,:)
  tvar2 = grb_file->CFNLF_129_HYBL(2,:,:)
  tvar = tvar2 - tvar
  hgt = grb_file->GH_129_ISBL(level,:,:) / 10. ; use decameters
  hybhgt = grb_file->GH_129_HYBL(:,:,:) * 3.28084 ; convert to feet
  wstrid = 20
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
  tvar = grb_file->CFNLF_236_HYBL(0,:,:)
  tvar2 = grb_file->CFNLF_236_HYBL(2,:,:)
  tvar = tvar2 - tvar
  hgt = grb_file->GH_236_ISBL(level,:,:) / 10. ; use decameters
  hybhgt = grb_file->GH_236_HYBL(:,:,:) * 3.28084 ; convert to feet
  wstrid = 3
;
  uvar = grb_file->UW_236_HYBL(0,:,:) * 1.9429 ; convert to knots
  vvar = grb_file->VW_236_HYBL(0,:,:) * 1.9429 ; convert to knots
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
  tvar = grb_file->CFNLF_GDS0_HYBL(0,:,:)
  tvar2 = grb_file->CFNLF_GDS0_HYBL(2,:,:)
  tvar = tvar2 - tvar
  hgt = grb_file->GH_GDS0_ISBL(level,:,:) / 10. ; use decameters
  hybhgt = grb_file->GH_GDS0_HYBL(:,:,:) * 3.28084 ; convert to feet
  wstrid = 5
  if (domainname .eq. "floating") then
    wstrid = 2
  end if
  if (domainname .eq. "w_pacific" .or. domainname .eq. "e_pacific") then
    wstrid = 7
  end if
;
  uvar = grb_file->UW_GDS0_HYBL(0,:,:) * 1.9429 ; convert to knots
  vvar = grb_file->VW_GDS0_HYBL(0,:,:) * 1.9429 ; convert to knots
;
end if
;
if (domainname .eq. "201") then
;
  fimsizes = filevardimsizes(grb_file,"gridlat_201")
  ny = fimsizes(0)
  nx = fimsizes(1)
;
  intm = grb_file->T_201_SFC@initial_time
  fctm = grb_file->T_201_SFC@forecast_time
  lat = grb_file->gridlat_201(:,:)
  lon = grb_file->gridlon_201(:,:)
  tvar = grb_file->CFNLF_201_HYBL(0,:,:)
  tvar2 = grb_file->CFNLF_201_HYBL(2,:,:)
  tvar = tvar2 - tvar
  hgt = grb_file->GH_201_ISBL(level,:,:) / 10. ; use decameters
  hybhgt = grb_file->GH_201_HYBL(:,:,:) * 3.28084 ; convert to feet
  wstrid = 10
;
  uvar = grb_file->UW_201_HYBL(0,:,:) * 1.9429 ; convert to knots
  vvar = grb_file->VW_201_HYBL(0,:,:) * 1.9429 ; convert to knots
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
  tvar = grb_file->CFNLF_28_HYBL(0,:,:)
  tvar2 = grb_file->CFNLF_28_HYBL(2,:,:)
  tvar = tvar2 - tvar
  hgt = grb_file->GH_28_ISBL(level,:,:) / 10. ; use decameters
  hybhgt = grb_file->GH_28_HYBL(:,:,:) * 3.28084 ; convert to feet
  wstrid = 2
;
  uvar = grb_file->UW_28_HYBL(0,:,:) * 1.9429 ; convert to knots
  vvar = grb_file->VW_28_HYBL(0,:,:) * 1.9429 ; convert to knots
  uvar = -uvar
  vvar = -vvar
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
  tvar = grb_file->CFNLF_244_HYBL(0,:,:)
  tvar2 = grb_file->CFNLF_244_HYBL(2,:,:)
  tvar = tvar2 - tvar
  hgt = grb_file->GH_244_ISBL(level,:,:) / 10. ; use decameters
  hybhgt = grb_file->GH_244_HYBL(:,:,:) * 3.28084 ; convert to feet
  wstrid = 7
;
  uvar = grb_file->UW_244_HYBL(0,:,:) * 1.9429 ; convert to knots
  vvar = grb_file->VW_244_HYBL(0,:,:) * 1.9429 ; convert to knots
;
end if
;
;print("min tvar = "+min(tvar))
;print("max tvar = "+max(tvar))
;
;tvar = tvar * 10.
;do ii = 0,nx-1
;  do jj = 0,ny-1
;    if (.not.ismissing(tvar(jj,ii))) then
;      if (tvar(jj,ii) .lt. 2.) then
;        tvar(jj,ii) = 1.
;      end if
;    end if
;  end do
;end do
;tvar = log10(tvar)
;
;print("min tvar after log10 = "+min(tvar))
;print("max tvar after log10 = "+max(tvar))
;
load "$NCL_HOME/obj/get_date_str_FIM.ncl"
;
print(date_str) 
print("Hyb level 2 CO2 - Sfc CO2")
print(min(tvar))
print(max(tvar))
;print("500 mb hgt")
;print(min(hgt))
;print(max(hgt))
;
if (max(tvar) .eq. 0) then  ; constant zero field
  tvar(0,0) = 0.11 ; set the ll corner point so it makes a plot anyway
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
           "wkMetaName" : "./3m1_CO2.ncgm"
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
           "wkPSFileName" : "./3m1_CO2.ps"
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
if (domainname .eq. "floating") then
  setvalues mpid
    "mpDataSetName" : "Earth..4"
  end setvalues
end if
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
  hfield = create "ScalarField" scalarFieldClass appid
      "sfDataArray" : hgt
      "sfXCStartV" : 0
      "sfYCStartV" : ydim
      "sfXCEndV" : xdim
      "sfYCEndV" : 0
      "sfMissingValueV" : hgt@_FillValue
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
      "vfXCStride" : wstrid       ; i.e., plot a vector at every (n+1)th point.
      "vfYCStride" : wstrid
  end create
;
end if
;
if (domainname .eq. "236" .or. domainname .eq. "africa" .or. \
    domainname .eq. "201" .or. domainname .eq. "w_pacific" .or. \
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
  vfield = create "VectorField" vectorFieldClass appid
      "vfUDataArray" : uvar
      "vfVDataArray" : vvar
      "vfXCStartV" : 0
      "vfYCStartV" : 0
      "vfXCEndV" : xdim
      "vfYCEndV" : ydim
      "vfMissingVValueV" : vvar@_FillValue
      "vfMissingUValueV" : uvar@_FillValue
      "vfXCStride" : wstrid        ; i.e., plot a vector at every (n+1)th point.
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
      "sfDataArray" : hgt
      "sfXCStartV" : xdim
      "sfYCStartV" : ydim
      "sfXCEndV" : 0
      "sfYCEndV" : 0
      "sfMissingValueV" : hgt@_FillValue
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
    "vcWindBarbTickSpacingF" : 0.1875
    "pmTickMarkDisplayMode" : "NoCreate"
    "vpXF"      : vpx        ; set viewport
    "vpYF"      : vpy       
    "vpWidthF"  : vpwidth
    "vpHeightF" : vpheight
end create
;
setvalues vcid
  "vcMinMagnitudeF" : 0
end setvalues
;
; Create a ContourPlot object.
;
; hybl 0 CO2
;
load "$NCL_HOME/obj/diff_-8-8by1.ncl"
;
if (domainname .eq. "236" .or. domainname .eq. "africa" .or. \
    domainname .eq. "244" .or. domainname .eq. "w_pacific" .or. \
    domainname .eq. "europe" .or. domainname .eq. "28" .or. \
    domainname .eq. "floating" .or. domainname .eq. "e_pacific" .or. \
    domainname .eq. "83" .or. domainname .eq. "201") then
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
    "cnMinLevelValF" : 486
    "cnMaxLevelValF" : 610
    "cnLevelSpacingF" : 6
    "cnLinesOn" : True
    "cnMonoLineColor" : True
    "cnLineColor" : 12
    "cnLineThicknessF" : 1.25
    "cnLineLabelsOn" : True
    "cnLineLabelFontHeightF" : 0.01
    "cnLineLabelPlacementMode" : "Computed"
;    "cnLineLabelDensityF" : 0.2
    "cnInfoLabelOn" : False
    "cnLineLabelInterval" : 1
    "cnLineLabelFont" : 21
    "cnLineLabelFontColor" : 12
    "cnLineLabelBackgroundColor" : -1
    "pmTickMarkDisplayMode" : "NoCreate"
    "cnFillOn" : False
    "cnGridBoundPerimOn" : False
    "vpXF"      : vpx        ; set viewport
    "vpYF"      : vpy
    "vpWidthF"  : vpwidth
    "vpHeightF" : vpheight
end create
;
if (domainname .eq. "201") then
;
  setvalues cnid2
    "cnLineLabelDensityF" : 1.0
  end setvalues
;
end if
;
titlestring = "Vertical diff CO2 conc, level 3 - level 1 (ppm)"
load "$NCL_HOME/obj/add_titles.ncl"
load "$NCL_HOME/obj/add_explabel.ncl"
;
; Draw map
;
draw(cnid)
draw(mpid)
;draw(cnid2)
;draw(vcid)
;
gsres = True
if (domainname .eq. "236") then
  do i = 0,dimsizes(airport_id)-1
    gsres@gsMarkerIndex = 16
    gsres@gsMarkerColor = 36
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

