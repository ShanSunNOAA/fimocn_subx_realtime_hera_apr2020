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
; Get the 850mb level index
;
levnam = grb_file->lv_ISBL3
do ii = 0, dimsizes(levnam)-1
  if (levnam(ii) .eq. 850) then
    level = ii
  end if
end do
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
  uvar = grb_file->UW_4_ISBL(level,:,:) * 1.9429 ; convert to knots
  vvar = grb_file->VW_4_ISBL(level,:,:) * 1.9429 ; convert to knots
  tvar = grb_file->RH_4_ISBL(level,:,:)
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
  uvar = grb_file->UW_83_ISBL(level,:,:) * 1.9429 ; convert to knots
  vvar = grb_file->VW_83_ISBL(level,:,:) * 1.9429 ; convert to knots
  tvar = grb_file->RH_83_ISBL(level,:,:)
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
  uvar = grb_file->UW_130_ISBL(level,:,:) * 1.9429 ; convert to knots
  vvar = grb_file->VW_130_ISBL(level,:,:) * 1.9429 ; convert to knots
  tvar = grb_file->RH_130_ISBL(level,:,:)
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
  uvar = grb_file->UW_236_ISBL(level,:,:) * 1.9429 ; convert to knots
  vvar = grb_file->VW_236_ISBL(level,:,:) * 1.9429 ; convert to knots
  tvar = grb_file->RH_236_ISBL(level,:,:)
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
  uvar = grb_file->UW_GDS0_ISBL(level,:,:) * 1.9429 ; convert to knots
  vvar = grb_file->VW_GDS0_ISBL(level,:,:) * 1.9429 ; convert to knots
  tvar = grb_file->RH_GDS0_ISBL(level,:,:)
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
  uvar = grb_file->UW_219_ISBL(level,:,:) * 1.9429 ; convert to knots
  vvar = grb_file->VW_219_ISBL(level,:,:) * 1.9429 ; convert to knots
  tvar = grb_file->RH_219_ISBL(level,:,:)
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
  uvar = grb_file->UW_28_ISBL(level,:,:) * 1.9429 ; convert to knots
  vvar = grb_file->VW_28_ISBL(level,:,:) * 1.9429 ; convert to knots
  uvar = -uvar
  vvar = -vvar
  tvar = grb_file->RH_28_ISBL(level,:,:)
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
  uvar = grb_file->UW_244_ISBL(level,:,:) * 1.9429 ; convert to knots
  vvar = grb_file->VW_244_ISBL(level,:,:) * 1.9429 ; convert to knots
  tvar = grb_file->RH_244_ISBL(level,:,:)
;
end if
;
load "$NCL_HOME/obj/get_date_str_FIM.ncl"
;
print(date_str)
print("RH 850")
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
           "wkMetaName" : "./ua_rh8.ncgm"
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
           "wkPSFileName" : "./ua_rh8.ps"
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
  load "$NCL_HOME/obj/rh.ncl"
;
  setvalues cnid
    "cnRasterModeOn" : False
  end setvalues
;
  titlestring = "850mb Rel Hum (%)"
  load "$NCL_HOME/obj/add_titles.ncl"
  load "$NCL_HOME/obj/add_explabel.ncl"
;
; Draw map
;
  draw(cnid)
  draw(mpid)
;
  gsres = True
  if (iq .eq. 0) then
    do i = 0,dimsizes(airport_id)-1
      gsres@gsMarkerIndex = 16
      gsres@gsMarkerColor = 29
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
      gsres@gsMarkerColor = 29
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
else  ; domain is not 130
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
  end if
;
; Create a ContourPlot object.
;
  load "$NCL_HOME/obj/rh.ncl"
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
  titlestring = "850mb Rel Hum (%)"
  load "$NCL_HOME/obj/add_titles.ncl"
  load "$NCL_HOME/obj/add_explabel.ncl"
;
; Draw map
;
  draw(cnid)
  draw(mpid)
;
  gsres = True
  if (domainname .eq. "236") then
    do i = 0,dimsizes(airport_id)-1
      gsres@gsMarkerIndex = 16
      gsres@gsMarkerColor = 29
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
end if  ; domainname .eq. "130"
;
end
;
