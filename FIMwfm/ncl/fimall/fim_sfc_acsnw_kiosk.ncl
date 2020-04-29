;
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_LIB/ncarg/nclscripts/csm/contributed.ncl"
modl = getenv("MODL")
interval = getenv("FCST_INTERVAL")
interval = 6
;
load "$NCL_HOME/obj/cmap1.ncl"
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
if (domainname .ne. "fim") then
  exit  ; only create these for the global domain
end if
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
  if (fctm .gt. 0) then
    if (isfilevar(grb_file,WEASDACC)) then
      zz = WEASDACC
      tvar = grb_file->$zz$(:,:) * 0.03937; convert to inches  
    else
      zz = WEASDACC + interval + "h"
      if (isfilevar(grb_file,zz)) then
        tvar = grb_file->$zz$(:,:) * 0.03937; convert to inches  
      else
        zz = WEASDACC + fctm + "h" ; use total snowfall
        tvar = grb_file->$zz$(:,:) * 0.03937; convert to inches  
        back6fname = asciiread("back6_file.txt",-1,"string")
        back6_grb_file = addfile(back6fname,"r")
        back6fctm = back6_grb_file->$TSFC$@forecast_time
        back6fctmunits = back6_grb_file->$TSFC$@forecast_time_units
        if (back6fctmunits .ne. "hours") then
          back6fctm = back6fctm * tointeger(back6fctmunits)
        end if
        back6fcst1 = back6fctm
        zz = WEASDACC + back6fcst1 + "h"
        back6tvar = back6_grb_file->$zz$(:,:) * 0.03937 ; convert to inches
        tvar = tvar - back6tvar
      end if
    end if
  else
    exit
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
  if (isfilevar(grb_file,WEASDACC)) then
    zz = WEASDACC
  else
    zz = WEASDACC + interval + "h"
  end if
  tvar = grb_file->$zz$(:,:) * 0.03937; convert to inches
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
  if (fctm .gt. 0) then
    if (isfilevar(grb_file,WEASDACC)) then
      zz = WEASDACC
      tvar = grb_file->$zz$(:,:) * 0.03937; convert to inches
    else
      zz = WEASDACC + interval + "h"
      if (isfilevar(grb_file,zz)) then
        tvar = grb_file->$zz$(:,:) * 0.03937; convert to inches
      else
        zz = WEASDACC + fctm + "h" ; use total snowfall
        tvar = grb_file->$zz$(:,:) * 0.03937; convert to inches
        back6fname = asciiread("back6_file.txt",-1,"string")
        back6_grb_file = addfile(back6fname,"r")
        back6fctm = back6_grb_file->$TSFC$@forecast_time
        back6fctmunits = back6_grb_file->$TSFC$@forecast_time_units
        if (back6fctmunits .ne. "hours") then
          back6fctm = back6fctm * tointeger(back6fctmunits)
        end if
        back6fcst1 = back6fctm
        zz = WEASDACC + back6fcst1 + "h"
        back6tvar = back6_grb_file->$zz$(:,:) * 0.03937 ; convert to inches
        tvar = tvar - back6tvar
      end if
    end if
  else
    exit
  end if
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
  if (fctm .gt. 0) then
    if (isfilevar(grb_file,WEASDACC)) then
      zz = WEASDACC
      tvar = grb_file->$zz$(:,:) * 0.03937; convert to inches
    else
      zz = WEASDACC + interval + "h"
      if (isfilevar(grb_file,zz)) then
        tvar = grb_file->$zz$(:,:) * 0.03937; convert to inches
      else
        zz = WEASDACC + fctm + "h" ; use total snowfall
        tvar = grb_file->$zz$(:,:) * 0.03937; convert to inches
        back6fname = asciiread("back6_file.txt",-1,"string")
        back6_grb_file = addfile(back6fname,"r")
        back6fctm = back6_grb_file->$TSFC$@forecast_time
        back6fctm = back6_grb_file->$TSFC$@forecast_time
        back6fctmunits = back6_grb_file->$TSFC$@forecast_time_units
        if (back6fctmunits .ne. "hours") then
          back6fctm = back6fctm * tointeger(back6fctmunits)
        end if
        back6fcst1 = back6fctm
        zz = WEASDACC + back6fcst1 + "h"
        back6tvar = back6_grb_file->$zz$(:,:) * 0.03937 ; convert to inches
        tvar = tvar - back6tvar
      end if
    end if
  else
    exit
  end if
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
  if (fctm .gt. 0) then
    if (isfilevar(grb_file,WEASDACC)) then
      zz = WEASDACC
      tvar = grb_file->$zz$(::-1,:) * 0.03937; convert to inches
    else
      zz = WEASDACC + interval + "h"
      if (isfilevar(grb_file,zz)) then
        tvar = grb_file->$zz$(::-1,:) * 0.03937; convert to inches
      else
        zz = WEASDACC + fctm + "h" ; use total snowfall
        tvar = grb_file->$zz$(::-1,:) * 0.03937; convert to inches
        back6fname = asciiread("back6_file.txt",-1,"string")
        back6_grb_file = addfile(back6fname,"r")
        back6fctm = back6_grb_file->$TSFC$@forecast_time
        back6fctmunits = back6_grb_file->$TSFC$@forecast_time_units
        if (back6fctmunits .ne. "hours") then
          back6fctm = back6fctm * tointeger(back6fctmunits)
        end if
        back6fcst1 = back6fctm
        zz = WEASDACC + back6fcst1 + "h"
        back6tvar = back6_grb_file->$zz$(::-1,:) * 0.03937 ; convert to inches
        tvar = tvar - back6tvar
      end if
    end if
  else
    exit
  end if
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
  if (fctm .gt. 0) then
    if (isfilevar(grb_file,WEASDACC)) then
      zz = WEASDACC
      tvar = grb_file->$zz$(:,:) * 0.03937; convert to inches
    else
      zz = WEASDACC + interval + "h"
      if (isfilevar(grb_file,zz)) then
        tvar = grb_file->$zz$(:,:) * 0.03937; convert to inches
      else
        zz = WEASDACC + fctm + "h" ; use total snowfall
        tvar = grb_file->$zz$(:,:) * 0.03937; convert to inches
        back6fname = asciiread("back6_file.txt",-1,"string")
        back6_grb_file = addfile(back6fname,"r")
        back6fctm = back6_grb_file->$TSFC$@forecast_time
        back6fctmunits = back6_grb_file->$TSFC$@forecast_time_units
        if (back6fctmunits .ne. "hours") then
          back6fctm = back6fctm * tointeger(back6fctmunits)
        end if
        back6fcst1 = back6fctm
        zz = WEASDACC + back6fcst1 + "h"
        back6tvar = back6_grb_file->$zz$(:,:) * 0.03937 ; convert to inches
        tvar = tvar - back6tvar
      end if
    end if
  else
    exit
  end if
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
  if (fctm .gt. 0) then
    if (isfilevar(grb_file,WEASDACC)) then
      zz = WEASDACC
      tvar = grb_file->$zz$(:,:) * 0.03937; convert to inches
    else
      zz = WEASDACC + interval + "h"
      if (isfilevar(grb_file,zz)) then
        tvar = grb_file->$zz$(:,:) * 0.03937; convert to inches
      else
        zz = WEASDACC + fctm + "h" ; use total snowfall
        tvar = grb_file->$zz$(:,:) * 0.03937; convert to inches
        back6fname = asciiread("back6_file.txt",-1,"string")
        back6_grb_file = addfile(back6fname,"r")
        back6fctm = back6_grb_file->$TSFC$@forecast_time
        back6fctmunits = back6_grb_file->$TSFC$@forecast_time_units
        if (back6fctmunits .ne. "hours") then
          back6fctm = back6fctm * tointeger(back6fctmunits)
        end if
        back6fcst1 = back6fctm
        zz = WEASDACC + back6fcst1 + "h"
        back6tvar = back6_grb_file->$zz$(:,:) * 0.03937 ; convert to inches
        tvar = tvar - back6tvar
      end if
    end if
  else
    exit
  end if
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
  if (fctm .gt. 0) then
    if (isfilevar(grb_file,WEASDACC)) then
      zz = WEASDACC
      tvar = grb_file->$zz$(:,:) * 0.03937; convert to inches
    else
      zz = WEASDACC + interval + "h"
      if (isfilevar(grb_file,zz)) then
        tvar = grb_file->$zz$(:,:) * 0.03937; convert to inches
      else
        zz = WEASDACC + fctm + "h" ; use total snowfall
        tvar = grb_file->$zz$(:,:) * 0.03937; convert to inches
        back6fname = asciiread("back6_file.txt",-1,"string")
        back6_grb_file = addfile(back6fname,"r")
        back6fctm = back6_grb_file->$TSFC$@forecast_time
        back6fctmunits = back6_grb_file->$TSFC$@forecast_time_units
        if (back6fctmunits .ne. "hours") then
          back6fctm = back6fctm * tointeger(back6fctmunits)
        end if
        back6fcst1 = back6fctm
        zz = WEASDACC + back6fcst1 + "h"
        back6tvar = back6_grb_file->$zz$(:,:) * 0.03937 ; convert to inches
        tvar = tvar - back6tvar
      end if
    end if
  else
    exit
  end if
;
end if
;
tvar = tvar * 10.  ; scaled for snow
;
load "$NCL_HOME/obj/get_date_str_FIM.ncl"
;
print(date_str)
print("acsnw")
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
wks_type = "png"
wks_type@wkWidth = 1600
wks_type@wkHeight = 1600
wid = gsn_open_wks(wks_type,"sfc_acsnw_kiosk")
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
      "mpCountyLineColor" : "background"
    end setvalues
  end if
;
  if (max(tvar(lly:ury,llx:urx)) .le. 0.01) then  ; set corner pt to avoid blank plot
    tvar(lly,llx) = 0.1
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
  load "$NCL_HOME/obj/acsnw.ncl"
;
  setvalues cnid
    "cnRasterModeOn" : False
  end setvalues
;
  titlestring = interval + "hr Accumulated Snow/Sleet Using 10:1 Ratio (in)"
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
  if (max(tvar(lly:ury,llx:urx)) .le. 0.01) then  ; set corner pt to avoid blank plot
    tvar(lly,llx) = 0.1
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
  load "$NCL_HOME/obj/weasd.ncl"
;
  titlestring = interval + "hr Accumulated Snow/Sleet Using 10:1 Ratio (in)"
  load "$NCL_HOME/obj/add_titles.ncl"
  load "$NCL_HOME/obj/add_explabel.ncl"
;
; Draw map
;
  draw(cnid)
  draw(mpid)
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
  end if
;
  load "$NCL_HOME/obj/acsnw.ncl"
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
  titlestring = interval + "hr Accumulated Snow/Sleet Using 10:1 Ratio (in)"
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
end if  ; domainname .ne. "130" .and. domainname .ne. "174"
;
end
