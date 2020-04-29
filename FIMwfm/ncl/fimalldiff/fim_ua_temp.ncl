;
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_LIB/ncarg/nclscripts/csm/contributed.ncl"
theos = systemfunc("uname")
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
if (domainname .eq. "africa" .or. domainname .eq. "europe") then
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
if (domainname .eq. "w_pacific") then
;
  ny = filevardimsizes(grb_file,"g0_lat_1")
  nx = filevardimsizes(grb_file,"g0_lon_2")
;
  intm = grb_file->T_GDS0_SFC@initial_time
  fctm = grb_file->T_GDS0_SFC@forecast_time
  lat = grb_file->g0_lat_1(:)
  lon = grb_file->g0_lon_2(:)
;
end if
;
lev = (/ 33, 9 /)
levnam = (/ 850, 250 /);
zmin = (/ 1100, 9000 /)
zmax = (/ 1800, 12000 /)
space = (/ 30, 300 /)
tmin = (/ -30, -70 /)
tmax = (/ 30, -10 /)
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
print(date_str)
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
           "wkMetaName" : "./ua_temp.ncgm"
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
           "wkPSFileName" : "./ua_temp.ps"
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
; Increase workspace memory
setvalues NhlGetWorkspaceObjectId()
   "wsMaximumSize": 247108864
end setvalues
;
load "$NCL_HOME/obj/fim_maps_blk.ncl"
;
do i = 0,1
;
  if (domainname .eq. "fim") then
;
    tvar = grb_file->T_4_ISBL(lev(i),:,:) - 273.15
    tvar2 = grb_file2->T_4_ISBL(lev(i),:,:) - 273.15
;
  end if
;
  if (domainname .eq. "236") then
;
    tvar = grb_file->T_236_ISBL(lev(i),:,:) - 273.15
    tvar2 = grb_file2->T_236_ISBL(lev(i),:,:) - 273.15
;
  end if
;
  if (domainname .eq. "africa" .or. domainname .eq. "europe") then
;
    tvar = grb_file->T_GDS0_ISBL(lev(i),:,:) - 273.15
    tvar2 = grb_file2->T_GDS0_ISBL(lev(i),:,:) - 273.15
;
  end if
;
  if (domainname .eq. "219") then
;
    tvar = grb_file->T_219_ISBL(lev(i),:,:) - 273.15
    tvar2 = grb_file2->T_219_ISBL(lev(i),:,:) - 273.15
;
  end if
;
  if (domainname .eq. "28") then
;
    tvar = grb_file->T_28_ISBL(lev(i),:,:) - 273.15
    tvar2 = grb_file2->T_28_ISBL(lev(i),:,:) - 273.15
;
  end if
;
  if (domainname .eq. "244") then
;
    tvar = grb_file->T_244_ISBL(lev(i),:,:) - 273.15
    tvar2 = grb_file2->T_244_ISBL(lev(i),:,:) - 273.15
;
  end if
;
  if (domainname .eq. "w_pacific") then
;
    tvar = grb_file->T_GDS0_ISBL(lev(i),:,:) - 273.15
    tvar2 = grb_file2->T_GDS0_ISBL(lev(i),:,:) - 273.15
;
  end if
;
  tvar = tvar - tvar2
;
  print("domain " + domainname)
  print(levnam(i))
  print("temp diff")
  print(min(tvar))
  print(max(tvar))
;
; Create ScalarField and VectorField objects.
;
if (domainname .eq. "fim" .or. domainname .eq. "244") then
;
; Temperature
;
    tfield = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : tvar
        "sfXCStartV" : 0
        "sfYCStartV" : ydim
        "sfXCEndV" : xdim
        "sfYCEndV" : 0
        "sfMissingValueV" : -999.0
    end create
;
  end if
;
  if (domainname .eq. "236" .or. domainname .eq. "africa" .or. \
      domainname .eq. "w_pacific" .or. domainname .eq. "europe") then
;
; Temperature
;
    tfield = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : tvar
        "sfXCStartV" : 0
        "sfYCStartV" : 0
        "sfXCEndV" : xdim
        "sfYCEndV" : ydim
        "sfMissingValueV" : -999.0
    end create
;
  end if
;
  if (domainname .eq. "219") then
;
; Temperature
;
    tfield = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : tvar
        "sfXCStartV" : 0
        "sfYCStartV" : 0
        "sfXCEndV" : xdim
        "sfYCEndV" : ydim
        "sfMissingValueV" : -999.0
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
      "sfMissingValueV" : -999.0
  end create
;
end if
;
; Create ContourPlot objects.
;
; Temperature
;
  if (i .eq. 0) then
;
    load "$NCL_HOME/obj/diff_-8-8by1.ncl"
;
  else
;
    load "$NCL_HOME/obj/diff_-8-8by1.ncl"
;
  end if
;
  if (domainname .eq. "236" .or. domainname .eq. "africa" .or. \
      domainname .eq. "244" .or. domainname .eq. "w_pacific" .or. \
      domainname .eq. "europe" .or. domainname .eq. "28") then
;
    setvalues cnid
;      "cnRasterModeOn" : False
      "cnRasterSmoothingOn" : True
    end setvalues
;
  end if
;
titleid = create "TextItems" textItemClass wid
  "txFuncCode" : "~"
  "txString" : levnam(i) + "mb Temperature (C), Height (m), Wind (kt)"
  "txFont" : 22
  "txFontHeightF" : 0.02
  "txFontColor" : 1
end create
;
txid = create "TextItems" textItemClass wid
  "txFuncCode" : "~"
  "txString" : modl+" "+intm+" "+fctm+" hr fcst"
  "txFont" : 22
  "txFontHeightF" : 0.015
  "txFontColor" : 1
end create
;
txid2 = create "TextItems" textItemClass wid
  "txFuncCode" : "~"
  "txString" : date_str
  "txFont" : 22
  "txFontHeightF" : 0.015
  "txFontColor" : 1
end create
;
iannotitle = NhlAddAnnotation(cnid,titleid)
setvalues iannotitle
  "amZone":           10
  "amSide":           "top"
  "amJust":           "bottomright"
  "amParallelPosF":    1.0
  "amOrthogonalPosF":  0.0
end setvalues
;
if (vpheight .ge. vpwidth) then
  opos = 0.032
else
  opos = 0.1 - ((vpheight/vpwidth)*0.07)
end if
;
ianno = NhlAddAnnotation(cnid,txid)
setvalues ianno
  "amZone":           10
  "amSide":           "top"
  "amJust":           "bottomleft"
  "amParallelPosF":    0.0
  "amOrthogonalPosF":  opos
end setvalues
;
ianno2 = NhlAddAnnotation(cnid,txid2)
setvalues ianno2
  "amZone":           10
  "amSide":           "top"
  "amJust":           "bottomright"
  "amParallelPosF":    1.0
  "amOrthogonalPosF":  opos
end setvalues
;
; Draw map of <levnam>mb Temp, Height, Winds
;
  draw(cnid)
  draw(mpid)
  frame(wid)
;
end do
;
end
;
