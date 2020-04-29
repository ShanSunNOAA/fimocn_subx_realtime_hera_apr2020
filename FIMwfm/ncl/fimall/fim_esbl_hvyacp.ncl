;
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_LIB/ncarg/nclscripts/csm/contributed.ncl"
modl = getenv("MODL")
interval = getenv("FCST_INTERVAL")
interval = 6
;
load "$NCL_HOME/obj/cmap1_esbl.ncl"
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
load "$NCL_HOME/obj/fim_grib_names.ncl"
;
load "$NCL_HOME/obj/esblfilecheck.ncl"
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
      tvar = grb_file->ACPCP_GDS0_SFC(::-1,:) + grb_file->NCPCP_GDS0_SFC(::-1,:)
      tvar = tvar * 0.03937 ; convert to inches
      tvar = lonFlip(tvar)
    else
      exit
    end if
;
  else
;
    ny = filevardimsizes(grb_file,LAT)
    nx = filevardimsizes(grb_file,LON)
;
    lat = grb_file->$LAT$(:)
    lon = grb_file->$LON$(:)
    if ( fctm .gt. 0 ) then
      if (isfilevar(grb_file,TPACC)) then
        zz = TPACC
      else
        zz = TPACC + interval + "h"
      end if
      tvar = grb_file->$zz$(:,:) * 0.03937 ; convert to inches
    else
      exit
    end if
;
    tvar = where (.not. ismissing(tvar) .and. tvar .lt. 0, 0, tvar)
    tvarback1 = tvar * 0.
    tvarback2 = tvar * 0.
    if (fnameback1test .eq. 1) then
      if ( isfilevar(grb_file_back1,TPACC) ) then
        zz = TPACC
      else
        zz = TPACC + interval + "h"
      end if
      tvarback1 = grb_file_back1->$zz$(:,:) * 0.03937 ; convert to inches
      tvarback1 = where (.not. ismissing(tvarback1) .and. tvarback1 .lt. 0, 0, tvarback1)
    end if
    if (fnameback2test .eq. 1) then
      if ( isfilevar(grb_file_back2,TPACC) ) then
        zz = TPACC
      else
        zz = TPACC + interval + "h"
      end if
      tvarback2 = grb_file_back2->$zz$(:,:) * 0.03937 ; convert to inches
      tvarback2 = where (.not. ismissing(tvarback2) .and. tvarback2 .lt. 0, 0, tvarback2)
    end if
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
    if (isfilevar(grb_file,TPACC)) then
      zz = TPACC
    else
      zz = TPACC + interval + "h"
    end if
    tvar = grb_file->$zz$(:,:) * 0.03937 ; convert to inches
  else
    exit
  end if
;
  tvar = where (.not. ismissing(tvar) .and. tvar .lt. 0, 0, tvar)
  tvarback1 = tvar * 0.
  tvarback2 = tvar * 0.
  if (fnameback1test .eq. 1) then
    if ( isfilevar(grb_file_back1,TPACC) ) then
      zz = TPACC
    else
      zz = TPACC + interval + "h"
    end if
    tvarback1 = grb_file_back1->$zz$(:,:) * 0.03937 ; convert to inches
    tvarback1 = where (.not. ismissing(tvarback1) .and. tvarback1 .lt. 0, 0, tvarback1)
  end if
  if (fnameback2test .eq. 1) then
    if ( isfilevar(grb_file_back2,TPACC) ) then
      zz = TPACC
    else
      zz = TPACC + interval + "h"
    end if
    tvarback2 = grb_file_back2->$zz$(:,:) * 0.03937 ; convert to inches
    tvarback2 = where (.not. ismissing(tvarback2) .and. tvarback2 .lt. 0, 0, tvarback2)
  end if
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
    if (isfilevar(grb_file,TPACC)) then
      zz = TPACC
    else
      zz = TPACC + interval + "h"
    end if
    tvar = grb_file->$zz$(:,:) * 0.03937 ; convert to inches
  else
    exit
  end if
;
  tvar = where (.not. ismissing(tvar) .and. tvar .lt. 0, 0, tvar)
  tvarback1 = tvar * 0.
  tvarback2 = tvar * 0.
  if (fnameback1test .eq. 1) then
    if ( isfilevar(grb_file_back1,TPACC) ) then
      zz = TPACC
    else
      zz = TPACC + interval + "h"
    end if
    tvarback1 = grb_file_back1->$zz$(:,:) * 0.03937 ; convert to inches
    tvarback1 = where (.not. ismissing(tvarback1) .and. tvarback1 .lt. 0, 0, tvarback1)
  end if
  if (fnameback2test .eq. 1) then
    if ( isfilevar(grb_file_back2,TPACC) ) then
      zz = TPACC
    else
      zz = TPACC + interval + "h"
    end if
    tvarback2 = grb_file_back2->$zz$(:,:) * 0.03937 ; convert to inches
    tvarback2 = where (.not. ismissing(tvarback2) .and. tvarback2 .lt. 0, 0, tvarback2)
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
  if ( fctm .gt. 0 ) then
    if (isfilevar(grb_file,TPACC)) then
      zz = TPACC
    else
      zz = TPACC + interval + "h"
    end if
    tvar = grb_file->$zz$(:,:) * 0.03937 ; convert to inches
  else
    exit
  end if
;
  tvar = where (.not. ismissing(tvar) .and. tvar .lt. 0, 0, tvar)
  tvarback1 = tvar * 0.
  tvarback2 = tvar * 0.
  if (fnameback1test .eq. 1) then
    if ( isfilevar(grb_file_back1,TPACC) ) then
      zz = TPACC
    else
      zz = TPACC + interval + "h"
    end if
    tvarback1 = grb_file_back1->$zz$(:,:) * 0.03937 ; convert to inches
    tvarback1 = where (.not. ismissing(tvarback1) .and. tvarback1 .lt. 0, 0, tvarback1)
  end if
  if (fnameback2test .eq. 1) then
    if ( isfilevar(grb_file_back2,TPACC) ) then
      zz = TPACC
    else
      zz = TPACC + interval + "h"
    end if
    tvarback2 = grb_file_back2->$zz$(:,:) * 0.03937 ; convert to inches
    tvarback2 = where (.not. ismissing(tvarback2) .and. tvarback2 .lt. 0, 0, tvarback2)
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
  if ( fctm .gt. 0 ) then
    if (isfilevar(grb_file,TPACC)) then
      zz = TPACC
    else
      zz = TPACC + interval + "h"
    end if
    tvar = grb_file->$zz$(::-1,:) * 0.03937 ; convert to inches
  else
    exit
  end if
;
  tvar = where (.not. ismissing(tvar) .and. tvar .lt. 0, 0, tvar)
  tvarback1 = tvar * 0.
  tvarback2 = tvar * 0.
  if (fnameback1test .eq. 1) then
    if ( isfilevar(grb_file_back1,TPACC) ) then
      zz = TPACC
    else
      zz = TPACC + interval + "h"
    end if
    tvarback1 = grb_file_back1->$zz$(::-1,:) * 0.03937 ; convert to inches
    tvarback1 = where (.not. ismissing(tvarback1) .and. tvarback1 .lt. 0, 0, tvarback1)
  end if
  if (fnameback2test .eq. 1) then
    if ( isfilevar(grb_file_back2,TPACC) ) then
      zz = TPACC
    else
      zz = TPACC + interval + "h"
    end if
    tvarback2 = grb_file_back2->$zz$(::-1,:) * 0.03937 ; convert to inches
    tvarback2 = where (.not. ismissing(tvarback2) .and. tvarback2 .lt. 0, 0, tvarback2)
  end if
;
  tvarorig = tvar
  fliptvar = lonFlip(tvar)
  tvarback1orig = tvarback1
  fliptvarback1 = lonFlip(tvarback1)
  tvarback2orig = tvarback2
  fliptvarback2 = lonFlip(tvarback2)
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
    if (isfilevar(grb_file,TPACC)) then
      zz = TPACC
    else
      zz = TPACC + interval + "h"
    end if
    tvar = grb_file->$zz$(:,:) * 0.03937 ; convert to inches
  else
    exit
  end if
;
  tvar = where (.not. ismissing(tvar) .and. tvar .lt. 0, 0, tvar)
  tvarback1 = tvar * 0.
  tvarback2 = tvar * 0.
  if (fnameback1test .eq. 1) then
    if ( isfilevar(grb_file_back1,TPACC) ) then
      zz = TPACC
    else
      zz = TPACC + interval + "h"
    end if
    tvarback1 = grb_file_back1->$zz$(:,:) * 0.03937 ; convert to inches
    tvarback1 = where (.not. ismissing(tvarback1) .and. tvarback1 .lt. 0, 0, tvarback1)
  end if
  if (fnameback2test .eq. 1) then
    if ( isfilevar(grb_file_back2,TPACC) ) then
      zz = TPACC
    else
      zz = TPACC + interval + "h"
    end if
    tvarback2 = grb_file_back2->$zz$(:,:) * 0.03937 ; convert to inches
    tvarback2 = where (.not. ismissing(tvarback2) .and. tvarback2 .lt. 0, 0, tvarback2)
  end if
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
    zz = TPACC + interval + "h"
    if (isfilevar(grb_file,TPACC)) then
      zz = TPACC
    else
      zz = TPACC + interval + "h"
    end if
    tvar = grb_file->$zz$(:,:) * 0.03937 ; convert to inches
  else
    exit
  end if
;
  tvar = where (.not. ismissing(tvar) .and. tvar .lt. 0, 0, tvar)
  tvarback1 = tvar * 0.
  tvarback2 = tvar * 0.
  if (fnameback1test .eq. 1) then
    if ( isfilevar(grb_file_back1,TPACC) ) then
      zz = TPACC
    else
      zz = TPACC + interval + "h"
    end if
    tvarback1 = grb_file_back1->$zz$(:,:) * 0.03937 ; convert to inches
    tvarback1 = where (.not. ismissing(tvarback1) .and. tvarback1 .lt. 0, 0, tvarback1)
  end if
  if (fnameback2test .eq. 1) then
    if ( isfilevar(grb_file_back2,TPACC) ) then
      zz = TPACC
    else
      zz = TPACC + interval + "h"
    end if
    tvarback2 = grb_file_back2->$zz$(:,:) * 0.03937 ; convert to inches
    tvarback2 = where (.not. ismissing(tvarback2) .and. tvarback2 .lt. 0, 0, tvarback2)
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
  if ( fctm .gt. 0 ) then
    if (isfilevar(grb_file,TPACC)) then
      zz = TPACC
    else
      zz = TPACC + interval + "h"
    end if
    tvar = grb_file->$zz$(:,:) * 0.03937 ; convert to inches
  else
    exit
  end if
;
  tvar = where (.not. ismissing(tvar) .and. tvar .lt. 0, 0, tvar)
  tvarback1 = tvar * 0.
  tvarback2 = tvar * 0.
  if (fnameback1test .eq. 1) then
    if ( isfilevar(grb_file_back1,TPACC) ) then
      zz = TPACC
    else
      zz = TPACC + interval + "h"
    end if
    tvarback1 = grb_file_back1->$zz$(:,:) * 0.03937 ; convert to inches
    tvarback1 = where (.not. ismissing(tvarback1) .and. tvarback1 .lt. 0, 0, tvarback1)
  end if
  if (fnameback2test .eq. 1) then
    if ( isfilevar(grb_file_back2,TPACC) ) then
      zz = TPACC
    else
      zz = TPACC + interval + "h"
    end if
    tvarback2 = grb_file_back2->$zz$(:,:) * 0.03937 ; convert to inches
    tvarback2 = where (.not. ismissing(tvarback2) .and. tvarback2 .lt. 0, 0, tvarback2)
  end if
;
end if
;
load "$NCL_HOME/obj/get_date_str_FIM.ncl"
;
print(date_str)
print("acp")
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
           "wkMetaName" : "./esbl_hvyacp.ncgm"
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
           "wkPSFileName" : "./esbl_hvyacp.ps"
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
  tvarbkg = tvar * 0.
  tvarbkg(lly,llx) = .1
;
  tfieldbkg = create "ScalarField" scalarFieldClass appid
      "sfDataArray" : tvarbkg
      "sfXCStartSubsetV" : llx
      "sfYCStartSubsetV" : lly
      "sfXCEndSubsetV" : urx
      "sfYCEndSubsetV" : ury
      "sfMissingValueV" : tvarbkg@_FillValue
  end create
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
  if (fnameback1test .eq. 1) then
    tfieldback1 = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : tvarback1
        "sfXCStartSubsetV" : llx
        "sfYCStartSubsetV" : lly
        "sfXCEndSubsetV" : urx
        "sfYCEndSubsetV" : ury
        "sfMissingValueV" : tvarback1@_FillValue
    end create
  end if
;
  if (fnameback2test .eq. 1) then
    tfieldback2 = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : tvarback2
        "sfXCStartSubsetV" : llx
        "sfYCStartSubsetV" : lly
        "sfXCEndSubsetV" : urx
        "sfYCEndSubsetV" : ury
        "sfMissingValueV" : tvarback2@_FillValue
    end create
  end if
;
  cl0 = .1
  cl1 = .25
  cl2 = .5
  cl3 = 1.0
  cl4 = 2.0
;
  load "$NCL_HOME/obj/esbl.ncl"
;
  setvalues cnid
    "cnRasterModeOn" : False
  end setvalues
;
  titlestring = "TL-Ensemble " + interval + "h Precip (in)"
  load "$NCL_HOME/obj/add_titles.ncl"
  load "$NCL_HOME/obj/add_explabel.ncl"
;
; Draw map
;
  draw(cnidbkg)
;
  if (fnameback2test .eq. 1) then
    draw(cnidback2abv5)
  end if
  if (fnameback1test .eq. 1) then
    draw(cnidback1abv5)
  end if
  draw(cnidcurrentabv5)
;
  if (fnameback2test .eq. 1) then
    draw(cnidback2)
  end if
  if (fnameback1test .eq. 1) then
    draw(cnidback1)
  end if
  draw(cnid)
  draw(mpid)
  draw(lbid)
;
  load "$NCL_HOME/obj/esbl2.ncl"
;
  draw(lbid)
;
  gsres = True
  if (iq .eq. 0) then
    do i = 0,dimsizes(airport_id)-1
      gsres@gsMarkerIndex = 16
      gsres@gsMarkerColor = 55
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
      gsres@gsMarkerColor = 55
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
  tvarbkg = tvar * 0.
  tvarbkg(lly,llx) = .1
;
  tfieldbkg = create "ScalarField" scalarFieldClass appid
      "sfDataArray" : tvarbkg
      "sfXCStartSubsetV" : llx
      "sfYCStartSubsetV" : lly
      "sfXCEndSubsetV" : urx
      "sfYCEndSubsetV" : ury
      "sfMissingValueV" : tvarbkg@_FillValue
  end create
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
  if (fnameback1test .eq. 1) then
    tfieldback1 = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : tvarback1
        "sfXCStartSubsetV" : llx
        "sfYCStartSubsetV" : lly
        "sfXCEndSubsetV" : urx
        "sfYCEndSubsetV" : ury
        "sfMissingValueV" : tvarback1@_FillValue
    end create
  end if
;
  if (fnameback2test .eq. 1) then
    tfieldback2 = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : tvarback2
        "sfXCStartSubsetV" : llx
        "sfYCStartSubsetV" : lly
        "sfXCEndSubsetV" : urx
        "sfYCEndSubsetV" : ury
        "sfMissingValueV" : tvarback2@_FillValue
    end create
  end if
;
  cl0 = .1
  cl1 = .25
  cl2 = .5
  cl3 = 1.0
  cl4 = 2.0
;
  load "$NCL_HOME/obj/esbl.ncl"
;
  setvalues cnid
    "cnRasterModeOn" : False
  end setvalues
;
  titlestring = "TL-Ensemble " + interval + "h Precip (in)"
  load "$NCL_HOME/obj/add_titles.ncl"
  load "$NCL_HOME/obj/add_explabel.ncl"
;
; Draw map
;
  draw(cnidbkg)
;
  if (fnameback2test .eq. 1) then
    draw(cnidback2abv5)
  end if
  if (fnameback1test .eq. 1) then
    draw(cnidback1abv5)
  end if
  draw(cnidcurrentabv5)
;
  if (fnameback2test .eq. 1) then
    draw(cnidback2)
  end if
  if (fnameback1test .eq. 1) then
    draw(cnidback1)
  end if
  draw(cnid)
  draw(mpid)
  draw(lbid)
;
  load "$NCL_HOME/obj/esbl2.ncl"
;
  draw(lbid)
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
  tvarbkg = tvar * 0.
  tvarbkg(lly,llx) = .1
;
  tfieldbkg = create "ScalarField" scalarFieldClass appid
      "sfDataArray" : tvarbkg
      "sfXCStartSubsetV" : llx
      "sfYCStartSubsetV" : lly
      "sfXCEndSubsetV" : urx
      "sfYCEndSubsetV" : ury
      "sfMissingValueV" : tvarbkg@_FillValue
  end create
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
  if (fnameback1test .eq. 1) then
    tfieldback1 = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : tvarback1
        "sfXCStartSubsetV" : llx
        "sfYCStartSubsetV" : lly
        "sfXCEndSubsetV" : urx
        "sfYCEndSubsetV" : ury
        "sfMissingValueV" : tvarback1@_FillValue
    end create
  end if
;
  if (fnameback2test .eq. 1) then
    tfieldback2 = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : tvarback2
        "sfXCStartSubsetV" : llx
        "sfYCStartSubsetV" : lly
        "sfXCEndSubsetV" : urx
        "sfYCEndSubsetV" : ury
        "sfMissingValueV" : tvarback2@_FillValue
    end create
  end if
;
  cl0 = .1
  cl1 = .25
  cl2 = .5
  cl3 = 1.0
  cl4 = 2.0
;
  load "$NCL_HOME/obj/esbl.ncl"
;
  setvalues cnid
    "cnRasterModeOn" : False
  end setvalues
;
  titlestring = "TL-Ensemble " + interval + "h Precip (in)"
  load "$NCL_HOME/obj/add_titles.ncl"
  load "$NCL_HOME/obj/add_explabel.ncl"
;
; Draw map
;
  draw(cnidbkg)
;
  if (fnameback2test .eq. 1) then
    draw(cnidback2abv5)
  end if
  if (fnameback1test .eq. 1) then
    draw(cnidback1abv5)
  end if
  draw(cnidcurrentabv5)
;
  if (fnameback2test .eq. 1) then
    draw(cnidback2)
  end if
  if (fnameback1test .eq. 1) then
    draw(cnidback1)
  end if
  draw(cnid)
  draw(mpid)
  draw(lbid)
;
  load "$NCL_HOME/obj/esbl2.ncl"
;
  draw(lbid)
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
  tvarbkg = tvar * 0.
  tvarbkg(0,0) = .1
;
; Create ScalarField objects.
;
  if (domainname .eq. "fim" .or. domainname .eq. "244" .or. domainname .eq. "129") then
;
    tfieldbkg = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : tvarbkg
        "sfXCStartV" : 0
        "sfYCStartV" : ydim
        "sfXCEndV" : xdim
        "sfYCEndV" : 0
        "sfMissingValueV" : tvarbkg@_FillValue
    end create
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
    if (fnameback1test .eq. 1) then
      tfieldback1 = create "ScalarField" scalarFieldClass appid
          "sfDataArray" : tvarback1
          "sfXCStartV" : 0
          "sfYCStartV" : ydim
          "sfXCEndV" : xdim
          "sfYCEndV" : 0
          "sfMissingValueV" : tvarback1@_FillValue
      end create
    end if
;
    if (fnameback2test .eq. 1) then
      tfieldback2 = create "ScalarField" scalarFieldClass appid
          "sfDataArray" : tvarback2
          "sfXCStartV" : 0
          "sfYCStartV" : ydim
          "sfXCEndV" : xdim
          "sfYCEndV" : 0
          "sfMissingValueV" : tvarback2@_FillValue
      end create
    end if
;
  end if
;
  if (domainname .eq. "236" .or. domainname .eq. "201" .or. domainname .eq. "224" .or. \
      domainname .eq. "83") then
;
    tfieldbkg = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : tvarbkg
        "sfXCStartV" : 0
        "sfYCStartV" : 0
        "sfXCEndV" : xdim
        "sfYCEndV" : ydim
        "sfMissingValueV" : tvarbkg@_FillValue
    end create
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
    if (fnameback1test .eq. 1) then
      tfieldback1 = create "ScalarField" scalarFieldClass appid
          "sfDataArray" : tvarback1
          "sfXCStartV" : 0
          "sfYCStartV" : 0
          "sfXCEndV" : xdim
          "sfYCEndV" : ydim
          "sfMissingValueV" : tvarback1@_FillValue
      end create
    end if
;
    if (fnameback2test .eq. 1) then
      tfieldback2 = create "ScalarField" scalarFieldClass appid
          "sfDataArray" : tvarback2
          "sfXCStartV" : 0
          "sfYCStartV" : 0
          "sfXCEndV" : xdim
          "sfYCEndV" : ydim
          "sfMissingValueV" : tvarback2@_FillValue
      end create
    end if
;
  end if
;
  cl0 = .1
  cl1 = .25
  cl2 = .5
  cl3 = 1.0
  cl4 = 2.0
;
  load "$NCL_HOME/obj/esbl.ncl"
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
  titlestring = "TL-Ensemble " + interval + "h Precip (in)"
  load "$NCL_HOME/obj/add_titles.ncl"
  load "$NCL_HOME/obj/add_explabel.ncl"
;
; Draw map
;
  draw(cnidbkg)
;
  if (fnameback2test .eq. 1) then
    draw(cnidback2abv5)
  end if
  if (fnameback1test .eq. 1) then
    draw(cnidback1abv5)
  end if
  draw(cnidcurrentabv5)
;
  if (fnameback2test .eq. 1) then
    draw(cnidback2)
  end if
  if (fnameback1test .eq. 1) then
    draw(cnidback1)
  end if
  draw(cnid)
  draw(mpid)
  draw(lbid)
;
  load "$NCL_HOME/obj/esbl2.ncl"
;
  draw(lbid)
;
  gsres = True
  if (domainname .eq. "236") then
    do i = 0,dimsizes(airport_id)-1
      gsres@gsMarkerIndex = 16
      gsres@gsMarkerColor = 55
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

