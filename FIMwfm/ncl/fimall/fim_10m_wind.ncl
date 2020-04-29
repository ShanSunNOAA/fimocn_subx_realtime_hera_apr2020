;
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_LIB/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
modl = getenv("MODL")
;
load "$NCL_HOME/obj/cmap1_track.ncl"
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
    wstrid = 10
;
    uvar = grb_file->U_GRD_GDS0_HTGL(::-1,:) * 1.9429 ; convert to knots
    vvar = grb_file->V_GRD_GDS0_HTGL(::-1,:) * 1.9429 ; convert to knots
;
    uvar = lonFlip(uvar)
    vvar = lonFlip(vvar)
;
  else
;
    ny = filevardimsizes(grb_file,LAT)
    nx = filevardimsizes(grb_file,LON)
;
    lat = grb_file->$LAT$(:)
    lon = grb_file->$LON$(:)
    wstrid = 10
;
    uvar = grb_file->$U10M$(:,:) * 1.9429 ; convert to knots
    vvar = grb_file->$V10M$(:,:) * 1.9429 ; convert to knots
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
  wstrid = 20
;
  uvar = grb_file->$U10M$(:,:) * 1.9429 ; convert to knots
  vvar = grb_file->$V10M$(:,:) * 1.9429 ; convert to knots
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
  wstrid = 10
;
  uvar = grb_file->$U10M$(:,:) * 1.9429 ; convert to knots
  vvar = grb_file->$V10M$(:,:) * 1.9429 ; convert to knots
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
  wstrid = 7
;
  uvar = grb_file->$U10M$(:,:) * 1.9429 ; convert to knots
  vvar = grb_file->$V10M$(:,:) * 1.9429 ; convert to knots
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
  wstrid = 50
;
  uvarorig = grb_file->$U10M$(::-1,:) * 1.9429 ; convert to knots
  uvar = uvarorig
  flipuvar = lonFlip(uvar)
  vvarorig = grb_file->$V10M$(::-1,:) * 1.9429 ; convert to knots
  vvar = vvarorig
  flipvvar = lonFlip(vvar)
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
  wstrid = 3
;
  uvar = grb_file->$U10M$(:,:) * 1.9429 ; convert to knots
  vvar = grb_file->$V10M$(:,:) * 1.9429 ; convert to knots
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
  uvartmp = grb_file->$U10M$(:,:)
  vvartmp = grb_file->$V10M$(:,:)
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
  wstrid = 5
;
  uvar = grb_file->$U10M$(:,:) * 1.9429 ; convert to knots
  vvar = grb_file->$V10M$(:,:) * 1.9429 ; convert to knots
;
end if
;
load "$NCL_HOME/obj/get_date_str_FIM.ncl"
;
uv2 = uvar * uvar
vv2 = vvar * vvar
wind = sqrt(uv2 + vv2) + .00001
;
if (domainname .eq. "174") then
  windorig = wind
  flipwind = lonFlip(wind)
end if
;
print(date_str)
print("wind spd")
print(min(wind))
print(max(wind))
;
if (max(wind) .le. 5) then ; set a corner pt to make a plot anyway
  wind(0,0) = 5.001
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
           "wkMetaName" : "./10m_wind.ncgm"
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
           "wkPSFileName" : "./10m_wind.ps"
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
; Create ScalarField and VectorField objects.
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
      "mpNationalLineColor" : "foreground"
    end setvalues
  end if
;
  tfield = create "ScalarField" scalarFieldClass appid
      "sfDataArray" : wind
      "sfXCStartSubsetV" : llx
      "sfYCStartSubsetV" : lly
      "sfXCEndSubsetV" : urx
      "sfYCEndSubsetV" : ury
      "sfMissingValueV" : wind@_FillValue
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
; Create a VectorPlot object.
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

  load "$NCL_HOME/obj/wind_5-90by5.ncl"
;
  titlestring = "10m Wind (kt)"
  load "$NCL_HOME/obj/add_titles.ncl"
  load "$NCL_HOME/obj/add_explabel_vcid.ncl"
;
  if (fileexists("track_file.txt")) then
    trkfname = asciiread("track_file.txt",-1,"string")
    print("trkfname = "+trkfname)
    trkpts = asciiread(trkfname,-1,"string")
;    print("trkpts = "+trkpts)
;    print(trkpts(0:40))
    trklat = new((/500,dimsizes(trkpts)/),float)
    tmplat = new((/dimsizes(trkpts)/),float)
    trklon = new((/500,dimsizes(trkpts)/),float)
    tmplon = new((/dimsizes(trkpts)/),float)
    latptcnt = new(500,integer)
    lonptcnt = new(500,integer)
    do i = 0, 499
      latptcnt(i) = 0
      lonptcnt(i) = 0
    end do
    print("dimsizes(trkpts) = "+dimsizes(trkpts))
    element = 0
    do while (element .lt. dimsizes(trkpts))
      storm = stringtoint(getSubString(trkpts(element),3,5))
      trk_region = getSubString(trkpts(element),0,1)
      if (trk_region .eq. "EP") then
        storm = storm + 100
      end if
      if (trk_region .eq. "AL") then
        storm = storm + 200
      end if
      if (trk_region .eq. "WP") then
        storm = storm + 300
      end if
;      print("storm "+storm)
      fhour = stringtoint(getSubString(trkpts(element),30,32))
      if (fhour .le. fctm) then
        trklat(storm,latptcnt(storm)) = stringtoint(getSubString(trkpts(element),34,37)) / 10.
        trklat_ns = getSubString(trkpts(element),38,38)
        if (trklat_ns .eq. "S") then
          trklat(storm,latptcnt(storm)) = -trklat(storm,latptcnt(storm))
        end if
;        print("lat = "+trklat(storm,latptcnt(storm)))
        latptcnt(storm) = latptcnt(storm) + 1
        trklon(storm,lonptcnt(storm)) = stringtoint(getSubString(trkpts(element),41,44)) / 10.
        trklon_ew = getSubString(trkpts(element),45,45)
        if (trklon_ew .eq. "W") then
          trklon(storm,lonptcnt(storm)) = 360. - trklon(storm,lonptcnt(storm))
        end if
;        print("lon = "+trklon(storm,lonptcnt(storm)))
        lonptcnt(storm) = lonptcnt(storm) + 1
      end if
      element = element + 1 ; skip to next line
    end do
    end if
;
; Draw map
;
  draw(cnid)
;
  if (fileexists("track_file.txt")) then
    gsres = True
    gsres@gsMarkerIndex = 1
    gsres@gsMarkerSizeF = 0.04
    gsres@gsMarkerColor = 23
    gsres@gsLineColor = 23
    gsres@gsLineThicknessF = 3
    txres = True
    txres@txFont = 37
    txres@txFontColor = 1
    txres@txFontHeightF = 0.015
    txres@txFontThicknessF = 1.5
;
    do i = 0, 399
      if (.not.ismissing(trklon(i,0))) then
        firstlat = trklat(i,0)
        firstlon = trklon(i,0)
        gsn_text(wid,mpid,"p",firstlon,firstlat,txres)
        delete(firstlat)
        delete(firstlon)
        tmplat = trklat(i,:)
        tmplon = trklon(i,:)
        gsn_polyline(wid,mpid,tmplon,tmplat,gsres)
        delete(tmplat)
        delete(tmplon)
      end if
    end do
  end if
;
  draw(mpid)
  draw(vcid)
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
  uv2 = uvar * uvar
  vv2 = vvar * vvar
  wind = sqrt(uv2 + vv2) + .00001
;
  if (max(wind) .le. 5) then ; set a corner pt to make a plot anyway
    wind(0,0) = 5.001
  end if
;
  if (idom .eq. 0 .or. idom .eq. 2) then
    load "$NCL_HOME/obj/fim_maps_cent0_blk.ncl"
  else
    load "$NCL_HOME/obj/fim_maps_cent180_blk.ncl"
  end if
;
  tfield = create "ScalarField" scalarFieldClass appid
      "sfDataArray" : wind
      "sfXCStartSubsetV" : llx
      "sfYCStartSubsetV" : lly
      "sfXCEndSubsetV" : urx
      "sfYCEndSubsetV" : ury
      "sfMissingValueV" : wind@_FillValue
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
; Create a VectorPlot object.
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
  load "$NCL_HOME/obj/wind_5-90by5.ncl"
;
  titlestring = "10m Wind (kt)"
  load "$NCL_HOME/obj/add_titles.ncl"
  load "$NCL_HOME/obj/add_explabel_vcid.ncl"
;
  if (fileexists("track_file.txt")) then
    trkfname = asciiread("track_file.txt",-1,"string")
    print("trkfname = "+trkfname)
    trkpts = asciiread(trkfname,-1,"string")
;    print("trkpts = "+trkpts)
;    print(trkpts(0:40))
    trklat = new((/500,dimsizes(trkpts)/),float)
    tmplat = new((/dimsizes(trkpts)/),float)
    trklon = new((/500,dimsizes(trkpts)/),float)
    tmplon = new((/dimsizes(trkpts)/),float)
    latptcnt = new(500,integer)
    lonptcnt = new(500,integer)
    do i = 0, 499
      latptcnt(i) = 0
      lonptcnt(i) = 0
    end do
    print("dimsizes(trkpts) = "+dimsizes(trkpts))
    element = 0
    do while (element .lt. dimsizes(trkpts))
      storm = stringtoint(getSubString(trkpts(element),3,5))
      trk_region = getSubString(trkpts(element),0,1)
      if (trk_region .eq. "EP") then
        storm = storm + 100
      end if
      if (trk_region .eq. "AL") then
        storm = storm + 200
      end if
      if (trk_region .eq. "WP") then
        storm = storm + 300
      end if
;      print("storm "+storm)
      fhour = stringtoint(getSubString(trkpts(element),30,32))
      if (fhour .le. fctm) then
        trklat(storm,latptcnt(storm)) = stringtoint(getSubString(trkpts(element),34,37)) / 10.
        trklat_ns = getSubString(trkpts(element),38,38)
        if (trklat_ns .eq. "S") then
          trklat(storm,latptcnt(storm)) = -trklat(storm,latptcnt(storm))
        end if
;        print("lat = "+trklat(storm,latptcnt(storm)))
        latptcnt(storm) = latptcnt(storm) + 1
        trklon(storm,lonptcnt(storm)) = stringtoint(getSubString(trkpts(element),41,44)) / 10.
        trklon_ew = getSubString(trkpts(element),45,45)
        if (trklon_ew .eq. "W") then
          trklon(storm,lonptcnt(storm)) = 360. - trklon(storm,lonptcnt(storm))
        end if
;         print("lon = "+trklon(storm,lonptcnt(storm)))
        lonptcnt(storm) = lonptcnt(storm) + 1
      end if
      element = element + 1 ; skip to next line
    end do
    end if
;
; Draw map
;
  draw(cnid)
;
  if (fileexists("track_file.txt")) then
    gsres = True
    gsres@gsMarkerIndex = 1
    gsres@gsMarkerSizeF = 0.04
    gsres@gsMarkerColor = 23
    gsres@gsLineColor = 23
    gsres@gsLineThicknessF = 3
    txres = True
    txres@txFont = 37
    txres@txFontColor = 1
    txres@txFontHeightF = 0.015
    txres@txFontThicknessF = 1.5
;
    do i = 0, 399
      if (.not.ismissing(trklon(i,0))) then
        firstlat = trklat(i,0)
        firstlon = trklon(i,0)
        gsn_text(wid,mpid,"p",firstlon,firstlat,txres)
        delete(firstlat)
        delete(firstlon)
        tmplat = trklat(i,:)
        tmplon = trklon(i,:)
        gsn_polyline(wid,mpid,tmplon,tmplat,gsres)
        delete(tmplat)
        delete(tmplon)
      end if
    end do
  end if
;
  draw(mpid)
  draw(vcid)
;
  gsres = True
  if (idom .eq. 1) then
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
;  else
;    do i = 0,dimsizes(airport_id)-1
;      gsres@gsMarkerIndex = 16
;      gsres@gsMarkerColor = 20
;      gsres@gsMarkerSizeF = 0.01
;      gsn_polymarker(wid,mpid,airport_loc(i,1),airport_loc(i,0),gsres)
;      gsres@gsMarkerIndex = 16
;      gsres@gsMarkerColor = 1
;      gsres@gsMarkerSizeF = 0.005
;      gsn_polymarker(wid,mpid,airport_loc(i,1),airport_loc(i,0),gsres)
;    end do
  end if
;
  frame(wid)
;
  end do  ; loop begins in fim_174subdomain_loop.ncl
;
end if  ; domainname .eq. "174"
;
if (domainname .eq. "242") then
;
  uv2 = uvar * uvar
  vv2 = vvar * vvar
  wind = sqrt(uv2 + vv2) + .00001
  windsave = wind
;
  load "$NCL_HOME/obj/242_tiles_loop.ncl"
;
  print("LcLat = "+LcLat)
  print("LcLon = "+LcLon)
  print("RcLat = "+RcLat)
  print("RcLon = "+RcLon)
  wind = windsave
  if (max(wind(lly:ury,llx:urx)) .le. 5) then ; set a corner pt to make a plot anyway
    wind(lly,llx) = 5.001
  end if
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
      "sfDataArray" : wind
      "sfXCStartSubsetV" : llx
      "sfYCStartSubsetV" : lly
      "sfXCEndSubsetV" : urx
      "sfYCEndSubsetV" : ury
      "sfMissingValueV" : wind@_FillValue
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
  setvalues vcid
    "vcMinMagnitudeF" : 0
  end setvalues
;
  load "$NCL_HOME/obj/wind_5-90by5.ncl"
;
  if (iq .ne. 0) then
    setvalues cnid
      "cnFillMode" : "AreaFill"
    end setvalues
  else
    setvalues cnid
      "cnFillMode" : "RasterFill"
    end setvalues
  end if
;
  titlestring = "10m Wind (kt)"
  load "$NCL_HOME/obj/add_titles.ncl"
  load "$NCL_HOME/obj/add_explabel_vcid.ncl"
;
; Draw map
;
  draw(cnid)
  draw(mpid)
  draw(vcid)
;
  frame(wid)
;
  end do  ; loop begins in 242_tiles_loop.ncl
;
end if  ; domainname .eq. "242"
;
if (domainname .ne. "130" .and. domainname .ne. "174") then
;
  load "$NCL_HOME/obj/fim_maps_blk.ncl"
;
; Create ScalarField and VectorField objects.
;
  if (domainname .eq. "fim" .or. domainname .eq. "244" .or. domainname .eq. "129") then
;
    tfield = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : wind
        "sfXCStartV" : 0
        "sfYCStartV" : ydim
        "sfXCEndV" : xdim
        "sfYCEndV" : 0
        "sfMissingValueV" : wind@_FillValue
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
  if (domainname .eq. "236" .or. domainname .eq. "201" .or. domainname .eq. "224" .or. \
      domainname .eq. "83") then
;
    tfield = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : wind
        "sfXCStartV" : 0
        "sfYCStartV" : 0
        "sfXCEndV" : xdim
        "sfYCEndV" : ydim
        "sfMissingValueV" : wind@_FillValue
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
  load "$NCL_HOME/obj/wind_5-90by5.ncl"
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
  titlestring = "10m Wind (kt)"
  load "$NCL_HOME/obj/add_titles.ncl"
  load "$NCL_HOME/obj/add_explabel.ncl"
;
  if (domainname .eq. "fim" .or. domainname .eq. "244" .or. \
      domainname .eq. "83" .or. domainname .eq. "129") then
    if (fileexists("track_file.txt")) then
      trkfname = asciiread("track_file.txt",-1,"string")
      print("trkfname = "+trkfname)
      trkpts = asciiread(trkfname,-1,"string")
;      print("trkpts = "+trkpts)
;      print(trkpts(0:40))
      trklat = new((/500,dimsizes(trkpts)/),float)
      tmplat = new((/dimsizes(trkpts)/),float)
      trklon = new((/500,dimsizes(trkpts)/),float)
      tmplon = new((/dimsizes(trkpts)/),float)
      latptcnt = new(500,integer)
      lonptcnt = new(500,integer)
      do i = 0, 499
        latptcnt(i) = 0
        lonptcnt(i) = 0
      end do
      print("dimsizes(trkpts) = "+dimsizes(trkpts))
      element = 0
      do while (element .lt. dimsizes(trkpts))
        storm = stringtoint(getSubString(trkpts(element),3,5))
        trk_region = getSubString(trkpts(element),0,1)
        if (trk_region .eq. "EP") then
          storm = storm + 100
        end if
        if (trk_region .eq. "AL") then
          storm = storm + 200
        end if
        if (trk_region .eq. "WP") then
          storm = storm + 300
        end if
;        print("storm "+storm)
        fhour = stringtoint(getSubString(trkpts(element),30,32))
        if (fhour .le. fctm) then
          trklat(storm,latptcnt(storm)) = stringtoint(getSubString(trkpts(element),34,37)) / 10.
          trklat_ns = getSubString(trkpts(element),38,38)
          if (trklat_ns .eq. "S") then
            trklat(storm,latptcnt(storm)) = -trklat(storm,latptcnt(storm))
          end if
;          print("lat = "+trklat(storm,latptcnt(storm)))
          latptcnt(storm) = latptcnt(storm) + 1
          trklon(storm,lonptcnt(storm)) = stringtoint(getSubString(trkpts(element),41,44)) / 10.
          trklon_ew = getSubString(trkpts(element),45,45)
          if (trklon_ew .eq. "W") then
            trklon(storm,lonptcnt(storm)) = 360. - trklon(storm,lonptcnt(storm))
          end if
;          print("lon = "+trklon(storm,lonptcnt(storm)))
          lonptcnt(storm) = lonptcnt(storm) + 1
        end if
        element = element + 1 ; skip to next line
      end do
    end if
  end if
;
; Draw map
;
  draw(cnid)
;
  if (domainname .eq. "fim" .or. domainname .eq. "129") then
    if (fileexists("track_file.txt")) then
      gsres = True
      gsres@gsMarkerIndex = 1
      gsres@gsMarkerSizeF = 0.005
      gsres@gsMarkerColor = 1
      gsres@gsLineColor = 23
      gsres@gsLineThicknessF = 2
      txres = True
      txres@txFont = 37
      txres@txFontColor = 1
      txres@txFontHeightF = 0.01
      txres@txFontThicknessF = 1
;
      do i = 0, 399
        if (.not.ismissing(trklon(i,0))) then
          firstlat = trklat(i,0)
          firstlon = trklon(i,0)
          gsn_text(wid,mpid,"p",firstlon,firstlat,txres)
          delete(firstlat)
          delete(firstlon)
          tmplat = trklat(i,:)
          tmplon = trklon(i,:)
          gsn_polyline(wid,mpid,tmplon,tmplat,gsres)
          delete(tmplat)
          delete(tmplon)
        end if
      end do
    end if
  end if
;
  if (domainname .eq. "244" .or. domainname .eq. "83") then
    if (fileexists("track_file.txt")) then
      gsres = True
      gsres@gsMarkerIndex = 1
      gsres@gsMarkerSizeF = 0.04
      gsres@gsMarkerColor = 23
      gsres@gsLineColor = 23
      gsres@gsLineThicknessF = 3
      txres = True
      txres@txFont = 37
      txres@txFontColor = 1
      txres@txFontHeightF = 0.015
      txres@txFontThicknessF = 1.5
;
      do i = 0, 399
        if (.not.ismissing(trklon(i,0))) then
          firstlat = trklat(i,0)
          firstlon = trklon(i,0)
          gsn_text(wid,mpid,"p",firstlon,firstlat,txres)
          delete(firstlat)
          delete(firstlon)
          tmplat = trklat(i,:)
          tmplon = trklon(i,:)
          gsn_polyline(wid,mpid,tmplon,tmplat,gsres)
          delete(tmplat)
          delete(tmplon)
        end if
      end do
    end if
  end if
;
  draw(vcid)
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
