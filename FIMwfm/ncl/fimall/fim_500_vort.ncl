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
;
; Get the 500mb level index
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
  uvar = grb_file->$UISBL$(level,:,:) 
  vvar = grb_file->$VISBL$(level,:,:)
  vort = new((/ny,nx/),float)
  flipuvar = new((/ny,nx/),float)
  flipvvar = new((/ny,nx/),float)
  fliplat = new((/ny/),float)
  do j = 0,ny-1
    do k = 0,nx-1
      fj = (ny-1) - j
      fliplat(j) = lat(fj)
      flipuvar(j,k) = uvar(fj,k)
      flipvvar(j,k) = vvar(fj,k)
    end do
  end do
  lat = fliplat
  uvar = flipuvar
  vvar = flipvvar
;
  uv2vrf(uvar,vvar,vort)
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
  uvar = grb_file->$UISBL$(level,:,:)
  vvar = grb_file->$VISBL$(level,:,:)
  vort = new((/ny,nx/),float)
  flipuvar = new((/ny,nx/),float)
  flipvvar = new((/ny,nx/),float)
  fliplat = new((/ny/),float)
  do j = 0,ny-1
    do k = 0,nx-1
      fj = (ny-1) - j
      fliplat(j) = lat(fj)
      flipuvar(j,k) = uvar(fj,k)
      flipvvar(j,k) = vvar(fj,k)
    end do
  end do
  lat = fliplat
  uvar = flipuvar
  vvar = flipvvar
;
  uv2vrf(uvar,vvar,vort)
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
  uvar = grb_file->$UISBL$(level,:,:)
  vvar = grb_file->$VISBL$(level,:,:)
  vort = new((/ny,nx/),float)
;
  uv2vrf(uvar,vvar,vort)
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
  rot = grb_file->$ROT$(:,:) ; use to rotate winds
  latin1 = grb_file->$LON$@Latin1
  latin2 = grb_file->$LON$@Latin2
  lov = grb_file->$LON$@Lov
  hvar = grb_file->$GHISBL$(level,:,:) / 10. ; use decameters
;  uvar = grb_file->$UISBL$(level,:,:)
;  vvar = grb_file->$VISBL$(level,:,:)
  uvartmp = grb_file->$UISBL$(level,:,:) ; use to rotate winds
  vvartmp = grb_file->$VISBL$(level,:,:) ; use to rotate winds
  vort = new((/ny,nx/),float)
;
; Rotate winds to earth-relative
;
  uearth = sin(rot)*vvartmp + cos(rot)*uvartmp
  vearth = cos(rot)*vvartmp - sin(rot)*uvartmp
  uvar = uearth
  vvar = vearth
;
  do j = 1,ny-2
    do k = 1,nx-2
      dx = 13545.087
      dy = 13545.087
      jm1 = j - 1
      jp1 = j + 1
      km1 = k - 1
      kp1 = k + 1
      dv = vvar(j,kp1) - vvar(j,km1)
      du = uvar(jp1,k) - uvar(jm1,k)
      vort(j,k) = ((dv/(2.*dx)) - (du/(2.*dy)))
    end do
  end do
  do k = 0,nx-1
    vort(0,k) = vort(1,k)
    vort(ny-1,k) = vort(ny-2,k)
  end do
  do j = 0,ny-1
    vort(j,0) = vort(j,1)
    vort(j,nx-1) = vort(j,nx-2)
  end do
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
  uvar = grb_file->$UISBL$(level,::-1,:)
  vvar = grb_file->$VISBL$(level,::-1,:)
;  uvar = grb_file->$UISBL$(level,::-1,:)
;  vvar = grb_file->$VISBL$(level,::-1,:)
;  flipuvar = new((/ny,nx/),float)
;  flipvvar = new((/ny,nx/),float)
;  fliplat = new((/ny/),float)
;  do j = 0,ny-1
;    do k = 0,nx-1
;      fj = (ny-1) - j
;      fliplat(j) = lat(fj)
;      flipuvar(j,k) = uvar(fj,k)
;      flipvvar(j,k) = vvar(fj,k)
;    end do
;  end do
;  lat = fliplat
;  uvar = flipuvar
;  vvar = flipvvar
;
  vort = new((/ny,nx/),float)
  uv2vrf(uvar,vvar,vort)
  delete(uvar)
  delete(vvar)
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
;  rot = grb_file->$ROT$(:,:) ; use to rotate winds
  latin1 = grb_file->$LON$@Latin1
  latin2 = grb_file->$LON$@Latin2
  lov = grb_file->$LON$@Lov
  hvar = grb_file->$GHISBL$(level,:,:) / 10. ; use decameters
  uvar = grb_file->$UISBL$(level,:,:)
  vvar = grb_file->$VISBL$(level,:,:)
;  uvartmp = grb_file->$UISBL$(level,:,:) ; use to rotate winds
;  vvartmp = grb_file->$VISBL$(level,:,:) ; use to rotate winds
  vort = new((/ny,nx/),float)
;;
;; Rotate winds to earth-relative
;;
;  uearth = sin(rot)*vvartmp + cos(rot)*uvartmp
;  vearth = cos(rot)*vvartmp - sin(rot)*uvartmp
;  do jj = 0,ny-1
;   do ii = 0,nx-1
;      print("uvar ("+jj+","+ii+") ="+uvar(jj,ii)+"  uearth ("+jj+","+ii+") ="+uearth(jj,ii))
;    end do
;  end do
;  uvar = uearth
;  vvar = vearth
;
  do j = 1,ny-2
    do k = 1,nx-2
      dx = 40635.3
      dy = 40635.3
      jm1 = j - 1
      jp1 = j + 1
      km1 = k - 1
      kp1 = k + 1
      dv = vvar(j,kp1) - vvar(j,km1)
      du = uvar(jp1,k) - uvar(jm1,k)
      vort(j,k) = ((dv/(2.*dx)) - (du/(2.*dy)))
    end do
  end do
  do k = 0,nx-1
    vort(0,k) = vort(1,k)
    vort(ny-1,k) = vort(ny-2,k)
  end do
  do j = 0,ny-1
    vort(j,0) = vort(j,1)
    vort(j,nx-1) = vort(j,nx-2)
  end do
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
  uvar = cos(rot)*uvartmp - sin(rot)*vvartmp
  vvar = sin(rot)*uvartmp + cos(rot)*vvartmp
;
  vort = new((/ny,nx/),float)
  do j = 1,ny-2
    do k = 1,nx-2
      dx = 94512.0
      dy = 94512.0
      jm1 = j - 1
      jp1 = j + 1
      km1 = k - 1
      kp1 = k + 1
      dv = vvar(j,kp1) - vvar(j,km1)
      du = uvar(jp1,k) - uvar(jm1,k)
      vort(j,k) = ((dv/(2.*dx)) - (du/(2.*dy)))
    end do
  end do
  do k = 0,nx-1
    vort(0,k) = vort(1,k)
    vort(ny-1,k) = vort(ny-2,k)
  end do
  do j = 0,ny-1
    vort(j,0) = vort(j,1)
    vort(j,nx-1) = vort(j,nx-2)
  end do
;
end if
;
if (domainname .eq. "244") then
;
  ny = filevardimsizes(grb_file,LAT)
  nx = filevardimsizes(grb_file,LON)
;
  lat = grb_file->$LAT$(:)
  lat = lat - 50.955
  lon = grb_file->$LON$(:)
  hvar = grb_file->$GHISBL$(level,:,:) / 10. ; use decameters
  uvar = grb_file->$UISBL$(level,:,:)
  vvar = grb_file->$VISBL$(level,:,:)
  vort = new((/ny,nx/),float)
  flipuvar = new((/ny,nx/),float)
  flipvvar = new((/ny,nx/),float)
  fliplat = new((/ny/),float)
  do j = 0,ny-1
    do k = 0,nx-1
      fj = (ny-1) - j
      fliplat(j) = lat(fj)
      flipuvar(j,k) = uvar(fj,k)
      flipvvar(j,k) = vvar(fj,k)
    end do
  end do
  lat = fliplat
;
  vort = uv2vr_cfd(uvar,vvar,lat,lon,2)
  vort = -vort
;
end if
;
load "$NCL_HOME/obj/get_date_str_FIM.ncl"
;
print(date_str)
vort = vort * 1.e05
vortorig = vort
if (domainname .eq. "174") then
  flipvort = lonFlip(vort)
end if
print("500 mb vort")
print(min(vort))
print(max(vort))
print("gh")
print(min(hvar))
print(max(hvar))
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
           "wkMetaName" : "./500_vort.ncgm"
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
           "wkPSFileName" : "./500_vort.ps"
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
   "wsMaximumSize": 494217728
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
      "sfDataArray" : vort
      "sfXCStartSubsetV" : llx
      "sfYCStartSubsetV" : lly
      "sfXCEndSubsetV" : urx
      "sfYCEndSubsetV" : ury
      "sfMissingValueV" : vort@_FillValue
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
  load "$NCL_HOME/obj/diff_-16-16by2.ncl"
;
  setvalues cnid
    "cnFillColors" : (/2,3,4,5,6,7,8,9,20,11,12,13,14,15,16,17,18/)
    "cnRasterModeOn" : False
  end setvalues
;
  load "$NCL_HOME/obj/height_60m.ncl"
;
  setvalues cnid4
    "cnLineColor" : 1
  end setvalues
;
  if (iq .ne. 0) then
    setvalues cnid4
      "cnLineThicknessF" : 2.5
    end setvalues
  end if
;
  titlestring = "500mb Rel Vort (10E-5/s), Height (dm)"
  load "$NCL_HOME/obj/add_titles.ncl"
  load "$NCL_HOME/obj/add_explabel_cnid4.ncl"
;
; Draw map
;
  draw(cnid)
  draw(mpid)
  draw(cnid4)
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
      "sfDataArray" : vort
      "sfXCStartSubsetV" : llx
      "sfYCStartSubsetV" : lly
      "sfXCEndSubsetV" : urx
      "sfYCEndSubsetV" : ury
      "sfMissingValueV" : vort@_FillValue
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
  load "$NCL_HOME/obj/diff_-16-16by2.ncl"
;
  setvalues cnid
    "cnFillColors" : (/2,3,4,5,6,7,8,9,20,11,12,13,14,15,16,17,18/)
  end setvalues
;
  load "$NCL_HOME/obj/height_60m.ncl"
;
  setvalues cnid4
    "cnLineColor" : 1
  end setvalues
;
  titlestring = "500mb Rel Vort (10E-5/s), Height (dm)"
  load "$NCL_HOME/obj/add_titles.ncl"
  load "$NCL_HOME/obj/add_explabel_cnid4.ncl"
;
; Draw map
;
  draw(cnid)
  draw(mpid)
  draw(cnid4)
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
      "sfDataArray" : vort
      "sfXCStartSubsetV" : llx
      "sfYCStartSubsetV" : lly
      "sfXCEndSubsetV" : urx
      "sfYCEndSubsetV" : ury
      "sfMissingValueV" : vort@_FillValue
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
  load "$NCL_HOME/obj/diff_-16-16by2.ncl"
;
  setvalues cnid
    "cnFillColors" : (/2,3,4,5,6,7,8,9,20,11,12,13,14,15,16,17,18/)
  end setvalues
;
  load "$NCL_HOME/obj/height_60m.ncl"
;
  setvalues cnid
    "cnFillMode" : "AreaFill"
  end setvalues
;
  setvalues cnid4
    "cnLineColor" : 1
  end setvalues
;
  titlestring = "500mb Rel Vort (10E-5/s), Height (dm)"
  load "$NCL_HOME/obj/add_titles.ncl"
  load "$NCL_HOME/obj/add_explabel_cnid4.ncl"
;
; Draw map
;
  draw(cnid)
  draw(mpid)
  draw(cnid4)
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
  if (domainname .eq. "fim" .or. domainname .eq. "129") then
;
    tfield = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : vort
        "sfXCStartV" : 0
        "sfYCStartV" : 0
        "sfXCEndV" : xdim
        "sfYCEndV" : ydim
        "sfMissingValueV" : vort@_FillValue
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
  end if
;
  if (domainname .eq. "244") then
;
    tfield = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : vort
        "sfXCStartV" : 0
        "sfYCStartV" : ydim
        "sfXCEndV" : xdim
        "sfYCEndV" : 0
        "sfMissingValueV" : vort@_FillValue
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
  end if
;
  if (domainname .eq. "236" .or. domainname .eq. "201" .or. domainname .eq. "224" .or. \
      domainname .eq. "83") then
;
    tfield = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : vort
        "sfXCStartV" : 0
        "sfYCStartV" : 0
        "sfXCEndV" : xdim
        "sfYCEndV" : ydim
        "sfMissingValueV" : vort@_FillValue
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
  end if
;
  load "$NCL_HOME/obj/diff_-16-16by2.ncl"
;
  if (domainname .eq. "236" .or. domainname .eq. "244" .or. domainname .eq. "83" .or. \
      domainname .eq. "201" .or. domainname .eq. "224") then
;
    setvalues cnid
      "cnFillColors" : (/2,3,4,5,6,7,8,9,20,11,12,13,14,15,16,17,18/)
      "cnRasterModeOn" : False
    end setvalues
;
  else
;
    setvalues cnid
      "cnFillColors" : (/2,3,4,5,6,7,8,9,20,11,12,13,14,15,16,17,18/)
    end setvalues
;
  end if
;
  load "$NCL_HOME/obj/height_60m.ncl"
;
  setvalues cnid4
    "cnLineColor" : 1
  end setvalues
;
  if (domainname .eq. "fim" .or. domainname .eq. "129") then
    setvalues cnid4
      "cnLineThicknessF" : 1
      "cnLineLabelDensityF" : 0.001
    end setvalues
  end if
;
  if (domainname .eq. "201") then
    setvalues cnid4
      "cnLineLabelDensityF" : 0.001
    end setvalues
  end if
;
  titlestring = "500mb Rel Vort (10E-5/s), Height (dm)"
  load "$NCL_HOME/obj/add_titles.ncl"
  load "$NCL_HOME/obj/add_explabel_cnid4.ncl"
;
; Draw map
;
  draw(cnid)
  draw(mpid)
  draw(cnid4)
;
  gsres = True
  if (domainname .eq. "236") then
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
  end if
;
  frame(wid)
;
end if  ; domainname .ne. "130" .and. domainname .ne. "174" .and. domainname .ne. "242"
;
end
;
