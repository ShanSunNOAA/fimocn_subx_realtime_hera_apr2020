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
domainname = asciiread("domain.txt",-1,"string")
;
; Get the 500mb level index
;
levnam = grb_file->lv_ISBL3
do ii = 0, dimsizes(levnam)-1
  if (levnam(ii) .eq. 500) then
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
  hvar = grb_file->GH_4_ISBL(level,:,:) / 10. ; use decameters
  uvar = grb_file->UW_4_ISBL(level,:,:) 
  vvar = grb_file->VW_4_ISBL(level,:,:)
  vort = new((/ny,nx/),float)
  flipuvar = new((/ny,nx/),float)
  flipvvar = new((/ny,nx/),float)
  fliplat = new((/ny/),float)
  do j = 0,ny-1
    do k = 0,nx-1
      fj = 360 - j
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
  ny = filevardimsizes(grb_file,"gridlat_83")
  nx = filevardimsizes(grb_file,"gridlon_83")
;
  intm = grb_file->T_83_SFC@initial_time
  fctm = grb_file->T_83_SFC@forecast_time
  lat = grb_file->gridlat_83(:,:)
  lon = grb_file->gridlon_83(:,:)
  hvar = grb_file->GH_83_ISBL(level,:,:) / 10. ; use decameters
  uvar = grb_file->UW_83_ISBL(level,:,:)
  vvar = grb_file->VW_83_ISBL(level,:,:)
  vort = new((/ny,nx/),float)
;
  uv2vrf(uvar,vvar,vort)
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
  rot = grb_file->gridrot_130(:,:) ; use to rotate winds
  latin1 = grb_file->gridlon_130@Latin1
  latin2 = grb_file->gridlon_130@Latin2
  lov = grb_file->gridlon_130@Lov
  hvar = grb_file->GH_130_ISBL(level,:,:) / 10. ; use decameters
;  uvar = grb_file->UW_130_ISBL(level,:,:)
;  vvar = grb_file->VW_130_ISBL(level,:,:)
  uvartmp = grb_file->UW_130_ISBL(level,:,:) ; use to rotate winds
  vvartmp = grb_file->VW_130_ISBL(level,:,:) ; use to rotate winds
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
;  rot = grb_file->gridrot_236(:,:) ; use to rotate winds
  latin1 = grb_file->gridlon_236@Latin1
  latin2 = grb_file->gridlon_236@Latin2
  lov = grb_file->gridlon_236@Lov
  hvar = grb_file->GH_236_ISBL(level,:,:) / 10. ; use decameters
  uvar = grb_file->UW_236_ISBL(level,:,:)
  vvar = grb_file->VW_236_ISBL(level,:,:)
;  uvartmp = grb_file->UW_236_ISBL(level,:,:) ; use to rotate winds
;  vvartmp = grb_file->VW_236_ISBL(level,:,:) ; use to rotate winds
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
if (domainname .eq. "africa" .or. domainname .eq. "europe" .or. \
    domainname .eq. "floating" .or. domainname .eq. "e_pacific") then
;
  nx = filevardimsizes(grb_file,"g0_lon_2")
  ny = filevardimsizes(grb_file,"g0_lat_1")
;
  intm = grb_file->T_GDS0_SFC@initial_time
  fctm = grb_file->T_GDS0_SFC@forecast_time
  lat = grb_file->g0_lat_1(:)
  lon = grb_file->g0_lon_2(:)
  hvar = grb_file->GH_GDS0_ISBL(level,:,:) / 10. ; use decameters
  uvar = grb_file->UW_GDS0_ISBL(level,:,:) 
  vvar = grb_file->VW_GDS0_ISBL(level,:,:)
  vort = new((/ny,nx/),float)
;
  vort = uv2vr_cfd(uvar,vvar,lat,lon,2)
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
  hvar = grb_file->GH_219_ISBL(level,:,:) / 10. ; use decameters
  uvar = grb_file->UW_219_ISBL(level,:,:)
  vvar = grb_file->VW_219_ISBL(level,:,:)
  vort = new((/ny,nx/),float)
  do j = 1,ny-2
    do k = 1,nx-2
      dx = 25400.0
      dy = 25400.0
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
  hvar = grb_file->GH_28_ISBL(level,:,:) / 10. ; use decameters
  uvar = grb_file->UW_28_ISBL(level,:,:) * 1.9429 ; convert to knots
  vvar = grb_file->VW_28_ISBL(level,:,:) * 1.9429 ; convert to knots
  vort = new((/ny,nx/),float)
;
  do j = 1,ny-2
    do k = 1,nx-2
      dx = 381000.
      dy = dx
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
  ny = filevardimsizes(grb_file,"lat_244")
  nx = filevardimsizes(grb_file,"lon_244")
;
  intm = grb_file->T_244_SFC@initial_time
  fctm = grb_file->T_244_SFC@forecast_time
  lat = grb_file->lat_244(:)
  lat = lat - 50.955
  lon = grb_file->lon_244(:)
  hvar = grb_file->GH_244_ISBL(level,:,:) / 10. ; use decameters
  uvar = grb_file->UW_244_ISBL(level,:,:)
  vvar = grb_file->VW_244_ISBL(level,:,:)
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
  hvar = grb_file->GH_GDS0_ISBL(level,:,:) / 10. ; use decameters
  uvar = grb_file->UW_GDS0_ISBL(level,:,:)
  vvar = grb_file->VW_GDS0_ISBL(level,:,:)
  vort = new((/ny,nx/),float)
;
  vort = uv2vr_cfd(uvar,vvar,lat,lon,2)
;
end if
;
load "$NCL_HOME/obj/get_date_str_FIM.ncl"
;
print(date_str)
vort = vort * 1.e05
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
           "wkMetaName" : "./ua_vort.ncgm"
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
           "wkPSFileName" : "./ua_vort.ps"
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
    "cnRasterModeOn" : False
  end setvalues
;
  load "$NCL_HOME/obj/height_60m.ncl"
;
  setvalues cnid4
    "cnLineColor" : 20
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
else  ; domain is not 130
;
  load "$NCL_HOME/obj/fim_maps_blk.ncl"
;
; Create ScalarField objects.
;
  if (domainname .eq. "fim") then
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
  if (domainname .eq. "236" .or. domainname .eq. "africa" .or. \
      domainname .eq. "219" .or. domainname .eq. "w_pacific" .or. \
      domainname .eq. "europe" .or. domainname .eq. "floating" .or. \
      domainname .eq. "e_pacific" .or. domainname .eq. "83") then
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
  if (domainname .eq. "28") then
;
    tfield = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : vort
        "sfXCStartV" : xdim
        "sfYCStartV" : ydim
        "sfXCEndV" : 0
        "sfYCEndV" : 0
        "sfMissingValueV" : vort@_FillValue
    end create
;
    hfield = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : hvar
        "sfXCStartV" : xdim
        "sfYCStartV" : ydim
        "sfXCEndV" : 0
        "sfYCEndV" : 0
        "sfMissingValueV" : hvar@_FillValue
    end create
;
  end if
;
; Create a ContourPlot object.
;
  load "$NCL_HOME/obj/diff_-16-16by2.ncl"
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
  load "$NCL_HOME/obj/height_60m.ncl"
;
  setvalues cnid4
    "cnLineColor" : 20
  end setvalues
;
  if (domainname .eq. "fim") then
    setvalues cnid4
      "cnLineThicknessF" : 1
      "cnLineLabelDensityF" : 0.001
    end setvalues
  end if
;
  if (domainname .eq. "219") then
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
end if  ; domainname .eq. "130"
;
end
;
