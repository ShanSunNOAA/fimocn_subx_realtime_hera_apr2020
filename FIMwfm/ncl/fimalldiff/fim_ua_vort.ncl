;
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_LIB/ncarg/nclscripts/csm/contributed.ncl"
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
;
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
  uvar_2 = grb_file2->UW_4_ISBL(level,:,:)
  vvar_2 = grb_file2->VW_4_ISBL(level,:,:)
  vort_2 = new((/ny,nx/),float)
  flipuvar_2 = new((/ny,nx/),float)
  flipvvar_2 = new((/ny,nx/),float)
  fliplat_2 = new((/ny/),float)
  do j = 0,ny-1
    do k = 0,nx-1
      fj = 360 - j
      fliplat_2(j) = lat(fj)
      flipuvar_2(j,k) = uvar_2(fj,k)
      flipvvar_2(j,k) = vvar_2(fj,k)
    end do
  end do
  lat = fliplat_2
  uvar_2 = flipuvar_2
  vvar_2 = flipvvar_2
;
  uv2vrf(uvar_2,vvar_2,vort_2)
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
  uvar = grb_file->UW_236_ISBL(level,:,:)
  vvar = grb_file->VW_236_ISBL(level,:,:)
  vort = new((/ny,nx/),float)
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
  uvar_2 = grb_file2->UW_236_ISBL(level,:,:)
  vvar_2 = grb_file2->VW_236_ISBL(level,:,:)
  vort_2 = new((/ny,nx/),float)
;
  do j = 1,ny-2
    do k = 1,nx-2
      dx = 40635.3
      dy = 40635.3
      jm1 = j - 1
      jp1 = j + 1
      km1 = k - 1
      kp1 = k + 1
      dv_2 = vvar_2(j,kp1) - vvar_2(j,km1)
      du_2 = uvar_2(jp1,k) - uvar_2(jm1,k)
      vort_2(j,k) = ((dv_2/(2.*dx)) - (du_2/(2.*dy)))
    end do
  end do
  do k = 0,nx-1
    vort_2(0,k) = vort_2(1,k)
    vort_2(ny-1,k) = vort_2(ny-2,k)
  end do
  do j = 0,ny-1
    vort_2(j,0) = vort_2(j,1)
    vort_2(j,nx-1) = vort_2(j,nx-2)
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
;
  uvar = grb_file->UW_GDS0_ISBL(level,:,:) 
  vvar = grb_file->VW_GDS0_ISBL(level,:,:)
  vort = new((/ny,nx/),float)
;
  vort = uv2vr_cfd(uvar,vvar,lat,lon,2)
;
  uvar_2 = grb_file2->UW_GDS0_ISBL(level,:,:)
  vvar_2 = grb_file2->VW_GDS0_ISBL(level,:,:)
  vort_2 = new((/ny,nx/),float)
;
  vort_2 = uv2vr_cfd(uvar_2,vvar_2,lat,lon,2)
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
  uvar = grb_file->UW_219_ISBL(level,:,:)
  vvar = grb_file->VW_219_ISBL(level,:,:)
  vort = new((/ny,nx/),float)
;
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
  uvar_2 = grb_file2->UW_219_ISBL(level,:,:)
  vvar_2 = grb_file2->VW_219_ISBL(level,:,:)
  vort_2 = new((/ny,nx/),float)
;
  do j = 1,ny-2
    do k = 1,nx-2
      dx = 25400.0
      dy = 25400.0
      jm1 = j - 1
      jp1 = j + 1
      km1 = k - 1
      kp1 = k + 1
      dv_2 = vvar_2(j,kp1) - vvar_2(j,km1)
      du_2 = uvar_2(jp1,k) - uvar_2(jm1,k)
      vort_2(j,k) = ((dv_2/(2.*dx)) - (du_2/(2.*dy)))
    end do
  end do
  do k = 0,nx-1
    vort_2(0,k) = vort_2(1,k)
    vort_2(ny-1,k) = vort_2(ny-2,k)
  end do
  do j = 0,ny-1
    vort_2(j,0) = vort_2(j,1)
    vort_2(j,nx-1) = vort_2(j,nx-2)
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
;
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
  uvar_2 = grb_file2->UW_28_ISBL(level,:,:) * 1.9429 ; convert to knots
  vvar_2 = grb_file2->VW_28_ISBL(level,:,:) * 1.9429 ; convert to knots
  vort_2 = new((/ny,nx/),float)
;
  do j = 1,ny-2
    do k = 1,nx-2
      dx = 381000.
      dy = dx
      jm1 = j - 1
      jp1 = j + 1
      km1 = k - 1
      kp1 = k + 1
      dv_2 = vvar_2(j,kp1) - vvar_2(j,km1)
      du_2 = uvar_2(jp1,k) - uvar_2(jm1,k)
      vort_2(j,k) = ((dv_2/(2.*dx)) - (du_2/(2.*dy)))
    end do
  end do
  do k = 0,nx-1
    vort_2(0,k) = vort_2(1,k)
    vort_2(ny-1,k) = vort_2(ny-2,k)
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
;
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
  lat_2 = grb_file->lat_244(:)
  lat_2 = lat_2 - 50.955
  lon_2 = grb_file->lon_244(:)
  uvar_2 = grb_file2->UW_244_ISBL(level,:,:)
  vvar_2 = grb_file2->VW_244_ISBL(level,:,:)
  vort_2 = new((/ny,nx/),float)
  flipuvar_2 = new((/ny,nx/),float)
  flipvvar_2 = new((/ny,nx/),float)
  fliplat_2 = new((/ny/),float)
  do j = 0,ny-1
    do k = 0,nx-1
      fj = (ny-1) - j
      fliplat_2(j) = lat_2(fj)
      flipuvar_2(j,k) = uvar_2(fj,k)
      flipvvar_2(j,k) = vvar_2(fj,k)
    end do
  end do
  lat_2 = fliplat_2
;
  vort_2 = uv2vr_cfd(uvar_2,vvar_2,lat_2,lon_2,2)
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
  uvar = grb_file->UW_GDS0_ISBL(level,:,:)
  vvar = grb_file->VW_GDS0_ISBL(level,:,:)
  vort = new((/ny,nx/),float)
;
  vort = uv2vr_cfd(uvar,vvar,lat,lon,2)
;
  uvar_2 = grb_file2->UW_GDS0_ISBL(level,:,:)
  vvar_2 = grb_file2->VW_GDS0_ISBL(level,:,:)
  vort_2 = new((/ny,nx/),float)
;
  vort_2 = uv2vr_cfd(uvar_2,vvar_2,lat,lon,2)
;
end if
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
vort = vort * 1.e05
vort_2 = vort_2 * 1.e05
vort = vort - vort_2
print("500 mb vort diff")
print(min(vort))
print(max(vort))
;
if (min(vort) .eq. 0 .and. max(vort) .eq. 0) then
  print("constant field of zero - exiting")
  exit
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
end if
;
if (domainname .eq. "236" .or. domainname .eq. "africa" .or. \
    domainname .eq. "219" .or. domainname .eq. "w_pacific" .or. \
    domainname .eq. "europe" .or. domainname .eq. "floating" .or. \
    domainname .eq. "e_pacific") then
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
end if
;
; Create a ContourPlot object.
;
; Relative Vorticity
;
load "$NCL_HOME/obj/diff_-16-16by2.ncl"
;
if (domainname .eq. "236" .or. domainname .eq. "africa" .or. \
    domainname .eq. "244" .or. domainname .eq. "w_pacific" .or. \
    domainname .eq. "europe" .or. domainname .eq. "28" .or. \
    domainname .eq. "floating" .or. domainname .eq. "e_pacific") then
;
  setvalues cnid
;    "cnRasterModeOn" : False
    "cnRasterSmoothingOn" : True
  end setvalues
;
end if
;
titlestring = "500mb Rel Vort (10E-5/s)"
load "$NCL_HOME/obj/add_titles.ncl"
;
; Draw map
;
draw(cnid)
draw(mpid)
frame(wid)
;
end
;
