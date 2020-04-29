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
; Get the 850-500mb level indices
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
if (modl .eq. "MPAS") then
  do ii = 0, dimsizes(levnam)-1
    if (levnam(ii) .eq. 850) then
      level850 = ii
    end if
;    if (levnam(ii) .eq. 825) then
;      level825 = ii
;    end if
;    if (levnam(ii) .eq. 800) then
;      level800 = ii
;    end if
;    if (levnam(ii) .eq. 775) then
;      level775 = ii
;    end if
;    if (levnam(ii) .eq. 750) then
;      level750 = ii
;    end if
;    if (levnam(ii) .eq. 725) then
;      level725 = ii
;    end if
    if (levnam(ii) .eq. 700) then
      level700 = ii
    end if
;    if (levnam(ii) .eq. 675) then
;      level675 = ii
;    end if
;    if (levnam(ii) .eq. 650) then
;      level650 = ii
;    end if
;    if (levnam(ii) .eq. 625) then
;      level625 = ii
;    end if
;    if (levnam(ii) .eq. 600) then
;      level600 = ii
;    end if
;    if (levnam(ii) .eq. 575) then
;      level575 = ii
;    end if
;    if (levnam(ii) .eq. 550) then
;      level550 = ii
;    end if
;    if (levnam(ii) .eq. 525) then
;      level525 = ii
;    end if
    if (levnam(ii) .eq. 500) then
      level500 = ii
    end if
  end do
else
  do ii = 0, dimsizes(levnam)-1
    if (levnam(ii) .eq. 850) then
      level850 = ii
    end if
;    if (levnam(ii) .eq. 825) then
;      level825 = ii
;    end if
;    if (levnam(ii) .eq. 800) then
;      level800 = ii
;    end if
;    if (levnam(ii) .eq. 775) then
;      level775 = ii
;    end if
;    if (levnam(ii) .eq. 750) then
;      level750 = ii
;    end if
;    if (levnam(ii) .eq. 725) then
;      level725 = ii
;    end if
    if (levnam(ii) .eq. 700) then
      level700 = ii
    end if
;    if (levnam(ii) .eq. 675) then
;      level675 = ii
;    end if
;    if (levnam(ii) .eq. 650) then
;      level650 = ii
;    end if
;    if (levnam(ii) .eq. 625) then
;      level625 = ii
;    end if
;    if (levnam(ii) .eq. 600) then
;      level600 = ii
;    end if
;    if (levnam(ii) .eq. 575) then
;      level575 = ii
;    end if
;    if (levnam(ii) .eq. 550) then
;      level550 = ii
;    end if
;    if (levnam(ii) .eq. 525) then
;      level525 = ii
;    end if
    if (levnam(ii) .eq. 500) then
      level500 = ii
    end if
  end do
end if
;
if (domainname .eq. "fim") then
;
  if (modl .eq. "MPAS") then
    ny = filevardimsizes(grb_file,"g0_lat_0")
    nx = filevardimsizes(grb_file,"g0_lon_1")
;
    intm = grb_file->PRES_GDS0_SFC@initial_time
    fctm = grb_file->PRES_GDS0_SFC@forecast_time
    lat = grb_file->g0_lat_0(:)
    lon = grb_file->g0_lon_1(:)
    wstrid = 10
    uvar = grb_file->U_GRD_GDS0_ISBL(level700,::-1,:) * 1.9429 ; convert to knots
    vvar = grb_file->V_GRD_GDS0_ISBL(level700,::-1,:) * 1.9429 ; convert to knots
    rh85 = grb_file->R_H_GDS0_ISBL(level850,::-1,:) * 100.
    rh70 = grb_file->R_H_GDS0_ISBL(level700,::-1,:) * 100.
    rh60 = grb_file->R_H_GDS0_ISBL(level600,::-1,:) * 100.
    rh50 = grb_file->R_H_GDS0_ISBL(level500,::-1,:) * 100.
    tvar = (rh85 + rh70 + rh60 + rh50) / 4.
;
    uvar = lonFlip(uvar)
    vvar = lonFlip(vvar)
    tvar = lonFlip(tvar)
  else
    ny = filevardimsizes(grb_file,LAT)
    nx = filevardimsizes(grb_file,LON)
;
    lat = grb_file->$LAT$(:)
    lon = grb_file->$LON$(:)
    wstrid = 10
;
    uvar = grb_file->$UISBL$(level700,:,:) * 1.9429 ; convert to knots
    vvar = grb_file->$VISBL$(level700,:,:) * 1.9429 ; convert to knots
    rh85 = grb_file->$RHISBL$(level850,:,:)
;    rh82 = grb_file->$RHISBL$(level825,:,:)
;    rh80 = grb_file->$RHISBL$(level800,:,:)
;    rh77 = grb_file->$RHISBL$(level775,:,:)
;    rh75 = grb_file->$RHISBL$(level750,:,:)
;    rh72 = grb_file->$RHISBL$(level725,:,:)
    rh70 = grb_file->$RHISBL$(level700,:,:)
;    rh67 = grb_file->$RHISBL$(level675,:,:)
;    rh65 = grb_file->$RHISBL$(level650,:,:)
;    rh62 = grb_file->$RHISBL$(level625,:,:)
;    rh60 = grb_file->$RHISBL$(level600,:,:)
;    rh57 = grb_file->$RHISBL$(level575,:,:)
;    rh55 = grb_file->$RHISBL$(level550,:,:)
;    rh52 = grb_file->$RHISBL$(level525,:,:)
    rh50 = grb_file->$RHISBL$(level500,:,:)
;    tvar = (rh85 + rh82 + rh80 + rh77 + rh75 + rh72 + rh70 + \
;            rh67 + rh65 + rh62 + rh60 + rh57 + rh55 + rh52 + rh50) / 15
    tvar = (rh85 + rh70 + rh50) / 3
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
  uvar = grb_file->$UISBL$(level700,:,:) * 1.9429 ; convert to knots
  vvar = grb_file->$VISBL$(level700,:,:) * 1.9429 ; convert to knots
  rh85 = grb_file->$RHISBL$(level850,:,:)
;  rh82 = grb_file->$RHISBL$(level825,:,:)
;  rh80 = grb_file->$RHISBL$(level800,:,:)
;  rh77 = grb_file->$RHISBL$(level775,:,:)
;  rh75 = grb_file->$RHISBL$(level750,:,:)
;  rh72 = grb_file->$RHISBL$(level725,:,:)
  rh70 = grb_file->$RHISBL$(level700,:,:)
;  rh67 = grb_file->$RHISBL$(level675,:,:)
;  rh65 = grb_file->$RHISBL$(level650,:,:)
;  rh62 = grb_file->$RHISBL$(level625,:,:)
;  rh60 = grb_file->$RHISBL$(level600,:,:)
;  rh57 = grb_file->$RHISBL$(level575,:,:)
;  rh55 = grb_file->$RHISBL$(level550,:,:)
;  rh52 = grb_file->$RHISBL$(level525,:,:)
  rh50 = grb_file->$RHISBL$(level500,:,:)
;  tvar = (rh85 + rh82 + rh80 + rh77 + rh75 + rh72 + rh70 + \
;          rh67 + rh65 + rh62 + rh60 + rh57 + rh55 + rh52 + rh50) / 15
  tvar = (rh85 + rh70 + rh50) / 3
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
  uvar = grb_file->$UISBL$(level700,:,:) * 1.9429 ; convert to knots
  vvar = grb_file->$VISBL$(level700,:,:) * 1.9429 ; convert to knots
  rh85 = grb_file->$RHISBL$(level850,:,:)
;  rh82 = grb_file->$RHISBL$(level825,:,:)
;  rh80 = grb_file->$RHISBL$(level800,:,:)
;  rh77 = grb_file->$RHISBL$(level775,:,:)
;  rh75 = grb_file->$RHISBL$(level750,:,:)
;  rh72 = grb_file->$RHISBL$(level725,:,:)
  rh70 = grb_file->$RHISBL$(level700,:,:)
;  rh67 = grb_file->$RHISBL$(level675,:,:)
;  rh65 = grb_file->$RHISBL$(level650,:,:)
;  rh62 = grb_file->$RHISBL$(level625,:,:)
;  rh60 = grb_file->$RHISBL$(level600,:,:)
;  rh57 = grb_file->$RHISBL$(level575,:,:)
;  rh55 = grb_file->$RHISBL$(level550,:,:)
;  rh52 = grb_file->$RHISBL$(level525,:,:)
  rh50 = grb_file->$RHISBL$(level500,:,:)
;  tvar = (rh85 + rh82 + rh80 + rh77 + rh75 + rh72 + rh70 + \
;          rh67 + rh65 + rh62 + rh60 + rh57 + rh55 + rh52 + rh50) / 15
  tvar = (rh85 + rh70 + rh50) / 3
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
  wstrid = 6
;
  uvar = grb_file->$UISBL$(level700,:,:) * 1.9429 ; convert to knots
  vvar = grb_file->$VISBL$(level700,:,:) * 1.9429 ; convert to knots
  rh85 = grb_file->$RHISBL$(level850,:,:)
;  rh82 = grb_file->$RHISBL$(level825,:,:)
;  rh80 = grb_file->$RHISBL$(level800,:,:)
;  rh77 = grb_file->$RHISBL$(level775,:,:)
;  rh75 = grb_file->$RHISBL$(level750,:,:)
;  rh72 = grb_file->$RHISBL$(level725,:,:)
  rh70 = grb_file->$RHISBL$(level700,:,:)
;  rh67 = grb_file->$RHISBL$(level675,:,:)
;  rh65 = grb_file->$RHISBL$(level650,:,:)
;  rh62 = grb_file->$RHISBL$(level625,:,:)
;  rh60 = grb_file->$RHISBL$(level600,:,:)
;  rh57 = grb_file->$RHISBL$(level575,:,:)
;  rh55 = grb_file->$RHISBL$(level550,:,:)
;  rh52 = grb_file->$RHISBL$(level525,:,:)
  rh50 = grb_file->$RHISBL$(level500,:,:)
;  tvar = (rh85 + rh82 + rh80 + rh77 + rh75 + rh72 + rh70 + \
;          rh67 + rh65 + rh62 + rh60 + rh57 + rh55 + rh52 + rh50) / 15
  tvar = (rh85 + rh70 + rh50) / 3
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
  uvar = grb_file->$UISBL$(level700,::-1,:) * 1.9429 ; convert to knots
  uvarorig = uvar
  flipuvar = lonFlip(uvar)
  vvar = grb_file->$VISBL$(level700,::-1,:) * 1.9429 ; convert to knots
  vvarorig = vvar
  flipvvar = lonFlip(vvar)
  rh85 = grb_file->$RHISBL$(level850,::-1,:)
;  rh82 = grb_file->$RHISBL$(level825,::-1,:)
;  rh80 = grb_file->$RHISBL$(level800,::-1,:)
;  rh77 = grb_file->$RHISBL$(level775,::-1,:)
;  rh75 = grb_file->$RHISBL$(level750,::-1,:)
;  rh72 = grb_file->$RHISBL$(level725,::-1,:)
  rh70 = grb_file->$RHISBL$(level700,::-1,:)
;  rh67 = grb_file->$RHISBL$(level675,::-1,:)
;  rh65 = grb_file->$RHISBL$(level650,::-1,:)
;  rh62 = grb_file->$RHISBL$(level625,::-1,:)
;  rh60 = grb_file->$RHISBL$(level600,::-1,:)
;  rh57 = grb_file->$RHISBL$(level575,::-1,:)
;  rh55 = grb_file->$RHISBL$(level550,::-1,:)
;  rh52 = grb_file->$RHISBL$(level525,::-1,:)
  rh50 = grb_file->$RHISBL$(level500,::-1,:)
;  tvar = (rh85 + rh82 + rh80 + rh77 + rh75 + rh72 + rh70 + \
;          rh67 + rh65 + rh62 + rh60 + rh57 + rh55 + rh52 + rh50) / 15
  tvar = (rh85 + rh70 + rh50) / 3
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
  wstrid = 3
;
  uvar = grb_file->$UISBL$(level700,:,:) * 1.9429 ; convert to knots
  vvar = grb_file->$VISBL$(level700,:,:) * 1.9429 ; convert to knots
  rh85 = grb_file->$RHISBL$(level850,:,:)
;  rh82 = grb_file->$RHISBL$(level825,:,:)
;  rh80 = grb_file->$RHISBL$(level800,:,:)
;  rh77 = grb_file->$RHISBL$(level775,:,:)
;  rh75 = grb_file->$RHISBL$(level750,:,:)
;  rh72 = grb_file->$RHISBL$(level725,:,:)
  rh70 = grb_file->$RHISBL$(level700,:,:)
;  rh67 = grb_file->$RHISBL$(level675,:,:)
;  rh65 = grb_file->$RHISBL$(level650,:,:)
;  rh62 = grb_file->$RHISBL$(level625,:,:)
;  rh60 = grb_file->$RHISBL$(level600,:,:)
;  rh57 = grb_file->$RHISBL$(level575,:,:)
;  rh55 = grb_file->$RHISBL$(level550,:,:)
;  rh52 = grb_file->$RHISBL$(level525,:,:)
  rh50 = grb_file->$RHISBL$(level500,:,:)
;  tvar = (rh85 + rh82 + rh80 + rh77 + rh75 + rh72 + rh70 + \
;          rh67 + rh65 + rh62 + rh60 + rh57 + rh55 + rh52 + rh50) / 15
  tvar = (rh85 + rh70 + rh50) / 3
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
  wstrid = 4
;
  uvartmp = grb_file->$UISBL$(level700,:,:)
  vvartmp = grb_file->$VISBL$(level700,:,:)
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
  rh85 = grb_file->$RHISBL$(level850,:,:)
;  rh82 = grb_file->$RHISBL$(level825,:,:)
;  rh80 = grb_file->$RHISBL$(level800,:,:)
;  rh77 = grb_file->$RHISBL$(level775,:,:)
;  rh75 = grb_file->$RHISBL$(level750,:,:)
;  rh72 = grb_file->$RHISBL$(level725,:,:)
  rh70 = grb_file->$RHISBL$(level700,:,:)
;  rh67 = grb_file->$RHISBL$(level675,:,:)
;  rh65 = grb_file->$RHISBL$(level650,:,:)
;  rh62 = grb_file->$RHISBL$(level625,:,:)
;  rh60 = grb_file->$RHISBL$(level600,:,:)
;  rh57 = grb_file->$RHISBL$(level575,:,:)
;  rh55 = grb_file->$RHISBL$(level550,:,:)
;  rh52 = grb_file->$RHISBL$(level525,:,:)
  rh50 = grb_file->$RHISBL$(level500,:,:)
;  tvar = (rh85 + rh82 + rh80 + rh77 + rh75 + rh72 + rh70 + \
;          rh67 + rh65 + rh62 + rh60 + rh57 + rh55 + rh52 + rh50) / 15
  tvar = (rh85 + rh70 + rh50) / 3
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
  wstrid = 7
;
  uvar = grb_file->$UISBL$(level700,:,:) * 1.9429 ; convert to knots
  vvar = grb_file->$VISBL$(level700,:,:) * 1.9429 ; convert to knots
  rh85 = grb_file->$RHISBL$(level850,:,:)
;  rh82 = grb_file->$RHISBL$(level825,:,:)
;  rh80 = grb_file->$RHISBL$(level800,:,:)
;  rh77 = grb_file->$RHISBL$(level775,:,:)
;  rh75 = grb_file->$RHISBL$(level750,:,:)
;  rh72 = grb_file->$RHISBL$(level725,:,:)
  rh70 = grb_file->$RHISBL$(level700,:,:)
;  rh67 = grb_file->$RHISBL$(level675,:,:)
;  rh65 = grb_file->$RHISBL$(level650,:,:)
;  rh62 = grb_file->$RHISBL$(level625,:,:)
;  rh60 = grb_file->$RHISBL$(level600,:,:)
;  rh57 = grb_file->$RHISBL$(level575,:,:)
;  rh55 = grb_file->$RHISBL$(level550,:,:)
;  rh52 = grb_file->$RHISBL$(level525,:,:)
  rh50 = grb_file->$RHISBL$(level500,:,:)
;  tvar = (rh85 + rh82 + rh80 + rh77 + rh75 + rh72 + rh70 + \
;          rh67 + rh65 + rh62 + rh60 + rh57 + rh55 + rh52 + rh50) / 15
  tvar = (rh85 + rh70 + rh50) / 3
;
end if
;
load "$NCL_HOME/obj/get_date_str_FIM.ncl"
;
print(date_str)
print("Mean 850-500 mb RH")
print(min(tvar))
print(max(tvar))
;
if (min(tvar) .eq. max(tvar)) then
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
           "wkMetaName" : "./ua_rh.ncgm"
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
           "wkPSFileName" : "./ua_rh.ps"
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
    setvalues mpid
      "mpUSStateLineColor" : "background"
      "mpUSStateLineThicknessF" : 2
      "mpNationalLineColor" : "background"
      "mpNationalLineThicknessF" : 2
      "mpGeophysicalLineColor" : "background"
      "mpGeophysicalLineThicknessF" : 2
    end setvalues
  else
    load "$NCL_HOME/obj/lc_map_counties.ncl"
    setvalues mpid
      "mpUSStateLineColor" : "background"
      "mpUSStateLineThicknessF" : 2
      "mpNationalLineColor" : "foreground"
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
;
  load "$NCL_HOME/obj/rh.ncl"
;
  setvalues cnid
    "cnRasterModeOn" : False
  end setvalues
;
  titlestring = "850-500mb mean Rel Hum (%), 700mb Wind (kt)"
  load "$NCL_HOME/obj/add_titles.ncl"
  load "$NCL_HOME/obj/add_explabel_vcid.ncl"
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
  draw(vcid)
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
      "sfDataArray" : tvar
      "sfXCStartSubsetV" : llx
      "sfYCStartSubsetV" : lly
      "sfXCEndSubsetV" : urx
      "sfYCEndSubsetV" : ury
      "sfMissingValueV" : tvar@_FillValue
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
  load "$NCL_HOME/obj/rh.ncl"
;
  titlestring = "850-500mb mean Rel Hum (%), 700mb Wind (kt)"
  load "$NCL_HOME/obj/add_titles.ncl"
  load "$NCL_HOME/obj/add_explabel_vcid.ncl"
;
  draw(cnid)
  draw(mpid)
  draw(vcid)
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
      "sfDataArray" : tvar
      "sfXCStartSubsetV" : llx
      "sfYCStartSubsetV" : lly
      "sfXCEndSubsetV" : urx
      "sfYCEndSubsetV" : ury
      "sfMissingValueV" : tvar@_FillValue
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
  load "$NCL_HOME/obj/rh.ncl"
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
  titlestring = "850-500mb mean Rel Hum (%), 700mb Wind (kt)"
  load "$NCL_HOME/obj/add_titles.ncl"
  load "$NCL_HOME/obj/add_explabel_vcid.ncl"
;
  draw(cnid)
  draw(mpid)
  draw(vcid)
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
; Create ScalarField and VectorField objects.
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
        "sfDataArray" : tvar
        "sfXCStartV" : 0
        "sfYCStartV" : 0
        "sfXCEndV" : xdim
        "sfYCEndV" : ydim
        "sfMissingValueV" : tvar@_FillValue
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
        "vfXCStride" : wstrid      ; i.e., plot a vector at every (n+1)th point.
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
  load "$NCL_HOME/obj/rh.ncl"
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
  titlestring = "850-500mb mean Rel Hum (%), 700mb Wind (kt)"
  load "$NCL_HOME/obj/add_titles.ncl"
  load "$NCL_HOME/obj/add_explabel_vcid.ncl"
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
  draw(vcid)
  frame(wid)
;
end if  ; domainname .ne. "130" .and. domainname .ne. "174" .and. domainname .ne. "242"
;
end
;
