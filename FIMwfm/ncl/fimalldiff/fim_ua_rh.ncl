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
fname2 = asciiread("arw_file2.txt",-1,"string")
grb_file2 = addfile(fname2,"r")
;
domainname = asciiread("domain.txt",-1,"string")
;
load "$NCL_HOME/obj/fim_grib_names.ncl"
;
; Get all level indexes between 850mb and 500mb
;
do ll = 0,5
  lnam = "lv_ISBL"+ll
  if (isfilevar(grb_file,lnam)) then
    break
  end if
end do
;
lcount = 0
levnam = grb_file->$lnam$
levunits = grb_file->$lnam$@units
if (levunits .eq. "Pa") then
  levnam = levnam / 100
end if
do ii = 0, dimsizes(levnam)-1
  if (levnam(ii) .eq. 850) then
    level850 = ii
    lcount = lcount + 1
  end if
;  if (levnam(ii) .eq. 825) then
;    level825 = ii
;    lcount = lcount + 1
;  end if
;  if (levnam(ii) .eq. 800) then
;    level800 = ii
;    lcount = lcount + 1
;  end if
;  if (levnam(ii) .eq. 775) then
;    level775 = ii
;    lcount = lcount + 1
;  end if
;  if (levnam(ii) .eq. 750) then
;    level750 = ii
;    lcount = lcount + 1
;  end if
;  if (levnam(ii) .eq. 725) then
;    level725 = ii
;    lcount = lcount + 1
;  end if
  if (levnam(ii) .eq. 700) then
    level700 = ii
    lcount = lcount + 1
  end if
;  if (levnam(ii) .eq. 675) then
;    level675 = ii
;    lcount = lcount + 1
;  end if
;  if (levnam(ii) .eq. 650) then
;    level650 = ii
;    lcount = lcount + 1
;  end if
;  if (levnam(ii) .eq. 625) then
;    level625 = ii
;    lcount = lcount + 1
;  end if
;  if (levnam(ii) .eq. 600) then
;    level600 = ii
;    lcount = lcount + 1
;  end if
;  if (levnam(ii) .eq. 575) then
;    level575 = ii
;    lcount = lcount + 1
;  end if
;  if (levnam(ii) .eq. 550) then
;    level550 = ii
;    lcount = lcount + 1
;  end if
;  if (levnam(ii) .eq. 525) then
;    level525 = ii
;    lcount = lcount + 1
;  end if
  if (levnam(ii) .eq. 500) then
    level500 = ii
    lcount = lcount + 1
  end if
end do
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
  rhsum = new((/ny,nx/),float)
  rhsum = 0.
  rhsum2 = new((/ny,nx/),float)
  rhsum2 = 0.
  if (isvar("level850")) then
    rh85 = grb_file->$RHISBL$(level850,:,:)
    rhsum = rhsum + rh85
    rh85_2 = grb_file2->$RHISBL$(level850,:,:)
    rhsum2 = rhsum2 + rh85_2
  end if
  if (isvar("level825")) then
    rh82 = grb_file->$RHISBL$(level825,:,:)
    rhsum = rhsum + rh82
    rh82_2 = grb_file2->$RHISBL$(level825,:,:)
    rhsum2 = rhsum2 + rh82_2
  end if
  if (isvar("level800")) then
    rh80 = grb_file->$RHISBL$(level800,:,:)
    rhsum = rhsum + rh80
    rh80_2 = grb_file2->$RHISBL$(level800,:,:)
    rhsum2 = rhsum2 + rh80_2
  end if
  if (isvar("level775")) then
    rh77 = grb_file->$RHISBL$(level775,:,:)
    rhsum = rhsum + rh77
    rh77_2 = grb_file2->$RHISBL$(level775,:,:)
    rhsum2 = rhsum2 + rh77_2
  end if
  if (isvar("level750")) then
    rh75 = grb_file->$RHISBL$(level750,:,:)
    rhsum = rhsum + rh75
    rh75_2 = grb_file2->$RHISBL$(level750,:,:)
    rhsum2 = rhsum2 + rh75_2
  end if
  if (isvar("level725")) then
    rh72 = grb_file->$RHISBL$(level725,:,:)
    rhsum = rhsum + rh72
    rh72_2 = grb_file2->$RHISBL$(level725,:,:)
    rhsum2 = rhsum2 + rh72_2
  end if
  if (isvar("level700")) then
    rh70 = grb_file->$RHISBL$(level700,:,:)
    rhsum = rhsum + rh70
    rh70_2 = grb_file2->$RHISBL$(level700,:,:)
    rhsum2 = rhsum2 + rh70_2
  end if
  if (isvar("level675")) then
    rh67 = grb_file->$RHISBL$(level675,:,:)
    rhsum = rhsum + rh67
    rh67_2 = grb_file2->$RHISBL$(level675,:,:)
    rhsum2 = rhsum2 + rh67_2
  end if
  if (isvar("level650")) then
    rh65 = grb_file->$RHISBL$(level650,:,:)
    rhsum = rhsum + rh65
    rh65_2 = grb_file2->$RHISBL$(level650,:,:)
    rhsum2 = rhsum2 + rh65_2
  end if
  if (isvar("level625")) then
    rh62 = grb_file->$RHISBL$(level625,:,:)
    rhsum = rhsum + rh62
    rh62_2 = grb_file2->$RHISBL$(level625,:,:)
    rhsum2 = rhsum2 + rh62_2
  end if
  if (isvar("level600")) then
    rh60 = grb_file->$RHISBL$(level600,:,:)
    rhsum = rhsum + rh60
    rh60_2 = grb_file2->$RHISBL$(level600,:,:)
    rhsum2 = rhsum2 + rh60_2
  end if
  if (isvar("level575")) then
    rh57 = grb_file->$RHISBL$(level575,:,:)
    rhsum = rhsum + rh57
    rh57_2 = grb_file2->$RHISBL$(level575,:,:)
    rhsum2 = rhsum2 + rh57_2
  end if
  if (isvar("level550")) then
    rh55 = grb_file->$RHISBL$(level550,:,:)
    rhsum = rhsum + rh55
    rh55_2 = grb_file2->$RHISBL$(level550,:,:)
    rhsum2 = rhsum2 + rh55_2
  end if
  if (isvar("level525")) then
    rh52 = grb_file->$RHISBL$(level525,:,:)
    rhsum = rhsum + rh52
    rh52_2 = grb_file2->$RHISBL$(level525,:,:)
    rhsum2 = rhsum2 + rh52_2
  end if
  if (isvar("level500")) then
    rh50 = grb_file->$RHISBL$(level500,:,:)
    rhsum = rhsum + rh50
    rh50_2 = grb_file2->$RHISBL$(level500,:,:)
    rhsum2 = rhsum2 + rh50_2
  end if
;
  tvar = rhsum / lcount
;
  tvar2 = rhsum2 / lcount
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
  rhsum = new((/ny,nx/),float)
  rhsum = 0.
  rhsum2 = new((/ny,nx/),float)
  rhsum2 = 0.
  if (isvar("level850")) then
    rh85 = grb_file->$RHISBL$(level850,:,:)
    rhsum = rhsum + rh85
    rh85_2 = grb_file2->$RHISBL$(level850,:,:)
    rhsum2 = rhsum2 + rh85_2
  end if
  if (isvar("level825")) then
    rh82 = grb_file->$RHISBL$(level825,:,:)
    rhsum = rhsum + rh82
    rh82_2 = grb_file2->$RHISBL$(level825,:,:)
    rhsum2 = rhsum2 + rh82_2
  end if
  if (isvar("level800")) then
    rh80 = grb_file->$RHISBL$(level800,:,:)
    rhsum = rhsum + rh80
    rh80_2 = grb_file2->$RHISBL$(level800,:,:)
    rhsum2 = rhsum2 + rh80_2
  end if
  if (isvar("level775")) then
    rh77 = grb_file->$RHISBL$(level775,:,:)
    rhsum = rhsum + rh77
    rh77_2 = grb_file2->$RHISBL$(level775,:,:)
    rhsum2 = rhsum2 + rh77_2
  end if
  if (isvar("level750")) then
    rh75 = grb_file->$RHISBL$(level750,:,:)
    rhsum = rhsum + rh75
    rh75_2 = grb_file2->$RHISBL$(level750,:,:)
    rhsum2 = rhsum2 + rh75_2
  end if
  if (isvar("level725")) then
    rh72 = grb_file->$RHISBL$(level725,:,:)
    rhsum = rhsum + rh72
    rh72_2 = grb_file2->$RHISBL$(level725,:,:)
    rhsum2 = rhsum2 + rh72_2
  end if
  if (isvar("level700")) then
    rh70 = grb_file->$RHISBL$(level700,:,:)
    rhsum = rhsum + rh70
    rh70_2 = grb_file2->$RHISBL$(level700,:,:)
    rhsum2 = rhsum2 + rh70_2
  end if
  if (isvar("level675")) then
    rh67 = grb_file->$RHISBL$(level675,:,:)
    rhsum = rhsum + rh67
    rh67_2 = grb_file2->$RHISBL$(level675,:,:)
    rhsum2 = rhsum2 + rh67_2
  end if
  if (isvar("level650")) then
    rh65 = grb_file->$RHISBL$(level650,:,:)
    rhsum = rhsum + rh65
    rh65_2 = grb_file2->$RHISBL$(level650,:,:)
    rhsum2 = rhsum2 + rh65_2
  end if
  if (isvar("level625")) then
    rh62 = grb_file->$RHISBL$(level625,:,:)
    rhsum = rhsum + rh62
    rh62_2 = grb_file2->$RHISBL$(level625,:,:)
    rhsum2 = rhsum2 + rh62_2
  end if
  if (isvar("level600")) then
    rh60 = grb_file->$RHISBL$(level600,:,:)
    rhsum = rhsum + rh60
    rh60_2 = grb_file2->$RHISBL$(level600,:,:)
    rhsum2 = rhsum2 + rh60_2
  end if
  if (isvar("level575")) then
    rh57 = grb_file->$RHISBL$(level575,:,:)
    rhsum = rhsum + rh57
    rh57_2 = grb_file2->$RHISBL$(level575,:,:)
    rhsum2 = rhsum2 + rh57_2
  end if
  if (isvar("level550")) then
    rh55 = grb_file->$RHISBL$(level550,:,:)
    rhsum = rhsum + rh55
    rh55_2 = grb_file2->$RHISBL$(level550,:,:)
    rhsum2 = rhsum2 + rh55_2
  end if
  if (isvar("level525")) then
    rh52 = grb_file->$RHISBL$(level525,:,:)
    rhsum = rhsum + rh52
    rh52_2 = grb_file2->$RHISBL$(level525,:,:)
    rhsum2 = rhsum2 + rh52_2
  end if
  if (isvar("level500")) then
    rh50 = grb_file->$RHISBL$(level500,:,:)
    rhsum = rhsum + rh50
    rh50_2 = grb_file2->$RHISBL$(level500,:,:)
    rhsum2 = rhsum2 + rh50_2
  end if
;
  tvar = rhsum / lcount
;
  tvar2 = rhsum2 / lcount
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
  rhsum = new((/ny,nx/),float)
  rhsum = 0.
  rhsum2 = new((/ny,nx/),float)
  rhsum2 = 0.
  if (isvar("level850")) then
    rh85 = grb_file->$RHISBL$(level850,:,:)
    rhsum = rhsum + rh85
    rh85_2 = grb_file2->$RHISBL$(level850,:,:)
    rhsum2 = rhsum2 + rh85_2
  end if
  if (isvar("level825")) then
    rh82 = grb_file->$RHISBL$(level825,:,:)
    rhsum = rhsum + rh82
    rh82_2 = grb_file2->$RHISBL$(level825,:,:)
    rhsum2 = rhsum2 + rh82_2
  end if
  if (isvar("level800")) then
    rh80 = grb_file->$RHISBL$(level800,:,:)
    rhsum = rhsum + rh80
    rh80_2 = grb_file2->$RHISBL$(level800,:,:)
    rhsum2 = rhsum2 + rh80_2
  end if
  if (isvar("level775")) then
    rh77 = grb_file->$RHISBL$(level775,:,:)
    rhsum = rhsum + rh77
    rh77_2 = grb_file2->$RHISBL$(level775,:,:)
    rhsum2 = rhsum2 + rh77_2
  end if
  if (isvar("level750")) then
    rh75 = grb_file->$RHISBL$(level750,:,:)
    rhsum = rhsum + rh75
    rh75_2 = grb_file2->$RHISBL$(level750,:,:)
    rhsum2 = rhsum2 + rh75_2
  end if
  if (isvar("level725")) then
    rh72 = grb_file->$RHISBL$(level725,:,:)
    rhsum = rhsum + rh72
    rh72_2 = grb_file2->$RHISBL$(level725,:,:)
    rhsum2 = rhsum2 + rh72_2
  end if
  if (isvar("level700")) then
    rh70 = grb_file->$RHISBL$(level700,:,:)
    rhsum = rhsum + rh70
    rh70_2 = grb_file2->$RHISBL$(level700,:,:)
    rhsum2 = rhsum2 + rh70_2
  end if
  if (isvar("level675")) then
    rh67 = grb_file->$RHISBL$(level675,:,:)
    rhsum = rhsum + rh67
    rh67_2 = grb_file2->$RHISBL$(level675,:,:)
    rhsum2 = rhsum2 + rh67_2
  end if
  if (isvar("level650")) then
    rh65 = grb_file->$RHISBL$(level650,:,:)
    rhsum = rhsum + rh65
    rh65_2 = grb_file2->$RHISBL$(level650,:,:)
    rhsum2 = rhsum2 + rh65_2
  end if
  if (isvar("level625")) then
    rh62 = grb_file->$RHISBL$(level625,:,:)
    rhsum = rhsum + rh62
    rh62_2 = grb_file2->$RHISBL$(level625,:,:)
    rhsum2 = rhsum2 + rh62_2
  end if
  if (isvar("level600")) then
    rh60 = grb_file->$RHISBL$(level600,:,:)
    rhsum = rhsum + rh60
    rh60_2 = grb_file2->$RHISBL$(level600,:,:)
    rhsum2 = rhsum2 + rh60_2
  end if
  if (isvar("level575")) then
    rh57 = grb_file->$RHISBL$(level575,:,:)
    rhsum = rhsum + rh57
    rh57_2 = grb_file2->$RHISBL$(level575,:,:)
    rhsum2 = rhsum2 + rh57_2
  end if
  if (isvar("level550")) then
    rh55 = grb_file->$RHISBL$(level550,:,:)
    rhsum = rhsum + rh55
    rh55_2 = grb_file2->$RHISBL$(level550,:,:)
    rhsum2 = rhsum2 + rh55_2
  end if
  if (isvar("level525")) then
    rh52 = grb_file->$RHISBL$(level525,:,:)
    rhsum = rhsum + rh52
    rh52_2 = grb_file2->$RHISBL$(level525,:,:)
    rhsum2 = rhsum2 + rh52_2
  end if
  if (isvar("level500")) then
    rh50 = grb_file->$RHISBL$(level500,:,:)
    rhsum = rhsum + rh50
    rh50_2 = grb_file2->$RHISBL$(level500,:,:)
    rhsum2 = rhsum2 + rh50_2
  end if
;
  tvar = rhsum / lcount
;
  tvar2 = rhsum2 / lcount
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
  rhsum = new((/ny,nx/),float)
  rhsum = 0.
  rhsum2 = new((/ny,nx/),float)
  rhsum2 = 0.
  if (isvar("level850")) then
    rh85 = grb_file->$RHISBL$(level850,:,:)
    rhsum = rhsum + rh85
    rh85_2 = grb_file2->$RHISBL$(level850,:,:)
    rhsum2 = rhsum2 + rh85_2
  end if
  if (isvar("level825")) then
    rh82 = grb_file->$RHISBL$(level825,:,:)
    rhsum = rhsum + rh82
    rh82_2 = grb_file2->$RHISBL$(level825,:,:)
    rhsum2 = rhsum2 + rh82_2
  end if
  if (isvar("level800")) then
    rh80 = grb_file->$RHISBL$(level800,:,:)
    rhsum = rhsum + rh80
    rh80_2 = grb_file2->$RHISBL$(level800,:,:)
    rhsum2 = rhsum2 + rh80_2
  end if
  if (isvar("level775")) then
    rh77 = grb_file->$RHISBL$(level775,:,:)
    rhsum = rhsum + rh77
    rh77_2 = grb_file2->$RHISBL$(level775,:,:)
    rhsum2 = rhsum2 + rh77_2
  end if
  if (isvar("level750")) then
    rh75 = grb_file->$RHISBL$(level750,:,:)
    rhsum = rhsum + rh75
    rh75_2 = grb_file2->$RHISBL$(level750,:,:)
    rhsum2 = rhsum2 + rh75_2
  end if
  if (isvar("level725")) then
    rh72 = grb_file->$RHISBL$(level725,:,:)
    rhsum = rhsum + rh72
    rh72_2 = grb_file2->$RHISBL$(level725,:,:)
    rhsum2 = rhsum2 + rh72_2
  end if
  if (isvar("level700")) then
    rh70 = grb_file->$RHISBL$(level700,:,:)
    rhsum = rhsum + rh70
    rh70_2 = grb_file2->$RHISBL$(level700,:,:)
    rhsum2 = rhsum2 + rh70_2
  end if
  if (isvar("level675")) then
    rh67 = grb_file->$RHISBL$(level675,:,:)
    rhsum = rhsum + rh67
    rh67_2 = grb_file2->$RHISBL$(level675,:,:)
    rhsum2 = rhsum2 + rh67_2
  end if
  if (isvar("level650")) then
    rh65 = grb_file->$RHISBL$(level650,:,:)
    rhsum = rhsum + rh65
    rh65_2 = grb_file2->$RHISBL$(level650,:,:)
    rhsum2 = rhsum2 + rh65_2
  end if
  if (isvar("level625")) then
    rh62 = grb_file->$RHISBL$(level625,:,:)
    rhsum = rhsum + rh62
    rh62_2 = grb_file2->$RHISBL$(level625,:,:)
    rhsum2 = rhsum2 + rh62_2
  end if
  if (isvar("level600")) then
    rh60 = grb_file->$RHISBL$(level600,:,:)
    rhsum = rhsum + rh60
    rh60_2 = grb_file2->$RHISBL$(level600,:,:)
    rhsum2 = rhsum2 + rh60_2
  end if
  if (isvar("level575")) then
    rh57 = grb_file->$RHISBL$(level575,:,:)
    rhsum = rhsum + rh57
    rh57_2 = grb_file2->$RHISBL$(level575,:,:)
    rhsum2 = rhsum2 + rh57_2
  end if
  if (isvar("level550")) then
    rh55 = grb_file->$RHISBL$(level550,:,:)
    rhsum = rhsum + rh55
    rh55_2 = grb_file2->$RHISBL$(level550,:,:)
    rhsum2 = rhsum2 + rh55_2
  end if
  if (isvar("level525")) then
    rh52 = grb_file->$RHISBL$(level525,:,:)
    rhsum = rhsum + rh52
    rh52_2 = grb_file2->$RHISBL$(level525,:,:)
    rhsum2 = rhsum2 + rh52_2
  end if
  if (isvar("level500")) then
    rh50 = grb_file->$RHISBL$(level500,:,:)
    rhsum = rhsum + rh50
    rh50_2 = grb_file2->$RHISBL$(level500,:,:)
    rhsum2 = rhsum2 + rh50_2
  end if
;
  tvar = rhsum / lcount
;
  tvar2 = rhsum2 / lcount
;
end if
;
if (domainname .eq. "174") then
;
  ny = filevardimsizes(grb_file,LAT)
  nx = filevardimsizes(grb_file,LON)
;
  lat = grb_file->$LAT$(:)
  lon = grb_file->$LON$(:)
  lonorig = lon
  fliplon = lonFlip(lon)
  rhsum = new((/ny,nx/),float)
  rhsum = 0.
  rhsum2 = new((/ny,nx/),float)
  rhsum2 = 0.
  if (isvar("level850")) then
    rh85 = grb_file->$RHISBL$(level850,:,:)
    rhsum = rhsum + rh85
    rh85_2 = grb_file2->$RHISBL$(level850,:,:)
    rhsum2 = rhsum2 + rh85_2
  end if
  if (isvar("level825")) then
    rh82 = grb_file->$RHISBL$(level825,:,:)
    rhsum = rhsum + rh82
    rh82_2 = grb_file2->$RHISBL$(level825,:,:)
    rhsum2 = rhsum2 + rh82_2
  end if
  if (isvar("level800")) then
    rh80 = grb_file->$RHISBL$(level800,:,:)
    rhsum = rhsum + rh80
    rh80_2 = grb_file2->$RHISBL$(level800,:,:)
    rhsum2 = rhsum2 + rh80_2
  end if
  if (isvar("level775")) then
    rh77 = grb_file->$RHISBL$(level775,:,:)
    rhsum = rhsum + rh77
    rh77_2 = grb_file2->$RHISBL$(level775,:,:)
    rhsum2 = rhsum2 + rh77_2
  end if
  if (isvar("level750")) then
    rh75 = grb_file->$RHISBL$(level750,:,:)
    rhsum = rhsum + rh75
    rh75_2 = grb_file2->$RHISBL$(level750,:,:)
    rhsum2 = rhsum2 + rh75_2
  end if
  if (isvar("level725")) then
    rh72 = grb_file->$RHISBL$(level725,:,:)
    rhsum = rhsum + rh72
    rh72_2 = grb_file2->$RHISBL$(level725,:,:)
    rhsum2 = rhsum2 + rh72_2
  end if
  if (isvar("level700")) then
    rh70 = grb_file->$RHISBL$(level700,:,:)
    rhsum = rhsum + rh70
    rh70_2 = grb_file2->$RHISBL$(level700,:,:)
    rhsum2 = rhsum2 + rh70_2
  end if
  if (isvar("level675")) then
    rh67 = grb_file->$RHISBL$(level675,:,:)
    rhsum = rhsum + rh67
    rh67_2 = grb_file2->$RHISBL$(level675,:,:)
    rhsum2 = rhsum2 + rh67_2
  end if
  if (isvar("level650")) then
    rh65 = grb_file->$RHISBL$(level650,:,:)
    rhsum = rhsum + rh65
    rh65_2 = grb_file2->$RHISBL$(level650,:,:)
    rhsum2 = rhsum2 + rh65_2
  end if
  if (isvar("level625")) then
    rh62 = grb_file->$RHISBL$(level625,:,:)
    rhsum = rhsum + rh62
    rh62_2 = grb_file2->$RHISBL$(level625,:,:)
    rhsum2 = rhsum2 + rh62_2
  end if
  if (isvar("level600")) then
    rh60 = grb_file->$RHISBL$(level600,:,:)
    rhsum = rhsum + rh60
    rh60_2 = grb_file2->$RHISBL$(level600,:,:)
    rhsum2 = rhsum2 + rh60_2
  end if
  if (isvar("level575")) then
    rh57 = grb_file->$RHISBL$(level575,:,:)
    rhsum = rhsum + rh57
    rh57_2 = grb_file2->$RHISBL$(level575,:,:)
    rhsum2 = rhsum2 + rh57_2
  end if
  if (isvar("level550")) then
    rh55 = grb_file->$RHISBL$(level550,:,:)
    rhsum = rhsum + rh55
    rh55_2 = grb_file2->$RHISBL$(level550,:,:)
    rhsum2 = rhsum2 + rh55_2
  end if
  if (isvar("level525")) then
    rh52 = grb_file->$RHISBL$(level525,:,:)
    rhsum = rhsum + rh52
    rh52_2 = grb_file2->$RHISBL$(level525,:,:)
    rhsum2 = rhsum2 + rh52_2
  end if
  if (isvar("level500")) then
    rh50 = grb_file->$RHISBL$(level500,:,:)
    rhsum = rhsum + rh50
    rh50_2 = grb_file2->$RHISBL$(level500,:,:)
    rhsum2 = rhsum2 + rh50_2
  end if
;
  tvar = rhsum / lcount
;
  tvar2 = rhsum2 / lcount
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
  rhsum = new((/ny,nx/),float)
  rhsum = 0.
  rhsum2 = new((/ny,nx/),float)
  rhsum2 = 0.
  if (isvar("level850")) then
    rh85 = grb_file->$RHISBL$(level850,:,:)
    rhsum = rhsum + rh85
    rh85_2 = grb_file2->$RHISBL$(level850,:,:)
    rhsum2 = rhsum2 + rh85_2
  end if
  if (isvar("level825")) then
    rh82 = grb_file->$RHISBL$(level825,:,:)
    rhsum = rhsum + rh82
    rh82_2 = grb_file2->$RHISBL$(level825,:,:)
    rhsum2 = rhsum2 + rh82_2
  end if
  if (isvar("level800")) then
    rh80 = grb_file->$RHISBL$(level800,:,:)
    rhsum = rhsum + rh80
    rh80_2 = grb_file2->$RHISBL$(level800,:,:)
    rhsum2 = rhsum2 + rh80_2
  end if
  if (isvar("level775")) then
    rh77 = grb_file->$RHISBL$(level775,:,:)
    rhsum = rhsum + rh77
    rh77_2 = grb_file2->$RHISBL$(level775,:,:)
    rhsum2 = rhsum2 + rh77_2
  end if
  if (isvar("level750")) then
    rh75 = grb_file->$RHISBL$(level750,:,:)
    rhsum = rhsum + rh75
    rh75_2 = grb_file2->$RHISBL$(level750,:,:)
    rhsum2 = rhsum2 + rh75_2
  end if
  if (isvar("level725")) then
    rh72 = grb_file->$RHISBL$(level725,:,:)
    rhsum = rhsum + rh72
    rh72_2 = grb_file2->$RHISBL$(level725,:,:)
    rhsum2 = rhsum2 + rh72_2
  end if
  if (isvar("level700")) then
    rh70 = grb_file->$RHISBL$(level700,:,:)
    rhsum = rhsum + rh70
    rh70_2 = grb_file2->$RHISBL$(level700,:,:)
    rhsum2 = rhsum2 + rh70_2
  end if
  if (isvar("level675")) then
    rh67 = grb_file->$RHISBL$(level675,:,:)
    rhsum = rhsum + rh67
    rh67_2 = grb_file2->$RHISBL$(level675,:,:)
    rhsum2 = rhsum2 + rh67_2
  end if
  if (isvar("level650")) then
    rh65 = grb_file->$RHISBL$(level650,:,:)
    rhsum = rhsum + rh65
    rh65_2 = grb_file2->$RHISBL$(level650,:,:)
    rhsum2 = rhsum2 + rh65_2
  end if
  if (isvar("level625")) then
    rh62 = grb_file->$RHISBL$(level625,:,:)
    rhsum = rhsum + rh62
    rh62_2 = grb_file2->$RHISBL$(level625,:,:)
    rhsum2 = rhsum2 + rh62_2
  end if
  if (isvar("level600")) then
    rh60 = grb_file->$RHISBL$(level600,:,:)
    rhsum = rhsum + rh60
    rh60_2 = grb_file2->$RHISBL$(level600,:,:)
    rhsum2 = rhsum2 + rh60_2
  end if
  if (isvar("level575")) then
    rh57 = grb_file->$RHISBL$(level575,:,:)
    rhsum = rhsum + rh57
    rh57_2 = grb_file2->$RHISBL$(level575,:,:)
    rhsum2 = rhsum2 + rh57_2
  end if
  if (isvar("level550")) then
    rh55 = grb_file->$RHISBL$(level550,:,:)
    rhsum = rhsum + rh55
    rh55_2 = grb_file2->$RHISBL$(level550,:,:)
    rhsum2 = rhsum2 + rh55_2
  end if
  if (isvar("level525")) then
    rh52 = grb_file->$RHISBL$(level525,:,:)
    rhsum = rhsum + rh52
    rh52_2 = grb_file2->$RHISBL$(level525,:,:)
    rhsum2 = rhsum2 + rh52_2
  end if
  if (isvar("level500")) then
    rh50 = grb_file->$RHISBL$(level500,:,:)
    rhsum = rhsum + rh50
    rh50_2 = grb_file2->$RHISBL$(level500,:,:)
    rhsum2 = rhsum2 + rh50_2
  end if
;
  tvar = rhsum / lcount
;
  tvar2 = rhsum2 / lcount
;
end if
;
if (domainname .eq. "africa" .or. domainname .eq. "europe" .or. \
    domainname .eq. "w_pacific" .or. domainname .eq. "floating" .or. \
    domainname .eq. "e_pacific") then
;
  nx = filevardimsizes(grb_file,LON)
  ny = filevardimsizes(grb_file,LAT)
;
  lat = grb_file->$LAT$(:)
  lon = grb_file->$LON$(:)
  rhsum = new((/ny,nx/),float)
  rhsum = 0.
  rhsum2 = new((/ny,nx/),float)
  rhsum2 = 0.
  if (isvar("level850")) then
    rh85 = grb_file->$RHISBL$(level850,:,:)
    rhsum = rhsum + rh85
    rh85_2 = grb_file2->$RHISBL$(level850,:,:)
    rhsum2 = rhsum2 + rh85_2
  end if
  if (isvar("level825")) then
    rh82 = grb_file->$RHISBL$(level825,:,:)
    rhsum = rhsum + rh82
    rh82_2 = grb_file2->$RHISBL$(level825,:,:)
    rhsum2 = rhsum2 + rh82_2
  end if
  if (isvar("level800")) then
    rh80 = grb_file->$RHISBL$(level800,:,:)
    rhsum = rhsum + rh80
    rh80_2 = grb_file2->$RHISBL$(level800,:,:)
    rhsum2 = rhsum2 + rh80_2
  end if
  if (isvar("level775")) then
    rh77 = grb_file->$RHISBL$(level775,:,:)
    rhsum = rhsum + rh77
    rh77_2 = grb_file2->$RHISBL$(level775,:,:)
    rhsum2 = rhsum2 + rh77_2
  end if
  if (isvar("level750")) then
    rh75 = grb_file->$RHISBL$(level750,:,:)
    rhsum = rhsum + rh75
    rh75_2 = grb_file2->$RHISBL$(level750,:,:)
    rhsum2 = rhsum2 + rh75_2
  end if
  if (isvar("level725")) then
    rh72 = grb_file->$RHISBL$(level725,:,:)
    rhsum = rhsum + rh72
    rh72_2 = grb_file2->$RHISBL$(level725,:,:)
    rhsum2 = rhsum2 + rh72_2
  end if
  if (isvar("level700")) then
    rh70 = grb_file->$RHISBL$(level700,:,:)
    rhsum = rhsum + rh70
    rh70_2 = grb_file2->$RHISBL$(level700,:,:)
    rhsum2 = rhsum2 + rh70_2
  end if
  if (isvar("level675")) then
    rh67 = grb_file->$RHISBL$(level675,:,:)
    rhsum = rhsum + rh67
    rh67_2 = grb_file2->$RHISBL$(level675,:,:)
    rhsum2 = rhsum2 + rh67_2
  end if
  if (isvar("level650")) then
    rh65 = grb_file->$RHISBL$(level650,:,:)
    rhsum = rhsum + rh65
    rh65_2 = grb_file2->$RHISBL$(level650,:,:)
    rhsum2 = rhsum2 + rh65_2
  end if
  if (isvar("level625")) then
    rh62 = grb_file->$RHISBL$(level625,:,:)
    rhsum = rhsum + rh62
    rh62_2 = grb_file2->$RHISBL$(level625,:,:)
    rhsum2 = rhsum2 + rh62_2
  end if
  if (isvar("level600")) then
    rh60 = grb_file->$RHISBL$(level600,:,:)
    rhsum = rhsum + rh60
    rh60_2 = grb_file2->$RHISBL$(level600,:,:)
    rhsum2 = rhsum2 + rh60_2
  end if
  if (isvar("level575")) then
    rh57 = grb_file->$RHISBL$(level575,:,:)
    rhsum = rhsum + rh57
    rh57_2 = grb_file2->$RHISBL$(level575,:,:)
    rhsum2 = rhsum2 + rh57_2
  end if
  if (isvar("level550")) then
    rh55 = grb_file->$RHISBL$(level550,:,:)
    rhsum = rhsum + rh55
    rh55_2 = grb_file2->$RHISBL$(level550,:,:)
    rhsum2 = rhsum2 + rh55_2
  end if
  if (isvar("level525")) then
    rh52 = grb_file->$RHISBL$(level525,:,:)
    rhsum = rhsum + rh52
    rh52_2 = grb_file2->$RHISBL$(level525,:,:)
    rhsum2 = rhsum2 + rh52_2
  end if
  if (isvar("level500")) then
    rh50 = grb_file->$RHISBL$(level500,:,:)
    rhsum = rhsum + rh50
    rh50_2 = grb_file2->$RHISBL$(level500,:,:)
    rhsum2 = rhsum2 + rh50_2
  end if
;
  tvar = rhsum / lcount
;
  tvar2 = rhsum2 / lcount
;
end if
;
if (domainname .eq. "201" .or domainname .eq. "224") then
;
  fimsizes = filevardimsizes(grb_file,LAT)
  ny = fimsizes(0)
  nx = fimsizes(1)
;
  lat = grb_file->$LAT$(:,:)
  lon = grb_file->$LON$(:,:)
  rhsum = new((/ny,nx/),float)
  rhsum = 0.
  rhsum2 = new((/ny,nx/),float)
  rhsum2 = 0.
  if (isvar("level850")) then
    rh85 = grb_file->$RHISBL$(level850,:,:)
    rhsum = rhsum + rh85
    rh85_2 = grb_file2->$RHISBL$(level850,:,:)
    rhsum2 = rhsum2 + rh85_2
  end if
  if (isvar("level825")) then
    rh82 = grb_file->$RHISBL$(level825,:,:)
    rhsum = rhsum + rh82
    rh82_2 = grb_file2->$RHISBL$(level825,:,:)
    rhsum2 = rhsum2 + rh82_2
  end if
  if (isvar("level800")) then
    rh80 = grb_file->$RHISBL$(level800,:,:)
    rhsum = rhsum + rh80
    rh80_2 = grb_file2->$RHISBL$(level800,:,:)
    rhsum2 = rhsum2 + rh80_2
  end if
  if (isvar("level775")) then
    rh77 = grb_file->$RHISBL$(level775,:,:)
    rhsum = rhsum + rh77
    rh77_2 = grb_file2->$RHISBL$(level775,:,:)
    rhsum2 = rhsum2 + rh77_2
  end if
  if (isvar("level750")) then
    rh75 = grb_file->$RHISBL$(level750,:,:)
    rhsum = rhsum + rh75
    rh75_2 = grb_file2->$RHISBL$(level750,:,:)
    rhsum2 = rhsum2 + rh75_2
  end if
  if (isvar("level725")) then
    rh72 = grb_file->$RHISBL$(level725,:,:)
    rhsum = rhsum + rh72
    rh72_2 = grb_file2->$RHISBL$(level725,:,:)
    rhsum2 = rhsum2 + rh72_2
  end if
  if (isvar("level700")) then
    rh70 = grb_file->$RHISBL$(level700,:,:)
    rhsum = rhsum + rh70
    rh70_2 = grb_file2->$RHISBL$(level700,:,:)
    rhsum2 = rhsum2 + rh70_2
  end if
  if (isvar("level675")) then
    rh67 = grb_file->$RHISBL$(level675,:,:)
    rhsum = rhsum + rh67
    rh67_2 = grb_file2->$RHISBL$(level675,:,:)
    rhsum2 = rhsum2 + rh67_2
  end if
  if (isvar("level650")) then
    rh65 = grb_file->$RHISBL$(level650,:,:)
    rhsum = rhsum + rh65
    rh65_2 = grb_file2->$RHISBL$(level650,:,:)
    rhsum2 = rhsum2 + rh65_2
  end if
  if (isvar("level625")) then
    rh62 = grb_file->$RHISBL$(level625,:,:)
    rhsum = rhsum + rh62
    rh62_2 = grb_file2->$RHISBL$(level625,:,:)
    rhsum2 = rhsum2 + rh62_2
  end if
  if (isvar("level600")) then
    rh60 = grb_file->$RHISBL$(level600,:,:)
    rhsum = rhsum + rh60
    rh60_2 = grb_file2->$RHISBL$(level600,:,:)
    rhsum2 = rhsum2 + rh60_2
  end if
  if (isvar("level575")) then
    rh57 = grb_file->$RHISBL$(level575,:,:)
    rhsum = rhsum + rh57
    rh57_2 = grb_file2->$RHISBL$(level575,:,:)
    rhsum2 = rhsum2 + rh57_2
  end if
  if (isvar("level550")) then
    rh55 = grb_file->$RHISBL$(level550,:,:)
    rhsum = rhsum + rh55
    rh55_2 = grb_file2->$RHISBL$(level550,:,:)
    rhsum2 = rhsum2 + rh55_2
  end if
  if (isvar("level525")) then
    rh52 = grb_file->$RHISBL$(level525,:,:)
    rhsum = rhsum + rh52
    rh52_2 = grb_file2->$RHISBL$(level525,:,:)
    rhsum2 = rhsum2 + rh52_2
  end if
  if (isvar("level500")) then
    rh50 = grb_file->$RHISBL$(level500,:,:)
    rhsum = rhsum + rh50
    rh50_2 = grb_file2->$RHISBL$(level500,:,:)
    rhsum2 = rhsum2 + rh50_2
  end if
;
  tvar = rhsum / lcount
;
  tvar2 = rhsum2 / lcount
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
  rhsum = new((/ny,nx/),float)
  rhsum = 0.
  rhsum2 = new((/ny,nx/),float)
  rhsum2 = 0.
  if (isvar("level850")) then
    rh85 = grb_file->$RHISBL$(level850,:,:)
    rhsum = rhsum + rh85
    rh85_2 = grb_file2->$RHISBL$(level850,:,:)
    rhsum2 = rhsum2 + rh85_2
  end if
  if (isvar("level825")) then
    rh82 = grb_file->$RHISBL$(level825,:,:)
    rhsum = rhsum + rh82
    rh82_2 = grb_file2->$RHISBL$(level825,:,:)
    rhsum2 = rhsum2 + rh82_2
  end if
  if (isvar("level800")) then
    rh80 = grb_file->$RHISBL$(level800,:,:)
    rhsum = rhsum + rh80
    rh80_2 = grb_file2->$RHISBL$(level800,:,:)
    rhsum2 = rhsum2 + rh80_2
  end if
  if (isvar("level775")) then
    rh77 = grb_file->$RHISBL$(level775,:,:)
    rhsum = rhsum + rh77
    rh77_2 = grb_file2->$RHISBL$(level775,:,:)
    rhsum2 = rhsum2 + rh77_2
  end if
  if (isvar("level750")) then
    rh75 = grb_file->$RHISBL$(level750,:,:)
    rhsum = rhsum + rh75
    rh75_2 = grb_file2->$RHISBL$(level750,:,:)
    rhsum2 = rhsum2 + rh75_2
  end if
  if (isvar("level725")) then
    rh72 = grb_file->$RHISBL$(level725,:,:)
    rhsum = rhsum + rh72
    rh72_2 = grb_file2->$RHISBL$(level725,:,:)
    rhsum2 = rhsum2 + rh72_2
  end if
  if (isvar("level700")) then
    rh70 = grb_file->$RHISBL$(level700,:,:)
    rhsum = rhsum + rh70
    rh70_2 = grb_file2->$RHISBL$(level700,:,:)
    rhsum2 = rhsum2 + rh70_2
  end if
  if (isvar("level675")) then
    rh67 = grb_file->$RHISBL$(level675,:,:)
    rhsum = rhsum + rh67
    rh67_2 = grb_file2->$RHISBL$(level675,:,:)
    rhsum2 = rhsum2 + rh67_2
  end if
  if (isvar("level650")) then
    rh65 = grb_file->$RHISBL$(level650,:,:)
    rhsum = rhsum + rh65
    rh65_2 = grb_file2->$RHISBL$(level650,:,:)
    rhsum2 = rhsum2 + rh65_2
  end if
  if (isvar("level625")) then
    rh62 = grb_file->$RHISBL$(level625,:,:)
    rhsum = rhsum + rh62
    rh62_2 = grb_file2->$RHISBL$(level625,:,:)
    rhsum2 = rhsum2 + rh62_2
  end if
  if (isvar("level600")) then
    rh60 = grb_file->$RHISBL$(level600,:,:)
    rhsum = rhsum + rh60
    rh60_2 = grb_file2->$RHISBL$(level600,:,:)
    rhsum2 = rhsum2 + rh60_2
  end if
  if (isvar("level575")) then
    rh57 = grb_file->$RHISBL$(level575,:,:)
    rhsum = rhsum + rh57
    rh57_2 = grb_file2->$RHISBL$(level575,:,:)
    rhsum2 = rhsum2 + rh57_2
  end if
  if (isvar("level550")) then
    rh55 = grb_file->$RHISBL$(level550,:,:)
    rhsum = rhsum + rh55
    rh55_2 = grb_file2->$RHISBL$(level550,:,:)
    rhsum2 = rhsum2 + rh55_2
  end if
  if (isvar("level525")) then
    rh52 = grb_file->$RHISBL$(level525,:,:)
    rhsum = rhsum + rh52
    rh52_2 = grb_file2->$RHISBL$(level525,:,:)
    rhsum2 = rhsum2 + rh52_2
  end if
  if (isvar("level500")) then
    rh50 = grb_file->$RHISBL$(level500,:,:)
    rhsum = rhsum + rh50
    rh50_2 = grb_file2->$RHISBL$(level500,:,:)
    rhsum2 = rhsum2 + rh50_2
  end if
;
  tvar = rhsum / lcount
;
  tvar2 = rhsum2 / lcount
;
end if
;
load "$NCL_HOME/obj/get_date_str_FIM.ncl"
;
tempvar = tvar
tvar = tvar - tvar2
print(date_str)
print("850-500 mean rh")
print(min(tvar))
print(max(tvar))
if (domainname .eq. "174") then
  tvarorig = tvar
  fliptvar = lonFlip(tvar)
end if
;
if (min(tvar) .eq. 0 .and. max(tvar) .eq. 0) then
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
if (domainname .eq. "130") then
;
  load "$NCL_HOME/obj/ruc_tiles_loop.ncl"
;
  if (iq .eq. 0) then
    load "$NCL_HOME/obj/lambert_conformal_map_blk.ncl"
    setvalues mpid
      "mpGeophysicalLineColor" : "background"
      "mpGeophysicalLineThicknessF" : 2
      "mpNationalLineColor" : "background"
      "mpNationalLineThicknessF" : 2
      "mpUSStateLineColor" : "background"
      "mpUSStateLineThicknessF" : 2
    end setvalues
  else
    load "$NCL_HOME/obj/lc_map_counties.ncl"
    setvalues mpid
      "mpGeophysicalLineColor" : "background"
      "mpGeophysicalLineThicknessF" : 2
      "mpNationalLineColor" : "background"
      "mpNationalLineThicknessF" : 2
      "mpUSStateLineColor" : "background"
      "mpUSStateLineThicknessF" : 2
      "mpCountyLineColor" : "background"
      "mpCountyLineThicknessF" : 1
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
  load "$NCL_HOME/obj/diff_-40-40by5.ncl"
;
  setvalues cnid
    "cnRasterModeOn" : False
  end setvalues
;
  titlestring = "850-500mb mean Rel Hum (%)"
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
      "sfDataArray" : tvar
      "sfXCStartSubsetV" : llx
      "sfYCStartSubsetV" : lly
      "sfXCEndSubsetV" : urx
      "sfYCEndSubsetV" : ury
      "sfMissingValueV" : tvar@_FillValue
  end create
;
  load "$NCL_HOME/obj/diff_-40-40by5.ncl"
;
  titlestring = "850-500mb mean Rel Hum (%)"
  load "$NCL_HOME/obj/add_titles.ncl"
  load "$NCL_HOME/obj/add_explabel.ncl"
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
  setvalues mpid
    "mpGeophysicalLineColor" : "background"
    "mpGeophysicalLineThicknessF" : 2
    "mpNationalLineColor" : "background"
    "mpNationalLineThicknessF" : 2
    "mpUSStateLineColor" : "background"
    "mpUSStateLineThicknessF" : 2
  end setvalues
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
    tempfield = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : tempvar
        "sfXCStartV" : 0
        "sfYCStartV" : ydim
        "sfXCEndV" : xdim
        "sfYCEndV" : 0
        "sfMissingValueV" : tempvar@_FillValue
    end create
;
  end if
;
  if (domainname .eq. "236" .or. domainname .eq. "83" .or. \
      domainname .eq. "201" .or. domainname .eq. "224") then
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
    tempfield = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : tempvar
        "sfXCStartV" : 0
        "sfYCStartV" : 0
        "sfXCEndV" : xdim
        "sfYCEndV" : ydim
        "sfMissingValueV" : tempvar@_FillValue
    end create
;
  end if
;
  load "$NCL_HOME/obj/diff_-40-40by5.ncl"
;
  if (domainname .eq. "236" .or. domainname .eq. "244" .or. domainname .eq. "83" .or. \
      domainname .eq. "201" .or. domainname .eq. "224") then
;
    setvalues cnid
;      "cnRasterModeOn" : False
      "cnRasterSmoothingOn" : True
    end setvalues
;
  end if
;
  titlestring = "850-500mb mean Rel Hum (%)"
  load "$NCL_HOME/obj/add_titles.ncl"
;
; Draw map
;
  draw(cnid)
  draw(mpid)
  frame(wid)
;
end if  ; domainname .ne. "130" .and. domainname .ne. "174"
;
end
;
