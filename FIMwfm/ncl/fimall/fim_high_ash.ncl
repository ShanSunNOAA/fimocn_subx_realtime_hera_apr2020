;
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_LIB/ncarg/nclscripts/csm/contributed.ncl"
modl = getenv("MODL")
;
load "$NCL_HOME/obj/cmap3.ncl"
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
  ash1 = grb_file->$CIN$(:,:,:)
  ash2 = grb_file->$GFLUX$(:,:,:)
  hgt500 = grb_file->$GHISBL$(level,:,:) / 10. ; use decameters
  isohgt = grb_file->$GHISBL$(:,:,:) * 3.28084 ; convert to feet
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
  ash1 = grb_file->$CIN$(:,:,:)
  ash2 = grb_file->$GFLUX$(:,:,:)
  hgt500 = grb_file->$GHISBL$(level,:,:) / 10. ; use decameters
  isohgt = grb_file->$GHISBL$(:,:,:) * 3.28084 ; convert to feet
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
  ash1 = grb_file->$CIN$(:,:,:)
  ash2 = grb_file->$GFLUX$(:,:,:)
  hgt500 = grb_file->$GHISBL$(level,:,:) / 10. ; use decameters
  isohgt = grb_file->$GHISBL$(:,:,:) * 3.28084 ; convert to feet
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
  ash1 = grb_file->$CIN$(:,:,:)
  ash2 = grb_file->$GFLUX$(:,:,:)
  hgt500 = grb_file->$GHISBL$(level,:,:) / 10. ; use decameters
  isohgt = grb_file->$GHISBL$(:,:,:) * 3.28084 ; convert to feet
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
  ash1 = grb_file->$CIN$(:,::-1,:)
  ash1orig = ash1
  flipash1 = lonFlip(ash1)
  ash2 = grb_file->$GFLUX$(:,::-1,:)
  ash2orig = ash2
  flipash2 = lonFlip(ash2)
  hgt500 = grb_file->$GHISBL$(level,::-1,:) / 10. ; use decameters
  hgt500orig = hgt500
  fliphgt500 = lonFlip(hgt500)
  isohgt = grb_file->$GHISBL$(:,::-1,:) * 3.28084 ; convert to feet
  isohgtorig = isohgt
  flipisohgt = lonFlip(isohgt)
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
  ash1 = grb_file->$CIN$(:,:,:)
  ash2 = grb_file->$GFLUX$(:,:,:)
  hgt500 = grb_file->$GHISBL$(level,:,:) / 10. ; use decameters
  isohgt = grb_file->$GHISBL$(:,:,:) * 3.28084 ; convert to feet
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
  ash1 = grb_file->$CIN$(:,:,:)
  ash2 = grb_file->$GFLUX$(:,:,:)
;  ash3 = grb_file->CSULF_201_HYBL(:,:,:)
;  ash4 = grb_file->CSDLF_201_HYBL(:,:,:)
  hgt500 = grb_file->$GHISBL$(level,:,:) / 10. ; use decameters
  isohgt = grb_file->$GHISBL$(:,:,:) * 3.28084 ; convert to feet
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
  ash1 = grb_file->$CIN$(:,:,:)
  ash2 = grb_file->$GFLUX$(:,:,:)
;  ash3 = grb_file->CSULF_244_HYBL(:,:,:)
;  ash4 = grb_file->CSDLF_244_HYBL(:,:,:)
  hgt500 = grb_file->$GHISBL$(level,:,:) / 10. ; use decameters
  isohgt = grb_file->$GHISBL$(:,:,:) * 3.28084 ; convert to feet
;
end if
;
numisolevs = grb_file->lv_ISBL3
ashlow = hgt500 * 0.
ashmid = hgt500 * 0.
ashhigh = hgt500 * 0.
do jj = 0,ny-1
  do ii = 0,nx-1
    lastlowhgt = 0.
    lastmidhgt = 0.
    lasthighhgt = 0.
    do ihlev = 1,dimsizes(numisolevs)-1
      ashtmp = 0.
      if (isohgt(ihlev,jj,ii) .le. 20000) then
        ashtmp = ash1(ihlev,jj,ii) + ash2(ihlev,jj,ii)
        dz = isohgt(ihlev,jj,ii) - isohgt(ihlev-1,jj,ii)
        ashlow(jj,ii) = ashlow(jj,ii) + (ashtmp * dz)
        lastlowhgt = isohgt(ihlev,jj,ii)
        lastlowhgtindex = ihlev
      end if
      if (isohgt(ihlev,jj,ii) .gt. 20000 .and. \
          isohgt(ihlev,jj,ii) .le. 35000) then
        ashtmp = ash1(ihlev,jj,ii) + ash2(ihlev,jj,ii)
        dz = isohgt(ihlev,jj,ii) - isohgt(ihlev-1,jj,ii)
        ashmid(jj,ii) = ashmid(jj,ii) + (ashtmp * dz)
        lastmidhgt = isohgt(ihlev,jj,ii)
        lastmidhgtindex = ihlev
      end if
      if (isohgt(ihlev,jj,ii) .gt. 35000 .and. \
          isohgt(ihlev,jj,ii) .le. 50000) then
        ashtmp = ash1(ihlev,jj,ii) + ash2(ihlev,jj,ii)
        dz = isohgt(ihlev,jj,ii) - isohgt(ihlev-1,jj,ii)
        ashmid(jj,ii) = ashmid(jj,ii) + (ashtmp * dz)
        lasthighhgt = isohgt(ihlev,jj,ii)
        lasthighhgtindex = ihlev
      end if
    end do
;
;    print("jj = "+jj+"   ii = "+ii)
;    print("ashlow before = "+ashlow(jj,ii))
;    print("lastlowhgt = "+lastlowhgt+"    hybhgt(0,jj,ii) = "+hybhgt(0,jj,ii))
    ashlow(jj,ii) = ashlow(jj,ii) / (lastlowhgt - isohgt(0,jj,ii))
;    print("ashlow after = "+ashlow(jj,ii))
;    print("ashmid before = "+ashmid(jj,ii))
;    print("lastmidhgt = "+lastmidhgt+"    hybhgt(lastlowhgtindex,jj,ii) = "+ \
;          hybhgt(lastlowhgtindex,jj,ii))
    ashmid(jj,ii) = ashmid(jj,ii) / (lastmidhgt - isohgt(lastlowhgtindex,jj,ii))
;    print("ashmid after = "+ashmid(jj,ii))
;    print("ashhigh before = "+ashhigh(jj,ii))
;    print("lasthighhgt = "+lasthighhgt+"    hybhgt(lastmidhgtindex,jj,ii) = "+ \
;          hybhgt(lastmidhgtindex,jj,ii))
    ashhigh(jj,ii) = ashhigh(jj,ii) / (lasthighhgt - isohgt(lastmidhgtindex,jj,ii))
;    print("ashhigh after = "+ashhigh(jj,ii))
;
;    if (ashlow(jj,ii) .lt. 1.) then
;      ashlow(jj,ii) = 1.
;    end if
;    if (ashmid(jj,ii) .lt. 1.) then
;      ashmid(jj,ii) = 1.
;    end if
;    if (ashhigh(jj,ii) .lt. 1.) then
;      ashhigh(jj,ii) = 1.
;    end if
  end do
end do
;
load "$NCL_HOME/obj/get_date_str_FIM.ncl"
;
print(date_str) 
print("ashhigh")
;ashhigh = log10(ashhigh)
;;ashhigh = log(ashhigh + sqrt((ashhigh*ashhigh) + 1))  ; inv hyp sin
if (min(ashhigh) .le. 0.1 .and. max(ashhigh) .le. 0.1) then
  ashhigh(0,0) = .001
end if
print(min(ashhigh))
print(max(ashhigh))
print("500 mb hgt")
print(min(hgt500))
print(max(hgt500))
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
           "wkMetaName" : "./high_ash.ncgm"
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
           "wkPSFileName" : "./high_ash.ps"
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
; Create ScalarField objects.
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
      "sfDataArray" : ashhigh
      "sfXCStartSubsetV" : llx
      "sfYCStartSubsetV" : lly
      "sfXCEndSubsetV" : urx
      "sfYCEndSubsetV" : ury
      "sfMissingValueV" : ashhigh@_FillValue
  end create
;
  hfield = create "ScalarField" scalarFieldClass appid
      "sfDataArray" : hgt500
      "sfXCStartSubsetV" : llx
      "sfYCStartSubsetV" : lly
      "sfXCEndSubsetV" : urx
      "sfYCEndSubsetV" : ury
      "sfMissingValueV" : hgt500@_FillValue
  end create
;
  load "$NCL_HOME/obj/uaash.ncl"
;
  setvalues cnid
    "cnRasterModeOn" : False
  end setvalues
;
  load "$NCL_HOME/obj/height_60m.ncl"
;
  titlestring = "Vert Avg Unspeciated Aerosols & Volcanic Ash, 35kft-50kft (~F33~m~F22~g/kg)"
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
      gsres@gsMarkerColor = 36
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
      gsres@gsMarkerColor = 36
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
      "sfDataArray" : ashhigh
      "sfXCStartSubsetV" : llx
      "sfYCStartSubsetV" : lly
      "sfXCEndSubsetV" : urx
      "sfYCEndSubsetV" : ury
      "sfMissingValueV" : ashhigh@_FillValue
  end create
;
  hfield = create "ScalarField" scalarFieldClass appid
      "sfDataArray" : hgt500
      "sfXCStartSubsetV" : llx
      "sfYCStartSubsetV" : lly
      "sfXCEndSubsetV" : urx
      "sfYCEndSubsetV" : ury
      "sfMissingValueV" : hgt500@_FillValue
  end create
;
  load "$NCL_HOME/obj/uaash.ncl"
;
  setvalues cnid
    "cnRasterModeOn" : False
  end setvalues
;
  load "$NCL_HOME/obj/height_60m.ncl"
;
  titlestring = "Vert Avg Unspeciated Aerosols & Volcanic Ash, 35kft-50kft (~F33~m~F22~g/kg)"
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
      "sfDataArray" : ashhigh
      "sfXCStartSubsetV" : llx
      "sfYCStartSubsetV" : lly
      "sfXCEndSubsetV" : urx
      "sfYCEndSubsetV" : ury
      "sfMissingValueV" : ashhigh@_FillValue
  end create
;
  hfield = create "ScalarField" scalarFieldClass appid
      "sfDataArray" : hgt500
      "sfXCStartSubsetV" : llx
      "sfYCStartSubsetV" : lly
      "sfXCEndSubsetV" : urx
      "sfYCEndSubsetV" : ury
      "sfMissingValueV" : hgt500@_FillValue
  end create
;
  load "$NCL_HOME/obj/uaash.ncl"
;
  setvalues cnid
    "cnFillMode" : "AreaFill"
  end setvalues
;
  load "$NCL_HOME/obj/height_60m.ncl"
;
  titlestring = "Vert Avg Unspeciated Aerosols & Volcanic Ash, 35kft-50kft (~F33~m~F22~g/kg)"
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
  if (domainname .eq. "fim" .or. domainname .eq. "244" .or. domainname .eq. "129") then
;
    tfield = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : ashhigh
        "sfXCStartV" : 0
        "sfYCStartV" : ydim
        "sfXCEndV" : xdim
        "sfYCEndV" : 0
        "sfMissingValueV" : ashhigh@_FillValue
    end create
;
    hfield = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : hgt500
        "sfXCStartV" : 0
        "sfYCStartV" : ydim
        "sfXCEndV" : xdim
        "sfYCEndV" : 0
        "sfMissingValueV" : hgt500@_FillValue
    end create
;
  end if
;
  if (domainname .eq. "236" .or. domainname .eq. "201" .or. domainname .eq. "224" .or. \
      domainname .eq. "83") then
;
    tfield = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : ashhigh
        "sfXCStartV" : 0
        "sfYCStartV" : 0
        "sfXCEndV" : xdim
        "sfYCEndV" : ydim
        "sfMissingValueV" : ashhigh@_FillValue
    end create
;
    hfield = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : hgt500
        "sfXCStartV" : 0
        "sfYCStartV" : 0
        "sfXCEndV" : xdim
        "sfYCEndV" : ydim
        "sfMissingValueV" : hgt500@_FillValue
    end create
;
  end if
;
  load "$NCL_HOME/obj/uaash.ncl"
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
  load "$NCL_HOME/obj/height_60m.ncl"
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
  titlestring = "Vert Avg Unspeciated Aerosols & Volcanic Ash, 35kft-50kft (~F33~m~F22~g/kg)"
  load "$NCL_HOME/obj/add_titles.ncl"
  load "$NCL_HOME/obj/add_explabel.ncl"
;
; Draw map
;
  draw(cnid)
  draw(mpid)
;  draw(cnid4)
;
  gsres = True
  if (domainname .eq. "236") then
    do i = 0,dimsizes(airport_id)-1
      gsres@gsMarkerIndex = 16
      gsres@gsMarkerColor = 36
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
;
end

