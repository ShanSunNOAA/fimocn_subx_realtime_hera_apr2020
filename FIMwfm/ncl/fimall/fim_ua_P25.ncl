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
  ash1 = grb_file->CIN_4_ISBL(:,:,:)
  ash2 = grb_file->GFLUX_4_ISBL(:,:,:)
  hgt500 = grb_file->GH_4_ISBL(level,:,:) / 10. ; use decameters
  isohgt = grb_file->GH_4_ISBL(:,:,:) * 3.28084 ; convert to feet
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
  ash1 = grb_file->CIN_83_ISBL(:,:,:)
  ash2 = grb_file->GFLUX_83_ISBL(:,:,:)
  hgt500 = grb_file->GH_83_ISBL(level,:,:) / 10. ; use decameters
  isohgt = grb_file->GH_83_ISBL(:,:,:) * 3.28084 ; convert to feet
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
  latin1 = grb_file->gridlon_130@Latin1
  latin2 = grb_file->gridlon_130@Latin2
  lov = grb_file->gridlon_130@Lov
  ash1 = grb_file->CIN_130_HYBL(:,:,:)
  ash2 = grb_file->GFLUX_130_HYBL(:,:,:)
  hgt500 = grb_file->GH_130_ISBL(level,:,:) / 10. ; use decameters
  isohgt = grb_file->GH_130_ISBL(:,:,:) * 3.28084 ; convert to feet
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
  ash1 = grb_file->CIN_236_HYBL(:,:,:)
  ash2 = grb_file->GFLUX_236_HYBL(:,:,:)
  hgt500 = grb_file->GH_236_ISBL(level,:,:) / 10. ; use decameters
  isohgt = grb_file->GH_236_ISBL(:,:,:) * 3.28084 ; convert to feet
;
end if
;
if (domainname .eq. "africa" .or. domainname .eq. "europe" .or. \
    domainname .eq. "w_pacific" .or. domainname .eq. "floating" .or. \
    domainname .eq. "e_pacific") then
;
  nx = filevardimsizes(grb_file,"g0_lon_2")
  ny = filevardimsizes(grb_file,"g0_lat_1")
;
  intm = grb_file->T_GDS0_SFC@initial_time
  fctm = grb_file->T_GDS0_SFC@forecast_time
  lat = grb_file->g0_lat_1(:)
  lon = grb_file->g0_lon_2(:)
  if (isfilevar(grb_file,"CIN_GDS0_ISBL")) then
    ash1 = grb_file->CIN_GDS0_ISBL(:,:,:)
  else
    print("variable not found")
    exit
  end if
  if (isfilevar(grb_file,"GFLUX_GDS0_ISBL")) then
    ash2 = grb_file->GFLUX_GDS0_ISBL(:,:,:)
  else
    print("variable not found")
    exit
  end if
  hgt500 = grb_file->GH_GDS0_ISBL(level,:,:) / 10. ; use decameters
  isohgt = grb_file->GH_GDS0_ISBL(:,:,:) * 3.28084 ; convert to feet
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
  ash1 = grb_file->CIN_219_HYBL(:,:,:)
  ash2 = grb_file->GFLUX_219_HYBL(:,:,:)
;  ash3 = grb_file->CSULF_219_HYBL(:,:,:)
;  ash4 = grb_file->CSDLF_219_HYBL(:,:,:)
  hgt500 = grb_file->GH_219_ISBL(level,:,:) / 10. ; use decameters
  isohgt = grb_file->GH_219_ISBL(:,:,:) * 3.28084 ; convert to feet
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
  ash1 = grb_file->CIN_28_HYBL(:,:,:)
  ash2 = grb_file->GFLUX_28_HYBL(:,:,:)
;  ash3 = grb_file->CSULF_28_HYBL(:,:,:)
;  ash4 = grb_file->CSDLF_28_HYBL(:,:,:)
  hgt500 = grb_file->GH_28_ISBL(level,:,:) / 10. ; use decameters
  isohgt = grb_file->GH_28_ISBL(:,:,:) * 3.28084 ; convert to feet
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
  lon = grb_file->lon_244(:)
  ash1 = grb_file->CIN_244_HYBL(:,:,:)
  ash2 = grb_file->GFLUX_244_HYBL(:,:,:)
;  ash3 = grb_file->CSULF_244_HYBL(:,:,:)
;  ash4 = grb_file->CSDLF_244_HYBL(:,:,:)
  hgt500 = grb_file->GH_244_ISBL(level,:,:) / 10. ; use decameters
  isohgt = grb_file->GH_244_ISBL(:,:,:) * 3.28084 ; convert to feet
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
print("ashlow")
;ashlow = log10(ashlow)
;;ashlow = log(ashlow + sqrt((ashlow*ashlow) + 1))  ; inv hyp sin
if (min(ashlow) .le. 0.1 .and. max(ashlow) .le. 0.1) then
  ashlow(0,0) = .001
end if
ashlow = (ashlow / 100.) * 100.
print(min(ashlow))
print(max(ashlow))
print("ashmid")
;ashmid = log10(ashmid)
;;ashmid = log(ashmid + sqrt((ashmid*ashmid) + 1))  ; inv hyp sin
if (min(ashmid) .le. 0.1 .and. max(ashmid) .le. 0.1) then
  ashmid(0,0) = .001
end if
print(min(ashmid))
print(max(ashmid))
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
           "wkMetaName" : "./ua_ash.ncgm"
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
           "wkPSFileName" : "./ua_ash.ps"
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
  tfield1 = create "ScalarField" scalarFieldClass appid
      "sfDataArray" : ashlow
      "sfXCStartSubsetV" : llx
      "sfYCStartSubsetV" : lly
      "sfXCEndSubsetV" : urx
      "sfYCEndSubsetV" : ury
      "sfMissingValueV" : ashlow@_FillValue
  end create
;
  tfield2 = create "ScalarField" scalarFieldClass appid
      "sfDataArray" : ashmid
      "sfXCStartSubsetV" : llx
      "sfYCStartSubsetV" : lly
      "sfXCEndSubsetV" : urx
      "sfYCEndSubsetV" : ury
      "sfMissingValueV" : ashmid@_FillValue
  end create
;
  tfield3 = create "ScalarField" scalarFieldClass appid
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
  do iplot = 1,3
;
    if (iplot .eq. 1) then
      tfield = tfield1
    end if
    if (iplot .eq. 2) then
      tfield = tfield2
    end if
    if (iplot .eq. 3) then
      tfield = tfield3
    end if
;
    load "$NCL_HOME/obj/uaash.ncl"
;
    setvalues cnid
      "cnRasterModeOn" : False
    end setvalues
;
    load "$NCL_HOME/obj/height_60m.ncl"
;
    if (iplot .eq. 1) then
      titlestring = "Vert Avg Unspeciated Aerosols & Volcanic Ash, Sfc-20kft (~F33~m~F22~g/kg)"
    end if
    if (iplot .eq. 2) then
      titlestring = "Vert Avg Unspeciated Aerosols & Volcanic Ash, 20kft-35kft (~F33~m~F22~g/kg)"
    end if
    if (iplot .eq. 3) then
      titlestring = "Vert Avg Unspeciated Aerosols & Volcanic Ash, 35kft-50kft (~F33~m~F22~g/kg)"
    end if
    load "$NCL_HOME/obj/add_titles.ncl"
    load "$NCL_HOME/obj/add_explabel.ncl"
;
; Draw map
;
    draw(cnid)
    draw(mpid)
    draw(vcid)
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
  end do  ; end of do iplot loop
;
  end do  ; loop begins in ruc_tiles_loop.ncl
;
else  ; domain is not 130
;
  load "$NCL_HOME/obj/fim_maps_blk.ncl"
;
; Create ScalarField objects.
;
  if (domainname .eq. "fim" .or. domainname .eq. "244") then
;
    tfield1 = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : ashlow
        "sfXCStartV" : 0
        "sfYCStartV" : ydim
        "sfXCEndV" : xdim
        "sfYCEndV" : 0
        "sfMissingValueV" : ashlow@_FillValue
    end create
;
    tfield2 = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : ashmid
        "sfXCStartV" : 0
        "sfYCStartV" : ydim
        "sfXCEndV" : xdim
        "sfYCEndV" : 0
        "sfMissingValueV" : ashmid@_FillValue
    end create
;
    tfield3 = create "ScalarField" scalarFieldClass appid
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
  if (domainname .eq. "236" .or. domainname .eq. "africa" .or. \
      domainname .eq. "219" .or. domainname .eq. "w_pacific" .or. \
      domainname .eq. "europe" .or. domainname .eq. "floating" .or. \
      domainname .eq. "e_pacific" .or. domainname .eq. "83") then
;
    tfield1 = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : ashlow
        "sfXCStartV" : 0
        "sfYCStartV" : 0
        "sfXCEndV" : xdim
        "sfYCEndV" : ydim
        "sfMissingValueV" : ashlow@_FillValue
    end create
;
    tfield2 = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : ashmid
        "sfXCStartV" : 0
        "sfYCStartV" : 0
        "sfXCEndV" : xdim
        "sfYCEndV" : ydim
        "sfMissingValueV" : ashmid@_FillValue
    end create
;
    tfield3 = create "ScalarField" scalarFieldClass appid
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
  if (domainname .eq. "28") then
;
    tfield1 = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : ashlow
        "sfXCStartV" : xdim
        "sfYCStartV" : ydim
        "sfXCEndV" : 0
        "sfYCEndV" : 0
        "sfMissingValueV" : ashlow@_FillValue
    end create
;
    tfield2 = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : ashmid
        "sfXCStartV" : xdim
        "sfYCStartV" : ydim
        "sfXCEndV" : 0
        "sfYCEndV" : 0
        "sfMissingValueV" : ashmid@_FillValue
    end create
;
    tfield3 = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : ashhigh
        "sfXCStartV" : xdim
        "sfYCStartV" : ydim
        "sfXCEndV" : 0
        "sfYCEndV" : 0
        "sfMissingValueV" : ashhigh@_FillValue
    end create
;
    hfield = create "ScalarField" scalarFieldClass appid
        "sfDataArray" : hgt500
        "sfXCStartV" : xdim
        "sfYCStartV" : ydim
        "sfXCEndV" : 0
        "sfYCEndV" : 0
        "sfMissingValueV" : hgt500@_FillValue
    end create
;
  end if
;
  do iplot = 1,3
;
    if (iplot .eq. 1) then
      tfield = tfield1
    end if
    if (iplot .eq. 2) then
      tfield = tfield2
    end if
    if (iplot .eq. 3) then
      tfield = tfield3
    end if
;
    load "$NCL_HOME/obj/uaash.ncl"
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
    if (iplot .eq. 1) then
      titlestring = "Vert Avg Unspeciated Aerosols & Volcanic Ash, Sfc-20kft (~F33~m~F22~g/kg)"
    end if
    if (iplot .eq. 2) then
      titlestring = "Vert Avg Unspeciated Aerosols & Volcanic Ash, 20kft-35kft (~F33~m~F22~g/kg)"
    end if
    if (iplot .eq. 3) then
      titlestring = "Vert Avg Unspeciated Aerosols & Volcanic Ash, 35kft-50kft (~F33~m~F22~g/kg)"
    end if
    load "$NCL_HOME/obj/add_titles.ncl"
    load "$NCL_HOME/obj/add_explabel.ncl"
;
; Draw map
;
    draw(cnid)
    draw(mpid)
;    draw(cnid4)
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
  end do  ; end of do iplot loop
;
end if  ; domainname .eq. "130"
;
;
end

