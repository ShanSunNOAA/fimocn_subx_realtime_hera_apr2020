;
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_LIB/ncarg/nclscripts/csm/contributed.ncl"
modl = getenv("MODL")
interval = getenv("FCST_INTERVAL")
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
  if ( fctm .gt. 0 ) then
    tvar = grb_file->$TP$(:,:) * 0.03937 ; convert to inches
    tvar2 = grb_file2->$TP$(:,:) * 0.03937 ; convert to inches
  else
    exit
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
  tvar = grb_file->$TP$(:,:) * 0.03937 ; convert to inches
  tvar2 = grb_file2->$TP$(:,:) * 0.03937 ; convert to inches
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
  tvar = grb_file->$TP$(:,:) * 0.03937 ; convert to inches
  tvar2 = grb_file2->$TP$(:,:) * 0.03937 ; convert to inches
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
  tvar = grb_file->$TP$(:,:) * 0.03937 ; convert to inches
  tvar2 = grb_file2->$TP$(:,:) * 0.03937 ; convert to inches
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
  tvar = grb_file->$TP$(::-1,:) * 0.03937 ; convert to inches
  tvar2 = grb_file2->$TP$(::-1,:) * 0.03937 ; convert to inches
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
;
  if ( fctm .gt. 0 ) then
    tvar = grb_file->$TP$(:,:) * 0.03937 ; convert to inches
    tvar2 = grb_file2->$TP$(:,:) * 0.03937 ; convert to inches
  else
    exit
  end if
;
end if
;
if (domainname .eq. "201" .or. domainname .eq. "224") then
;
  fimsizes = filevardimsizes(grb_file,LAT)
  ny = fimsizes(0)
  nx = fimsizes(1)
;
  lat = grb_file->$LAT$(:,:)
  lon = grb_file->$LON$(:,:)
  if ( fctm .gt. 0 ) then
    tvar = grb_file->$TP$(:,:) * 0.03937 ; convert to inches
    tvar2 = grb_file2->$TP$(:,:) * 0.03937 ; convert to inches
  else
    exit
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
    tvar = grb_file->$TP$(:,:) * 0.03937 ; convert to inches
    tvar2 = grb_file2->$TP$(:,:) * 0.03937 ; convert to inches
  else
    exit
  end if
;
end if
;
load "$NCL_HOME/obj/get_date_str_FIM.ncl"
;
tempvar = tvar
tvar = tvar - tvar2
;
print(date_str)
print("acp diff")
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
           "wkMetaName" : "./sfc_acp.ncgm"
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
           "wkPSFileName" : "./sfc_acp.ps"
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
  load "$NCL_HOME/obj/diff_-2-2.ncl"
;
  setvalues cnid
    "cnRasterModeOn" : False
  end setvalues
;
  titlestring = interval+"h Acc Precip (in)"
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
  load "$NCL_HOME/obj/diff_-2-2.ncl"
;
  titlestring = interval+"h Acc Precip (in)"
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
  load "$NCL_HOME/obj/diff_-2-2.ncl"
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
  load "$NCL_HOME/obj/precip_0-2_contours.ncl"
;
  titlestring = interval+"h Acc Precip (in)"
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
