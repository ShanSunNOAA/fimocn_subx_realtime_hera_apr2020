; do loop to make subdomain tiles for domain 242 (Alaska)
;
CENLAT = 90 ; not given in grib file
;
do iq = 0,2
;
if (iq .eq. 0) then ; full domain
  llx = 0
  lly = 0
  urx = nx - 1
  ury = ny - 1
  wstride = 10
end if
;
if (iq .eq. 1) then ; zoomed mainland tile
  llx = 220
  lly = 180
  urx = 365
  ury = 325
  wstride = 3
end if
;
if (iq .eq. 2) then ; zoom including Kamchatka
  llx = 10
  lly = 116
  urx = 417
  ury = 407
  wstride = 6
end if
;
tmpnx = urx - llx
tmpny = ury - lly
;xdim = nx - 1
;ydim = ny - 1
xdim = tmpnx
ydim = tmpny
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
LcLon = lon(lly,llx)
LcLat = lat(lly,llx)
RcLon = lon(ury,urx)
RcLat = lat(ury,urx)

