PRO overplot_fcstens, ensfcst, nmembers, colortableno, colortb, fcst_len, fcst_output_int

; --- first calculate the ensmemble-mean track
fcst_len = FIX (fcst_len)
fcst_output_int = FIX (fcst_output_int)
ensmean_lon = fltarr(29)
ensmean_lat = fltarr(29)
iflag_pos = lonarr(29)
for ilead = 0, fcst_len, fcst_output_int do begin
  idx = ilead/fcst_output_int
  nen = 0L
  for imem = 0, nmembers-1 do begin
     if (ensfcst.ctr_lat(idx,imem) gt -999) then begin
        ensmean_lon(idx) = ensmean_lon(idx) + ensfcst.ctr_lon(idx,imem)
        ensmean_lat(idx) = ensmean_lat(idx) + ensfcst.ctr_lat(idx,imem)
        iflag_pos(idx) = 1
        nen = nen+1 
     endif
  endfor
  ensmean_lon(idx) = ensmean_lon(idx) / FLOAT(nen)
  ensmean_lat(idx) = ensmean_lat(idx) / FLOAT(nen)
  print, 'fcst: ',idx, ' mean lon: ',ensmean_lon(idx)
  print, 'fcst: ',idx, ' mean lat: ',ensmean_lat(idx)
endfor

loadct,0

; ---- now plot each member and then mean track

for ilead = fcst_output_int, fcst_len, fcst_output_int do begin
  idx = ilead/fcst_output_int
  for imem = 0, nmembers-1 do begin
    enscolor=75
    if (imem LT 10) then enscolor=0
    if (imem LT 10) then thick=2
    if (imem GE 10) then thick=1
    if (ensfcst.ctr_lat(idx,imem) gt -999 and ensfcst.ctr_lat(idx-1,imem) gt -999) then begin
    ; print, 'mem: ',imem, ' TRACK lon: ',ensfcst.ctr_lon(idx,imem),' ', ensfcst.ctr_lon(idx-1,imem)
    ; print, 'mem: ',imem, ' TRACK lat: ', ensfcst.ctr_lat(idx,imem),' ',ensfcst.ctr_lat(idx-1,imem)
      oplot,[ensfcst.ctr_lon(idx,imem), ensfcst.ctr_lon(idx-1,imem)],$
        [ensfcst.ctr_lat(idx,imem),ensfcst.ctr_lat(idx-1,imem)],color=enscolor,thick=thick
    endif
 endfor
endfor

for ilead = fcst_output_int, fcst_len, fcst_output_int do begin
  idx = ilead/fcst_output_int
  if (iflag_pos(idx) NE 0 AND iflag_pos(idx-1) NE 0) then begin
    oplot,[ensmean_lon(idx), ensmean_lon(idx-1)],$
        [ensmean_lat(idx),ensmean_lat(idx-1)],color=0,thick=7
  endif
endfor

loadct,colortableno

for ilead = 0, fcst_len, 24 do begin

  ;print,'------------- LEAD = ',ilead

  idx = ilead/fcst_output_int
  cidx = string(ilead/24,format='(i1)')

  ; --- plot dots of the appropriate color at each of the ensemble positions

  for imem = 0, nmembers-1 do begin
    if (ensfcst.ctr_lat(idx,imem) gt -999) then begin
      oplot,[ensfcst.ctr_lon(idx,imem), ensfcst.ctr_lon(idx,imem)],$
        [ensfcst.ctr_lat(idx,imem),ensfcst.ctr_lat(idx,imem)],psym=8,symsize=0.3,$
        color=colortb(cidx)
    endif
  endfor

  ; ---- now calculate ellipse (bivariate normal) that 
  ;      encompasses 90% of the probability

  xlats = [ensfcst.ctr_lat(idx,*)]
  xlons = [ensfcst.ctr_lon(idx,*)]

  a = where (xlats GE -90 AND xlats LE 90.)
  IF (a(0) ne -1 AND n_elements(a) GE 3) THEN BEGIN
    xlatmean = TOTAL(xlats(a))/FLOAT(n_elements(a))
    xlonmean = TOTAL(xlons(a))/FLOAT(n_elements(a))
    xlatdiffs = [xlats(a) - xlatmean]
    xlondiffs = [xlons(a) - xlonmean]
    Pb = fltarr(2,2)
    Pb(0,0) = total(xlondiffs^2) / FLOAT(n_elements(a)-1)
    Pb(1,1) = total(xlatdiffs^2) / FLOAT(n_elements(a)-1)
    Pb(1,0) = total(xlatdiffs*xlondiffs) / FLOAT(n_elements(a)-1)
    Pb(0,1) = Pb(1,0)
    rho = Pb(1,0) / (SQRT(pb(0,0)) * SQRT(pb(1,1)))
    sigmax = SQRT(Pb(0,0))
    sigmay = SQRT(Pb(1,1))
;
;    ; loop around a circle
;
    xcontour = fltarr(361)
    ycontour = fltarr(361)
;
    for idegree = 0, 360 do begin
      radians = FLOAT(idegree)*!pi/180.
      xstart = COS(radians)
      ystart = SIN(radians)
      for rdistance = 0,1600,1 do begin
        xloc = xstart*rdistance/40.
        yloc = ystart*rdistance/40.
        dist = sqrt(xloc^2 + yloc^2)
        fac = 1./(2.*(1-rho^2))
        prob =  exp(-1.0 * fac*($
          (xloc/sigmax)^2 + (yloc/sigmay)^2 - $
          2.*rho* (xloc/sigmax)*(yloc/sigmay)))
        if (prob lt 0.256) then begin
          xcontour(idegree) = xloc + xlonmean
          ycontour(idegree) = yloc + xlatmean
          goto, after
        endif
      endfor  ; -- rdistance
      after:
   endfor     ; -- idegree
;
;   ; ---- actually plot the contour 
;
   oplot,xcontour, ycontour, color=colortb(cidx), thick = 4 
;
;   ; ---- annotate with a number indicating the forecast lead
;
   xmax = max(xcontour,minss)
   ymax = ycontour(minss)
   cday = string(ilead/24,format='(i1)')
   print,cday,cidx
   xyouts,xmax+0.5,ymax,cday,color=colortb(cidx),charsize=1.4

  ENDIF  ; a(0) ne -1) 

  ; ----- plot a bigger dot at the location of the ensemble mean, and 
  ;       connect mean locations

  IF (a(0) ne -1) THEN BEGIN
    xlonold = xlonmean
    xlatold = xlatmean
    oplot,[xlonmean,xlonmean],[xlatmean,xlatmean],psym=8,$
      symsize=1.,color=colortb(cidx)
  ENDIF

endfor

return
end
