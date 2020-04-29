PRO plot_latlon_labels,ilonmin,ilonmax,ilatmin,ilatmax


ilonmin2 = FIX((ilonmin)/5.0)*5 +5
ilonmax2 = FIX((ilonmax+5.)/5.0)*5-5 
ilatmin2 = FIX((ilatmin)/5.0)*5 +5
ilatmax2 = FIX((ilatmax+5.)/5.0)*5-5 

rlat = FLOAT(ilatmin)-0.05*FLOAT(ilatmax-ilatmin)
rlon = FLOAT(ilonmin)-0.07*FLOAT(ilonmax-ilonmin)

for ilon = ilonmin2,ilonmax2,5 do begin
  if (ilon gt 180) then begin
    clon = string(ABS(ilon-360),format='(i3)')+'W'
  endif else begin
    clon = string(ilon,format='(i3)')+'E'
  endelse
  xyouts,FLOAT(ilon),rlat,clon,align=0.5,charsize=0.7
endfor

for ilat = ilatmin2, ilatmax2, 5 do begin
  clat = string(ilat,format='(i2)')+'N'
  xyouts,rlon,FLOAT(ilat)-0.01*FLOAT(ilonmax-ilonmin),clat,align=0,charsize=0.7
endfor

return
end
