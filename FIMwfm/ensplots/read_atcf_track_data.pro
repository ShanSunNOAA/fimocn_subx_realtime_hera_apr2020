PRO read_atcf_track_data,infile,ensfcst,nmembers,istat,fcst_len, fcst_output_int

; ----- read forecast track data from file

istat = 0
print,'beginning subroutine read_atcf_track_data'
print,'in read_atcf_track_data: nmembers:', nmembers

ileadidx_in = 0L
imem_in = 0L
ctr_lat_in = 0.0
ctr_lon_in = 0.0
centralpressure_in = 0.0
windspeed_in = 0.0
iflead_in = 0L
nleads = 29
ensfcst = {$
  iflead:replicate(-999999L, nleads),$
  ctr_lat: replicate(-9999.99, nleads, nmembers),$
  ctr_lon: replicate(-9999.99, nleads, nmembers),$
  centralpressure: replicate(-9999.99, nleads, nmembers),$
  windspeed: replicate(-9999.99, nleads, nmembers)}

close,1
PRINT, '*** infile: ',infile
openr,1, infile
FOR ilead  = 0, fcst_len, fcst_output_int do begin
  ileadidx = ilead/fcst_output_int
  ;PRINT, '*** ileadidx: ',ileadidx
  readf, 1, ileadidx_in, format='(i3)'

  IF (ileadidx_in EQ '666') THEN BEGIN
    PRINT,'no data available due to corrupt atcf files.'
    PRINT,'sorry, ellipses will not be generated. Terminating'
    GOTO, after
  ENDIF

  istat = 1
  FOR imem = 0, nmembers-1 do begin
    print, 'imem: ',imem
    readf,1,imem_in,iflead_in, $
      ctr_lat_in, ctr_lon_in, centralpressure_in,$
      windspeed_in,$
      format='(i2,1x,i3,1x,4(f8.2,1x))'
      print, 'imem_in:',imem_in, ' iflead_in: ',iflead_in, ' ctr_lat_in: ',ctr_lat_in, $
             'ctr_lon_in: ', ctr_lon_in, 'centrappre: ',centralpressure_in, 'windspeed: ',windspeed_in

    ensfcst.iflead(ileadidx) = iflead_in
    ensfcst.ctr_lat(ileadidx,imem) = ctr_lat_in
    ensfcst.ctr_lon(ileadidx,imem) = ctr_lon_in
    ensfcst.centralpressure(ileadidx,imem) = centralpressure_in
    ensfcst.windspeed(ileadidx,imem) = windspeed_in
  ENDFOR
ENDFOR 

print,'subroutine read_atcf_track_data:  data read in'

after:

return
end
