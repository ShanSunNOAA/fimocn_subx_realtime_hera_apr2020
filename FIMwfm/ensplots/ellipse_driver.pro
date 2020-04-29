pro ellipse_driver
;
; idl program to generate ensemble track forecasts overlaid with
; ellipses.  these ellipses come from a bivariate normal distribution
; fit to the ensemble forecasts, and the ellipses contain 90% of the 
; fitted probability.
;
; assumption here is that this driver program will read a file 
; named 'ensemble-storm-data.txt' and each line will contain data for
; a particular storm.  The line of text will have the date/time,
; the basin (AL, EP or WP), and the storm number, all separated 
; by a space.  So, for example, the file may contain the following:
;
; 2010090112 AL 06
; 2010090112 AL 07
; 2010090112 AL 08
;
; The output will be a postscript file for each storm that can be
; converted to a gif file, or other format.  The filename convention
; will be YYYYMMDDHH-basin-stormno.ps, e.g., 2010090112-AL-06.ps
;

print,'idl program ellipse_driver entered....'

; ---- read in the forecast information

tcvitalsPath = GETENV('tcvitalsPath')
glvl_in = GETENV('glvl_in')
fcst_len_in = GETENV('fcst_len_in')
fcst_output_int_in = GETENV('fcst_output_int_in')
; fimensPath = GETENV('fimensPath')
plotPath = GETENV('plotPath')
nmembers = GETENV('nmembers')
fim_bin = GETENV('fim_bin')

print,'in ellipse_driver.pro glvl_in = ',glvl_in
print,'in ellipse_driver.pro fcst_len_in = ',fcst_len_in
print,'in ellipse_driver.pro fcst_output_int_in = ',fcst_output_int_in
print,'in ellipse_driver.pro plotPath = ',plotPath
print,'start'
;input_file = '/Users/thamill/hfip/ellipse/ensemble-storm-data.txt'
input_file = plotPath+'/ensemble-storm-data.txt'
print,'input file = ',input_file

files = file_search(input_file)
print,'file search = ',files
if (files EQ '') then begin
   print,'ellipse_driver.pro could not find the expected input file ',input_file
   print,'stopping.'
   stop
endif else begin
   close,8
   print,' reading ', input_file
   openr,8,input_file
   while  not eof(8) do begin
      linein = ''
      readf,8,linein
      print,' ---------- PROCESSING  ------- ', linein
      stringparts = STR_SEP(linein,' ')
      cyyyymmddhh = stringparts(0)
      cbasin = stringparts(1)
      cstormno = stringparts(2)
      cstormcode = stringparts(3)
      ;outfile = '/Users/thamill/hfip/ellipse/ellipses_'+cyyyymmddhh+'-'+cbasin+'-'+cstormno+'.txt'
      outfile = plotPath+'/'+cbasin+'-'+cstormno+'.txt'
  
      ; ---- make the plot

      ; plot_ellipses, istat, cyyyymmddhh, cbasin, cstormno, outfile, fcst_len_in, fcst_output_int_in, fim_bin
      plot_ellipses, istat, cyyyymmddhh, cbasin, cstormno, outfile, fcst_len_in, fcst_output_int_in, newPlotPath
      if (istat EQ 1) then $
         print,outfile,' created with ellipse data plotting information ***** '
         print,'NEW PLOT PATH: ',newPlotPath
       ;outfilepng = newPlotPath+'/'+cbasin+'-'+cstormcode+'.png'
       outfilepng = newPlotPath+'/ellipses.png'
       cmd = 'convert -density 150x150 -background white '+outfile+' '+outfilepng
      spawn, cmd
      print,'PNG file ',outfilepng,' also created.'
      cmd = '/bin/rm '+outfile
      spawn, cmd
   endwhile
   close, 8
   cmd = '/bin/rm '+input_file
   spawn, cmd
   print,'idl program ellipse_driver finished'
endelse

end
