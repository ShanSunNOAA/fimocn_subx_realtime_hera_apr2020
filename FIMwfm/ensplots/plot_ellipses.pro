pro plot_ellipses, istat, cyyyymmddhh, cbasin, cstormno, output_file, fcst_len, fcst_output_int, plotPath

; ---- controls the actual plotting of the ellipses.
;
; inputs: cyyyymmddhh - year, month, day, hour as character string
;         cbasin - WP, EP, or AL
;         cstormno - storm number identifier for the atcf file
;         outfile - the name of the output file to store data to
; outputs: none through calling list, though this program is expected
;    to create and fill 'outfile'
;
print,'begin subroutine plot_ellipses'

; ----- spawn a fortran program to read FIMens ensemble forecast
;       track data and store it to the working file "idl_fcst.dat"

tcvitalsPath = GETENV('tcvitalsPath') 
glvl_in = GETENV('glvl_in')
fcst_len_in = GETENV('fcst_len_in')
fcst_output_int = GETENV('fcst_output_int_in')
rundir_in = GETENV('rundir_in')
nmembers = GETENV('nmembers')
fimensPath = GETENV('fimensPath')
plotPath = GETENV('plotPath')
fim_bin = GETENV('fim_bin')
gribPath = GETENV('gribPath')

print,'in plot_ellipses.pro plotPath: ',plotPath
print,'in plot_ellipses.pro fimensPath: ',fimensPath
print,'in plot_ellipses.pro output_file: ',output_file

print,'calling extract_atcf'
cglvl = string(glvl_in,format='(i2)')
cfcst_len = string(fcst_len_in,format='(i3)')
cfcst_output_int = string(fcst_output_int,format='(i3)')

; print,'BEFORE 1st MKDIR plotPath: ',plotPath
FILE_MKDIR, plotPath

cmd = string(fim_bin+'/extract_atcf_FIM '+cyyyymmddhh+' '+cstormno+' '+cbasin+' '+gribPath+' '+plotPath+' '+cglvl+' '+cfcst_len+' '+cfcst_output_int)
print, 'in PLOT_ELLIPSES extract cmd: ', cmd
spawn, cmd, n_avail
; print, '************** after PLOT_ELLIPSES extract avail: ', n_avail, '****'
; print, "HELP"
; print, 'plotPath: ',plotPath
help, n_avail
infile = plotPath+'/idl_fcst.dat'
print, 'before read: infile: ',infile

read_atcf_track_data,infile,ensfcst,nmembers,istat,fcst_len_in,fcst_output_int
if (istat EQ 0) then begin
  print, 'ERROR in read_atcf_track_data - SHUTDOWN'
  goto, closedown
endif

; ----- compute the plot domain boundary and the lat/lon label locations

rlons = [ensfcst.ctr_lon(*,*)]
rlats = [ensfcst.ctr_lat(*,*)]
a = where (rlons GT -999,numa)
if (numa GE 1) THEN BEGIN
   rlon2 = rlons(a)
   a = where (rlats GT -999)
   rlat2 = rlats(a)
   rlonmax = max(rlon2)
   rlonmin = min(rlon2)
   rlatmax = max(rlat2)
   rlatmin = min(rlat2)
   print,rlonmax,rlonmin,rlatmax,rlatmin

   ilonmin = FIX((rlonmin)/5.0)*5 - 1
   ilonmax = FIX((rlonmax+5.)/5.0)*5 + 1
   ilatmin = FIX((rlatmin)/5.0)*5 - 1
   ilatmax = FIX((rlatmax+5.)/5.0)*5 + 1

   print,'before ',ilonmin,ilonmax,ilatmin,ilatmax
   dx=abs(ilonmax-ilonmin)
   dy=abs(ilatmax-ilatmin)
   xinches = 5.0
yinches = xinches/1.15
  
corr=1.15/(FLOAT(dx) / FLOAT(dy))
IF (dy*1.15 GT dx) THEN BEGIN
   print,'X'
   xmid=(ilonmin+ilonmax)/2
   dxnew=dx*corr
   ilonmin=xmid-dxnew/2
   ilonmax=xmid+dxnew/2
ENDIF
IF (dy*1.15 LT dx) THEN BEGIN
   print,'Y'
   ymid=(ilatmax+ilatmin)/2
   dynew=dy/corr
   ilatmin=ymid-dynew/2
   ilatmax=ymid+dynew/2
ENDIF


; ----- plot the data

!p.multi=0.
!p.color=0
!p.background=!d.n_colors-1
colortableno = 5
set_plot,'Z'
loadct,colortableno
colortb = [0, 160, 70, 110, 170, 50, 30, 130]

basins_c = ['E','W','C','L','B','S','P']
basins_cc = ['EP','WP','CP','AL','IO','SI','SP']

tcvitalsFile = tcvitalsPath +  "/tcvitals."+cyyyymmddhh+".txt"
read_tcvitals_data,tcvitalsFile,stormNos,stormNames
print, 'in plot_ellipses: stormNos: ',stormNos
print, 'in plot_ellipses: stormNames: ',stormNames
print, 'cstormno: ',cstormno, 'cbasin: ',cbasin
help,basins_cc
help,cbasin
index = where(strmatch(basins_cc,cbasin,/fold_case) eq 1)
help,index
print,index

print, 'index[0]: ',index[0]
if (index[0] NE '-1') then begin
   s = strtrim(string(cstormno+basins_c[index[0]]))
endif else begin
   s = strtrim(string(cstormno+cbasin))
endelse
print, "*********s: ",s	

; help,stormNos
; help,s
print, 's[0]: ',s[0]
s2=s[0]
help,s2
index = where(strmatch(stormNos,s2,/fold_case) eq 1)
; help,index
; print, "index:",index
if (index[0] NE '-1' ) then begin
;    print,'index ne -1'
   stormDir = stormNames[index[0]]+'_'+s2
endif else begin
;   print,'index ge 0'
   stormDir = s2
endelse

; print,'stormDir: ',stormDir

infile = plotPath+'/idl_fcst.dat'
print, 'infile: ',infile

plotPath=plotPath+'/'+stormDir
print,'plotPath: ',plotPath

; print,'BEFORE 2nd MKDIR plotPath: ',plotPath
FILE_MKDIR, plotPath

limit = [FLOAT(ilatmin),FLOAT(ilonmin),FLOAT(ilatmax),FLOAT(ilonmax)]
print,'limit: ',limit

output_file = plotPath+'/ellipses_'+cyyyymmddhh+'_'+cstormno+'_'+cbasin+'.ps'
print,'writing output to ',output_file
a1=findgen(17)*(!pi*2./16.)
usersym,cos(a1),sin(a1),/fill

set_post_filename,output_file,xinches,yinches

;!p.position=[.08,.07,.97,.87]
!p.position=[.08,.10,.97,.87]
title = '!5FIM-40km '+$ 
        '!cIC='+cyyyymmddhh+' for storm number '+cstormno+' in the '+cbasin+' basin'
iusa = 1
if (cbasin EQ 'WP' or cbasin EQ 'EP') then iusa=0
if (ilonmin LT 180) THEN BEGIN
   map_set,/CYLINDRICAL, 0,180,$
   title=title,/continents,charsize=0.75,usa=iusa,$
   limit=limit,/grid,latdel=5,londel=5,mlinethick=1,/hires
   xyouts,0.27,0.03,'NOAA/ESRL Global Systems Division',/normal
ENDIF ELSE BEGIN
   map_set,/CYLINDRICAL, 0,0,$
   title=title,/continents,charsize=0.75,usa=iusa,$
   limit=limit,/grid,latdel=5,londel=5,mlinethick=1,/hires
   xyouts,0.27,0.03,'NOAA/ESRL Global Systems Division',/normal
ENDELSE

loadct,8
map_continents,color=235,/fill
loadct,colortableno

plot_latlon_labels,ilonmin,ilonmax,ilatmin,ilatmax

overplot_fcstens, ensfcst, nmembers, colortableno, colortb, fcst_len, fcst_output_int
if (ilonmin LT 180) THEN BEGIN
   map_set,/CYLINDRICAL, 0,180,$
   title=title,/continents,charsize=0.75,usa=iusa,$
   limit=limit,/grid,latdel=5,londel=5,mlinethick=1,/hires,/noerase
ENDIF ELSE BEGIN
   map_set,/CYLINDRICAL, 0,0,$
   title=title,/continents,charsize=0.75,usa=iusa,$
   limit=limit,/grid,latdel=5,londel=5,mlinethick=1,/hires,/noerase
ENDELSE
device,/close
device,/portrait

;cmd = '/bin/rm '+infile
;spawn,cmd

closedown:
ENDIF

cmd = '/bin/rm '+infile
spawn,cmd

return
end
