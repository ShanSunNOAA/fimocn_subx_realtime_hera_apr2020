check = ""
check = systemfunc("ls back1_run.txt")
fnameback1test = 0
fnameback2test = 0
if (.not.ismissing(check.eq."back1_run.txt")) then
   fnameback1 = asciiread("back1_run.txt",-1,"string")
   back1exists = isfilepresent(fnameback1)
   if (back1exists) then
     print("reading "+fnameback1)
     fnameback1test = 1
     grb_file_back1 = addfile(fnameback1,"r")
   else
     print("error - missing file "+fnameback1)
   end if
end if
check2 = ""
check2 = systemfunc("ls back2_run.txt")
if (.not.ismissing(check2.eq."back2_run.txt")) then
   fnameback2 = asciiread("back2_run.txt",-1,"string")
   back2exists = isfilepresent(fnameback2)
   if (back2exists) then
     print("reading "+fnameback2)
     fnameback2test = 1
     grb_file_back2 = addfile(fnameback2,"r")
   else
     print("error - missing file "+fnameback2)
   end if
end if
;
if (fnameback1test .eq. 1) then
  print("fnameback1 ="+fnameback1)
  intmback1 = grb_file_back1->$T2M$@initial_time
  fctmback1 = grb_file_back1->$T2M$@forecast_time
else
  intmback1 = -999
  fctmback1 = -999
end if
if (fnameback2test .eq. 1) then
  print("fnameback2 = "+fnameback2)
  intmback2 = grb_file_back2->$T2M$@initial_time
  fctmback2 = grb_file_back2->$T2M$@forecast_time
else
  intmback2 = -999
  fctmback2 = -999
end if
;

