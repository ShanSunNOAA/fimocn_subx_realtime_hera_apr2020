  load "$NCL_HOME/obj/cref_bkg.ncl"
  intmdate = getSubString(intm,0,9)
  intminit = getSubString(intm,12,13)
  intmcurrentruntime = intmdate+" "+intminit+"Z"
  currentstring = intmcurrentruntime+" "+fctm+"hr fcst"
  if (fnameback2test .eq. 1) then
    load "$NCL_HOME/obj/cref_back2_abv5.ncl"
    setvalues cnidback2abv5
      "cnFillMode" : "AreaFill"
      "cnLevels" : (/ cl0 /)
      "cnFillColors" : (/ -1,53 /)
    end setvalues
    load "$NCL_HOME/obj/cref_back2.ncl"
    setvalues cnidback2
      "cnFillMode" : "AreaFill"
      "cnLevels" : (/ cl1,cl2,cl3,cl4 /)
      "cnFillColors" : (/ -1,35,32,29,26 /)
    end setvalues
    intmback2date = getSubString(intmback2,0,9)
    intmback2init = getSubString(intmback2,12,13)
    intmback2runtime = intmback2date+" "+intmback2init+"Z"
    back2string = intmback2runtime+" "+fctmback2+"hr fcst"
  else
    back2string = "MISSING"
  end if
  if (fnameback1test .eq. 1) then
    load "$NCL_HOME/obj/cref_back1_abv5.ncl"
    setvalues cnidback1abv5
      "cnFillMode" : "AreaFill"
      "cnLevels" : (/ cl0 /)
      "cnFillColors" : (/ -1,53 /)
    end setvalues
    load "$NCL_HOME/obj/cref_back1.ncl"
    setvalues cnidback1
      "cnFillMode" : "AreaFill"
      "cnLevels" : (/ cl1,cl2,cl3,cl4 /)
      "cnFillColors" : (/ -1,23,20,17,14 /)
    end setvalues
    intmback1date = getSubString(intmback1,0,9)
    intmback1init = getSubString(intmback1,12,13)
    intmback1runtime = intmback1date+" "+intmback1init+"Z"
    back1string = intmback1runtime+" "+fctmback1+"hr fcst"
  else
    back1string = "MISSING"
  end if
;
  load "$NCL_HOME/obj/cref_current_abv5.ncl"
  setvalues cnidcurrentabv5
    "cnFillMode" : "AreaFill"
    "cnLevels" : (/ cl0 /)
    "cnFillColors" : (/ -1,53 /)
  end setvalues
  load "$NCL_HOME/obj/cref_current.ncl"
  setvalues cnid
    "cnFillMode" : "AreaFill"
    "cnLevels" : (/ cl1,cl2,cl3,cl4 /)
    "cnFillColors" : (/ -1,11,8,5,2 /)
  end setvalues
;
  scl0 = tostring_with_format(cl0,"%4.2f")
  scl1 = tostring_with_format(cl1,"%4.2f")
  scl2 = tostring_with_format(cl2,"%4.2f")
  scl3 = tostring_with_format(cl3,"%4.2f")
  scl4 = tostring_with_format(cl4,"%4.2f")
  print("scl0 = "+scl0)
  print("scl1 = "+scl1)
  print("scl2 = "+scl2)
  print("scl3 = "+scl3)
  print("scl4 = "+scl4)
;
  lbid = create "LabelBar" labelBarClass wid
      "lbOrientation" : "Horizontal"
      "lbAutoManage" : False
      "lbLabelAutoStride" : False
      "lbLabelFont" : 4
      "lbLabelFontHeightF" : .011
      "lbPerimOn" : False
      "lbBoxCount" : 20
      "lbBoxLinesOn" : False
      "lbLabelFont" : 21
      "lbLabelFontColor" : 1
      "lbLabelFontThicknessF" : 2
;      "lbFillColors" : (/ 52,53,35,32,29,26,0,52,53,23,20,17,14,0,52,53,11,8,5,2 /)
      "lbFillColors" : (/ 1,53,35,32,29,26,0,1,53,23,20,17,14,0,1,53,11,8,5,2 /)
      "lbMonoFillPattern" : True
      "lbLabelAlignment" : "InteriorEdges"
      "lbLabelPosition" : "Bottom"
      "lbLabelStrings" : (/ scl0,scl1,scl2,scl3,scl4,"","",scl0,scl1,scl2,scl3,scl4, \
                            "","",scl0,scl1,scl2,scl3,scl4 /)
;      "vpXF"      : 0.18
      "vpXF"      : vpx
      "vpYF"      : (vpy - vpheight) - 0.008
;      "vpYF"      : (vpy - vpheight)
;      "vpWidthF"  : vpwidth-0.26
      "vpWidthF"  : vpwidth
;      "vpHeightF" : 0.1
      "vpHeightF" : 0.06
  end create
;
; be sure to include 2nd part of lbid (esbl2.ncl) after drawing 1st part
;
