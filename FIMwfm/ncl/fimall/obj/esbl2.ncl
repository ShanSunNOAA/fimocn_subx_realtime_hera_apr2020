  setvalues lbid
    "lbLabelAlignment" : "InteriorEdges"
    "lbLabelPosition" : "Top"
;    "lbFillColors" : (/ 52,53,35,32,29,26,0,52,53,23,20,17,14,0,52,53,11,8,5,2 /)
    "lbFillColors" : (/ 1,53,35,32,29,26,0,1,53,23,20,17,14,0,1,53,11,8,5,2 /)
    "lbLabelFont" : 21
    "lbLabelFontColor" : 1
    "lbLabelFontThicknessF" : 2
    "lbLabelStrings" : (/ "","",back2string,"","","","","","",back1string, \
                          "","","","","","",currentstring,"","" /)

;    "vpXF"      : 0.18
    "vpXF"      : vpx
    "vpYF"      : (vpy - vpheight) + 0.011
;    "vpYF"      : (vpy - vpheight)
;    "vpWidthF"  : vpwidth-0.26
    "vpWidthF"  : vpwidth
;    "vpHeightF" : 0.1
    "vpHeightF" : 0.06
  end setvalues
