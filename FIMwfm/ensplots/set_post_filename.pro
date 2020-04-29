pro set_post_filename,outfile,xinch,yinch

set_plot,'ps'
!p.font=0
device,/close
device,/portrait
device,file=outfile,bits=8,/inches,encapsulated=0,$
  xs=xinch,ys=yinch,xoff=1.0,yoff=1.0
!x.thick=2
!y.thick=2
!p.charthick=2
!p.font=-1
device,/color

return
end
