#!/bin/csh
set date=2011082912
while ($date <= 2011082912) 
   csh plot_ellipse.csh $date
   /bin/mv -f ellipses_${date}*gif $date
 set date=`incdate $date 12`
end
