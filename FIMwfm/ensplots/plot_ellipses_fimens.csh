#!/bin/csh
set python=/lfs2/projects/fim/whitaker/bin/python
set date=$analdate
cd /lfs2/projects/gfsenkf/hurrplots/fimens
mkdir ${date}_tst
csh plot_ellipse.csh $date
/bin/mv -f ellipses_${date}*gif ${date}_tst
/bin/rm -f ellipses_${date}*ps
cd /lfs2/projects/gfsenkf/hurrplots/mme
mkdir ${date}_tst
csh plot_ellipse.csh $date
/bin/mv -f ellipses_${date}*gif ${date}_tst
/bin/rm -f ellipses_${date}*ps
