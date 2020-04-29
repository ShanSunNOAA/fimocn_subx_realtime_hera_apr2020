#!/usr/bin/python

import sys
import os
import time
import datetime
from datetime import timedelta
from time import gmtime, strftime
now = datetime.datetime.now()
bdates = []
bdates.append("2013080512")
bdates.append("2013080600")
bdates.append("2013080612")
bdates.append("2013080700")
bdates.append("2013080712")
bdates.append("2013080800")
bdates.append("2013080812")
bdates.append("2013080900")
bdates.append("2013080912")
bdates.append("2013081000")
bdates.append("2013081012")
bdates.append("2013081100")
bdates.append("2013081112")
bdates.append("2013081200")
bdates.append("2013081212")
bdates.append("2013081300")
bdates.append("2013081312")
bdates.append("2013081400")
bdates.append("2013081412")
bdates.append("2013081500")
bdates.append("2013081512")
bdates.append("2013081600")
bdates.append("2013081612")
bdates.append("2013081700")
bdates.append("2013081712")
bdates.append("2013081800")
bdates.append("2013081812")
bdates.append("2013081900")
bdates.append("2013081912")
bdates.append("2013082000")
bdates.append("2013082012")
bdates.append("2013082100")
bdates.append("2013082112")
bdates.append("2013082200")
bdates.append("2013082212")
bdates.append("2013082300")
bdates.append("2013082312")
bdates.append("2013082400")
bdates.append("2013082412")
bdates.append("2013082500")
bdates.append("2013082512")
bdates.append("2013082600")
bdates.append("2013082612")

print str(now)
input_dir = '/pan2/projects/nrtrr/amb-verif/blocking/gfs/'
output_dir = '/data/amb/graphics/gfs_blocking_jet/'
lockfiles_dir = '/home/rtfim/lockfiles/'

for bdate in bdates :
  if os.path.isfile(input_dir+bdate+'00/gfsblock.png') :
    if not os.path.isfile(lockfiles_dir+'mv_gfs_blocking_jet.lock') :
      open(lockfiles_dir+'mv_gfs_blocking_jet.lock', 'w')
#  os.utime(lockfiles_dir+'mv_hrrr_test.lock', None)
      mkdir_cmd  = 'ssh ratchet mkdir -m 777 -p '+output_dir+bdate+'/full'
      print mkdir_cmd 
      os.system(mkdir_cmd)
      sync_cmd = 'rsync -avz '+input_dir+bdate+'00/gfsblock.png ratchet:' + output_dir+bdate+'/full/block_f00.png --delete --rsh=ssh'
      print sync_cmd
      os.system(sync_cmd)
      makeimageslist_cmd = 'ssh ratchet.fsl.noaa.gov \'(cd ' + output_dir+bdate+'/full; /bin/ls -1 > images.txt; chmod 777 images.txt)\''
      print makeimageslist_cmd
      os.system(makeimageslist_cmd)
      chm_cmd  = 'ssh ratchet.fsl.noaa.gov chmod -R a=rwx '+output_dir+bdate
      print chm_cmd
      os.system(chm_cmd)
      makedatelist_cmd = 'ssh ratchet.fsl.noaa.gov \'(cd ' + output_dir+'; /bin/ls -1 > dates.txt; chmod 777 dates.txt)\''
      print makedatelist_cmd
      os.system(makedatelist_cmd)
      os.remove(lockfiles_dir+'mv_gfs_blocking_jet.lock')
      now2 = datetime.datetime.now()
      print str(now2)
