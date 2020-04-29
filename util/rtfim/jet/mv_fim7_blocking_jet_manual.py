#!/usr/bin/python

import sys
import os
import time
import datetime
from datetime import timedelta
from time import gmtime, strftime
now = datetime.datetime.now()
bdates = []
# bdates.append("2013110100")
# bdates.append("2013110112")
# bdates.append("2013110200")
# bdates.append("2013110212")
# bdates.append("2013110300")
# bdates.append("2013110312")
# bdates.append("2013110400")
# bdates.append("2013110412")
# bdates.append("2013110500")
# bdates.append("2013110512")
# bdates.append("2013110600")
bdates.append("2014041000")
bdates.append("2014041012")
bdates.append("2014041100")
bdates.append("2014041112")
bdates.append("2014041200")
bdates.append("2014041212")
bdates.append("2014041300")
bdates.append("2014041312")
bdates.append("2014041400")
bdates.append("2014041412")
bdates.append("2014041500")
bdates.append("2014041512")
bdates.append("2014041600")
bdates.append("2014041612")
bdates.append("2014041700")
bdates.append("2014041712")

print str(now)
input_dir = '/pan2/projects/nrtrr/amb-verif/blocking/fim7/'
output_dir = '/data/amb/graphics/fim7_blocking_jet/'
lockfiles_dir = '/home/rtfim/lockfiles/'

for bdate in bdates :
  if os.path.isfile(input_dir+bdate+'00/fim7block.png') :
    if not os.path.isfile(lockfiles_dir+'mv_fim7_blocking_jet.lock') :
      open(lockfiles_dir+'mv_fim7_blocking_jet.lock', 'w')
#  os.utime(lockfiles_dir+'mv_hrrr_test.lock', None)
      mkdir_cmd  = 'ssh ratchet mkdir -m 777 -p '+output_dir+bdate+'/full'
      print mkdir_cmd 
      os.system(mkdir_cmd)
      sync_cmd = 'rsync -avz '+input_dir+bdate+'00/fim7block.png ratchet:' + output_dir+bdate+'/full/block_f00.png --delete --rsh=ssh'
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
      os.remove(lockfiles_dir+'mv_fim7_blocking_jet.lock')
      now2 = datetime.datetime.now()
      print str(now2)
