#!/usr/bin/python

import sys
import os
import time
import datetime
from datetime import timedelta
from time import gmtime, strftime

now = datetime.datetime.now()
print str(now)
input_dir = '/home/rtfim/'
output_dir = '/w3/rapb/fim/from_jet/'

if os.path.isfile(input_dir+'statusAll') :
  sync_cmd = 'rsync -avz '+input_dir+'statusAll ratchet.fsl.noaa.gov:'+output_dir+'statusAll_alt --rsh=ssh'
  print sync_cmd
  os.system(sync_cmd)
  chm_cmd  = 'ssh ratchet.fsl.noaa.gov chmod -R a=rwx '+output_dir+'statusAll_alt'
  print chm_cmd
  os.system(chm_cmd)

  now2 = datetime.datetime.now()
  print str(now2)

