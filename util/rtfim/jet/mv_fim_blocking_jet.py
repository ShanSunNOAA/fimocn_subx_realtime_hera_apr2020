#!/usr/bin/python

import sys
import os
import time
import datetime
from datetime import timedelta
from time import gmtime, strftime
now = datetime.datetime.now()
back15 = now - timedelta(hours=15)
b15 = back15.strftime("%Y%m%d%H")
print b15

print str(now)
input_dir = '/pan2/projects/nrtrr/amb-verif/blocking/fim/'
output_dir = '/data/amb/graphics/fim_blocking_jet/'
lockfiles_dir = '/home/rtfim/lockfiles/'

if os.path.isfile(input_dir+b15+'00/fim8block.png') :
  if not os.path.isfile(lockfiles_dir+'mv_fim_blocking_jet.lock') :
    open(lockfiles_dir+'mv_fim_blocking_jet.lock', 'w')
#  os.utime(lockfiles_dir+'mv_hrrr_test.lock', None)
    mkdir_cmd  = 'ssh ratchet mkdir -m 777 -p '+output_dir+b15+'/full'
    print mkdir_cmd 
    os.system(mkdir_cmd)
    sync_cmd = 'rsync -avz '+input_dir+b15+'00/fim8block.png ratchet:' + output_dir+b15+'/full/block_f00.png --delete --rsh=ssh'
    print sync_cmd
    os.system(sync_cmd)
    makeimageslist_cmd = 'ssh ratchet.fsl.noaa.gov \'(cd ' + output_dir+b15+'/full; /bin/ls -1 > images.txt; chmod 777 images.txt)\''
    print makeimageslist_cmd
    os.system(makeimageslist_cmd)
    chm_cmd  = 'ssh ratchet.fsl.noaa.gov chmod -R a=rwx '+output_dir+b15
    print chm_cmd
    os.system(chm_cmd)
    makedatelist_cmd = 'ssh ratchet.fsl.noaa.gov \'(cd ' + output_dir+'; /bin/ls -1 > dates.txt; chmod 777 dates.txt)\''
    print makedatelist_cmd
    os.system(makedatelist_cmd)
    os.remove(lockfiles_dir+'mv_fim_blocking_jet.lock')
    now2 = datetime.datetime.now()
    print str(now2)
