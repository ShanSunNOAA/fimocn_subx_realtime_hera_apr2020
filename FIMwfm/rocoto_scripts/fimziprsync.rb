#!/usr/bin/ruby

BASEDIR=ARGV[0]  # directory from which images will be transferred
DSKEY=ARGV[1]    # directory name on OUTPUTDIR (e.g. "fim_jet")
DOMAIN=ARGV[2]   # domain directory name below DSKEY
OUTHOST=ARGV[3]  # host machine address for transfer (e.g. "clank.fsl.noaa.gov")
INTERVAL=ARGV[4] # run interval (not forecast interval)
KEEP=ARGV[5]     # how many runs to look back
PTYPE=ARGV[6]    # set this to "ncl_C" or "ncldiff_C"
# RSYNC="/usr/bin/rsync -av -e ssh"
RSYNC="/usr/bin/rsync -a"
OUTPUTDIR="/data/amb/graphics"

DEBUG=true

# Get most recent runs (set by KEEP)
interv=INTERVAL.to_i
kp=KEEP.to_i
now=Time.now.to_i
last_run=Time.at(now - now % (3600*interv))

puts "now: #{now}"
puts "last_run #{last_run}"

time1 = Time.new
puts "current time1 : " + time1.inspect

if FileTest.exists?("#{BASEDIR}/.testlock#{DOMAIN}") then
   fileMtime=File.mtime("#{BASEDIR}/.testlock#{DOMAIN}").to_i
   overwriteTime=fileMtime+(2*3600)
   puts "fileMtime #{fileMtime}"
   puts "overwriteTime #{overwriteTime}"
end

0.upto(kp-1) { |i|

  yyyymmddhh=(last_run - (i*3600*interv)).strftime("%Y%m%d%H")
  yyyymmddhhmm=(last_run - (i*3600*interv)).strftime("%Y%m%d%H%M")
  puts "yyyymmddhhmm #{yyyymmddhhmm}"
  hh=(last_run - (i*3600*interv)).strftime("%H").to_i
  puts "hh #{hh}"

  Dir["#{BASEDIR}/fim_[0-9]*_[0-9]*_[0-9]*_#{yyyymmddhhmm}/#{PTYPE}/*"].each { |dir|
    if Dir["#{dir}/*"].empty? then
       next
    end
    if FileTest.exists?("#{dir}/files.zip") == false then
       next
    end
    if FileTest.exists?("#{dir}/.testlock#{DOMAIN}") == false then
      `touch #{dir}/.testlock#{DOMAIN}`
      dirs= dir.split("/")
      dName=dirs.at(-1)
      if dName == DOMAIN
        puts "dName #{dName} *** \n"

        cmd="ssh -v #{OUTHOST}  mkdir -p #{OUTPUTDIR}/#{DSKEY}/#{yyyymmddhh}/#{dName}"
        puts "cmd #{cmd} *** \n"
        output=`#{cmd} 2>&1`
        error=$?.exitstatus
        puts output if DEBUG

        if error != 0
          puts "ERROR: #{cmd} failed!  Exit status=#{error}"
          puts output
          `rm -f #{dir}/.testlock#{DOMAIN}`
          puts "exiting ...."
          Process.exit
        end

        cmd="ssh -v #{OUTHOST}  'ls -1 #{OUTPUTDIR}/#{DSKEY}/ > #{OUTPUTDIR}/#{DSKEY}/dates.txt; chmod 775 #{OUTPUTDIR}/#{DSKEY}/dates.txt'"
        puts "cmd #{cmd} *** \n"
        output=`#{cmd} 2>&1`
        error=$?.exitstatus
        puts output if DEBUG

        if error != 0
          puts "ERROR: #{cmd} failed!  Exit status=#{error}"
          puts output
          `rm -f #{dir}/.testlock#{DOMAIN}`
          puts "exiting ...."
          Process.exit
        end

        time2 = Time.new
        puts "current time2 : " + time2.inspect
        # thisdir=Dir.pwd
        # cmd="cd #{dir}; zip -n .png files.zip * -i \*.png; cd #{thisdir}"
        # puts "cmd #{cmd} **** \n"
        # output=`#{cmd} 2>&1`
        # error=$?.exitstatus
        # puts output if DEBUG
        # time3 = Time.new
        # puts "current time3 : " + time3.inspect
        # cmd="#{RSYNC}  #{dir}/files.zip #{OUTHOST}:#{OUTPUTDIR}/#{DSKEY}/#{yyyymmddhh}/#{dName}"
        cmd="#{RSYNC}  #{dir}/files.zip #{OUTHOST}:#{OUTPUTDIR}/#{DSKEY}/#{yyyymmddhh}/#{dName} --rsh=ssh"
        puts "cmd #{cmd} **** \n"
        output=`#{cmd} 2>&1`
        error=$?.exitstatus
        puts output if DEBUG
        time4 = Time.new
        puts "current time4 : " + time4.inspect

        if error != 0
          puts "ERROR: #{cmd} failed!  Exit status=#{error}"
          puts output
          `rm -f #{dir}/.testlock#{DOMAIN}`
          puts "exiting ...."
          Process.exit
        end
      end
    else
      puts "#{dir}/.testlock#{DOMAIN} exists"
      puts "exiting ...."
      Process.exit
    end
    `rm -f #{dir}/.testlock#{DOMAIN}`
  }
}
