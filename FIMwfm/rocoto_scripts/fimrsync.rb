#!/usr/bin/ruby

RSYNC="/usr/bin/rsync -vrpgoDu -e ssh"
BASEDIR="/pan2/projects/fim-njet/FIM/FIMrun_sjet_p"
DEBUG=true
KEEP=2

# Get most recent run
now=Time.now.to_i
last_run=Time.at(now - now % (3600*12))

puts "*******************************"
0.upto(KEEP-1) { |i|

  yyyymmddhhmm=(last_run - (i*3600*12)).strftime("%Y%m%d%H%M")
  hh=(last_run - (i*3600*12)).strftime("%H").to_i

  # rsync grib1 files 
  puts "*GRIB1*"
  Dir["#{BASEDIR}/fim_[0-9]*_[0-9]*_[0-9]*_#{yyyymmddhhmm}/post_C/fim/NAT/grib1"].each { |dir|
    cmd="#{RSYNC}  #{dir}/[0-9]*[0-9] /rt/fim/nat/grib1"
    output=`#{cmd} 2>&1`
    error=$?.exitstatus
    puts output if DEBUG
    if error != 0
      puts "ERROR: #{cmd} failed!  Exit status=#{error}"
      puts output
    end

  }

  # Only rsync grib2 files that are valid at 00Z and 12Z
  puts "*GRIB2*"
  Dir["#{BASEDIR}/fim_[0-9]*_[0-9]*_[0-9]*_#{yyyymmddhhmm}/post_C/fim/NAT/grib2"].each { |dir|
    cmd="#{RSYNC}  #{dir}/[0-9]*[0-9] /rt/fim/nat/grib2"
    output=`#{cmd} 2>&1`
    error=$?.exitstatus
    puts output if DEBUG
    if error != 0
      puts "ERROR: #{cmd} failed!  Exit status=#{error}"
      puts output
      next
    end

  }

  # sync tracker files
  #   NOTE:  to use wildcard in directory name, need to use Dir.glob
  puts "*TRACKER*"
  puts "#{BASEDIR}/fim_8_64_*_#{yyyymmddhhmm}/tracker_C/240"
  if File.exists?("#{BASEDIR}/fim_8_64_240_#{yyyymmddhhmm}/tracker_C/240")
    cmd="#{RSYNC} #{BASEDIR}/fim_8_64_240_#{yyyymmddhhmm}/tracker_C/240/*track* /rt/fim/tracker"
    output=`#{cmd} 2>&1`
    error=$?.exitstatus
    puts output if DEBUG
    if error != 0
      puts "ERROR: #{cmd} failed!  Exit status=#{error}"
      puts output
    end
  end

}
