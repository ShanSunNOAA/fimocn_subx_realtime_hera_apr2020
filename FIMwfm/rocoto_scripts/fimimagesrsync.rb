#!/usr/bin/ruby

RSYNC="/usr/bin/rsync -vrpgoDu -e ssh"
TAR="tar --append --file="
BASEDIR="/pan2/projects/fim-njet/FIM/FIMrun_sjet_p"
OUTDIR="/rt/fim/images"
DEBUG=true
KEEP=2

# Get most recent run
now=Time.now.to_i
last_run=Time.at(now - now % (3600*12))

puts "*******************************"
0.upto(KEEP-1) { |i|

  yyyymmddhhmm=(last_run - (i*3600*12)).strftime("%Y%m%d%H%M")
  hh=(last_run - (i*3600*12)).strftime("%H").to_i

  # tar and sync image files 
  bdir=yyyymmddhhmm[0..9]
  if not File.exists?("#{OUTDIR}/#{bdir}.tar.gz")
    puts "*IMAGES*"
    Dir.chdir("#{BASEDIR}")
    Dir["#{bdir}/ncl_C/*/files.zip"].each { |dir|
      puts "#{dir}"
      tarcmd="#{TAR}#{BASEDIR}/#{yyyymmddhhmm}.tar #{dir}"
      output=`#{tarcmd} 2>&1`
      error=$?.exitstatus
      puts output if DEBUG
      if error != 0
        puts "ERROR: #{tarcmd} failed!  Exit status=#{error}"
        puts output
      end
    }

    gzipcmd="gzip #{BASEDIR}/#{yyyymmddhhmm}.tar"
    output=`#{gzipcmd} 2>&1`
    error=$?.exitstatus
    puts output if DEBUG
    if error != 0
      puts "ERROR: #{gzipcmd} failed!  Exit status=#{error}"
      puts output
    end
    rsynccmd="#{RSYNC} #{BASEDIR}/#{yyyymmddhhmm}.tar.gz #{OUTDIR}"
    output=`#{rsynccmd} 2>&1`
    error=$?.exitstatus
    puts output if DEBUG
    if error != 0
      puts "ERROR: #{rsynccmd} failed!  Exit status=#{error}"
      puts output
    end
    File.delete("#{BASEDIR}/#{yyyymmddhhmm}.tar.gz")
  end

}
