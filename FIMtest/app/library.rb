module Library

  # REQUIRED METHODS (CALLED BY DRIVER)

  def lib_build(env,prepkit)

    # When running inside a suite, env.qmon will already point to an initialized
    # Qmon object. If not, this is a single-run invocation: Create a Qmon object
    # and make note that we need to shut it down later.

    unless env.qmon
      env.qmon=Qmon.new(60)
      env.must_stop_qmon=true
    end

    logd "*   (Also see #{env.build.log})"
    cmd="cd #{env.build.builddir} && qsub #{env.build.script}"
    queued=false
    output=[]

    # Try twice to queue the build job, a sad but necessary response to issues
    # with a flaky qsub (on jet, IIRC) that randomly fails.

    2.times do |n|
      props={:msg=>"Job submission failed, retrying...",:die=>false}
      output,status=Thread.exclusive { ext(cmd,props) }
      if status==0
        queued=true
        break
      end
    end

    # Make a 3rd (and final) attempt to queue the job, if necessary.

    unless queued
      props={:msg=>"Job submission failed, see #{logfile}",:die=>true}
      output,status=Thread.exclusive { ext(cmd,props) }
    end

    # Collect and verify the job id reported by qsub.

    re=Regexp.new('^(\d+)\..*')
    jobid=nil
    output.each do |e|
      m=re.match(e)
      jobid=m[1] if m
    end
    die "Job ID not found in batch-submit output" unless jobid

    # Register this job id for deletion in case of a test-suite failure.

    job_activate(jobid,self)

    # Monitor the job until completion.

    jobstatus=lib_wait_for_job(env,jobid)

    # De-register the job from the list of active jobs.

    job_deactivate(jobid)

    # Check for job success and return result.

    if jobstatus and not jobstatus==0
      die "Build #{env.build.name} failed. "+
        "See #{env.build.log} (if it exists) for details."
    end

    File.join(env.build.builddir,"FIMrun#{lib_suffix(env)}")
  end

  def lib_build_post(env,postkit)
    postkit
  end

  def lib_build_prep_jet(env)
    lib_build_prep_std(env)
    script=<<SCRIPT
#!/bin/sh
#PBS -A fim
#PBS -d #{env.build.builddir}
#PBS -j oe
#PBS -l procs=16
#PBS -l partition=vjet
#PBS -l walltime=00:20:00
#PBS -o #{env.build.log}
#PBS -q batch
#PBS -W umask=022
#{env.build.copy_src_cmd}
#{env.build.copy_run_cmd}
cd #{env.build.src}
./makefim #{env.build.args}
SCRIPT
    File.open(env.build.script,'w') { |f| f.write(script) }
    FileUtils.chmod(0755,env.build.script)
    logd "Wrote batch-build script #{env.build.script}"
    nil
  end

  def lib_build_prep_std(env)

    # Some setup.

    makefimargs=env.build.args.squeeze
    buildname=makefimargs.gsub(/  */,'_')
    env.build.name=buildname
    logd "Set build name: #{buildname}"

    # Create build location.

    builddir=(env.build.builddir=env.build.ddts_root)
    unless Dir.exist?(builddir)
      FileUtils.mkdir_p(builddir)
      logd "Made directory: #{builddir}"
    end

    # Use knowledge of this app's relationship to model directory hierarchy.

    topdir=File.expand_path("..")

    # Prepare copy commands.

    src='FIMsrc'
    srcdir=valid_dir(File.join(topdir,src))
    dstdir=builddir
    env.build.copy_src_cmd="rsync -a #{srcdir} #{dstdir}"

    run='FIMrun'
    srcdir=valid_dir(File.join(topdir,run))
    dstdir=File.join(builddir,run)
    env.build.copy_run_cmd="rsync -a --exclude '*/' #{srcdir}/ #{dstdir}/"

    # More setup.

    env.build.log=File.join(builddir,"build.log")
    env.build.script=File.join(builddir,"build.#{buildname}")
    env.build.src=File.join(builddir,src)

    nil
  end

  def lib_build_prep_theia(env)
    lib_build_prep_std(env)
    script=<<SCRIPT
#!/bin/sh
#PBS -A gsd-fv3-test
#PBS -d #{env.build.builddir}
#PBS -j oe
#PBS -l procs=24
#PBS -l walltime=00:20:00
#PBS -o #{env.build.log}
#PBS -q batch
#PBS -W umask=022
#{env.build.copy_src_cmd}
#{env.build.copy_run_cmd}
cd #{env.build.src}
./makefim #{env.build.args}
SCRIPT
    File.open(env.build.script,'w') { |f| f.write(script) }
    FileUtils.chmod(0755,env.build.script)
    logd "Wrote batch-build script #{env.build.script}"
    nil
  end

  def lib_data_jet(env)
    lib_link_data("/lfs2/projects/fim/test-suite-data/ocean/2016-02-11")
  end

  def lib_data_theia(env)
    lib_link_data("/scratch3/BMC/fim/test-suite-data/ocean/2016-02-11")
  end

  def lib_outfiles(env,path)

    # See README for a description of what this method returns.

    if env.run.namelists.oceannamelist.coupled==".true."       # Coupled atm/ocn
      restrs=['(.*/)(fim/fim_out_.*)','(.*/)(fim/ocn_out_.*)','(.*/)(post/[0-9]+)']
    else
      restrs=['(.*/)(fim/fim_out_.*)','(.*/)(post/[0-9]+)']
    end
    res=restrs.map { |e| Regexp.new(e) }
    outs=[]
    cmd="find -L #{path} -type f"
    props={:msg=>"Error executing: #{cmd}",:out=>false}
    output,status=ext(cmd,props)
    allfiles=output
    allfiles.each do |e|
      res.each do |re|
        m=re.match(e)
        unless m.nil?
          outs << [m[1],m[2]]
          break
        end
      end
    end
    outs
  end

  def lib_queue_del_cmd(env)
    'qdel'
  end

  def lib_run_check(env,postkit)
    fim_stdout,rundir=postkit
    fim_ok=job_check(fim_stdout,lib_re_str_success)
    logi "Run failed: fim stdout is #{fim_stdout}" unless fim_ok
    pop_ok=true
    unless env.run.ignore_pop
      fimdir=File.dirname(fim_stdout)
      pop_stdout=File.expand_path(File.join(fimdir,"..","post","stdout"))
      pop_ok=job_check(pop_stdout,"pop completed successfully")
      logi "Run failed: pop stdout is #{pop_stdout}" unless pop_ok
    end
    (fim_ok and pop_ok)?(rundir):(nil)
  end

  def lib_run_check_rocoto(env,postkit)
    logfile,subdir=postkit
    fim_ok=job_check(logfile,"This cycle is complete: Success")
    logi "Run failed: Rocoto log is #{logfile}" unless fim_ok
    (fim_ok)?(subdir):(nil)
  end

  def lib_run_compare_var(env,prepkit)

    # Compare_Var runs are the same as standard runs, except that their output
    # directory name is modified by the external run automation. Reproduce and
    # return that directory name here.

    rundir=prepkit
    lib_run_std(env,rundir)
    nls=env.run.namelists
    g=nls.cntlnamelist.glvl
    k=nls.cntlnamelist.nvl
    pes1=env.run.pes1
    pes2=env.run.pes2
    stdout=File.join(rundir,"fim#{g}_#{k}_cv.#{pes1}.vs.#{pes2}","fim","stdout")
    [stdout,rundir]
  end

  def lib_run_post(env,runkit)
    env.qmon.stop if env.must_stop_qmon
    runkit
  end

  def lib_run_prep_compare_var(env)

    # In addition to doing the standard run prep, set up SMSnamelist for
    # Compare_Var operation.

    rundir=lib_run_prep_std(env)
    totalpes=env.run.namelists.queuenamelist.computetasks.to_i
    pes1=1
    pes2=totalpes-pes1
    smsenv=OpenStruct.new
    smsenv.smsnamelist=OpenStruct.new
    smsnl=smsenv.smsnamelist
    smsnl.compare_var_on=YAML_Unquoted.new('.true.')
    smsnl.compare_var_ntasks_1=pes1
    smsnl.compare_var_ntasks_2=pes2
    smsnlfile=valid_file(File.join(rundir,"SMSnamelist"))
    lib_mod_namelist_file(env,smsnlfile,smsenv)
    env.run.pes1=pes1
    env.run.pes2=pes2
    rundir
  end

  def lib_run_prep_enkf(env)

    # Like a restart run, the EnKF run happens in two parts. Do the standard run
    # prep first, then set up a second namelist file for the second part of the
    # EnKF run.

    rundir=lib_run_prep_std(env)
    nlfile1=File.join(rundir,"FIMnamelist")
    nlfile2=File.join(rundir,"fimenkf_second_run")
    die "ERROR: #{nlfile2} already exists" if File.exist?(nlfile2)
    FileUtils.cp(nlfile1,nlfile2)
    logd "Copied #{nlfile1} to #{nlfile2}"
    datadir=valid_dir(File.join(tmp_dir,"data"))
    nl=env.run.namelists.modelnamelist
    nl.enkfio_in=YAML_Unquoted.new('.true.')
    nl.enkfio_out=YAML_Unquoted.new('.false.')
    nl.incr_fname=File.join(datadir,'gincr.b')
    nl=env.run.namelists.outputnamelist
    lib_mod_namelist_file(env,nlfile2,env.run.namelists)
    rundir
  end

  def lib_run_prep_rocoto(env)

    # Do the standard run prep first.

    rundir=lib_run_prep_std(env)

    # Copy FIMwfm to run directory.

    rocoto_srcdir=valid_dir(File.join(File.expand_path(".."),"FIMwfm"))
    rocoto_rundir=File.join(rundir,"..","FIMwfm")
    logd "Copying #{rocoto_srcdir} -> #{rocoto_rundir}"
    FileUtils.cp_r(rocoto_srcdir,rocoto_rundir)

    # Prepare for XML modifications.

    require 'rexml/document'
    rocoto_xmldir=File.join(rocoto_rundir,"rocoto_xml")
    env.run.rocoto_xml=valid_dir(rocoto_xmldir)

    # Modify prep XML file: Disable retries.

    prepxml=valid_file(File.join(rocoto_xmldir,"FIM_Prep.xml"))
    prepdoc=REXML::Document.new(File.new(prepxml),{:raw=>:all})
    prepdoc.delete_element("//task/dependency")
    REXML::XPath.first(prepdoc,"//task").add_attribute("maxtries","1")
    File.open(prepxml,"w") { |f| f.puts prepdoc }
    logd "Modified prep XML file #{prepxml}"

    # Modify fim XML file: Disable retries, and remove dependency on the
    # "spectral" task, which normally copied initialization data for the model,
    # but which we do not run here because the test suite's initializatino data
    # is already available.

    fimxml=valid_file(File.join(rocoto_xmldir,"FIM_FIM.xml"))
    fimdoc=REXML::Document.new(File.new(fimxml),{:raw=>:all})
    REXML::XPath.first(fimdoc,"//task").add_attribute("maxtries","1")
    File.open(fimxml,"w") { |f| f.puts fimdoc }
    logd "Modified fim XML file #{fimxml}"

    # Modify post XML file: Disable retries, and set the GRID_NAMES and
    # GRID_SPECS environment-variable settings to use the values defined in the
    # workflow XML file.

    postxml=valid_file(File.join(rocoto_xmldir,"FIM_PostFimoutTrue.xml"))
    postdoc=REXML::Document.new(File.new(postxml),{:raw=>:all})
    REXML::XPath.first(postdoc,"//task").add_attribute("maxtries","1")
    REXML::XPath.match(postdoc,"//task/envar").each do |e|
      e[1].text="&GRID_NAMES;" if e[0].text=="GRID_NAMES"
      e[1].text="&GRID_SPECS;" if e[0].text=="GRID_SPECS"
    end
    File.open(postxml,"w") { |f| f.puts postdoc }
    logd "Modified post XML file #{postxml}"

    # Modify workflow XML file: Set some environment-variable values as
    # appropriate to this run.

    wfxml=valid_file(File.join(rocoto_xmldir,env.run.wfxml))
    wfdoc=REXML::Document.new(File.new(wfxml))
    datadir=valid_dir(File.join(tmp_dir,"data"))
    entity_map={
      "DATADIR"=>datadir,
      "DATADR2"=>datadir,
      "FIM_HOME"=>File.expand_path(File.join(rundir,"..")),
      "FIM_RUN"=>rundir
    }
    wfdoc.doctype.each do |e|
      if e.is_a?(REXML::Entity) and entity_map.has_key?(e.name)
        e.replace_with(REXML::Entity.new(e.name,entity_map[e.name]))
      end
    end
    File.open(wfxml,"w") { |f| f.puts wfdoc }
    logd "Modified workflow XML file #{wfxml}"

    # Symlink "fim" and "post" directories so that output can be found in a
    # directory structure similar to that created for non-Rocoto runs.

    g=env.run.namelists.cntlnamelist.glvl
    k=env.run.namelists.cntlnamelist.nvl
    p=env.run.namelists.queuenamelist.computetasks
    t=env.run.namelists.timenamelist.yyyymmddhhmm
    subdir=File.expand_path(File.join(rundir,"fim_#{g}_#{k}_#{p}_#{t}"))
    FileUtils.mkdir_p(subdir)
    fimdir_target=File.join(subdir,"fim_C")
    fimdir_linkname=File.join(subdir,"fim")
    postdir_target=File.join(subdir,"post_C","228","NAT","grib1")
    postdir_linkname=File.join(subdir,"post")
    FileUtils.ln_s(fimdir_target,fimdir_linkname)
    logd "Linked #{fimdir_linkname} -> #{fimdir_target}"
    FileUtils.ln_s(postdir_target,postdir_linkname)
    logd "Linked #{postdir_linkname} -> #{postdir_target}"

    # Remember useful paths for later.

    env.run.subdir=subdir
    env.run.wfxml=wfxml

    # Return the path to the unique run directory.

    rocoto_rundir
  end

  def lib_run_prep_std(env)
    rundir=env.run.ddts_root
    runfiles=env.build.ddts_result
    logd "Copying #{runfiles} -> #{rundir}"
    FileUtils.cp_r(runfiles,rundir)
    rundir=File.join(rundir,File.basename(runfiles))
    datadir=valid_dir(File.join(tmp_dir,"data"))
    nl=env.run.namelists.queuenamelist
    nl.chem_datadir=File.join(datadir,"chem")
    nl.datadir=datadir
    nl.datadr2=datadir
    unless (nl=env.run.namelists.landnamelist)
      nl=(env.run.namelists.landnamelist=OpenStruct.new)
    end
    nl.landdatdir=datadir
    unless (nl=env.run.namelists.toponamelist)
      nl=(env.run.namelists.toponamelist=OpenStruct.new)
    end
    nl.topodatfile=File.join(datadir,"wrf5mintopo.dat")
# Added for ocean data path
    unless (nl=env.run.namelists.modelnamelist)
      nl=(env.run.namelists.modelnamelist=OpenStruct.new)
    end
    nl.atm_trnsecdir=File.join(datadir,"ocean/")
    unless (nl=env.run.namelists.oceannamelist)
      nl=(env.run.namelists.oceannamelist=OpenStruct.new)
    end
    nl.ocn_trnsecdir=File.join(datadir,"ocean/")
    nl.inicondir=File.join(datadir,"ocean/")

    nlfile=File.join(rundir,"FIMnamelist")
    lib_mod_namelist_file(env,nlfile,env.run.namelists)
    rundir
  end

  def lib_run_restart(env,prepkit)

    # Run once, check the run, then set 'restart' flag and run again.

    fim_stdout,rundir=lib_run_std(env,prepkit)
    return nil if fim_stdout.nil?
    unless job_check(fim_stdout,lib_re_str_success)
      die "Restart first half failed"
    end
    env.run.restart=true
    lib_run_std(env,prepkit)
  end

  def lib_run_rocoto(env,prepkit)

    # Do some setup.

    rocoto_rundir=prepkit
    ts=env.run.namelists.timenamelist.yyyymmddhhmm[0..-3]
    storefile=File.join(rocoto_rundir,"fimts.store")
    logfile=File.join(rocoto_rundir,"log","workflow","workflow_#{ts}.log")
    FileUtils.touch(logfile)
    logd "Rocoto log file is #{logfile}"
    logd "Rocoto store file is #{storefile}"
    logd "Rocoto workflow file is #{env.run.wfxml}"

    # Iterate Rocoto. For unknown reasons, ext() hangs when running this shell
    # command (IO::read appears to be the source of the hang), so run it with
    # IO.popen directly instead.

    cmd="module purge && module load rocoto && rocotorun -d #{storefile} -w #{env.run.wfxml}"
    logd "Iterating Rocto with command: #{cmd}"
    loop do
      IO.popen(cmd)
      str1="This cycle is complete" # All workflow tasks completed
      str2="in state DEAD"          # Some workflow task failed
      break if job_check(logfile,str1) or job_check(logfile,str2)
      sleep 60
    end
    [logfile,env.run.subdir]
  end

  def lib_run_std(env,prepkit)
    rundir=prepkit
    jobid=nil
    subdir=nil
    re1=Regexp.new(lib_re_str_job_id)
    re2=Regexp.new(lib_re_str_run_dir)
    ss='subfim'
    mach_name=invoke(:lib_submit_script,:run,env)
    restart=env.run.restart
    ss+=".restart" if restart
    cmd="#{File.join(rundir,ss)} #{mach_name} #{rundir}"
    logd "Submitting job with command: #{cmd}"
    output,status=Thread.exclusive { ext(cmd,{:msg=>"Job submission failed"}) }
    output.each do |e|
      e.chomp!
      logd e unless e=~/^\s*$/
      jobid=e.gsub(re1,'\1') if re1.match(e)
      subdir=e.gsub(re2,'\1') if re2.match(e)
    end
    if jobid.nil?
      logi "ERROR: Job ID not found in #{ss} output"
      return nil
    end
    job_activate(jobid,self)
    if subdir.nil? and not restart
      logi "ERROR: Run directory not found in #{ss} output"
      return nil
    end
    rundir.replace(File.join(rundir,subdir)) unless restart
    qs="Queued with job ID #{jobid}"
    qs+=" (restart)" if restart
    logi qs
    lib_wait_for_job(env,jobid)
    job_deactivate(jobid)
    nls=env.run.namelists
    g=nls.cntlnamelist.glvl
    k=nls.cntlnamelist.nvl
    p=nls.queuenamelist.computetasks.delete('"')
    stdout=File.join(rundir,"fim#{g}_#{k}_#{p}","fim","stdout")
    [stdout,rundir]
  end

  def lib_suite_post(env)
    env.qmon.stop
  end

  def lib_suite_prep(env)
    env.qmon=Qmon.new(60)
  end

  # CUSTOM METHODS (NOT CALLED BY DRIVER)

  def lib_link_data(target)
    linkname=File.join(tmp_dir,"data")
    FileUtils.rm_f(linkname)
    FileUtils.ln_s(target,linkname)
    lib_validate_data(valid_dir(linkname))
  end

  def lib_mod_namelist_file(env,nlfile,nlenv)
    h=convert_o2h(nlenv)
    sets=h.reduce([]) do |m0,(n,kv)|
      inner=kv.reduce([]) do |m1,(k,v)|
        v="\"#{quote_string(v)}\""
        logd "Set namelist #{n}:#{k}=#{v}"
        m1.push("-s #{n}:#{k}=#{v}")
      end
      m0.concat(inner)
    end
    nml=valid_file(File.join(env.build.ddts_result,"nml"))
    cmd="#{nml} -i #{nlfile} -o #{nlfile} #{sets.join(" ")}"
    Thread.exclusive { ext(cmd,{:msg=>"Failed to edit #{nlfile}"}) }
  end

  def lib_re_str_job_id
    'The job (\d+).* has been submitted.'
  end

  def lib_re_str_run_dir
    'Made directory (.*)'
  end

  def lib_re_str_success
    '(Program exited normally)|(PROGRAM nems *HAS ENDED)'
  end

  def lib_submit_script_jet(env)
    'jet'
  end

  def lib_submit_script_theia(env)
    'theia'
  end

  def lib_suffix(env)
    p_or_s=(env.build.args.include?("serial"))?("s"):("p")
    omp=(env.build.args.include?("omp"))?("_omp"):("")
    debug=(env.build.args.include?("debug"))?("_debug"):("")
    "_"+env.build.args.sub(/^\s+/,'').sub(/\s+.*$/,'')+"_#{p_or_s}#{omp}#{debug}"
  end

  def lib_validate_data(dir)
    logd "Validating data..."
    expected={
      '110010000.gfs.t00z.sanl'             => 'f0f0438321bae6bed3a1ac5fef6f5bb2',
      '110010000.gfs.t00z.sfcanl'           => '767938b50aed1e69ce63fb946f67b553',
      '142460000.gfs.t00z.sanl'             => '289a94ab0db93c6065d58b6b312cd546',
      '142460000.gfs.t00z.sfcanl'           => 'e6a39b2d6ad7c77fb58482b8ace9f7e8',
      'chem/anthro_binary'                  => '1667e3c46700f8016fabc7e89e66bb4d',
      'chem/chemltln.dat'                   => '1c9b7062dbab9600a02b600da57a9237',
      'chem/clay.dat'                       => 'cbe104bdddfefd2413aea634ae32337b',
      'chem/current.dat'                    => '5c0015327f4229236a3c742ab798c2a1',
      'chem/dm0_binary'                     => 'f08eaf8cdcc128356b211a015f7fc8c8',
      'chem/e_bc.dat'                       => '43bfe2b9055e306010d78b2b78d37924',
      'chem/e_oc.dat'                       => 'a7c38366b9067951f13be8ad9b6d1a51',
      'chem/e_pm_10.dat'                    => 'b394b2ec3cc1b813cfc53cf34bfa04c5',
      'chem/e_pm_25.dat'                    => '92c2b8062a324be4e002af6dfe4a07ff',
      'chem/e_so2.dat'                      => 'dcd5e1314566eeb2be462b9e3fabfe33',
      'chem/e_sulf.dat'                     => 'bd5a040b3f34a3a8d9ad87811161f5ec',
      'chem/erod1.dat'                      => '469a1cc0ce9b0221362e453619462f21',
      'chem/erod2.dat'                      => 'febc03d214dc6e9fa36b5b21eeaa005b',
      'chem/erod3.dat'                      => 'febc03d214dc6e9fa36b5b21eeaa005b',
      'chem/erod_binary'                    => '4c32c20c9b77084f6c55c102627e14ab',
      'chem/gocart_backgd_littlee'          => 'aa4b0410959c1c803ad67d1ea599371e',
      'chem/MODIS55N155W'                   => 'fcb9e86f1ba0f58af3eb4425e8123ed7',
      'chem/MODIS60N150W'                   => '4112914e2f7042d87e32f9fc3fe28b40',
      'chem/MODISHEADER'                    => '47470b2a003eb72723aa6861c777936b',
      'chem/OGE55N155W'                     => '7b0483a4c42a7a5c7a5f0d87c4feeb2d',
      'chem/OGE60N150W'                     => '00a1c0722687c1e808c3fe238fd06aa2',
      'chem/OGEHEADER'                      => '47470b2a003eb72723aa6861c777936b',
      'chem/prep_chem_sources'              => '6d45c40b52eefbd5d1ebbbc9f4edb133',
      'chem/prep_chem_sources_template.inp' => 'd327bf954b03723e45f4429784c0c588',
      'chem/sand.dat'                       => '6cb851c1cfe05ef135b2d53a7ccbc6b1',
      'chem/volcanic.dat'                   => '486a1b8d37795669d71c62602dc2eb9a',
      'climaeropac_global.txt'              => '0056ee94cc520ff9ba759e9660362042',
      'co2historicaldata_2007.txt'          => '4dbae00d7c84d4158659c7bb42d5ce52',
      'co2historicaldata_2008.txt'          => '2706b78df491d9311b93dbb35193e2d1',
      'co2historicaldata_2009.txt'          => 'ee945c7df4e3e9d61adac649a9df5aea',
      'co2historicaldata_2010.txt'          => '51c88e0bb297f00e7e0db1ae08e85c9f',
      'co2historicaldata_2011.txt'          => '29f21c1480e41a598eb999d48fdb0062',
      'co2historicaldata_2012.txt'          => '4f58a15d5d3096ef0af6aae28141ef99',
      'co2historicaldata_glob.txt'          => 'd0b1e1403f7726c5fc5745f9312c9ac3',
      'dm0_binary'                          => 'f08eaf8cdcc128356b211a015f7fc8c8',
      'geo_em.d01.nc'                       => '95af13ee9fa08e684d3b7b3187e83852',
      'gfsltln_t1534.dat'                   => '0f99cf12ec9adcbed091f652fe24c17a',
      'gfsltln_t382.dat'                    => 'e811632fc303d44c6ecb78827c1614cb',
      'gfsltln_t574.dat'                    => '348172e9f4d23a590482f0d501ff11c9',
      'gincr.b'                             => '7b70a0f9cf28bfc3bf2543aa751eeade',
      'global_mtnvar.t1534'                 => '49152f2339d0fa3c5582940c3c292422',
      'global_mtnvar.t382'                  => 'da9c6a8e15eebd3588e6b69e9b7baa86',
      'global_mtnvar.t574'                  => '7b34cdde6e35dc56c93e71f9a3dec716',
      'global_o3prdlos.f77'                 => '79fc6d61a66eb2512b5b95e84198fdb2',
      'grid_spec_05.nc.T382'                => '44528ca574770eaf89505aabad3b811e',
      'HADISST_MONTHLY.1991-present'        => 'cb0abf1ee3f99fd40b8cec1ccbfe702c',
      'HADISST_MONTHLY_CLIM.1981-2010'      => '43533f846ba98a3690c740d6cdab034d',
      'ocean/G4_Bering_thrufl.dat'          => 'b1ddf91334b8279921712d3a524383b1',
      'ocean/G4_Drake_thrufl.dat'           => '4d0e0508d7e460d4cad779ec87124439',
      'ocean/G4_Indo_thrufl.dat'            => '76589932780c3cf178102d33fdf8399d',
      'ocean/G4_05N_latcirc.dat'            => 'f43f35883c4bd80b634f8b3f916d7ad4',
      'ocean/G4_05S_latcirc.dat'            => 'd690219ad924b8ffb45ba42ddf208d62',
      'ocean/G4_15N_latcirc.dat'            => '71a606ef3adbe9a86cbe09c419b2e037',
      'ocean/G4_15S_latcirc.dat'            => '782ae081b76ac48a6d369fbcd8ee8dfe',
      'ocean/G5_05N_latcirc.dat'            => '93f1346f5d98327f4c330338ac291449',
      'ocean/G5_20N_latcirc.dat'            => '0a0cc0ec2fbb30eae0890caaa8361acf',
      'ocean/G5_35N_latcirc.dat'            => '71a961a73213294025d5e2a1a8f98861',
      'ocean/G5_50N_latcirc.dat'            => '60e9e90bc61c5ca87a38a5eb2599ee9a',
      'ocean/G5_65N_latcirc.dat'            => 'cee3c4f32616b51e4447a986e648bd70',
      'ocean/G4_25N_latcirc.dat'            => '3c6ae4366790d44b9b4856893176369d',
      'ocean/G4_25S_latcirc.dat'            => '1f8e699bf5aa5564c75d96ba1c367dd1',
      'ocean/G4_35N_latcirc.dat'            => '96944ec86f30af1f406222051dbf826f',
      'ocean/G4_35S_latcirc.dat'            => 'e5243d2dd94bb643404576f8dcc135d8',
      'ocean/G4_45N_latcirc.dat'            => '923156721e98a471421dcfe64787113a',
      'ocean/G4_45S_latcirc.dat'            => 'e24aeff016a8b307355ff0db98fa08f1',
      'ocean/G4_55N_latcirc.dat'            => '167aa0d25a5d2f745103c31dc244f421',
      'ocean/G4_55S_latcirc.dat'            => '67ebcbb925668be0bf666ae41576399b',
      'ocean/G4_65N_latcirc.dat'            => '2dcc1c136f13206c46bed7e72c943fa7',
      'ocean/G4_65S_latcirc.dat'            => 'e6fb45ef7e8793bc867c6f40b00e2fa1',
      'ocean/G5_05S_latcirc.dat'            => '2d058a3a39b36ba1417a3f16cc9b7be4',
      'ocean/G5_15N_latcirc.dat'            => '704672af3ac3d69ea63b8d396d6c1498',
      'ocean/G5_15S_latcirc.dat'            => 'f15a66db5951b3de080f36ce8814faea',
      'ocean/G5_25N_latcirc.dat'            => '8446b3f2485a8d2f6fc8b55021d7d385',
      'ocean/G5_25S_latcirc.dat'            => '937a56bbed399990bef9170a20b0fb0e',
      'ocean/G5_35S_latcirc.dat'            => 'b8bdf06288d6db27ec48033b2bb34ad6',
      'ocean/G5_45N_latcirc.dat'            => 'fd2184dd7f56dbe5e6d08853c734f524',
      'ocean/G5_45S_latcirc.dat'            => '82538a94bc73017e840cf76599a0ed46',
      'ocean/G5_55N_latcirc.dat'            => 'bc32cc01a86ab30d6350c4322d67e5ae',
      'ocean/G5_55S_latcirc.dat'            => '5cdce36c92d2e5cd025b393deafbf26d',
      'ocean/G5_65S_latcirc.dat'            => 'eae26398d0bf363ac74605ae4f5e0b0e',
      'ocean/G5_30N_Atl.dat'                => '9425618203a67f1563e91736a33501e0',
      'ocean/G5_35N_Atl.dat'                => '726446c42b387c27a4dadb7839d5c9bb',
      'ocean/G5_40N_Atl.dat'                => 'ed9b9b9ba984f4a19c37a4ca657778b6',
      'ocean/G5_45N_Atl.dat'                => 'c75b36c9868fe93e974c56ad4fbad02f',
      'ocean/G5_40N_Pac.dat'                => '5b00a385eed86dc37780ca9efd13855b',
      'ocean/G5_35N_Pac.dat'                => '63a8a7b4080a02ad94386eca3bc2bdbc',
      'ocean/G5_Indo_thrufl.dat'            => '36d660c801e9cd7aa6941cd0b58056d2',
      'ocean/G5_Bering_thrufl.dat'          => 'f8a46c67bcb5aaf2e6886dd180384097',
      'ocean/G5_Drake_thrufl.dat'           => 'c868ad664a503b9fd92616aa85dfb4f6',
      'ocean/transects_g5.nml'              => '8c8760c8f987be878d5bfa9a5dc4761a',
      'ocean/transects_atm_G4.nml'          => '436d076a5a82f0325778d56f45e58174',
      'ocean/transects_atm_G5.nml'          => '4daee1d92ac93403faca1d271a9d2ae1',
      'ocean/slantbounds_g4.nml'            => '1b47d9d5e4cef7363384eeb93d331986',
      'ocean/slantbounds_g5.nml'            => '43818ed3a14bc6078adbb189c53138bf',
      'ocean/transects_g4.nml'              => '4e1de49c5830cf0a65743400c922835d',
      'ocean/input_g4/topo1_g4.4bin'        => '1fa2d68fe19ff485bf0084e4afe0c279',
      'ocean/input_g4/lwdn_core1_g4.nc'     => '4e508deb96c07ce7518a1243ee9168b8',
      'ocean/input_g4/pout_g4_k26_sig1b.4bin' => '105d29d12cf8e4e97bc1fcd2c4f9ef8d',
      'ocean/input_g4/precip_core1_g4.nc'   => '1d68d8865f490fa1d607b53c518ce484',
      'ocean/input_g4/q_10_core1_g4.nc'     => '03c50cfcd7623df42bacdbd36d9d26ab',
      'ocean/input_g4/runoff_core1_g4.nc'   => 'e4046d9840adb6bb29933efcd014d46f',
      'ocean/input_g4/saln_g4_k26_sig1b.4bin' => '286b46ce881c7c3d88a8670a4cac99ba',
      'ocean/input_g4/sss_core1_g4.nc'      => 'd2de46bee883d4444f79c5725ca92d13',
      'ocean/input_g4/swdn_core1_g4.nc'     => '378a18d1a746197bc7cf8bd9db3e3e10',
      'ocean/input_g4/t_10_core1_g4.nc'     => '770e2bd49d3208b0d27f89ec65cfd9d3',
      'ocean/input_g4/temp_g4_k26_sig1b.4bin' => 'bd4c084762039c54576a1b47359ffeaa',
      'ocean/input_g4/u_10_core1_g4.nc'     => '603b086cae9363b103a72b117c8d90d2',
      'ocean/input_g4/v_10_core1_g4.nc'     => '184bbc8de58af3137d984abe8ca084ab',
      'ocean/input_g5/lwdn_core1_g5.nc'     => '0022ee08350f687eff375b4870ea29d7',
      'ocean/input_g5/pout_g5_k26_sig1b.4bin' => '720e49590d68425ae85b3d021df96dcb',
      'ocean/input_g5/precip_core1_g5.nc'   => '3f2c889a8dfca68380b9f28c45339875',
      'ocean/input_g5/q_10_core1_g5.nc'     => '8c56d64ff94bd01f04d50c39f581a335',
      'ocean/input_g5/runoff_core1_g5.nc'   => '1e5eba5ff2162eef9007681d564530ab',
      'ocean/input_g5/saln_g5_k26_sig1b.4bin' => 'cbc5686175645a65f833af702a76ee6c',
      'ocean/input_g5/sss_core1_g5.nc'      => '0d9cdaa387dcda15f88d5f375c8483e1',
      'ocean/input_g5/swdn_core1_g5.nc'     => '67f03752d46349c576fd450734977554',
      'ocean/input_g5/t_10_core1_g5.nc'     => '64f214221977a79fa72bc3e720fac5e9',
      'ocean/input_g5/temp_g5_k26_sig1b.4bin' => 'bae3539fb89d5ea3e08ed5aad87d5153',
      'ocean/input_g5/topo1_g5.4bin'        => 'd830ffb81003e2d349ee99d97334839e',
      'ocean/input_g5/u_10_core1_g5.nc'     => '83d98befbb10d7c55bdc5f23ece804f1',
      'ocean/input_g5/v_10_core1_g5.nc'     => '47e272725bd967da23389d4f148f2ce0',
      'ocean_bcs_ltln.360x180.dat'          => '079b9bd8513429c467c68570030dbdf2',
      'ocnanl.gdas.2011010100.ice_model.res.nc' => '84f7c37ec6c5303a802810cbb9dda787',
      'ocnanl.gdas.2011010100.ocean_temp_salt.res.nc' => '796ea2967e191fec59a3c504466aa284',
      'ocnanl.gdas.2014090300.ice_model.res.nc' => 'da94fc86c6badf9041bc94456221ced4',
      'ocnanl.gdas.2014090300.ocean_temp_salt.res.nc' => '72ad634c2b634fdee7e9672c649efa2f',
      'rucgrid'                             => '521d14c526e6500fdfeac1b413ca0945',
      'sfc_emissivity_idx.txt'              => '8f64bd9e63cc049909ec287decec1c44',
      'sfcanl.gdas.2011010100'              => '68218a6e48a92fad64a8e76c9e657162',
      'siganl.gdas.2011010100'              => '8d722e344e559c8532d715017d4877ac',
      'solarconstant_noaa_an.txt'           => '1f065539ad66484fd1be236132167760',
      'volcanic_aerosols_1980-1989.txt'     => 'aff7c065dd6b7b0ee4327fa93f018731',
      'volcanic_aerosols_1990-1999.txt'     => 'b7c758cfca095a0743a7dbeda8ee8f08',
      'wrf5mintopo.dat'                     => 'ab846b2154b57c003c0dceb8e0c6d5ee'
    }
    actual=Dir.glob("#{dir}/**/*")
    actual=actual.delete_if { |e| File.directory?(e) }
    actual=actual.reduce({}) do
      |m,e| m.merge!({e.sub(/#{dir}\/?/,'')=>Digest::MD5.file(e).to_s})
    end
    actual.keys.sort.each do |k|
      die "Unexpected data file: #{k}" unless expected[k]
      unless actual[k]==expected[k]
        logd "Checksum validation failed for #{dir}/#{k}"
        logd "  Expected #{expected[k]}"
        logd "    Actual #{actual[k]}"
        die "Error validating test-suite data, see #{logfile}"
      end
      logd "  #{k}: OK"
    end
    logd "Validating data: OK"
  end

  def lib_wait_for_job(env,jobid)
    live=%w[E H Q R T W S]
    interval=env.qmon.interval
    tolerance=interval*5
    job_update_time=nil
    job_clock=nil
    while true
      sleep interval*1.25
      result=env.qmon.query(jobid)
      if result[:query_status]==0
        if (state=result[:job_state])
          job_update_time=Time.now
          break unless live.include?(state)
        else
          if job_update_time
            if (age=Time.now-job_update_time) > tolerance
              logd "No news on job #{jobid} in #{age} seconds, assuming complete"
              break
            end
          else
            job_clock||=Time.now # clock is ticking...
            if (age=Time.now-job_clock) > tolerance
              die "Job #{jobid} unknown after #{age} seconds, aborting"
            end
          end
        end
      else
        logd "Batch-system query error:"
        result[:error].split("\n").each { |e| logd e }
        if (age=Time.now-result[:updated]) > tolerance
          die "No response from batch system in #{age} seconds, aborting"
        end
      end
    end
    result[:exit_status]
  end

  class Qmon

    attr_reader :interval

    def initialize(interval)
      require "rexml/document"
      require "thread"
      @error=nil
      @interval=interval
      @jobs={}
      @lock=Mutex.new
      @running=false
      @thread=nil
      start
    end

    def query(jobid)
      @lock.synchronize do
        {
          :error=>@error,
          :exit_status=>((x=@jobs[jobid])?(x[:exit_status]):(nil)),
          :job_state=>((x=@jobs[jobid])?(x[:job_state]):(nil)),
          :query_status=>@status,
          :updated=>@updated
        }
      end
    end

    def start
      return if @running
      update
      @updated=Time.now
      @running=true
      sleep @interval
      @thread=Thread.new do
        while @running
          t0=Time.now
          update
          elapsed=Time.now-t0
          if (adjusted_interval=@interval-elapsed) > 0
            sleep adjusted_interval
          end
        end
      end
    end

    def stop
      return unless @running
      @running=false
      @thread.join if @thread
      @thread=nil
    end

    private

    def update
      jobs={}
      output=IO.popen("qstat -x 2>&1") { |io| io.read }
      status=$?.exitstatus
      if status==0
        error=nil
        doc=REXML::Document.new(output)
        jobs=REXML::XPath.each(doc,"//Data/Job").reduce({}) do |m,e|
          jobid=REXML::XPath.first(e,"Job_Id").text.sub(/^([0-9]+).*/,'\1')
          job_state=REXML::XPath.first(e,"job_state").text
          exit_status=(x=REXML::XPath.first(e,"exit_status"))?(x.text.to_i):(nil)
          m[jobid]={:job_state=>job_state,:exit_status=>exit_status}
          m
        end
      else
        error=output
      end
      @lock.synchronize do
        @error=error
        @status=status
        unless error
          @jobs=jobs
          @updated=Time.now
        end
      end
    end

  end # class Qmon

end
