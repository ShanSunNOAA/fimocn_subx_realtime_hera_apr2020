 &PlotIcos
  ginfofile ="/scratch1/portfolios/BMC/nim/Ning.Wang/FIM_20130822_1833/FIMrun/zeussubfim_9900/fim8_64_128/prep/icos_grid_info_level.dat"
  datafile ="/scratch1/portfolios/BMC/nim/Ning.Wang/FIM_20130822_1833/FIMrun/zeussubfim_17918/fim8_64_128/fim/fim_out_th3D000066hr"
  grid_level = 8
  var_name = "th3D"    ! name of the variable
  nvlvls = 64                ! number of vertical levels of the dataset
  level = 25                 ! the level of the data to plot
  proj = 'CE'
  latc = 27.344             ! US
  lonc = -80.0              ! US
  extent = 30.0 35.0        ! extent of the domain
!  extent = 180.0 180.0      ! extent of the domain
  map_vis = 1               ! plot map
  cell_vis = 0              ! plot voronoi cell
  ll_vis = 1                ! plot lat lon lines
  ipn_label = 0             ! plot ipn index label
  print_version = .false./        ! create graphic file for printing
/

