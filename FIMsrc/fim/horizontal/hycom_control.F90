module hycom_control

integer            :: numtr=1		! # of passive tracers
integer, parameter :: subcycle_ocn=1	! ocn_timestep / atm_timestep
logical, parameter :: modesplit=.false.	! barotropic/clinic mode split [y/n]
logical, parameter :: buoyfor=.true.	! if true, add thermal to wind forcing
logical, parameter :: flxcor=.true.	! if true, nudge surface fluxes to maintain global salt and heat
integer, parameter :: trcfrq=30		! trcr transport every 'trcfrq' steps
integer, parameter :: bclin_frq=30	! baclin=batrop*bclin_frq
integer, parameter :: iocnmx=2		! 1 = full column; 2 = mxlyr only KPP
integer, parameter :: kgap=5		! k spacing for lyr-by-lyr diagnostics
logical, parameter :: latdiw=.true.	! K-PROF:  activate lat.dep. int.wave mixing
logical            :: diag_sst=.false.	! perform SST diagnostics
integer, parameter :: nsecs=10		! max.# of massflx diagno.transects
character*40       :: thruflsecs(nsecs)	! transects for thruflow diagnostics
character*40       :: bdcurrsecs(nsecs)	! transects for bdry current diagno.
character*40       :: bylayrsecs(nsecs)	! transects for massflux by layer
end module hycom_control
