!*********************************************************************
!	Stop program for icosahedral flow-following global model
!	Alexander E. MacDonald  12/27/2004
!*********************************************************************

subroutine finalize

  use module_core_setup            ,only: iam_fim_task, my_comm
  use module_fim_chem_finalize     ,only: chem_finalize
  use module_fim_cpl_finalize      ,only: cpl_finalize
  use module_fim_dyn_finalize      ,only: dyn_finalize
  use module_fim_phy_finalize      ,only: phy_finalize
  use module_fim_wrf_phy_finalize  ,only: wrf_phy_finalize
  use module_initial_chem_namelists,only: chem_opt

  implicit none

#include <gptl.inc>

  integer :: mype = 0               ! mpi rank (use 0 for serial)
  integer :: ret                    ! return code from gptl functions
  real*8  :: tmain                  ! time taken by 'main'
  integer :: procsiz, rss, share, text, datastack ! output from gptlget_memusage
  real*8  :: rssmin, rssmax         ! min, max process sizes from GPTL
  character(len=64) :: name         ! region name

!sms$comm_rank(mype)

  if (mype >= 0) then
    ret = gptlpr (mype) ! print indented timing info to timing.<rank>
  end if

  if (iam_fim_task) then
    ! Print timing summary info across all tasks and (if applicable) threads
!SMS$insert ret = gptlpr_summary (my_comm)
    ret = gptlget_memusage (procsiz, rss, share, text, datastack)
    rssmin = rss
    rssmax = rss
!SMS$REDUCE(rssmin,min)
!SMS$REDUCE(rssmax,max)

    write(6,*)'If on zeus, dont trust process size numbers!'
    write(6,'(a,1p,2e9.2)') 'min,max process sizes=', rssmin, rssmax
    
    ret = gptlget_wallclock ('main', 0, tmain)  ! The "0" is thread number
    print*,'Total time =', tmain
    print*,' '

    call cpl_finalize ()
    call dyn_finalize ()
    call phy_finalize ()
    if (chem_opt > 0) then
      call wrf_phy_finalize ()
      call chem_finalize ()
    end if
    call datetime ()
  end if

  return
end subroutine finalize
