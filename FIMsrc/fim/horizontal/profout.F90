module module_profout

contains

!*********************************************************************
! profout
! Profile output program for fim global model
! Jacques Middlecoff      April 2007
!*********************************************************************

  subroutine profout(     &
    its,PrintIpnDiag,     & ! index time step
    us3d,vs3d,dp3d,       & ! west wind, south wind, delta pres
    pr3d,th3d,mp3d,tk3d,  & ! pressure, theta, mont pot, temp (k)
    qv3d,rh3d,            & ! specific and relative humidity
    ph3d,                 &
    ts2d,us2d,hf2d,qf2d   ) ! skin temperature, ustar, sensible heat flux, water vapor flux

    use module_control  ,only: nvlp1,nip,ArchvStep
    use fimnamelist     ,only: nvl
    use module_constants,only: rd,cp

    implicit none

    integer,intent(IN) :: its,PrintIpnDiag
!sms$distribute (dh,1) begin
    real   ,intent(IN) :: ts2d(nip),us2d(nip),hf2d(nip),qf2d(nip)
!sms$distribute end
!sms$distribute (dh,2) begin
    real   ,intent(IN) :: us3d(nvl,nip),vs3d(nvl,nip),dp3d(nvl,nip)
    real   ,intent(IN) :: pr3d(nvlp1,nip),th3d(nvl,nip)
    real   ,intent(IN) :: mp3d(nvl,nip),tk3d(nvl,nip)
    real   ,intent(IN) :: ph3d(nvlp1,nip)
    real   ,intent(IN) :: qv3d(nvl,nip),rh3d(nvl,nip)
!sms$distribute end
    integer :: k

101 format (/' '/i7,1x,a/' '/'(ustar, T-skin, shflx & lhflx:',f9.2,2x,f11.6,2x,f11.6,2x,f11.6,')'/ &
      ' '/(2x,i3,2x,f9.3,2x,f9.2,2x,f12.5,2x,f11.3))

    if (PrintIpnDiag>0.and.(its.eq.0.or.mod(its,ArchvStep).eq.0)) then
!sms$serial (<us2d,ts2d,hf2d,qf2d,th3d,pr3d,qv3d,pr3d,in>:default=ignore) begin
      write (6,101)                                                                        &
        its-1,                                                                             &
        '    theta        t        specific humidity    press',                            &
        us2d(PrintIpnDiag),                                                                &
        ts2d(PrintIpnDiag),                                                                &
        hf2d(PrintIpnDiag),                                                                &
        qf2d(PrintIpnDiag),                                                                &
        (k,th3d(k,PrintIpnDiag),                                                           &
        th3d(k,PrintIpnDiag)*(5.e-6*(pr3d(k,PrintIpnDiag)+pr3d(k,PrintIpnDiag)))**(rd/cp), &
        qv3d(k,PrintIpnDiag),                                                              &
        pr3d(k,PrintIpnDiag),k=1,nvl)
!sms$serial end
    endif

  end subroutine profout

end module module_profout
