module ssfc2icos_mtnvar

contains

  subroutine read_mtnvar(mdrag3d,imax,jmax,mvar,filename)

    integer,intent(in)::imax,jmax,mvar
    real,intent(out)::mdrag3d(imax,jmax,mvar)
    character(len=*),intent(in)::filename
    integer::i, j, v_idx

    open(21,file=trim(filename),form="unformatted",status='old',iostat=i)
    if (i.ne.0) then
      write (*,'(a,a,a)') 'ERROR in read_mtnvar: Could not open ',trim(filename),'.'
      stop
    endif
    read(21,iostat=i) mdrag3d
    if (i.ne.0) then
      write (*,'(a)') 'ERROR in read_mtnvar: Could not read mdrag3d.'
      stop
    endif
    close(21)

!    do v_idx = 7, 10
!      do j = 1, jmax
!        do i = 1, imax
!          mdrag3d(i,j,v_idx) = max(0.0, mdrag3d(i,j,v_idx))
!        end do
!      end do
!    end do

  end subroutine read_mtnvar

end module ssfc2icos_mtnvar
