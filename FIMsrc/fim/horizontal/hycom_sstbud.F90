   subroutine sstbud(nproc,text,tem,thk,leap)
!
   use module_control  ,only: nip
   use module_constants,only: perm
   use fimnamelist    ,only: kdm,itest
   use hycom_constants ,only: wet,onem,nsstpro,sstproc
   use hycom_variables ,only: tdpold,sstndcy
   use stencilprint
!
! --- record time-averaged SST tendencies caused by various physical processes
!
! --- nproc - process counter (normally 1,...,nsstpro, except:)
!             nproc = 0 sets tdpold = tem*dp (initialization)
! --- text  - character string identifying process associated with 'nproc'
! --- tem   - temperature profile after completion of process 'nproc'
! --- thk   - layer thickness profile after completion of process 'nproc'
! --- leap  - leapfrog time slot
!
   implicit none
   integer  ,intent(IN) :: nproc,leap
   character,intent(IN) :: text*16
!SMS$DISTRIBUTE (dh,2) BEGIN
   real    ,intent(IN)  :: thk (kdm,nip,2)		! layer thickness
   real    ,intent(IN)  :: tem (kdm,nip)		! temperature
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE (dh,1) BEGIN
   real prs(nip),tdpnew(nip),work(nip)
!SMS$DISTRIBUTE END
   integer              :: i,k
   real                 :: q,delp,r_delp
   character            :: flnm*120
   logical              :: vrbos
   real,parameter       :: dpth = 30.	!  control volume depth (m)
!
   delp=dpth*onem			!  control volume depth in Pa
   r_delp=1./delp
!
   if (nproc.eq.0) then

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE(vrbos)
    do i=1,nip
     tdpold(i)=0.
     prs(i)=0.
     if (wet(i) > 0) then
      do k=1,kdm
       tdpold(i)=tdpold(i)+tem(k,i)*(min(prs(i)+thk(k,i,leap),delp)	&
                                    -min(prs(i)              ,delp))
       prs(i)=prs(i)+thk(k,i,leap)
       if (prs(i).gt.delp) exit
      end do
     end if
    end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END
!
   else				!  nproc > 0
!
    if (nproc.gt.nsstpro) stop 'sstbud error: nproc > nsstpro'
    sstproc(nproc)=text
!
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE(vrbos)
    do i=1,nip
     tdpnew(i)=0.
     prs(i)=0.
     if (wet(i) > 0) then
      vrbos=i.eq.itest
      do k=1,kdm
       tdpnew(i)=tdpnew(i)+tem(k,i)*(min(prs(i)+thk(k,i,leap),delp)	&
                                    -min(prs(i)              ,delp))
       prs(i)=prs(i)+thk(k,i,leap)
       if (prs(i).gt.delp) exit
      end do
      if (vrbos) print '(i8,a,i3,a,2f9.3)',perm(i),			&
        ' (sstbud) old/new T after process',nproc,sstproc(nproc),	&
        tdpold(i)*r_delp,tdpnew(i)*r_delp

! --- build up time integral of t*dpth change attributed to process 'nproc'
      sstndcy(nproc,i)=sstndcy(nproc,i)+(tdpnew(i)-tdpold(i))*r_delp
!
! --- 'new' t*dpth becomes 'old' t*dpth for next process
      tdpold(i)=tdpnew(i)
     end if
    end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END
!
    work(:)=sstndcy(nproc,:)
    call stencl(work,1,100.,'SST change (.01 C) due to '//sstproc(nproc))
   end if			! nproc = or > 0
!
   return
   end
