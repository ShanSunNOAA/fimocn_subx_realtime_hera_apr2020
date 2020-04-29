module module_outFMTed
contains
subroutine outFMTed(its,pr3d,ph3d,us3d,vs3d,dp3d,mp3d,vor,tr,rh3d,tk3d,ws3d,st3d,sm3d,rn2d,pw2d,pq2d,ts2d,us2d,hf2d,qf2d,rsds,rlds,time)
use module_control,only: nvlp1,nip,ntra,dt
use fimnamelist   ,only: nvl

implicit none
! External variable declarations:
integer,intent(IN)    :: its
!sms$distribute (dh,1) begin
real   ,intent(IN)    :: rn2d(nip)
real   ,intent(INOUT) :: pw2d(nip)
real   ,intent(INOUT) :: pq2d(nip)
real   ,intent(IN)    :: ts2d(nip),us2d(nip),hf2d(nip),qf2d(nip),rsds(nip),rlds(nip)
!sms$distribute end
!sms$distribute (dh,2) begin
real   ,intent(IN)    :: us3d(nvl  ,nip),vs3d(nvl,nip),dp3d(nvl,nip)
real   ,intent(IN)    :: pr3d(nvlp1,nip)
real   ,intent(IN)    :: mp3d(nvl  ,nip)
real   ,intent(IN)    :: vor (nvl  ,nip)
real   ,intent(IN)    :: ws3d(nvl  ,nip)
real   ,intent(IN)    :: ph3d(nvlp1,nip)
real   ,intent(IN)    :: tr(nvl    ,nip,ntra)
real   ,intent(INOUT) :: rh3d(nvl  ,nip)
real   ,intent(INOUT) :: tk3d(nvl  ,nip)
real   ,intent(IN)    :: st3d(4,nip),sm3d(4,nip)
!sms$distribute end
integer,intent(in)    :: time
integer               :: ivl,ipn,lunout=40
character(80)         :: filename

write(filename,"('fim_out_',i6.6)") time
open (lunout,file=filename,form="formatted")
write(lunout,*) 'its =',its
write(lunout,100) 'pr3d',((ivl,ipn,pr3d(ivl,ipn),ipn=1,nip),ivl=1,nvlp1)
write(lunout,100) 'ph3d',((ivl,ipn,ph3d(ivl,ipn),ipn=1,nip),ivl=1,nvlp1)
write(lunout,100) 'us3d',((ivl,ipn,us3d(ivl,ipn),ipn=1,nip),ivl=1,nvl)
write(lunout,100) 'vs3d',((ivl,ipn,vs3d(ivl,ipn),ipn=1,nip),ivl=1,nvl)
write(lunout,100) 'dp3d',((ivl,ipn,dp3d(ivl,ipn),ipn=1,nip),ivl=1,nvl)
write(lunout,100) 'mp3d',((ivl,ipn,mp3d(ivl,ipn),ipn=1,nip),ivl=1,nvl)
write(lunout,100) 'vor ',((ivl,ipn,vor (ivl,ipn),ipn=1,nip),ivl=1,nvl)
write(lunout,100) 'th3d',((ivl,ipn,tr(ivl,ipn,1),ipn=1,nip),ivl=1,nvl)
write(lunout,100) 'qv3d',((ivl,ipn,tr(ivl,ipn,2),ipn=1,nip),ivl=1,nvl)
write(lunout,100) 'qw3d',((ivl,ipn,tr(ivl,ipn,3),ipn=1,nip),ivl=1,nvl)
write(lunout,100) 'rh3d',((ivl,ipn,rh3d(ivl,ipn),ipn=1,nip),ivl=1,nvl)
write(lunout,100) 'tk3d',((ivl,ipn,rh3d(ivl,ipn),ipn=1,nip),ivl=1,nvl)
write(lunout,100) 'ws3d',((ivl,ipn,ws3d(ivl,ipn),ipn=1,nip),ivl=1,nvl)
write(lunout,100) 'st3d',((ivl,ipn,st3d(ivl,ipn),ipn=1,nip),ivl=1,4  )
write(lunout,100) 'sm3d',((ivl,ipn,sm3d(ivl,ipn),ipn=1,nip),ivl=1,4  )
write(lunout,101) 'rn2d',(     ipn,rn2d(    ipn),ipn=1,nip)
write(lunout,101) 'pw2d',(     ipn,pw2d(    ipn),ipn=1,nip)
write(lunout,101) 'pq2d',(     ipn,pw2d(    ipn),ipn=1,nip)
write(lunout,101) 'ts2d',(     ipn,ts2d(    ipn),ipn=1,nip)
write(lunout,101) 'us2d',(     ipn,us2d(    ipn),ipn=1,nip)
write(lunout,101) 'hf2d',(     ipn,hf2d(    ipn),ipn=1,nip)
write(lunout,101) 'qf2d',(     ipn,qf2d(    ipn),ipn=1,nip)
write(lunout,101) 'rsds',(     ipn,rsds(    ipn),ipn=1,nip)
write(lunout,101) 'rlds',(     ipn,rlds(    ipn),ipn=1,nip)
close(lunout)
100 format(' Variable ',a5/(2i10,1pe20.7))
101 format(' Variable ',a5/( i20,1pe20.7))
end subroutine outFMTed
end module module_outFMTed
