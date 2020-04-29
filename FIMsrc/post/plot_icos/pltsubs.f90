      subroutine dfclrs
      dimension rgbv(3,15)
      dimension rgbv1(3,100)
! 15 basic color defined as NCAR GRAPHICS GUIDE TO NEW UTILITYS 
      data rgbv /0.00,0.00,0.00,       &        !black
                 0.70,0.70,0.70,       &
                 0.75,0.50,1.00,       &
                 0.50,0.00,1.00,       &
                 0.00,0.00,1.00,       &
                 0.00,0.50,1.00,       &
                 0.00,1.00,1.00,       &
                 0.00,1.00,0.60,       &
                 0.00,1.00,0.00,       &
                 0.70,1.00,0.00,       &
                 1.00,1.00,0.00,       &
                 1.00,0.75,0.00,       &
                 1.00,0.38,0.38,       &
                 1.00,0.00,0.38,       &
                 1.00,0.00,0.00/  
!
      call gscr(1,0,0.,0.,0.)  ! 0 black
      call gscr(1,1,1.,1.,1.)  ! 1 white
!     do  i=1,15
      do  i=2,15
      call gscr(1,i,rgbv(1,i),rgbv(2,i),rgbv(3,i))
      enddo
!
! cseri-1  color(16-20)
!
      do k=16,20
      rgbv1(1,k)=0.5 
      rgbv1(2,k)=1.-1./4.*(k-16) 
      rgbv1(3,k)=1. 
      enddo
!
! cseri-2  color(21-25)
      do k=21,25
      rgbv1(1,k)=0.5 
      rgbv1(2,k)=1.-1./4.*(k-21) 
      rgbv1(3,k)=0.63 
      enddo
!
! cseri-3  color(26-30)
      do k=26,30
      rgbv1(1,k)=0.5 
      rgbv1(2,k)=1.-1./4.*(k-26) 
      rgbv1(3,k)=0.0
      enddo
!
! cseri-4  color(31-35)
      do k=31,35
      rgbv1(1,k)=0. 
      rgbv1(2,k)=1.-1./4.*(k-31) 
      rgbv1(3,k)=1. 
      enddo
!
! cseri-5  color(36-40)
      do k=36,40
      rgbv1(1,k)=0. 
      rgbv1(2,k)=1.-1./4.*(k-36) 
      rgbv1(3,k)=0.63 
      enddo
!
! cseri-6  color(41-45) 
      do k=41,45
      rgbv1(1,k)=0.0 
      rgbv1(2,k)=1.-1./4.*(k-41) 
      rgbv1(3,k)=0.0
      enddo
!
! cseri-7  color(46-50)
!
      do k=46,50
      rgbv1(1,k)=0.25 
      rgbv1(2,k)=1.-1./4.*(k-46) 
      rgbv1(3,k)=0.0
      enddo
!
! cseri-8  color(51-55)
!
      do k=51,55
      rgbv1(1,k)=0.25 
      rgbv1(2,k)=1.-1./4.*(k-51) 
      rgbv1(3,k)=0.63
      enddo
!
! cseri-9  color(56-60)
!
      do k=56,60
      rgbv1(1,k)=0.25 
      rgbv1(2,k)=1.-1./4.*(k-56) 
      rgbv1(3,k)=1.0
      enddo
!
! cseri-10  color(61-65) 
      do k=61,65
      rgbv1(1,k)=0.75 
      rgbv1(2,k)=0.88  
      rgbv1(3,k)=1.-1./4.*(k-61) 
      enddo
!
! cseri-11  color(66-70) 
      do k=66,70
      rgbv1(1,k)=0.75 
      rgbv1(2,k)=0.63 
      rgbv1(3,k)=1.-1./4.*(k-66) 
      enddo
!
! cseri-12  color(71-75)
!
      do k=71,75
      rgbv1(1,k)=0.75 
      rgbv1(2,k)=0.
      rgbv1(3,k)=1.-1./4.*(k-71) 
      enddo
!
! cseri-13  color(76-80)
!
      do k=76,80
      rgbv1(1,k)=1
      rgbv1(2,k)=1.-1./4.*(k-76) 
      rgbv1(3,k)=0.
      enddo
!
! cseri-14  color(81-85)
!
      do k=86,90
      rgbv1(1,k)=1. 
      rgbv1(2,k)=1.-1./4.*(k-86) 
      rgbv1(3,k)=0.5
      enddo
!
! cseri-15  color(86-90)
!
      do k=81,85
      rgbv1(1,k)=0.9 
      rgbv1(2,k)=1.-1./4.*(k-81) 
      rgbv1(3,k)=0.9
      enddo
!
      do 100 i=16,90
      call gscr(1,i,rgbv1(1,i),rgbv1(2,i),rgbv1(3,i))
  100 continue
!
      call gscr(1,91,0.66,0.66,0.66)     ! index=91 gray
      call gscr(1,92,0.4,0.4,0.4)        ! index=92 dark gray
      call gscr(1,93,0.25,0.25,0.25)     ! index=93 dark gray
      call gscr(1,94,0.25,0.12,0.12)     ! index=94 dark gray
      call gscr(1,95,0.86,0.58,0.44)     ! index=95 tan
      call gscr(1,96,0.65,0.16,0.16)     ! index=96 brown
      call gscr(1,97,1.,0.,0.)           ! index=97 red 
      call gscr(1,98,0.1,0.3,0.1)        ! index=98 dark gray
      call gscr(1,99,0.14,0.56,0.14)     ! index=99 forest green
      call gscr(1,100,0.2,0.56,0.8)      ! index=100 sky blue 
      call gscr(1,100,0.1,0.28,0.4)      ! index=100 sky blue 
      return
      end
