!	----------------------------------------------------------------------
!	File name: GenerateSeries.f90
!
!	
!	----------------------------------------------------------------------

subroutine GenerateSeries()

  use Globals
  implicit none

 
  integer t
  real(8), parameter:: tiny = 1.0D-15

  ! Assume that the initial state is good
  ! THIS HAS TO BE CHANGED IF z(1)=B or D
  ! p(i,t) is the probability of state i (1 bad, 2 good) at period t
  
  pi(1,1) = 1.0D0-tiny
  
  pi(2,1) = 1.0D0-pi(1,1)

  do t = 2, N1+1
    
	pi(2,t) =  pi(1,t-1)*trZ(1,2)+pi(2,t-1)*trZ(2,2)
    pi(1,t) = 1.0D0-pi(2,t)

  end do
 
  do t = 1, N1
    
	zseries(t) =  pi(2,t)*zgrid(2)+pi(1,t)*zgrid(1)
    useries(t) =  pi(2,t)*mu(2)+pi(1,t)*mu(1)

  end do
   open(1, file='Output\pi2.txt', status='unknown')
   open(2, file='Output\pi1.txt', status='unknown')
   open(3, file='Output\zseries.txt', status='unknown')
   open(4, file='Output\useries.txt', status='unknown')
   write(1, '(F15.6)') pi(2,:)
   write(2, '(F15.6)') pi(1,:)
   write(3, '(F15.6)') zseries
   write(4, '(F15.6)') useries
  end subroutine
