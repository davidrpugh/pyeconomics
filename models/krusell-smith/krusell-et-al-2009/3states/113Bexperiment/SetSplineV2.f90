!--------------------------------------
!File name SetSplineV.f90
!--------------------------------------

subroutine setsplineV2()


	use Globals
	use Numerical_Libraries	


	
	
	implicit none
	integer, parameter :: n=na
	integer i
	double precision a(n), b(n), c(n), r(n), xx(n), yy(n), uu(n)


	a(1)=0.0D0
	a(n)=1.0D0
	b(1)=1.0D0
	b(n)=-1.05D0
	c(1)=-1.05D0
    !c(1)=-1.2D0
	c(n)=0.0D0
	r(1)=0.0D0
	r(n)=0.0D0

	
	do iPge = 1, nP	
	do iPbe = 1, nP	
	do iPgs = 1, nP	
	do iPbs = 1, nP	
	do iPbf = 1, nP	
	do ibta = 1, nbta
	
	  	  


	do i=1,n
		xx(i)=agrid(i)
		yy(i)=EV2(i,iPge,iPbe,iPgs,iPbs,iPbf,ibta)
	end do

	do i=2,n-1
		a(i)=(xx(i)-xx(i-1))/6.0D0
		b(i)=(xx(i+1)-xx(i-1))/3.0D0
		c(i)=(xx(i+1)-xx(i))/6.0D0
		r(i)=(yy(i+1)-yy(i))/(xx(i+1)-xx(i))-(yy(i)-yy(i-1))/(xx(i)-xx(i-1))
	end do

	call tridiag(a,b,c,r,uu)

	do i=1,n
		EV2pp(i,iPge,iPbe,iPgs,iPbs,iPbf,ibta)=uu(i)
	end do

	end do
	end do
	end do
	end do
	end do
	end do

end subroutine

