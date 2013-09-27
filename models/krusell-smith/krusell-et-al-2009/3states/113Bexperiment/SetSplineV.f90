!--------------------------------------
!File name SetSplineV.f90
!--------------------------------------

subroutine setsplineV()


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

	do ik = 1, nk
	do iPge = 1, nP	
	do iPbe = 1, nP	
	do iPgs = 1, nP	
	do iPbs = 1, nP	
	do iPbf = 1, nP	
	do ibta = 1, nbta
	
	do i=1,n
		xx(i)=agrid(i)
		if (ibta==1) then
		yy(i)=EVl(i,iPge,iPbe,iPgs,iPbs,iPbf,ik)
		else if (ibta==2) then
		yy(i)=EVm(i,iPge,iPbe,iPgs,iPbs,iPbf,ik)
		else
		yy(i)=EVh(i,iPge,iPbe,iPgs,iPbs,iPbf,ik)
		end if
	end do

	do i=2,n-1
		a(i)=(xx(i)-xx(i-1))/6.0D0
		b(i)=(xx(i+1)-xx(i-1))/3.0D0
		c(i)=(xx(i+1)-xx(i))/6.0D0
		r(i)=(yy(i+1)-yy(i))/(xx(i+1)-xx(i))-(yy(i)-yy(i-1))/(xx(i)-xx(i-1))
	end do

	call tridiag(a,b,c,r,uu)

	do i=1,n
		if (ibta==1) then
		EV2pl(i,iPge,iPbe,iPgs,iPbs,iPbf,ik)=uu(i)
		else if (ibta==2) then
		EV2pm(i,iPge,iPbe,iPgs,iPbs,iPbf,ik)=uu(i)
		else 		
		EV2ph(i,iPge,iPbe,iPgs,iPbs,iPbf,ik)=uu(i)
		end if

	end do

	end do
	end do
	end do
	end do
	end do
	end do
	end do

end subroutine

