	SUBROUTINE tridiag(a,b,c,r,uu)

	use Globals
	IMPLICIT NONE
	integer, parameter :: n=na
	double precision a(n),b(n),c(n),r(n),uu(n),gam(n),bet
	
	integer j



	bet=b(1)

	uu(1)=r(1)/bet
	do j=2,n
		gam(j)=c(j-1)/bet
		bet=b(j)-a(j-1)*gam(j)
		uu(j)=(r(j)-a(j-1)*uu(j-1))/bet
	end do
	do j=n-1,1,-1
		uu(j)=uu(j)-gam(j+1)*uu(j+1)
	end do
	END SUBROUTINE tridiag
