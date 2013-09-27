!--------------------------------------
!File name splint.f90
!--------------------------------------

subroutine splint(xa,ya,y2a,xx,yy,yyp,yydp,n)


	use Globals
	use Numerical_Libraries	


	
	
	implicit none
	
	integer n
	double precision xa(n), ya(n), y2a(n), xx, yy, yyp, yydp

	integer j, jhi, jlo
	double precision a, b, h, asq, bsq
	

	jlo=1
	jhi=n
	do 
		if (jhi-jlo<=1) then
		exit
		else
		j=(jhi+jlo)/2
		if (xa(j)>xx) then
		jhi=j
		else
		jlo=j
		end if
		end if

	end do
		
	h=xa(jhi)-xa(jlo)
	a=(xa(jhi)-xx)/h
	b=(xx-xa(jlo))/h

	asq=a*a
	bsq=b*b

	yy=a*ya(jlo)+b*ya(jhi)+((asq*a-a)*y2a(jlo)+(bsq*b-b)*y2a(jhi))*(h*h)/6.0D0
	yyp=(ya(jhi)-ya(jlo))/h-(3.0D0*asq-1.0D0)*h*y2a(jlo)/6.0D0+(3.0D0*bsq-1.0D0)*h*y2a(jhi)/6.0D0
	yydp=a*y2a(jlo)+b*y2a(jhi)




end subroutine

