double precision function FCNU2(X)

!note: ia has to be in Globals. 

	use Globals
	use Numerical_Libraries	
	use ValueModule


	
	
	implicit none
	real(8) X

	
	FCNU2=-ValueF2(X)
	

end function

