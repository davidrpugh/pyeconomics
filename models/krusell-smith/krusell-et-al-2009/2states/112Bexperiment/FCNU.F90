double precision function FCNU(X)

!note: ia has to be in Globals. 

	use Globals
	use Numerical_Libraries	
	use ValueModule


	
	
	implicit none
	real(8) X

	
	FCNU=-ValueF(X)
	

end function

