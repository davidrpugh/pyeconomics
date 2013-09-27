!	----------------------------------------------------------------------
!	File name: BisectU.f90
!	
!	----------------------------------------------------------------------


subroutine BisectU2(ap2, vf2)

	use Globals

	use Numerical_Libraries
	use ValueModule
	
	implicit none

	real(8), parameter:: tol=1.0e-6
    REAL(8) ap2, vf2
	real(8) x0, x1, x3, ff0, ff1, ff3
    logical ifstop
      	
	x0= amin

	x3 = min(agrid(ia)+10.0,amax)

	x1 = (x0+x3)/2.0D0


	ff0 = ValueF2p(x0)
	ff1 = ValueF2p(x1)
	ff3 = ValueF2p(x3)

	if (ff0<0.0D0) then !corner
		ap2 = x0
		vf2=ValueF2(ap2)
	else if (ff3>0.0D0) then
		ap2 = x3
!
		vf2=ValueF2(ap2)
	else ! start the bisection loop.
		ifstop = .true.
		do while (ifstop)
		   if (ff1>0.0D0) then
		   x0=x1
		   x1 = (x0+x3)/2.0D0
		   ff1 = ValueF2p(x1)
		   else  
		   x3=x1
		   x1 = (x0+x3)/2.0D0
		   ff1 = ValueF2p(x1)
		   end if
		   if (abs(x3-x0)<TOL) then
		   ap2=x1
		   vf2=ValueF2(ap2)
		   ifstop = .false.
		   end if

		end do
	end if


end subroutine

