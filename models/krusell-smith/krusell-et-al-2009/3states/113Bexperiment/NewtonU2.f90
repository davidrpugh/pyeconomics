!	----------------------------------------------------------------------
!	File name: NewtonU2.f90
!	
!	----------------------------------------------------------------------


subroutine NewtonU2(ap2, vf2)

	use Globals

	use Numerical_Libraries
	use ValueModule
	
	implicit none

	real(8), parameter:: tol=1.0e-12
    REAL(8) ap2, vf2
	integer IT
	integer, parameter:: ITMAX=20
	real(8) x0, x1, fp, fp2
    logical ifstop
   
   
	if (ibta==1) then   	
	x0 = AS2l(time,ia,iPge,iPbe,iPgs,iPbs,iPbf)
	else if (ibta==2) then
	x0 = AS2m(time,ia,iPge,iPbe,iPgs,iPbs,iPbf)
	else
	x0 = AS2h(time,ia,iPge,iPbe,iPgs,iPbs,iPbf)
	end if

	ifstop = .true.
	IT=0
	do while (ifstop)
	IT=IT+1

	if (IT>ITMAX) then
!	print*, "warning at Newton"
	ap2=x0
	vf2=ValueF2(x0)
	ifstop = .false.
	end if

	fp = ValueF2p(x0)
	
	if (abs(fp)<tol) then
	ap2=x0
	vf2=ValueF2(x0)
	ifstop = .false.
	else
	fp2 = min(ValueF2p2(x0),-0.000000000001)
	x1=x0-fp/fp2
	x0=max(amin,min(amax,x1))
	end if

	end do

end subroutine

