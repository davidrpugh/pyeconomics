!	----------------------------------------------------------------------
!	File name: NewtonU.f90
!	
!	----------------------------------------------------------------------


subroutine NewtonU(ap2, vf2)

	use Globals

	use Numerical_Libraries
	use ValueModule
	
	implicit none

	real(8), parameter:: tol=1.0e-10
    REAL(8) ap2, vf2
	integer IT
	integer, parameter:: ITMAX=20
	real(8) x0, x1, fp, fp2
    logical ifstop
      	
	x0 = AS(ia,ieps,ibta,ik,iz)


	ifstop = .true.
	IT=0
	do while (ifstop)
	IT=IT+1

	if (IT>ITMAX) then
!	print*, "warning at Newton"
	ap2=x0
	vf2=ValueF(x0)
	ifstop = .false.
	end if

	fp = ValueFp(x0)
	
	if (abs(fp)<tol) then
	ap2=x0
	vf2=ValueF(x0)
	ifstop = .false.
	else
	fp2 = min(ValueFp2(x0),-0.0000001)
	x1=x0-fp/fp2
	x0=max(amin,min(amax,x1))
	end if

	end do

end subroutine

