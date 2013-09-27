!	----------------------------------------------------------------------
!	File name: SolveValueFunction.f90
!	
!	Sovles individual's maximization problem by iterating value function, 
!	generates decision rule	for asset holdings.
!	----------------------------------------------------------------------


subroutine SolveValueFunction()

	use Globals

	implicit none

	real(8) errV, errA
	integer iter
	real(8) NewV(na,neps,nbta,nk,nz)
	real(8) NewAS(na,neps,nbta,nk,nz)

	if (iterLOM == 1) then
	MAXITER = 1000
	else
	MAXITER = 1000
	end if

	! start iteration

	do iter = 1, MAXITER 


		!	evaluate expected value function

		call ExpectedValueFunction()


		!	maximization procedure

		do iz = 1, nz
		do ik = 1, nk
	
		irate = zgrid(iz)*alpha*(kgrid(ik)/((1.0-mu(iz))*hour))**(alpha-1)-delta
		!print*, irate
		wage =  zgrid(iz)*(1-alpha)*(kgrid(ik)/((1.0-mu(iz))*hour))**alpha
	    !print*, "wage=", wage
			do ia = 1, na
			do ieps = 1, neps
			! **** Calculate income  and maximize ****
				income = (1.0+irate)*agrid(ia) + wage*epsgrid(ieps)*hour+ui*(1.0-epsgrid(ieps))
				do ibta = 1, nbta
				call BrentU(NewAS(ia,ieps,ibta,ik,iz), NewV(ia,ieps,ibta,ik,iz))
!			
				end do



            end do
			end do

		end do
		end do


		!	calculate errors

		errV = maxval(dabs(NewV - V))
		errA = maxval(abs(NewAS - AS))
 
		! update value function and decision rules

		V  = NewV
		AS  = NewAS



		! convergence check

		if (errV < tol_value .and. errA < tol_value) then
			print*, "Individual Optimization is done at iteration " ,iter
			print*, " "
			exit 
		else if (mod(iter,showerr) == 0) then
			print*, "---------SolveValueFunction iteration------ ", iter
			print*, "errV = ", errV, "   errA = ", errA

			print*, " "

		end if

	end do

	
	!	save value functions and decision rules if final simulation

!	if (final) then
		call SaveValueFunction()
!	end if

end subroutine
