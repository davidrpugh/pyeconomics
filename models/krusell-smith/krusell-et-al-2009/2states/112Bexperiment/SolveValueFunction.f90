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
	real(8) NewV(na,nP,nP,nbta,nk)
	real(8) NewAS(na,nP,nP,nbta,nk)

	if (iterLOM == 1) then
	MAXITER = 2000
	else
	MAXITER = 1000
	end if

	! start iteration

	do iter = 1, MAXITER 


		!	evaluate expected value function

		call ExpectedValueFunction()


		!	maximization procedure

		
		do ik = 1, nk

		irate(ik) = alpha*(kgrid(ik)/((1.0-mubar)*hour))**(alpha-1)-delta
		wage(ik) =  (1-alpha)*(kgrid(ik)/((1.0-mubar)*hour))**alpha
	   
			do iPg = 1, nP
			do iPb = 1, nP
			do ia = 1, na
		
			! **** Calculate income  and maximize ****
				income = (1.0+irate(ik))*agrid(ia) + wage(ik)*(0.5*PGgrid(iPg)+0.5*PBgrid(iPb))*hour+ui*(1.0-(0.5*PGgrid(iPg))-(0.5*PBgrid(iPb)))
				do ibta = 1, nbta
				 
				  call BrentU(NewAS(ia,iPg,iPb,ibta,ik), NewV(ia,iPg,iPb,ibta,ik))
!
				end do

            end do
			end do
            end do
		
        end do
	

		!	calculate errors

		errV = maxval(dabs(NewV - V))
		errA = maxval(dabs(NewAS - AS))
 
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

	if (final) then
		call SaveValueFunction()
	end if

end subroutine
