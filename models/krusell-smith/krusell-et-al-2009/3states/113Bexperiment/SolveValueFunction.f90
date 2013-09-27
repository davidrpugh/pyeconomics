!	----------------------------------------------------------------------
!	File name: SolveValueFunction.f90
!	
!	Sovles individual's maximization problem by iterating value function, 
!	generates decision rule	for asset holdings.
!	----------------------------------------------------------------------


subroutine SolveValueFunction()

	use Globals

	implicit none

	real(8) errVl, errVm, errVh, errAlo, errAmi, errAhi
	integer iter
	real(8) NewVl(na,nP,nP,nP,nP,nP,nk), NewVm(na,nP,nP,nP,nP,nP,nk), NewVh(na,nP,nP,nP,nP,nP,nk)
	real(8) NewASlo(na,nP,nP,nP,nP,nP,nk), NewASmi(na,nP,nP,nP,nP,nP,nk), NewAShi(na,nP,nP,nP,nP,nP,nk)

	if (iterLOM <4) then
	MAXITER = 300
	else
	MAXITER = 200
	end if

	! start iteration

	do iter = 1, MAXITER 


		!	evaluate expected value function

		call ExpectedValueFunction()


		!	maximization procedure

		
		do ik = 1, nk



		irate(ik) = alpha*(kgrid(ik)/((1.0-mubar)*hour))**(alpha-1)-delta
		wage(ik) =  (1-alpha)*(kgrid(ik)/((1.0-mubar)*hour))**alpha
	   
			do iPge = 1, nP
			do iPbe = 1, nP
			do iPgs = 1, nP
			do iPbs = 1, nP
			do iPbf = 1, nP
			do ia = 1, na
		
			! **** Calculate income  and maximize ****
				income = (1.0+irate(ik))*agrid(ia) + wage(ik)*(0.5*Pgrid(iPge)+0.5*Pgrid(iPbe))*hour &
				         +uih*(0.5*(1.0D0-Pgrid(iPge))*Pgrid(iPgs)+0.5*(1.0D0-Pgrid(iPbe))*Pgrid(iPbs)+0.5*(1.0D0-Pgrid(iPbe))*(1.0D0-Pgrid(iPbs))*Pgrid(iPbf)) &
						 +uil*(1.0-0.5*Pgrid(iPge)-0.5*Pgrid(iPbe)-0.5*(1.0D0-Pgrid(iPge))*Pgrid(iPgs)-0.5*(1.0D0-Pgrid(iPbe))*Pgrid(iPbs)-0.5*(1.0D0-Pgrid(iPbe))*(1.0D0-Pgrid(iPbs))*Pgrid(iPbf))

				do ibta = 1, nbta
				
					if (ibta==1) then
						if (iterLOM==1) then
						if (agrid(ia)>cutvalue .and. iter>50) then
						call NewtonU(NewASlo(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik), NewVl(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik))
						else
						call BisectU(NewASlo(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik), NewVl(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik))
						end if
						else
						if (agrid(ia)>cutvalue) then
						call NewtonU(NewASlo(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik), NewVl(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik))
						else
						call BisectU(NewASlo(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik), NewVl(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik))
						end if
						end if
					else if (ibta==2) then
						if (iterLOM==1) then
						if (agrid(ia)>cutvalue .and. iter>50) then
						call NewtonU(NewASmi(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik), NewVm(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik))
						else
						call BisectU(NewASmi(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik), NewVm(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik))
						end if
						else
						if (agrid(ia)>cutvalue) then
						call NewtonU(NewASmi(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik), NewVm(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik))
						else
						call BisectU(NewASmi(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik), NewVm(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik))
						end if
						end if
					else 
						if (iterLOM==1) then
						if (agrid(ia)>cutvalue .and. iter>50) then
						call NewtonU(NewAShi(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik), NewVh(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik))
						else
						call BisectU(NewAShi(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik), NewVh(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik))
						end if
						else
						if (agrid(ia)>cutvalue) then
						call NewtonU(NewAShi(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik), NewVh(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik))
						else
						call BisectU(NewAShi(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik), NewVh(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik))
						end if
						end if
											
					end if

				end do

            end do
			end do
            end do
            end do
			end do
            end do	
			
        end do
	

		!	calculate errors

		errVl = maxval(dabs(NewVl - Vl))
		errVm = maxval(dabs(NewVm - Vm))
		errVh = maxval(dabs(NewVh - Vh))
		errAlo = maxval(dabs(NewASlo - ASlo))
 		errAmi = maxval(dabs(NewASmi - ASmi))
		errAhi = maxval(dabs(NewAShi - AShi))
		! update value function and decision rules

		Vl  = NewVl
		Vm  = NewVm
		Vh  = NewVh
		ASlo  = NewASlo
		ASmi  = NewASmi		
		AShi  = NewAShi


!	if (iter==50) then
!		call SaveValueFunction()
		
!	end if


		! convergence check

		if (errVl+errVm+errVh < tol_value .and. errAlo+errAmi+errAhi < tol_value) then
			print*, "Individual Optimization is done at iteration " ,iter
			print*, " "
			exit 
		else if (mod(iter,showerr) == 0) then
			print*, "---------SolveValueFunction iteration------ ", iter
			print*, "iterLOM = ", iterLOM
			print*, "errVl = ", errVl, "   errAlo = ", errAlo
			print*, "errVm = ", errVm, "   errAmi = ", errAmi
			print*, "errVh = ", errVh, "   errAhi = ", errAhi
			print*, " "

		end if

	end do

	
	!	save value functions and decision rules if final simulation

	if (final) then
		call SaveValueFunction()
	end if

end subroutine
