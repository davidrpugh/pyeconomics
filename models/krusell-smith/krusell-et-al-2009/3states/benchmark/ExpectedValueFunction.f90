!	----------------------------------------------------------------------
!	File name: ExpectedValueFunction.f90
!
!	Calculates expected value function given current state variables.
!	----------------------------------------------------------------------



subroutine ExpectedValueFunction()

	use Globals
	use PolyInterpModule

	implicit none

	real(8) kp, VP(na,neps,nbta,nz)
	integer iepsp, izp, ibtap


	EV = 0
	do ik = 1, nk
	do iz = 1, nz


		!	next period aggregate capital

		if (iz==3) then
		kp = dexp(KcoefG(1) + KcoefG(2)*dlog(kgrid(ik)))
		else if (iz==2) then
		kp = dexp(KcoefD(1) + KcoefD(2)*dlog(kgrid(ik)))
		else
		kp = dexp(KcoefB(1) + KcoefB(2)*dlog(kgrid(ik)))
		end if

		do izp = 1, nz 


			!	next period value function: VP(ia,iepsp,izp) with kp(ik,iz)

            do iepsp = 1, neps
			do ia = 1, na
            do ibtap = 1, nbta
				VP(ia,iepsp,ibtap,izp) = polyinterp1(kp, kgrid, V(ia,iepsp,ibtap,:,izp))
            end do
			end do
			end do
           
		end do
			
		!	calculate expected value function

		do ieps = 1, neps
		do ia = 1, na
        do ibta = 1, nbta
			do iepsp = 1, neps
			do izp = 1, nz
            do ibtap = 1, nbta
				EV(ia,ieps,ibta,ik,iz) = EV(ia,ieps,ibta,ik,iz)		&
								 + VP(ia,iepsp,ibtap,izp)*trZ(iz,izp)*trEps(iz,izp,ieps,iepsp)*trBta(ibta,ibtap)
            end do
			end do
			end do
        end do
		end do
        end do

	end do
	end do
		
	call setsplineV()

end subroutine
