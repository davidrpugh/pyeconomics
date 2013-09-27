subroutine ComputePprimes2

	use Globals
	use Numerical_Libraries

	implicit none

	integer izp
  


	!calculate the values of Pprime based on gbars and the probability
	Pprime=0.0D0
	do izp=1, nz
		do iinterv=1, ninterv
			do iz = 1, nz
			do ieps = 1, neps
				if (ieps==2) then
				Pprime(iinterv, izp) = min(max(Pprime(iinterv, izp) + gbar(iinterv,iz,izp,ieps)*trZ(iz,izp)*real(pi(iz,time)/pi(izp,time+1))*Pvalue(iz),0.0),1.0)
				else
				Pprime(iinterv, izp) = min(max(Pprime(iinterv, izp) + gbar(iinterv,iz,izp,ieps)*trZ(iz,izp)*real(pi(iz,time)/pi(izp,time+1))*(1-Pvalue(iz)),0.0),1.0)
				end if
			end do
			end do
		end do
	end do

  Pprimeb = Pprime(:, 1)
  Pprimeg = Pprime(:, 2)

end subroutine