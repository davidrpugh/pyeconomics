subroutine ComputePprimes2

	use Globals
	use Numerical_Libraries

	implicit none

	integer izp, i
    real(8), parameter:: tiny = 1.0D-15


	!calculate the values of Pprime based on gbars and the probability

	Pprimege=0.0D0

		do iinterv=1, ninterv
			do iz = 1, nz
			do ieps = 1, neps
				
				Pprimege(iinterv) = Pprimege(iinterv) + (1.0D0/max(pi(2,time+1),tiny))*gbar(iinterv,iz,2,ieps,4)*trZ(iz,2)*Pvalue(iz,ieps)

		
			end do
			end do
		end do
	
	Pprimebe=0.0D0

		do iinterv=1, ninterv
			do iz = 1, nz
			do ieps = 1, neps
				
				Pprimebe(iinterv) = Pprimebe(iinterv) + (1.0D0/max(pi(1,time+1),tiny))*gbar(iinterv,iz,1,ieps,4)*trZ(iz,1)*Pvalue(iz,ieps)

			
			end do
			end do
		end do

	Pprimegs=0.0D0

		do iinterv=1, ninterv
			do iz = 1, nz
			do ieps = 1, neps
				
				Pprimegs(iinterv) = Pprimegs(iinterv) + (1.0D0/max(pi(2,time+1),tiny))*(1.0D0/max(1.0D0-Pprimege(iinterv),tiny))*gbar(iinterv,iz,2,ieps,3)*trZ(iz,2)*Pvalue(iz,ieps)

			end do
			end do
		end do

	Pprimebs=0.0D0

		do iinterv=1, ninterv
			do iz = 1, nz
			do ieps = 1, neps
				
				Pprimebs(iinterv) = Pprimebs(iinterv) + (1.0D0/max(pi(1,time+1),tiny))*(1.0D0/max(1.0D0-Pprimebe(iinterv),tiny))*gbar(iinterv,iz,1,ieps,3)*trZ(iz,1)*Pvalue(iz,ieps)

	
			end do
			end do
		end do

	Pprimebf=0.0D0

		do iinterv=1, ninterv
			do iz = 1, nz
			do ieps = 1, neps
				
				Pprimebf(iinterv) = Pprimebf(iinterv) + (1.0D0/max(pi(1,time+1),tiny))*(1.0D0/max(1.0D0-Pprimebe(iinterv),tiny))*(1.0D0/max(1.0D0-Pprimebs(iinterv),tiny))*gbar(iinterv,iz,1,ieps,2)*trZ(iz,1)*Pvalue(iz,ieps)


			end do
			end do
		end do

    do i=1, ninterv

      Pprimege(i)=max(0.0D0, min(1.0D0,Pprimege(i)))
      Pprimebe(i)=max(0.0D0, min(1.0D0,Pprimebe(i)))
      Pprimegs(i)=max(0.0D0, min(1.0D0,Pprimegs(i)))
      Pprimebs(i)=max(0.0D0, min(1.0D0,Pprimebs(i)))
      Pprimebf(i)=max(0.0D0, min(1.0D0,Pprimebf(i)))

    end do


end subroutine