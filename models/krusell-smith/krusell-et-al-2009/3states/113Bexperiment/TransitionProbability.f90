!	----------------------------------------------------------------------
!	File name: TransitionProbability.f90
!
!	calculates transition probabilities for idiosyncratic shocks and 
!	aggregate productivity shock.
!	----------------------------------------------------------------------



subroutine IdiosyncraticTransition()

	use Globals

	implicit none
	
	integer  i, j
	integer  izp
	! transition probability matrix for idiosyncratic shocks

	!trEps(z,z',eps,eps')
	trEps = 0.0
	
	! (z,z') = (b,b)
	trEps(1,1,1,1) = 0.9875D0
	trEps(1,1,1,2) = 0.0
    trEps(1,1,1,3) = 0.0
	trEps(1,1,1,4) = 1.0D0-0.9875D0

	trEps(1,1,2,1) = 0.9875D0
	trEps(1,1,2,2) = 0.0
    trEps(1,1,2,3) = 0.0
	trEps(1,1,2,4) = 1.0D0-0.9875D0

	trEps(1,1,3,1) = 0.033796296D0
	trEps(1,1,3,2) = 0.0
    trEps(1,1,3,3) = 0.0
	trEps(1,1,3,4) = 1.0D0-0.033796296D0

	trEps(1,1,4,1) = 0.0
	trEps(1,1,4,2) = 0.0
    trEps(1,1,4,3) = 0.03D0
	trEps(1,1,4,4) = 0.97D0


	! (z,z') = (g,b)
    trEps(2,1,1,1) = 0.93939393939D0
	trEps(2,1,1,2) = 0.0
	trEps(2,1,1,3) = 0.0
	trEps(2,1,1,4) = 1.0D0-0.93939393939D0

    trEps(2,1,2,1) = 0.93939393939D0
	trEps(2,1,2,2) = 0.0
	trEps(2,1,2,3) = 0.0
	trEps(2,1,2,4) = 1.0D0-0.93939393939D0

	trEps(2,1,3,1) = 0.75D0
	trEps(2,1,3,2) = 0.0
	trEps(2,1,3,3) = 0.0
	trEps(2,1,3,4) = 0.25D0

    trEps(2,1,4,1) = 0.0
	trEps(2,1,4,2) = 0.0421875D0
	trEps(2,1,4,3) = 0.028125D0
	trEps(2,1,4,4) = 1.0D0-0.0421875D0-0.028125D0


    ! (z,z') = (b,g)
    trEps(1,2,1,1) = 0.169726D0
	trEps(1,2,1,2) = 0.0
	trEps(1,2,1,3) = 0.0
	trEps(1,2,1,4) = 1.0D0-0.169726D0

    trEps(1,2,2,1) = 0.169726D0
	trEps(1,2,2,2) = 0.0
	trEps(1,2,2,3) = 0.0
	trEps(1,2,2,4) = 1.0D0-0.169726D0
 
    trEps(1,2,3,1) = 0.03D0
	trEps(1,2,3,2) = 0.0
	trEps(1,2,3,3) = 0.0
	trEps(1,2,3,4) = 0.97D0

    trEps(1,2,4,1) = 0.0
	trEps(1,2,4,2) = 0.0
	trEps(1,2,4,3) = 0.02977777778D0
	trEps(1,2,4,4) = 1.0D0-0.02977777778D0


	! (z,z') = (g,g)
	trEps(2,2,1,1) = 0.50D0
	trEps(2,2,1,2) = 0.0
	trEps(2,2,1,3) = 0.0
	trEps(2,2,1,4) = 0.50D0

	trEps(2,2,2,1) = 0.50D0
	trEps(2,2,2,2) = 0.0
	trEps(2,2,2,3) = 0.0
	trEps(2,2,2,4) = 0.50D0

	trEps(2,2,3,1) = 0.2462686567D0
	trEps(2,2,3,2) = 0.0
	trEps(2,2,3,3) = 0.0
	trEps(2,2,3,4) = 1.0D0-0.2462686567D0

	trEps(2,2,4,1) = 0.0
	trEps(2,2,4,2) = 0.0
	trEps(2,2,4,3) = 0.02791666666667D0
	trEps(2,2,4,4) = 1.0D0-0.02791666666667D0









		 	   

	!	cumulative transition matrix of idiosyncratic shocks for random drawing

	CtrEps(:,:,:,1) = trEps(:,:,:,1)
	do ieps = 2, neps
		CtrEps(:,:,:,ieps) = CtrEps(:,:,:,ieps-1) + trEps(:,:,:,ieps)
	end do
    
	

	!	Cutoff1(nz,nz',neps,1) : first cutoff

	do iz = 1, nz
	do izp = 1, nz
	do ieps = 1, neps

	Cutoff(iz,izp,ieps,1)=trEps(iz,izp,ieps,4)

	end do
	end do
	end do

	!	Cutoff(nz,nz',neps,2) : second cutoff

	do iz = 1, nz
	do izp = 1, nz
	do ieps = 1, neps

	Cutoff(iz,izp,ieps,2)=trEps(iz,izp,ieps,4)+trEps(iz,izp,ieps,3)

	end do
	end do
	end do

	!	Cutoff(nz,nz',neps,3) : third cutoff

	do iz = 1, nz
	do izp = 1, nz
	do ieps = 1, neps

	Cutoff(iz,izp,ieps,3)=trEps(iz,izp,ieps,4)+trEps(iz,izp,ieps,3)+trEps(iz,izp,ieps,2)

	end do
	end do
	end do
    
end subroutine


subroutine AggregateTransition()

	use Globals
	!use TauchenMod

	implicit none
	
	integer i, j


	! transition probability matrix for aggregate shock: 2 states case

	trZ(1,1) = 0.875D0
	trZ(1,2) = 1.0D0 - 0.875D0
	trZ(2,1) = 1.0D0 - 0.875D0
	trZ(2,2) = 0.875D0

	open(1, file='Output\trZ.txt', status='unknown')
	write(1, '(<nz>f12.6)') ((trZ(i,j),j=1,nz),i=1,nz)


	!	cumulative transition matrix of aggregate shock for random drawing

	CtrZ(:,1) = trZ(:,1)
	do iz = 2, nz
		CtrZ(:,iz) = CtrZ(:,iz-1) + trZ(:,iz)
	end do

	open(2, file='Output\CtrZ.txt', status='unknown')
	write(2, '(<nz>f12.6)') ((CtrZ(i,j),j=1,nz),i=1,nz)


end subroutine

subroutine BetaTransition()

	use Globals
	

	implicit none
	
	integer i, j


	! transition probability matrix for beta: 3-state stochastic beta


	trBta(1,1) = 199.0D0/200.0D0
	trBta(1,2) = 1.0D0/200D0
	trBta(1,3) = 0.0D0
	trBta(2,1) = 1.0D0/1600.0D0
    trBta(2,2) = 1598.0D0/1600.0D0
	trBta(2,3) = 1.0D0/1600.0D0
	trBta(3,1) = 0.0D0
	trBta(3,2) = 1.0D0/200.0D0
    trBta(3,3) = 199.0D0/200.0D0

	open(1, file='Output\trBta.txt', status='unknown')
	write(1, '(<nbta>f12.6)') ((trBta(i,j),j=1,nbta),i=1,nbta)


	!	cumulative transition matrix of aggregate shock for random drawing

	CtrBta(:,1) = trBta(:,1)
	do ibta = 2, nbta
		CtrBta(:,ibta) = CtrBta(:,ibta-1) + trBta(:,ibta)
	end do

	open(2, file='Output\CtrBta.txt', status='unknown')
	write(2, '(<nbta>f12.6)') ((CtrBta(i,j),j=1,nbta),i=1,nbta)


end subroutine
