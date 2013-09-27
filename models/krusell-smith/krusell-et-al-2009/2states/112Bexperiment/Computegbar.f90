subroutine Computegbar

	use Globals
	use Numerical_Libraries

	implicit none

	integer izp
  

!   sort the cutoff vector
!   We have eight cutoff values. We sort them here since it is not clear how they rank. 


!   First record where they come from
 
	icut=1
	do iz = 1, nz
	do izp = 1, nz
	do ieps = 1, neps
	
	indcut(icut)=icut !just indexing
	cutvecttemp(icut)=Cutoff(iz,izp,ieps)
	zindtemp(icut)=iz
	zpindtemp(icut)=izp
	epsindtemp(icut)=ieps
	icut=icut+1

	end do
	end do
	end do

!   sort by using DSVRGP. The sorting is by an increasing order.
!   cutvect is the ordered vector and it is the permutation of cutvecttemp by indcut

	call DSVRGP (ncut, cutvecttemp, cutvect, indcut)

!   record the values of z,zp and eps that they came from

	do icut=1,ncut

	  zind(icut)=zindtemp(indcut(icut))
	  zpind(icut)=zpindtemp(indcut(icut))
	  epsind(icut)=epsindtemp(indcut(icut))

	end do


!   for each interval, calculate the probability and the value of the P's.

!   probability and evaluating points


	probint(1)=cutvect(1)
	evalpoint(1)=cutvect(1)-smallval
	
	do iinterv=2, ninterv-1
	 probint(iinterv)=cutvect(iinterv)-cutvect(iinterv-1)
	 evalpoint(iinterv)=cutvect(iinterv)-smallval
	end do
	
	probint(ninterv)=1.0D0-cutvect(ninterv-1)
	evalpoint(ninterv)=1.0D0-smallval

	!construct gbar at each interval, by comparing the evaluating point with each cutoff.
	do iinterv=1, ninterv
	do iz = 1, nz
	do izp = 1, nz
	do ieps = 1, neps
	
	if (evalpoint(iinterv)<Cutoff(iz,izp,ieps)) then 
	gbar(iinterv,iz,izp,ieps)=1.0D0
	else
	gbar(iinterv,iz,izp,ieps)=0.0D0
	end if
	
	end do
	end do
	end do
	end do

	end subroutine