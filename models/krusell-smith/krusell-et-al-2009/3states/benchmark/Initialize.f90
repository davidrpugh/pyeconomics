!	----------------------------------------------------------------------
!	File name: Initialize.f90
!	----------------------------------------------------------------------



subroutine InitializeParameters()


	use Globals

	implicit none
	
!    **** The unemployment rate for bad and good times. First row is bad times, second row is good times ****

	mu(1) = 0.10
	mu(2) = 0.10
	mu(3) = 0.04



end subroutine



subroutine InitializeCoefficients()

	use Globals

	implicit none

	KcoefG(1)=0.104264 
	KcoefG(2)=0.958536
	KcoefB(1)=0.117081
	KcoefB(2)=0.951348
	KcoefD(1)=0.090291
	KcoefD(2)=0.962583


	
end subroutine

subroutine InitializeValue()

	use Globals
	use ValueModule

	implicit none
    real ll

   	!initialize value function and decision rules

	do iz = 1, nz
	do ik = 1, nk

		irate = alpha*zgrid(iz)*(kgrid(ik)/((1.0-mu(iz))*hour))**(alpha-1)-delta
		wage =  (1-alpha)*zgrid(iz)*(kgrid(ik)/((1.0-mu(iz))*hour))**alpha

		do ia = 1, na
		do ieps = 1, neps
			income = (1.0+irate)*agrid(ia) + wage*epsgrid(ieps)*hour+ui(ieps)*(1.0-epsgrid(ieps))
		do ibta = 1, nbta
			AS(ia,ieps,ibta,ik,iz) = min(agrid(ia)+ui(ieps),amax)
			V(ia,ieps,ibta,ik,iz) = Utility(AS(ia,ieps,ibta,ik,iz))/(1.0-btagrid(ibta))
			EV2p(ia,ieps,ibta,ik,iz) = 0.0
			EV(ia,ieps,ibta,ik,iz) = 0.0
		end do				
        end do  
		end do
       
	end do
	end do


    open(1, file='Output\Vinitial.txt', status='unknown')
	open(2, file='Output\ASinitial.txt', status='unknown')
	

    do ia = 1, na
	do ieps = 1, neps	
    do ibta = 1, nbta
	do ik = 1, nk
	do iz = 1, nz
	

		write(1, '(F15.6)') V(ia,ieps,ibta,ik,iz)
		write(2, '(F15.6)') AS(ia,ieps,ibta,ik,iz)
		
    end do
	end do
	end do
	end do
	end do 
	print*, "initialize is done"




end subroutine
