!	----------------------------------------------------------------------
!	File name: Initialize.f90
!	----------------------------------------------------------------------



subroutine InitializeParameters()


	use Globals

	implicit none
	
!    **** The unemployment rate for bad and good times. First row is bad times, second row is good times ****

	mu(1) = .10
	mu(2) = .04



end subroutine



subroutine InitializeCoefficients()

	use Globals

	implicit none

	KcoefG(1)=0.10207 
	KcoefG(2)=0.95963
	KcoefB(1)=0.09397
	KcoefB(2)=0.96096

	
end subroutine

subroutine InitializeValue()

	use Globals
	use ValueModule

	implicit none
    real ll

   	!initialize value function and decision rules

	do iz = 1, nz
	do ik = 1, nk

		irate = alpha*zgrid(iz)*(3.0*kgrid(ik)/(1.0-mu(iz)))**(alpha-1)-delta
		wage =  (1-alpha)*zgrid(iz)*(3.0*kgrid(ik)/(1.0-mu(iz)))**alpha

		do ia = 1, na
		do ieps = 1, neps
			income = (1.0+irate)*agrid(ia) + wage*epsgrid(ieps)*hour+ui*(1.0-epsgrid(ieps))
		do ibta = 1, nbta
			AS(ia,ieps,ibta,ik,iz) = min(agrid(ia)+ui,amax)
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
