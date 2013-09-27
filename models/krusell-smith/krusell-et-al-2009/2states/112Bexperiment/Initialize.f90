!	----------------------------------------------------------------------
!	File name: Initialize.f90
!	----------------------------------------------------------------------



subroutine InitializeParameters()


	use Globals

	implicit none
	
!   The unemployment rate for bad (1) and good (2) states 

	mu(1) = 0.10
	mu(2) = 0.04



end subroutine



subroutine InitializeK()

	use Globals

	implicit none

	Kcoef(1) = 	0.02456 !(0.10205+0.09423)/2.0
	Kcoef(2) =  0.99003 !(0.95963+0.96084)/2.0
   
	K(1) = 11.20

	do ik = 2, Nperiod+1
   
     K(ik) = dexp(Kcoef(1)+Kcoef(2)*dlog(K(ik-1)))

    end do

    open(1, file='Output\K.txt', status='unknown')
	write(1, '(F15.6)') K



end subroutine


subroutine InitializeValue()

	use Globals
	use ValueModule

	implicit none
    real ll


   	!initialize value function and decision rules

    mubar = (mu(1)+mu(2))/2.0

	do ik = 1, nk

		irate(ik) = alpha*(3.0*kgrid(ik)/(1.0-mubar))**(alpha-1)-delta
		wage(ik) =  (1-alpha)*(3.0*kgrid(ik)/(1.0-mubar))**alpha

		do ia = 1, na
		do iPg = 1, nP
        do iPb = 1, nP
			income = (1.0+irate(ik))*agrid(ia) + wage(ik)*(PGgrid(iPg)+PBgrid(iPb))+ui*(1.0-PGgrid(iPg)-PBgrid(iPb))
		do ibta = 1, nbta
			AS(ia,iPg,iPb,ibta,ik) = min(agrid(ia)+ui,amax)
			V(ia,iPg,iPb,ibta,ik) = Utility(AS(ia,iPg,iPb,ibta,ik))/(1.0-btagrid(ibta))
		end do				
        end do  
		end do
       
	end do
	end do


    open(1, file='Output\Vinitial.txt', status='unknown')
	open(2, file='Output\ASinitial.txt', status='unknown')
	

    do ia = 1, na
	do iPg = 1, nP
	do iPb = 1, nP		
    do ibta = 1, nbta
	do ik = 1, nk

	

		write(1, '(F15.6)') V(ia,iPg,iPb,ibta,ik)
		write(2, '(F15.6)') AS(ia,iPg,iPb,ibta,ik) 
		
    end do
	end do
	end do
	end do
	end do 




end subroutine


subroutine InitializeValue2()

  use Globals
  use ValueModule
  
  implicit none
 

  do time = 1, N1
     
	 irate2(time) = zseries(time)*alpha*(real(K(time)/((1-useries(time))*hour))**(alpha-1))-delta
	 wage2(time) = zseries(time)*(1-alpha)*(real(K(time)/((1.0-useries(time))*hour))**(alpha))
   
	 do ia = 1, na
	 do iPg = 1, nP
	 do iPb = 1, nP
 
		do ibta = 1, nbta


		AS2(time,ia,iPg,iPb,ibta) = (AS(ia,iPg,iPb,ibta,1)+AS(ia,iPg,iPb,ibta,2)+AS(ia,iPg,iPb,ibta,3)+AS(ia,iPg,iPb,ibta,4))/4.0   
		V2(time,ia,iPg,iPb,ibta) = (V(ia,iPg,iPb,ibta,1)+V(ia,iPg,iPb,ibta,2)+V(ia,iPg,iPb,ibta,3)+V(ia,iPg,iPb,ibta,4))/4.0

     end do
	 end do
	 end do
	 end do
end do
end subroutine

subroutine InitializeValueRead()

	use Globals
	use ValueModule

	implicit none
    real(8) ll
	
    open(1, file='Input\V.txt', status='old')
	open(2, file='Input\AS.txt', status='old')

    do ia = 1, na
	do ieps = 1, neps		
    do ibta = 1, nbta
	do ik = 1, nk
	do iz = 1, nz
	

		 read(1, *) V1(ia,ieps,ibta,ik,iz)
		 read(2, *) AS1(ia,ieps,ibta,ik,iz)
		
    end do
	end do
	end do
	end do
	end do 

	
   	!initialize value function and decision rules

    mubar = (mu(1)+mu(2))/2.0

	do ik = 1, nk

		do ia = 1, na
		do iPg = 1, nP
        do iPb = 1, nP
		do ibta = 1, nbta

			AS(ia,iPg,iPb,ibta,ik) = (AS1(ia,1,ibta,ik,1)+AS1(ia,1,ibta,ik,2)+AS1(ia,2,ibta,ik,1)+AS1(ia,2,ibta,ik,2))/4.0
			V(ia,iPg,iPb,ibta,ik) = (V1(ia,1,ibta,ik,1)+V1(ia,1,ibta,ik,2)+V1(ia,2,ibta,ik,1)+V1(ia,2,ibta,ik,2))/4.0
		
		end do				
        end do  
		end do
       
	end do
	end do


    open(1, file='Output\Vinitial.txt', status='unknown')
	open(2, file='Output\ASinitial.txt', status='unknown')
	

    do ia = 1, na
	do iPg = 1, nP
	do iPb = 1, nP		
    do ibta = 1, nbta
	do ik = 1, nk

	

		write(1, '(F15.6)') V(ia,iPg,iPb,ibta,ik)
		write(2, '(F15.6)') AS(ia,iPg,iPb,ibta,ik) 
		
    end do
	end do
	end do
	end do
	end do 
	



end subroutine

