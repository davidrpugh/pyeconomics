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

	Kcoef(1) = 0.04058 
	Kcoef(2) = 0.98339 
   


    open(1, file='Input\K.txt', status='old')


	do ik = 1, Nperiod+1

	

		 read(1, *) K(ik)

		
	end do 



    open(2, file='Output\K.txt', status='unknown')
	write(2, '(F15.6)') K



end subroutine




subroutine InitializeValue2()

  use Globals
  use ValueModule
  
  implicit none
 
  ik=2
  do time = 1, N1
        
     
	 irate2(time) = zseries(time)*alpha*(real(K(time)/((1-useries(time))*hour))**(alpha-1))-delta
	 wage2(time) = zseries(time)*(1-alpha)*(real(K(time)/((1.0-useries(time))*hour))**(alpha))
   


		do ia = 1, na
		do iPge = 1, nP
        do iPbe = 1, nP
		do iPgs = 1, nP
        do iPbs = 1, nP	
		do iPbf = 1, nP

			AS2l(time,ia,iPge,iPbe,iPgs,iPbs,iPbf) = (AS1(ia,1,1,ik,3)+AS1(ia,1,1,ik,2)+AS1(ia,4,1,ik,3)+AS1(ia,4,1,ik,2))/4.0
			AS2m(time,ia,iPge,iPbe,iPgs,iPbs,iPbf) = (AS1(ia,1,2,ik,3)+AS1(ia,1,2,ik,2)+AS1(ia,4,2,ik,3)+AS1(ia,4,2,ik,2))/4.0
			AS2h(time,ia,iPge,iPbe,iPgs,iPbs,iPbf) = (AS1(ia,1,3,ik,3)+AS1(ia,1,3,ik,2)+AS1(ia,4,3,ik,3)+AS1(ia,4,3,ik,2))/4.0

			V2l(time,ia,iPge,iPbe,iPgs,iPbs,iPbf) = (V1(ia,1,1,ik,3)+V1(ia,1,1,ik,2)+V1(ia,4,1,ik,3)+V1(ia,4,1,ik,2))/4.0
			V2m(time,ia,iPge,iPbe,iPgs,iPbs,iPbf) = (V1(ia,1,2,ik,3)+V1(ia,1,2,ik,2)+V1(ia,4,2,ik,3)+V1(ia,4,2,ik,2))/4.0
			V2h(time,ia,iPge,iPbe,iPgs,iPbs,iPbf) = (V1(ia,1,3,ik,3)+V1(ia,1,3,ik,2)+V1(ia,4,3,ik,3)+V1(ia,4,3,ik,2))/4.0
		
        end do  
		end do
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
	do iz = 1, nz+1
	

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
		do iPge = 1, nP
        do iPbe = 1, nP
		do iPgs = 1, nP
        do iPbs = 1, nP	
		do iPbf = 1, nP

			ASlo(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik) = (AS1(ia,1,1,ik,3)+AS1(ia,1,1,ik,2))*0.1+(AS1(ia,4,1,ik,3)+AS1(ia,4,1,ik,2))*0.4
			ASmi(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik) = (AS1(ia,1,2,ik,3)+AS1(ia,1,2,ik,2))*0.1+(AS1(ia,4,2,ik,3)+AS1(ia,4,2,ik,2))*0.4
			AShi(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik) = (AS1(ia,1,3,ik,3)+AS1(ia,1,3,ik,2))*0.1+(AS1(ia,4,3,ik,3)+AS1(ia,4,3,ik,2))*0.4

			Vl(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik) = (V1(ia,1,1,ik,3)+V1(ia,1,1,ik,2))*0.1+(V1(ia,4,1,ik,3)+V1(ia,4,1,ik,2))*0.4
			Vm(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik) = (V1(ia,1,2,ik,3)+V1(ia,1,2,ik,2))*0.1+(V1(ia,4,2,ik,3)+V1(ia,4,2,ik,2))*0.4
			Vh(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik) = (V1(ia,1,3,ik,3)+V1(ia,1,3,ik,2))*0.1+(V1(ia,4,3,ik,3)+V1(ia,4,3,ik,2))*0.4


				
        end do  
		end do
       	end do
        end do  
		end do
       	end do
	end do


    open(1, file='Output\Vlinitial.txt', status='unknown')
    open(2, file='Output\Vminitial.txt', status='unknown')
	open(3, file='Output\Vhinitial.txt', status='unknown')

	open(4, file='Output\ASloinitial.txt', status='unknown')
	open(5, file='Output\ASmiinitial.txt', status='unknown')	
	open(6, file='Output\AShiinitial.txt', status='unknown')


    do ia = 1, na
	do iPge = 1, nP
	do iPbe = 1, nP		
	do iPgs = 1, nP
	do iPbs = 1, nP	
	do iPbf = 1, nP
	do ik = 1, nk

		write(1, '(F15.6)') Vl(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik)
		write(2, '(F15.6)') Vm(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik)
		write(3, '(F15.6)') Vh(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik)
		write(4, '(F15.6)') ASlo(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik) 
		write(5, '(F15.6)') ASmi(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik) 
		write(6, '(F15.6)') AShi(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik) 
				
    end do
	end do
	end do
	end do
	end do 
	end do
	end do 
		



end subroutine

