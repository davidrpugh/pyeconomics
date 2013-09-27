!	----------------------------------------------------------------------
!	File name: SimulateData.f90
!
!	generate artificial data through simulating individuals' decision rules. These data
!	will be used in RegressLOM().
!	----------------------------------------------------------------------


subroutine SimulateData()

	use Globals
	use LinInterpModule
	use PolyInterpModule
	use Numerical_Libraries
	
	implicit none

	integer, parameter:: xseed = 2

	real(8), dimension(:,:,:), allocatable:: ASkz
	real(8), dimension(:,:), allocatable:: CrA
	integer, dimension(:,:), allocatable:: CrEps, CrBta

   
	!	initialize time series data

	Kdata = 0.0D0
	Ldata = 0.0D0


	!	set the seed number for random number generator

	call rnset(xseed)


	!	allocate space for cross section data generation

	allocate(CrA(Nindiv,2), CrEps(Nindiv,2), CrBta(Nindiv,2),		&
			 ASkz(na,neps,nbta))


	!	intialize distribution of workers

	do indiv = 1, Nindiv
		CrA(indiv,1) = kss
		CrBta(indiv,1) = 2
		if (indiv<(Nindiv)*(1.0-mu(iZdata(1)))) then
	    	CrEps(indiv,1) = 4
		else if (indiv>(Nindiv)*(1.0-mu(iZdata(1)))+(Nindiv)*mu(iZdata(1))*(1/3)) then
	    	CrEps(indiv,1) = 3
		else 
            CrEps(indiv,1) = 1
		end if
        if (mod(indiv,10)==1) then 
		    CrBta(indiv,1) = 1
        else if (mod(indiv,10)==2) then 
	    	CrBta(indiv,1) = 3
		end if  
	end do

	!	initial aggregate capital

	Kdata(1) = sum(CrA(:,1))/Nindiv
	

	!	start generating artificial time series data

	do time = 1, Nperiod


		! interpolate value function and asset decision rule at (Kapital, z)
        ! **** compute the value function and decision rules for each a,eps and eta ****

		do ia = 1, na
        do ieps = 1, neps
        do ibta = 1, nbta
		
			ASkz(ia,ieps,ibta) = polyinterp1(Kdata(time), kgrid, AS(ia,ieps,ibta,:,iZdata(time)))
		
		end do	
		end do
        end do

		!	cross-section data for asset and hours given (Kapital, z)

		do indiv = 1, Nindiv

				CrA(indiv,2) = lininterp1(CrA(indiv,1), agrid, ASkz(:, CrEps(indiv,1),CrBta(indiv,1)))
				CrA(indiv,2) = min(max(CrA(indiv,2), agrid(1)), agrid(na))
				Ldata(time) = Ldata(time) + epsgrid(CrEps(indiv,1))*hour

		end do
		

	
		
		!	generate time series data

		Ldata(time) = Ldata(time)/Nindiv
		Rdata(time) = Zdata(time)*alpha*(Kdata(time)/Ldata(time))**(alpha-1.0D0)
		Wdata(time) = Zdata(time)*(1.0D0-alpha)*(Kdata(time)/Ldata(time))**(alpha)

		!	next period variables

		Kdata(time+1) = min(max(sum(CrA(:,2))/Nindiv, kgrid(1)), kgrid(nk))

		call NextBeta(CrBta)

		call NextIdiosyncraticShock(CrEps)

		
		!	lag cross section data one period

		CrA(:,1) = CrA(:,2)
		CrEps(:,1) = CrEps(:,2)
        CrBta(:,1) = CrBta(:,2)

		!	show how simulation is going on

		if (mod(time, 100) == 0) then
			write(*,'(A,I6)')   "time = ", time	
			write(*,'(A,F8.5)') "z = ", Zdata(time)
			write(*,'(A,F12.6)') "K = ", Kdata(time)
			write(*,'(A,F12.6)') "Kp = ", Kdata(time+1)
			write(*,'(A,F8.5)') "L = ", Ldata(time)
			write(*,'(/)')
		end if

	end do
	

	!	deallocate memory

	deallocate(CrA, CrEps, ASkz, CrBta)


	!	save time sries data

	open(1, file='Output\Kdata.txt', status='unknown')
	open(2, file='Output\Ldata.txt', status='unknown')
	
	close(1, status='delete')
	close(2, status='delete')

	open(1, file='Output\Kdata.txt', status='unknown')
	open(2, file='Output\Ldata.txt', status='unknown')

	write(1, '(f12.6)') Kdata
	write(2, '(f12.6)') Ldata


end subroutine
