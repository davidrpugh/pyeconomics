!	----------------------------------------------------------------------
!	File name: FinalSimulation.f90
!
!	generate artificial data for all aggregate variables, which will be 
!	used in ModelStatistics().
!	----------------------------------------------------------------------


subroutine FinalSimulation()

	use Globals
	use LinInterpModule
	use PolyInterpModule
	use Numerical_Libraries

	implicit none
    integer, parameter:: xseed = 2
	! **** Allocate space for Akz and Vkz ****
	! **** allocate space for CrA, CrEps ****
	real(8), dimension(:,:,:), allocatable:: ASkz
	real(8), dimension(:,:), allocatable:: CrA
	integer, dimension(:,:), allocatable:: CrEps, CrBta
	real(8) avgG, avgB, varG, varB, errorG, errorB, sumN, sumD, sum1, sum5, sum10, sum20, sum30
	integer Negsum

	!	save final output

	final = .true.


	!	write the cnverged coefficients of equations

	open(0, file='Output\Logfile.txt', status='unknown')

	write(0, '(//A)') "Converged Coefficients of Equations:"
	write(0, '(/A, 2F10.5)') "Law of Motion for KapitalG:", KcoefG
	write(0, '(/A, 2F10.5)') "Law of Motion for KapitalB:", KcoefB

	call SaveCoefficients()


	!	solve value function
	
	call SolveValueFunction()


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
	    	CrEps(indiv,1) = 2
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


	!calculate statistics



	call DSVRGN(Nindiv,CrA(:,1),Adist)

	sumN=0.0D0
	sumD=0.0D0
	do indiv=1, Nindiv
	sumN=sumN+(real(Nindiv)+1.0-real(indiv))*Adist(indiv)
	sumD=sumD+Adist(indiv)
	end do

	Gini(time)=(real(Nindiv)+1.0-2.0*(sumN/sumD))/real(Nindiv)
	


	
	sum1=0.0D0
	do indiv=1, int(0.01*real(Nindiv))
	sum1=sum1+Adist(Nindiv-indiv+1)
	end do
	wealth1(time)=sum1/sumD

	sum5=0.0D0
	do indiv=1, int(0.05*real(Nindiv))
	sum5=sum5+Adist(Nindiv-indiv+1)
	end do
	wealth5(time)=sum5/sumD
	

	sum10=0.0D0
	do indiv=1, int(0.1*real(Nindiv))
	sum10=sum10+Adist(Nindiv-indiv+1)
	end do
	wealth10(time)=sum10/sumD
	

	sum20=0.0D0
	do indiv=1, int(0.2*real(Nindiv))
	sum20=sum20+Adist(Nindiv-indiv+1)
	end do
	wealth20(time)=sum20/sumD

	sum30=0.0D0
	do indiv=1, int(0.3*real(Nindiv))
	sum30=sum30+Adist(Nindiv-indiv+1)
	end do
	wealth30(time)=sum30/sumD

	Negsum=0
	do indiv=1, Nindiv
	if (Adist(indiv)<0) then
	Negsum=Negsum+1
	end if
	end do
	Negratio(time)=real(Negsum)/real(Nindiv)


	if (time==10941) then !k=11.2, bad  10941
	PanelA112B = CrA(:,1)
	PanelEps112B = CrEps(:,1)
	PanelBta112B = CrBta(:,1)
	end if

	if (time==9606) then !k=11.2, good 9606
	PanelA112G = CrA(:,1)
	PanelEps112G = CrEps(:,1)
	PanelBta112G = CrBta(:,1)
	end if

	if (time==7562) then !k=12.3, bad 7562
	PanelA123B = CrA(:,1)
	PanelEps123B = CrEps(:,1)
	PanelBta123B = CrBta(:,1)
	end if

	if (time==9080) then !k=12.3, good 9080
	PanelA123G = CrA(:,1)
	PanelEps123G = CrEps(:,1)
	PanelBta123G = CrBta(:,1)
	end if

	end do
	


	!	deallocate memory

	deallocate(CrA, CrEps, ASkz, CrBta)




	!	save simulated time series data

	call SaveSimulatedTimeSeriesData()


	avgG=0.0D0
	avgB=0.0D0

	do time=Nskip+1, Nperiod
	if (iZdata(time)==2) then
	avgG=avgG+dlog(Kdata(time+1))/real(goodtime)
	else
	avgB=avgB+dlog(Kdata(time+1))/real(badtime)
	end if
	end do

	varG=0.0D0
	varB=0.0D0

	do time=Nskip+1, Nperiod
	if (iZdata(time)==2) then
	varG=varG+(dlog(Kdata(time+1))-avgG)**2
	else
	varB=varB+(dlog(Kdata(time+1))-avgB)**2
	end if
	end do

	errorG=0.0D0
	errorB=0.0D0

	do time=Nskip+1, Nperiod
	if (iZdata(time)==2) then
	errorG=errorG+(dlog(Kdata(time+1))-KcoefG(1)-KcoefG(2)*dlog(Kdata(time)))**2
	else
	errorB=errorB+(dlog(Kdata(time+1))-KcoefB(1)-KcoefB(2)*dlog(Kdata(time)))**2
	end if
	end do


	open(1, file='Output\RsqG.txt', status='unknown')
 
	write(1, '(F15.10)') 1.0-errorG/varG

	open(2, file='Output\RsqB.txt', status='unknown')
 
	write(2, '(F15.10)') 1.0-errorB/varB

end subroutine
