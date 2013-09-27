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
	integer, dimension(:,:), allocatable:: CrEps, CrBta, CrDur
	real(8) avgG, avgB, avgD, varG, varB, varD, errorG, errorB, errorD, sumN, sumD, sum1, sum5, sum10, sum20, sum30
	integer Negsum, Numunemp, NumDurnonzero, Numover2

	!	save final output

	final = .true.


	!	write the cnverged coefficients of equations

	open(0, file='Output\Logfile.txt', status='unknown')

	write(0, '(//A)') "Converged Coefficients of Equations:"
	write(0, '(/A, 2F10.5)') "Law of Motion for KapitalG:", KcoefG
	write(0, '(/A, 2F10.5)') "Law of Motion for KapitalB:", KcoefB
	write(0, '(/A, 2F10.5)') "Law of Motion for KapitalD:", KcoefD	

	call SaveCoefficients()


	!	solve value function
	
	call SolveValueFunction()


	!	initialize time series data

	Kdata = 0.0D0
	Ldata = 0.0D0
	AvgDur = 0.0D0
	AvgDurNonzero = 0.0D0
	ratioover2 = 0.0D0
	NumUnemp=0
	NumDurNonzero = 0
	Numover2 = 0

	!	set the seed number for random number generator

	call rnset(xseed)


	!	allocate space for cross section data generation

	allocate(CrA(Nindiv,2), CrEps(Nindiv,2), CrBta(Nindiv,2),		&
			 CrDur(Nindiv,2), ASkz(na,neps,nbta))


	!	intialize distribution of workers

	do indiv = 1, Nindiv
		CrA(indiv,1) = kss
		CrBta(indiv,1) = 2
		if (indiv<(Nindiv)*(1.0-mu(iZdata(1)))) then
	    	CrEps(indiv,1) = 4
			CrDur(indiv,1) = 0
		else if (indiv>(Nindiv)*(1.0-mu(iZdata(1)))+(Nindiv)*mu(iZdata(1))*(1/3)) then
	    	CrEps(indiv,1) = 3
			CrDur(indiv,1) = 1
		else 
            CrEps(indiv,1) = 1
			CrDur(indiv,1) = 1
		end if
        if (mod(indiv,10)==1) then 
		    CrBta(indiv,1) = 1
        else if (mod(indiv,10)==2) then 
	    	CrBta(indiv,1) = 3
        end if  
	end do		



	!	initial aggregate capital

	Kdata(1) = sum(CrA(:,1))/Nindiv
	Empdata=0.0
	Shodata=0.0
	Firdata=0.0
	Londata=0.0

	!	start generating artificial time series data

	do time = 1, Nperiod

	NumUnemp=0
	NumDurNonzero = 0
	Numover2 = 0
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

			if (CrEps(indiv,1)==4) then
				Empdata(time)=Empdata(time)+1.0/real(Nindiv)
			else if (CrEps(indiv,1)==3) then
				Shodata(time)=Shodata(time)+1.0/real(Nindiv)
			else if (CrEps(indiv,1)==2) then
				Firdata(time)=Firdata(time)+1.0/real(Nindiv)
			else
				Londata(time)=Londata(time)+1.0/real(Nindiv)
			end if


		end do
		

	
		
		!	generate time series data

		Ldata(time) = Ldata(time)/Nindiv
		Rdata(time) = Zdata(time)*alpha*(Kdata(time)/Ldata(time))**(alpha-1.0D0)
		Wdata(time) = Zdata(time)*(1.0D0-alpha)*(Kdata(time)/Ldata(time))**(alpha)


		!	next period variables

		Kdata(time+1) = min(max(sum(CrA(:,2))/Nindiv, kgrid(1)), kgrid(nk))

		call NextBeta(CrBta)

		call NextIdiosyncraticShock(CrEps)
		do indiv = 1, Nindiv

			if (CrEps(indiv,2)==4) then
				CrDur(indiv,2)=0
			else
				CrDur(indiv,2)=CrDur(indiv,1)+1
				
			end if

			if (CrEps(indiv,1)==4) then
				Numunemp=Numunemp
			else
				Numunemp=Numunemp+1
			end if

		AvgDur(time)=AvgDur(time)+real(CrDur(indiv,1))/real(Nindiv)
		NumDurNonzero=NumDurNonzero+CrDur(indiv,1)
		if (CrDur(indiv,1)>2) then
		Numover2=Numover2+1
		end if
		end do

		AvgDurNonzero(time)=real(NumDurNonzero)/real(NumUnemp)
		ratioover2(time)=real(Numover2)/real(NumUnemp)

	if (time==10912) then !k=11.3, B
	PanelA113B = CrA(:,1)
	PanelEps113B = CrEps(:,1)
	PanelBta113B = CrBta(:,1)
	PanelDur113B = CrDur(:,1)

	end if

	if (time==10252) then !k=11.3, D
	PanelA113D = CrA(:,1)
	PanelEps113D = CrEps(:,1)
	PanelBta113D = CrBta(:,1)
	PanelDur113D = CrDur(:,1)

	end if

	if (time==10155) then !k=11.3, G
	PanelA113G = CrA(:,1)
	PanelEps113G = CrEps(:,1)
	PanelBta113G = CrBta(:,1)
	PanelDur113G = CrDur(:,1)

	end if

	if (time==9316) then !k=12.1, B
	PanelA121B = CrA(:,1)
	PanelEps121B = CrEps(:,1)
	PanelBta121B = CrBta(:,1)
	PanelDur121B = CrDur(:,1)
	end if


	if (time==7045) then !k=12.1, D
	PanelA121D = CrA(:,1)
	PanelEps121D = CrEps(:,1)
	PanelBta121D = CrBta(:,1)
	PanelDur121D = CrDur(:,1)
	end if

	if (time==10471) then !k=12.1, G
	PanelA121G = CrA(:,1)
	PanelEps121G = CrEps(:,1)
	PanelBta121G = CrBta(:,1)
	PanelDur121G = CrDur(:,1)
	end if

		
		!	lag cross section data one period

		CrA(:,1) = CrA(:,2)
		CrEps(:,1) = CrEps(:,2)
        CrBta(:,1) = CrBta(:,2)
        CrDur(:,1) = CrDur(:,2)

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











	end do
	


	!	deallocate memory

	deallocate(CrA, CrEps, ASkz, CrBta, CrDur)




	!	save simulated time series data

	call SaveSimulatedTimeSeriesData()


	avgG=0.0D0
	avgB=0.0D0
	avgB=0.0D0

	do time=Nskip+1, Nperiod
	if (iZdata(time)==3) then
	avgG=avgG+dlog(Kdata(time+1))/real(goodtime)
	else if (iZdata(time)==2) then
	avgD=avgD+dlog(Kdata(time+1))/real(downtime)
	else
	avgB=avgB+dlog(Kdata(time+1))/real(badtime)
	end if
	end do

	varG=0.0D0
	varB=0.0D0
	varD=0.0D0

	do time=Nskip+1, Nperiod
	if (iZdata(time)==3) then
	varG=varG+(dlog(Kdata(time+1))-avgG)**2
	else if (iZdata(time)==2) then
	varD=varD+(dlog(Kdata(time+1))-avgD)**2
	else
	varB=varB+(dlog(Kdata(time+1))-avgB)**2
	end if
	end do

	errorG=0.0D0
	errorB=0.0D0
	errorD=0.0D0

	do time=Nskip+1, Nperiod
	if (iZdata(time)==3) then
	errorG=errorG+(dlog(Kdata(time+1))-KcoefG(1)-KcoefG(2)*dlog(Kdata(time)))**2
	else if (iZdata(time)==2) then
	errorD=errorD+(dlog(Kdata(time+1))-KcoefD(1)-KcoefD(2)*dlog(Kdata(time)))**2
	else
	errorB=errorB+(dlog(Kdata(time+1))-KcoefB(1)-KcoefB(2)*dlog(Kdata(time)))**2
	end if
	end do


	open(1, file='Output\RsqG.txt', status='unknown')
 
	write(1, '(F15.10)') 1.0-errorG/varG

	open(2, file='Output\RsqB.txt', status='unknown')
 
	write(2, '(F15.10)') 1.0-errorB/varB

	open(3, file='Output\RsqD.txt', status='unknown')
 
	write(3, '(F15.10)') 1.0-errorD/varD

end subroutine
