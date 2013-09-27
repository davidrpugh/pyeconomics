!	----------------------------------------------------------------------
!	File name: FinalSimulation.f90
!
!	generate artificial data through simulating individuals' decision rules. These data
!	will be used in RegressLOM().
!	----------------------------------------------------------------------


subroutine FinalSimulation()

	use Globals
	use LinInterpModule
	use PolyInterpModule
	use Numerical_Libraries
		 include 'link_fnl_static.h'
   !DEC$ OBJCOMMENT LIB:'libguide.lib'

	
	implicit none

	integer, parameter:: xseed = 2
	real(8), parameter:: tiny = 1.0D-20
	integer Z1
	real(8) avgG, avgB, varG, varB, errorG, errorB, sumN, sumD, sum1, sum5, sum10, sum20, sum30
	integer Negsum

	real(8), dimension(:,:,:,:,:,:,:), allocatable:: ASkz
	real(8), dimension(:,:), allocatable:: CrA, CrPge, CrPbe, CrPgs, CrPbs, CrPbf
	integer, dimension(:,:), allocatable:: CrBta, CrEps


	!	save final output

	final = .true.


	!	write the cnverged coefficients of equations

	open(0, file='Output\Logfile.txt', status='unknown')

	write(0, '(//A)') "Converged Coefficients of Equations:"
	write(0, '(/A, 2F10.5)') "Law of Motion for Kapital(1):", 	Kcoef(1) 
	write(0, '(/A, 2F10.5)') "Law of Motion for Kapital(2):", 	Kcoef(2) 

	call SaveCoefficients()


	!	solve value function
	
	call SolveValueFunction()

   
	!	initialize time series data

	Kdata = 0.0D0
	Ldata = 0.0D0


	!	set the seed number for random number generator

	call rnset(xseed)


	!	allocate space for cross section data generation

	allocate(CrA(Nindiv,2), CrEps(Nindiv,2), CrPge(Nindiv,2), CrPbe(Nindiv,2), CrPgs(Nindiv,2), CrPbs(Nindiv,2), CrPbf(Nindiv,2), CrBta(Nindiv,2),		&
			 ASkz(na,nP,nP,nP,nP,nP,nbta))
      
	! Initial state assumed to be good. Also change this in generateseries if the initial state is a bad state
	! CHANGE TO B STATE
	Zdata(1) = 1
	
	!   Read the initial distribution from the file 
	open(1, file='Input\PanelA113B.txt', status = 'old')
    do indiv = 1, Nindiv
     read(1, *) CrA(indiv,1)
      CrA(indiv,2) = 0.0D0
    end do


    open(1, file='Input\PanelEps113B.txt', status = 'old')
    do indiv = 1, Nindiv
      read(1, *) CrEps(indiv,1)
      CrEps(indiv,2) = 0.0D0
    end do
 
    open(1, file='Input\PanelBta113B.txt', status = 'old')
    do indiv = 1, Nindiv
      read(1, *) CrBta(indiv,1)
      CrBta(indiv,2) = 0.0D0
    end do
	

 
	CrPge = 0.0
	CrPbe = 0.0
	CrPgs = 0.0
	CrPbs = 0.0
	CrPbf = 0.0

    do indiv = 1, Nindiv

     if (CrEps(indiv,1) == 4) then 
	     CrPbe(indiv,1) = 1.0
	 else if (CrEps(indiv,1)== 3) then
		 CrPbs(indiv,1) = 1.0
	 end if
   
    end do

	!	initial aggregate capital

	Kdata(1) = real(sum(CrA(:,1))/Nindiv)



	!	start generating artificial time series data

	do time = 1, N1

		
		do indiv = 1, Nindiv
        if (CrBta(indiv,1)==1) then
		  CrA(indiv,2) = lininterp6(CrA(indiv,1), CrPge(indiv,1), CrPbe(indiv,1), CrPgs(indiv,1), CrPbs(indiv,1), CrPbf(indiv,1), agrid, Pgrid, Pgrid, Pgrid, Pgrid, Pgrid, AS2l(time,:,:,:,:,:,:))
		  CrA(indiv,2) = min(max(CrA(indiv,2), agrid(1)), agrid(na))
		else if (CrBta(indiv,1)==2) then
		  CrA(indiv,2) = lininterp6(CrA(indiv,1), CrPge(indiv,1), CrPbe(indiv,1), CrPgs(indiv,1), CrPbs(indiv,1), CrPbf(indiv,1), agrid, Pgrid, Pgrid, Pgrid, Pgrid, Pgrid, AS2m(time,:,:,:,:,:,:))
		  CrA(indiv,2) = min(max(CrA(indiv,2), agrid(1)), agrid(na))		
		else
		  CrA(indiv,2) = lininterp6(CrA(indiv,1), CrPge(indiv,1), CrPbe(indiv,1), CrPgs(indiv,1), CrPbs(indiv,1), CrPbf(indiv,1), agrid, Pgrid, Pgrid, Pgrid, Pgrid, Pgrid, AS2h(time,:,:,:,:,:,:))
		  CrA(indiv,2) = min(max(CrA(indiv,2), agrid(1)), agrid(na))
		end if
				 
		end do
     
		!	next period variables

		Kdata(time+1) = sum(CrA(:,2))/Nindiv

		
		call NextBeta(CrBta)
		call NextP2(CrPge,CrPbe,CrPgs,CrPbs,CrPbf)
       if (time==1) then
	   	open(1, file='Output\CrA.txt', status='unknown')
	    open(2, file='Output\CrPge.txt', status='unknown')
	    open(3, file='Output\CrPbe.txt', status='unknown')
	    open(4, file='Output\CrPgs.txt', status='unknown')
	    open(5, file='Output\CrPbs.txt', status='unknown')
	    open(6, file='Output\CrPbf.txt', status='unknown')


		 write(1, '(F15.6)') CrA(:,2)
		 write(2, '(F15.6)') CrPge(:,2)
		 write(3, '(F15.6)') CrPbe(:,2)
		 write(4, '(F15.6)') CrPgs(:,2)
		 write(5, '(F15.6)') CrPbs(:,2)		
		 write(6, '(F15.6)') CrPbf(:,2)
	   
	   end if  
	
		!	lag cross section data one period

		CrA(:,1) = CrA(:,2)
        CrBta(:,1) = CrBta(:,2)
		CrPge(:,1) = CrPge(:,2)
		CrPbe(:,1) = CrPbe(:,2)
		CrPgs(:,1) = CrPgs(:,2)
		CrPbs(:,1) = CrPbs(:,2)
		CrPbf(:,1) = CrPbf(:,2)


		!	show how simulation is going on

		if (mod(time, 100) == 0) then
			write(*,'(A,I6)')   "time = ", time	
			write(*,'(A,F12.6)') "K = ", Kdata(time)
			write(*,'(A,F12.6)') "Kp = ", Kdata(time+1)
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

	!	start generating artificial time series data

	do time = N1+1, Nperiod


	! interpolate value function and asset decision rule at (Kapital, z)
       

		do ia = 1, na
        do iPge = 1, nP
        do iPbe = 1, nP
		do iPgs = 1, nP
        do iPbs = 1, nP
		do iPbf = 1, nP
			ASkz(ia,iPge,iPbe,iPgs,iPbs,iPbf,1) = polyinterp1(Kdata(time), kgrid, ASlo(ia,iPge,iPbe,iPgs,iPbs,iPbf,:))
			ASkz(ia,iPge,iPbe,iPgs,iPbs,iPbf,2) = polyinterp1(Kdata(time), kgrid, ASmi(ia,iPge,iPbe,iPgs,iPbs,iPbf,:))
			ASkz(ia,iPge,iPbe,iPgs,iPbs,iPbf,3) = polyinterp1(Kdata(time), kgrid, AShi(ia,iPge,iPbe,iPgs,iPbs,iPbf,:))		
		end do
        end do
		end do	
		end do
        end do
		end do

		!	cross-section data for asset and hours given (Kapital)

		do indiv = 1, Nindiv

		  CrA(indiv,2) = lininterp6(CrA(indiv,1), CrPge(indiv,1), CrPbe(indiv,1), CrPgs(indiv,1), CrPbs(indiv,1), CrPbf(indiv,1), agrid, Pgrid, Pgrid, Pgrid, Pgrid, Pgrid, ASkz(:,:,:,:,:,:,CrBta(indiv,1)))
		  CrA(indiv,2) = min(max(CrA(indiv,2), agrid(1)), agrid(na))  

		end do
		
		Kdata(time+1) = min(max(sum(CrA(:,2))/Nindiv, kgrid(1)), kgrid(nk))
        
        
		call NextBeta(CrBta)
		call NextP(CrPge,CrPbe,CrPgs,CrPbs,CrPbf)

		
		!	lag cross section data one period

		CrA(:,1) = CrA(:,2)
        CrBta(:,1) = CrBta(:,2)
		CrPge(:,1) = CrPge(:,2)
		CrPbe(:,1) = CrPbe(:,2)
		CrPgs(:,1) = CrPgs(:,2)
		CrPbs(:,1) = CrPbs(:,2)
		CrPbf(:,1) = CrPbf(:,2)
		
		
		!	show how simulation is going on

		if (mod(time, 100) == 0) then
			write(*,'(A,I6)')   "time = ", time	
			write(*,'(A,F12.6)') "K = ", Kdata(time)
			write(*,'(A,F12.6)') "Kp = ", Kdata(time+1)
			write(*,'(/)')
		end if
  	!record panel data

		if (time == 600) then
		PanelA = CrA(:,1)
		PanelPge = CrPge(:,1)
		PanelPbe = CrPbe(:,1)
		PanelPgs = CrPgs(:,1)
		PanelPbs = CrPbs(:,1)
		PanelPbf = CrPbf(:,1)

		PanelBta = CrBta(:,1)
	
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

	deallocate(CrA, CrEps, CrPge, CrPbe, CrPgs, CrPbs, CrPbf, ASkz, CrBta)



	!	save simulated time series data

	call SaveSimulatedTimeSeriesData()



end subroutine
