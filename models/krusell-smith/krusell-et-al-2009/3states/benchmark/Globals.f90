!	----------------------------------------------------------------------
!	File name : Globals.f90
!
!	Defines Global Variables
!	----------------------------------------------------------------------



module Globals

	
implicit none

	
	!	bounds for individual and aggregate capital

	real(8), parameter:: amin = -2.75D0,		&
						 amax = 3000.0D0		


	real(8) kss, kmax, kmin
	real(8), parameter:: hour=0.3271, cutoff=-1.0


	!	number of grids

	integer, parameter:: na = 50, neps = 4, nk = 4, nz = 3, nbta = 3



	!	maximum number of iterations

	integer, parameter:: maxiterLOM = 30, showerr = 50
	integer iterLOM, MAXITER, iter

	!	tolerence for convergence check

	real(8), parameter:: tol_value = 1.0D-5, tol_LOM = 1.0D-4


	!	number of periods, individuals for simulation

	integer, parameter:: Nperiod = 11000, Nindiv = 90000, Nskip = 1000


	!	adjustment factor

	real(8) sfac 

	
	!	transition probability of aggregate shock
	! ****   expected duration is 2 years. 8=(1/(1-p)) gives us p=7/8=0.875  ****
  
!	real(8), parameter:: trProb = 0.875D0


	!	model parameters

	real(8), parameter:: alpha = 0.36d0, delta = 0.025d0
	real(8) bta(nbta)
    
	! **** Unemployment rate for good and bad times ****
    real(8)  mu(nz), ui(neps)
    
	 
	!	auxiliary variables
	real(8) irate, wage, income, employ
    integer goodtime, badtime, downtime

	!	grid vectors
	real(8) agrid(na), epsgrid(neps), kgrid(nk), zgrid(nz), btagrid(nbta)
	

	!	transition probability
	real(8) trEps(nz,nz,neps,neps), CtrEps(nz,nz,neps,neps), &
	        trZ(nz,nz), CtrZ(nz,nz), trBta(nbta,nbta), CtrBta(nbta,nbta)


	!	simulated time series data
	real(8) Kdata(Nperiod+1), Ldata(Nperiod),  Zdata(Nperiod+1), Ydata(Nperiod), AvgDur(Nperiod), ratioover2(Nperiod), AvgDurNonzero(Nperiod),&
	        Idata(Nperiod), Cdata(Nperiod), Adata(Nperiod), Edata(Nperiod), Rdata(Nperiod), Wdata(Nperiod),&
			Empdata(Nperiod), Shodata(Nperiod), Firdata(Nperiod), Londata(Nperiod)
	integer iZdata(Nperiod+1)

	!  cross section
	real(8) Adist(Nindiv), Gini(Nperiod), wealth1(Nperiod), wealth5(Nperiod), wealth10(Nperiod), &
	        wealth20(Nperiod), wealth30(Nperiod), Negratio(Nperiod)


	!	coefficients of equations
	real(8)	KcoefG(2), NewKcoefG(2), KcoefB(2), NewKcoefB(2), KcoefD(2), NewKcoefD(2)


	!	indexing variables
	 !**** iz=1 is bad times and iz=2 is good times ****	
	integer ia, iz, ik, ieps, ibta, indiv, time
	logical final


	!	spaces for value functions and decision rules
	
	real(8) V(na,neps,nbta,nk,nz), EV(na,neps,nbta,nk,nz), AS(na,neps,nbta,nk,nz), EV2p(na,neps,nbta,nk,nz)

	real(8) PanelA113G(Nindiv), PanelA113B(Nindiv), PanelA113D(Nindiv), PanelA121G(Nindiv), PanelA121B(Nindiv), PanelA121D(Nindiv)
	integer PanelEps113G(Nindiv), PanelEps113B(Nindiv), PanelEps113D(Nindiv), PanelEps121G(Nindiv), PanelEps121B(Nindiv), PanelEps121D(Nindiv)
	integer PanelBta113G(Nindiv), PanelBta113B(Nindiv), PanelBta113D(Nindiv), PanelBta121G(Nindiv), PanelBta121B(Nindiv), PanelBta121D(Nindiv)
	integer PanelDur113G(Nindiv), PanelDur113B(Nindiv), PanelDur113D(Nindiv), PanelDur121G(Nindiv), PanelDur121B(Nindiv), PanelDur121D(Nindiv)





end module
