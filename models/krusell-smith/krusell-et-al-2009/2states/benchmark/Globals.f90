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
	real(8), parameter:: hour=0.3271
	!	number of grids

	integer, parameter:: na = 45, neps = 2, nk = 4, nz = 2, nbta = 3



	!	maximum number of iterations

	integer, parameter:: maxiterLOM = 30, showerr = 50
	integer iterLOM, MAXITER

	!	tolerence for convergence check

	real(8), parameter:: tol_value = 1.0D-5, tol_LOM = 1.0D-4


	!	number of periods, individuals for simulation

	integer, parameter:: Nperiod = 11000, Nindiv = 90000, Nskip = 1000


	!	adjustment factor

	real(8) sfac 

	
	!	transition probability of aggregate shock
	! ****   expected duration is 2 years. 8=(1/(1-p)) gives us p=7/8=0.875  ****
  
	real(8), parameter:: trProb = 0.875D0


	!	model parameters

	real(8), parameter:: alpha = 0.36d0, delta = 0.025d0, ui=0.0334
	real(8) bta(nbta)
    
	! **** Unemployment rate for good and bad times ****
    real(8)  mu(nz)
    
	 
	!	auxiliary variables
	real(8) irate, wage, income, employ
    integer goodtime, badtime

	!	grid vectors
	real(8) agrid(na), epsgrid(neps), kgrid(nk), zgrid(nz), btagrid(nbta)
	

	!	transition probability
	real(8) trEps(nz,nz,neps,neps), CtrEps(nz,nz,neps,neps), &
	        trZ(nz,nz), CtrZ(nz,nz), trBta(nbta,nbta), CtrBta(nbta,nbta)


	!	simulated time series data
	real(8) Kdata(Nperiod+1), Ldata(Nperiod),  Zdata(Nperiod+1), Ydata(Nperiod), &
	        Idata(Nperiod), Cdata(Nperiod), Adata(Nperiod), Edata(Nperiod), Rdata(Nperiod), Wdata(Nperiod)
	integer iZdata(Nperiod+1)

	!  cross section
	real(8) Adist(Nindiv), Gini(Nperiod), wealth1(Nperiod), wealth5(Nperiod), wealth10(Nperiod), &
	        wealth20(Nperiod), wealth30(Nperiod), Negratio(Nperiod)


	!	coefficients of equations
	real(8)	KcoefG(2), NewKcoefG(2), KcoefB(2), NewKcoefB(2)


	!	indexing variables
	 !**** iz=1 is bad times and iz=2 is good times ****	
	integer ia, iz, ik, ieps, ibta, indiv, time
	logical final


	!	spaces for value functions and decision rules
	
	real(8) V(na,neps,nbta,nk,nz), EV(na,neps,nbta,nk,nz), AS(na,neps,nbta,nk,nz), EV2p(na,neps,nbta,nk,nz)

	real(8) PanelA112G(Nindiv), PanelA112B(Nindiv), PanelA123G(Nindiv), PanelA123B(Nindiv)
	integer PanelEps112G(Nindiv), PanelEps112B(Nindiv), PanelEps123G(Nindiv), PanelEps123B(Nindiv)
	integer PanelBta112G(Nindiv), PanelBta112B(Nindiv), PanelBta123G(Nindiv), PanelBta123B(Nindiv)




end module
