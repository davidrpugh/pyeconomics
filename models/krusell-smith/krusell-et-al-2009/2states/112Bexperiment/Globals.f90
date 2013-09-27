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
    integer, parameter:: ncut=nz*nz*neps, nP=6,  ninterv=ncut+1


	!	maximum number of iterations

	integer, parameter:: maxiterLOM = 50, showerr = 50, maxiterK = 200
	integer iterLOM, MAXITER

	!	tolerence for convergence check

	real(8), parameter:: tol_value = 1.0D-5, tol_LOM = 1.0D-4, TolErrK = 1.0D-5
	real(8), parameter:: smallval=0.000000001



	!	number of periods, individuals for simulation

	integer, parameter:: Nperiod = 600, Nindiv = 90000, N1 = 125


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
	real(8) irate(nk), wage(nk), income, employ
    integer goodtime, badtime

	!	grid vectors
	real(8) agrid(na), epsgrid(neps), kgrid(nk), zgrid(nz), btagrid(nbta), 	PGgrid (nP), PBgrid (nP)
	

	!	transition probability
	real(8) trEps(nz,nz,neps,neps), CtrEps(nz,nz,neps,neps), &
	        trZ(nz,nz), CtrZ(nz,nz), trBta(nbta,nbta), CtrBta(nbta,nbta)


	!	simulated time series data
	real(8) Kdata(Nperiod+1), Ldata(Nperiod), Ydata(Nperiod), K(Nperiod+1)
	real(8) Idata(Nperiod), Cdata(Nperiod), Adata(Nperiod), Edata(Nperiod), Rdata(Nperiod), Wdata(Nperiod)
	integer iZdata(Nperiod+1), Zdata(Nperiod+1)
	real(8) errKdata(Nperiod+1)


	!	coefficients of equations
	real(8)	Kcoef(2), NewKcoef(2)


	!	indexing variables
	 !**** iz=1 is bad times and iz=2 is good times ****	
	integer ia, iz, ik, ieps, ibta, indiv, time
	logical final
	integer icut, iinterv
	integer iPg, iPb


	!	spaces for value functions and decision rules
	
	real(8) V(na,nP, nP, nbta,nk), EV(na,nP,nP,nbta,nk), AS(na,nP,nP,nbta,nk)
	real(8)  EV2p(na,nP,nP,nbta,nk),EV2pp(na,nP,nP,nbta)


!   cutoff values
    real(8)  Cutoff(nz,nz,neps), cutvecttemp(ncut), cutvect(ncut)
	integer  zindtemp(ncut), zpindtemp(ncut), epsindtemp(ncut)
	integer  zind(ncut), zpind(ncut), epsind(ncut), indcut(ncut)
	
!   P computation variables
	real(8)  probint(ninterv), Pprime(ninterv,nz), Pvalue(nz), gbar(ninterv,nz,nz,neps), evalpoint(ninterv)
    real(8) mubar, Pprimeb(ninterv), Pprimeg(ninterv)


    real(8) zseries(N1), useries(N1)
	real(8) irate2(N1), wage2(N1)
	real(8) AS2(N1,na,nP,nP,nbta), V2(N1,na,nP,nP,nbta), EV2(na,nP,nP,nbta)
	real(8) pi(2,N1+1)
	integer Backtime

!   Panel data
     real(8) PanelA(Nindiv), PanelPg(Nindiv),PanelPb(Nindiv)
	 integer PanelBta(Nindiv)

!  cross section

	real(8) Adist(Nindiv), Gini(Nperiod), wealth1(Nperiod), wealth5(Nperiod), wealth10(Nperiod)
	real(8) wealth20(Nperiod), wealth30(Nperiod), Negratio(Nperiod)

	real(8) V1(na,neps,nbta,nk,nz),AS1(na,neps,nbta,nk,nz)

end module
