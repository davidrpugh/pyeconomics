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
	real(8), parameter:: hour=0.3271, cutvalue=-1.0


!	number of grids
	integer, parameter:: na = 50, neps = 4, nk = 4, nz = 2, nbta = 3
    integer, parameter:: ncut=nz*nz*neps*(neps-1), nP=5,  ninterv=ncut+1


	!	maximum number of iterations

	integer, parameter:: maxiterLOM = 60, showerr = 2, maxiterK = 200
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

	real(8), parameter:: alpha = 0.36d0, delta = 0.025d0, uih=0.391, uil=0.038
	real(8) bta(nbta)
    
	! **** Unemployment rate for good and bad times ****
    real(8)  mu(nz)
    
	 
	!	auxiliary variables
	real(8) irate(nk), wage(nk), income, employ
    integer goodtime, badtime

	!	grid vectors
	real(8) agrid(na), epsgrid(neps), kgrid(nk), zgrid(nz), btagrid(nbta), 	Pgrid(nP)
	

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
	integer ia, iz, ik, ieps, ibta, indiv, time, iepsp
	logical final
	integer icut, iinterv
	integer iPge, iPbe, iPgs, iPbs, iPbf


	!	spaces for value functions and decision rules
	
	real(8) Vl(na,nP,nP,nP,nP,nP,nk), EVl(na,nP,nP,nP,nP,nP,nk), ASlo(na,nP,nP,nP,nP,nP,nk)
	real(8) Vm(na,nP,nP,nP,nP,nP,nk), EVm(na,nP,nP,nP,nP,nP,nk), ASmi(na,nP,nP,nP,nP,nP,nk)
	real(8) Vh(na,nP,nP,nP,nP,nP,nk), EVh(na,nP,nP,nP,nP,nP,nk), AShi(na,nP,nP,nP,nP,nP,nk)

	real(8)  EV2pl(na,nP,nP,nP,nP,nP,nk)
	real(8)  EV2pm(na,nP,nP,nP,nP,nP,nk)
	real(8)  EV2ph(na,nP,nP,nP,nP,nP,nk)


!   cutoff values
    real(8)  Cutoff(nz,nz,neps,neps), cutvecttemp(ncut), cutvect(ncut)
	integer  zindtemp(ncut), zpindtemp(ncut), epsindtemp(ncut), epspindtemp(ncut)
	integer  zind(ncut), zpind(ncut), epsind(ncut), indcut(ncut), epspind(ncut)
	
!   P computation variables
	real(8)  probint(ninterv), Pvalue(nz,neps), gbar(ninterv,nz,nz,neps,neps), evalpoint(ninterv)
    real(8)  mubar, Pprimege(ninterv), Pprimebe(ninterv), Pprimegs(ninterv), Pprimebs(ninterv), Pprimebf(ninterv)


    real(8) zseries(N1), useries(N1)
	real(8) irate2(N1), wage2(N1)
	real(8) AS2l(N1,na,nP,nP,nP,nP,nP), V2l(N1,na,nP,nP,nP,nP,nP)
	real(8) AS2m(N1,na,nP,nP,nP,nP,nP), V2m(N1,na,nP,nP,nP,nP,nP)
	real(8) AS2h(N1,na,nP,nP,nP,nP,nP), V2h(N1,na,nP,nP,nP,nP,nP)
    real(8) EV2(na,nP,nP,nP,nP,nP,nbta), EV2pp(na,nP,nP,nP,nP,nP,nbta)

	real(8) pi(2,N1+1)
	integer Backtime

!   Panel data
     real(8) PanelA(Nindiv), PanelPge(Nindiv), PanelPbe(Nindiv), PanelPgs(Nindiv), PanelPbs(Nindiv), PanelPbf(Nindiv)
	 integer PanelBta(Nindiv)

!  cross section

	real(8) Adist(Nindiv), Gini(Nperiod), wealth1(Nperiod), wealth5(Nperiod), wealth10(Nperiod)
	real(8) wealth20(Nperiod), wealth30(Nperiod), Negratio(Nperiod)

	real(8) V1(na,neps,nbta,nk,nz+1),AS1(na,neps,nbta,nk,nz+1)

end module
