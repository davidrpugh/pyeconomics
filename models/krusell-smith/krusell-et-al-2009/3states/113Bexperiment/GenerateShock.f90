!	----------------------------------------------------------------------
!	File name: GenerateShock.f90
!
!	Calculates expected value function given current state variables.
!	----------------------------------------------------------------------



subroutine NextBeta(CrBta)

	use Numerical_Libraries
	use Globals

	integer j, CrBta(Nindiv,2)
	real(8) ushock


	do indiv = 1, Nindiv

		ushock = drnunf()
        
		do j = 1, nbta

			if (ushock <= CtrBta(CrBta(indiv,1),j)) then
				CrBta(indiv,2) = j
				exit
			end if
		
		end do

	end do
end subroutine

subroutine NextP(CrPge,CrPbe,CrPgs,CrPbs,CrPbf)

	use Numerical_Libraries
	use Globals

	integer j
	real(8) ushock, CrPge(Nindiv,2),CrPbe(Nindiv,2),CrPgs(Nindiv,2),CrPbs(Nindiv,2),CrPbf(Nindiv,2)


	do indiv = 1, Nindiv

		ushock = drnunf()

		 Pvalue(1,4) = 0.5D0*CrPbe(indiv,1)  !bad, e
		 Pvalue(1,3) = 0.5D0*(1.0D0-CrPbe(indiv,1))*CrPbs(indiv,1)  !bad, s
		 Pvalue(1,2) = 0.5D0*(1.0D0-CrPbe(indiv,1))*(1.0D0-CrPbs(indiv,1))*CrPbf(indiv,1)  !bad, f
		 Pvalue(1,1) = 0.5D0-Pvalue(1,4)-Pvalue(1,3)-Pvalue(1,2) !bad, l
		 Pvalue(2,4) = 0.5D0*CrPge(indiv,1)  !good, e
		 Pvalue(2,3) = 0.5D0*(1.0D0-CrPge(indiv,1))*CrPgs(indiv,1)  !good, s
		 Pvalue(2,2) = 0.0D0  !good, f
		 Pvalue(2,1) = 0.5D0-Pvalue(2,4)-Pvalue(2,3) !good, l

		call ComputePprimes()
        CrPge(indiv,2) = Pprimege(ninterv)
		CrPbe(indiv,2) = Pprimebe(ninterv)
        CrPgs(indiv,2) = Pprimegs(ninterv)
		CrPbs(indiv,2) = Pprimebs(ninterv)
		CrPbf(indiv,2) = Pprimebf(ninterv)
        
		do j = 1, ncut      	
			if (ushock <= cutvect(j)) then
				CrPge(indiv,2) = Pprimege(j)
				CrPbe(indiv,2) = Pprimebe(j)
				CrPgs(indiv,2) = Pprimegs(j)
				CrPbs(indiv,2) = Pprimebs(j)
				CrPbf(indiv,2) = Pprimebf(j)


			exit
			end if		
		end do
	end do
end subroutine


subroutine NextP2(CrPge,CrPbe,CrPgs,CrPbs,CrPbf)

	use Numerical_Libraries
	use Globals

	integer j
	real(8) ushock, CrPge(Nindiv,2),CrPbe(Nindiv,2),CrPgs(Nindiv,2),CrPbs(Nindiv,2),CrPbf(Nindiv,2)


	do indiv = 1, Nindiv

		ushock = drnunf()

		 Pvalue(1,4) = pi(1,time)*CrPbe(indiv,1)  !bad, e
		 Pvalue(1,3) = pi(1,time)*(1.0D0-CrPbe(indiv,1))*CrPbs(indiv,1)  !bad, s
		 Pvalue(1,2) = pi(1,time)*(1.0D0-CrPbe(indiv,1))*(1.0D0-CrPbs(indiv,1))*CrPbf(indiv,1)  !bad, f
		 Pvalue(1,1) = pi(1,time)-Pvalue(1,4)-Pvalue(1,3)-Pvalue(1,2) !bad, l
		 Pvalue(2,4) = pi(2,time)*CrPge(indiv,1)  !good, e
		 Pvalue(2,3) = pi(2,time)*(1.0D0-CrPge(indiv,1))*CrPgs(indiv,1)  !good, s
		 Pvalue(2,2) = 0.0D0  !good, f
		 Pvalue(2,1) = pi(2,time)-Pvalue(2,4)-Pvalue(2,3) !good, l

		call ComputePprimes2()
        CrPge(indiv,2) = Pprimege(ninterv)
		CrPbe(indiv,2) = Pprimebe(ninterv)
        CrPgs(indiv,2) = Pprimegs(ninterv)
		CrPbs(indiv,2) = Pprimebs(ninterv)
		CrPbf(indiv,2) = Pprimebf(ninterv)
        
		do j = 1, ncut      	
			if (ushock <= cutvect(j)) then
				CrPge(indiv,2) = Pprimege(j)
				CrPbe(indiv,2) = Pprimebe(j)
				CrPgs(indiv,2) = Pprimegs(j)
				CrPbs(indiv,2) = Pprimebs(j)
				CrPbf(indiv,2) = Pprimebf(j)


			exit
			end if		
		end do
	end do

end subroutine

