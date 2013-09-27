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

subroutine NextP(CrPg,CrPb)

	use Numerical_Libraries
	use Globals

	integer j
	real(8) ushock, CrPg(Nindiv,2),CrPb(Nindiv,2)


	do indiv = 1, Nindiv

		ushock = drnunf()
        Pvalue(1) = CrPb(indiv,1)
		Pvalue(2) = CrPg(indiv,1)
		call ComputePprimes()
        CrPg(indiv,2) = Pprimeg(ninterv)
		CrPb(indiv,2) = Pprimeb(ninterv)

		do j = 1, ncut
			if (ushock <= cutvect(j)) then
				CrPg(indiv,2) = Pprimeg(j)
				CrPb(indiv,2) = Pprimeb(j)
			exit
			end if
		
		end do

	end do

end subroutine


subroutine NextP2(CrPg,CrPb)

	use Numerical_Libraries
	use Globals

	integer j
	real(8) ushock, CrPg(Nindiv,2),CrPb(Nindiv,2)


	do indiv = 1, Nindiv

		ushock = drnunf()
        Pvalue(1) = CrPb(indiv,1)
		Pvalue(2) = CrPg(indiv,1)
		call ComputePprimes2()
        CrPg(indiv,2) = Pprimeg(ninterv)
		CrPb(indiv,2) = Pprimeb(ninterv)
        
		do j = 1, ncut      	
			if (ushock <= cutvect(j)) then
				CrPg(indiv,2) = Pprimeg(j)
				CrPb(indiv,2) = Pprimeb(j)
			exit
			end if		
		end do
	end do

end subroutine

