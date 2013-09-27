!	----------------------------------------------------------------------
!	File name: SaveValueFunction.f90
!	----------------------------------------------------------------------
	
	

subroutine SaveCoefficients()

	use Globals

	open(1, file='Output\CoefficientsOUTG.txt', status='unknown')
	write(1, '(2F10.6)') KcoefG(1), KcoefG(2)
	open(2, file='Output\CoefficientsOUTB.txt', status='unknown')
	write(2, '(2F10.6)') KcoefB(1), KcoefB(2)
	open(3, file='Output\CoefficientsOUTD.txt', status='unknown')
	write(3, '(2F10.6)') KcoefD(1), KcoefD(2)



end subroutine



subroutine SaveSimulatedTimeSeriesData()

	use Globals

	open(1, file='Output\Kdata.txt', status='unknown')
	open(2, file='Output\Rdata.txt', status='unknown')
	open(3, file='Output\Wdata.txt', status='unknown')
	open(4, file='Output\Ldata.txt', status='unknown')
	open(5, file='Output\Adata.txt', status='unknown')
	open(6, file='Output\Gini.txt', status='unknown')
	open(7, file='Output\wealth1.txt', status='unknown')
	open(8, file='Output\wealth5.txt', status='unknown')
	open(9, file='Output\wealth10.txt', status='unknown')
	open(10, file='Output\wealth20.txt', status='unknown')  
	open(11, file='Output\wealth30.txt', status='unknown')
	open(12, file='Output\Negratio.txt', status='unknown')  

	open(31, file='Output\AvgDur.txt', status='unknown')

	open(13, file='Output\PanelA113B.txt', status='unknown')
	open(14, file='Output\PanelEps113B.txt', status='unknown')
    open(15, file='Output\PanelBta113B.txt', status='unknown')
	open(16, file='Output\PanelA113D.txt', status='unknown')
	open(17, file='Output\PanelEps113D.txt', status='unknown')
    open(18, file='Output\PanelBta113D.txt', status='unknown')
	open(19, file='Output\PanelA113G.txt', status='unknown')
	open(20, file='Output\PanelEps113G.txt', status='unknown')
    open(21, file='Output\PanelBta113G.txt', status='unknown')
	open(22, file='Output\PanelA121B.txt', status='unknown')
	open(23, file='Output\PanelEps121B.txt', status='unknown')
    open(24, file='Output\PanelBta121B.txt', status='unknown')
	open(25, file='Output\PanelA121D.txt', status='unknown')
	open(26, file='Output\PanelEps121D.txt', status='unknown')
    open(27, file='Output\PanelBta121D.txt', status='unknown')
	open(28, file='Output\PanelA121G.txt', status='unknown')
	open(29, file='Output\PanelEps121G.txt', status='unknown')
    open(30, file='Output\PanelBta121G.txt', status='unknown')

    open(32, file='Output\PanelDur113B.txt', status='unknown')
    open(33, file='Output\PanelDur113D.txt', status='unknown')
    open(34, file='Output\PanelDur113G.txt', status='unknown')
    open(35, file='Output\PanelDur121B.txt', status='unknown')
    open(36, file='Output\PanelDur121D.txt', status='unknown')
    open(37, file='Output\PanelDur121G.txt', status='unknown')

    open(38, file='Output\AvgDurNonzero.txt', status='unknown')
    open(39, file='Output\ratioover2.txt', status='unknown')

	open(40, file='Output\Empdata.txt', status='unknown')
	open(41, file='Output\Shodata.txt', status='unknown')
	open(42, file='Output\Firdata.txt', status='unknown')
	open(43, file='Output\Londata.txt', status='unknown')



	write(1, '(F12.6)') Kdata
	write(2, '(F12.6)') Rdata
	write(3, '(F12.6)') Wdata
	write(4, '(F12.6)') Ldata
	write(5, '(F12.6)') Adist
	write(6, '(F12.6)') Gini
	write(7, '(F12.6)') wealth1
	write(8, '(F12.6)') wealth5
	write(9, '(F12.6)') wealth10
	write(10, '(F12.6)') wealth20
	write(11, '(F12.6)') wealth30
	write(12, '(F12.6)') Negratio

	write(13, '(F14.8)') PanelA113B
	write(14, '(I1)') PanelEps113B
	write(15, '(I1)') PanelBta113B
	write(16, '(F14.8)') PanelA113D
	write(17, '(I1)') PanelEps113D
	write(18, '(I1)') PanelBta113D
	write(19, '(F14.8)') PanelA113G
	write(20, '(I1)') PanelEps113G
	write(21, '(I1)') PanelBta113G
	write(22, '(F14.8)') PanelA121B
	write(23, '(I1)') PanelEps121B
	write(24, '(I1)') PanelBta121B
	write(25, '(F14.8)') PanelA121D
	write(26, '(I1)') PanelEps121D
	write(27, '(I1)') PanelBta121D
	write(28, '(F14.8)') PanelA121G
	write(29, '(I1)') PanelEps121G
	write(30, '(I1)') PanelBta121G

	write(31, '(F12.6)') AvgDur
	write(38, '(F12.6)') AvgDurNonzero
	write(39, '(F12.6)') ratioover2

	write(32, '(I6)') PanelDur113B
	write(33, '(I6)') PanelDur113D
	write(34, '(I6)') PanelDur113G
	write(35, '(I6)') PanelDur121B
	write(36, '(I6)') PanelDur121D
	write(37, '(I6)') PanelDur121G


	write(40, '(F12.6)') Empdata
	write(41, '(F12.6)') Shodata
	write(42, '(F12.6)') Firdata
	write(43, '(F12.6)') Londata

end subroutine



subroutine SaveValueFunction()

	use Globals


	open(1, file='Output\V.txt', status='unknown')
	open(2, file='Output\AS.txt', status='unknown')

    do ia = 1, na
	do ieps = 1, neps		
    do ibta = 1, nbta
	do ik = 1, nk
	do iz = 1, nz
	

		 write(1, '(F15.6)') V(ia,ieps,ibta,ik,iz)
		 write(2, '(F15.6)') AS(ia,ieps,ibta,ik,iz)
		
    end do
	end do
	end do
	end do
	end do 


end subroutine

