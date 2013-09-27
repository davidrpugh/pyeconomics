!	----------------------------------------------------------------------
!	File name: SaveValueFunction.f90
!	----------------------------------------------------------------------
	
	

subroutine SaveCoefficients()

	use Globals

	open(1, file='Output\CoefficientsOUTG.txt', status='unknown')
	write(1, '(2F10.6)') KcoefG(1), KcoefG(2)
	open(2, file='Output\CoefficientsOUTB.txt', status='unknown')
	write(2, '(2F10.6)') KcoefB(1), KcoefB(2)



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

	open(13, file='Output\PanelA112G.txt', status='unknown')
	open(14, file='Output\PanelEps112G.txt', status='unknown')
    open(15, file='Output\PanelBta112G.txt', status='unknown')
	open(16, file='Output\PanelA112B.txt', status='unknown')
	open(17, file='Output\PanelEps112B.txt', status='unknown')
    open(18, file='Output\PanelBta112B.txt', status='unknown')
	open(19, file='Output\PanelA123G.txt', status='unknown')
	open(20, file='Output\PanelEps123G.txt', status='unknown')
    open(21, file='Output\PanelBta123G.txt', status='unknown')
	open(22, file='Output\PanelA123B.txt', status='unknown')
	open(23, file='Output\PanelEps123B.txt', status='unknown')
    open(24, file='Output\PanelBta123B.txt', status='unknown')


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

	write(13, '(F14.8)') PanelA112G
	write(14, '(I1)') PanelEps112G
	write(15, '(I1)') PanelBta112G
	write(16, '(F14.8)') PanelA112B
	write(17, '(I1)') PanelEps112B
	write(18, '(I1)') PanelBta112B
	write(19, '(F14.8)') PanelA123G
	write(20, '(I1)') PanelEps123G
	write(21, '(I1)') PanelBta123G
	write(22, '(F14.8)') PanelA123B
	write(23, '(I1)') PanelEps123B
	write(24, '(I1)') PanelBta123B

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

