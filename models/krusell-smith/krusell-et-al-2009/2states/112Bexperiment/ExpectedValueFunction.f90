!	----------------------------------------------------------------------
!	File name: ExpectedValueFunction.f90
!
!	Calculates expected value function given current state variables.
!	----------------------------------------------------------------------



subroutine ExpectedValueFunction()

	use Globals
	use PolyInterpModule
	use LinInterpModule

	implicit none

	real(8) kp, VP(na,nP,nP,nbta), VPP(na,nbta)
	integer ibtap, iPgp,iPbp, i



	EV = 0D0

	do ik = 1, nk

    	kp = dexp(Kcoef(1) + Kcoef(2)*dlog(kgrid(ik)))
	        
			do ia = 1, na
            do iPgp = 1, nP
			do iPbp = 1, nP
            do ibtap = 1, nbta

				VP(ia,iPgp,iPbp,ibtap) = polyinterp1(kp, kgrid, V(ia,iPgp,iPbp,ibtap,:))
            
			end do
			end do
			end do
            end do
			
		!	calculate expected value function	
		
		do iPg = 1, nP
		do iPb = 1, nP
         
		 ! Compute next period's Pg and Pb 
	

		 Pvalue(1) = PBgrid(iPb)
		 Pvalue(2) = PGgrid(iPg)
		 call ComputePprimes()
		
		! I treat the first and last interval differently. These intervals are defined as (0,cutvect(1))
		! and (cutvect(ninterv-1),1). There are eight cutoff points but nine different values. 
	

         do ia = 1, na
		 do ibta = 1, nbta   
		      
			  do i = 1,  ninterv
			   if (i==1) then  
			   	do ibtap = 1, nbta
			    VPP(ia,ibtap) = lininterp2(Pprimeg(i),Pprimeb(i), PGgrid, PBgrid, VP(ia,:,:,ibtap))           
  	            EV(ia,iPg,iPb,ibta,ik) =  EV(ia,iPg,iPb,ibta,ik)+(cutvect(1))*VPP(ia,ibtap)*trBta(ibta,ibtap)
				end do
			   else if (i==ninterv) then
			   do ibtap = 1, nbta
			    VPP(ia,ibtap) = lininterp2(Pprimeg(i),Pprimeb(i), PGgrid, PBgrid, VP(ia,:,:,ibtap))    		   
  	            EV(ia,iPg,iPb,ibta,ik) =  EV(ia,iPg,iPb,ibta,ik)+(1.0-cutvect(ninterv-1))*VPP(ia,ibtap)*trBta(ibta,ibtap)
				end do
			   else 
			    do ibtap = 1, nbta
			    VPP(ia,ibtap) = lininterp2(Pprimeg(i),Pprimeb(i), PGgrid, PBgrid, VP(ia,:,:,ibtap))    
  	            EV(ia,iPg,iPb,ibta,ik) =  EV(ia,iPg,iPb,ibta,ik)+(cutvect(i)-cutvect(i-1))*VPP(ia,ibtap)*trBta(ibta,ibtap)
				end do
			  end if 
             
			 end do
		   

		 end do
        end do
        end do
		end do
	end do	

	call setsplineV()
  
end subroutine
