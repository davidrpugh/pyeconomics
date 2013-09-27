!	----------------------------------------------------------------------
!	File name: IterateBack.f90
!	
!	----------------------------------------------------------------------


subroutine IterateBack()

	use Globals
	use PolyInterpModule
	use LinInterpModule

	implicit none

	real(8) kp, VP(na,nP,nP,nbta), VPP(na,nbta)
	integer  ibtap, i, iPgp, iPbp




	kp = K(N1+1)
	EV2 = 0.0D0
	
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
		
         do ia = 1, na
		 do ibta = 1, nbta   
		      
			  do i = 1,  ninterv
			   if (i==1) then  
			   	do ibtap = 1, nbta
			    VPP(ia,ibtap) = lininterp2(Pprimeg(i),Pprimeb(i), PGgrid, PBgrid, VP(ia,:,:,ibtap))           
  	            EV2(ia,iPg,iPb,ibta) =  EV2(ia,iPg,iPb,ibta)+(cutvect(1))*VPP(ia,ibtap)*trBta(ibta,ibtap)
				end do
			   else if (i==ninterv) then
			   do ibtap = 1, nbta
			    VPP(ia,ibtap) = lininterp2(Pprimeg(i),Pprimeb(i), PGgrid, PBgrid, VP(ia,:,:,ibtap))    		   
  	            EV2(ia,iPg,iPb,ibta) =  EV2(ia,iPg,iPb,ibta)+(1.0-cutvect(ninterv-1))*VPP(ia,ibtap)*trBta(ibta,ibtap)
				end do
			   else 
			    do ibtap = 1, nbta
			    VPP(ia,ibtap) = lininterp2(Pprimeg(i),Pprimeb(i), PGgrid, PBgrid, VP(ia,:,:,ibtap))    
  	            EV2(ia,iPg,iPb,ibta) =  EV2(ia,iPg,iPb,ibta)+(cutvect(i)-cutvect(i-1))*VPP(ia,ibtap)*trBta(ibta,ibtap)
				end do
			   end if      
			  end do
		   
       
		end do
        end do
        end do
		end do
	

   call setsplineV2()

    do ia = 1, na
	do iPg = 1, nP
	do iPb = 1, nP
	  
    	 income = (1.0+irate2(N1))*agrid(ia) + wage2(N1)*(0.5*PGgrid(iPg)+0.5*PBgrid(iPb))*hour+ui*(1.0-(0.5*PGgrid(iPg))-(0.5*PBgrid(iPb)))
	     do ibta = 1,nbta
         
		  time = N1
	      call BrentU2(AS2(N1,ia,iPg,iPb,ibta), V2(N1,ia,iPg,iPb,ibta))

	     end do
	end do
	end do
	end do
	

	do Backtime = 1, N1-1
     time = N1-Backtime
	 kp = K(N1+1-Backtime)
     EV2 = 0.0D0

	!	calculate expected value function	
		
		do iPg = 1, nP
		do iPb = 1, nP
         
		 ! Compute next period's Pg and Pb 
		 Pvalue(1) = PBgrid(iPb)
		 Pvalue(2) = PGgrid(iPg)
		 call ComputePprimes2()
		
         do ia = 1, na
		 do ibta = 1, nbta   
		      
			  do i = 1,  ninterv
			   if (i==1) then  
			   	do ibtap = 1, nbta
				
			    VPP(ia,ibtap) = lininterp2(Pprimeg(i),Pprimeb(i), PGgrid, PBgrid, V2(N1+1-Backtime,ia,:,:,ibtap))           
  	            EV2(ia,iPg,iPb,ibta) =  EV2(ia,iPg,iPb,ibta)+(cutvect(1))*VPP(ia,ibtap)*trBta(ibta,ibtap)
				end do
			   else if (i==ninterv) then
			   do ibtap = 1, nbta
			    VPP(ia,ibtap) = lininterp2(Pprimeg(i),Pprimeb(i), PGgrid, PBgrid, V2(N1+1-Backtime,ia,:,:,ibtap))    		   
  	            EV2(ia,iPg,iPb,ibta) =  EV2(ia,iPg,iPb,ibta)+(1.0-cutvect(ninterv-1))*VPP(ia,ibtap)*trBta(ibta,ibtap)
				end do
			   else 
			    do ibtap = 1, nbta
			    VPP(ia,ibtap) = lininterp2(Pprimeg(i),Pprimeb(i), PGgrid, PBgrid, V2(N1+1-Backtime,ia,:,:,ibtap))    
  	            EV2(ia,iPg,iPb,ibta) =  EV2(ia,iPg,iPb,ibta)+(cutvect(i)-cutvect(i-1))*VPP(ia,ibtap)*trBta(ibta,ibtap)
				end do
			  end if 
             
			 end do
		   

		 end do
        end do
        end do
		end do

    call setsplineV2()

    do ia = 1, na
	do iPg = 1, nP
	do iPb = 1, nP
	  
    	 income = (1.0+irate2(N1-Backtime))*agrid(ia)+wage2(N1-Backtime)*(pi(2,time)*PGgrid(iPg)+pi(1,time)*PBgrid(iPb))*hour+ui*(1.0-(pi(2,time)*PGgrid(iPg))-(pi(1,time)*PBgrid(iPb)))
	     do ibta = 1,nbta
         
	       call BrentU2(AS2(N1-Backtime,ia,iPg,iPb,ibta), V2(N1-Backtime,ia,iPg,iPb,ibta))
           
	     end do
	end do
	end do
	end do
	
    end do

	!call SaveValueFunction()

 
end subroutine

	    






