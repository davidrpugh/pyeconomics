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

	real(8) kp, VP(na,nP,nP,nP,nP,nP,nbta), VPP(na,nbta)
	integer ibtap, iPgep, iPbep, iPgsp, iPbsp, iPbfp, i




	EVl = 0.0D0
	EVm = 0.0D0
    EVh = 0.0D0

	do ik = 1, nk


    	kp = dexp(Kcoef(1) + Kcoef(2)*dlog(kgrid(ik)))
	        
			do ia = 1, na
            do iPgep = 1, nP
			do iPbep = 1, nP
            do iPgsp = 1, nP
			do iPbsp = 1, nP
			do iPbfp = 1, nP

            
				VP(ia,iPgep,iPbep,iPgsp,iPbsp,iPbfp,1) = polyinterp1(kp, kgrid, Vl(ia,iPgep,iPbep,iPgsp,iPbsp,iPbfp,:))
				VP(ia,iPgep,iPbep,iPgsp,iPbsp,iPbfp,2) = polyinterp1(kp, kgrid, Vm(ia,iPgep,iPbep,iPgsp,iPbsp,iPbfp,:))				
				VP(ia,iPgep,iPbep,iPgsp,iPbsp,iPbfp,3) = polyinterp1(kp, kgrid, Vh(ia,iPgep,iPbep,iPgsp,iPbsp,iPbfp,:))
				
			end do
			end do
			end do
			end do
			end do
            end do
			
		!	calculate expected value function	
		
		do iPge = 1, nP
		do iPbe = 1, nP
		do iPgs = 1, nP
		do iPbs = 1, nP
		do iPbf = 1, nP
         
		 ! Compute next period's Ps
				!first give today's P values.  THESE ARE UNCONDITIONAL VALUES!!!!!!
		 Pvalue(1,4) = 0.5D0*Pgrid(iPbe)  !bad, e
		 Pvalue(1,3) = 0.5D0*(1.0D0-Pgrid(iPbe))*Pgrid(iPbs)  !bad, s
		 Pvalue(1,2) = 0.5D0*(1.0D0-Pgrid(iPbe))*(1.0D0-Pgrid(iPbs))*Pgrid(iPbf)  !bad, f
		 Pvalue(1,1) = 0.5D0-Pvalue(1,4)-Pvalue(1,3)-Pvalue(1,2) !bad, l
		 Pvalue(2,4) = 0.5D0*Pgrid(iPge)  !good, e
		 Pvalue(2,3) = 0.5D0*(1.0D0-Pgrid(iPge))*Pgrid(iPgs)  !good, s
		 Pvalue(2,2) = 0.0D0  !good, f
		 Pvalue(2,1) = 0.5D0-Pvalue(2,4)-Pvalue(2,3) !good, l

		 call ComputePprimes()
		
		! I treat the first and last interval differently. These intervals are defined as (0,cutvect(1))
		! and (cutvect(ninterv-1),1). There are eight cutoff points but nine different values. 

         do ia = 1, na
	  
		      
			  do i = 1,  ninterv
			   if (i==1) then  
			   	do ibtap = 1, nbta
			    VPP(ia,ibtap) = lininterp5(Pprimege(i),Pprimebe(i),Pprimegs(i),Pprimebs(i),Pprimebf(i),Pgrid,Pgrid,Pgrid,Pgrid,Pgrid,VP(ia,:,:,:,:,:,ibtap))           
  	            EVl(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik) =  EVl(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik)+(cutvect(1))*VPP(ia,ibtap)*trBta(1,ibtap)
				EVm(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik) =  EVm(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik)+(cutvect(1))*VPP(ia,ibtap)*trBta(2,ibtap)
				EVh(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik) =  EVh(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik)+(cutvect(1))*VPP(ia,ibtap)*trBta(3,ibtap)


				end do
			   else if (i==ninterv) then
			   do ibtap = 1, nbta
			    VPP(ia,ibtap) = lininterp5(Pprimege(i),Pprimebe(i),Pprimegs(i),Pprimebs(i),Pprimebf(i),Pgrid,Pgrid,Pgrid,Pgrid,Pgrid,VP(ia,:,:,:,:,:,ibtap))           
  	            EVl(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik) =  EVl(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik)+(1.0D0-cutvect(ninterv-1))*VPP(ia,ibtap)*trBta(1,ibtap)
				EVm(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik) =  EVm(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik)+(1.0D0-cutvect(ninterv-1))*VPP(ia,ibtap)*trBta(2,ibtap)
				EVh(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik) =  EVh(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik)+(1.0D0-cutvect(ninterv-1))*VPP(ia,ibtap)*trBta(3,ibtap)
				end do
			   else 
			    do ibtap = 1, nbta
			    VPP(ia,ibtap) = lininterp5(Pprimege(i),Pprimebe(i),Pprimegs(i),Pprimebs(i),Pprimebf(i),Pgrid,Pgrid,Pgrid,Pgrid,Pgrid,VP(ia,:,:,:,:,:,ibtap))           
  	            EVl(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik) =  EVl(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik)+(cutvect(i)-cutvect(i-1))*VPP(ia,ibtap)*trBta(1,ibtap)
				EVm(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik) =  EVm(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik)+(cutvect(i)-cutvect(i-1))*VPP(ia,ibtap)*trBta(2,ibtap)
				EVh(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik) =  EVh(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik)+(cutvect(i)-cutvect(i-1))*VPP(ia,ibtap)*trBta(3,ibtap)
				end do
			  end if 
             
			 end do
		   

		
        end do
        
		
		end do
		end do
		end do
		end do
		end do

end do	



	call setsplineV()
  
end subroutine
