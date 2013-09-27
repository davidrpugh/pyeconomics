!	----------------------------------------------------------------------
!	File name: IterateBack.f90
!	
!	----------------------------------------------------------------------


subroutine IterateBack()

	use Globals
	use PolyInterpModule
	use LinInterpModule

	implicit none

	real(8) kp, VP(na,nP,nP,nP,nP,nP,nbta), VPP(na,nbta)
	integer  ibtap, i, iPgep, iPbep, iPgsp, iPbsp, iPbfp




	kp = K(N1+1)
	EV2 = 0.0D0
	
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
		 do ibta = 1, nbta   
		      
			  do i = 1,  ninterv
			   if (i==1) then  
			   	do ibtap = 1, nbta
			    VPP(ia,ibtap) = lininterp5(Pprimege(i),Pprimebe(i),Pprimegs(i),Pprimebs(i),Pprimebf(i),Pgrid,Pgrid,Pgrid,Pgrid,Pgrid,VP(ia,:,:,:,:,:,ibtap))           
  	            EV2(ia,iPge,iPbe,iPgs,iPbs,iPbf,ibta) =  EV2(ia,iPge,iPbe,iPgs,iPbs,iPbf,ibta)+(cutvect(1))*VPP(ia,ibtap)*trBta(ibta,ibtap)
	     	    end do
			   else if (i==ninterv) then
			   do ibtap = 1, nbta
			    VPP(ia,ibtap) = lininterp5(Pprimege(i),Pprimebe(i),Pprimegs(i),Pprimebs(i),Pprimebf(i),Pgrid,Pgrid,Pgrid,Pgrid,Pgrid,VP(ia,:,:,:,:,:,ibtap))           
  	            EV2(ia,iPge,iPbe,iPgs,iPbs,iPbf,ibta) =  EV2(ia,iPge,iPbe,iPgs,iPbs,iPbf,ibta)+(1.0D0-cutvect(ninterv-1))*VPP(ia,ibtap)*trBta(ibta,ibtap)
			   end do
			   else 
			    do ibtap = 1, nbta
			    VPP(ia,ibtap) = lininterp5(Pprimege(i),Pprimebe(i),Pprimegs(i),Pprimebs(i),Pprimebf(i),Pgrid,Pgrid,Pgrid,Pgrid,Pgrid,VP(ia,:,:,:,:,:,ibtap))           
  	            EV2(ia,iPge,iPbe,iPgs,iPbs,iPbf,ibta) =  EV2(ia,iPge,iPbe,iPgs,iPbs,iPbf,ibta)+(cutvect(i)-cutvect(i-1))*VPP(ia,ibtap)*trBta(ibta,ibtap)
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

   call setsplineV2()

    do ia = 1, na
	do iPge = 1, nP
    do iPbe = 1, nP
	do iPgs = 1, nP
	do iPbs = 1, nP
	do iPbf = 1, nP

		income = (1.0+irate2(N1))*agrid(ia) + wage2(N1)*(0.5D0*Pgrid(iPge)+0.5D0*Pgrid(iPbe))*hour &
				         +uih*(0.5D0*(1.0D0-Pgrid(iPge))*Pgrid(iPgs)+0.5D0*(1.0D0-Pgrid(iPbe))*Pgrid(iPbs)+0.5D0*(1.0D0-Pgrid(iPbe))*(1.0D0-Pgrid(iPbs))*Pgrid(iPbf)) &
						 +uil*(1.0-0.5D0*Pgrid(iPge)-0.5D0*Pgrid(iPbe)-0.5D0*(1.0D0-Pgrid(iPge))*Pgrid(iPgs)-0.5D0*(1.0D0-Pgrid(iPbe))*Pgrid(iPbs)-0.5D0*(1.0D0-Pgrid(iPbe))*(1.0D0-Pgrid(iPbs))*Pgrid(iPbf))


	     do ibta = 1,nbta
         
		  time = N1

					if (ibta==1) then

						if (agrid(ia)>cutvalue) then
						call NewtonU2(AS2l(N1,ia,iPge,iPbe,iPgs,iPbs,iPbf), V2l(N1,ia,iPge,iPbe,iPgs,iPbs,iPbf))
						else
						call BisectU2(AS2l(N1,ia,iPge,iPbe,iPgs,iPbs,iPbf), V2l(N1,ia,iPge,iPbe,iPgs,iPbs,iPbf))
						end if

					else if (ibta==2) then

						if (agrid(ia)>cutvalue) then
						call NewtonU2(AS2m(N1,ia,iPge,iPbe,iPgs,iPbs,iPbf), V2m(N1,ia,iPge,iPbe,iPgs,iPbs,iPbf))
						else
						call BisectU2(AS2m(N1,ia,iPge,iPbe,iPgs,iPbs,iPbf), V2m(N1,ia,iPge,iPbe,iPgs,iPbs,iPbf))
						end if
					
					else 

						if (agrid(ia)>cutvalue) then
						call NewtonU2(AS2h(N1,ia,iPge,iPbe,iPgs,iPbs,iPbf), V2h(N1,ia,iPge,iPbe,iPgs,iPbs,iPbf))
						else
						call BisectU2(AS2h(N1,ia,iPge,iPbe,iPgs,iPbs,iPbf), V2h(N1,ia,iPge,iPbe,iPgs,iPbs,iPbf))						
						end if
						
											
					end if


	     end do
	end do
	end do
	end do
	end do
	end do
	end do	

	do Backtime = 1, N1-1
     time = N1-Backtime
	 kp = K(N1+1-Backtime)
     EV2 = 0.0D0
print*, time


		!	calculate expected value function	
		
		do iPge = 1, nP
		do iPbe = 1, nP
		do iPgs = 1, nP
		do iPbs = 1, nP
		do iPbf = 1, nP
         
		 ! Compute next period's Ps
				!first give today's P values.  THESE ARE UNCONDITIONAL VALUES!!!!!!
		 Pvalue(1,4) = pi(1,time)*Pgrid(iPbe)  !bad, e
		 Pvalue(1,3) = pi(1,time)*(1.0D0-Pgrid(iPbe))*Pgrid(iPbs)  !bad, s
		 Pvalue(1,2) = pi(1,time)*(1.0D0-Pgrid(iPbe))*(1.0D0-Pgrid(iPbs))*Pgrid(iPbf)  !bad, f
		 Pvalue(1,1) = pi(1,time)-Pvalue(1,4)-Pvalue(1,3)-Pvalue(1,2) !bad, l
		 Pvalue(2,4) = pi(2,time)*Pgrid(iPge)  !good, e
		 Pvalue(2,3) = pi(2,time)*(1.0D0-Pgrid(iPge))*Pgrid(iPgs)  !good, s
		 Pvalue(2,2) = 0.0D0  !good, f
		 Pvalue(2,1) = pi(2,time)-Pvalue(2,4)-Pvalue(2,3) !good, l

		 call ComputePprimes2()
		
		! I treat the first and last interval differently. These intervals are defined as (0,cutvect(1))
		! and (cutvect(ninterv-1),1). There are eight cutoff points but nine different values. 

         do ia = 1, na
		 do ibta = 1, nbta   
		      
			  do i = 1,  ninterv
			   if (i==1) then  
			   	do ibtap = 1, nbta
				if (ibtap==1) then
			    VPP(ia,ibtap) = lininterp5(Pprimege(i),Pprimebe(i),Pprimegs(i),Pprimebs(i),Pprimebf(i),Pgrid,Pgrid,Pgrid,Pgrid,Pgrid,V2l(N1+1-Backtime,ia,:,:,:,:,:))           
  	            else if (ibtap==2) then
				VPP(ia,ibtap) = lininterp5(Pprimege(i),Pprimebe(i),Pprimegs(i),Pprimebs(i),Pprimebf(i),Pgrid,Pgrid,Pgrid,Pgrid,Pgrid,V2m(N1+1-Backtime,ia,:,:,:,:,:))		
				else
			    VPP(ia,ibtap) = lininterp5(Pprimege(i),Pprimebe(i),Pprimegs(i),Pprimebs(i),Pprimebf(i),Pgrid,Pgrid,Pgrid,Pgrid,Pgrid,V2h(N1+1-Backtime,ia,:,:,:,:,:))
				end if

				EV2(ia,iPge,iPbe,iPgs,iPbs,iPbf,ibta) =  EV2(ia,iPge,iPbe,iPgs,iPbs,iPbf,ibta)+(cutvect(1))*VPP(ia,ibtap)*trBta(ibta,ibtap)
	     	   end do
			   else if (i==ninterv) then
			   do ibtap = 1, nbta
				if (ibtap==1) then
			    VPP(ia,ibtap) = lininterp5(Pprimege(i),Pprimebe(i),Pprimegs(i),Pprimebs(i),Pprimebf(i),Pgrid,Pgrid,Pgrid,Pgrid,Pgrid,V2l(N1+1-Backtime,ia,:,:,:,:,:))           
  	            else if (ibtap==2) then
				VPP(ia,ibtap) = lininterp5(Pprimege(i),Pprimebe(i),Pprimegs(i),Pprimebs(i),Pprimebf(i),Pgrid,Pgrid,Pgrid,Pgrid,Pgrid,V2m(N1+1-Backtime,ia,:,:,:,:,:))		
				else
			    VPP(ia,ibtap) = lininterp5(Pprimege(i),Pprimebe(i),Pprimegs(i),Pprimebs(i),Pprimebf(i),Pgrid,Pgrid,Pgrid,Pgrid,Pgrid,V2h(N1+1-Backtime,ia,:,:,:,:,:))
				end if

  	            EV2(ia,iPge,iPbe,iPgs,iPbs,iPbf,ibta) =  EV2(ia,iPge,iPbe,iPgs,iPbs,iPbf,ibta)+(1.0D0-cutvect(ninterv-1))*VPP(ia,ibtap)*trBta(ibta,ibtap)
				
			   end do
			   else 
			    do ibtap = 1, nbta
				if (ibtap==1) then
				VPP(ia,ibtap) = lininterp5(Pprimege(i),Pprimebe(i),Pprimegs(i),Pprimebs(i),Pprimebf(i),Pgrid,Pgrid,Pgrid,Pgrid,Pgrid,V2l(N1+1-Backtime,ia,:,:,:,:,:))           
  	            else if (ibtap==2) then
				VPP(ia,ibtap) = lininterp5(Pprimege(i),Pprimebe(i),Pprimegs(i),Pprimebs(i),Pprimebf(i),Pgrid,Pgrid,Pgrid,Pgrid,Pgrid,V2m(N1+1-Backtime,ia,:,:,:,:,:))		
				else
				VPP(ia,ibtap) = lininterp5(Pprimege(i),Pprimebe(i),Pprimegs(i),Pprimebs(i),Pprimebf(i),Pgrid,Pgrid,Pgrid,Pgrid,Pgrid,V2h(N1+1-Backtime,ia,:,:,:,:,:))
				end if
  	            EV2(ia,iPge,iPbe,iPgs,iPbs,iPbf,ibta) =  EV2(ia,iPge,iPbe,iPgs,iPbs,iPbf,ibta)+(cutvect(i)-cutvect(i-1))*VPP(ia,ibtap)*trBta(ibta,ibtap)
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


    call setsplineV2()

    do ia = 1, na
	do iPge = 1, nP
    do iPbe = 1, nP
	do iPgs = 1, nP
	do iPbs = 1, nP
	do iPbf = 1, nP


		income = (1.0+irate2(N1-Backtime))*agrid(ia) + wage2(N1-Backtime)*(pi(2,time)*Pgrid(iPge)+pi(1,time)*Pgrid(iPbe))*hour &
				         +uih*(pi(2,time)*(1.0D0-Pgrid(iPge))*Pgrid(iPgs)+pi(1,time)*(1.0D0-Pgrid(iPbe))*Pgrid(iPbs)+pi(1,time)*(1.0D0-Pgrid(iPbe))*(1.0D0-Pgrid(iPbs))*Pgrid(iPbf)) &
						 +uil*(1.0-pi(2,time)*Pgrid(iPge)-pi(1,time)*Pgrid(iPbe)-pi(2,time)*(1.0D0-Pgrid(iPge))*Pgrid(iPgs)-pi(1,time)*(1.0D0-Pgrid(iPbe))*Pgrid(iPbs)-pi(1,time)*(1.0D0-Pgrid(iPbe))*(1.0D0-Pgrid(iPbs))*Pgrid(iPbf))



	     do ibta = 1,nbta
         
					if (ibta==1) then

						if (agrid(ia)>cutvalue) then
						call NewtonU2(AS2l(N1-Backtime,ia,iPge,iPbe,iPgs,iPbs,iPbf), V2l(N1-Backtime,ia,iPge,iPbe,iPgs,iPbs,iPbf))
						else
						call BisectU2(AS2l(N1-Backtime,ia,iPge,iPbe,iPgs,iPbs,iPbf), V2l(N1-Backtime,ia,iPge,iPbe,iPgs,iPbs,iPbf))
						end if
						
					else if (ibta==2) then

						if (agrid(ia)>cutvalue) then
						call NewtonU2(AS2m(N1-Backtime,ia,iPge,iPbe,iPgs,iPbs,iPbf), V2m(N1-Backtime,ia,iPge,iPbe,iPgs,iPbs,iPbf))
						else
						call BisectU2(AS2m(N1-Backtime,ia,iPge,iPbe,iPgs,iPbs,iPbf), V2m(N1-Backtime,ia,iPge,iPbe,iPgs,iPbs,iPbf))
						end if
	
					else 

						if (agrid(ia)>cutvalue) then
						call NewtonU2(AS2h(N1-Backtime,ia,iPge,iPbe,iPgs,iPbs,iPbf), V2h(N1-Backtime,ia,iPge,iPbe,iPgs,iPbs,iPbf))
						else
						call BisectU2(AS2h(N1-Backtime,ia,iPge,iPbe,iPgs,iPbs,iPbf), V2h(N1-Backtime,ia,iPge,iPbe,iPgs,iPbs,iPbf))						
						end if
						
											
					end if
           
	     end do
	end do
	end do
	end do
	end do
	end do
	end do

	
    end do

	if (final) then
		call SaveValueFunction()
	end if

end subroutine

	    






