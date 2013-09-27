!	----------------------------------------------------------------------
!	File name: ComputeWelfare.f90
!	----------------------------------------------------------------------


subroutine ComputeWelfare()

	use Globals
	use LinInterpModule
	use PolyInterpModule
	use Numerical_Libraries
	
	implicit none

	integer, parameter:: xseed = 2
	real(8), parameter:: tiny = 1.0D-20
	
    real(8) ll, V0(na,neps,nbta,nk,nz)
	real(8) AS0(na,neps,nbta,nk,nz)
	real(8) VA(na), oldmul(Nindiv), newmul(Nindiv), Bsum(Nindiv), B(3), lambda(Nindiv), sumlambda
	real(8) lambdaeps1, lambdaeps0, lambdabeta1, lambdabeta2, lambdabeta3
	real(8) Asorted(Nindiv)
	integer Num1,Num05, Num5, Num25, Num50, Num75, Num95, Num995, Num99, Num999
	integer  Nindiveps0,  Nindiveps1,  Nindivbeta1, Nindivbeta2, Nindivbeta3
	real(8) lambdasum1, lambdasum15, lambdasum525, lambdasum2550
	real(8) lambdasum5075, lambdasum7595, lambdasum9599, lambdasum99 
	integer IPERM(Nindiv)
	real(8) lambdasum1eps0, lambdasum15eps0, lambdasum525eps0, lambdasum2550eps0,lambdasum5075eps0,lambdasum7595eps0
    real(8) lambdasum9599eps0, lambdasum99eps0 
    real(8) lambdasum1eps1, lambdasum15eps1, lambdasum525eps1, lambdasum2550eps1, lambdasum5075eps1, lambdasum7595eps1
    real(8) lambdasum9599eps1, lambdasum99eps1 

	integer Num1eps0, Num15eps0, Num525eps0, Num2550eps0, Num5075eps0, Num7595eps0, Num9599eps0, Num99eps0 
    integer Num1eps1, Num15eps1, Num525eps1, Num2550eps1, Num5075eps1, Num7595eps1, Num9599eps1, Num99eps1 

	real(8) lambdaCeps1beta1, lambdaCeps1beta2, lambdaCeps1beta3, lambdaCeps0beta1, lambdaCeps0beta2, lambdaCeps0beta3
    real(8) lambda05eps1beta1, lambda05eps1beta2, lambda05eps1beta3, lambda05eps0beta1, lambda05eps0beta2, lambda05eps0beta3
	real(8) lambda5eps1beta1, lambda5eps1beta2, lambda5eps1beta3, lambda5eps0beta1, lambda5eps0beta2, lambda5eps0beta3
	real(8) lambda50eps1beta1, lambda50eps1beta2, lambda50eps1beta3, lambda50eps0beta1, lambda50eps0beta2, lambda50eps0beta3
	real(8) lambda95eps1beta1, lambda95eps1beta2, lambda95eps1beta3, lambda95eps0beta1, lambda95eps0beta2, lambda95eps0beta3
	real(8) lambda995eps1beta1, lambda995eps1beta2, lambda995eps1beta3, lambda995eps0beta1, lambda995eps0beta2, lambda995eps0beta3
	real(8) lambda999eps1beta1, lambda999eps1beta2, lambda999eps1beta3, lambda999eps0beta1, lambda999eps0beta2, lambda999eps0beta3

	real(8) VC_f, VC_s 

	real(8), dimension(:), allocatable:: V_s, V_f
	real(8), dimension(:,:), allocatable:: CrA
	integer, dimension(:,:), allocatable:: CrBta, CrEps, CrPg, CrPb




	!	allocate space for cross section data generation

	allocate(CrA(Nindiv,2), CrEps(Nindiv,2), CrPg(Nindiv,2), CrPb(Nindiv,2), CrBta(Nindiv,2),		&
			  V_s(Nindiv),  V_f(Nindiv))
    
		! Initial state assumed to be good. Also change this in generateseries if the initial state is a bad state
	call rnset(xseed)

	!   Read the initial distribution from the file 
	
	open(1, file='Input\PanelA112B.txt', status = 'old')
    do indiv = 1, Nindiv
     read(1, *) CrA(indiv,1)
      CrA(indiv,2) = 0.0D0
    end do


    open(1, file='Input\PanelEps112B.txt', status = 'old')
    do indiv = 1, Nindiv
      read(1, *) CrEps(indiv,1)
      CrEps(indiv,2) = 0
    end do
 
    open(1, file='Input\PanelBta112B.txt', status = 'old')
    do indiv = 1, Nindiv
      read(1, *) CrBta(indiv,1)
      CrBta(indiv,2) = 0
    end do


 ! Initial aggeraget states-CHANGE FOR B STATE
 	Zdata(1) = 1
	Kdata(1) =  real(sum(CrA(:,1))/Nindiv)
	
    CrPg = 1
	CrPb = 1

    do indiv = 1, Nindiv
     
     if (CrEps(indiv,1) == 2 .and. Zdata(1) == 2) then 
	     CrPg(indiv,1) = nP
		 CrPb(indiv,1) = 1
	 end if
	 
	 if (CrEps(indiv,1) == 2 .and. Zdata(1) == 1) then 
	     CrPg(indiv,1) = 1
		 CrPb(indiv,1) = nP
     end if
   
    end do

    
    do indiv = 1, Nindiv
     IPERM(indiv)  = indiv
	end do

	CALL DSVRGP (Nindiv, CrA(:,1),Asorted, IPERM)

    Num1 = int(0.01*real(Nindiv))
	Num05 = int(0.005*real(Nindiv))
	Num5 = int(0.05*real(Nindiv))
	Num25 = int(0.25*real(Nindiv))
	Num50 = int(0.5*real(Nindiv))
	Num75 = int(0.75*real(Nindiv))
	Num95 = int(0.95*real(Nindiv))
	Num995 = int(0.995*real(Nindiv))
	Num99 = int(0.99*real(Nindiv))
	Num999 = int(0.999*real(Nindiv))



   B(1) = 75.73359878
   B(2) = 93.91519803
   B(3) = 121.96313901
   
   sumlambda = 0.0
   lambdaeps1 = 0.0
   lambdaeps0 = 0.0
   lambdabeta1 = 0.0
   lambdabeta2 = 0.0
   lambdabeta3 = 0.0
   
   Nindiveps0 = 0
   Nindiveps1 = 0  
   Nindivbeta1 = 0
   Nindivbeta2 = 0
   Nindivbeta3 = 0
   
   lambdasum1 = 0.0D0
   lambdasum15 = 0.0D0
   lambdasum525 = 0.0D0
   lambdasum2550 = 0.0D0
   lambdasum5075 = 0.0D0
   lambdasum7595 = 0.0D0
   lambdasum9599 = 0.0D0
   lambdasum99 = 0.0D0

   lambdasum1eps0 = 0.0D0
   lambdasum15eps0 = 0.0D0
   lambdasum525eps0 = 0.0D0
   lambdasum2550eps0 = 0.0D0
   lambdasum5075eps0 = 0.0D0
   lambdasum7595eps0 = 0.0D0
   lambdasum9599eps0 = 0.0D0
   lambdasum99eps0 = 0.0D0

   lambdasum1eps1 = 0.0D0
   lambdasum15eps1 = 0.0D0
   lambdasum525eps1 = 0.0D0
   lambdasum2550eps1 = 0.0D0
   lambdasum5075eps1 = 0.0D0
   lambdasum7595eps1 = 0.0D0
   lambdasum9599eps1 = 0.0D0
   lambdasum99eps1 = 0.0D0

   Num1eps0 = 0
   Num15eps0 = 0
   Num525eps0 = 0
   Num2550eps0 = 0
   Num5075eps0 = 0
   Num7595eps0 = 0
   Num9599eps0 = 0
   Num99eps0 = 0

   Num1eps1 = 0
   Num15eps1 = 0
   Num525eps1 = 0
   Num2550eps1 = 0
   Num5075eps1 = 0
   Num7595eps1 = 0
   Num9599eps1 = 0
   Num99eps1 = 0

   lambdaCeps1beta1 = 0.0
   lambdaCeps1beta2 = 0.0 
   lambdaCeps1beta3 = 0.0
   lambdaCeps0beta1 = 0.0 
   lambdaCeps0beta2 = 0.0
   lambdaCeps0beta3 = 0.0

   

   lambda05eps1beta1 = 0.0
   lambda05eps1beta2 = 0.0 
   lambda05eps1beta3 = 0.0
   lambda05eps0beta1 = 0.0 
   lambda05eps0beta2 = 0.0
   lambda05eps0beta3 = 0.0

  

   lambda5eps1beta1 = 0.0
   lambda5eps1beta2 = 0.0 
   lambda5eps1beta3 = 0.0
   lambda5eps0beta1 = 0.0 
   lambda5eps0beta2 = 0.0
   lambda5eps0beta3 = 0.0

   lambda50eps1beta1 = 0.0
   lambda50eps1beta2 = 0.0 
   lambda50eps1beta3 = 0.0
   lambda50eps0beta1 = 0.0 
   lambda50eps0beta2 = 0.0
   lambda50eps0beta3 = 0.0

   
   lambda95eps1beta1 = 0.0
   lambda95eps1beta2 = 0.0 
   lambda95eps1beta3 = 0.0
   lambda95eps0beta1 = 0.0 
   lambda95eps0beta2 = 0.0
   lambda95eps0beta3 = 0.0

  

   lambda995eps1beta1 = 0.0
   lambda995eps1beta2 = 0.0 
   lambda995eps1beta3 = 0.0
   lambda995eps0beta1 = 0.0 
   lambda995eps0beta2 = 0.0
   lambda995eps0beta3 = 0.0


   lambda999eps1beta1 = 0.0
   lambda999eps1beta2 = 0.0 
   lambda999eps1beta3 = 0.0
   lambda999eps0beta1 = 0.0 
   lambda999eps0beta2 = 0.0
   lambda999eps0beta3 = 0.0

  

!   Now we have the initial distributions and the initial value functions

	
	do indiv = 1, Nindiv
          do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,CrEps(IPERM(indiv),1),CrBta(IPERM(indiv),1),:, Zdata(1)))
		  end do
         
		  V_f(indiv) = lininterp1(Asorted(indiv), agrid, VA(:))
		  V_s(indiv) = lininterp1(Asorted(indiv), agrid, V2(1,:,CrPg(IPERM(indiv),1),CrPb(IPERM(indiv),1),CrBta(IPERM(indiv),1)))
		  lambda(indiv) = (dexp((V_s(indiv)-V_f(indiv))/B(CrBta(IPERM(indiv),1)))-1.0)*100	 
		  sumlambda = sumlambda+lambda(indiv)
		  
		  ! Table 1
		  if (CrEps(IPERM(indiv),1) == 1) then 
		   lambdaeps0 = lambdaeps0+lambda(indiv)
		   Nindiveps0 =  Nindiveps0+1
		  else
		   lambdaeps1 = lambdaeps1+lambda(indiv) 
		   Nindiveps1 =  Nindiveps1+1
		  end if
          
		  if (CrBta(IPERM(indiv),1) == 1) then 
		   lambdabeta1 = lambdabeta1+lambda(indiv)
           Nindivbeta1 = Nindivbeta1+1
		  elseif (CrBta(IPERM(indiv),1) == 2) then 
		   lambdabeta2 = lambdabeta2+lambda(indiv)
		   Nindivbeta2 = Nindivbeta2+1
		  else 
           lambdabeta3 = lambdabeta3+lambda(indiv)
		   Nindivbeta3 = Nindivbeta3+1
		  end if
          
		  ! Table 2
		  if (indiv<Num1) then 
		   lambdasum1 = lambdasum1+lambda(indiv)
		   if (CrEps(IPERM(indiv),1) ==1) then 
		    lambdasum1eps0 = lambdasum1eps0+lambda(indiv)
			Num1eps0 = Num1eps0+1
		   else 
            lambdasum1eps1 = lambdasum1eps1+lambda(indiv)
			Num1eps1 = Num1eps1+1
		   endif 
		  
		  else if (indiv>=Num1 .and. indiv<Num5) then 
		   lambdasum15 = lambdasum15+lambda(indiv)
		   if (CrEps(IPERM(indiv),1) ==1) then 
		    lambdasum15eps0 = lambdasum15eps0+lambda(indiv)
			Num15eps0 = Num15eps0+1
		   else 
            lambdasum15eps1 = lambdasum15eps1+lambda(indiv)
			Num15eps1 = Num15eps1+1
		   endif 
		  
		  elseif (indiv>=Num5 .and. indiv<Num25) then 
		   lambdasum525 = lambdasum525+lambda(indiv)
		   if (CrEps(IPERM(indiv),1) ==1) then 
		    lambdasum525eps0 = lambdasum525eps0+lambda(indiv)
			Num525eps0 = Num525eps0+1
		   else 
            lambdasum525eps1 = lambdasum525eps1+lambda(indiv)
			Num525eps1 = Num525eps1+1
		   endif 

		  elseif (indiv>=Num25 .and. indiv<Num50) then 
		   lambdasum2550 = lambdasum2550+lambda(indiv)
            if (CrEps(IPERM(indiv),1) ==1) then 
		    lambdasum2550eps0 = lambdasum2550eps0+lambda(indiv)
			Num2550eps0 = Num2550eps0+1
		   else 
            lambdasum2550eps1 = lambdasum2550eps1+lambda(indiv)
			Num2550eps1 = Num2550eps1+1
		   endif 

		  elseif (indiv>=Num50 .and. indiv<Num75) then 
		   lambdasum5075 = lambdasum5075+lambda(indiv)
		   if (CrEps(IPERM(indiv),1) ==1) then 
		    lambdasum5075eps0 = lambdasum5075eps0+lambda(indiv)
			Num5075eps0 = Num5075eps0+1
		   else 
            lambdasum5075eps1 = lambdasum5075eps1+lambda(indiv)
			Num5075eps1 = Num5075eps1+1
		   endif 

		  elseif (indiv>=Num75 .and. indiv<Num95) then 
		   	lambdasum7595 = lambdasum7595+lambda(indiv)
			if (CrEps(IPERM(indiv),1) ==1) then 
		    lambdasum7595eps0 = lambdasum7595eps0+lambda(indiv)
			Num7595eps0 = Num7595eps0+1
		   else 
            lambdasum7595eps1 = lambdasum7595eps1+lambda(indiv)
			Num7595eps1 = Num7595eps1+1
		   endif 

		  elseif (indiv>=Num95 .and. indiv<Num99) then 
           lambdasum9599 = lambdasum9599+lambda(indiv)
		   if (CrEps(IPERM(indiv),1) ==1) then 
		    lambdasum9599eps0 = lambdasum9599eps0+lambda(indiv)
			Num9599eps0 = Num9599eps0+1
		   else 
            lambdasum9599eps1 = lambdasum9599eps1+lambda(indiv)
			Num9599eps1 = Num9599eps1+1
		   endif 
		  else 
		   lambdasum99 = lambdasum99+lambda(indiv)
		   if (CrEps(IPERM(indiv),1) ==1) then 
		    lambdasum99eps0 = lambdasum99eps0+lambda(indiv)
			Num99eps0 = Num99eps0+1
		   else 
            lambdasum99eps1 = lambdasum99eps1+lambda(indiv)
			Num99eps1 = Num99eps1+1
		   endif 
		  endif 
     end do     
		  ! Table 3, V1(ia,ieps,ibta,ik,iz) and V2(time,ia,iPg,iPb,ibta)
          
        ! Constrained agents
		! eps=1, bta=1
                 
		  VC_f = polyinterp1(Kdata(1), kgrid, V1(1,2,1,:, Zdata(1)))
		  VC_s =  V2(1,1,nP,1,1)
          lambdaCeps1beta1 = (dexp((VC_s-VC_f)/B(1))-1.0)*100	
		
		! eps=1, bta=2          
		  VC_f = polyinterp1(Kdata(1), kgrid, V1(1,2,2,:, Zdata(1)))
		  VC_s =  V2(1,1,nP,1,2)
          lambdaCeps1beta2 = (dexp((VC_s-VC_f)/B(2))-1.0)*100	
		  
		! eps=1, bta=3             
		  VC_f = polyinterp1(Kdata(1), kgrid, V1(1,2,3,:, Zdata(1)))
		  VC_s =  V2(1,1,nP,1,3)
          lambdaCeps1beta3 = (dexp((VC_s-VC_f)/B(3))-1.0)*100  
		
		! eps=0, bta=1          
		  VC_f = polyinterp1(Kdata(1), kgrid, V1(1,1,1,:, Zdata(1)))
		  VC_s =  V2(1,1,1,1,1)
          lambdaCeps0beta1 = (dexp((VC_s-VC_f)/B(1))-1.0)*100	
		
		! eps=0, bta=2         
		  VC_f = polyinterp1(Kdata(1), kgrid, V1(1,1,2,:, Zdata(1)))
		  VC_s =  V2(1,1,1,1,2)
          lambdaCeps0beta2 = (dexp((VC_s-VC_f)/B(2))-1.0)*100	
		  
		! eps=0, bta=3           
		  VC_f = polyinterp1(Kdata(1), kgrid, V1(1,1,3,:, Zdata(1)))
		  VC_s =  V2(1,1,1,1,3)
          lambdaCeps0beta3 = (dexp((VC_s-VC_f)/B(3))-1.0)*100    
		  
		 
		! Agents at 0.005, V1(ia,ieps,ibta,ik,iz) and V2(time,ia,iPg,iPb,ibta)
	    
		! eps=1, bta=1
         do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,2,1,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num05), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num05), agrid, V2(1,:,nP,1,1))
          lambda05eps1beta1 = (dexp((VC_s-VC_f)/B(1))-1.0)*100
         
       	! eps=1, bta=2
         do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,2,2,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num05), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num05), agrid, V2(1,:,nP,1,2))
          lambda05eps1beta2 = (dexp((VC_s-VC_f)/B(2))-1.0)*100

		! eps=1, bta=3
          do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,2,3,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num05), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num05), agrid, V2(1,:,nP,1,3))
          lambda05eps1beta3 = (dexp((VC_s-VC_f)/B(3))-1.0)*100


        ! eps=0, bta=1
          do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,1,1,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num05), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num05), agrid, V2(1,:,1,1,1))
          lambda05eps0beta1 = (dexp((VC_s-VC_f)/B(1))-1.0)*100

        ! eps=0, bta=2
          do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,1,2,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num05), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num05), agrid, V2(1,:,1,1,2))
          lambda05eps0beta2 = (dexp((VC_s-VC_f)/B(2))-1.0)*100

        ! eps=0, bta=3
          do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,1,3,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num05), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num05), agrid, V2(1,:,1,1,3))
          lambda05eps0beta3 = (dexp((VC_s-VC_f)/B(3))-1.0)*100



        ! Agents at 0.05
	    
		! eps=1, bta=1
         do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,2,1,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num5), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num5), agrid, V2(1,:,nP,1,1))
          lambda5eps1beta1 = (dexp((VC_s-VC_f)/B(1))-1.0)*100
         
       	! eps=1, bta=2
         do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,2,2,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num5), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num5), agrid, V2(1,:,nP,1,2))
          lambda5eps1beta2 = (dexp((VC_s-VC_f)/B(2))-1.0)*100

		! eps=1, bta=3
          do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,2,3,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num5), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num5), agrid, V2(1,:,nP,1,3))
          lambda5eps1beta3 = (dexp((VC_s-VC_f)/B(3))-1.0)*100


        ! eps=0, bta=1
          do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,1,1,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num5), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num5), agrid, V2(1,:,1,1,1))
          lambda5eps0beta1 = (dexp((VC_s-VC_f)/B(1))-1.0)*100

        ! eps=0, bta=2
          do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,1,2,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num5), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num5), agrid, V2(1,:,1,1,2))
          lambda5eps0beta2 = (dexp((VC_s-VC_f)/B(2))-1.0)*100

        ! eps=0, bta=3
          do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,1,3,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num5), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num5), agrid, V2(1,:,1,1,3))
          lambda5eps0beta3 = (dexp((VC_s-VC_f)/B(3))-1.0)*100



        ! Agents at 0.5
	    
		! eps=1, bta=1
         do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,2,1,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num50), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num50), agrid, V2(1,:,nP,1,1))
          lambda50eps1beta1 = (dexp((VC_s-VC_f)/B(1))-1.0)*100
         
       	! eps=1, bta=2
         do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,2,2,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num50), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num50), agrid, V2(1,:,nP,1,2))
          lambda50eps1beta2 = (dexp((VC_s-VC_f)/B(2))-1.0)*100

		! eps=1, bta=3
          do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,2,3,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num50), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num50), agrid, V2(1,:,nP,1,3))
          lambda50eps1beta3 = (dexp((VC_s-VC_f)/B(3))-1.0)*100


        ! eps=0, bta=1
          do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,1,1,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num50), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num50), agrid, V2(1,:,1,1,1))
          lambda50eps0beta1 = (dexp((VC_s-VC_f)/B(1))-1.0)*100

        ! eps=0, bta=2
          do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,1,2,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num50), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num50), agrid, V2(1,:,1,1,2))
          lambda50eps0beta2 = (dexp((VC_s-VC_f)/B(2))-1.0)*100

        ! eps=0, bta=3
          do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,1,3,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num50), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num50), agrid, V2(1,:,1,1,3))
          lambda50eps0beta3 = (dexp((VC_s-VC_f)/B(3))-1.0)*100


		! Agents at 0.95
	    
		! eps=1, bta=1
         do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,2,1,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num95), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num95), agrid, V2(1,:,nP,1,1))
          lambda95eps1beta1 = (dexp((VC_s-VC_f)/B(1))-1.0)*100
         
       	! eps=1, bta=2
         do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,2,2,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num95), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num95), agrid, V2(1,:,nP,1,2))
          lambda95eps1beta2 = (dexp((VC_s-VC_f)/B(2))-1.0)*100

		! eps=1, bta=3
          do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,2,3,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num95), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num95), agrid, V2(1,:,nP,1,3))
          lambda95eps1beta3 = (dexp((VC_s-VC_f)/B(3))-1.0)*100


        ! eps=0, bta=1
          do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,1,1,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num95), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num95), agrid, V2(1,:,1,1,1))
          lambda95eps0beta1 = (dexp((VC_s-VC_f)/B(1))-1.0)*100

        ! eps=0, bta=2
          do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,1,2,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num95), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num95), agrid, V2(1,:,1,1,2))
          lambda95eps0beta2 = (dexp((VC_s-VC_f)/B(2))-1.0)*100

        ! eps=0, bta=3
          do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,1,3,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num95), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num95), agrid, V2(1,:,1,1,3))
          lambda95eps0beta3 = (dexp((VC_s-VC_f)/B(3))-1.0)*100

 
        ! Agents at 0.995
	    
		! eps=1, bta=1
         do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,2,1,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num995), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num995), agrid, V2(1,:,nP,1,1))
          lambda995eps1beta1 = (dexp((VC_s-VC_f)/B(1))-1.0)*100
         
       	! eps=1, bta=2
         do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,2,2,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num995), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num995), agrid, V2(1,:,nP,1,2))
          lambda995eps1beta2 = (dexp((VC_s-VC_f)/B(2))-1.0)*100

		! eps=1, bta=3
          do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,2,3,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num995), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num995), agrid, V2(1,:,nP,1,3))
          lambda995eps1beta3 = (dexp((VC_s-VC_f)/B(3))-1.0)*100


        ! eps=0, bta=1
          do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,1,1,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num995), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num995), agrid, V2(1,:,1,1,1))
          lambda995eps0beta1 = (dexp((VC_s-VC_f)/B(1))-1.0)*100

        ! eps=0, bta=2
          do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,1,2,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num995), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num995), agrid, V2(1,:,1,1,2))
          lambda995eps0beta2 = (dexp((VC_s-VC_f)/B(2))-1.0)*100

        ! eps=0, bta=2
          do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,1,3,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num995), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num995), agrid, V2(1,:,1,1,3))
          lambda995eps0beta3 = (dexp((VC_s-VC_f)/B(3))-1.0)*100

		  
       ! Agents at 0.999
	    
		! eps=1, bta=1
         do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,2,1,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num999), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num999), agrid, V2(1,:,nP,1,1))
          lambda999eps1beta1 = (dexp((VC_s-VC_f)/B(1))-1.0)*100
         
       	! eps=1, bta=2
         do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,2,2,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num999), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num999), agrid, V2(1,:,nP,1,2))
          lambda999eps1beta2 = (dexp((VC_s-VC_f)/B(2))-1.0)*100

		! eps=1, bta=3
          do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,2,3,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num999), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num999), agrid, V2(1,:,nP,1,3))
          lambda999eps1beta3 = (dexp((VC_s-VC_f)/B(3))-1.0)*100


        ! eps=0, bta=1
          do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,1,1,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num999), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num999), agrid, V2(1,:,1,1,1))
          lambda999eps0beta1 = (dexp((VC_s-VC_f)/B(1))-1.0)*100

        ! eps=0, bta=2
          do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,1,2,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num999), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num999), agrid, V2(1,:,1,1,2))
          lambda999eps0beta2 = (dexp((VC_s-VC_f)/B(2))-1.0)*100

        ! eps=0, bta=3
          do ia = 1, na               
		   VA(ia) = polyinterp1(Kdata(1), kgrid, V1(ia,1,3,:, Zdata(1)))
		  end do
         
		  VC_f = lininterp1(Asorted(Num999), agrid, VA(:))
		  VC_s = lininterp1(Asorted(Num999), agrid, V2(1,:,1,1,3))
          lambda999eps0beta3 = (dexp((VC_s-VC_f)/B(3))-1.0)*100

      
   	
	

    open(1, file='Output\lambda.txt', status='unknown')
	open(2, file='Output\sumlambda.txt', status='unknown')
	open(3, file='Output\lambdaepsbeta.txt', status='unknown')
	
    
    do indiv = 1, Nindiv
	
		 write(1, '(F15.6)') lambda(indiv)
         write(3,'(I5,F15.6,I2,I2)') indiv, lambda(indiv), CrEps(IPERM(indiv),1), CrBta(IPERM(indiv),1)

	end do 
   
	write(2, '(/A,F10.6)')"lambda= ", sumlambda/Nindiv


	write(2, '(/A,F10.6)')"lambdaeps0= ", lambdaeps0/Nindiveps0
	write(2, '(/A,F10.6)')"lambdaeps1= ", lambdaeps1/Nindiveps1
	write(2, '(/A,F10.6)')"lambdabeta1= ", lambdabeta1/Nindivbeta1
	write(2, '(/A,F10.6)')"lambdabeta2= ", lambdabeta2/Nindivbeta2
	write(2, '(/A,F10.6)')"lambdabeta3= ", lambdabeta3/Nindivbeta3
	write(2, '(/A,F10.6)')"lambdasum1= ", lambdasum1/Num1
	write(2, '(/A,F10.6)')"lambdasum1_5= ", lambdasum15/(Num5-Num1)
	write(2, '(/A,F10.6)')"lambdasum5_25= ", lambdasum525/(Num25-Num5)
	write(2, '(/A,F10.6)')"lambdasum25_50= ", lambdasum2550/(Num50-Num25)
	write(2, '(/A,F10.6)')"lambdasum50_75= ", lambdasum5075/(Num75-Num50)
	write(2, '(/A,F10.6)')"lambdasum75_95= ", lambdasum7595/(Num95-Num75)
	write(2, '(/A,F10.6)')"lambdasum95_99= ", lambdasum9599/(Num99-Num95)
	write(2, '(/A,F10.6)')"lambdasum99= ", lambdasum99/(Nindiv-Num99)

	write(2, '(/A,F10.6)')"lambdasum1eps0= ", lambdasum1eps0/Num1eps0
	write(2, '(/A,F10.6)')"lambdasum1_5eps0= ", lambdasum15eps0/Num15eps0
	write(2, '(/A,F10.6)')"lambdasum5_25eps0= ", lambdasum525eps0/Num525eps0
	write(2, '(/A,F10.6)')"lambdasum25_50eps0= ", lambdasum2550eps0/Num2550eps0
	write(2, '(/A,F10.6)')"lambdasum50_75eps0= ", lambdasum5075eps0/Num5075eps0
	write(2, '(/A,F10.6)')"lambdasum75_95eps0= ", lambdasum7595eps0/Num5075eps0
	write(2, '(/A,F10.6)')"lambdasum95_99eps0= ", lambdasum9599eps0/Num9599eps0
	write(2, '(/A,F10.6)')"lambdasum99eps0= ", lambdasum99eps0/Num99eps0 	

	write(2, '(/A,F10.6)')"lambdasum1eps1= ", lambdasum1eps1/Num1eps1
	write(2, '(/A,F10.6)')"lambdasum1_5eps1= ", lambdasum15eps1/Num15eps1
	write(2, '(/A,F10.6)')"lambdasum5_25eps1= ", lambdasum525eps1/Num525eps1
	write(2, '(/A,F10.6)')"lambdasum25_50eps1= ", lambdasum2550eps1/Num2550eps1
	write(2, '(/A,F10.6)')"lambdasum50_75eps1= ", lambdasum5075eps1/Num5075eps1
	write(2, '(/A,F10.6)')"lambdasum75_95eps1= ", lambdasum7595eps1/Num5075eps1
	write(2, '(/A,F10.6)')"lambdasum95_99eps1= ", lambdasum9599eps1/Num9599eps1
	write(2, '(/A,F10.6)')"lambdasum99eps1= ", lambdasum99eps1/Num99eps1	

    write(2, '(/A,F10.6)')"lambdaCeps1beta1= ", lambdaCeps1beta1
	write(2, '(/A,F10.6)')"lambdaCeps1beta2= ", lambdaCeps1beta2
	write(2, '(/A,F10.6)')"lambdaCeps1beta3= ", lambdaCeps1beta3
	write(2, '(/A,F10.6)')"lambdaCeps0beta1= ", lambdaCeps0beta1
	write(2, '(/A,F10.6)')"lambdaCeps0beta2= ", lambdaCeps0beta2
	write(2, '(/A,F10.6)')"lambdaCeps0beta3= ", lambdaCeps0beta3

	write(2, '(/A,F10.6)')"lambda05eps1beta1= ", lambda05eps1beta1
	write(2, '(/A,F10.6)')"lambda05eps1beta2= ", lambda05eps1beta2
	write(2, '(/A,F10.6)')"lambda05eps1beta3= ", lambda05eps1beta3
	write(2, '(/A,F10.6)')"lambda05eps0beta1= ", lambda05eps0beta1
	write(2, '(/A,F10.6)')"lambda05eps0beta2= ", lambda05eps0beta2
	write(2, '(/A,F10.6)')"lambda05eps0beta3= ", lambda05eps0beta3
	
	write(2, '(/A,F10.6)')"lambda5eps1beta1= ", lambda5eps1beta1
	write(2, '(/A,F10.6)')"lambda5eps1beta2= ", lambda5eps1beta2
	write(2, '(/A,F10.6)')"lambda5eps1beta3= ", lambda5eps1beta3
	write(2, '(/A,F10.6)')"lambda5eps0beta1= ", lambda5eps0beta1
	write(2, '(/A,F10.6)')"lambda5eps0beta2= ", lambda5eps0beta2
	write(2, '(/A,F10.6)')"lambda5eps0beta3= ", lambda5eps0beta3

	write(2, '(/A,F10.6)')"lambda5eps1beta1= ", lambda50eps1beta1
	write(2, '(/A,F10.6)')"lambda5eps1beta2= ", lambda50eps1beta2
	write(2, '(/A,F10.6)')"lambda5eps1beta3= ", lambda50eps1beta3
	write(2, '(/A,F10.6)')"lambda5eps0beta1= ", lambda50eps0beta1
	write(2, '(/A,F10.6)')"lambda5eps0beta2= ", lambda50eps0beta2
	write(2, '(/A,F10.6)')"lambda5eps0beta3= ", lambda50eps0beta3


	write(2, '(/A,F10.6)')"lambda95eps1beta1= ", lambda95eps1beta1
	write(2, '(/A,F10.6)')"lambda95eps1beta2= ", lambda95eps1beta2
	write(2, '(/A,F10.6)')"lambda95eps1beta3= ", lambda95eps1beta3
	write(2, '(/A,F10.6)')"lambda95eps0beta1= ", lambda95eps0beta1
	write(2, '(/A,F10.6)')"lambda95eps0beta2= ", lambda95eps0beta2
	write(2, '(/A,F10.6)')"lambda95eps0beta3= ", lambda95eps0beta3

	write(2, '(/A,F10.6)')"lambda995eps1beta1= ", lambda995eps1beta1
	write(2, '(/A,F10.6)')"lambda995eps1beta2= ", lambda995eps1beta2
	write(2, '(/A,F10.6)')"lambda995eps1beta3= ", lambda995eps1beta3
	write(2, '(/A,F10.6)')"lambda995eps0beta1= ", lambda995eps0beta1
	write(2, '(/A,F10.6)')"lambda995eps0beta2= ", lambda995eps0beta2
	write(2, '(/A,F10.6)')"lambda995eps0beta3= ", lambda995eps0beta3

	write(2, '(/A,F10.6)')"lambda999eps1beta1= ", lambda999eps1beta1
	write(2, '(/A,F10.6)')"lambda999eps1beta2= ", lambda999eps1beta2
	write(2, '(/A,F10.6)')"lambda999eps1beta3= ", lambda999eps1beta3
	write(2, '(/A,F10.6)')"lambda999eps0beta1= ", lambda999eps0beta1
	write(2, '(/A,F10.6)')"lambda999eps0beta2= ", lambda999eps0beta2
	write(2, '(/A,F10.6)')"lambda999eps0beta3= ", lambda999eps0beta3
	

	!	deallocate memory

	deallocate(CrA, CrEps, CrPg, CrPb, CrBta, V_s,  V_f)





end subroutine
