!	----------------------------------------------------------------------
!	File name: BrentU.f90
!	
!	----------------------------------------------------------------------


subroutine BrentU2(ap2, vf2)

	use Globals

	use Numerical_Libraries
	implicit none

	double precision FCNU2
	real(8), parameter:: CGOLD=.3819660,ZEPS=1.0e-12, tol=1.0e-12, RR = 0.61803399D0, C = 1.0D0-RR
	integer, parameter:: ITMAX=100
    REAL(8) ap2, vf2, ax, bx, cx
    INTEGER iterB
    REAL(8) a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,uu,vv,ww,xx,xm
	real(8) x0, x1, x2, x3, ff0, ff1, ff2, ff3
    logical ifstop, nobracket


    x0= max(amin,AS2(time,ia,iPg,iPb,ibta)-0.5)
	x3= min(amax,AS2(time,ia,iPg,iPb,ibta)+0.5)

	x1 = (x0+x3)/2.0D0


!
	ff0 = FCNU2(x0)
	ff1 = FCNU2(x1)
	ff3 = FCNU2(x3)



	nobracket= .true.
	if (ff1>max(ff0,ff3)) then !corner
		if (ff0>ff3) then
		ap2 = x3
!
		vf2=-FCNU2(ap2)
		else
		ap2 = x0
!
		vf2=-FCNU2(ap2)
		end if
	else ! start the goldensecion loop. Once bracketed, switch to brent 
		ifstop = .true.
		x2 = RR*x1+C*x3
!
		ff2 = FCNU2(x2)
		do while (ifstop .and. nobracket)
		   if (ff1>ff2) then
		   x0=x1
		   x1=x2
		   x2=RR*x1+C*x3
		   ff1=ff2
!
		   ff2=FCNU2(x2)
		   	if (ff1<min(ff0,ff2)) then ! bracketed
			nobracket = .false.
			else
			nobracket = .true.
			end if
		   else
		   x3=x2
		   x2=x1
		   x1=RR*x2+C*x0
		   ff2=ff1
!
		   ff1=FCNU2(x1)
		    if (ff1<min(ff0,ff2)) then ! bracketed
			nobracket = .false.
			else
			nobracket = .true.
			end if
		   end if
		
		   if (abs(x3-x0)<TOL) then
		       if (ff1<ff2) then
				ap2 = x1
				vf2 = -ff1
			   else
			    ap2 = x2
				vf2 = -ff2
			   end if
		    
		   ifstop = .false.
		   end if

		end do

	  if (ifstop) then !start brent
	  ax=x0
	  cx=x2
	  
	  a=min(ax,cx)
      b=max(ax,cx)
      vv=x1
      ww=vv
      xx=vv
      e=0.

!
      fx=FCNU2(xx)
      fv=fx
      fw=fx
      do 11 iterB=1,ITMAX
        xm=0.5*(a+b)
        tol1=tol*abs(xx)+ZEPS
        tol2=2.*tol1
        if(abs(xx-xm).le.(tol2-.5*(b-a))) goto 3
        if(abs(e).gt.tol1) then
          r=(xx-ww)*(fx-fv)
          q=(xx-vv)*(fx-fw)
          p=(xx-vv)*q-(xx-ww)*r
          q=2.*(q-r)
          if(q.gt.0.) p=-p
          q=abs(q)
          etemp=e
          e=d
          if(abs(p).ge.abs(.5*q*etemp).or.p.le.q*(a-xx).or.p.ge.q*(b-xx)) goto 1
          d=p/q
          uu=xx+d
          if(uu-a.lt.tol2 .or. b-uu.lt.tol2) d=sign(tol1,xm-xx)
          goto 2
        endif
1       if(xx.ge.xm) then
          e=a-xx
        else
          e=b-xx
        endif
        d=CGOLD*e
2       if(abs(d).ge.tol1) then
          uu=xx+d
        else
          uu=xx+sign(tol1,d)
        endif
!
        fu=FCNU2(uu)
        if(fu.le.fx) then
          if(uu.ge.xx) then
            a=xx
          else
            b=xx
          endif
          vv=ww
          fv=fw
          ww=xx
          fw=fx
          xx=uu
          fx=fu
        else
          if(uu.lt.xx) then
            a=uu
          else
            b=uu
          endif
          if(fu.le.fw .or. ww.eq.xx) then
            vv=ww
            fv=fw
            ww=uu
            fw=fu
          else if(fu.le.fv .or. vv.eq.xx .or. vv.eq.ww) then
            vv=uu
            fv=fu
          endif
        endif
11    continue
      print*, 'brent exceed maximum iterations'
3     ap2=xx
      vf2=-fx
 
 end if
 end if

end subroutine

