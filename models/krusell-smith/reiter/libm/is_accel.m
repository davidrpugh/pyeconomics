function i = is_accel(iter,nAccel,nStart,distia)
  if(iter<=nStart)
    i = false;
  else
    i = mod(iter,nAccel+1)~=0;
    if(nargin>3 & i==true & iter>10 & nAccel<100)
      allfull = find(distia(1:iter-2,2)==0);
      j = allfull(end);  %last full step
      if(distia(j,1)>1.05*max(distia(j-1,1),distia(j+1,1)) & iter-j<10*nAccel)
	i = false;
      end;
    end;
  end;