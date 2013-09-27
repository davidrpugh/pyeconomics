% 1-dimensional interpolation, same calling syntax as Matlab's interp1
% does no range check on x0!
% uses lookup from Miranda/Fackler
function [z,dz] = interp1lin(x,Y,x0)
  loc = lookup(x,x0,3);
  xl = x(loc);
  xr = x(loc+1);
  fac = (x0-xl)./(xr-xl);
  nc = size(Y,2);
  if(nc==1)
    z = (1-fac).*Y(loc) + fac.*Y(loc+1);
    if(nargout>1)
      dY = diff(Y)./diff(x);
      dz = dY(loc);
    end
  else
    assert(size(x0,2)==1);
    z = zeros(length(x0),nc);
    dz = z;
    for j=1:nc
      y = Y(:,j);
      z(:,j) = (1-fac).*y(loc) + fac.*y(loc+1);
      if(nargout>1)
	dy = diff(y)./diff(x);
	dz(:,j) = dy(loc);
      end
    end
  end
    