% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% interpolates function gives by (x,v,dv) at points x0
function [v0,dv0,c0] = interp_grid1d(x0,x,v,dv)
  global Params;

  % init:
  v0 = zeros(size(x0));
  dv0 = v0;
  c0 = zeros(size(x0));
  nx = length(x);

  xcrit = x(1);
  assert(abs(xcrit-invmargutil(dv(1)))<1e-8);
  
  % interpol. or extrapol.:
  ilo = x0 < xcrit;
  ihi = x0 > x(nx);
  iin = ~(ilo | ihi);
  
  % interp:
  [v0(iin),dv0(iin)] = schumaker(x0(iin),x,v,dv);
  %iin2 = iin & x0>20;
  %dv0(iin2) = cubsplinedv(x0(iin2),x,dv);


  if(any(ilo))
    c0low = x0(ilo);
    v0(ilo) = v(1) - util(xcrit) + util(c0low);
    dv0(ilo) = margutil(c0low);
    c0(ilo) = c0low;
  end;
  if(any(ihi))
    cmax = invmargutil(dv(nx-1:nx));
    cslope = (cmax(2)-cmax(1)) / (x(nx)-x(nx-1));
    c0hi = cmax(2) + cslope*(x0(ihi)-x(nx));
    dv0(ihi) = margutil(c0hi);
    v0(ihi) = v(nx) + (util(c0hi) - util(cmax(2)))/cslope;
    c0(ihi) = c0hi;
  end;

  if(nargout>2)
    c0(iin) = invmargutil(dv0(iin));
  end


  %  v(x) = v(x0) + int_x0^x dv/dx
  %       = v(x0) + int_x0^x u'(c(x))
  %       = v(x0) + int_x0^x u'(c(x0)+s(x-x0))
  %       = v(x0) + int_x0^x (c(x0)+s(x-x0))^{-gam}  dx
  %       = v(x0) + int_x0^x (c(x0)+s(x-x0))^{-gam}

