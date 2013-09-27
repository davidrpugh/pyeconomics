% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% distribution function_ of D at x;
% separately for_ each column of D;
%   inverse of dst_quantile if_ iNormalize=1;
function p = dst_distrfunc(D,x,iNormalize)

  if(nargin<3)
    iNormalize = 0;
  end

  n = size(D.x,1);
  for i=1:size(D.p,2)
    d = D.p(:,i);
    if(iNormalize)
      d = d/sum(d);
    end
    d = cumsum([0;d]);
    jlos = lookup(D.x(:,i),x,0);

    ilo = x<=D.x(1,i);
    ihi = x>=D.x(n,i);
    iin = ~(ilo|ihi);
    p(ilo,i) = 0;
    p(ihi,i) = d(n-1);
    jlo = jlos(iin);
    fac = (x(iin)-D.x(jlo,i)) ./ (D.x(jlo+1,i)-D.x(jlo,i));
    p(iin,i) = fac.*d(jlo+1) + (1-fac).*d(jlo);

  end

    
  
