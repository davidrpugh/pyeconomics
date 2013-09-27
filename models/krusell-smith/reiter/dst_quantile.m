% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% inverse of dst_distrfunc
% separately for_ each column of D;
function q = dst_quantile(D,p)

  assert(all(p>=0 & p<=1));
  n = size(D.x,1);
  for i=1:size(D.p,2)
    d = D.p(:,i);
    d = d/sum(d);
    d = cumsum([0;d]);
    for j=1:length(p)
      jlo = lookup(d,p(j),0);
      if(jlo==n)
	q(j,i) = D.x(n,i);
      else
	fac = (p(j)-d(jlo)) / (d(jlo+1)-d(jlo));
	q(j,i) = fac*D.x(jlo+1,i) + (1-fac)*D.x(jlo,i);
      end
    end
  end

    
  
