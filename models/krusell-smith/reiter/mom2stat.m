% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% transforms moments into the statistics on which the grid is defined
function s = mom2stat(m,nMom)
  global Params;
  if(nargin<2)
    nMom = Params.nMoments;
  end
  s = m;
  if(nMom==2)
    if(Params.useC)
      k = m(3);
      if(m(2)>0.07)
	ia = 1;
      else
	ia = 2;
      end
      s(4) = m(4) - [1 k k^2]*Params.CofK{ia};
    else
      s(4) = sqrt(m(4) - m(3).^2) / m(3) ;  %second noncentral moment to rel. stdev.
    end
  end;

  if(Params.nMoments>2)
    error('not yet implemented');
  end
