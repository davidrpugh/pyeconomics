% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
function m = stat2mom(s)
  global Params;
  m = s;
  if(Params.nMoments==2)
    if(Params.useC)
      k = m(3);
      if(m(2)>0.07)
	ia = 1;
      else
	ia = 2;
      end
      m(4) = s(4) + [1 k k^2]*Params.CofK{ia};
    else
      m(4) = s(3).^2 + (m(3)*s(4)).^2;  %variance to second noncentral moment;
    end
  end;

  if(Params.nMoments>2)
    error('not yet implemented');
  end
