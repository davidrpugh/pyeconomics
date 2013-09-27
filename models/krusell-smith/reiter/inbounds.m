% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% projects grid point into the rectangle defined by Params.lowB,Params.uppb
% returns column vector
function m2 = inbounds(m)
  global Params;
  if(size(m,2)==1)  % column vector
    m2 = min(max(m',Params.lowB+1e-8),Params.uppB-1e-8)';
  else
    m2 = min(max(m,Params.lowB+1e-8),Params.uppB-1e-8);
  end