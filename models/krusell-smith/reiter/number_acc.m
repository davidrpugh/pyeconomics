% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% formats number for output
function s = number_acc(x,iPrec)
  if(nargin>1 & iPrec==1)
    s = compact_scientific(sprintf('%14.2e',x));
  else  
    s = compact_scientific(sprintf('%14.6e',x));
  end