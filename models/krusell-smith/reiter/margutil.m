% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% marginal utility
function mu = margutil(c)
  global Params;
  assert(all(c>0),'c not positive');
  mu = c.^(-Params.gam);
