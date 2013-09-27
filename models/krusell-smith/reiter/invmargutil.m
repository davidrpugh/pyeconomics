% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% inverse marginal utility
function c = invmargutil(mu)
  global Params;
  assert(all(mu>0),'mu not positive');
  c = mu.^(-1/Params.gam);

