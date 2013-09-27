% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% forward distribution of end-of-period capital to distribution of cash-on-hand
function Dnext = advancea2x(Daac,TransProbInd,Rnext,wnexti)
  global Params;
  Add = wnexti;
  Fac = Rnext*ones(1,Params.nInd);
  Dnext = d_mixtypes(Daac,TransProbInd,Add,Fac);
