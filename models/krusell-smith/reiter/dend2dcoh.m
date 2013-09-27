% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% forward distribution of cash-on-hand to distribution end-of-period capital
%   with exogenous variable ia (this period) and ianext (next period)
function [DCoh,wnext,Rnext] = dbeg2dcoh(Dbeg,ia,ianext)
  global Params;

  K = dst_moment(1,Dbeg);
  L = Params.Lvec(ianext);
  z = Params.zAggr(ianext);
  [Y,wnext,Rnext] = prodfunc(K,L,z);


  DCoh = advancea2x(Dbeg,Params.ProbInd{ia,ianext},Rnext,Params.YInd{ianext}'*wnext);
