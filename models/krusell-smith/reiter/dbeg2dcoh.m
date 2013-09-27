% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% forward distribution of beginning-of-period capital to distribution of cash-on-hand,
%   with exogenous variable ia
function [DCoh,wnext,Rnext] = dbeg2dcoh(Dbeg,ia)
  global Params;

  K = dst_moment(1,Dbeg);
  L = Params.Lvec(ia);
  z = Params.zAggr(ia);
  [Y,wnext,Rnext] = prodfunc(K,L,z);


  DCoh = advancea2x(Dbeg,eye(2),Rnext,Params.YInd{ia}'*wnext);
