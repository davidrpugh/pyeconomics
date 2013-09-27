% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% forwards end-of-period distribution Daaclag into next period's distribution Daac
% Inputs:
%   VAac:  value function at end of period
%   xequ:  gues of equilibrium next period's moments
%   Daaclag: distribution to start from
%   ialag: last period's exogenous state
%   ia: this period's exogenous state

function [Daac,VCoh,wage,R] = forward(VAac,xequ,Daaclag,ialag,ia)
  global Params;

  Dbeg = d_mixtypes(Daaclag,Params.ProbInd{ialag,ia});
  Dbeg.p = adjustd(Dbeg.p,ia,1);  %adjust each period, other wise the error in U is of order 1e-6;
  [DCoh,wage,R] = dbeg2dcoh(Dbeg,ia);
  xequ = solverecu(xequ,DCoh,[],[],ia,VAac,1);  % solve for_ equilibrium moments
  Daac = Params.DaacTry;  %was set_ in last call to reaction.m;
  VCoh = Params.VcohEndog;

