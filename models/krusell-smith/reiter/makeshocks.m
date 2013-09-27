% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% creates shock series, used for simulation to update reference distribution
rand('seed',412347);
shockalld = rand(10000,2);
shocktest = rand(4000,2);
shockalld2 = rand(10000,2);
shocktest2 = rand(4000,2);
save shocks shockalld shocktest shockalld2 shocktest2;