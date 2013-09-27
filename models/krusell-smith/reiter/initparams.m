% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% initializes parameters
global Params;
Params = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  WHICH MODEL?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Params.nAssets = 1;  % number of assets>1 not handled in these files
Params.approxtype = 'lin';
Params.nAccel = 7;   % acceleration steps in backward solution
Params.nStart = 20;  % first iterations no acceleration steps
Params.useC = 1;

% size grid individual k:
Params.n = 500;
% number of intervals in k-distribution
Params.nDistr = 1500;

% size grid aggregate K:
Params.nK = 8;

Params.iSaveDistr = 0;
Params.iMakeForec = 0;
Params.iEulerres = 0;  %DON'T COMPUTE EULER RESIDUALS IN SIMUL (WILL BE CHANGED LATER)

% print and save results each n iterations in solution:
Params.iPrint = 10;
Params.iSave = 10;


% WILL ALL BE CHANGED LATER:
Params.useDsel = 0;
Params.nMoments = 1;
Params.BadStSt = 0;
