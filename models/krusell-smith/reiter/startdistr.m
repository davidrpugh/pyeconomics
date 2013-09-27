% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% provides distribution to start the simulations
% xAdd: is added to steady state distribution
function D = startdistr(ia)
  global Params;

  if(Params.iSaveDistr)  %used to generate reference distribution
    load(resfilename('stst',maketag(Params,1)),'DKstst');
    % adjust for_ unemployment rate:
    D = DKstst;
    D.p = adjustd(D.p,ia);
  else  %use distribution given by DJJ
    D = loaddistr(Params.U(ia));
  end

