% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% loads steady state solution for later use
  load(resfilename('stst',maketag(Params,1)),'DXstst','DKstst','Pstst','VCohStSt','VAacStSt');
  % check we have the right stst file:
  assert(Params.gam==Pstst.gam & Params.mu==Pstst.mu);
  Params.nGridAggr = Params.nGridMom*Params.nAggr;
  % will be used in momstart for starting value in simulation:
  Params.Kstst = dst_moment(1,DKstst);
  Params.Xstst = dst_moment(1,DXstst);
  disp(sprintf('Kstst RA: %f; HA: %f; Coh: %f',Params.Knorm,Params.Kstst,Params.Xstst));
