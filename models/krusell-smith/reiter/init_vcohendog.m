% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% initializes value function
function VcohEndog = init_vcohendog
  global Params;
  VcohEndog = cell(Params.nInd,1);
  for i=1:prod(size(VcohEndog));
    grid_coh = grid1d_def_interp(Params.KIndGrid*Params.Rnorm + Params.YIndStst(i)*Params.wnorm); 
    VcohEndog{i} = grid_coh;
  end;
