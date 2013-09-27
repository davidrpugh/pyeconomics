% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% SOLVE FOR NEXT PERIOD'S MOMENTS
% solvelevel: 
%   0: don't solve for_ moments, just use given input xequTry
%   other values: call broydn to find equilibrium moments
function xequFound = solverecu(xequTry,D,iM,ia,ianext, V,solvelevel)
  global Params;
  global broydn_init_B;

  xequFound = xequTry;  %default, used in acceleration step;

  if(solvelevel==0) 
    %don't solve for_ equilibrium, but calling fp_z has side effect: sets Params.DaacTry;
    resid = fp_z(xequTry,D,ianext, V);
    return;
  end
  if(isempty(iM) | mod(Params.t,30)==1)
    broydn_init_B = [];
  else
    broydn_init_B = Params.Bbroydn{iM,ia,ianext};
  end
  [xequFound,check,B] = broydn(@fp_z, xequTry, [1e-7,0,0],D,ianext, V);
  if(~isempty(iM))
    Params.Bbroydn{iM,ia,ianext} = B;
  end
  if(check)
    disp('Warning: broydn NOT converged');
  else
    % disp('Warning: broydn converged!!');    
  end



function z = fp_z(momnext,varargin)
  global Params;
  momnext2 = feval(@reaction,momnext,varargin{:});
  z = momnext2 - momnext; %gives right answer when momnext is at lower bound
  
