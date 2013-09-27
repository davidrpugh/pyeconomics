% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% solves dynamic programming problem of the household in steady state with capital K
function [Vcoh,VAac2] = dynpro(K)
  global Params;
  
  setparam;

  [dummy,Params.wstst,Params.Rstst] = prodfunc(K,Params.Lstst,1);

  Vcoh = cell(Params.nInd,1);
  VAac = cell(Params.nInd,1);
  for i=1:Params.nInd
    % x given only at start, later endogenous:
    Vcoh{i} = grid1d_def_interp(Params.KIndGrid*Params.Rstst + Params.YIndStst(i)*Params.wstst); 
    Vcoh{i}.v = util(Vcoh{i}.x);
    Vcoh{i}.dv = margutil(Vcoh{i}.x);
    VAac{i} = grid1d_def(Params.KIndGrid);
  end;

  for iter=1:2000  %set iter to 1000, don't stop early, to allow differentiation
    for id=1:Params.nInd
      % grid for Vcoh is after interest, but before labor income:
      VAac{id}.v = 0;
      VAac{id}.dv = 0;
      for idnext=1:Params.nInd
	wage_i = Params.wstst*Params.YIndStst(idnext);
	CoH = VAac{id}.x *Params.Rstst + wage_i;
	prob = Params.TransProbStst(id,idnext);
	[v,dv] = interp_grid1d(CoH,Vcoh{idnext}.x,Vcoh{idnext}.v,Vcoh{idnext}.dv);
	VAac{id}.v  = VAac{id}.v  + (Params.betta*prob)*v;
	VAac{id}.dv = VAac{id}.dv + (Params.betta*prob*Params.Rstst)*dv;
      end;
    end;
    Vold = Vcoh;
    dist = 0;
    for id=1:Params.nInd
      [Vcoh{id}.x,Vcoh{id}.v, Vcoh{id}.dv] = v_optimc(VAac{id}.x,VAac{id}.v,VAac{id}.dv);
      dist = max(dist,norm(Vcoh{id}.v-Vold{id}.v));
    end;
    % if(mod(iter,100)==0), disp(sprintf('iter %d: dist = %e',iter,dist)); end;
  end;

  % Change the format of VAac:
  VAac2.x = VAac{1}.x;  %same x for VAac{1} and VAac{2};
  VAac2.v = [VAac{1}.v VAac{2}.v];
  VAac2.dv = [VAac{1}.dv VAac{2}.dv];
