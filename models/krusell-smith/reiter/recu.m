% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% GRIDS:
% VAac: cash on hand, after production, including wages
% Vaac: assets after consumption, lower bound =0
%
% VAac and xEqu are guess on input;
%   if no argument given or empty, defaults are used
function [VAac,xEqu] = recu(VAac,xEqu)
  global Params;  
  tic;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % INITIALIZATIONS:
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % load steady state results, including distribution;
  loadstst;
  
  Params.nReaction = 0;
  Params.nCallMom = 0;
  Params.cpu = cputime;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % CREATE DISTRIBUTIONS:
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % create distribution at each point in the state space:
  Params.AllDistribs = cell(Params.nGridMom,Params.nAggr);

  m = [1;0;0];
  for ia=1:Params.nAggr
    for iK=1:Params.nGridMom
      fprintf(1,'%d  %d\n',ia,iK);
      s = [1;Params.U(ia);Params.GridMom(iK,:)'];
      D = dsf_s(s,ia);
      Params.AllDistribs{iK,ia} = D;
    end;
  end;
  disp('dists done');

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % INITIALIZE GRIDS:
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if(nargin==0 | isempty(VAac))
    VAac = vbeg_init(VAacStSt);
    VAac.approxmom = Params.approxmom;    
    [isc,notc] = v_concave(VAac);
    assert(isc);
  end;


  % initialize equilibrium outcome:
  if(nargin==2 & ~isempty(xEqu))
    Params.xequ_given=1;
  else
    Params.xequ_given=0;
    xEqu = xequ_init;
  end;



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PREPARE ITERATIONS:
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Vneu = VAac;
  vtol = 1e-6;
  % maximum number of iterations in time:
  maxiter = 2000;

  % information on acceleration steps:
  distia = zeros(maxiter,2);
  % initialize distance in k-space:
  distk = 0;
  % store information on Jacobian for_ root finding:
  Params.Bbroydn = cell(Params.nGridMom,Params.nAggr);

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % DO ITERATIONS:
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for iter=1:maxiter
    Params.t = iter;  %used in solverecu

    % find out whether this step is an acceleration step, and store the information:
    isacc = is_accel(iter,Params.nAccel,Params.nStart,distia);
    distia(iter,2) = isacc;

    
    VAac = Vneu;
    Vneu.approxmom = Params.approxmom;    
    VAac = prepare_mom(VAac,Params.approxmom);
    [isc,notc] = v_concave(VAac);
    assert(isc);

    dist = 1;  % distance between new and old value function_, will be measured from iteration 2 on
    xEquOld = xEqu;

    if(isacc)
      solvelevel = 0;
    else
      solvelevel = 1;
    end;
    Vneu.v = zeros(size(VAac.v)); Vneu.dv = Vneu.v;
    for iK=1:Params.nGridMom
      for ia=1:Params.nAggr
	D = Params.AllDistribs{iK,ia};
	for ianext=1:Params.nAggr
	  % compute distribution of cash-on-hand
	  [DCoh,wnext,Rnext] = dend2dcoh(D,ia,ianext);
	  % solve for equilibrium and update value functions at CoH:
	  xequTry = xEqu(iK,:,ia,ianext);
	  xequTry2 = solverecu(xequTry(:),DCoh,iK,ia,ianext,VAac,solvelevel);
	  xEqu(iK,:,ia,ianext) = xequTry2';
	  % update value function at end-of-period:
	  for ii2=1:Params.nInd  %next period's individual state
	    jnext = (ianext-1)*2 + ii2;
	    x = VAac.x*Rnext + Params.YInd{ianext}(ii2)*wnext;
	    [v,dv] = interp_grid1d(x,Params.VcohEndog{ii2}.x,Params.VcohEndog{ii2}.v,Params.VcohEndog{ii2}.dv);
	    if(max(diff(dv))>=0)
	      error('not concave here')
	    end;
	    for ii=1:Params.nInd
	      j = (ia-1)*2 + ii;
	      fac = Params.TransProb(j,jnext) * Params.betta;
	      Vneu.v(iK,:,ia,ii) = Vneu.v(iK,:,ia,ii) + fac*v';
	      Vneu.dv(iK,:,ia,ii) = Vneu.dv(iK,:,ia,ii) + (fac*Rnext)*dv';
	    end
	  end;
	end;
      end;
    end;
    if(iter>1)  %maximum difference in marginal utility:
      dist = max(abs(VAac.dv(:)./Vneu.dv(:) - 1));
    end
    distia(iter,1) = dist;
    if(~isacc)
      distk = max(abs(xEquOld(:)-xEqu(:)));
      % disp(sprintf('iter %d: dist policy = %e',iter,distk));
    end;
    if(Params.iPrint>0 & mod(iter,Params.iPrint)==0)
      disp(sprintf('iter %d: dist = %e; distpol = %e;',iter,dist,distk));
      toc;tic;
    end;
    if(~isacc)
      if(dist<vtol & iter>2)
	disp(sprintf('done after %d iter with dist=%e',iter,dist));
	break;
      end;
    end;
    Params.xequ_given = 1;
    if(Params.iSave>0 & mod(iter,Params.iSave)==0)
      Params2 = Params;
      Params2.AllDistribs = [];
      save recures VAac xEqu Params2;
    end;
  end;
  fn = resfilename('V',maketag(Params),Params.iter);
  Params2 = Params;
  Params2.AllDistribs = [];
  save(fn,'VAac','xEqu','Params2');



