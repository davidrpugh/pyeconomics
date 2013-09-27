% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% simulate all models; must first run allsolve
global Params;
initparams;



for uj=1:3 
  Params.whichShocks = uj;

  Params.gam = 1;
  Params.iEulerres = 100;
  Params.iSaveDistr = 0;

  

  n = 1000;
  if(Params.whichShocks==2)
    n = 200;
  end
  if(Params.whichShocks==3)
    n = 400;
  end


  atypes = {'lin','pwcub'};
  aletter = 'LC';
  nmoms = [1 2];
  its = [0 1];
  ds = [0 1];
  
  setrefdist;
  for i=1:length(nmoms);
    Params.nMoments = nmoms(i);
    for j=1:length(atypes);
      Params.approxtype = atypes{j};
      for l=1:length(its);  
	Params.iter=its(l);
	for m=1:length(ds);
	  setparam;
	  
	  Params.useDsel = ds(m);
	  tagSim = sprintf('%s%d%d%d_%d',...
			   aletter(j),Params.whichShocks,Params.nMoments,Params.iter,Params.useDsel);
	  
	  [ser,id] = dosimul(maketag(Params),tagSim,n,0);
	end
      end
    end
  end
  

end
