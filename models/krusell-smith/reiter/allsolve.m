% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% computes solution for all model variants
global Params;
% initialize parameters:
initparams;

% risk aversion parameter:
Params.gam = 1;

% also solve for steady state? (if not, load earlier results)
doStst = 0;
% do first round of results?
doFirst = 1;


% FIND STEADY STATE:
Params.useC = 0;
if(doStst)
  % compute steady state for_ comparison, with high unemployment rate:
  Params.BadStSt = 1;
  setparam;
  findkstst;

  % compute steady state with average unemployment rate, used for_ reference distribution:
  Params.BadStSt = 0;
  setparam;
  tcpu = cputime;
  findkstst;
  fprintf(1,'cpu time for_ steady state: %0.2f seconds\n',cputime-tcpu);
end
% use consumption vector as second moment, if_ more_ than one moment used:
Params.useC = 1;

% FIRST SOLUTION:
% piecewise linear approximation in aggregate variables:
Params.approxtype = 'lin';
% use 1 moment:
Params.nMoments = 1;
% use 7 acceleration steps in solution:
Params.nAccel = 7;
% use steady state for_ reference distribution:
Params.iter = 0;
if(doFirst)
  tcpu = cputime;
  doone;
  fprintf(1,'cpu time for_ solution: %0.2f seconds\n',cputime-tcpu);

  % simulate and save distributions:
  Params.iSaveDistr=1;
  Params.whichShocks = 0;  %draw random numbers
  Params.useDsel = 0;
  n = 7000;

  % write file which_ indicates where reference distribution is stored:
  FF = fopen('setrefdist.m','w');
  fprintf(FF,'Params.RefDistrib = ''%s'';\n',resfilename('D',maketag(Params,0),1));
  fclose(FF);
  % simulate and store the resulting distributions:
  [ser1_0,id1] = dosimul(maketag(Params),'Ref',n,0);
  % analyze the distribution and store the results:
  setrefdist;
  analdistr;
end
setparam;
setrefdist;


atypes = {'lin','pwcub'};
nmoms = [1 2];
its = [0 1];

for i=1:length(nmoms);
  Params.nMoments = nmoms(i);
  for j=1:length(atypes);
    Params.approxtype = atypes{j};
    for l=1:length(its);
      Params.iter=its(l);

      tcpu = cputime;
      doone;
      fprintf(1,'cpu time for_ solution: %0.2f seconds\n',cputime-tcpu);
    end
  end
end

