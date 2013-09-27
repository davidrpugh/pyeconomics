% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% generates structure of global parameters Params;
function setparam;
  global Params;
  global broydn_step_jacob;

  % set 1 to check consistency of stochastic solution with StSt
  Params.checkStSt=0;


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  PARAMETERS OF SOLUTION METHOD:
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % number of moments to be used:
  % assert(isfield(Params,'nMoments'));
  % spline or piecewiese linear:
  Params.ispline = 0;
  % assert(isfield(Params,'approxtype'));
  % number of acceleration steps:
  % assert(isfield(Params,'nAccel'));
  % initial number of iter in which no accel. steps are taken:
  % assert(isfield(Params,'nStart'));
  % how to determine reference distribution (for recursive method only)

  % aggregate productivity states:
  Params.nAggr = 2;
  % size grid 2nd moment (in case it applies):
  Params.nMom2 = 8;

  % persistent individual income states;
  Params.nInd = 2;  

  % make conditional forecasts up to 40 periods:
  Params.nPred = 0;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  MODEL PARAMETERS:
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  Params.betta = 0.99; % disount factor
  
  Params.mu = 0.15;  % income if_ unemployed
  
  Delta_a = 0.01;  % aggregate productivity fluctuations 
  Params.deltavec = 0;  %depreciation in excess of delta;
  ndv = length(Params.deltavec);
  Params.deltaprob = ones(ndv,1)/ndv;  %depreciation in excess of delta;
 
  % R = 1 + Params.alpha*A*(K/L)^(Params.alpha-1) - Params.delta;
  Params.alpha = 0.36;  % production function
  Params.delta = 0.025; % deprecitation, quarterly
  % TransProb(i,j) is prob. to go from i to j:
  Params.TransProb = [  0.525 0.35 0.03125 0.09375;
		      0.038889 0.836111 0.002083 0.122917;
		      0.09375 0.03125 0.291667 0.583333;
		      0.009115 0.115885 0.024306 0.850694];


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % DERIVED PARAMETERS:

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % invariant distribution of aggregate productivity:
  Params.InvD = invdistr_gen(Params.TransProb);
  % Steady state unemployment:
  Params.Ustst = Params.InvD(1)+Params.InvD(3);
  % Transition matrix steady state:
  nA = Params.nAggr;
  nI = Params.nInd;
  % Stst transition matrix of individual states:
  np = size(Params.TransProb,1);
  offs = 0:nI:np-1;
  for i=1:nI
    for j=1:nI
      pst(i,j) = sum(sum(Params.TransProb(i+offs,j+offs)))/nA;
    end;
  end;
  
  Params.TransProbBadStst = [0.6 0.4;0.044445 0.955555];
  assert(all(abs(sum(Params.TransProbBadStst')-1)<1e-14));
  % Steady state unemployment:
  invdBad = invdistr_gen(Params.TransProbBadStst);
  Params.UststBad = invdBad(1);
  Params.lbar = 1/(1-Params.UststBad);  %normalization, DJJ
  if(Params.BadStSt)  %the case_ required by DJJ
    if(0) %take probabilities from Table 2
      Params.TransProbStst = Params.TransProb(1:2,1:2);
      Params.TransProbStst = Params.TransProbStst ./ repmat(sum(Params.TransProbStst,2),1,2);
    else  %take probabilities from Table 3
      Params.TransProbStst = Params.TransProbBadStst;
    end
  else % the case_ with mean unemployment, used for_ reference distribution:
    p11 = pst(1,1);
    % find transition matrix that generates right unemployment rate:
    p22 = fzero(@findpstst_z,[.001 .999],[],p11);
    Params.TransProbStst = [p11 1-p11;1-p22 p22];
  end
  if(Params.checkStSt==1)
    Params.TransProb = kron(0.5*ones(2,2),Params.TransProbStst);
    Delta_a = 0.0;
  end


  % REPEAT, in case_that TransProb was modified:
  Params.InvD = invdistr_gen(Params.TransProb);
  % Steady state unemployment:
  if(Params.BadStSt)  %the case_ required by DJJ
    Params.Ustst = Params.UststBad;
  else
    Params.Ustst = Params.InvD(1)+Params.InvD(3);
  end
  % Steady state labor input:
  Params.Lstst = Params.lbar*(1-Params.Ustst);
  Params.Emplvec = zeros(Params.nAggr,1);
  Params.Lvec = zeros(Params.nAggr,1);
  for i=1:Params.nAggr
    Params.Emplvec(i) = 1 - Params.InvD(2*i-1)/sum(sum(Params.InvD(2*i-1:2*i)));
    Params.Lvec(i) = Params.Emplvec(i)*Params.lbar;
  end;
  Params.Rnorm = 1/Params.betta;
  % R = 1 + alpha*A*(K/L)^(alpha-1) - delta;
  Params.A = 1;  %following Equ.(2) of HJJ
  KL = ((Params.Rnorm + Params.delta - 1)/Params.alpha)^(1/(Params.alpha-1));
  Params.Knorm = Params.Lstst*KL;
  [dummy,Params.wnorm,R] = prodfunc(Params.Knorm,Params.Lstst,1);
  % only to check StSt calculations:
  assert(abs(R-Params.Rnorm)<1e-8);
  if(Params.BadStSt==0)
    assert(abs(0.07-Params.Ustst)<1e-5);
  end
  if(Params.nAssets ==2)
    broydn_step_jacob = 1e-7;
  else
    broydn_step_jacob = 1e-6;
  end;

  Params.ProbInd = cell(nA,nA);
  Params.ProbAggr = zeros(nA,nA);
  for i=1:nA
    offsi = (i-1)*nI;
    for j=1:nA
      offsj = (j-1)*nI;
      Params.ProbInd{i,j} = Params.TransProb(offsi+(1:nI),offsj+(1:nI));
      pa = mean(sum(Params.ProbInd{i,j}'));
      Params.ProbAggr(i,j) = pa;
      Params.ProbInd{i,j} = Params.ProbInd{i,j}/pa;
    end;
  end;
  Params.U = zeros(nA,1);
  % Notice: slightly deviates from 0.04 and 0.1!!!
  for i=1:nA
    Params.U(i) = Params.InvD(i*nI-1) / sum(Params.InvD(i*nI-1:i*nI));
    % vector of income realizations:
    Params.YInd{i} = make_yind(Params.U(i));
  end;
  [Params.YIndStst,Params.taustst] = make_yind(Params.Ustst);
  Params.zAggr = [1-Delta_a;1+Delta_a;];


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Define grids and approximations:
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Params.lowB = 0.92*Params.Knorm;
  Params.uppB = 1.08*Params.Knorm;
  if(Params.checkStSt==1)
    load(resfilename('stst',maketag(Params,1)),'DXstst');
    xstst = dst_moment(1,DXstst);
    Params.lowB = 0.9999*xstst;    
    Params.uppB = 1.0001*xstst;
  end
  % Grid aggregate K:
  Params.KAggrGrid = linspace(Params.lowB(1),Params.uppB(1),Params.nK)';
  %   take logs of xEqu?



  % Construct nodes for individual k (value function and distribution)
  % grid will be equidistant in log(k+kshift):
  kmin = 0;  %liquidity constraint:
  kmax = 10*Params.Knorm;
  kshift = 0.1*Params.wnorm;  %shift by 1 quarter of wages:
  Params.KIndGrid = logspaceshift(kmin,kmax*1.001,Params.n,kshift)';

  % compute nodes for cross-sectional distribution:
  assert(kmin==0);
  node1 = (0:0.1:99.99)';
  % node2 = linspace(100,kmax*1.01,Params.nDistr-length(node1))';
  node2 = linspace(100,kmax,Params.nDistr-length(node1))';
  Params.DNodes = [0;1e-8;node1(2:end);node2];

  Params.nEqu = Params.nMoments;  % 1 equ. variable to solve for at each grid point;


  if(Params.nMoments==1)
    % for_ regression of xEqu on Moments,
    %   take logs of Moments?
    Params.inlogs = 1;
    Params.GridMom = Params.KAggrGrid;
    [Params.GridMat,Params.nGrid] = gridmat({Params.KAggrGrid});
    Params.approxmom = momdef(Params.KAggrGrid,Params.approxtype);
  elseif(Params.nMoments==2)
    Params.inlogs = [1;0];
    Params.lowB(2) = -0.004;
    Params.uppB(2) =  0.0025;
    Params.Mom2Grid = linspace(Params.lowB(2),Params.uppB(2),Params.nMom2)';
    Params.GridMom = gridmake(Params.KAggrGrid,Params.Mom2Grid);
    [Params.GridMat,Params.nGrid] = gridmat({Params.KAggrGrid,Params.Mom2Grid});
    Params.approxmom = momdef({Params.KAggrGrid,Params.Mom2Grid},Params.approxtype);
  end;
  % prepare for_ case_ useC:
  if(Params.useC)
    load(resfilename('stst',maketag(Params,1)),'DKstst','Cstst','Caggrstst');
    Params.Cstst = Cstst;
    h = make_hvec(repmat(Params.DNodes,1,2),2);
    Params.Caggrstst = h(:,4)'*DKstst.p(:);
  end

  Params.nGridMom = size(Params.GridMom,1);

  if(isfield(Params,'RefDistrib'))
    load(Params.RefDistrib,'CofK');
    Params.CofK = CofK;
  end

  Params.VcohEndog = init_vcohendog;


function z = findpstst_z(p22,p11)
  global Params;
  p = [p11 1-p11;1-p22 p22];
  InvD = invdistr_gen(p);
  z = InvD(1) - Params.Ustst;
  
function [Y,tau] = make_yind(U)
  global Params;
  tau = Params.mu/Params.lbar * U/(1-U);
  Y = [Params.mu;Params.lbar*(1-tau)];
