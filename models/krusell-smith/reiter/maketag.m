% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% GENERATES FILE NAME DEPENDING ON PARAMETER CONFIGURATION:
function tag = maketag(Params,iStSt,apprType)
  gam = sprintf('g%1.0f',Params.gam);
  if(nargin>1 & iStSt)
    tag = gam;
    return;
  end
  gam = sprintf('%sd%d',gam,length(Params.deltavec));
  if(0)  %NOT USED ANY LONGER
    b0 = 'r';  % 'Regression'
  else
    b0 = 'n';  % 'Neighboring points'
  end;    
  base = sprintf('%s%d',b0,Params.nAssets);

  mom = sprintf('m%d',Params.nMoments);
  suffix = [];
  if(Params.nMoments>1)
    if(Params.useC)
      suffix = [suffix 'C'];
    else
      suffix = [suffix 'M'];  %stands for "moments"
    end;
  end;

  suff2 = sprintf('nn%d_%d_%d',Params.nK,Params.nDistr,Params.n);
  suffix = [suffix suff2];

  if(nargin<3)
    apprType = Params.approxtype;
  end
  switch(apprType)
    case 'lin'
      suffix = [suffix 'L'];
    case 'pwcub'
      suffix = [suffix 'C'];
    case 'spli'
      suffix = [suffix 'B'];  %B-splines
    case 'complpoly'
      suffix = [suffix 'P'];
    otherwise
      error('wrong approximation type');
  end;

  tag = [base gam mom suffix];