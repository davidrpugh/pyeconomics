% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% gives reference distribution for moms mom and aggregate state ia
function D = finddref(mom,ia,iterUse)
  global Params;
  if(nargin<3 | isempty(iterUse))
    iterUse = Params.iter;
  end
  if(iterUse==0)  % use steady state distribution
    % UPDATED 29.11.2007
    load(resfilename('stst',maketag(Params,1)),'DKstst');
    U = Params.U(ia);
    D = DKstst;
    mD = sum(D.p);
    D.p(:,1) = D.p(:,1) * U/mD(1);
    D.p(:,2) = D.p(:,2) * (1-U)/mD(2);
  elseif(iterUse==1)  % use steady state distribution
    load(Params.RefDistrib,'Dmean');
    D = Dmean{ia};
  else
    assert(0);
  end
  mRef(1) = dst_moment(1,D);
  % WARNING: NOT YET IMPLEMENTED FOR BORROWOING CONSTRAINT ~=0: 
  D.x = D.x*(mom(1)/mRef(1));
