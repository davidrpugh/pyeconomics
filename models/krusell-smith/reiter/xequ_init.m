% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
function xEqu = xequ_init
  global Params;
  xEqu = zeros(Params.nGridMom,Params.nEqu,Params.nAggr,Params.nAggr);

  load(resfilename('stst',maketag(Params,1)),'DKstst');
  h = make_hvec(DKstst.x);
  m = mom2stat(h'*DKstst.p(:));
  m = m(Params.iMomK:end);
  if(Params.nMoments>=4)
    m = [m;h2'*DKstst.p(:)];
  end

  for iK=1:Params.nGridMom
    for ia=1:Params.nAggr
      for ianext=1:Params.nAggr
	% initialize with steady state:
	xEqu(iK,:,ia,ianext) = m';
      end
    end
  end;


