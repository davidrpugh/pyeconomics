% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
function V = vbeg_init(VAacStSt)
  global Params;
  % V = grid1d_def_interp(GridK,[Params.nGridMom,Params.n,Params.nAggr,Params.nInd]);

  if(Params.nMoments>=4)
    P2 = Params;
    P2.nMoments = 1;
    P2.approxtype = 'lin';
    rf = resfilename('V',maketag(P2),0);
    load(rf,'Vcoh');
  end

  V.x = VAacStSt.x;
  V.v = zeros([Params.nGridMom,Params.n,Params.nAggr,Params.nInd]);
  V.dv = V.v;
  for iM=1:Params.nGridMom
    for ia=1:Params.nAggr
      for ii=1:Params.nInd
	V.v(iM,:,ia,ii) = VAacStSt.v(:,ii);
	V.dv(iM,:,ia,ii) = VAacStSt.dv(:,ii);
      end
    end;
  end;
