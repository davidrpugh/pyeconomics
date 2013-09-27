% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
function outp = transmatx2a(gridCoH,gridA,VaacInt,D)
  [n,m] = size(gridCoH);
  [n2,m2] = size(gridA);
  assert(m==m2);

  iD = nargin>3;  %forward dstribution, otherwise_ give transition matrix

  if(iD)
    D2.x = gridA;
    D2.p = zeros(n2-1,2);
  end
  for j=1:m
    a = gridA(:,j);
    [v_unused,dv] = schumaker(a,VaacInt.x,VaacInt.v(:,j),VaacInt.dv(:,j));
    %iBig = a>20;
    %dv(iBig) = cubsplinedv(a(iBig),VaacInt.x,VaacInt.dv(:,j));
    c = invmargutil(dv);
    x = a+c;  %CoH
    X = gridCoH(:,j);
    x = min(max(x,X(1)),X(n));
    pos = lookup(X,x,3);
    frac = (x-X(pos))./(X(pos+1)-X(pos));
    % tails of D will be matched into first and last element of D2:
    pos(1) = 1;
    frac(1) = 0;
    pos(n2) = n-1;
    frac(n2) = 1;

    if(iD)
      D2.p(:,j) = forwardx2a(pos,frac,D.p(:,j));
    else
      [iFr{j},iTo{j},Val{j}] = forwmatx2a(pos,frac);
    end
  end

  if(iD)
    outp = D2;
  else
    n21 = n2-1;
    nn = 2*n21;
    outp = sparse([iTo{1};iTo{2}+n21],[iFr{1};iFr{2}+n21],[Val{1};Val{2}],nn,nn);
  end


