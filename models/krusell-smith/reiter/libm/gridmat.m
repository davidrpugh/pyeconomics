function [mat,nn] = gridmat(x)
  assert(iscell(x));
  nc = length(x);
  nr = 0;
  for i=1:nc
    nr = max(nr,length(x{i}));
  end;
  mat = NaN*ones(nr,nc); 
  for i=1:nc
    n = length(x{i});
    nn(i) = n;
    mat(1:n,i) = x{i};
  end;
  
  