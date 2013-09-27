% P is transition matrix: P(i,j) = prob(i->j)
function InvD = invdistr_gen(P)
  if(any(abs(sum(P')-1)>1e-10))
    error('no prob matrix in invdistr_gen');
  end
  n = size(P,1);
  m = eye(n) - P';
  [U,w,V] = svd(m);
  % m = U*w*V';
  w = diag(w);
  InvD = [];
  for i=1:n
    if(abs(w(i))<1e-14)
      D = V(:,i);
      D = D*(1/sum(D));
      InvD = [InvD D];
    end;
  end;

