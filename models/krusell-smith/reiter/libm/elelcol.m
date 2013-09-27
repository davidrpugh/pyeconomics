function z = elelcol(x,v)
  if(size(v,1)~=1)
    v = v';
  end;
  z = x .* repmat(v,size(x,1),1);
