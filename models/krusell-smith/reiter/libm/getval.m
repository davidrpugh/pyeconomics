function v = getval(x)
  if(isa(x,'deriv1'))
    v = x.v;
  else
    v = x;
  end
