function V = prepare_mom(V,fapprox)
  global Params;
  if(strcmp(fapprox.type,'lin') | strcmp(fapprox.type,'pwcub'))
    % do nothing;
  elseif(strcmp(fapprox.type,'spli'))
    nn = size(V.v);
    assert(length(nn)==4);
    nn2 = [nn(1) prod(nn(2:4))];
    vv = reshape(V.v,nn2);
    dvv = reshape(V.dv,nn2);
    V.coeff_v = fapprox.B\vv;
    V.coeff_dv = fapprox.B\dvv;
  else
    error('not yet impl');
  end;
