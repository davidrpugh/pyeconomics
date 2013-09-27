function z = is_concave(x,f,df)
  tol = 1e-10;
  z = 1;
  n = length(x);
  df1 = df(1:n-1);
  df2 = df(2:n);
  if(~all(df2-df1<=0))  %also catches NaN
    z = 0;
    return;
  end;
  step = diff(x);

  slope = diff(f) ./ step;
  if(~all( (df1-slope)>=-tol & (df2-slope<=tol)))   %also catches NaN
    z = 0;
  end;