function inext = changestaterandomly(icurr,randu, transprob)
  assert(randu>=0 & randu<=1, 'invalid random number');
  sum = 0;
  n = size(transprob,2);
  inext=n;
  for(i=1:n-1)
    sum = sum + transprob(icurr,i);  
    if(randu<=sum)
      inext = i;
      return;
    end;
  end;

