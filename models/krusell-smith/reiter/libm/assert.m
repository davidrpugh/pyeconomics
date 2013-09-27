function assert(cond,txt)
  if(~cond)
    if(nargin<2)
      txt = 'assertion wrong';
    end;
    myerror(txt);
  end;