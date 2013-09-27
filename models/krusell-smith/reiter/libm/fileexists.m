function ex = fileexists(fname)
  ex = 0;  %default: does not exist;
  fid = fopen(fname);
  if(fid>=0)  % could open for reading
    ex = 1;
    fclose(fid);
  end

%  try
%    kk = fileattrib(fname);
%  catch
%    ex = 0;  %if error: doesn't exist;
%  end;