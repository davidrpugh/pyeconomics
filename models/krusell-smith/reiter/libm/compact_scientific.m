function s = compact_scientific(s)
  s = strrep(s,'e+0','e+');
  s = strrep(s,'e+0','e+');
  s = strrep(s,'e-0','e-');
  s = strrep(s,'e-0','e-');
