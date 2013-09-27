% a version where wealth comes first is in directory svi2!!
function [ic,notc] = v_concave(V)
  notc = [];
  ic = 1;
  n = size(V.v);
  for i=1:n(1)
    for j=1:n(3)
      for k=1:n(4)
	v = V.v(i,:,j,k);
	dv = V.dv(i,:,j,k);
	isc = is_concave(V.x,v(:),dv(:));
	if(max(diff(dv(:)))>=0)
	  error('in v_con');
	end;
	if(isc==0)
	  ic=0;
	  notc = [notc [i;j;k]];
	end;
      end;
    end;
  end