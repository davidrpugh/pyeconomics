function [m,mi,m_cond] = dst_moment(k,D,iCentral)
  if(nargin>2 & iCentral & k>1)
    meank = dst_moment(1,D);
    Db = D;
    Db.x = D.x - meank;
    m = dst_moment(k,Db);
  else
    H = dst_h_moment(k,D.x);
    m = H'*D.p(:);
    if(nargout>1)
      for i=1:size(D.p,2)
	h = dst_h_moment(k,D.x(:,i));
	mi(i) = h'*D.p(:,i);
	m_cond(i) = mi(i) / sum(D.p(:,i));
      end
    end
  end

