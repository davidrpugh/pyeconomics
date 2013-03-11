function [residual, g1, g2] = diamond_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 11, 1);

%
% Model equations
%

T52 = params(2)^((-1)/params(5));
T58 = 1+T52*(1+y(6)-params(3))^((-(1-params(5)))/params(5));
lhs =y(5);
rhs =y(2)+y(7);
residual(1)= lhs-rhs;
lhs =y(1);
rhs =y(2)+1/((1+params(6))*(1+params(7)))*y(3);
residual(2)= lhs-rhs;
lhs =y(9);
rhs =y(1)+y(8);
residual(3)= lhs-rhs;
lhs =y(3);
rhs =y(7)*(1+y(6)-params(3));
residual(4)= lhs-rhs;
lhs =y(4);
rhs =y(8)+y(4)*(1-params(3))/((1+params(6))*(1+params(7)));
residual(5)= lhs-rhs;
lhs =y(9);
rhs =y(4)^params(1);
residual(6)= lhs-rhs;
lhs =y(5);
rhs =y(9)*(1-params(1));
residual(7)= lhs-rhs;
lhs =y(6);
rhs =y(9)*params(1)/y(4);
residual(8)= lhs-rhs;
lhs =y(7);
rhs =y(5)/T58;
residual(9)= lhs-rhs;
lhs =y(10);
rhs =y(9)-(y(4)-y(4)*(1-params(3))/((1+params(6))*(1+params(7))))-y(1);
residual(10)= lhs-rhs;
lhs =y(11);
rhs =y(3)-y(4)*(1+params(6))*(1+params(7))*(1+y(6)-params(3));
residual(11)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(11, 11);

%
% Jacobian matrix
%

  g1(1,2)=(-1);
  g1(1,5)=1;
  g1(1,7)=(-1);
  g1(2,1)=1;
  g1(2,2)=(-1);
  g1(2,3)=(-(1/((1+params(6))*(1+params(7)))));
  g1(3,1)=(-1);
  g1(3,8)=(-1);
  g1(3,9)=1;
  g1(4,3)=1;
  g1(4,6)=(-y(7));
  g1(4,7)=(-(1+y(6)-params(3)));
  g1(5,4)=1-(1-params(3))/((1+params(6))*(1+params(7)));
  g1(5,8)=(-1);
  g1(6,4)=(-(getPowerDeriv(y(4),params(1),1)));
  g1(6,9)=1;
  g1(7,5)=1;
  g1(7,9)=(-(1-params(1)));
  g1(8,4)=(-((-(y(9)*params(1)))/(y(4)*y(4))));
  g1(8,6)=1;
  g1(8,9)=(-(params(1)/y(4)));
  g1(9,5)=(-(1/T58));
  g1(9,6)=(-((-(y(5)*T52*getPowerDeriv(1+y(6)-params(3),(-(1-params(5)))/params(5),1)))/(T58*T58)));
  g1(9,7)=1;
  g1(10,1)=1;
  g1(10,4)=1-(1-params(3))/((1+params(6))*(1+params(7)));
  g1(10,9)=(-1);
  g1(10,10)=1;
  g1(11,3)=(-1);
  g1(11,4)=(1+params(6))*(1+params(7))*(1+y(6)-params(3));
  g1(11,6)=(1+params(6))*(1+params(7))*y(4);
  g1(11,11)=1;
  if ~isreal(g1)
    g1 = real(g1)+imag(g1).^2;
  end
end
if nargout >= 3,
%
% Hessian matrix
%

  g2 = sparse([],[],[],11,121);
end
end
