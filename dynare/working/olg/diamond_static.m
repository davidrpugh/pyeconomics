function [residual, g1, g2] = diamond_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 10, 1);

%
% Model equations
%

T19 = 1/((1+params(6))*(1+params(7)));
T55 = params(2)^((-1)/params(5));
T60 = 1+T55*(1+y(6)-params(3))^((params(5)-1)/params(5));
lhs =y(8);
rhs =y(1)+y(9);
residual(1)= lhs-rhs;
lhs =y(1);
rhs =y(2)+T19*y(3);
residual(2)= lhs-rhs;
lhs =y(8);
rhs =y(5)+y(6)*y(4);
residual(3)= lhs-rhs;
lhs =y(4);
rhs =T19*(y(9)+y(4)*(1-params(3)));
residual(4)= lhs-rhs;
lhs =y(6);
rhs =params(1)*y(4)^(params(1)-1);
residual(5)= lhs-rhs;
lhs =y(5);
rhs =y(8)*(1-params(1));
residual(6)= lhs-rhs;
lhs =y(2);
rhs =y(5)-y(7);
residual(7)= lhs-rhs;
lhs =y(3);
rhs =y(4)*(1+params(6))*(1+params(7))*(1+y(6)-params(3));
residual(8)= lhs-rhs;
lhs =y(7);
rhs =y(5)/T60;
residual(9)= lhs-rhs;
lhs =y(10);
rhs =y(4)-T19*y(7);
residual(10)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(10, 10);

%
% Jacobian matrix
%

  g1(1,1)=(-1);
  g1(1,8)=1;
  g1(1,9)=(-1);
  g1(2,1)=1;
  g1(2,2)=(-1);
  g1(2,3)=(-T19);
  g1(3,4)=(-y(6));
  g1(3,5)=(-1);
  g1(3,6)=(-y(4));
  g1(3,8)=1;
  g1(4,4)=1-T19*(1-params(3));
  g1(4,9)=(-T19);
  g1(5,4)=(-(params(1)*getPowerDeriv(y(4),params(1)-1,1)));
  g1(5,6)=1;
  g1(6,5)=1;
  g1(6,8)=(-(1-params(1)));
  g1(7,2)=1;
  g1(7,5)=(-1);
  g1(7,7)=1;
  g1(8,3)=1;
  g1(8,4)=(-((1+params(6))*(1+params(7))*(1+y(6)-params(3))));
  g1(8,6)=(-((1+params(6))*(1+params(7))*y(4)));
  g1(9,5)=(-(1/T60));
  g1(9,6)=(-((-(y(5)*T55*getPowerDeriv(1+y(6)-params(3),(params(5)-1)/params(5),1)))/(T60*T60)));
  g1(9,7)=1;
  g1(10,4)=(-1);
  g1(10,7)=T19;
  g1(10,10)=1;
  if ~isreal(g1)
    g1 = real(g1)+imag(g1).^2;
  end
end
if nargout >= 3,
%
% Hessian matrix
%

  g2 = sparse([],[],[],10,100);
end
end
