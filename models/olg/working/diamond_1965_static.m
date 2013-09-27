function [residual, g1, g2] = diamond_1965_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 10, 1);

%
% Model equations
%

T14 = 1/((1+params(7))*(1+params(6)));
T47 = params(2)^((-1)/params(5));
T52 = 1+T47*(1+y(2))^((params(5)-1)/params(5));
lhs =y(1);
rhs =T14^params(1)*y(8)^params(1);
residual(1)= lhs-rhs;
lhs =y(2);
rhs =params(1)*T14^(params(1)-1)*y(8)^(params(1)-1)-params(3);
residual(2)= lhs-rhs;
lhs =y(3);
rhs =y(1)*(1-params(1));
residual(3)= lhs-rhs;
lhs =y(5)+y(7);
rhs =y(3);
residual(4)= lhs-rhs;
lhs =y(6);
rhs =y(7)*1/(1+params(7))*(1+y(2));
residual(5)= lhs-rhs;
lhs =y(7);
rhs =y(3)/T52;
residual(6)= lhs-rhs;
lhs =y(4);
rhs =y(5)+y(6)*1/(1+params(6));
residual(7)= lhs-rhs;
lhs =y(8);
rhs =y(7)+y(8)*(1-params(3))/((1+params(7))*(1+params(6)));
residual(8)= lhs-rhs;
lhs =y(9);
rhs =y(1)-y(8)*T14*(y(2)+params(3))-y(3);
residual(9)= lhs-rhs;
lhs =y(10);
rhs =y(5)^(-params(5))-(1+y(2))*params(2)*(1+params(7))^(-params(5))*y(6)^(-params(5));
residual(10)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(10, 10);

%
% Jacobian matrix
%

  g1(1,1)=1;
  g1(1,8)=(-(T14^params(1)*getPowerDeriv(y(8),params(1),1)));
  g1(2,2)=1;
  g1(2,8)=(-(params(1)*T14^(params(1)-1)*getPowerDeriv(y(8),params(1)-1,1)));
  g1(3,1)=(-(1-params(1)));
  g1(3,3)=1;
  g1(4,3)=(-1);
  g1(4,5)=1;
  g1(4,7)=1;
  g1(5,2)=(-(y(7)*1/(1+params(7))));
  g1(5,6)=1;
  g1(5,7)=(-(1/(1+params(7))*(1+y(2))));
  g1(6,2)=(-((-(y(3)*T47*getPowerDeriv(1+y(2),(params(5)-1)/params(5),1)))/(T52*T52)));
  g1(6,3)=(-(1/T52));
  g1(6,7)=1;
  g1(7,4)=1;
  g1(7,5)=(-1);
  g1(7,6)=(-(1/(1+params(6))));
  g1(8,7)=(-1);
  g1(8,8)=1-(1-params(3))/((1+params(7))*(1+params(6)));
  g1(9,1)=(-1);
  g1(9,2)=T14*y(8);
  g1(9,3)=1;
  g1(9,8)=T14*(y(2)+params(3));
  g1(9,9)=1;
  g1(10,2)=params(2)*(1+params(7))^(-params(5))*y(6)^(-params(5));
  g1(10,5)=(-(getPowerDeriv(y(5),(-params(5)),1)));
  g1(10,6)=(1+y(2))*params(2)*(1+params(7))^(-params(5))*getPowerDeriv(y(6),(-params(5)),1);
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
