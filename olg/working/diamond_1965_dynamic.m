function [residual, g1, g2, g3] = diamond_1965_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(10, 1);
T14 = 1/((1+params(7))*(1+params(6)));
T47 = params(2)^((-1)/params(5));
T54 = 1+T47*(1+y(12))^((params(5)-1)/params(5));
lhs =y(2);
rhs =T14^params(1)*y(1)^params(1);
residual(1)= lhs-rhs;
lhs =y(3);
rhs =params(1)*T14^(params(1)-1)*y(1)^(params(1)-1)-params(3);
residual(2)= lhs-rhs;
lhs =y(4);
rhs =y(2)*(1-params(1));
residual(3)= lhs-rhs;
lhs =y(6)+y(8);
rhs =y(4);
residual(4)= lhs-rhs;
lhs =y(7);
rhs =y(8)*1/(1+params(7))*(1+y(3));
residual(5)= lhs-rhs;
lhs =y(8);
rhs =y(4)/T54;
residual(6)= lhs-rhs;
lhs =y(5);
rhs =y(6)+y(7)*1/(1+params(6));
residual(7)= lhs-rhs;
lhs =y(9);
rhs =y(8)+y(1)*(1-params(3))/((1+params(7))*(1+params(6)));
residual(8)= lhs-rhs;
lhs =y(10);
rhs =y(2)-y(1)*T14*(y(3)+params(3))-y(4);
residual(9)= lhs-rhs;
lhs =y(11);
rhs =y(6)^(-params(5))-(1+y(12))*params(2)*(1+params(7))^(-params(5))*y(7)^(-params(5));
residual(10)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(10, 12);

%
% Jacobian matrix
%

g1(1,2)=1;
g1(1,1)=(-(T14^params(1)*getPowerDeriv(y(1),params(1),1)));
g1(2,3)=1;
g1(2,1)=(-(params(1)*T14^(params(1)-1)*getPowerDeriv(y(1),params(1)-1,1)));
g1(3,2)=(-(1-params(1)));
g1(3,4)=1;
g1(4,4)=(-1);
g1(4,6)=1;
g1(4,8)=1;
g1(5,3)=(-(y(8)*1/(1+params(7))));
g1(5,7)=1;
g1(5,8)=(-(1/(1+params(7))*(1+y(3))));
g1(6,12)=(-((-(y(4)*T47*getPowerDeriv(1+y(12),(params(5)-1)/params(5),1)))/(T54*T54)));
g1(6,4)=(-(1/T54));
g1(6,8)=1;
g1(7,5)=1;
g1(7,6)=(-1);
g1(7,7)=(-(1/(1+params(6))));
g1(8,8)=(-1);
g1(8,1)=(-((1-params(3))/((1+params(7))*(1+params(6)))));
g1(8,9)=1;
g1(9,2)=(-1);
g1(9,3)=T14*y(1);
g1(9,4)=1;
g1(9,1)=T14*(y(3)+params(3));
g1(9,10)=1;
g1(10,12)=params(2)*(1+params(7))^(-params(5))*y(7)^(-params(5));
g1(10,6)=(-(getPowerDeriv(y(6),(-params(5)),1)));
g1(10,7)=(1+y(12))*params(2)*(1+params(7))^(-params(5))*getPowerDeriv(y(7),(-params(5)),1);
g1(10,11)=1;
end
if nargout >= 3,
%
% Hessian matrix
%

  g2 = sparse([],[],[],10,144);
end
if nargout >= 4,
%
% Third order derivatives
%

  g3 = sparse([],[],[],10,1728);
end
end
