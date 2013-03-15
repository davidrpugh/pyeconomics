function [residual, g1, g2, g3] = diamond_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(10, 1);
T19 = 1/((1+params(6))*(1+params(7)));
T56 = params(2)^((-1)/params(5));
T64 = 1+T56*(1+y(12)-params(3))^((params(5)-1)/params(5));
T85 = T56*getPowerDeriv(1+y(12)-params(3),(params(5)-1)/params(5),1);
lhs =y(9);
rhs =y(2)+y(10);
residual(1)= lhs-rhs;
lhs =y(2);
rhs =y(3)+T19*y(4);
residual(2)= lhs-rhs;
lhs =y(9);
rhs =y(6)+y(7)*y(1);
residual(3)= lhs-rhs;
lhs =y(5);
rhs =T19*(y(10)+y(1)*(1-params(3)));
residual(4)= lhs-rhs;
lhs =y(7);
rhs =params(1)*y(1)^(params(1)-1);
residual(5)= lhs-rhs;
lhs =y(6);
rhs =y(9)*(1-params(1));
residual(6)= lhs-rhs;
lhs =y(3);
rhs =y(6)-y(8);
residual(7)= lhs-rhs;
lhs =y(4);
rhs =y(1)*(1+params(6))*(1+params(7))*(1+y(7)-params(3));
residual(8)= lhs-rhs;
lhs =y(8);
rhs =y(6)/T64;
residual(9)= lhs-rhs;
lhs =y(11);
rhs =y(5)-T19*y(8);
residual(10)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(10, 12);

%
% Jacobian matrix
%

g1(1,2)=(-1);
g1(1,9)=1;
g1(1,10)=(-1);
g1(2,2)=1;
g1(2,3)=(-1);
g1(2,4)=(-T19);
g1(3,1)=(-y(7));
g1(3,6)=(-1);
g1(3,7)=(-y(1));
g1(3,9)=1;
g1(4,1)=(-(T19*(1-params(3))));
g1(4,5)=1;
g1(4,10)=(-T19);
g1(5,1)=(-(params(1)*getPowerDeriv(y(1),params(1)-1,1)));
g1(5,7)=1;
g1(6,6)=1;
g1(6,9)=(-(1-params(1)));
g1(7,3)=1;
g1(7,6)=(-1);
g1(7,8)=1;
g1(8,4)=1;
g1(8,1)=(-((1+params(6))*(1+params(7))*(1+y(7)-params(3))));
g1(8,7)=(-((1+params(6))*(1+params(7))*y(1)));
g1(9,6)=(-(1/T64));
g1(9,12)=(-((-(y(6)*T85))/(T64*T64)));
g1(9,8)=1;
g1(10,5)=(-1);
g1(10,8)=T19;
g1(10,11)=1;
end
if nargout >= 3,
%
% Hessian matrix
%

  v2 = zeros(8,3);
v2(1,1)=3;
v2(1,2)=73;
v2(1,3)=(-1);
v2(2,1)=3;
v2(2,2)=7;
v2(2,3)=v2(1,3);
v2(3,1)=5;
v2(3,2)=1;
v2(3,3)=(-(params(1)*getPowerDeriv(y(1),params(1)-1,2)));
v2(4,1)=8;
v2(4,2)=73;
v2(4,3)=(-((1+params(6))*(1+params(7))));
v2(5,1)=8;
v2(5,2)=7;
v2(5,3)=v2(4,3);
v2(6,1)=9;
v2(6,2)=138;
v2(6,3)=(-((-T85)/(T64*T64)));
v2(7,1)=9;
v2(7,2)=72;
v2(7,3)=v2(6,3);
v2(8,1)=9;
v2(8,2)=144;
v2(8,3)=(-((T64*T64*(-(y(6)*T56*getPowerDeriv(1+y(12)-params(3),(params(5)-1)/params(5),2)))-(-(y(6)*T85))*(T64*T85+T64*T85))/(T64*T64*T64*T64)));
  g2 = sparse(v2(:,1),v2(:,2),v2(:,3),10,144);
end
if nargout >= 4,
%
% Third order derivatives
%

  g3 = sparse([],[],[],10,1728);
end
end
