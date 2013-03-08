function [residual, g1, g2, g3] = diamond_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(11, 1);
T55 = params(2)^((-1)/params(5));
T61 = 1+T55*(1+y(14)-params(3))^((-(1-params(5)))/params(5));
T91 = T55*getPowerDeriv(1+y(14)-params(3),(-(1-params(5)))/params(5),1);
lhs =y(6);
rhs =y(3)+y(8);
residual(1)= lhs-rhs;
lhs =y(2);
rhs =y(3)+1/((1+params(6))*(1+params(7)))*y(4);
residual(2)= lhs-rhs;
lhs =y(10);
rhs =y(2)+y(9);
residual(3)= lhs-rhs;
lhs =y(13);
rhs =y(8)*(1+y(14)-params(3));
residual(4)= lhs-rhs;
lhs =y(5);
rhs =y(9)+(1-params(3))/((1+params(6))*(1+params(7)))*y(1);
residual(5)= lhs-rhs;
lhs =y(10);
rhs =y(1)^params(1);
residual(6)= lhs-rhs;
lhs =y(6);
rhs =y(10)*(1-params(1));
residual(7)= lhs-rhs;
lhs =y(7);
rhs =y(10)*params(1)/y(5);
residual(8)= lhs-rhs;
lhs =y(8);
rhs =y(6)/T61;
residual(9)= lhs-rhs;
lhs =y(11);
rhs =y(10)-(y(5)-(1-params(3))/((1+params(6))*(1+params(7)))*y(1))-y(2);
residual(10)= lhs-rhs;
lhs =y(12);
rhs =y(4)-y(1)*(1+params(6))*(1+params(7))*(1+y(7)-params(3));
residual(11)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(11, 14);

%
% Jacobian matrix
%

g1(1,3)=(-1);
g1(1,6)=1;
g1(1,8)=(-1);
g1(2,2)=1;
g1(2,3)=(-1);
g1(2,4)=(-(1/((1+params(6))*(1+params(7)))));
g1(3,2)=(-1);
g1(3,9)=(-1);
g1(3,10)=1;
g1(4,13)=1;
g1(4,14)=(-y(8));
g1(4,8)=(-(1+y(14)-params(3)));
g1(5,1)=(-((1-params(3))/((1+params(6))*(1+params(7)))));
g1(5,5)=1;
g1(5,9)=(-1);
g1(6,1)=(-(getPowerDeriv(y(1),params(1),1)));
g1(6,10)=1;
g1(7,6)=1;
g1(7,10)=(-(1-params(1)));
g1(8,5)=(-((-(y(10)*params(1)))/(y(5)*y(5))));
g1(8,7)=1;
g1(8,10)=(-(params(1)/y(5)));
g1(9,6)=(-(1/T61));
g1(9,14)=(-((-(y(6)*T91))/(T61*T61)));
g1(9,8)=1;
g1(10,2)=1;
g1(10,1)=(-((1-params(3))/((1+params(6))*(1+params(7)))));
g1(10,5)=1;
g1(10,10)=(-1);
g1(10,11)=1;
g1(11,4)=(-1);
g1(11,1)=(1+params(6))*(1+params(7))*(1+y(7)-params(3));
g1(11,7)=(1+params(6))*(1+params(7))*y(1);
g1(11,12)=1;
end
if nargout >= 3,
%
% Hessian matrix
%

  v2 = zeros(11,3);
v2(1,1)=4;
v2(1,2)=112;
v2(1,3)=(-1);
v2(2,1)=4;
v2(2,2)=190;
v2(2,3)=v2(1,3);
v2(3,1)=6;
v2(3,2)=1;
v2(3,3)=(-(getPowerDeriv(y(1),params(1),2)));
v2(4,1)=8;
v2(4,2)=61;
v2(4,3)=(-((-((-(y(10)*params(1)))*(y(5)+y(5))))/(y(5)*y(5)*y(5)*y(5))));
v2(5,1)=8;
v2(5,2)=131;
v2(5,3)=(-((-params(1))/(y(5)*y(5))));
v2(6,1)=8;
v2(6,2)=66;
v2(6,3)=v2(5,3);
v2(7,1)=9;
v2(7,2)=188;
v2(7,3)=(-((-T91)/(T61*T61)));
v2(8,1)=9;
v2(8,2)=84;
v2(8,3)=v2(7,3);
v2(9,1)=9;
v2(9,2)=196;
v2(9,3)=(-((T61*T61*(-(y(6)*T55*getPowerDeriv(1+y(14)-params(3),(-(1-params(5)))/params(5),2)))-(-(y(6)*T91))*(T61*T91+T61*T91))/(T61*T61*T61*T61)));
v2(10,1)=11;
v2(10,2)=85;
v2(10,3)=(1+params(6))*(1+params(7));
v2(11,1)=11;
v2(11,2)=7;
v2(11,3)=v2(10,3);
  g2 = sparse(v2(:,1),v2(:,2),v2(:,3),11,196);
end
if nargout >= 4,
%
% Third order derivatives
%

  g3 = sparse([],[],[],11,2744);
end
end
