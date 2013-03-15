function [residual, g1, g2, g3] = diamondv2_dynamic(y, x, params, steady_state, it_)
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
T24 = params(1)*T14^(params(1)-1);
T48 = params(2)^((-1)/params(5));
T53 = 1+T48*(1+y(12))^((params(5)-1)/params(5));
T79 = params(2)*(1+params(7))^(-params(5));
T91 = T48*getPowerDeriv(1+y(12),(params(5)-1)/params(5),1);
T93 = (-(y(4)*T91));
lhs =y(2);
rhs =T14^params(1)*y(1)^params(1);
residual(1)= lhs-rhs;
lhs =y(3);
rhs =T24*y(1)^(params(1)-1)-params(3);
residual(2)= lhs-rhs;
lhs =y(4);
rhs =y(2)*(1-params(1));
residual(3)= lhs-rhs;
lhs =y(6)+y(8);
rhs =y(4);
residual(4)= lhs-rhs;
lhs =y(13);
rhs =y(8)*1/(1+params(7))*(1+y(12));
residual(5)= lhs-rhs;
lhs =y(8);
rhs =y(4)/T53;
residual(6)= lhs-rhs;
lhs =y(5);
rhs =y(6)+1/(1+params(6))*y(7);
residual(7)= lhs-rhs;
lhs =y(9);
rhs =y(8)+y(1)*(1-params(3))/((1+params(7))*(1+params(6)));
residual(8)= lhs-rhs;
lhs =y(10);
rhs =y(2)-y(1)*T14*(y(3)+params(3))-y(4);
residual(9)= lhs-rhs;
lhs =y(11);
rhs =y(6)^(-params(5))-(1+y(12))*T79*y(7)^(-params(5));
residual(10)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(10, 13);

%
% Jacobian matrix
%

g1(1,2)=1;
g1(1,1)=(-(T14^params(1)*getPowerDeriv(y(1),params(1),1)));
g1(2,3)=1;
g1(2,1)=(-(T24*getPowerDeriv(y(1),params(1)-1,1)));
g1(3,2)=(-(1-params(1)));
g1(3,4)=1;
g1(4,4)=(-1);
g1(4,6)=1;
g1(4,8)=1;
g1(5,12)=(-(y(8)*1/(1+params(7))));
g1(5,13)=1;
g1(5,8)=(-(1/(1+params(7))*(1+y(12))));
g1(6,12)=(-(T93/(T53*T53)));
g1(6,4)=(-(1/T53));
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
g1(10,12)=T79*y(7)^(-params(5));
g1(10,6)=(-(getPowerDeriv(y(6),(-params(5)),1)));
g1(10,7)=(1+y(12))*T79*getPowerDeriv(y(7),(-params(5)),1);
g1(10,11)=1;
end
if nargout >= 3,
%
% Hessian matrix
%

  v2 = zeros(13,3);
v2(1,1)=1;
v2(1,2)=1;
v2(1,3)=(-(T14^params(1)*getPowerDeriv(y(1),params(1),2)));
v2(2,1)=2;
v2(2,2)=1;
v2(2,3)=(-(T24*getPowerDeriv(y(1),params(1)-1,2)));
v2(3,1)=5;
v2(3,2)=103;
v2(3,3)=(-(1/(1+params(7))));
v2(4,1)=5;
v2(4,2)=151;
v2(4,3)=v2(3,3);
v2(5,1)=6;
v2(5,2)=155;
v2(5,3)=(-((T53*T53*(-(y(4)*T48*getPowerDeriv(1+y(12),(params(5)-1)/params(5),2)))-T93*(T53*T91+T53*T91))/(T53*T53*T53*T53)));
v2(6,1)=6;
v2(6,2)=51;
v2(6,3)=(-((-T91)/(T53*T53)));
v2(7,1)=6;
v2(7,2)=147;
v2(7,3)=v2(6,3);
v2(8,1)=9;
v2(8,2)=3;
v2(8,3)=T14;
v2(9,1)=9;
v2(9,2)=27;
v2(9,3)=v2(8,3);
v2(10,1)=10;
v2(10,2)=71;
v2(10,3)=(-(getPowerDeriv(y(6),(-params(5)),2)));
v2(11,1)=10;
v2(11,2)=90;
v2(11,3)=T79*getPowerDeriv(y(7),(-params(5)),1);
v2(12,1)=10;
v2(12,2)=150;
v2(12,3)=v2(11,3);
v2(13,1)=10;
v2(13,2)=85;
v2(13,3)=(1+y(12))*T79*getPowerDeriv(y(7),(-params(5)),2);
  g2 = sparse(v2(:,1),v2(:,2),v2(:,3),10,169);
end
if nargout >= 4,
%
% Third order derivatives
%

  g3 = sparse([],[],[],10,2197);
end
end
