function [residual, g1, g2, g3] = rbc_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(7, 1);
T13 = params(1)*1/y(10);
T18 = params(4)*y(5)^(params(4)-1);
T28 = 1+T18*(exp(y(12))*y(11))^(1-params(4))-params(3);
T37 = y(1)^params(4);
T42 = (1-params(4))*T37*exp(y(9))^(1-params(4));
T44 = y(7)^(-params(4));
lhs =1/y(4);
rhs =T13*T28;
residual(1)= lhs-rhs;
lhs =y(4)*params(2)/(1-y(7));
rhs =T42*T44;
residual(2)= lhs-rhs;
lhs =y(4)+y(6);
rhs =y(3);
residual(3)= lhs-rhs;
lhs =y(3);
rhs =T37*(y(7)*exp(y(9)))^(1-params(4));
residual(4)= lhs-rhs;
lhs =y(6);
rhs =y(5)-y(1)*(1-params(3));
residual(5)= lhs-rhs;
lhs =y(8);
rhs =y(3)/y(7);
residual(6)= lhs-rhs;
lhs =y(9);
rhs =params(5)*y(2)+x(it_, 1);
residual(7)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(7, 13);

%
% Jacobian matrix
%

g1(1,4)=(-1)/(y(4)*y(4));
g1(1,10)=(-(T28*params(1)*(-1)/(y(10)*y(10))));
g1(1,5)=(-(T13*(exp(y(12))*y(11))^(1-params(4))*params(4)*getPowerDeriv(y(5),params(4)-1,1)));
g1(1,11)=(-(T13*T18*exp(y(12))*getPowerDeriv(exp(y(12))*y(11),1-params(4),1)));
g1(1,12)=(-(T13*T18*exp(y(12))*y(11)*getPowerDeriv(exp(y(12))*y(11),1-params(4),1)));
g1(2,4)=params(2)/(1-y(7));
g1(2,1)=(-(T44*exp(y(9))^(1-params(4))*(1-params(4))*getPowerDeriv(y(1),params(4),1)));
g1(2,7)=y(4)*params(2)/((1-y(7))*(1-y(7)))-T42*getPowerDeriv(y(7),(-params(4)),1);
g1(2,9)=(-(T44*(1-params(4))*T37*exp(y(9))*getPowerDeriv(exp(y(9)),1-params(4),1)));
g1(3,3)=(-1);
g1(3,4)=1;
g1(3,6)=1;
g1(4,3)=1;
g1(4,1)=(-((y(7)*exp(y(9)))^(1-params(4))*getPowerDeriv(y(1),params(4),1)));
g1(4,7)=(-(T37*exp(y(9))*getPowerDeriv(y(7)*exp(y(9)),1-params(4),1)));
g1(4,9)=(-(T37*y(7)*exp(y(9))*getPowerDeriv(y(7)*exp(y(9)),1-params(4),1)));
g1(5,1)=1-params(3);
g1(5,5)=(-1);
g1(5,6)=1;
g1(6,3)=(-(1/y(7)));
g1(6,7)=(-((-y(3))/(y(7)*y(7))));
g1(6,8)=1;
g1(7,2)=(-params(5));
g1(7,9)=1;
g1(7,13)=(-1);
end
if nargout >= 3,
%
% Hessian matrix
%

  g2 = sparse([],[],[],7,169);
end
if nargout >= 4,
%
% Third order derivatives
%

  g3 = sparse([],[],[],7,2197);
end
end
