function [residual, g1, g2, g3] = ramsey_levels_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(4, 1);
T21 = y(2)*y(1)/((1+params(5))*(1+params(6))*y(5));
T24 = T21^params(3);
T33 = params(1)*(1+params(5))^(-params(2));
T40 = T33*(y(7)*y(8)/(y(5)*y(4)))^(-params(2));
T48 = 1+params(3)*(y(3)*y(5)/((1+params(5))*(1+params(6))*y(8)))^(params(3)-1)-params(4);
T64 = params(4)+params(3)*T21^(params(3)-1)-params(4);
T68 = y(2)/((1+params(5))*(1+params(6))*y(5));
T71 = T68*getPowerDeriv(T21,params(3),1);
T76 = getPowerDeriv(T21,params(3)-1,1);
T85 = getPowerDeriv(y(3)*y(5)/((1+params(5))*(1+params(6))*y(8)),params(3)-1,1);
T94 = getPowerDeriv(y(7)*y(8)/(y(5)*y(4)),(-params(2)),1);
T104 = y(1)/((1+params(5))*(1+params(6))*y(5));
T106 = getPowerDeriv(T21,params(3),1)*T104;
T124 = (-(y(2)*y(1)*(1+params(5))*(1+params(6))))/((1+params(5))*(1+params(6))*y(5)*(1+params(5))*(1+params(6))*y(5));
T126 = getPowerDeriv(T21,params(3),1)*T124;
lhs =y(3);
rhs =(1-params(4))*T21+T24-y(4);
residual(1)= lhs-rhs;
lhs =1;
rhs =T40*T48;
residual(2)= lhs-rhs;
lhs =y(5);
rhs =y(2)^params(7)*exp(x(it_, 1));
residual(3)= lhs-rhs;
lhs =y(6);
rhs =T24-T24*(1-params(3))-T21*T64;
residual(4)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(4, 9);

  %
  % Jacobian matrix
  %

  g1(1,1)=(-((1-params(4))*T68+T71));
  g1(1,3)=1;
  g1(1,4)=1;
  g1(1,2)=(-((1-params(4))*T104+T106));
  g1(1,5)=(-((1-params(4))*T124+T126));
  g1(2,3)=(-(T40*params(3)*y(5)/((1+params(5))*(1+params(6))*y(8))*T85));
  g1(2,4)=(-(T48*T33*(-(y(5)*y(7)*y(8)))/(y(5)*y(4)*y(5)*y(4))*T94));
  g1(2,7)=(-(T48*T33*T94*y(8)/(y(5)*y(4))));
  g1(2,5)=(-(T48*T33*T94*(-(y(4)*y(7)*y(8)))/(y(5)*y(4)*y(5)*y(4))+T40*params(3)*T85*y(3)/((1+params(5))*(1+params(6))*y(8))));
  g1(2,8)=(-(T48*T33*T94*y(7)/(y(5)*y(4))+T40*params(3)*T85*(-((1+params(5))*(1+params(6))*y(3)*y(5)))/((1+params(5))*(1+params(6))*y(8)*(1+params(5))*(1+params(6))*y(8))));
  g1(3,2)=(-(exp(x(it_, 1))*getPowerDeriv(y(2),params(7),1)));
  g1(3,5)=1;
  g1(3,9)=(-(y(2)^params(7)*exp(x(it_, 1))));
  g1(4,1)=(-(T71-(1-params(3))*T71-(T64*T68+T21*params(3)*T68*T76)));
  g1(4,2)=(-(T106-(1-params(3))*T106-(T64*T104+T21*params(3)*T76*T104)));
  g1(4,5)=(-(T126-(1-params(3))*T126-(T64*T124+T21*params(3)*T76*T124)));
  g1(4,6)=1;
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],4,81);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],4,729);
end
end
