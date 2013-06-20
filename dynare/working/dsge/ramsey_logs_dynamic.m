function [residual, g1, g2, g3] = ramsey_logs_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(3, 1);
T23 = exp(y(2)-y(5)+y(1))/((1+params(5))*(1+params(6)));
T36 = params(1)*(1+params(5))^(-params(2));
T44 = T36*exp((-params(2))*(y(6)+y(7)-y(4)-y(5)));
T48 = exp(y(3)+y(5)-y(7))/((1+params(5))*(1+params(6)));
T53 = 1+params(3)*T48^(params(3)-1)-params(4);
lhs =exp(y(3));
rhs =(1-params(4))*T23+T23^params(3)-exp(y(4));
residual(1)= lhs-rhs;
lhs =1;
rhs =T44*T53;
residual(2)= lhs-rhs;
lhs =exp(y(5));
rhs =exp(y(2)*params(7)+x(it_, 1));
residual(3)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(3, 8);

  %
  % Jacobian matrix
  %

  g1(1,1)=(-((1-params(4))*T23+T23*getPowerDeriv(T23,params(3),1)));
  g1(1,3)=exp(y(3));
  g1(1,4)=exp(y(4));
  g1(1,2)=(-((1-params(4))*T23+T23*getPowerDeriv(T23,params(3),1)));
  g1(1,5)=(-((1-params(4))*(-exp(y(2)-y(5)+y(1)))/((1+params(5))*(1+params(6)))+getPowerDeriv(T23,params(3),1)*(-exp(y(2)-y(5)+y(1)))/((1+params(5))*(1+params(6)))));
  g1(2,3)=(-(T44*params(3)*T48*getPowerDeriv(T48,params(3)-1,1)));
  g1(2,4)=(-(T53*T36*params(2)*exp((-params(2))*(y(6)+y(7)-y(4)-y(5)))));
  g1(2,6)=(-(T53*T36*(-params(2))*exp((-params(2))*(y(6)+y(7)-y(4)-y(5)))));
  g1(2,5)=(-(T44*params(3)*T48*getPowerDeriv(T48,params(3)-1,1)+T53*T36*params(2)*exp((-params(2))*(y(6)+y(7)-y(4)-y(5)))));
  g1(2,7)=(-(T53*T36*(-params(2))*exp((-params(2))*(y(6)+y(7)-y(4)-y(5)))+T44*params(3)*getPowerDeriv(T48,params(3)-1,1)*(-exp(y(3)+y(5)-y(7)))/((1+params(5))*(1+params(6)))));
  g1(3,2)=(-(params(7)*exp(y(2)*params(7)+x(it_, 1))));
  g1(3,5)=exp(y(5));
  g1(3,8)=(-exp(y(2)*params(7)+x(it_, 1)));
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],3,64);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],3,512);
end
end
