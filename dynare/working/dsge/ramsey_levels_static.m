function [residual, g1, g2] = ramsey_levels_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 4, 1);

%
% Model equations
%

T19 = y(1)*y(3)/(y(3)*(1+params(5))*(1+params(6)));
T31 = params(1)*(1+params(5))^(-params(2));
T51 = params(4)+params(3)*T19^(params(3)-1)-params(4);
T55 = y(3)/(y(3)*(1+params(5))*(1+params(6)));
T58 = T55*getPowerDeriv(T19,params(3),1);
T77 = (y(1)*y(3)*(1+params(5))*(1+params(6))-y(1)*y(3)*(1+params(5))*(1+params(6)))/(y(3)*(1+params(5))*(1+params(6))*y(3)*(1+params(5))*(1+params(6)));
T79 = getPowerDeriv(T19,params(3),1)*T77;
lhs =y(1);
rhs =(1-params(4))*T19+T19^params(3)-y(2);
residual(1)= lhs-rhs;
lhs =1;
rhs =T31*(1+params(3)*T19^(params(3)-1)-params(4));
residual(2)= lhs-rhs;
lhs =y(3);
rhs =y(3)^params(7)*exp(x(1));
residual(3)= lhs-rhs;
lhs =y(4);
rhs =T19^params(3)-T19^params(3)*(1-params(3))-T19*T51;
residual(4)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(4, 4);

  %
  % Jacobian matrix
  %

  g1(1,1)=1-((1-params(4))*T55+T58);
  g1(1,2)=1;
  g1(1,3)=(-((1-params(4))*T77+T79));
  g1(2,1)=(-(T31*params(3)*T55*getPowerDeriv(T19,params(3)-1,1)));
  g1(2,3)=(-(T31*params(3)*getPowerDeriv(T19,params(3)-1,1)*T77));
  g1(3,3)=1-exp(x(1))*getPowerDeriv(y(3),params(7),1);
  g1(4,1)=(-(T58-(1-params(3))*T58-(T51*T55+T19*params(3)*T55*getPowerDeriv(T19,params(3)-1,1))));
  g1(4,3)=(-(T79-(1-params(3))*T79-(T51*T77+T19*params(3)*getPowerDeriv(T19,params(3)-1,1)*T77)));
  g1(4,4)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],4,16);
end
end
