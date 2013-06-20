function [residual, g1, g2] = ramsey_logs_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 3, 1);

%
% Model equations
%

T18 = exp(y(1))/((1+params(5))*(1+params(6)));
T37 = params(1)*(1+params(5))^(-params(2))*exp((-params(2))*(y(3)+y(2)-y(2)-y(3)));
lhs =exp(y(1));
rhs =(1-params(4))*T18+T18^params(3)-exp(y(2));
residual(1)= lhs-rhs;
lhs =1;
rhs =T37*(1+params(3)*T18^(params(3)-1)-params(4));
residual(2)= lhs-rhs;
lhs =exp(y(3));
rhs =exp(y(3)*params(7)+x(1));
residual(3)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(3, 3);

  %
  % Jacobian matrix
  %

  g1(1,1)=exp(y(1))-((1-params(4))*T18+T18*getPowerDeriv(T18,params(3),1));
  g1(1,2)=exp(y(2));
  g1(2,1)=(-(T37*params(3)*T18*getPowerDeriv(T18,params(3)-1,1)));
  g1(3,3)=exp(y(3))-params(7)*exp(y(3)*params(7)+x(1));
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],3,9);
end
end
