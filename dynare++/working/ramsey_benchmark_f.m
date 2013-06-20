% Usage:
%       out = ramsey_benchmark_f(params, y)
%   where
%       out    is a (7,1) column vector of the residuals
%              of the static system
%       params is a (6,1) vector of parameter values
%              in the ordering as declared
%       y      is a (7,1) vector of endogenous variables
%              in the ordering as declared
%
% Created by Dynare++ v. 4.3.1

% params ordering
% =====================
% beta
% theta
% alpha
% delta
% rho
% sigma
%
% y ordering
% =====================
% k
% c
% z
% y
% i
% w
% r

function out = ramsey_benchmark_f(params, y)
if size(y) ~= [7,1]
	error('Wrong size of y, must be [7,1]');
end
if size(params) ~= [6,1]
	error('Wrong size of params, must be [6,1]');
end

% hardwired constants
a0 =            0;
a1 =            1;
a2 = NaN;
a3 =    1.1283792;
% numerical constants
a13 =            1;
% parameter values
a32 = params(1); % beta
a29 = params(2); % theta
a8 = params(3); % alpha
a22 = params(4); % delta
a45 = params(5); % rho
% sigma not used in the model
% exogenous variables to zeros
a48 = 0.0; % eps
% endogenous variables to y
a7 = y(1); % k
a40 = y(1); % k
a25 = y(2); % c
a33 = y(2); % c
a46 = y(3); % z
a5 = y(3); % z
a4 = y(4); % y
a26 = y(5); % i
a12 = y(6); % w
a17 = y(7); % r
a36 = y(7); % r

t6 = exp(a5);
t9 = a7 ^ a8;
t10 = t6 * t9;
t11 = a4 - t10;
t14 = a13 - a8;
t15 = a4 * t14;
t16 = a12 - t15;
t18 = t6 * a8;
t19 = a8 - a13;
t20 = a7 ^ t19;
t21 = t18 * t20;
t23 = t21 - a22;
t24 = a17 - t23;
t27 = a25 + a26;
t28 = a4 - t27;
t30 = -(a29);
t31 = a25 ^ t30;
t34 = a33 ^ t30;
t35 = a32 * t34;
t37 = a13 + a36;
t38 = t35 * t37;
t39 = t31 - t38;
t41 = a13 - a22;
t42 = a7 * t41;
t43 = a26 + t42;
t44 = a40 - t43;
t47 = a45 * a46;
t49 = t47 + a48;
t50 = a5 - t49;
% setting the output variable
out = zeros(7, 1);
out(1) = t11;
out(2) = t16;
out(3) = t24;
out(4) = t28;
out(5) = t39;
out(6) = t44;
out(7) = t50;
