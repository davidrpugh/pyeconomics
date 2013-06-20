% Usage:
%       out = ramsey_fixed_wages_f(params, y)
%   where
%       out    is a (7,1) column vector of the residuals
%              of the static system
%       params is a (7,1) vector of parameter values
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
% w
% rho
% sigma
%
% y ordering
% =====================
% k
% c
% l
% z
% y
% i
% r

function out = ramsey_fixed_wages_f(params, y)
if size(y) ~= [7,1]
	error('Wrong size of y, must be [7,1]');
end
if size(params) ~= [7,1]
	error('Wrong size of params, must be [7,1]');
end

% hardwired constants
a0 =            0;
a1 =            1;
a2 = NaN;
a3 =    1.1283792;
% numerical constants
a12 =            1;
% parameter values
a37 = params(1); % beta
a34 = params(2); % theta
a8 = params(3); % alpha
a27 = params(4); % delta
a17 = params(5); % w
a50 = params(6); % rho
% sigma not used in the model
% exogenous variables to zeros
a53 = 0.0; % eps
% endogenous variables to y
a7 = y(1); % k
a45 = y(1); % k
a30 = y(2); % c
a38 = y(2); % c
a11 = y(3); % l
a51 = y(4); % z
a5 = y(4); % z
a4 = y(5); % y
a31 = y(6); % i
a21 = y(7); % r
a41 = y(7); % r

t6 = exp(a5);
t9 = a7 ^ a8;
t10 = t6 * t9;
t13 = a12 - a8;
t14 = a11 ^ t13;
t15 = t10 * t14;
t16 = a4 - t15;
t18 = a4 / a11;
t19 = t13 * t18;
t20 = a17 - t19;
t22 = t6 * a8;
t23 = a8 - a12;
t24 = a7 ^ t23;
t25 = t22 * t24;
t26 = t14 * t25;
t28 = t26 - a27;
t29 = a21 - t28;
t32 = a30 + a31;
t33 = a4 - t32;
t35 = -(a34);
t36 = a30 ^ t35;
t39 = a38 ^ t35;
t40 = a37 * t39;
t42 = a12 + a41;
t43 = t40 * t42;
t44 = t36 - t43;
t46 = a12 - a27;
t47 = a7 * t46;
t48 = a31 + t47;
t49 = a45 - t48;
t52 = a50 * a51;
t54 = t52 + a53;
t55 = a5 - t54;
% setting the output variable
out = zeros(7, 1);
out(1) = t16;
out(2) = t20;
out(3) = t29;
out(4) = t33;
out(5) = t44;
out(6) = t49;
out(7) = t55;
