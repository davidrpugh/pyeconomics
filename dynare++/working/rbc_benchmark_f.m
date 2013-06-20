% Usage:
%       out = rbc_benchmark_f(params, y)
%   where
%       out    is a (8,1) column vector of the residuals
%              of the static system
%       params is a (9,1) vector of parameter values
%              in the ordering as declared
%       y      is a (8,1) vector of endogenous variables
%              in the ordering as declared
%
% Created by Dynare++ v. 4.3.1

% params ordering
% =====================
% beta
% theta
% omega
% alpha
% delta
% rho
% sigma
% kss
% lss
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
% l

function out = rbc_benchmark_f(params, y)
if size(y) ~= [8,1]
	error('Wrong size of y, must be [8,1]');
end
if size(params) ~= [9,1]
	error('Wrong size of params, must be [9,1]');
end

% hardwired constants
a0 =            0;
a1 =            1;
a2 = NaN;
a3 =    1.1283792;
% numerical constants
a12 =            1;
% parameter values
a44 = params(1); % beta
a40 = params(2); % theta
a31 = params(3); % omega
a8 = params(4); % alpha
a24 = params(5); % delta
a69 = params(6); % rho
% sigma not used in the model
% kss not used in the model
% lss not used in the model
% exogenous variables to zeros
a72 = 0.0; % eps
% endogenous variables to y
a7 = y(1); % k
a64 = y(1); % k
a27 = y(2); % c
a45 = y(2); % c
a70 = y(3); % z
a5 = y(3); % z
a4 = y(4); % y
a28 = y(5); % i
a17 = y(6); % w
a21 = y(7); % r
a56 = y(7); % r
a11 = y(8); % l
a48 = y(8); % l

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
t22 = a4 / a7;
t23 = a8 * t22;
t25 = t23 - a24;
t26 = a21 - t25;
t29 = a27 + a28;
t30 = a4 - t29;
t32 = a31 - a12;
t33 = a27 ^ t32;
t34 = a12 - a11;
t35 = a12 - a31;
t36 = t34 ^ t35;
t37 = t33 * t36;
t38 = a27 ^ a31;
t39 = t36 * t38;
t41 = -(a40);
t42 = t39 ^ t41;
t43 = t37 * t42;
t46 = a45 ^ t32;
t47 = a44 * t46;
t49 = a12 - a48;
t50 = t49 ^ t35;
t51 = t47 * t50;
t52 = a45 ^ a31;
t53 = t50 * t52;
t54 = t53 ^ t41;
t55 = t51 * t54;
t57 = a12 + a56;
t58 = t55 * t57;
t59 = t43 - t58;
t60 = a31 * t34;
t61 = a17 * t60;
t62 = a27 * t35;
t63 = t61 - t62;
t65 = a12 - a24;
t66 = a7 * t65;
t67 = a28 + t66;
t68 = a64 - t67;
t71 = a69 * a70;
t73 = t71 + a72;
t74 = a5 - t73;
% setting the output variable
out = zeros(8, 1);
out(1) = t16;
out(2) = t20;
out(3) = t26;
out(4) = t30;
out(5) = t59;
out(6) = t63;
out(7) = t68;
out(8) = t74;
