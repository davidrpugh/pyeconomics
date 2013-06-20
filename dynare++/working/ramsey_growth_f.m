% Usage:
%       out = ramsey_growth_f(params, y)
%   where
%       out    is a (6,1) column vector of the residuals
%              of the static system
%       params is a (8,1) vector of parameter values
%              in the ordering as declared
%       y      is a (6,1) vector of endogenous variables
%              in the ordering as declared
%
% Created by Dynare++ v. 4.3.1

% params ordering
% =====================
% beta
% theta
% alpha
% delta
% n
% g
% rho
% sigma
%
% y ordering
% =====================
% y
% k
% c
% i
% z
% check1

function out = ramsey_growth_f(params, y)
if size(y) ~= [6,1]
	error('Wrong size of y, must be [6,1]');
end
if size(params) ~= [8,1]
	error('Wrong size of params, must be [8,1]');
end

% hardwired constants
a0 =            0;
a1 =            1;
a2 = NaN;
a3 =    1.1283792;
% numerical constants
a8 =            1;
% parameter values
a30 = params(1); % beta
a31 = params(2); % theta
a17 = params(3); % alpha
a25 = params(4); % delta
a11 = params(5); % n
a9 = params(6); % g
a52 = params(7); % rho
% sigma not used in the model
% exogenous variables to zeros
a54 = 0.0; % eps
% endogenous variables to y
a4 = y(1); % y
a6 = y(2); % k
a24 = y(2); % k
a20 = y(3); % c
a35 = y(3); % c
a21 = y(4); % i
a5 = y(5); % z
a14 = y(5); % z
a36 = y(5); % z
a58 = y(6); % check1

t7 = a5 * a6;
t10 = a8 + a9;
t12 = a8 + a11;
t13 = t10 * t12;
t15 = t13 * a14;
t16 = t7 / t15;
t18 = t16 ^ a17;
t19 = a4 - t18;
t22 = a20 + a21;
t23 = a4 - t22;
t26 = a8 - a25;
t27 = t16 * t26;
t28 = a21 + t27;
t29 = a24 - t28;
t32 = -(a31);
t33 = t10 ^ t32;
t34 = a30 * t33;
t37 = a35 * a36;
t38 = a14 * a20;
t39 = t37 / t38;
t40 = t39 ^ t32;
t41 = t34 * t40;
t42 = a14 * a24;
t43 = t13 * a36;
t44 = t42 / t43;
t45 = a17 - a8;
t46 = t44 ^ t45;
t47 = a17 * t46;
t48 = a8 + t47;
t49 = t48 - a25;
t50 = t41 * t49;
t51 = a8 - t50;
t53 = a5 ^ a52;
t55 = exp(a54);
t56 = t53 * t55;
t57 = a14 - t56;
t59 = a8 - a17;
t60 = t18 * t59;
t61 = a4 - t60;
t62 = t16 ^ t45;
t63 = a17 * t62;
t64 = t63 - a25;
t65 = a25 + t64;
t66 = t16 * t65;
t67 = t61 - t66;
t68 = a58 - t67;
% setting the output variable
out = zeros(6, 1);
out(1) = t19;
out(2) = t23;
out(3) = t29;
out(4) = t51;
out(5) = t57;
out(6) = t68;
