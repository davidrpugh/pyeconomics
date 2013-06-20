% Usage:
%       out = ramsey_benchmark_ff(params, y)
%   where
%       out    is a (7,7) matrix of the first order
%              derivatives of the static system residuals
%              columns correspond to endo variables in
%              the ordering as declared
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

function out = ramsey_benchmark_ff(params, y)
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

t1 = a1;
t6 = exp(a5);
t9 = a7 ^ a8;
t10 = t6 * t9;
t51 = -(t10);
t54 = a8 - a1;
t55 = a7 ^ t54;
t56 = a8 * t55;
t57 = t6 * t56;
t58 = -(t57);
t14 = a13 - a8;
t93 = -(t14);
t18 = t6 * a8;
t19 = a8 - a13;
t20 = a7 ^ t19;
t21 = t18 * t20;
t94 = -(t21);
t96 = t19 - a1;
t97 = a7 ^ t96;
t98 = t19 * t97;
t99 = t18 * t98;
t100 = -(t99);
t135 = -(a1);
t37 = a13 + a36;
t30 = -(a29);
t138 = t30 - a1;
t139 = a33 ^ t138;
t140 = t30 * t139;
t141 = a32 * t140;
t142 = t37 * t141;
t143 = -(t142);
t34 = a33 ^ t30;
t35 = a32 * t34;
t144 = -(t35);
t147 = a25 ^ t138;
t148 = t30 * t147;
t41 = a13 - a22;
t213 = -(t41);
t214 = -(a45);
% setting the output variable
out = zeros(7, 7);
out(1,4) = out(1,4) + t1; % y(0)
out(1,3) = out(1,3) + t51; % z(0)
out(1,1) = out(1,1) + t58; % k(-1)
out(2,4) = out(2,4) + t93; % y(0)
out(2,6) = out(2,6) + t1; % w(0)
out(3,3) = out(3,3) + t94; % z(0)
out(3,7) = out(3,7) + t1; % r(0)
out(3,1) = out(3,1) + t100; % k(-1)
out(4,4) = out(4,4) + t1; % y(0)
out(4,5) = out(4,5) + t135; % i(0)
out(4,2) = out(4,2) + t135; % c(0)
out(5,2) = out(5,2) + t143; % c(1)
out(5,7) = out(5,7) + t144; % r(1)
out(5,2) = out(5,2) + t148; % c(0)
out(6,5) = out(6,5) + t135; % i(0)
out(6,1) = out(6,1) + t1; % k(0)
out(6,1) = out(6,1) + t213; % k(-1)
out(7,3) = out(7,3) + t1; % z(0)
out(7,3) = out(7,3) + t214; % z(-1)
