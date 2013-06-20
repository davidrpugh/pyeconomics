% Usage:
%       out = ramsey_fixed_wages_ff(params, y)
%   where
%       out    is a (7,7) matrix of the first order
%              derivatives of the static system residuals
%              columns correspond to endo variables in
%              the ordering as declared
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

function out = ramsey_fixed_wages_ff(params, y)
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
t58 = t13 - a1;
t59 = a11 ^ t58;
t60 = t13 * t59;
t61 = t10 * t60;
t62 = -(t61);
t1 = a1;
t14 = a11 ^ t13;
t15 = t10 * t14;
t63 = -(t15);
t66 = a8 - a1;
t67 = a7 ^ t66;
t68 = a8 * t67;
t69 = t6 * t68;
t70 = t14 * t69;
t71 = -(t70);
t164 = -(a4);
t165 = a11 * a11;
t166 = t164 / t165;
t167 = t13 * t166;
t168 = -(t167);
t169 = a1 / a11;
t170 = t13 * t169;
t171 = -(t170);
t22 = t6 * a8;
t23 = a8 - a12;
t24 = a7 ^ t23;
t25 = t22 * t24;
t263 = t25 * t60;
t264 = -(t263);
t26 = t14 * t25;
t265 = -(t26);
t267 = t23 - a1;
t268 = a7 ^ t267;
t269 = t23 * t268;
t270 = t22 * t269;
t271 = t14 * t270;
t272 = -(t271);
t179 = -(a1);
t42 = a12 + a41;
t35 = -(a34);
t341 = t35 - a1;
t342 = a38 ^ t341;
t343 = t35 * t342;
t344 = a37 * t343;
t345 = t42 * t344;
t346 = -(t345);
t39 = a38 ^ t35;
t40 = a37 * t39;
t347 = -(t40);
t350 = a30 ^ t341;
t351 = t35 * t350;
t46 = a12 - a27;
t416 = -(t46);
t417 = -(a50);
% setting the output variable
out = zeros(7, 7);
out(1,3) = out(1,3) + t62; % l(0)
out(1,5) = out(1,5) + t1; % y(0)
out(1,4) = out(1,4) + t63; % z(0)
out(1,1) = out(1,1) + t71; % k(-1)
out(2,3) = out(2,3) + t168; % l(0)
out(2,5) = out(2,5) + t171; % y(0)
out(3,3) = out(3,3) + t264; % l(0)
out(3,4) = out(3,4) + t265; % z(0)
out(3,7) = out(3,7) + t1; % r(0)
out(3,1) = out(3,1) + t272; % k(-1)
out(4,5) = out(4,5) + t1; % y(0)
out(4,6) = out(4,6) + t179; % i(0)
out(4,2) = out(4,2) + t179; % c(0)
out(5,2) = out(5,2) + t346; % c(1)
out(5,7) = out(5,7) + t347; % r(1)
out(5,2) = out(5,2) + t351; % c(0)
out(6,6) = out(6,6) + t179; % i(0)
out(6,1) = out(6,1) + t1; % k(0)
out(6,1) = out(6,1) + t416; % k(-1)
out(7,4) = out(7,4) + t1; % z(0)
out(7,4) = out(7,4) + t417; % z(-1)
