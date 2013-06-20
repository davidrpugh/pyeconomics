% Usage:
%       out = rbc_benchmark_ff(params, y)
%   where
%       out    is a (8,8) matrix of the first order
%              derivatives of the static system residuals
%              columns correspond to endo variables in
%              the ordering as declared
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

function out = rbc_benchmark_ff(params, y)
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

t1 = a1;
t6 = exp(a5);
t9 = a7 ^ a8;
t10 = t6 * t9;
t13 = a12 - a8;
t14 = a11 ^ t13;
t15 = t10 * t14;
t75 = -(t15);
t78 = t13 - a1;
t79 = a11 ^ t78;
t80 = t13 * t79;
t81 = t10 * t80;
t82 = -(t81);
t85 = a8 - a1;
t86 = a7 ^ t85;
t87 = a8 * t86;
t88 = t6 * t87;
t89 = t14 * t88;
t90 = -(t89);
t183 = a1 / a11;
t184 = t13 * t183;
t185 = -(t184);
t186 = -(a4);
t187 = a11 * a11;
t188 = t186 / t187;
t189 = t13 * t188;
t190 = -(t189);
t282 = a1 / a7;
t283 = a8 * t282;
t284 = -(t283);
t285 = a7 * a7;
t286 = t186 / t285;
t287 = a8 * t286;
t288 = -(t287);
t191 = -(a1);
t57 = a12 + a56;
t32 = a31 - a12;
t46 = a45 ^ t32;
t47 = a44 * t46;
t49 = a12 - a48;
t35 = a12 - a31;
t50 = t49 ^ t35;
t51 = t47 * t50;
t382 = a31 - a1;
t383 = a45 ^ t382;
t384 = a31 * t383;
t385 = t50 * t384;
t41 = -(a40);
t52 = a45 ^ a31;
t53 = t50 * t52;
t388 = t41 - a1;
t389 = t53 ^ t388;
t390 = t41 * t389;
t391 = t385 * t390;
t392 = t51 * t391;
t54 = t53 ^ t41;
t376 = t32 - a1;
t377 = a45 ^ t376;
t378 = t32 * t377;
t379 = a44 * t378;
t380 = t50 * t379;
t393 = t54 * t380;
t394 = t392 + t393;
t395 = t57 * t394;
t396 = -(t395);
t55 = t51 * t54;
t397 = -(t55);
t400 = t35 - a1;
t401 = t49 ^ t400;
t402 = t35 * t401;
t403 = t191 * t402;
t405 = t52 * t403;
t406 = t390 * t405;
t407 = t51 * t406;
t404 = t47 * t403;
t408 = t54 * t404;
t409 = t407 + t408;
t410 = t57 * t409;
t411 = -(t410);
t33 = a27 ^ t32;
t34 = a12 - a11;
t36 = t34 ^ t35;
t37 = t33 * t36;
t418 = a27 ^ t382;
t419 = a31 * t418;
t420 = t36 * t419;
t38 = a27 ^ a31;
t39 = t36 * t38;
t423 = t39 ^ t388;
t424 = t41 * t423;
t425 = t420 * t424;
t426 = t37 * t425;
t42 = t39 ^ t41;
t414 = a27 ^ t376;
t415 = t32 * t414;
t416 = t36 * t415;
t427 = t42 * t416;
t428 = t426 + t427;
t431 = t34 ^ t400;
t432 = t35 * t431;
t433 = t191 * t432;
t435 = t38 * t433;
t436 = t424 * t435;
t437 = t37 * t436;
t434 = t33 * t433;
t438 = t42 * t434;
t439 = t437 + t438;
t60 = a31 * t34;
t2124 = -(t35);
t2125 = a31 * t191;
t2126 = a17 * t2125;
t65 = a12 - a24;
t2127 = -(t65);
t2128 = -(a69);
% setting the output variable
out = zeros(8, 8);
out(1,4) = out(1,4) + t1; % y(0)
out(1,3) = out(1,3) + t75; % z(0)
out(1,8) = out(1,8) + t82; % l(0)
out(1,1) = out(1,1) + t90; % k(-1)
out(2,4) = out(2,4) + t185; % y(0)
out(2,6) = out(2,6) + t1; % w(0)
out(2,8) = out(2,8) + t190; % l(0)
out(3,4) = out(3,4) + t284; % y(0)
out(3,7) = out(3,7) + t1; % r(0)
out(3,1) = out(3,1) + t288; % k(-1)
out(4,4) = out(4,4) + t1; % y(0)
out(4,5) = out(4,5) + t191; % i(0)
out(4,2) = out(4,2) + t191; % c(0)
out(5,2) = out(5,2) + t396; % c(1)
out(5,7) = out(5,7) + t397; % r(1)
out(5,8) = out(5,8) + t411; % l(1)
out(5,2) = out(5,2) + t428; % c(0)
out(5,8) = out(5,8) + t439; % l(0)
out(6,6) = out(6,6) + t60; % w(0)
out(6,2) = out(6,2) + t2124; % c(0)
out(6,8) = out(6,8) + t2126; % l(0)
out(7,5) = out(7,5) + t191; % i(0)
out(7,1) = out(7,1) + t1; % k(0)
out(7,1) = out(7,1) + t2127; % k(-1)
out(8,3) = out(8,3) + t1; % z(0)
out(8,3) = out(8,3) + t2128; % z(-1)
