% Usage:
%       out = ramsey_growth_ff(params, y)
%   where
%       out    is a (6,6) matrix of the first order
%              derivatives of the static system residuals
%              columns correspond to endo variables in
%              the ordering as declared
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

function out = ramsey_growth_ff(params, y)
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

t1 = a1;
t7 = a5 * a6;
t10 = a8 + a9;
t12 = a8 + a11;
t13 = t10 * t12;
t69 = t7 * t13;
t70 = -(t69);
t15 = t13 * a14;
t71 = t15 * t15;
t72 = t70 / t71;
t16 = t7 / t15;
t75 = a17 - a1;
t76 = t16 ^ t75;
t77 = a17 * t76;
t78 = t72 * t77;
t79 = -(t78);
t80 = a5 / t15;
t81 = t77 * t80;
t82 = -(t81);
t83 = a6 / t15;
t84 = t77 * t83;
t85 = -(t84);
t970 = -(a1);
t26 = a8 - a25;
t971 = t26 * t72;
t972 = -(t971);
t973 = t26 * t80;
t974 = -(t973);
t975 = t26 * t83;
t976 = -(t975);
t32 = -(a31);
t33 = t10 ^ t32;
t34 = a30 * t33;
t37 = a35 * a36;
t38 = a14 * a20;
t39 = t37 / t38;
t40 = t39 ^ t32;
t41 = t34 * t40;
t42 = a14 * a24;
t1017 = t13 * t42;
t1018 = -(t1017);
t43 = t13 * a36;
t1019 = t43 * t43;
t1020 = t1018 / t1019;
t45 = a17 - a8;
t44 = t42 / t43;
t1023 = t45 - a1;
t1024 = t44 ^ t1023;
t1025 = t45 * t1024;
t1026 = t1020 * t1025;
t1027 = a17 * t1026;
t1028 = t41 * t1027;
t46 = t44 ^ t45;
t47 = a17 * t46;
t48 = a8 + t47;
t49 = t48 - a25;
t1009 = a35 / t38;
t1012 = t32 - a1;
t1013 = t39 ^ t1012;
t1014 = t32 * t1013;
t1015 = t1009 * t1014;
t1016 = t34 * t1015;
t1029 = t49 * t1016;
t1030 = t1028 + t1029;
t1031 = -(t1030);
t1032 = a36 / t38;
t1033 = t1014 * t1032;
t1034 = t34 * t1033;
t1035 = t49 * t1034;
t1036 = -(t1035);
t1037 = a14 / t43;
t1038 = t1025 * t1037;
t1039 = a17 * t1038;
t1040 = t41 * t1039;
t1041 = -(t1040);
t1048 = a24 / t43;
t1049 = t1025 * t1048;
t1050 = a17 * t1049;
t1051 = t41 * t1050;
t1042 = a20 * t37;
t1043 = -(t1042);
t1044 = t38 * t38;
t1045 = t1043 / t1044;
t1046 = t1014 * t1045;
t1047 = t34 * t1046;
t1052 = t49 * t1047;
t1053 = t1051 + t1052;
t1054 = -(t1053);
t1055 = a14 * t37;
t1056 = -(t1055);
t1057 = t1056 / t1044;
t1058 = t1014 * t1057;
t1059 = t34 * t1058;
t1060 = t49 * t1059;
t1061 = -(t1060);
t55 = exp(a54);
t6245 = a52 - a1;
t6246 = a5 ^ t6245;
t6247 = a52 * t6246;
t6248 = t55 * t6247;
t6249 = -(t6248);
t59 = a8 - a17;
t6285 = t59 * t78;
t6286 = -(t6285);
t6288 = t16 ^ t1023;
t6289 = t45 * t6288;
t6290 = t72 * t6289;
t6291 = a17 * t6290;
t6292 = t16 * t6291;
t62 = t16 ^ t45;
t63 = a17 * t62;
t64 = t63 - a25;
t65 = a25 + t64;
t6293 = t65 * t72;
t6294 = t6292 + t6293;
t6295 = t6286 - t6294;
t6296 = -(t6295);
t6297 = t59 * t81;
t6298 = -(t6297);
t6299 = t80 * t6289;
t6300 = a17 * t6299;
t6301 = t16 * t6300;
t6302 = t65 * t80;
t6303 = t6301 + t6302;
t6304 = t6298 - t6303;
t6305 = -(t6304);
t6306 = t59 * t84;
t6307 = -(t6306);
t6308 = t83 * t6289;
t6309 = a17 * t6308;
t6310 = t16 * t6309;
t6311 = t65 * t83;
t6312 = t6310 + t6311;
t6313 = t6307 - t6312;
t6314 = -(t6313);
% setting the output variable
out = zeros(6, 6);
out(1,1) = out(1,1) + t1; % y(0)
out(1,5) = out(1,5) + t79; % z(0)
out(1,2) = out(1,2) + t82; % k(-1)
out(1,5) = out(1,5) + t85; % z(-1)
out(2,1) = out(2,1) + t1; % y(0)
out(2,4) = out(2,4) + t970; % i(0)
out(2,3) = out(2,3) + t970; % c(0)
out(3,4) = out(3,4) + t970; % i(0)
out(3,2) = out(3,2) + t1; % k(0)
out(3,5) = out(3,5) + t972; % z(0)
out(3,2) = out(3,2) + t974; % k(-1)
out(3,5) = out(3,5) + t976; % z(-1)
out(4,5) = out(4,5) + t1031; % z(1)
out(4,3) = out(4,3) + t1036; % c(1)
out(4,2) = out(4,2) + t1041; % k(0)
out(4,5) = out(4,5) + t1054; % z(0)
out(4,3) = out(4,3) + t1061; % c(0)
out(5,5) = out(5,5) + t1; % z(0)
out(5,5) = out(5,5) + t6249; % z(-1)
out(6,1) = out(6,1) + t970; % y(0)
out(6,6) = out(6,6) + t1; % check1(0)
out(6,5) = out(6,5) + t6296; % z(0)
out(6,2) = out(6,2) + t6305; % k(-1)
out(6,5) = out(6,5) + t6314; % z(-1)
