%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

clear all
clear global
tic;
global M_ oo_ options_ ys0_ ex0_
options_ = [];
M_.fname = 'diamond_1965';
%
% Some global variables initialization
%
global_initialization;
diary off;
logname_ = 'diamond_1965.log';
if exist(logname_, 'file')
    delete(logname_)
end
diary(logname_)
M_.endo_names = 'y';
M_.endo_names_tex = 'y';
M_.endo_names = char(M_.endo_names, 'r');
M_.endo_names_tex = char(M_.endo_names_tex, 'r');
M_.endo_names = char(M_.endo_names, 'w');
M_.endo_names_tex = char(M_.endo_names_tex, 'w');
M_.endo_names = char(M_.endo_names, 'c');
M_.endo_names_tex = char(M_.endo_names_tex, 'c');
M_.endo_names = char(M_.endo_names, 'c1');
M_.endo_names_tex = char(M_.endo_names_tex, 'c1');
M_.endo_names = char(M_.endo_names, 'c2');
M_.endo_names_tex = char(M_.endo_names_tex, 'c2');
M_.endo_names = char(M_.endo_names, 's');
M_.endo_names_tex = char(M_.endo_names_tex, 's');
M_.endo_names = char(M_.endo_names, 'k');
M_.endo_names_tex = char(M_.endo_names_tex, 'k');
M_.endo_names = char(M_.endo_names, 'check1');
M_.endo_names_tex = char(M_.endo_names_tex, 'check1');
M_.endo_names = char(M_.endo_names, 'check2');
M_.endo_names_tex = char(M_.endo_names_tex, 'check2');
M_.param_names = 'alpha';
M_.param_names_tex = 'alpha';
M_.param_names = char(M_.param_names, 'beta');
M_.param_names_tex = char(M_.param_names_tex, 'beta');
M_.param_names = char(M_.param_names, 'delta');
M_.param_names_tex = char(M_.param_names_tex, 'delta');
M_.param_names = char(M_.param_names, 'rho');
M_.param_names_tex = char(M_.param_names_tex, 'rho');
M_.param_names = char(M_.param_names, 'theta');
M_.param_names_tex = char(M_.param_names_tex, 'theta');
M_.param_names = char(M_.param_names, 'n');
M_.param_names_tex = char(M_.param_names_tex, 'n');
M_.param_names = char(M_.param_names, 'g');
M_.param_names_tex = char(M_.param_names_tex, 'g');
M_.param_names = char(M_.param_names, 'T');
M_.param_names_tex = char(M_.param_names_tex, 'T');
M_.exo_det_nbr = 0;
M_.exo_nbr = 0;
M_.endo_nbr = 10;
M_.param_nbr = 8;
M_.orig_endo_nbr = 10;
M_.aux_vars = [];
M_.Sigma_e = zeros(0, 0);
M_.H = 0;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('diamond_1965_dynamic');
M_.lead_lag_incidence = [
 0 2 0;
 0 3 12;
 0 4 0;
 0 5 0;
 0 6 0;
 0 7 0;
 0 8 0;
 1 9 0;
 0 10 0;
 0 11 0;]';
M_.equations_tags = {
};
M_.exo_names_orig_ord = [1:0];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(10, 1);
M_.params = NaN(8, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 30;
M_.NNZDerivatives(2) = 0;
M_.NNZDerivatives(3) = -1;
M_.params( 8 ) = 30;
T = M_.params( 8 );
M_.params( 4 ) = 0.03;
rho = M_.params( 4 );
M_.params( 2 ) = (1-M_.params(4))^M_.params(8);
beta = M_.params( 2 );
M_.params( 6 ) = 1.02^M_.params(8)-1;
n = M_.params( 6 );
M_.params( 7 ) = 1.02^M_.params(8)-1;
g = M_.params( 7 );
M_.params( 1 ) = 0.33;
alpha = M_.params( 1 );
M_.params( 3 ) = 1-0.9^M_.params(8);
delta = M_.params( 3 );
M_.params( 5 ) = 1.0;
theta = M_.params( 5 );
kss  = 1.0;
yss  = (1 / ((1 + g) * (1 + n)))^alpha * kss^alpha;
wss  = (1 - alpha) * yss;
rss  = alpha * (1 / ((1 + g) * (1 + n)))^(alpha - 1) * kss^(alpha - 1) - delta; 
sss  = wss / (1 + beta^(- 1 / theta) * (1 + rss - delta)^((theta - 1) / theta));
c1ss = wss - sss;
c2ss = (1 / (1 + g)) * (1 + rss) * sss;
css  = c1ss + (1 / (1 + n)) * c2ss;  
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.steady_state( 8 ) = kss;
oo_.steady_state( 5 ) = c1ss;
oo_.steady_state( 6 ) = c2ss;
oo_.steady_state( 2 ) = rss;
oo_.steady_state( 3 ) = wss;
oo_.steady_state( 7 ) = sss;
oo_.steady_state( 4 ) = css;
oo_.steady_state( 1 ) = yss;
oo_.endo_simul=[oo_.steady_state*ones(1,M_.maximum_lag)];
if M_.exo_nbr > 0;
	oo_.exo_simul = [ones(M_.maximum_lag,1)*oo_.exo_steady_state'];
end;
if M_.exo_det_nbr > 0;
	oo_.exo_det_simul = [ones(M_.maximum_lag,1)*oo_.exo_det_steady_state'];
end;
steady;
check;
options_.periods = 100;
options_.datafile = 'diamond_1965_initvals.m';
simul();
save('diamond_1965_results.mat', 'oo_', 'M_', 'options_');
diary off

disp(['Total computing time : ' dynsec2hms(toc) ]);
