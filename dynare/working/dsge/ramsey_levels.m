%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

clear all
clear global
tic;
global M_ oo_ options_ ys0_ ex0_ estimation_info
options_ = [];
M_.fname = 'ramsey_levels';
%
% Some global variables initialization
%
global_initialization;
diary off;
logname_ = 'ramsey_levels.log';
if exist(logname_, 'file')
    delete(logname_)
end
diary(logname_)
M_.exo_names = 'eps';
M_.exo_names_tex = 'eps';
M_.endo_names = 'k';
M_.endo_names_tex = 'k';
M_.endo_names = char(M_.endo_names, 'c');
M_.endo_names_tex = char(M_.endo_names_tex, 'c');
M_.endo_names = char(M_.endo_names, 'z');
M_.endo_names_tex = char(M_.endo_names_tex, 'z');
M_.endo_names = char(M_.endo_names, 'check1');
M_.endo_names_tex = char(M_.endo_names_tex, 'check1');
M_.param_names = 'beta';
M_.param_names_tex = 'beta';
M_.param_names = char(M_.param_names, 'theta');
M_.param_names_tex = char(M_.param_names_tex, 'theta');
M_.param_names = char(M_.param_names, 'alpha');
M_.param_names_tex = char(M_.param_names_tex, 'alpha');
M_.param_names = char(M_.param_names, 'delta');
M_.param_names_tex = char(M_.param_names_tex, 'delta');
M_.param_names = char(M_.param_names, 'g');
M_.param_names_tex = char(M_.param_names_tex, 'g');
M_.param_names = char(M_.param_names, 'n');
M_.param_names_tex = char(M_.param_names_tex, 'n');
M_.param_names = char(M_.param_names, 'rho');
M_.param_names_tex = char(M_.param_names_tex, 'rho');
M_.param_names = char(M_.param_names, 'sigma');
M_.param_names_tex = char(M_.param_names_tex, 'sigma');
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 4;
M_.param_nbr = 8;
M_.orig_endo_nbr = 4;
M_.aux_vars = [];
M_.Sigma_e = zeros(1, 1);
M_.H = 0;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('ramsey_levels_static');
erase_compiled_function('ramsey_levels_dynamic');
M_.lead_lag_incidence = [
 1 3 0;
 0 4 7;
 2 5 8;
 0 6 0;]';
M_.nstatic = 1;
M_.nfwrd   = 1;
M_.npred   = 1;
M_.nboth   = 1;
M_.equations_tags = {
};
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(4, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(8, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 17;
M_.NNZDerivatives(2) = -1;
M_.NNZDerivatives(3) = -1;
M_.params( 1 ) = 0.9896;
beta = M_.params( 1 );
M_.params( 2 ) = 2.0;
theta = M_.params( 2 );
M_.params( 3 ) = 0.4;
alpha = M_.params( 3 );
M_.params( 4 ) = 0.0196;
delta = M_.params( 4 );
M_.params( 5 ) = 0.0031875459362321018;
g = M_.params( 5 );
M_.params( 6 ) = 0.0037722357644148748;
n = M_.params( 6 );
M_.params( 7 ) = 0.95;
rho = M_.params( 7 );
M_.params( 8 ) = 0.01;
sigma = M_.params( 8 );
steady;
oo_.dr.eigval = check(M_,options_,oo_);
%
% SHOCKS instructions
%
make_ex_;
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (M_.params(8))^2;
M_.sigma_e_is_diagonal = 1;
options_.irf = 40;
options_.order = 1;
var_list_=[];
info = stoch_simul(var_list_);
save('ramsey_levels_results.mat', 'oo_', 'M_', 'options_');


disp(['Total computing time : ' dynsec2hms(toc) ]);
diary off
