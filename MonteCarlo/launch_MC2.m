cd C:\Users\marti\PhD\MySimMatlab\MultiPAniso\nucl_barrier\MonteCarlo
addpath '../shape_factor/'
clear
% 
resultsfile_sf = 'Results_shape_factors.mat';
resultsfile_sf = 'Results_shape_factors3.mat';
load(resultsfile_sf)

% ind =39;
ind =60;
res = results(ind);
exp_homog = 20; % between ca 10-78

%% MC simulation specification 
MCspec = gen_MCspec(resultsfile_sf, ind);
MCspec.exp_homog = exp_homog; % assumed between 20-300

%% launch

plotting = true;
results = run_MCv2(MCspec,plotting);



