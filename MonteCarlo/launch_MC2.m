% COPYRIGHT NOTICE 
 % This file is part of a dataset <Minar, Martin (2023), “Influence of surface energy anisotropy on nucleation and crystallographic texture of polycrystalline deposits”, Mendeley Data, V1, doi: 10.17632/bsdff8shbz.1>, coupled to publication of the same name by Minar, Moelans submitted to Computational Materials Science in September 2023. 
 % Distributed under GPLv3 license.
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



