% COPYRIGHT NOTICE 
 % This file is part of a dataset <Minar, Martin (2023), “Influence of surface energy anisotropy on nucleation and crystallographic texture of polycrystalline deposits”, Mendeley Data, V1, doi: 10.17632/bsdff8shbz.1>, coupled to publication of the same name by Minar, Moelans submitted to Computational Materials Science in September 2023. 
 % Distributed under GPLv3 license.
%% launching different Monte Carlo simulations
% __ note that the current implementation assumes shape factor maps 50x50
%       - usage of different maps will result in wrong initial condition
%       and confusing results
%% single-run code usage demonstration 
clear
initialize_paths
% __ generate simulation spefification using gen_MCspec.m
help gen_MCspec.m
sf_resultsfile = 'shape_factor\results_SF.mat';
load(sf_resultsfile,'results_SF')
A = [(1:length(results_SF))' , [results_SF.n]', [results_SF.soaIE]', cellfun(@(x) size(x,1),{results_SF.S}','UniformOutput',true)];
results_overview = array2table(A,'VariableNames',{'index','n','anisotropy_strength','map_resolution'})

% __ example: 4-fold, 50x50, soaIE = 0.5
sf_ind = 5; 

% __ on-screen outpt, nucleation with anisotropic interface energy on,
% 300x200 system
MCspec = gen_MCspec(sf_resultsfile,sf_ind)

% __ run the simulation, only on-screen output
plotting = true;
run_MCv2(MCspec,plotting)

%% Use of launcher file for the published results
% __ the launcher allows reproducible launching of the three compared modes
% of nucleation in comparable conditions
%   - no nucleation
%   - nucleation with isotropic interface energy
%   - nucleation with anisotropic interface energy
clear
initialize_paths

% __ the launching the script runs 3x10 individual
% simulations in a system 50x50 (results in the paper based on 1500x200)
% and saves the output to 'MonteCarlo\demo_2D\results_MC2_demo.mat'
open launch_MC2_demo.m
launch_MC2_demo
% __ load the produced data
load('MonteCarlo\demo_2D\results_MC2_demo.mat','results')

% __ preview the anisotropic-nucleation deposit in results(3)
ind = 3;
figure(1)
ax = gca;
display_nucleation_sites = true;
plotMC2D_deposit(results,3,ax, display_nucleation_sites)
% __ see the through-thickness orientation evolution
fontsize = 11;
ylims_top = 0.25;
plotMC2D_ori_through_thicnkess(results, ind,fontsize, ylims_top)

%% Run and plot the 2-grain validation simulation
clear
initialize_paths
sf_resultsfile = 'shape_factor\results_SF.mat';
load(sf_resultsfile,'results_SF');
sf_ind = 6;
% MC simulation specification 
MCspec = gen_MCspec(sf_resultsfile,sf_ind);
MCspec = MCspec_2grains(MCspec);

MCspec.nucleation_on = false;
MCspec.reps = 1;
MCspec.exp_homog = 10;
MCspec.nucl_fluct.code = 'none';
resultsfile = '';
MCspec.outfile = resultsfile;
plotting = true;
% __ run the simulation
results = run_MCv2(MCspec,plotting);
