% __ the anisotropy function is f(\theta) = 1 + soaIE*cos(n*(theta-alpha))
% __    n ... order of symmetry
% __    soaIE ... strength of anisotropy in interface energy
% __    theta... interface normal angle
% __    alpha ... rotation

%% visualize the shape factor-orientation maps and Winterbottom construction for its individual points
initialize_paths
% __ load the shape factor maps presented in the paper
load('shape_factor\results_SF.mat')
% __ legend to the results
% ind = 1:7; % 4-fold, 50x50, soaIE = 0.1:.1:07
% ind = 8:14; % 4-fold, 180x180, soaIE = 0.1:.1:07
% ind = 15:21; % 6-fold, 50x50, soaIE = 0.1:.1:07
% ind = 22:28; % 6-fold, 180x180, soaIE = 0.1:.1:07
% __ see results_overview
A = [(1:length(results_SF))' , [results_SF.n]', [results_SF.soaIE]', cellfun(@(x) size(x,1),{results_SF.S}','UniformOutput',true)];
results_overview = array2table(A,'VariableNames',{'index','n','anisotropy_strength','map_resolution'})

% __ example: 4-fold, 180x180, soaIE = 0.6
ind = 13; 

% __view the shape factor-orientation map with the stability condition
% applied
figure
imagesc(results_SF(ind).S,'AlphaData',results_SF(ind).stab_cond)
set(gca,'YDir','normal','DataAspectRatio',[1,1,1])
title([num2str(results_SF(ind).n) '-fold, \delta=' num2str(results_SF(ind).soaIE) ])

% __ explore the Winterbottom construction in individual points in the map
% (click on the point of interest once the figure opens)
plot_shapefacetor_wrapper(results_SF(ind))


%% compute a shape factor-orientation map
initialize_paths

% __ open launching file with default input and adjust strength of anisotropy
% and/or order of symmetry if you like
% __ note that in 3-fold symmetry the inverted solution is not used (but should be)
open run_ori_prob_maps.m
% __ run the script to carry out the computation
run_ori_prob_maps

%% comment on the implementation of Winterbottom construction
% - the function 'find_stable_solution.m' carries out the actual Winterbottom
% construction (and area comparison in case of multiple solutions)
% - it is used in 
%   - run_ori_prob_maps > get_shape_factor_oridep > parser_get_shape_factor
%       > get_shape_factor > find_stable_solution
%   - plot_shapefacetor_wrapper > plot_shapefacetor > find_stable_solution
% - it uses a standard input produced by gen_input_stable_sol
initialize_paths
help gen_input_stable_sol
help find_stable_solution
