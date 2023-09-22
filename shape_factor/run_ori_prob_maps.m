% COPYRIGHT NOTICE 
 % This file is part of a dataset <Minar, Martin (2023), “Influence of surface energy anisotropy on nucleation and crystallographic texture of polycrystalline deposits”, Mendeley Data, V1, doi: 10.17632/bsdff8shbz.1>, coupled to publication of the same name by Minar, Moelans submitted to Computational Materials Science in September 2023. 
 % Distributed under GPLv3 license.
%%
 clear 
 
resultsfile = [mfilename('fullpath') '\dummy_results_SF.mat'];
% __ change to true to cumulatively add computed map
addtoresults = true;

% __ standard user input 
soaIE = 0.7; % strength of anisotropy 0-1
n =2; % nfold - choose 2, 3, 4, or 6
numoris = 50; % number of points in each direction in the map

% __ other parameters in the input
in.misor_0 = 15/180*pi; 
in.SLE = 0.3; % scalar solid-liquid interface energy
in.GBE = 0.3; % scalar grain boundary energy
in.n = n;
in.tori = linspace(0,2*pi/in.n,numoris+1)'; % tori ... top grain orientation
in.tori(end) = [];
in.bori = in.tori; % bori ... bottom grain orientation; take a symmetric map
in.Omega = soaIE*(in.n^2-1); % normalized strength of anisotropy
disp(['delta=' num2str(soaIE) ', n=' num2str(in.n)])
tic
results = get_shape_factor_oridep(in,{addtoresults,resultsfile});
toc





