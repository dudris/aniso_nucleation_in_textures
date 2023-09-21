%% MCspec = gen_MCspec(sf_resultsfile,sf_ind)
% - generate standard input to the Monte Carlo algorithm in 'run_MCv2'
% - loads existing shape-factor orientation map and assumes that anisotropy
% - specifies system dimensions, initial condition, nucleation behavior and simulation output handling
% - INPUT
%   - sf_resultsfile ... path to the *.mat file with the shape
%      factor-orientation maps as produced by 'run_ori_prob_maps.m'
%   - sf_ind                ... index of the shape factor-orientation map in the saved structure
% - OUTPUT
%   - MCspec            ... structure containing standard input to the Monte Carlo
%       algorithm in 'run_MCv2.m'

function MCspec = gen_MCspec(sf_resultsfile,sf_ind)

% __ INITIALIZE THE INPUT STRUCTURE
    MCspec = struct;
    
% __ DIMENSIONS OF THE SYSTEM  (growth in positive Y direction)
    MCspec.sx = 300;
    MCspec.sy = 200;
    
% __ENABLE/DISABLE NUCLEATION
    MCspec.nucleation_on = true; 
    
% __ NUCLEATION MODE: 
%       'uniform'... nucleation with ISOTROPIC interface energy (irrespective of the loaded shape factor-orientation map)
%       'aniso'     ... nucleation with ANISOTROPIC interface energy (using the shape loaded factor-orientation map)
    MCspec.nucl_mode = 'aniso'; % 'uniform' or 'aniso'

% LOAD THE ANISOTROPY FROM THE SHAPE FACTOR-ORIENTATION MAP AND 
% SET DEFAULT SETTINGS FOR THE MONTE CARLO ALGORITHM
    MCspec = set_default_specs_and_anisotropy(MCspec, sf_resultsfile,sf_ind);
    
% __ THE INITIAL CONDITION OF THE SEED ROW - choose one of the options
%   - each seed pixel to have uniformly distributed random value:
%       - 'rand' ...  1-50
%       - 'rand_31pts_center' ... 11-40
%       - 'rand_31pts_edges' ... 1-15 united with 36-50
%   - '2grains' ... 2 grains with maximal and minimal solid-liquid
%       interface energy (used for validation) 
    MCspec.init_code = 'rand'; 

% __ NON-DIMENSIONAL NUCLEATION BARRIER (\beta in the paper)
    MCspec.exp_homog = 50; % assumed between 20-300

% __ LOCAL ENERGY FLUCTUATIONS - turn off with 'none'
    MCspec.nucl_fluct.code = 'shifteduni'; % valid options: 'shifteduni', 'none'
    if strcmp(MCspec.nucl_fluct.code,'shifteduni')
        MCspec.nucl_fluct.shift = 0.25; 
        MCspec.nucl_fluct.ampl = MCspec.SLE*0.001;
    end
  
% __ NUMBER OF SIMULATION REPETITIONS
    MCspec.reps = 1; 

% __ OUTPUT CONTROL AND COMMENT
    MCspec.save_deposit = false; % to (not) save the simulated deposit to the output data 
    MCspec.outfile = ''; % when empty, no record of the simulation is saved ; 
%     MCspec.outfile = 'results_MCv2_dev_2grains.mat'; % save there
    MCspec.comment = ''; % comment to be saved to the simulation run in the output structure
    
end% func

function MCspec = set_default_specs_and_anisotropy(MCspec, sf_resultsfile,sf_ind)
% __SPECIFY THE FILE WITH THE SHAPE FACTOR-ORIENTATION MAP (irrelevant in case of nucleation with isotropic interface energy)
% __ ANISOTROPY IN INTERFACE ENERGY READ AS WELL
    MCspec.sf_resultsfile = sf_resultsfile;
    MCspec.sf_ind = sf_ind;
    load(sf_resultsfile,'results_SF')
%     load(sf_resultsfile)
%     load(which(sf_resultsfile),'results')
    res = results_SF(sf_ind);
    MCspec.n = res.n;
    MCspec.Omega = res.Omega;
    MCspec.soaIE = res.soaIE;
    stringaniso = ['@(ang,bori) 1+' num2str(MCspec.soaIE) '*cos(' num2str(MCspec.n) '*(ang-bori))'];
    MCspec.aniso = str2func(stringaniso);
    
    MCspec.GBE = res.GBE; % (J/m2), grain boundary energy
    MCspec.SLE = res.SLE; % (J/m2), solid-liquid interface energy

    %     MCspec.intf_normal_code = 'findslope';
    MCspec.intf_normal_code = 'fix';
    %
    %     MCspec.neighb_code = 'nn';
    MCspec.neighb_code = 'snn'; 
    %     MCspec.neighb_code = 'tnn';
    MCspec = gen_MCspec_assign_neighbor_weights(MCspec);
    %
    % MCspec.GBE_code = 'R-S';
    MCspec.GBE_code = 'iso';
    %
    %  __ set MCspec.nucl_mode_uniform_S
    MCspec = gen_MCspec_assign_uniform_S(MCspec); % in optional second argument the uniform shape factor can be specified
    MCspec.nucl_tres_nobarr = 1; % probability above which no barrier is assumed and other oris are disregarded
    %
    MCspec.BCs = 'no_interaction'; 
    % MCspec.BCs = 'periodic'; 
    %
    %     MCspec.growthrule = 'loc_DE<=0';
    MCspec.growthrule = 'loc_DE>=0';
    % MCspec.growthrule = 'loc_DE<0';
end



