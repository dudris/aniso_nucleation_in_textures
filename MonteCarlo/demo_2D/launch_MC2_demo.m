clear
% 
sf_resultsfile = 'shape_factor\results_SF.mat';
% __ specify shape factor-orientation map to be used (uncomment the
%   respective order of symmetry below)
SFIND = 1; % 0.1 ; 4-fold, 50x50
% SFIND = [1,3,5,7]; % 0.1,0.3,0.5,0.7 ; 4-fold, 50x50
% SFIND = [15,17,19,21]; % 0.1,0.3,0.5,0.7 ; 6-fold, 50x50

%__ uncomment below to check the shape factor data
% load(sf_resultsfile,'results_SF')
% A = [(1:length(results_SF))' , [results_SF.n]', [results_SF.soaIE]', cellfun(@(x) size(x,1),{results_SF.S}','UniformOutput',true)];
% results_overview = array2table(A,'VariableNames',{'index','n','anisotropy_strength','map_resolution'})

% __ homogeneous nucleation barriers to run the simulation for
EXPHOM = 10;
% EXPHOM = [10, 50, 90];


% __ runs through all combinations of SFIND and EXPHOM
% __ the particular settings are always locally modified by MCspec_01,
%    MCspec_02 and MCspec_03 depending on the case to have it right
% __ ADJUST THE BELOW SPECS TO CONTROL THE INITIAL CONDITION AND OUTPUT
%   - MCspec.init_code
%   - MCspec.outfile
for ee = 1:length(EXPHOM)
    exp_homog = EXPHOM(ee);
    for soa = 1:length(SFIND)
        sf_ind = SFIND(soa);
        % __ MC simulation specification 
        MCspec = gen_MCspec(sf_resultsfile,sf_ind);
        MCspec.sx = 1500;
        MCspec.sy = 200;
        MCspec.sx = 50;
        MCspec.sy = 50;
        SIM = {'no_nucleation', 'iso_nucleation','aniso_nucleation'};
        run_simtypes = logical([1,1,1]);
        simtypes = SIM(run_simtypes);
        reps  = 10;
        MCspec.save_deposit = true;
        MCspec.reps = reps; % number of repetitions
        MCspec.init_code = 'rand_31pts_edges';
        MCspec.outfile = 'MonteCarlo\demo_2D\results_MC2_demo.mat';
        is_simtype = @(simtype,simtypes) cellfun(@(x) strcmp(x,simtype),simtypes);

        plotting = false;

        % __ no_nucleation
        if run_simtypes(1) 
            simtype = SIM{1};
            MCspec.exp_homog = exp_homog;
            MCspec.reps = reps; % number of repetitions
            MCspec = MCspec_01_no_nucl(MCspec);
            MCspec.comment = [simtype ', seed: ' MCspec.init_code ', soaIE = ' num2str(MCspec.soaIE) ', ' num2str(MCspec.sx) 'x' num2str(MCspec.sy) ', exp hom=' num2str(exp_homog)];
            disp(['Launching ''' simtype ''' simulation with ' num2str(reps) ' repetitions.' ])
            results = run_MCv2(MCspec,plotting);
        end

        if run_simtypes(2) 
            simtype = SIM{2};
            MCspec.exp_homog = exp_homog;
            MCspec.reps = reps; % number of repetitions
            MCspec = MCspec_02_iso_nucl(MCspec);
            MCspec.comment = [simtype ', seed: ' MCspec.init_code ', soaIE = ' num2str(MCspec.soaIE) ', ' num2str(MCspec.sx) 'x' num2str(MCspec.sy) ', exp hom=' num2str(exp_homog) ', RS GB iso nucl.'];
            disp(['Launching ''' simtype ''' simulation with ' num2str(reps) ' repetitions.' ])
%             disp(MCspec.comment)
            results = run_MCv2(MCspec,plotting);
        end


        if run_simtypes(3) 
            simtype = SIM{3};
            MCspec.exp_homog = exp_homog;
            MCspec.reps = reps; % number of repetitions
            MCspec = MCspec_03_aniso_nucl(MCspec);
            MCspec.comment = [simtype ', seed: ' MCspec.init_code ', soaIE = ' num2str(MCspec.soaIE) ', ' num2str(MCspec.sx) 'x' num2str(MCspec.sy) ', exp hom=' num2str(exp_homog)];
            disp(['Launching ''' simtype ''' simulation with ' num2str(reps) ' repetitions.' ])
            results = run_MCv2(MCspec,plotting);
        end
    end% soa
end%  exp_homog




    

