% COPYRIGHT NOTICE 
 % This file is part of a dataset <Minar, Martin (2023), “Influence of surface energy anisotropy on nucleation and crystallographic texture of polycrystalline deposits”, Mendeley Data, V1, doi: 10.17632/bsdff8shbz.1>, coupled to publication of the same name by Minar, Moelans submitted to Computational Materials Science in September 2023. 
 % Distributed under GPLv3 license.
%% MC simulation specification 

function MCspec = MCspec_2grains(MCspec)
   
    MCspec.sx = 50;
    MCspec.sy = 100;
    %
    MCspec.BCs = 'no_interaction'; 
    %
    MCspec.nucleation_on = false; 
    %
%     MCspec.neighb_code = 'nn';
    MCspec.neighb_code = 'snn'; 
%     MCspec.neighb_code = 'tnn'; MCspec.neighb_snn_weight = 1/4; MCspec.neighb_tnn_weight = 1/10;
    %
    MCspec.init_code = '2grains';
    %
    MCspec.intf_normal_code = 'fix';
    %
    MCspec.GBE_code = 'iso';
    %
    MCspec.nucl_fluct.code = 'none';
    %
    MCspec.reps = 30; % number of repetitions
    %
    MCspec.nucl_mode = 'uniform'; % 'uniform' or 'aniso'
    MCspec.nucl_tres_nobarr = 1; % probability above which no barrier is assumed and other oris are disregarded
    %
    MCspec.save_deposit = false;
    %
end% func