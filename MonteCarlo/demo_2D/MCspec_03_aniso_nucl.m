% COPYRIGHT NOTICE 
 % This file is part of a dataset <Minar, Martin (2023), “Influence of surface energy anisotropy on nucleation and crystallographic texture of polycrystalline deposits”, Mendeley Data, V1, doi: 10.17632/bsdff8shbz.1>, coupled to publication of the same name by Minar, Moelans submitted to Computational Materials Science in September 2023. 
 % Distributed under GPLv3 license.
function MCspec = MCspec_03_aniso_nucl(MCspec)

MCspec.neighb_code = 'snn';
MCspec = gen_MCspec_assign_neighbor_weights(MCspec);
MCspec.nucleation_on = true; 
MCspec.GBE_code = 'iso';

MCspec.nucl_fluct.code = 'shifteduni';
% MCspec.nucl_fluct.shift = 1e-2; 
MCspec.nucl_fluct.shift = 0.25; 
MCspec.nucl_fluct.ampl = MCspec.SLE*0.001;

MCspec.nucl_mode = 'aniso'; % 'uniform' or 'aniso'
MCspec.nucl_tres_nobarr = 1; % probability above which no barrier is assumed and other oris are disregarded

if strcmp(MCspec.neighb_code,'snn') && isfield(MCspec,'neighb_tnn_weight')
    MCspec  = rmfield(MCspec,'neighb_tnn_weight');
end

end