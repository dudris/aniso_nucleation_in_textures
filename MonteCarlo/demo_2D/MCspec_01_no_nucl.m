% COPYRIGHT NOTICE 
 % This file is part of a dataset <Minar, Martin (2023), “Influence of surface energy anisotropy on nucleation and crystallographic texture of polycrystalline deposits”, Mendeley Data, V1, doi: 10.17632/bsdff8shbz.1>, coupled to publication of the same name by Minar, Moelans submitted to Computational Materials Science in September 2023. 
 % Distributed under GPLv3 license.
function MCspec = MCspec_01_no_nucl(MCspec)

MCspec.neighb_code = 'snn';
MCspec = gen_MCspec_assign_neighbor_weights(MCspec);
MCspec.nucleation_on = false; 
MCspec.GBE_code = 'iso';
MCspec.nucl_fluct.code = 'none';

if strcmp(MCspec.neighb_code,'snn') && isfield(MCspec,'neighb_tnn_weight')
    MCspec  = rmfield(MCspec,'neighb_tnn_weight');
end

end