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