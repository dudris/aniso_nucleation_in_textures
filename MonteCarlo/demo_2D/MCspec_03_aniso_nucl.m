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