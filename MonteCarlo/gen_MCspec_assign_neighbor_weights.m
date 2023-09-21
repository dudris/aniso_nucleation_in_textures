function MCspec = gen_MCspec_assign_neighbor_weights(MCspec)
    if MCspec.n == 3
        if strcmp(MCspec.neighb_code,'snn')
            MCspec.neighb_snn_weight = 1/4.3;
        elseif strcmp(MCspec.neighb_code,'nn')
        elseif strcmp(MCspec.neighb_code,'tnn')
            error('Weights for neighborcode tnn were not optimized' )
        end% if neighbor_code
        
    elseif MCspec.n == 4
        if strcmp(MCspec.neighb_code,'snn')
            MCspec.neighb_snn_weight = 1/4.4;
        elseif strcmp(MCspec.neighb_code,'nn')
        elseif strcmp(MCspec.neighb_code,'tnn')
            error('Weights for neighborcode tnn were not optimized' )
        end% if neighbor_code
        
    elseif MCspec.n == 6
         if strcmp(MCspec.neighb_code,'snn')
            MCspec.neighb_snn_weight = 1/4.4; % might be optimized but works well enough
        elseif strcmp(MCspec.neighb_code,'nn')
        elseif strcmp(MCspec.neighb_code,'tnn')
            error('Weights for neighborcode tnn in 6-fold were not optimized' )
        end% if neighbor_code
        
    end % if n
end% func
        