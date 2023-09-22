% COPYRIGHT NOTICE 
 % This file is part of a dataset <Minar, Martin (2023), “Influence of surface energy anisotropy on nucleation and crystallographic texture of polycrystalline deposits”, Mendeley Data, V1, doi: 10.17632/bsdff8shbz.1>, coupled to publication of the same name by Minar, Moelans submitted to Computational Materials Science in September 2023. 
 % Distributed under GPLv3 license.
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
        