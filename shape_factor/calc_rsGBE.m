% COPYRIGHT NOTICE 
 % This file is part of a dataset <Minar, Martin (2023), “Influence of surface energy anisotropy on nucleation and crystallographic texture of polycrystalline deposits”, Mendeley Data, V1, doi: 10.17632/bsdff8shbz.1>, coupled to publication of the same name by Minar, Moelans submitted to Computational Materials Science in September 2023. 
 % Distributed under GPLv3 license.
% rsGBE = calc_rsGBE(GBE,misor_0,misor)
% test
% misor = linspace(-1,1,120);
% rsGBE = calc_rsGBE(1,0.3,misor);
% plot(misor,rsGBE,'o')

function rsGBE = calc_rsGBE(GBE,misor_0,misor)

    assert(misor_0>=0)

    rsGBE = zeros(size(misor));
    cond1 = abs(misor)>=misor_0;
    rsGBE(cond1) = GBE;
    
    x = abs(misor)/misor_0;
    cond2 = abs(misor)<misor_0;
    rsGBE(cond2) = GBE*x(cond2).*(1-log(x(cond2)));
    
%     rsGBE(isnan(rsGBE)) = 0;
    
end%func

