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

