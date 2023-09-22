% COPYRIGHT NOTICE 
 % This file is part of a dataset <Minar, Martin (2023), “Influence of surface energy anisotropy on nucleation and crystallographic texture of polycrystalline deposits”, Mendeley Data, V1, doi: 10.17632/bsdff8shbz.1>, coupled to publication of the same name by Minar, Moelans submitted to Computational Materials Science in September 2023. 
 % Distributed under GPLv3 license.
%% normang = suppl_MC_findslope(MCspec,growth_sites,IND) 
% - take into accounr 'numng' points around the one where the slope is
% seeked (midpoint)
% - nearest neighbors are taken with weight 'weightnn'
% - weighted mean difference value on the left and right from the center point is
% taken and the values on the left and right are connected by straight line
% - that is taken to be the tangent to interface in the midpoint
% - assumes unit spacing in both x and y direction

function normang = suppl_MC_findslope(MCspec,growth_sites,IND)
numng = 2;
weightnn = 3;

numpts = length(growth_sites);

% numng points to the left and right
lind = (IND-numng):(IND-1);
rind = (IND+1):(IND+numng);

% is near system boundary?
beyond_LB = lind<1;
beyond_RB = rind>numpts;
% discard indices which are outside of the system
if any(beyond_LB) % near left boundary
    if strcmp(MCspec.BCs,'no_interaction')
        normang = pi/2;
        return
        
    elseif strcmp(MCspec.BCs,'periodic')
        sum_beyond_LB = sum(beyond_LB);
        lind(beyond_LB) = numpts- (sum_beyond_LB:-1:1) + 1; % 
    end % switch BC
    
elseif any(beyond_RB) % near right boundary
    if strcmp(MCspec.BCs,'no_interaction')
        normang = pi/2;
        return
        
    elseif strcmp(MCspec.BCs,'periodic')
        sum_beyond_RB = sum(beyond_RB);
        rind(beyond_RB) = 1:sum_beyond_RB; % 
    end % switch BC
    
end % if near boundary


if ~isempty(lind)
    lw = ones(size(lind'));
    lw(end) = weightnn;
    ldiffs = growth_sites(lind)-growth_sites(IND);
    % weighted mean 
    vl = ldiffs*lw/sum(lw);
    xvl = mean(lind);
else % redundant
    % center the left point in the edge growth site
    vl = 0; 
    xvl = IND;
end

if ~isempty(rind)
    rw = ones(size(rind'));
    rw(1) = weightnn;
    rdiffs = growth_sites(rind)-growth_sites(IND);
    % weighted mean 
    vr = rdiffs*rw/sum(rw);
    xvr = mean(rind);
else % redundant
    % center the right point in the edge growth site
    vr = 0;
    xvr = IND;
end

slope = (vr-vl)/(xvr-xvl);
tangang = atan2(vr-vl,xvr-xvl);
% tangang = atan(slope)*180/pi;
% if tangang<0
%     normang = tangang +pi/2;
% else 
%     normang = tangang +pi/2;
% end
normang = tangang+pi/2;

% if 
% normang = atan(slope)*180/pi+90;


plotting = false;
if plotting
    id = 1:numpts;
    sites = [lind,IND,rind];
    
    bar(id(sites),growth_sites(sites))
    grid on
    hold on
    plot(id(sites),slope*(id((sites))-IND)+growth_sites(IND))
    hold off
    axis equal
    pause
end % plotting

end %end func