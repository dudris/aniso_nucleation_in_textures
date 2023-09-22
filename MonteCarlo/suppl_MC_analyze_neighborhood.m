% COPYRIGHT NOTICE 
 % This file is part of a dataset <Minar, Martin (2023), “Influence of surface energy anisotropy on nucleation and crystallographic texture of polycrystalline deposits”, Mendeley Data, V1, doi: 10.17632/bsdff8shbz.1>, coupled to publication of the same name by Minar, Moelans submitted to Computational Materials Science in September 2023. 
 % Distributed under GPLv3 license.
function [s,loc_DEmin,growth_favoured] = suppl_MC_analyze_neighborhood(MCspec,res,sitex,sitey,s,normang,nucl_probability)
    
    oris = res.tori;
    loc_totE1 = 0;
    
    indlist = get_neighbor_indices(sitex,sitey,MCspec);
%     indlist = get_neighbor_indices(sitex,sitey,size(s,2),'nearest');
%     indlist = get_neighbor_indices(sitex,sitey,size(s,2),'second_nearest');
    
    siteval = s(sitey,sitex);
    bori1 = oris(s(sitey-1,sitex));
    [neighb(:,1),is_equal_val(:,1),is_surf(:,1),is_GB(:,1),loc_totE1] = ...
        check_nighborhood(MCspec,res,indlist,siteval,s,bori1,normang(1));
    
    % unique nearest neighbors
    uniquenbs = unique(neighb(indlist.nn,1));
    uniquenbs(uniquenbs==0)=[]; % leave the liquid out
        
    % energies when a neighbor is promoted to the site
    loc_totE2 = zeros(size(uniquenbs));
    
    % get energies when nearest-neighbor-value is added to site
    for unn = 1:length(uniquenbs)
        iind = unn+1;
        bori2 = oris(uniquenbs(unn));
        [~,is_equal_val(:,iind),is_surf(:,iind),is_GB(:,iind),loc_totE2(unn,1)] = ...
            check_nighborhood(MCspec,res,indlist,uniquenbs(unn),s,bori2,normang(2));
    end % for unique neighbors
    
    loc_DE = (loc_totE2-loc_totE1);
    
%     if any(loc_DE<0) 
%     if any(any(is_GB))
%         disp(['loc_totE1=[' num2str(loc_totE1,'%4.3e,') '], loc_totE2=[' num2str(loc_totE2','%4.3e,') ']'  ])
%         disp(num2str(indlist.ind))
%         imagesc(s,'AlphaData',s>0)
%         set(gca,'YDir','normal')
% %         ylim([0,10])
%         [MCspec.aniso(normang,oris(uniquenbs(1))), MCspec.aniso(normang,oris(uniquenbs(2)))]
%     end
    
    loc_DE = add_fluctuation(loc_DE,MCspec);

%     disp(['DeltaE=' num2str(loc_DE')])
    if strcmp(MCspec.growthrule,'loc_DE<=0')
        min_vals = loc_DE==max(loc_DE);
    elseif strcmp(MCspec.growthrule,'loc_DE>=0')
        min_vals = loc_DE==min(loc_DE);
    end
    
    if sum(min_vals)>1 % multiple options energetically equivalent
        % choose the one which is more frequent in second-nearest neighborhood
        % value 0 from North left out => 3 sites around
        
        counts = sum(neighb==uniquenbs',1); % row vector, counts of neighbors 
        mostfreq = counts'==max(counts(min_vals));
        subselected_ind = find(mostfreq);
        
        if length(subselected_ind)==1 % unambguously selected most frequent neighbor with minimal energy change
            min_val_ind = subselected_ind;
        else % must choose randomly
            randind = ceil(rand*length(subselected_ind));
            min_val_ind = subselected_ind(randind); % bool made into index (the smalles of the equivalent)
        end  % if 
            
    else
        min_val_ind = find(min_vals);
    end % if
    loc_DEmin = loc_DE(min_val_ind);
    
    growth_favoured = is_growth_favoured(MCspec.growthrule,loc_DEmin);
    
    if growth_favoured
        s(sitey,sitex) = uniquenbs(min_val_ind);
    else
        if MCspec.nucleation_on
            randnr = rand;
            s(sitey,sitex) = nucleation(nucl_probability,randnr,MCspec.nucl_mode,MCspec.nucl_tres_nobarr);
        else
            s(sitey,sitex) = 0;
        end% if nucleation is on

    end
    
    
end % function

%%
function [neighb,is_equal_val,is_surf,is_GB,loc_totE] = check_nighborhood(MCspec,res,indlist,siteval,s,bori,normang)
    
    oris = res.tori;
    GBE = res.GBE;
    SLE = res.SLE;
    
    aniso = MCspec.aniso;
    loc_totE = 0;
    
    % find linear indices of the sites in the matrix as given in indlist
    lininds = sub2ind(size(s),indlist.ind(:,1),indlist.ind(:,2));
    % get neighbor values
    neighb = s(lininds);
    % assess values
    ndiff = neighb-siteval;
    % assess type of bond
    is_equal_val = ndiff==0;
    is_surf = (ndiff~=0) & (neighb==0 | siteval==0); % 
    is_GB = (ndiff~=0) & (neighb~=0 & siteval~=0);
    
    if strcmp(MCspec.growthrule,'loc_DE<=0')
        sign_ = -1;
    elseif strcmp(MCspec.growthrule,'loc_DE>=0')
        sign_ = 1;
    end
    
    
    if any(is_GB) % nearest neighbor
        
        nnGB = indlist.nn(is_GB);
        
        if strcmp(MCspec.GBE_code,'R-S') % Read-shockley
            oris_neigh = oris(neighb(is_GB));
            disori = abs ( oris_neigh - oris(siteval) );
            % __ add periodicity to R-S
            ppp = disori>pi/in.n;
            disori(ppp) = -pi/in.n + (disori(ppp)-pi/in.n);
            loc_GB = calc_rsGBE(GBE,15/180*pi,disori);
            loc_totE = loc_totE + sign_*sum(loc_GB(nnGB)); % nearest neighbors
            if strcmp(MCspec.neighb_code,'snn') % nearest neighbors AND second nn
                snnGB = indlist.snn(is_GB);
                loc_totE = loc_totE + sign_*sum(loc_GB(snnGB))*MCspec.neighb_snn_weight; % second nearest neighbors
            end % if snn
            if strcmp(MCspec.neighb_code,'tnn') % nearest neighbors AND second nn
                tnnGB = indlist.tnn(is_GB);
                loc_totE = loc_totE + sign_*sum(loc_GB(tnnGB))*MCspec.neighb_tnn_weight; % third nearest neighbors
            end % if tnn
            
            
        elseif strcmp(MCspec.GBE_code,'iso') % isotropic
            loc_totE = loc_totE + sign_*sum(nnGB)*GBE;
            if strcmp(MCspec.neighb_code,'snn') % nearest neighbors AND second nn
                snnGB = indlist.snn(is_GB);
                loc_totE = loc_totE + sign_*sum(snnGB)*GBE*MCspec.neighb_snn_weight;
            end % if snn
            if strcmp(MCspec.neighb_code,'tnn') % nearest neighbors AND second nn
                tnnGB = indlist.tnn(is_GB);
                loc_totE = loc_totE + sign_*sum(tnnGB)*GBE*MCspec.neighb_tnn_weight;
            end % if tnn
        end %  spec GBE_code
        
    end % if is_GB
    
    if any(is_surf)
        % __ the orientation of the bottom pixel assumed in the anisotropy
        loc_anisoIE = SLE*aniso(normang,bori);
%         % __ option1 ) anisoIE only added once
        loc_totE = loc_totE + sign_*loc_anisoIE;
        % __ option2) every surface pixel around site adds a dangling bond
        % - NOT WORKING
%         nn_surf = indlist.nn(is_surf);
%         loc_totE = loc_totE + sign_*sum(nn_surf)*loc_anisoIE;
% %         if strcmp(MCspec.neighb_code,'snn') % nearest neighbors AND second nn
% %             snn_surf = indlist.snn(is_surf);
% %             loc_totE = loc_totE + sign_*sum(snn_surf)*loc_anisoIE*MCspec.neighb_snn_weight;
% %         end % if snn
% %         if strcmp(MCspec.neighb_code,'tnn') % nearest neighbors AND second nn
% %             tnn_surf = indlist.tnn(is_surf);
% %             loc_totE = loc_totE + sign_*sum(tnn_surf)*loc_anisoIE*MCspec.neighb_tnn_weight;
% %         end % if tnn
        
    end % if is_surf
    
end% function

%%
function indlist = get_neighbor_indices(sitex,sitey,MCspec)
    
    sizeS = MCspec.sx;
    switch MCspec.neighb_code
        case 'nn'
            neighbor_order = 1;
        case 'snn'
            neighbor_order = 2;
        case 'tnn'
            neighbor_order = 3;
    end
    
    if neighbor_order>=1 % nearest neighbors
        indlist.ind = [sitey,sitex+1; ... % East
                        sitey+1,sitex; ... % North
                        sitey,sitex-1; ... % West
                        sitey-1,sitex]; % South
        indlist.nn = true(4,1);
        
        is_neigh_off_W = indlist.ind(:,2)==0; % West neighbor off-domain
        is_neigh_off_E = indlist.ind(:,2)==sizeS+1; % East neighbor off-domain
    end % nearest neighbors
    
    if neighbor_order>=2 % second nearest neighbors
        indlist.ind = [indlist.ind;
                        sitey+1,sitex+1; ... % NorthEast
                        sitey+1,sitex-1; ... % NorthWest
                        sitey-1,sitex-1; ... % SouthWest
                        sitey-1,sitex+1; ... % SouthEast
                        ]; 
        indlist.nn = [indlist.nn; false(4,1)];
        indlist.snn = [false(4,1); true(4,1)];
        
        % recompute to include off-domain snn
        is_neigh_off_W = indlist.ind(:,2)==0; % West neighbor out of domain
        is_neigh_off_E = indlist.ind(:,2)==sizeS+1; % East neighbor out of domain
    end % nearest neighbors
    
    % treat off-domain neighors
    if any(is_neigh_off_W) || any(is_neigh_off_E)
        if strcmp(MCspec.BCs,'no_interaction' )
            is_neigh_off = is_neigh_off_W | is_neigh_off_E;
            indlist.ind(is_neigh_off,:) = [];
            indlist.nn(is_neigh_off,:) = [];
            if neighbor_order>=2
                indlist.snn(is_neigh_off,:) = [];
            end

        elseif strcmp(MCspec.BCs,'periodic' )
            indlist.ind(is_neigh_off_W,2) = sizeS;
            indlist.ind(is_neigh_off_E,2) = 1;
        end % switch BCs
    end % if some neighbor is off
    
    if neighbor_order>=3 % third nearest neighbors
        indlist.ind = [indlist.ind;
                        sitey,sitex+2; ... % EastEast
                        sitey+2,sitex; ... % NorthNorth
                        sitey,sitex-2; ... % WestWest
                        sitey-2,sitex; ... % SouthSouth
                        ]; 
        indlist.nn = [indlist.nn; false(4,1)];
        indlist.snn = [indlist.snn; false(4,1)];
        indlist.tnn = [false(8,1); true(4,1)];
        
        % E,W,S off-domain tnn
        is_neigh_off_W = indlist.ind(:,2)<1; % West off-domain neighbor 
        is_neigh_off_E = indlist.ind(:,2)>sizeS; % East off-domain neighbor 
        is_neigh_out_S = indlist.ind(:,1)==0; % South off-domain neighbor  - first row
        
        % exclude S off-domain neighbor
        indlist.ind(is_neigh_out_S,:) = [];
        indlist.nn(is_neigh_out_S,:) = [];
        indlist.snn(is_neigh_out_S,:) = [];
        indlist.tnn(is_neigh_out_S,:) = [];
        
        % treat oher off-domain tnn
        if any(is_neigh_off_W) || any(is_neigh_off_E)
            if strcmp(MCspec.BCs,'no_interaction' )
                is_neigh_off = is_neigh_off_W | is_neigh_off_E;
                indlist.ind(is_neigh_off,:) = [];
                indlist.nn(is_neigh_off,:) = [];
                indlist.snn(is_neigh_off,:) = [];
                indlist.tnn(is_neigh_off,:) = [];

            elseif strcmp(MCspec.BCs,'periodic' )
                indlist.ind(indlist.ind(:,2)==sizeS+1,2) = 1;
                indlist.ind(indlist.ind(:,2)==sizeS+2,2) = 2;
                indlist.ind(indlist.ind(:,2)==0,2) = sizeS;
                indlist.ind(indlist.ind(:,2)==-1,2) = sizeS-1;
            end % switch BCs
        end
    end % nearest neighbors
    
end % func

%%
function growth_favoured = is_growth_favoured(growthrule,loc_DE)

	if strcmp(growthrule,'loc_DE<=0')
        growth_favoured = any(loc_DE<=0);
        
    elseif strcmp(growthrule,'loc_DE<0')
        growth_favoured = any(loc_DE<0);
        
    elseif strcmp(growthrule,'loc_DE>=0')
        growth_favoured = any(loc_DE>=0);
        
    end % switch growthrule
    
end % func

%%
function loc_DE = add_fluctuation(loc_DE,MCspec)
    spec = MCspec.nucl_fluct.code;
    
    if strcmp(spec, 'none')
        return
    elseif strcmp(spec, 'shifteduni')
        if strcmp(MCspec.growthrule,'loc_DE<=0')
            loc_DE = loc_DE + MCspec.nucl_fluct.ampl*(rand-1+MCspec.nucl_fluct.shift);

        elseif strcmp(MCspec.growthrule,'loc_DE<0')
            loc_DE = loc_DE + MCspec.nucl_fluct.ampl*(rand-1+MCspec.nucl_fluct.shift);

        elseif strcmp(MCspec.growthrule,'loc_DE>=0')
            loc_DE = loc_DE + MCspec.nucl_fluct.ampl*(rand-MCspec.nucl_fluct.shift);

        end % switch growthrule
% %         loc_DE = (loc_totE2-loc_totE1) + res.SLE*0.001*(rand-MCspec.exp_homog/320-0.06);
%         loc_DE = loc_DE + MCspec.nucl_fluct.ampl*(rand-MCspec.nucl_fluct.shift);
% %         loc_DE = loc_totE2-loc_totE1 + res.SLE*0.01*(-2.5 + randn);
        
    end % switch spec

end % func
