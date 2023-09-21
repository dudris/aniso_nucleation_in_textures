clear
% 
resultsfile = 'Results_shape_factors.mat';
load(resultsfile)

% ind =39;
ind =61;
res = results(ind);
exp_homog = 100; % between ca 10-78
T = 300; % K, temperature
k = 1.3807e-23;% Boltzmann constant, (J/K)
rc2D = exp_homog*k*T/res.A_hom/res.SLE/2; 
rc3D = sqrt(exp_homog*k*T/res.SLE); 
% res.h*rc3D
rc3D/(1.35e-10) % rc3D w.r.t. atomic radius
% rc3D/(5e-9) % rc3D w.r.t. nm


figure(5)
subplot(121)
    imagesc(res.S,'AlphaData',~isnan(res.S))
%     imagesc(res.h)
    set(gca,'YDir','normal','DataAspectRatio',[1,1,1])
    colorbar
    colormap(parula)
%     title(['S(\alpha_1,\alpha_2), ID= ' num2str(ind) ', n = ' num2str(res.n) ', \Omega=' num2str(res.Omega)])
    title(['S(\alpha_1,\alpha_2), ID= ' num2str(ind) ', n = ' num2str(res.n) ', \delta=' num2str(res.soaIE)])
subplot(122)
    P = exp(-res.S*exp_homog);
    complete_wet = ( isnan(P) & res.Cshift'>0 ) | isnan(res.Cshift');
    complete_nonwet = isnan(P) & res.Cshift'<0;
%     P(complete_wet) = 1;
    P(complete_wet) = 0;
    P(complete_nonwet) = 0;
    imagesc(P,'AlphaData',~isnan(res.S))
%     imagesc(res.Cshift')
    set(gca,'YDir','normal','DataAspectRatio',[1,1,1])
    colorbar
    title(['P(\alpha_1,\alpha_2)=exp(-' num2str(exp_homog) 'S)'])

ID = 1:length(results);

figure(1)
    subplot(311)
        res_n = [results.n]';
        plot(ID,res_n,'o')
        hold on
        plot(ind,res_n(ind),'o','MarkerFaceColor','r')
        hold off
        ylabel('n')
        grid on
        xticks(ID)
        set(gca, 'XTickLabelRotation',45)
    subplot(312)
        res_O  = [results.Omega]';
        plot(ID,res_O,'o')
        hold on
        plot(ind,res_O(ind),'o','MarkerFaceColor','r')
        hold off
        ylabel('\Omega')
        grid on
        xticks(ID)
        set(gca, 'XTickLabelRotation',45)
    subplot(313)
        soaIE = [results.Omega]'./([results.n]'.^2-1);
        plot(ID,soaIE,'o')
        ylabel('\delta')
        grid on
        xticks(ID)
        set(gca, 'XTickLabelRotation',45)
        sumnan = arrayfun(@(x)(sum(isnan(x.S(:)))), results, 'UniformOutput', false);
        numell = arrayfun(@(x)(numel(x.S(:))), results, 'UniformOutput', false);
        sumnan = cellfun(@(x) x,sumnan);
        numell = cellfun(@(x) x,numell);
        failed = sumnan./numell>0.25;
        hold on
        plot(ID(failed),soaIE(failed),'xr')
        plot(ind,soaIE(ind),'o','MarkerFaceColor','r')
        hold off
        legend('\delta','>25% NaNs')
        

%% MC simulation specification 
MCspec.sx = 50;
MCspec.sy = 200;
MCspec.growthrule = 'loc_DE<=0';
% MCspec.growthrule = 'loc_DE<0';
MCspec.nucleation_on = true; 
MCspec.neighb_code = 'nn';
% MCspec.neighb_code = 'snn'; MCspec.neighb_snn_weight = 1/8;
% MCspec.init_code = '3pm2';
% MCspec.init_code = 'rand';
MCspec.init_code = '2grains';
MCspec.aniso = @(ang,bori) 1+res.Omega/(res.n^2-1)*cos(res.n*(ang-bori));
% MCspec.aniso = @(ang,bori) 1;
MCspec.GBE_code = 'R-S';
% MCspec.GBE_code = 'iso';
MCspec.exp_homog = exp_homog; % assumed between 20-300
MCspec.nucl_fluct.code = 'shifteduni';
MCspec.nucl_fluct.shift = exp_homog/320; 
MCspec.nucl_fluct.ampl = res.SLE*0.001;
% MCspec.nucl_fluct.code = 'none';
% 

%% initiallization
sx = MCspec.sx;
sy = MCspec.sy;
s = zeros(sy,sx,'int8');

oris = res.tori;
[numgr,s(1,:) ] = suppl_MC_initialize(oris,sx,MCspec);

gs = 2*ones(1,sx); % growth site
ind_tr = 0;
%% MC simulation
tic
run_next=true;
tt = 0;
while run_next
    tt=tt+1;
    
    % a bias to nucleate in sites with sudden change in height
    % the larger the change, the larger bias - promotes spikes
%     heightvar(1) = 1 + abs(gs(1)-gs(2));
%     heightvar(2:sx) = cumsum(1+ abs(diff(gs)))+heightvar(1);
%     growprob = heightvar/max(heightvar);
%     sitex = find(growprob>rand,1,'first');     
    
    %bias to nucleate at the lowest sites
    H = abs(gs-max(gs))+1; 
    growprob = cumsum(H);
    growprob = growprob/max(growprob);
    sitex = find(growprob>rand,1,'first');     
    
    % unbiased column selection
%     sitex = ceil(rand*sx); % random column
    
    % growth site y-position
    sitey =gs(sitex); 
    
    
    normang(1) = suppl_MC_findslope(gs,sitex);
    % nucleation PDF corresponding to bottom grain orientation, assumes that values of s are indices of res.bori
     nucl_probability = P(s(sitey-1,sitex),:); 
     % advance interface
    gs(sitex) = gs(sitex)+1;
    % recalculate local interface orientation after advance
    normang(2) = suppl_MC_findslope(gs,sitex);
    
    [s,loc_DE(tt), growth_favoured(tt)]  = suppl_MC_analyze_neighborhood(MCspec,res,sitex,sitey,s,normang,nucl_probability);
    
    % feedback on nucleation
    if ~growth_favoured(tt) % negate 
        ind_tr  = ind_tr+1;
        nucl_trial(ind_tr,1) = tt;
        nucl_trial(ind_tr,2) = sitex; 
        nucl_trial(ind_tr,3) = sitey; 
        if s(sitey,sitex)~=0
            nucl_trial(ind_tr,4) = 1;
        else
            nucl_trial(ind_tr,4) = 0;
        end
    end
        
    if s(sitey,sitex)==0 % the film has not progressed
        gs(sitex) = gs(sitex)-1;
    end
    
    run_next = gs(sitex)<=(sy-1) ;

%     if mod(tt,10)==0
%         figure(3)
%             imagesc(s,'AlphaData',s>0)
%             hold on
%             plot(gs,'s')
%             hold off
%             set(gca,'YDir','normal')
%         %         set(gca,'DataAspectRatio',[1,1,1],'YDir','normal')
%             colorbar
%         pause
%     end
%     bar(gs)
%     plot(diff(gs),'o')
end% for tt
calctime = toc;

disp(['nucleation success rate=' num2str(sum(nucl_trial(:,4)==1 )) '/' num2str(size(nucl_trial,1)) '/' num2str(tt) ])
nonedge = nucl_trial(:,2)~=1 & nucl_trial(:,2)~=sx;
disp(['non-boundary nucleation success rate=' num2str(sum(nucl_trial(:,4)==1 & nonedge)) '/' num2str(sum(nonedge)) '/' num2str(tt) ])

nucl_succ_rate = [sum(nucl_trial(:,4)==1) , size(nucl_trial,1) , tt ; ...
    sum(nucl_trial(:,4)==1 & nonedge), sum(nonedge), tt];
%%
numsims = size(results(ind).MC,1);
results(ind).MC{numsims+1,1} = datestr(clock);
results(ind).MC{numsims+1,2} = calctime;
results(ind).MC{numsims+1,3} = exp_homog;
results(ind).MC{numsims+1,4} = s;
results(ind).MC{numsims+1,5} = MCspec;
results(ind).MC{numsims+1,6} = nucl_succ_rate;


numoris = length(oris);

figure(4)
    subplot(121)
        plot(loc_DE,'o')
    subplot(122)
        nucl_event = nucl_trial(:,4)==1;
        ixn = nucl_trial(nucl_event,2);
        iyn = nucl_trial(nucl_event,3);
        ilin = sub2ind(size(s),iyn,ixn);
        on = oris(s(ilin));
        plot(nucl_trial(nucl_event,1),on,'o','MarkerFaceColor','r')

ss = s;
% s=ss;
% s(s<25&s>0) = s(s<25&s>0)+50;
% s(s~=0)=s(s~=0)-24;
figure(3)
% contour(s,50)
subplot(121)
    imagesc(s,'AlphaData',s>0)
    % set(gca,'YDir','normal')
    set(gca,'DataAspectRatio',[1,1,1],'YDir','normal')
    colorbar
subplot(122)
    imagesc(s,'AlphaData',s>0)
    % set(gca,'YDir','normal')
    set(gca,'DataAspectRatio',[1,1,1],'YDir','normal')
    colorbar
    hold on 
    plot(ixn,iyn,'xr','LineWidth',1.5)
    ix = nucl_trial(~nucl_event,2);
    iy = nucl_trial(~nucl_event,3);
    plot(ix,iy,'xk','LineWidth',1.5)
    hold off

figure(6)
%     oris = results(indMC(ind)).tori;
%     [fr_bins,fr_edges] = discretize(oris(s(1,:)),oris);
%     [lr_bins,lr_edges] = discretize(oris(s(min(gs)-1,:)),oris);
%     [fr_bins,fr_edges] = discretize(s(1,:),1:50);
%     [lr_bins,lr_edges] = discretize(s(min(gs)-1,:),1:50);
%     histogram(fr_bins,fr_edges)
%     plot(lr_bins)
    histogram(oris(s(1,:))*180/pi,oris*180/pi);
    hold on
%     histogram(lr_bins,lr_edges)
    histogram(oris(s((min(gs)-1),:))*180/pi,oris*180/pi)
    hold off
    legend('first row','last row','Location','northwest')
    set(gca,'FontSize',13)
%     xticks(0:10:50)
%     xticklabels({'0','0.2','0.4','0.6','0.8','1'})
    xlabel('orientation (deg)')
    ylabel('count')

figure(7)
    bdries = diff(s,[],2)~=0;
    h = 1:min(gs);
    numgrs_h = sum(bdries,2);
    numgrs_h = numgrs_h(h);
    plot(h,numgrs_h,'o')
    title('num. of grains through thickness')
    
    
%     imagesc(s,'AlphaData',0.6*double(s>0))
%     % set(gca,'YDir','normal')
%     set(gca,'DataAspectRatio',[1,1,1],'YDir','normal')
%     colorbar
%     hold on
%     contour(bdries,1)
%     hold off
    
    

%%
% save(resultsfile,'results')

%% test the bias to grow in deep growth sites
% gs = ceil(rand(50,1)*5);
% H = abs(gs-max(gs))+1; 
% growprob = cumsum(H);
% growprob = growprob/max(growprob);
% bar(gs);
% hold on
% plot(H,'k')
% hold off
