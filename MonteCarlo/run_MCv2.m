%% results = run_MCv2(MCspec,plotting) 
function results = run_MCv2(MCspec,plotting)

load(MCspec.sf_resultsfile,'results_SF')
res = results_SF(MCspec.sf_ind);

init_msg(MCspec)

if strcmp(MCspec.nucl_mode,'aniso')
    S = res.S;
    S(~res.stab_cond) = nan;
elseif strcmp(MCspec.nucl_mode,'uniform')
    if isnan(MCspec.nucl_mode_uniform_S)
        S = get_iso_SF(res);
    else
        S = MCspec.nucl_mode_uniform_S*ones(size(res.S));
    end
end
P = exp(-S*MCspec.exp_homog);
P = P - diag(nan(size(P,1),1),0);


%% plot shape factor and prob. matrix
if plotting
    figure(5)
    ang= res.tori*180/pi;
    subplot(121)
        imagesc(ang,ang,S,'AlphaData',~isnan(S))
    %     imagesc(res.h)
        set(gca,'YDir','normal','DataAspectRatio',[1,1,1])
        colorbar
        colormap(parula)
    %     title(['S(\alpha_1,\alpha_2), ID= ' num2str(ind) ', n = ' num2str(res.n) ', \Omega=' num2str(res.Omega)])
        title(['S(\alpha_1,\alpha_2), n = ' num2str(res.n) ', \delta=' num2str(res.soaIE)])
    subplot(122)
        imagesc(ang,ang,P,'AlphaData',~isnan(S))
    %     imagesc(res.Cshift')
        set(gca,'YDir','normal','DataAspectRatio',[1,1,1])
        colorbar
        title(['P(\alpha_1,\alpha_2)=exp(-' num2str(MCspec.exp_homog) 'S)'])
end % if plotting
%% initiallization
sx = MCspec.sx;
sy = MCspec.sy;
s = zeros(sy,sx,'int8');

oris = res.tori;

resentry = cell(0,1);

for rep=1:MCspec.reps
    disp(['Repetition ' num2str(rep) '/' num2str(MCspec.reps)])
    [numgr,s(1,:) ] = suppl_MC_initialize(oris,sx,MCspec);

    gs = 2*ones(1,sx); % growth site
    ind_tr = 0;
    nucl_trial = nan(1,4);
    % % MC simulation
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

        % nucleation PDF corresponding to bottom grain orientation, assumes that values of s are indices of res.bori
         nucl_probability = P(:,s(sitey-1,sitex)); 

        if strcmp(MCspec.intf_normal_code,'fix')
            normang = pi/2*[1, 1];
            % advance interface
            gs(sitex) = gs(sitex)+1;
        elseif strcmp(MCspec.intf_normal_code,'findslope')
            normang(1) = suppl_MC_findslope(MCspec,gs,sitex);
            % advance interface
            gs(sitex) = gs(sitex)+1;
            % recalculate local interface orientation after advance
            normang(2) = suppl_MC_findslope(MCspec,gs,sitex);
        end

        [s,loc_DE(tt), growth_favoured(tt)]  = suppl_MC_analyze_neighborhood(MCspec,res,sitex,sitey,s,normang,nucl_probability);

        if s(sitey,sitex)==0 % the film has not progressed
            gs(sitex) = gs(sitex)-1;
        end

        run_next = gs(sitex)<=(sy-2) ;

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

    resentry = process_results(resentry,rep,MCspec,s,nucl_trial,tt,calctime);
    
    if plotting
        plot_all(s,loc_DE,nucl_trial,oris)
    end
        
    clear nucl_trial calctime
    s = zeros(sy,sx,'int8');
end % repetitions

if ~isempty(MCspec.outfile)
%     addpath '..\..\..\results_processing\'
    results = AddResultsEntry(MCspec.outfile,resentry,MCspec.comment);
    save(MCspec.outfile,'results')
else
    results = resentry;
end
% nonedge = nucl_trial(:,2)~=1 & nucl_trial(:,2)~=sx;
% disp(['non-boundary nucleation success rate=' num2str(sum(nucl_trial(:,4)==1 & nonedge)) '/' num2str(sum(nonedge)) '/' num2str(tt) ])


% nucl_succ_rate = [sum(nucl_trial(:,4)==1) , size(nucl_trial,1) , tt ]


end % func

%%
function init_msg(MCspec)

disp('MCv2 launched')
disp(['nfold=' num2str(MCspec.n) ', soa=' num2str(MCspec.soaIE) ', exp_homog=' num2str(MCspec.exp_homog)])
disp(['IC: ' MCspec.init_code])
disp(['BCs: ' MCspec.BCs])
disp(['nucleation on: ' num2str(MCspec.nucleation_on) ', fluctuation code: ' MCspec.nucl_fluct.code])
disp(['interface normal code: ' MCspec.intf_normal_code])
end
%%
function resentry = process_results(resentry,kk,MCspec,s,nucl_trial,tt,calctime)
    
    if kk==1
        resentry = MCspec;
    end
    
    if MCspec.save_deposit
        resentry.deposit{kk} = s;
        resentry.nucl_summary.events_txy{kk} = nucl_trial(nucl_trial(:,4)==1,1:3) ;
    end
    
    resentry.nucl_summary.events(kk) = sum(nucl_trial(:,4)==1);
    resentry.nucl_summary.trials(kk) = size(nucl_trial,1);
    resentry.MCSs(kk) = tt;
    resentry.calctime(kk) = calctime;
    
    if strcmp(MCspec.init_code,'2grains') && ~MCspec.nucleation_on
        resentry.ori_inds = [s(1,1),s(1,end)];
        thickness = size(s,1) - max(sum(s==0));
        vol1 = sum(s==s(1,1),2);
        vol2 = sum(s==s(1,end),2);
        resentry.volfrac_h(:,kk) = vol1./(vol1+vol2);
        resentry.volfrac_h((thickness+1):end,kk)=nan;
        
    elseif strcmp(MCspec.init_code,'3pm2')
        
    else
%         warning('run_MCv2: undefined output settings.')
    end % if init_code

end % func

%% 
    function plot_all(s,loc_DE,nucl_trial,oris)
    numoris = length(oris);
    
%     figure(4)
%         subplot(121)
%             plot(loc_DE,'o')
%             title('loc \Delta E ')
%             xlabel('MCSs')
%         subplot(122)
            nucl_event = nucl_trial(:,4)==1;
            ixn = nucl_trial(nucl_event,2);
            iyn = nucl_trial(nucl_event,3);
            ilin = sub2ind(size(s),iyn,ixn);
            on = oris(s(ilin))*180/pi;
%             plot(nucl_trial(nucl_event,1),on,'o','MarkerFaceColor','r')
%             ylabel('tori (deg)')
%             xlabel('MCSs')

    figure(3)
    % contour(s,50)
    subplot(121)
        imagesc(s,'AlphaData',s>0)
        % set(gca,'YDir','normal')
        set(gca,'DataAspectRatio',[1,1,1],'YDir','normal')
        colorbar('Limits',[0,50])
        title('Deposit')
    subplot(122)
        imagesc(s,'AlphaData',s>0)
        % set(gca,'YDir','normal')
        set(gca,'DataAspectRatio',[1,1,1],'YDir','normal')
        colorbar
        colorbar('Limits',[0,50])
        hold on 
        plot(ixn,iyn,'o','LineWidth',1,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,0,0],'MarkerSize',5)
%         ix = nucl_trial(~nucl_event,2);
%         iy = nucl_trial(~nucl_event,3);
%         plot(ix,iy,'xk','LineWidth',1.5)
        hold off
        legend('successful')
        title('Deposit with successful nucleation events')
    
    gs = sum(s~=0);
        
%     figure(6)
%     %     oris = results(indMC(ind)).tori;
%     %     [fr_bins,fr_edges] = discretize(oris(s(1,:)),oris);
%     %     [lr_bins,lr_edges] = discretize(oris(s(min(gs)-1,:)),oris);
%     %     [fr_bins,fr_edges] = discretize(s(1,:),1:50);
%     %     [lr_bins,lr_edges] = discretize(s(min(gs)-1,:),1:50);
%     %     histogram(fr_bins,fr_edges)
%     %     plot(lr_bins)
%         histogram(oris(s(1,:))*180/pi,oris*180/pi);
%         hold on
%     %     histogram(lr_bins,lr_edges)
%         histogram(oris(s((min(gs)-1),:))*180/pi,oris*180/pi)
%         hold off
%         legend('first row','last row','Location','northwest')
%         set(gca,'FontSize',13)
%     %     xticks(0:10:50)
%     %     xticklabels({'0','0.2','0.4','0.6','0.8','1'})
%         xlabel('orientation (deg)')
%         ylabel('count')

    figure(7)
        bdries = diff(s,[],2)~=0;
        h = 1:min(gs);
        numgrs_h = sum(bdries,2);
        numgrs_h = numgrs_h(h);
        plot(h,numgrs_h,'o')
        xlabel('thickness')
        ylabel('num. of grains')
        title('num. of grains through thickness evolution')

end % func