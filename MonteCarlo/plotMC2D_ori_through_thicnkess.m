% convenient dimensions of a figure for single results plot with 4 levels
% through thickness: set(gcf,'pos',[20 250 860 165])
% currently it displays rows of indices idxs = [1,10,50,180]
function figname = plotMC2D_ori_through_thicnkess(results, ind,fontsize, ylims_top, varargin)

    resentry = results(ind);

    D = [];
    nucl_ev_y = [];
    for num_of_rep=1:resentry.reps
        d = resentry.deposit{num_of_rep};
        D = [ D , d];
        nucl_ev_y = [ nucl_ev_y ; resentry.nucl_summary.events_txy{num_of_rep}(:,2:3) ];
    end

    [cmap , simtype_short] = plotMC2D_get_colormap_and_annotation(resentry);

%     [fcn_dir,~, ~] = fileparts(which(mfilename)) ;
%     addpath([fcn_dir '\..\..\results_processing'],[fcn_dir '\..\shape_factor'])
    
    sf_res = load_results(which(resentry.sf_resultsfile),'results_SF');
    oris = sf_res(resentry.sf_ind).tori;
    ori_spacing = oris(2) - oris(1);
    ori_range = 2*pi/resentry.n;
%     idxs = [1,10,50,180];
    idxs = [1, (size(d,1)-15)];
    
    angs = [size(resentry.deposit{1},2),length(idxs)];

            for ll = 1:length(idxs)
                edges = (0.5:1:50.5)*ori_spacing+(ll-1)*ori_range;
                reg_oris = d(idxs(ll),:)>0;
                zero_oris = d(idxs(ll),:)==0;
                ori_idxs = d(idxs(ll),:);
                angs(reg_oris,ll) = oris(ori_idxs(reg_oris))+(ll-1)*ori_range;
                angs(zero_oris,ll) = nan;
                binvals = oris+(ll-1)*ori_range;
                [N(:,ll),~] = histcounts(angs(:,ll),edges,'Normalization','probability');
                bc = bar(binvals,N(:,ll),1,'facecolor', 'flat','LineStyle','none');
                bc.CData = cmap;
                hold on
                leg{ll} =['y = ' num2str(idxs(ll))];
            end
            ylims = get(gca,'YLim')';
%             ytxt = 0.9*ylims(2);
            ytxt = ylims_top - 0.04;
            for ll = 1:length(idxs)
                plot([1,1]*ori_range*ll,[0,1],'k-')
                text(0.04*2*pi/resentry.n+(ll-1)*pi/2,ytxt,leg{ll},...
                    'FontSize',fontsize, 'BackgroundColor', 'white')
    %             text(xfirsttxt+(ll-1)*xspacingtxt,ytxt,leg{ll})
            end
    %         XTICKS = 0:(oris(end)/4):oris(end)*4;
            XTICKS = 0:(ori_range/4):ori_range*4;
            xticks(XTICKS)
            xticklabs = {'0', '0.25', '0.5', '0.75','1&0' ,... 
                                      '0.25', '0.5', '0.75','1' };
            xticklabels(xticklabs)
               
%         ylim(ylims)
        ylim([0,ylims_top])
        yticks(0:0.05:0.4)
        hold off
%         legend(leg,'Location','northwest')
%         if isempty(varargin)
%             title([simtype_short ', \delta=' num2str(resentry.soaIE) ', \beta='  num2str(resentry.exp_homog) events_str])
%         end
        xlabel('orientation (2\pi/n)')
        ylabel('probability')
        grid on
        set(gca,'FontSize',fontsize)
        
        
        figname = ['ori_distrib_thru_thick_' simtype_short '_soa' num2str(resentry.soaIE,'%1.0e') '_exphom' num2str(resentry.exp_homog) '_ind' num2str(ind)];
        
        %% alternative plot - 3D bar chart
%         gap = 5;
%         numcols = (1+gap)*size(N,2);
%         NN = zeros(size(N,1),numcols);
% %         ind_res = (1:size(N,2) [1, gap*(1:size(N,2)-1)+1 ];
%         ind_res = 1:(gap+1):numcols;  
%         for pp = 1:size(N,2)
%             NN(:,ind_res(pp))=N(:,pp);
%             NN(:,ind_res(pp))=N(:,pp);
%         end
%         bc3 = bar3(oris*180/pi,NN,1);
%         colormap(cmap)
%         for ff = 1:size(NN,2)
%             if any(ff==ind_res)
%                 bc3(ff).CDataMapping = 'direct';
%                 bc3(ff).CData = kron((1:length(cmap))',ones(6,4));
%             else
%                 bc3(ff).FaceAlpha = 0;
%                 bc3(ff).EdgeAlpha = 0;
%             end
%         end
        
    end % func