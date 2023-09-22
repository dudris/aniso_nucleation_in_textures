% COPYRIGHT NOTICE 
 % This file is part of a dataset <Minar, Martin (2023), “Influence of surface energy anisotropy on nucleation and crystallographic texture of polycrystalline deposits”, Mendeley Data, V1, doi: 10.17632/bsdff8shbz.1>, coupled to publication of the same name by Minar, Moelans submitted to Computational Materials Science in September 2023. 
 % Distributed under GPLv3 license.
% plot_shapefacetor(top_bot_ori,resultsentry)

function plot_shapefacetor(top_bot_ori,resultsentry)
    res = resultsentry;
    
    ind_toprot = top_bot_ori(1);
    ind_botrot = top_bot_ori(2);

    toprot = res.tori(ind_toprot);
    spec.toprot = toprot;
    cycoord = res.Cshift(ind_toprot, ind_botrot);
    spec = gen_input_stable_sol(toprot,cycoord,res.n,res.Omega);
    
    % __ get area of isolated Wulff
    global A_hom
    A_hom = res.A_hom;
    
    [Aij, Sij, hij, contangij, sol_type_ij,th_wulff_cell,ind_min_area] = find_stable_solution(spec);
    
    figure(33)
    ba_forbb = border_angles_v2(spec.params,'forbb',false);
    %     ba_ncnvx = border_angles_v2(params,'forbb',false);
    ba_forbb = rotate_to_first_inetrval(ba_forbb,toprot);
    [th_allwd, ~] = GetAllowedAngles(ba_forbb,360,false);
    
    subplot(121)
        t=  res.tori*180/pi;
        b=  res.bori*180/pi;
        S = res.S;
        % __ for more comprehensive plots don't show the diagonal NaNs
        alphadata = (~isnan(S) | logical(eye(sqrt(numel(res.stab_cond)))) ) & res.stab_cond;
        alphadata =  res.stab_cond;
        alphadata = alphadata - diag(diag(alphadata)) + diag(diag(res.stab_cond)); % stability condition to have preference in plotting
%         imagesc(b,t,S)
        imagesc(b,t,S,'AlphaData',alphadata)
        set(gca,'DataAspectRatio',[1,1,1],'YDir','normal')
        colorbar
        hold on
        plot(b(ind_botrot),t(ind_toprot),'o','MarkerFaceColor','r')
%         plot(t(ind_toprot),b(ind_botrot),'o','MarkerFaceColor','r')
        hold off
        ylabel('top grain \alpha_2 (deg)')
        xlabel('bottom grain \alpha_1 (deg)')
        title(['n =' num2str(res.n) ', \delta =' num2str(res.soaIE) ', (\alpha_1,\alpha_2) = (' num2str([res.bori(ind_botrot),res.tori(ind_toprot)]*180/pi,'%3.1f ') ')'])
    
    subplot(122)
        lims = [-3,3];
        area(lims,[0,0],-1.5,'FaceColor',[.9,.9,.9])
        hold on
        plot(spec.wx(th_allwd),spec.wy_shftd(th_allwd),'b.','markersize',4)
        if ~isempty(th_wulff_cell)
            for hh = 1:length(th_wulff_cell)
                if hh~=ind_min_area
                    th_sol = th_wulff_cell{hh};
                    plot(spec.wx(th_sol),spec.wy_shftd(th_sol),'o','Color',0.6*[1,1,1])
                end% if
            end% for
            th_sol = th_wulff_cell{ind_min_area};
            plot(spec.wx(th_sol),spec.wy_shftd(th_sol),'go')
        end
        
%         hold on
%         plot(0,0,'``bx','LineWidth',1.5)
        TH([1,2]) = res.contangs(ind_toprot,ind_botrot,:);
        res_contangs_col = reshape(res.contangs(ind_toprot,ind_botrot,:),[2,1]);
        check_contangs = abs(res_contangs_col-contangij')./contangij' < 1e-3;
        
%         TH([1,2]) = contangij;
        if ~all(abs(TH)<1e-10) % NEGATE: ( both are around zero => no intersection )
            plot(spec.wx(TH),spec.wy_shftd(TH),'xr','linewidth',1)
            quiver(spec.wx(TH),spec.wy_shftd(TH),cos(TH),sin(TH),0,'Color','r')
        end
%         plot([-1.1,1.1],res.Cshift(ind_botrot,ind_toprot)*[1,1],'k-','linewidth',1.2)
%         plot(lims,(res.h(ind_toprot,ind_botrot)+res.Cshift(ind_toprot,ind_botrot))*[1,1],'k--')
        hold off
        axis  equal
        grid on
        ylim(lims)
        xlim(lims)
        
        if ~all(check_contangs)
            if all(~isnan([res_contangs_col;contangij']))
                msg = ['Contact angles re-computed for visualization are not consistent with what is in the results. Computed: ' num2str(res_contangs_col') ', saved: ' num2str(contangij) '.'];
                ME = MException('find_stable_solution:results_inconsistency',msg);
                throw(ME)
            end
        end
end % func