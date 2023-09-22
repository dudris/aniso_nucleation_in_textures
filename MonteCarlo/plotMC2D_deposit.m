% COPYRIGHT NOTICE 
 % This file is part of a dataset <Minar, Martin (2023), “Influence of surface energy anisotropy on nucleation and crystallographic texture of polycrystalline deposits”, Mendeley Data, V1, doi: 10.17632/bsdff8shbz.1>, coupled to publication of the same name by Minar, Moelans submitted to Computational Materials Science in September 2023. 
 % Distributed under GPLv3 license.
%% filename = plotMC2D_deposit(results, ind,ax, disp_nucl_sites)
function filename = plotMC2D_deposit(results, ind,ax, disp_nucl_sites)

resentry = results(ind);

num_of_rep = 2;
d = resentry.deposit{num_of_rep};
ev_txy = resentry.nucl_summary.events_txy{num_of_rep};

% __ plootting to current figure
    imagesc(d,'AlphaData',d>0)
    set(gca,'YDir','normal','DataAspectRatio',[1,1,1])
    
    [cmap , simtype_short] = plotMC2D_get_colormap_and_annotation(resentry);
    
    colormap(ax, cmap)
    cbr = colorbar;
    cbr.Ticks = [1,10:10:50];
%     cbr.
    set(gca,'FontSize',12,'TickDir','out')
%     title([simtype_short ', \delta=' num2str(resentry.soaIE) ', \beta='  num2str(resentry.exp_homog) ', ' resentry.init_code])
    if disp_nucl_sites
        hold on
        plot(ev_txy(:,2),ev_txy(:,3),'.k')
%         legend('nucl. sites')
        hold off
        nucl_sites_string = 'nucl-sites-on';
    else
        nucl_sites_string = 'nucl-sites-off';
    end
        
%     xlim([0,300])
%     yticks([10,50,180])

filename = ['deposit_ ' simtype_short '_soa' num2str(resentry.soaIE,'%1.0e') '_exphom' num2str(resentry.exp_homog) '_ind' num2str(ind) '_repnum' num2str(num_of_rep) '_' nucl_sites_string];

end % func