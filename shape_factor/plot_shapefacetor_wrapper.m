% COPYRIGHT NOTICE 
 % This file is part of a dataset <Minar, Martin (2023), “Influence of surface energy anisotropy on nucleation and crystallographic texture of polycrystalline deposits”, Mendeley Data, V1, doi: 10.17632/bsdff8shbz.1>, coupled to publication of the same name by Minar, Moelans submitted to Computational Materials Science in September 2023. 
 % Distributed under GPLv3 license.
function plot_shapefacetor_wrapper(resultsentry)
    actions_continue = true;
    initpoint = ceil(size(resultsentry.S)*0.3)+[-1,1];
    oristep = resultsentry.tori(2)-resultsentry.tori(1);
    plot_shapefacetor(initpoint,resultsentry)
    while actions_continue
        [x, y] = ginput(1);
%         vals = round([x,y],0);
        find_ind = @(val) find(abs(resultsentry.tori-val*pi/180 )<0.5*oristep);
        inds = [find_ind(y), find_ind(x)];
        try
            plot_shapefacetor(inds,resultsentry)
        catch errM
            if strcmp(errM.identifier,'find_stable_solution:results_inconsistency')
                actions_continue = false;
                errordlg([errM.message '\n Plotted saved contact angles. Further plotting aborted.'],'error_inconsistency_dialog')
            else
                helpdlg('Not valid point. Click again.')
            end
        end
    end
end
    
    
    
