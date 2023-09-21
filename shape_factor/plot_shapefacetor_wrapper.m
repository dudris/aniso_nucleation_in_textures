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
    
    
    
