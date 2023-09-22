% COPYRIGHT NOTICE 
 % This file is part of a dataset <Minar, Martin (2023), “Influence of surface energy anisotropy on nucleation and crystallographic texture of polycrystalline deposits”, Mendeley Data, V1, doi: 10.17632/bsdff8shbz.1>, coupled to publication of the same name by Minar, Moelans submitted to Computational Materials Science in September 2023. 
 % Distributed under GPLv3 license.
%%
function [Aij,Sij,contangij,hij,sol_type_ij] = get_shape_factor(toprot,cycoord,n,Omg,plotting)
    
    if isnan(cycoord)
       Aij= 0;
       Sij = 0;
       hij= nan;
       contangij = nan(1,2);
       sol_type_ij = -1; % no solution, submerged
       return
    end
    
    spec = gen_input_stable_sol(toprot,cycoord,n,Omg);
    
    [Aij, Sij, hij, contangij, sol_type_ij,th_wulff_cell,ind_min_area] = find_stable_solution(spec);
    
    if plotting 
        % allowed angles
        th = linspace(-pi,pi,1001)';
        th(1)=[];
        wx = spec.wx;
        wy_shftd = spec.wy_shftd;
        params = spec.params;
        tth = th((1-params.Omega*cos(params.nfold*(th-toprot))) >=0);
        
        figure(56)
            plot(wx(tth),wy_shftd(tth),'.'),
            hold on, 
            linespec = {'rd', 'gx'};
            for hh = 1:length(th_wulff_cell)
                ls = linespec{1};
                if hh==ind_min_area
                    ls = linespec{2};
                end
                plot(wx(th_wulff_cell{hh}),wy_shftd(th_wulff_cell{hh}),ls),
            end
%             plot(wx(th_wulff),wy_shftd(th_wulff),'go'),
            plot(wx(contangij),wy_shftd(contangij),'xk'),
            quiver(wx(contangij),wy_shftd(contangij),0.1*cos(contangij),0.1*sin(contangij),0.3,'Color','r')
            plot([-3,3],[0,0],'k-','linewidth',1.2)
            hold off
            axis equal
            grid on  
            pause(.1)
    end% plotting

end % func
