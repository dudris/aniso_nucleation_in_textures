% COPYRIGHT NOTICE 
 % This file is part of a dataset <Minar, Martin (2023), “Influence of surface energy anisotropy on nucleation and crystallographic texture of polycrystalline deposits”, Mendeley Data, V1, doi: 10.17632/bsdff8shbz.1>, coupled to publication of the same name by Minar, Moelans submitted to Computational Materials Science in September 2023. 
 % Distributed under GPLv3 license.
%% spec = gen_input_stable_sol(toprot,cycoord,n,Omg)
% - generate standard input to the Winterbottom construction solver in 'find_stable_solution.m'
% - INPUT
%   - toprot ... rotation of the top grain (nucleus) in radians
%   - cycoord ... vertical shift of the shape (the 'wetting parameter' \Gamma from the paper)
%   - n ... order of symmetry
%   - Omg ... normalized strength of anisotropy \Omega = \delta*(n^2-1)
% - OUTPUT
%   - spec ... structure with the input
function spec = gen_input_stable_sol(toprot,cycoord,n,Omg)

    soaIE = Omg/(n^2-1);
    f = @(th) 1+soaIE*cos(n*(th-toprot));
    df = @(th) -n*soaIE*sin(n*(th-toprot));
    
    params.Omega = Omg;
    params.nfold=  n;
    params.soaIE = soaIE;
    spec.params = params;
    spec.toprot = toprot;

    spec.ba_forb = border_angles_v2(params,'forbb',false);
    spec.ba_forb = rotate_to_first_inetrval(spec.ba_forb,toprot);
    [spec.th_allwd, spec.ba_allwd] = GetAllowedAngles(spec.ba_forb,360,false);
    spec.ba_allwd = rotate_to_first_inetrval(spec.ba_allwd,0);
    
    % __ first probe the equilibrium isolated Wulff
    spec.ba_earscross = border_angles_v2(params,'nconvex',false);
    spec.ba_earscross = rotate_to_first_inetrval(spec.ba_earscross,toprot);
    [spec.th_allwd_earscross{1}, spec.ba_allwd_earscross] = GetAllowedAngles(spec.ba_earscross,360,false);
    % __ higher order(s) of crossed ears
    ba_earscross2 = border_angles_v2(params,'earscross',false);
    if all(all(all(~isnan(ba_earscross2))))
        for kk = 1:size(ba_earscross2,3)
            spec.ba_earscross(:,:,1+kk) = rotate_to_first_inetrval(ba_earscross2(:,:,kk),toprot);
            [spec.th_allwd_earscross{kk+1} , spec.ba_allwd_earscross(:,:,kk+1)] = GetAllowedAngles(spec.ba_earscross(:,:,1+kk),360,false);
        end
    end
    % __ for backwards compatibility
    spec.ba_ncnvx = spec.ba_earscross(:,:,1);
    spec.ba_allwd_ncnvx = spec.ba_allwd_earscross(:,:,1);
    spec.th_allwd_ncnvx = spec.th_allwd_earscross{1};
    th_allwd_ncnvx = spec.th_allwd_ncnvx;
    
    wx = @(ang) f(ang).*cos(ang) - df(ang).*sin(ang);
    wy = @(ang) f(ang).*sin(ang) + df(ang).*cos(ang);
    
    % __ y coordinates of shifted Wulff, look for zeros
    wy_shftd = @(ang) (wy(ang)-cycoord);
    spec.wy_shftd = wy_shftd;
    spec.wx = wx;
    spec.cycoord = cycoord;
    
    %__ min max value on isolated Wulff
    spec.wy_isol_minmax = [min(wy_shftd(spec.th_allwd_ncnvx(:))) , max(wy_shftd(spec.th_allwd_ncnvx(:)))];
    % __ xy positions of the corners
    spec.w_isol_corners = [wx(spec.ba_ncnvx(:,1)),wy_shftd(spec.ba_ncnvx(:,1))];
    for kk = 1:size(spec.ba_earscross,3)
        spec.w_xy_cross(:,:,kk) = [wx(spec.ba_earscross(:,1,kk)),wy_shftd(spec.ba_earscross(:,1,kk))];
    end
%     plot(wx(th_allwd_ncnvx),wy_shftd(th_allwd_ncnvx),'.'), axis equal, hold on,plot(w_isol_corners(:,1),w_isol_corners(:,2),'xk')

    % __ Wulff type... categorical, either of: {'weak','strong','strong earcross'}
    spec.wulff_type = get_wulff_type(n,Omg);
    
    % __ find roots on spline through wy_shifted
    th = linspace(-pi-0.1,pi+0.1,1001)'; % larger than (-pi,pi), technicality
    splwy = spline(th,wy_shftd(th));
    rootss_all = fnzeros(splwy); % includes forbidden
    rootss_all = rotate_to_first_inetrval(rootss_all,0);
%     plot(th,wy_shftd(th),'.',th,ppval(th,splwy),'-'), grid on
    rootss_all = rootss_all(1,:); % wy is smooth, polynomial-like function, no touch-zeros
    spec.rootss_all = rootss_all;
    
    %__ leave out forbidden angles
    rootss_allwd = is_allowed_ang(rootss_all,n,Omg,toprot,'forbb');
    % __ roots on the isolated Wulff shape
    rootss_isol = is_allowed_ang(spec.rootss_all,n,Omg,toprot,'nconvex');
    spec.rootss_allwd = remove_doubled(rootss_allwd);
    spec.rootss_isol = remove_doubled(rootss_isol);
    
    % __ get area of isolated Wulff
    spec.A_hom = polyarea(spec.wx(spec.th_allwd_ncnvx),spec.wy_shftd(spec.th_allwd_ncnvx));
    
end