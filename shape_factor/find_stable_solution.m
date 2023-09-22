% COPYRIGHT NOTICE 
 % This file is part of a dataset <Minar, Martin (2023), “Influence of surface energy anisotropy on nucleation and crystallographic texture of polycrystalline deposits”, Mendeley Data, V1, doi: 10.17632/bsdff8shbz.1>, coupled to publication of the same name by Minar, Moelans submitted to Computational Materials Science in September 2023. 
 % Distributed under GPLv3 license.
%% [Aij, Sij, hij, contangij, sol_type_ij,th_wulff_cell,ind_min_area] = find_stable_solution(spec)
% - function to find the minimal-area solution to Winterbottom construction
% with arbitrary orientation of the nucleus (top grain)
% -INPUT
%   - generated using: spec = gen_input_stable_sol(toprot,cycoord,n,Omg)
% - OUTPUT (consider unit-radius Wulff shape)
%   - Aij ... area of the minimal-area solution (scalar)
%   - Sij ... shape factor of the minimal-area solution (scalar)
%   - hij ... largest distance between the truncating line and the points 
%   on the minimal-area solution (tallest point) (scalar)
%   - contangij ... a pair of contact angles in the minimal-area solution
%   (vector 1x2)
%   - sol_type_ij ... integer indicating the type of solution in the minimal-area solution
%       - -2: no solution, emerged
%       - -1: no solution, submerged
%       -   1: truncated isolated Wulff shape (regular Winterbottom construction) 
%       -  2: inverted shape solution
%       -  3: single-branch solution
%       -  4: crossed-ears solution
%       -  5: emerged shape solution
%   - th_wulff_cell ... normal angles along all the considered solutions (cell of vectors, a vector per solution)
%   - ind_min_area ... index of the minimal-area solution in th_wulff_cell

function [Aij, Sij, hij, contangij, sol_type_ij,th_wulff_cell,ind_min_area] = find_stable_solution(spec)

%     global A_hom
    wx = spec.wx;
    wy_shftd = spec.wy_shftd;
    
    tolerance_y = 1e-4;
    
    % __no solution 
    if numel(spec.rootss_allwd)==0 || numel(spec.rootss_allwd)==1 || ... 
       (strcmp(spec.wulff_type,'strong') && any(abs(spec.wy_isol_minmax))<=tolerance_y )
        
        ind_min_area = nan;
        if spec.cycoord < 0 % complete dewetting    
           Aij= 1;
           Sij = 1;
           hij= nan;
           contangij = nan(1,2);
           sol_type_ij = -2; % no solution, emerged
           th_wulff_cell = {};
        else % complete wetting
           Aij= 0;
           Sij = 0;
           hij= nan;
           contangij = nan(1,2);
           sol_type_ij = -1; % no solution, submerged
           th_wulff_cell = {};
        end % if shift up/down
        return
    
    elseif spec.params.nfold==3 && spec.wy_isol_minmax(2)<=1e-10 % 3fold has no inverted eq.shapes or crossed ears
        Aij= nan;
        Sij = nan;
        hij= nan;
        contangij = nan(1,2);
        sol_type_ij = -1; % no solution, submerged
        th_wulff_cell = {};
        ind_min_area = nan;
        return
        
%     elseif all(spec.wy_isol_minmax>0) % for strong and strong_earcross solution only larger than isolatd Wulff => disregarded for now
%         Aij= 1;
%         Sij = 1;
%         hij= nan;
%         contangij = nan(1,2);
%         sol_type_ij = -2; % no solution, emerged
%         th_wulff_cell = {};
%         ind_min_area = nan;
%         return
        
    end % if no solution 
    
    
    % __ select the class of solutions
    % __ when isolated Wulff truncated - 2 roots assumed
    if spec.wy_isol_minmax(1)<0 && spec.wy_isol_minmax(2)>0 
%             assert(numel(rootss_isol)==2,['number of roots on isolated Wulff = ' num2str(numel(rootss_isol))])
        
        th_allwd_ncnvx = spec.th_allwd_ncnvx;
        % __ add roots to the contour
        th_allwd_ncnvx = unique([th_allwd_ncnvx; spec.rootss_isol(:)]);
        % __ indices of theta where Wulff above the line
        gz2 = wy_shftd(th_allwd_ncnvx) > -1e-10;
        % __ theta angless of truncated Wulff
        th_wulff = th_allwd_ncnvx(gz2);
        % for plotting compatibility
        th_wulff_cell{1}=th_wulff;
        ind_min_area = 1;
        
        if (numel(spec.rootss_isol))==2 % the ordinary case
            contangij = spec.rootss_isol;
            % check that Left and Right contacts are there, negative goes first
            contangij = check_contangs(contangij);
            
        elseif (numel(spec.rootss_isol))~=2 % when corners are near 0, numerics gets tricky; numel can be 0, 1 or 3
            contangij = nan(1,2); % this case to be identified in results by having sol_type==1 and contangs==nan
        end
        Aij = polyarea(wx(th_wulff),wy_shftd(th_wulff));
        hij= spec.wy_isol_minmax(2);
        Sij = Aij/spec.A_hom;
        sol_type_ij = 1; % truncated isolated Wulff solution
    
    elseif spec.wy_isol_minmax(2)<=1e-10 % wetting, isolated shape submerged
        
        % __ only 4fold or 6fold; 3fold was ruled out above
        if strcmp(spec.wulff_type,'strong') % only inverted shape or single_branch solution possible
            
            [th_wulff_cell{1}, contangij_all(1,:)]= find_single_br_w(spec);
            
            [th_wulff_cell, contangij_all]= find_inverted_wulff(spec, th_wulff_cell, contangij_all);
            
            % __ compute areas of the possible solutions
            Aij_all=nan(length(th_wulff_cell),1);
            for kk=1:length(Aij_all)
                Aij_all(kk) = polyarea(wx(th_wulff_cell{kk}),wy_shftd(th_wulff_cell{kk}));
            end
                
            [Aij, ind_min_area ] = min(Aij_all);
            Sij = Aij/spec.A_hom;
            th_wulff = th_wulff_cell{ind_min_area};
            contangij = contangij_all(ind_min_area,:);
            contangij = check_contangs(contangij);
            
            % __ result depending of solution
            if any(size(th_wulff)==0) % can happen with 2-fold - no solution
               Aij= 0;
               Sij = 0;
               hij= nan;
               contangij = nan(1,2);
               sol_type_ij = -1; % no solution, submerged
               th_wulff_cell = {};
            elseif  ind_min_area==1 % single branch wulff
                hij= max(wy_shftd(th_wulff));
                sol_type_ij = 3; % single branch solution 
                
            elseif ind_min_area>1 % inverted wulff
                hij= -spec.wy_isol_minmax(2);
                sol_type_ij = 2; % inverted shape solution 
                
            else
                error(['unexpected solution in wulff_type=' wulff_type ])
            end
            
        
        elseif strcmp(spec.wulff_type,'strong earcross') % multiple shapes possible
            % __ single branch wulff
            [th_wulff_cell{1}, contangij_all(1,:)]= find_single_br_w(spec);
            % __ inverted wulff
            [th_wulff_cell, contangij_all]= find_inverted_wulff(spec, th_wulff_cell, contangij_all);
%             [th_wulff_cell{2}, contangij_all{2}]= find_inverted_wulff(spec);

            % __ cross-ear wulff
%             [th_wulff_cell, contangij_cell]= find_crossear_wulff(params,toprot,wx,wy_shftd,rootss_allwd,th_wulff_cell,contangij_cell);
            [th_wulff_cell, contangij_all] = find_crossear_wulff(spec, th_wulff_cell, contangij_all);
            
            Aij_all=nan(length(th_wulff_cell),1);
            Aij_all = cellfun(@(x) polyarea(wx(x),wy_shftd(x)), th_wulff_cell');
            
            [Aij, ind_min_area ] = min(Aij_all);
            Sij = Aij/spec.A_hom;
            th_wulff = th_wulff_cell{ind_min_area};
            contangij = contangij_all(ind_min_area,:);
%             contangij = check_contangs(contangij);
            
            max_cross_order = size(spec.ba_allwd_earscross,3);
            if ind_min_area==1
                hij= max(wy_shftd(th_wulff));
                sol_type_ij = 3; % single branch wulff
            elseif ind_min_area>1 && ind_min_area<=(max_cross_order+1)
                hij= -spec.wy_isol_minmax(2);
                sol_type_ij = 2; % inverted shape solution 
            elseif ind_min_area>(max_cross_order+1) 
                hij= max(wy_shftd(th_wulff));
                sol_type_ij = 4; % crossed ears solution 
            end 
            
%         elseif strcmp(spec.wulff_type,'general') 
            
        end
        
    elseif spec.wy_isol_minmax(1)>0% dewetting, isolated shape emerged
        [th_wulff_cell, contangij_all]= find_emerged_wulff(spec);
        Aij_all=nan(length(th_wulff_cell),1);
        Aij_all = cellfun(@(x) polyarea(wx(x),wy_shftd(x)), th_wulff_cell');

        [Aij, ind_min_area ] = min(Aij_all);
        Sij = Aij/spec.A_hom;
        th_wulff = th_wulff_cell{ind_min_area};
        contangij = contangij_all(ind_min_area,:);
        hij= max(wy_shftd(th_wulff));
        sol_type_ij = 5; % emerged Wulff solution
    end% class of solutions
    
end % func


%%
% function [th_invert_wulff,contangij] = find_inverted_wulff(ba_allwd_ncnvx,params,toprot,wy_isol_minmax,w_isol_corners,ba_ncnvx,wx,wy_shftd,rootss)
function [th_wulff_cell, contangij_all] = find_inverted_wulff(spec, th_wulff_cell, contangij_all)

params = spec.params;
wy_shftd = spec.wy_shftd ;

max_cross_order = size(spec.ba_allwd_earscross,3);
th_invert_wulff_cell = cell(1,max_cross_order);
contangij = nan(max_cross_order,2);

is_branch_roots = find_branch(spec.rootss_allwd,params,spec.ba_allwd);
    
% __ there may be 2 roots on 1 branch, one not present in inverted
% shape and thus to be omitted
br_double_cross = sum(is_branch_roots,2)==2;
if any(br_double_cross)  % do this only once
    % __ corner x coord
%     x_corner = spec.wx(angs_topcorn(1));
    % __ the two roots on 1 branch 
    ind_dcr = find(is_branch_roots(br_double_cross,:));
    double_cross_roots = spec.rootss_allwd(ind_dcr);
    % __ absolute of the x distance
    x_roots_dist = abs(spec.wx(double_cross_roots));
    % __ index of the farther
    [~,indmaxdist] = max(x_roots_dist);
    % __ remove the farther from the list
    is_branch_roots(br_double_cross,ind_dcr(indmaxdist))= false;
end % if branch double cross

for cross_order = 1:max_cross_order
    ba_allwd_subtract = spec.ba_allwd_earscross(:,:,cross_order);
    ba_mix_intrm = [spec.ba_forb;ba_allwd_subtract];
    [th_allwd_ears, ~] = GetAllowedAngles(ba_mix_intrm ,360,false);

    % __ take the topmost corner/cross of the corresponding order SMALLER
    % THAN ZERO
    cond_cross_lz = spec.w_xy_cross(:,2,cross_order) < -0.001;
    ind_topcorn = spec.w_xy_cross(:,2,cross_order) == max(spec.w_xy_cross(cond_cross_lz,2,cross_order)) ;
    if sum(ind_topcorn)==2 % unambiguous determination of the corner 
        ind_topcorn = find(ind_topcorn,1); % take the first
    end
    % __ the two normal angles at the corner = border angles
    angs_topcorn = spec.ba_earscross(ind_topcorn,:,cross_order);
    cy_topcorn = wy_shftd(angs_topcorn(1));
    %     angs_topcorn = spec.ba_ncnvx(ind_topcorn,:);

    from_corner_to_zero = (wy_shftd(th_allwd_ears) >= cy_topcorn )  & wy_shftd(th_allwd_ears)<=0;
%     from_corner_to_zero = (wy_shftd(th_allwd_ears) >= spec.wy_isol_minmax(2) )  & wy_shftd(th_allwd_ears)<=0;
    
    is_branch_corner = find_branch(angs_topcorn,params,spec.ba_allwd);
    
    br_is_corn_and_has_roots = all([any(is_branch_corner,2) , any(is_branch_roots,2) ],2);
    % __ there are 2 branches with roots
    two_br_roots_exist = sum(br_is_corn_and_has_roots)==2; 
    
    if two_br_roots_exist
        ind_br = find(br_is_corn_and_has_roots);
        ind_br_root = find(any(is_branch_roots(ind_br,:))); % 
        br_roots = spec.rootss_allwd(ind_br_root);
        br_roots_tang = rotate_to_first_inetrval(br_roots-pi/2,0);
        % __ the 2 roots correspond to L and R contat point (+ and -)
        are_left_right = min(br_roots_tang)<0 && max(br_roots_tang)>0;
        % __ right contact is on the left and vice versa as should be in inverted 
        are_inverted = spec.wx(min(br_roots)) < spec.wx(max(br_roots));
        is_wulff_inv = two_br_roots_exist & are_left_right & are_inverted;
    else
        is_wulff_inv = false;
    end
    
    if is_wulff_inv
        contangij(cross_order,:) = br_roots;
        
        th_ear_tozero = th_allwd_ears(from_corner_to_zero);
        
        % __ through the crossing branches 
        is_branch_all = find_branch(th_ear_tozero, params,spec.ba_allwd);
        is_on_either_of_two = is_branch_all(ind_br(1),:) | is_branch_all(ind_br(2),:);
        is_th_invert_wulff = is_on_either_of_two;
%         is_th_invert_wulff = false(size(th_ear_tozero));
%         for bb=1:2
% %             % __ set correct limits on the intervals (both from single branch)
%             th_limits(1) = angs_topcorn(1,is_branch_corner(ind_br(bb),:));
%             th_limits(2) = spec.rootss_allwd(ind_br_root(bb));
%             is_branch = ( th_ear_tozero >= min(th_limits)) & ( th_ear_tozero <= max(th_limits) );
%             is_th_invert_wulff = is_th_invert_wulff | is_branch;
%         end % for branch in single cross point
%             clear is_branch
        
        th_invert_wulff_cell{cross_order} = th_ear_tozero(is_th_invert_wulff);
        th_invert_wulff_cell{cross_order} = unique([th_invert_wulff_cell{cross_order} ; contangij(cross_order,:)' ; angs_topcorn(1)]);
        
%         th = linspace(-pi,pi,2001)';
%         th(1)=[];
%         tth = th((1-params.Omega*cos(params.nfold*(th-spec.toprot))) >=0);
%         figure(56)
%             plot(spec.wx(tth),spec.wy_shftd(tth),'.'),
%             hold on, 
%             plot(wx(th_invert_wulff),wy_shftd(th_invert_wulff),'o'),
%             plot(spec.wx(spec.rootss_allwd),spec.wy_shftd(spec.rootss_allwd),'xk'),
% %             plot(wx(th_invert_wulff),wy_shftd(th_invert_wulff),'k'),
%             plot(spec.wx(angs_topcorn(1)),spec.wy_shftd(angs_topcorn(1)),'sg')
%             hold off
%             axis equal
%             grid on

    else
        contangij(cross_order,:) = nan(1,2);
        th_invert_wulff_cell{cross_order} = nan;
    end % two branches at the topmost corner have roots
end

% __ append the found solutions to the others
th_wulff_cell = [th_wulff_cell , th_invert_wulff_cell];
contangij_all = [contangij_all ; contangij];

end% func

%% function
% function [th_wulff_cell, contangij_cell]= find_crossear_wulff(params,toprot,wx,wy_shftd,rootss,th_wulff_cell,contangij_cell)
function [th_wulff_cell, contangij_all]= find_crossear_wulff(spec, th_wulff_cell, contangij_all)
    
    params = spec.params;
    wx = spec.wx;
    wy_shftd = spec.wy_shftd;
    rootss = spec.rootss_allwd;
    
    th_wulff_cell_length = length(th_wulff_cell);
    
    th_allwd_gz = spec.th_allwd(wy_shftd(spec.th_allwd) >=0);
    
    % __ assign roots to branches
    % __ size(is_branch_root)=[n,numel(rootss)]
    is_branch_root =  find_branch(rootss,params,spec.ba_allwd);
   
    max_cross_order = size(spec.ba_allwd_earscross,3);
    
    % __ cycle through crossed ears orders
    % __ cross order = 1 are the isolated shape corners - irreleveant here
    for cross_order = 2:max_cross_order    
        % __ location and angles in cross-ear points
        % __ in th_crossear  1 row... theta angles in 1 cross point, ie. 1 row is NOT 1 branch
        % as in border angles
        th_crossear = spec.ba_earscross(:,:,cross_order);
        xy_crser = spec.w_xy_cross(:,:,cross_order);
        
        % __ consider only crossed ears above 0 and not too close to 0
        cond_crossd_ears_gz = xy_crser(:,2)>0.001;
        if any(cond_crossd_ears_gz)
            % __ consider only the one closest to zero - that implies smaller area
            closest_to_zero = min(xy_crser(cond_crossd_ears_gz,2));
            if any(size(closest_to_zero)>1) % in case two corners at the same height
                closest_to_zero = closest_to_zero(1);
            end
            cond_min_crossd_ears_gz =  xy_crser(:,2) ==  closest_to_zero;
            ind_candidates = find(cond_min_crossd_ears_gz);
        else
            ind_candidates = [];
        end % if any above 0

        % __ check 2 roots on 1 branch
        br_double_cross = sum(is_branch_root,2)==2;
        % __ then BOTH these roots are not at base of cross-ear solution and
        % one or both can be removed (single-branch treated elsewhere)
        if any(br_double_cross) && cross_order==2 % do this only once
            % __ to exclude the root of single branch which is farther in
            % x-coordinate on Wulff
            % __ find the two roots 
            ind_roots_dc_branch = abs(is_branch_root(br_double_cross,:)-1)<1e-16;
            % __ which of all roots is the one which in single-branch solution
            % and has the farther x-coordinate on Wulff
            ind_farther = abs(abs(wx(rootss))-max(abs(wx(rootss(ind_roots_dc_branch)))))<1e-10;
            % __ get rid of it
            rootss(ind_farther) = [];
            % __ recompute branch-root assignment
            is_branch_root =  find_branch(rootss,params,spec.ba_allwd);
        end
    
        % __ cycle left from when all crossed ears solutions probed
        % __ now assume that length(ind_candidates)=1
        for kk=1:length(ind_candidates)
            % __ index in CROSS POINTS (not branches)
            indgz = ind_candidates(kk); % greater than zero
            is_branch_crser{kk} =  find_branch(th_crossear(indgz,:),params,spec.ba_allwd);

            % __ true for branch containing both cross-ear point and a root
            br_is_crser_and_has_roots(:,kk) = all([any(is_branch_crser{kk},2) , any(is_branch_root,2) ],2);
            % __ flag is true when there are 2 branches crossing ears, each having a root
            is_wulff_crossear(kk,1) = sum(br_is_crser_and_has_roots(:,kk))==2; 

            if is_wulff_crossear(kk,1)
                % __ indices of branches meeting in this cross point
                ind_br(:,kk) = find(br_is_crser_and_has_roots(:,kk));

                % __ indices of the two roots in the rootss vetor
                contangij_all(th_wulff_cell_length+1,:) = rootss(any(is_branch_root(ind_br(:,kk),:))); % 
                contangij_all(th_wulff_cell_length+1,:) = sort(contangij_all(th_wulff_cell_length+1,:),'ascend');

                is_th_wulff_cell{kk} = false(size(th_allwd_gz));

                % __ through the crossing branches 
                for bb=1:2
                    % __ set correct limits on the intervals (both from single branch)
                    th_limits(1) = th_crossear(indgz,is_branch_crser{kk}(ind_br(bb,kk),:));
    %                 if length(is_branch_root(ind_br(bb,kk),:))>1
    %                     pause
    %                 end
                    th_limits(2) = rootss(is_branch_root(ind_br(bb,kk),:));

                    is_nojump_in_lim = abs(diff(th_limits))<pi;

                    if is_nojump_in_lim
                        is_branch = ( th_allwd_gz >= min(th_limits)) & ( th_allwd_gz <= max(th_limits) );
                    else
                        is_branch = ( ( th_allwd_gz <= min(th_limits)) & ( th_allwd_gz >= -pi) ) | ...
                                            ( ( th_allwd_gz >= max(th_limits) ) & ( th_allwd_gz <= pi ) );
    %                     is_branch = ( th_allwd_gz_cond >= min(th_limits)) & ( th_allwd_gz_cond <= max(th_limits) );
                    end

                    is_th_wulff_cell{kk} = is_th_wulff_cell{kk} | is_branch;
                end % for branch in single cross point
                    clear is_branch

    %             % right branch ... negative angles
    %             is_branchR = ( th_allwd_gz_cond <= th_crossear(indgz,1)  & th_allwd_gz_cond >= contangij_cell{kk}(1) );
    %             is_th_wulff_cell{kk} = is_th_wulff_cell{kk} |  is_branchR;
    %             
    % %             figure(56)
    % %             plot(wx(th_allwd),wy_shftd(th_allwd),'.')
    % %             hold on
    % %             plot([-3,3],[0,0],'k-')
    % %             plot(wx(th_allwd_gz_cond(is_th_wulff_cell{kk})),wy_shftd(th_allwd_gz_cond(is_th_wulff_cell{kk})),'ro')
    %             
    %             % left branch ... positive angles
    %             is_branchL = ( th_allwd_gz_cond >= th_crossear(indgz,1)  & th_allwd_gz_cond <= contangij_cell{kk}(2) );
    %             is_th_wulff_cell{kk} = is_th_wulff_cell{kk} |  is_branchL;

                % __ possibilities of 2 other solutions tried in th_wulff_cell and contangij_cell
                th_wulff_cell{th_wulff_cell_length+1} = th_allwd_gz(is_th_wulff_cell{kk});
                % __ add roots and cross point to the result
                th_wulff_cell{th_wulff_cell_length+1} = unique([th_wulff_cell{th_wulff_cell_length+1} ; contangij_all(th_wulff_cell_length+1,:)' ; th_crossear(indgz,1)  ]);
                th_wulff_cell_length = th_wulff_cell_length+1;

    %             figure(57)
    %             plot(wx(th_allwd),wy_shftd(th_allwd),'.'),axis equal,grid on, hold on, plot(wx(rootss),wy_shftd(rootss),'xk')
    %             hold on
    %             plot([-3,3],[0,0],'jk-')
    %             plot(wx(th_wulff_cell{kk+2}),wy_shftd(th_wulff_cell{kk+2}),'ro')
    % %             plot(wx(th_crossear(ind(1,kk),1)),wy_shftd(th_crossear(ind(1,kk),1)),'ks','linewidth',1.4)
    %             plot(wx(th_crossear(indgz,1)),wy_shftd(th_crossear(indgz,1)),'ks','linewidth',1.4)
    %             plot(wx(contangij_cell{kk+2}),wy_shftd(contangij_cell{kk+2}),'kx','linewidth',1.4)
    %             quiver(wx(contangij_cell{kk+2}),wy_shftd(contangij_cell{kk+2}),cos(contangij_cell{kk+2}),sin(contangij_cell{kk+2}))
    %             quiver(wx(contangij_cell{1}),wy_shftd(contangij_cell{1}),cos(contangij_cell{1}),sin(contangij_cell{1}))
    %             hold off
    % %             plot(wx(th_wulff_cell{1}),wy_shftd(th_wulff_cell{1}),'ro')
            else
                th_wulff_cell{th_wulff_cell_length+1} = nan;
                contangij_all(th_wulff_cell_length+1,:) = nan(1,2);
                th_wulff_cell_length = th_wulff_cell_length+1;
            end % if is solution

        end%for candidate cross-ear solutions
    end % for cross order
    
%     figure(56)
%     hold off
end % func

%% function
function [th_wulff_cell, contangij_all]= find_emerged_wulff(spec)
    
    params = spec.params;
    wx = spec.wx;
    wy_shftd = spec.wy_shftd;
    rootss = spec.rootss_allwd;
    
%     th_wulff_cell_length = length(th_wulff_cell);
    
    th_allwd_gz = spec.th_allwd(wy_shftd(spec.th_allwd) >=0);
    
    % __ assign roots to branches
    % __ size(is_branch_root)=[n,numel(rootss)]
    is_branch_root =  find_branch(rootss,params,spec.ba_allwd);
    
    % __ check 2 roots on 1 branch
    br_double_cross = sum(is_branch_root,2)==2;
    % __ then BOTH these roots are not at base of cross-ear solution and
    % one or both can be removed (single-branch treated elsewhere)
    if any(br_double_cross)
        % __ to exclude the root of single branch which is farther in
        % x-coordinate on Wulff
        % __ find the two roots 
        ind_roots_dc_branch = abs(is_branch_root(br_double_cross,:)-1)<1e-16;
        % __ which of all roots is the one which in single-branch solution
        % and has the farther x-coordinate on Wulff
        ind_farther = abs(abs(wx(rootss))-max(abs(wx(rootss(ind_roots_dc_branch)))))<1e-10;
        % __ get rid of it
        rootss(ind_farther) = [];
        % __ recompute branch-root assignment
        is_branch_root =  find_branch(rootss,params,spec.ba_allwd);
    end
    
    % __ form all possible pairs from the rootss
    combs_roots = nchoosek(rootss,2);
    % __ tangent angles
    combs_roots_tang = rotate_to_first_inetrval(combs_roots-pi/2,0);
    combs_roots_tang = sort(combs_roots_tang,2,'ascend');
    are_left_right = combs_roots_tang(:,1)<0 & combs_roots_tang(:,2)>0;
    % __ right contact is on the right and left on the left
    combs_roots_sorted = rotate_to_first_inetrval(combs_roots_tang+pi/2,0);
    are_not_inverted = spec.wx(combs_roots_sorted(:,1)) > spec.wx(combs_roots_sorted(:,2));
    % __ complying pais of contact angles
    combs_roots = combs_roots(are_left_right & are_not_inverted,:);
    
    is_wulff_emrgd = ~isempty(combs_roots);
    
    if is_wulff_emrgd
    
        if size(combs_roots,1)==1
            roots_base = combs_roots;
        elseif size(combs_roots,1)>1
            % __ the pair with shortest x distance corresponds to lowerst possible
            % order of nearest crossing and thus connects the  isolated Wulff with
            % substrate plane by the smallest possible area
        %     DOES NOT ALWAYS WORK - take 2 with shortest distance and compare
        %     areas
            xcoor_roots = spec.wx(combs_roots);
            xdist_roots = abs(diff(xcoor_roots,[],2));
            [~, ind_min_xdist] = sort( xdist_roots,'ascend') ;
            % __ take num_sols with smallest xdist
            num_sols = 2;
            roots_base = combs_roots(ind_min_xdist(1:num_sols),:);
    %         roots_base = combs_roots(ind_min_xdist,:);
        end
        
        nfold = spec.params.nfold;
        flat_side_angs = unique(rotate_to_first_inetrval((1:nfold)*2*pi/nfold+pi/nfold,spec.toprot));
        is_branch_flat_side = find_branch(flat_side_angs,params,spec.ba_allwd);
        
        for rr = 1:size(roots_base,1)
            rb = roots_base(rr,:);
            is_branch_root_base = find_branch(rb,params,spec.ba_allwd);

            ind_flat_side1 = all((is_branch_flat_side-is_branch_root_base(:,1))==0);
            base1_ends = [ flat_side_angs(ind_flat_side1),rb(1)];
            % __ jump
            if abs(diff(base1_ends))>pi
                base1 = [ linspace(min(base1_ends), -pi, 50) , linspace(pi, max(base1_ends), 50) ]';
                base1 = rotate_to_first_inetrval(base1,0);
            else
                base1 = linspace(base1_ends(1), base1_ends(2), 100)';
            end
            
            ind_flat_side2 = all((is_branch_flat_side-is_branch_root_base(:,2))==0);
            base2_ends = [ flat_side_angs(ind_flat_side2),rb(2)];
            % __ jump
            if abs(diff(base2_ends))>pi
                base2 = [ linspace(min(base2_ends), -pi, 50) , linspace(pi, max(base2_ends), 50) ]';
                base2 = rotate_to_first_inetrval(base2,0);
            else
                base2 = linspace(base2_ends(1), base2_ends(2), 100)';
            end
            
%             ind_flat_side2 = all((is_branch_flat_side-is_branch_root_base(:,2))==0);
%             base2 = linspace(flat_side_angs(ind_flat_side2), rb(2), 100)';
%             base2 = rotate_to_first_inetrval(base2,0);
            
            th_wulff = [spec.th_allwd_ncnvx ; base1 ; base2];
            k = convhull(spec.wx(th_wulff), spec.wy_shftd(th_wulff));
            th_wulff = th_wulff(k);
        
            th_wulff_cell{rr} = th_wulff;
            contangij_all(rr,:) = rb;

%             figure(56)
%                 plot(wx(spec.th_allwd),wy_shftd(spec.th_allwd),'.')
%                 axis equal
%                 hold on
%                 plot([-3,3],[0,0],'k-')
%     %             plot(wx(flat_side_angs),wy_shftd(flat_side_angs),'o')
% %                 plot(wx(th_wulff),wy_shftd(th_wulff),'o')
%                 plotpart = base2;
%                 plot(wx(plotpart),wy_shftd(plotpart),'o')
%                 hold off
%                 for kk = 1:size(ind_min_xdist)
%                     rr = combs_roots(ind_min_xdist(kk),:);
%                     plot(wx(spec.th_allwd),wy_shftd(spec.th_allwd),'.')
%                     axis equal
%                     hold on
%                     plot([-3,3],[0,0],'k-')
%                     plot(wx(rr),wy_shftd(rr'),'o')
%                     hold off
%                     pause
%                 end
%                 hold off
%                 plot(wx(th_allwd_gz_cond(is_th_wulff_cell{kk})),wy_shftd(th_allwd_gz_cond(is_th_wulff_cell{kk})),'ro')
        end
        
    else
        th_wulff_cell{1} = nan;
        contangij_all(1,:) = nan(1,2);
        
    end % if is_wulff_emrgd
end% function

%% 
% __ to decide on what branch angle(s) x lay
% is_on_branch ... bool, size=[n,length(x)]
%   branches in rows, probed angle values in column. True in the respective
%   branch-row
function is_on_branch = find_branch(x,params,ba_allwd)
    ba_allwd = sort(ba_allwd,2,'ascend');
    is_on_branch = false(params.nfold,numel(x));
    
    if size(x,1)>1 && size(x,2)==1
        x = x';
    elseif size(x,1)>1 && size(x,2)>1
        error(['Function ''find_branch'' operates only on vectors but matrix was given, size(x)=[' num2str(size(x)) ']'])
    end
    nojump = diff(ba_allwd,[],2)<pi;
    
    find_branch_nojump = @(x) ba_allwd(nojump,1) <= x & ba_allwd(nojump,2) >= x;
    is_on_branch(nojump,:) = find_branch_nojump(x);

    find_branch_jump = @(x) ( ba_allwd(~nojump,1) >= x & x >= -pi ) | (  ba_allwd(~nojump,2) <= x  &  x <= pi);
    is_on_branch(~nojump,:) = find_branch_jump(x);
    
end% func

%%
% function [th_wulff, contangij]= find_single_br_w(ba_allwd_ncnvx,params,toprot,wy_isol_minmax,w_isol_corners,ba_ncnvx,wy_shftd,rootss)
function [th_wulff, contangij]= find_single_br_w(spec)
    
    wy_shftd = spec.wy_shftd ;

    % __ necessary modification to obtain correct allowed angles 
    spec.ba_allwd_ncnvx = rotate_to_first_inetrval(spec.ba_allwd_ncnvx,0);
    % __ get border angles when only forbidden excluded
    ba_forb = border_angles_v2(spec.params,'forbb',false);
    ba_forb = rotate_to_first_inetrval(ba_forb,spec.toprot);
    % __ from only-forbidden exclude non-convex allowed angles
    %  __ then only 'allowed ears' remain
    ba_comb = [ba_forb;spec.ba_allwd_ncnvx];
    [th_allwd_ears, ~] = GetAllowedAngles(ba_comb,360,false);

    % __ take the topmost corner
    ind_topcorn = spec.w_isol_corners(:,2) == max(spec.w_isol_corners(:,2)) ;
    % __ the two normal angles at the corner = border angles
    angs_topcorn = spec.ba_ncnvx(ind_topcorn,:);

    corner_to_zero = (wy_shftd(th_allwd_ears) >= spec.wy_isol_minmax(2) )  & wy_shftd(th_allwd_ears)<=0;

    [~, ba_allwd] = GetAllowedAngles(ba_forb,360,false);
    ba_allwd = rotate_to_first_inetrval(ba_allwd,0);
    
    is_branch_roots = find_branch(spec.rootss_allwd,spec.params,ba_allwd);
    
    % __ look for 2 roots on 1 branch
    br_double_cross = sum(is_branch_roots,2)==2;
    if any(br_double_cross) 
        % __ the two roots on 1 branch 
        ind_dcr = is_branch_roots(br_double_cross,:);
        contangij = spec.rootss_allwd(ind_dcr);
        
        th_dc = th_allwd_ears ( th_allwd_ears >= min(contangij) &  th_allwd_ears <= max(contangij) );
        th_wulff = th_dc ( wy_shftd(th_dc) >= 0 );
        contangij = check_contangs(contangij);
        
    else
        
        th_wulff = nan;
        contangij = nan(1,2);
    end % if branch double cross
end%

%% 
% - check that contangij are 2 contact angles (Left and Right) and order ascendingly
% - contangij are NORMAL angles in contact points
% - tangent_ang = norm_ang - pi/2
% - tangent ang in (0,pi) corresponds to LEFT contact
% - tangent ang in (-pi,0) corresponds to RIGHT contact
function contangij = check_contangs(contangij)
    
    if any(isnan(contangij)) 
        return
    elseif isempty(contangij)
        contangij
    else
        assert(numel(contangij)==2,['not 2 contact angles in checkout. Num angles=' num2str(numel(contangij))])
        tang_ang = rotate_to_first_inetrval(contangij - pi/2,0);
        check(1) = min(tang_ang) <0;
        check(2) = max(tang_ang) >0;
        assert(all(check),['Not Left-Right pair. contangij=' num2str(tang_ang*180/pi) 'deg. Must be + and -'])

        contangij = [min(contangij), max(contangij)];
    end
end % func