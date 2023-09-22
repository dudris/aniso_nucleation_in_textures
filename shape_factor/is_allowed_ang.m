% COPYRIGHT NOTICE 
 % This file is part of a dataset <Minar, Martin (2023), “Influence of surface energy anisotropy on nucleation and crystallographic texture of polycrystalline deposits”, Mendeley Data, V1, doi: 10.17632/bsdff8shbz.1>, coupled to publication of the same name by Minar, Moelans submitted to Computational Materials Science in September 2023. 
 % Distributed under GPLv3 license.
%% leaves out non-allowed angles
% allwd_ang = is_allowed_ang(ang,n,Omg,toprot,spec)
%   spec ... 'nconvex', 'forbb'
% 
function allwd_ang = is_allowed_ang(ang,n,Omg,toprot,spec)

    if strcmp(spec,'nconvex')
        % __in a piece-wise fashon rule out all angles from forbidden intervals
        params.Omega = Omg;
        params.nfold = n;
        params.soaIE = Omg/(n^2-1);
        ba_ncnvx = border_angles_v2(params,'nconvex',false);
        ba_ncnvx = rotate_to_first_inetrval(ba_ncnvx,toprot);
        ba_ncnvx = sort(ba_ncnvx,2);
%         figure,polarplot(ba_ncnvx(:),ones(size(ba_ncnvx(:))),'o')
        nojump = diff(ba_ncnvx,[],2)<pi;
        nojump_ind = find(nojump);
        is_forb_cx = false(size(ang));
        for ll = 1:length(nojump_ind)
            kk = nojump_ind(ll);
            is_forb_cx = is_forb_cx | (ang > ba_ncnvx(kk,1) ) & (ang < ba_ncnvx(kk,2) );
        end
        
        if any (~nojump)
            % __ assumes single jump
            is_forb_cx = is_forb_cx | (ang < ba_ncnvx(~nojump,1) ) ;
            is_forb_cx = is_forb_cx | (ang > ba_ncnvx(~nojump,2) );
        end
        is_allwd_cx = ~is_forb_cx;
        
        if Omg >1
            lim_ang = 5e-4;
            for kk=1:n
                is_allwd_cx = is_allwd_cx | abs(ang-ba_ncnvx(kk,1))<lim_ang | abs(ang-ba_ncnvx(kk,2))<lim_ang ;
            end
        end
        
        allwd_ang = ang(is_allwd_cx);
        
    elseif strcmp(spec,'forbb')
        % __ interface stiffness to be positive
        allwd_ang = ang( (1-Omg*cos(n*(ang-toprot)))>=0 );
        
    end % switch spec
    
end% func