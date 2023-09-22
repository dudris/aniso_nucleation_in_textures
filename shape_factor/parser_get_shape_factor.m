% COPYRIGHT NOTICE 
 % This file is part of a dataset <Minar, Martin (2023), “Influence of surface energy anisotropy on nucleation and crystallographic texture of polycrystalline deposits”, Mendeley Data, V1, doi: 10.17632/bsdff8shbz.1>, coupled to publication of the same name by Minar, Moelans submitted to Computational Materials Science in September 2023. 
 % Distributed under GPLv3 license.
%% get_shape_factor
% IN
%   TOPROT ... top grain orientation in rad
%   CYCOORD ... shift of the Wulff center due to wetting condition 
%   n ... nfold
%   Omg ... normalized strength of anisotropy
% OUT
%   A ... area of the heterogeneous nucleus, size(A) = [length(TOPROT),length(CYCOORD)]
%   S ... shape factor, S = A/A_hom, A_hom is untruncated Wulff, size(S) = [length(TOPROT),length(CYCOORD)]
%   A_hom ... area of the homogeneous nucleus
%   contangs ... contact angles


% function [A,S,A_hom,contangs, h,sol_type] = parser_get_shape_factor(TOPROT,CYCOORD,n,Omg,stabcond,plotting)
function [A,S,contangs, h,sol_type] = parser_get_shape_factor(TOPROT,CYCOORD,n,Omg,plotting)
    
    disp('****************************************************************************')
    disp([mfilename ' function initiated' ])
    status = '';
    
    if any(size(TOPROT)==1) && any(size(CYCOORD)==1) % input as column vectors
        dims = [length(TOPROT),length(CYCOORD)];
        disp(['Dimensions of the result: ' num2str([dims(2), dims(1)])])
        disp('Progress:')
        A = zeros(dims);
        S = A;
        contangs(:,:,1) = A;
        contangs(:,:,2) = A;
        sol_type = A;
        for kk = 1:length(TOPROT)
            for kkk = 1:length(CYCOORD)
                toprot = TOPROT(kk);
                cycoord = CYCOORD(kkk);
                [A(kk,kkk),S(kk,kkk),contangs(kk,kkk,:),h(kk,kkk),sol_type(kk,kkk)] = get_shape_factor(toprot,cycoord,n,Omg,plotting);
%                 if stabcond(kk,kkk)
%                     [A(kk,kkk),S(kk,kkk),contangs(kk,kkk,:),h(kk,kkk),sol_type(kk,kkk)] = get_shape_factor(toprot,cycoord,n,Omg,plotting);
%                 else
%                     [A(kk,kkk),S(kk,kkk),contangs(kk,kkk,:),h(kk,kkk),sol_type(kk,kkk)] = no_solution_output;
%                 end
            end % for CYCOORD
            disp(['    ' num2str(kk/dims(1)*100) ' %'])
        end % for TOPROT

    elseif all(size(TOPROT)==size(CYCOORD))
        dims = size(TOPROT);
        disp(['Dimensions of the result: ' num2str([dims(2), dims(1)])])
        A = zeros(dims);
        S = A;
        contangs(:,:,1) = A;
        contangs(:,:,2) = A;
        sol_type = A;
        for kk = 1:dims(1)
            for kkk=1:dims(2)
                toprot = TOPROT(kk,kkk);
                cycoord = CYCOORD(kk,kkk);
                [A(kk,kkk),S(kk,kkk),contangs(kk,kkk,:),h(kk,kkk),sol_type(kk,kkk)] = get_shape_factor(toprot,cycoord,n,Omg,plotting);
%                 if stabcond(kk,kkk)
%                     [A(kk,kkk),S(kk,kkk),contangs(kk,kkk,:),h(kk,kkk),sol_type(kk,kkk)] = get_shape_factor(toprot,cycoord,n,Omg,plotting);
%                 else
%                     [A(kk,kkk),S(kk,kkk),contangs(kk,kkk,:),h(kk,kkk),sol_type(kk,kkk)] = no_solution_output;
%                 end
            end % for kkk
            revertstr = repmat('\b',1,length(status));
            status = sprintf('Percent done: %3.1f \n', kk/dims(1)*100);
            fprintf([revertstr, status  ])
        end% for kk
        
    end % if input
    
    
    
end% func


%%
function [Aij,Sij,contangij,hij,sol_type_ij] = no_solution_output
    Aij= nan;
    Sij = nan;
    hij= nan;
    contangij = nan(1,2);
    sol_type_ij = -3; % no solution, unstable shape
end