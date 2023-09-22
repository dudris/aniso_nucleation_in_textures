% COPYRIGHT NOTICE 
 % This file is part of a dataset <Minar, Martin (2023), “Influence of surface energy anisotropy on nucleation and crystallographic texture of polycrystalline deposits”, Mendeley Data, V1, doi: 10.17632/bsdff8shbz.1>, coupled to publication of the same name by Minar, Moelans submitted to Computational Materials Science in September 2023. 
 % Distributed under GPLv3 license.
%% function th_m = border_angles(input,plotting)
% th_m is not rotated, should there be some offset, it must be applied on
% the resulting th_m
% outspec ... 'forbb' , 'nconvex', 'earscross'
function th_m = border_angles_v2(input,outspec,plotting)
    
    if input.Omega <= 1
        th_m = 0;
        
    else
        
        load('dtb_border_angles.mat','BA')
        
        [th_m, entry_found] = lookup_in_dtb(input,outspec,BA);
%         entry_found = false;
        if ~entry_found % compute both forbb and nconvex angs and save
            entry = size(BA,1)+1;
            
            % __ forbidden angles
            intf_stiff = @(th) 1 - input.Omega*cos(input.nfold*th);
            % __ to eliminate angles th there intf_stiff(th)<0
            % __ min(intf_stiff) = 1-Omega at th=0
            % __ intf_stiff(pi/nfold)>0 always
            th_ba_forb(1,1) = fzero(intf_stiff,[0, pi/input.nfold]);
            th_ba_forb(1,2) = -th_ba_forb(1,1);
            th_ba_forb = ext_ba_to_other_segs(input,th_ba_forb);
            
            % __ non-convex angles
            th_ba_ncvx = border_angles_segm1(input,plotting);
            th_ba_ncvx = ext_ba_to_other_segs(input,th_ba_ncvx);
           
            % __ ears crossing angles
            th_ears_cross = find_ears_cross(input,th_ba_forb,plotting);
            
            % __ save to dtb
            
            BA = [BA; cell2table({input.soaIE, input.nfold, input.Omega, th_ba_forb, th_ba_ncvx,th_ears_cross},'VariableNames',{'soaIE','nfold','Omega','ba_forbb','ba_nconvex','ears_cross'})];
%             BA(end+1,:) = cell2table({input.soaIE, input.nfold, input.Omega, th_ba_forb, th_ba_ncvx,th_ears_cross})
%             BA(end+1,:) = {input.soaIE, input.nfold, input.Omega, th_ba_forb, th_ba_ncvx,th_ears_cross};
%             BA{end+1,:} = {input.soaIE, input.nfold, input.Omega, th_ba_forb, th_ba_ncvx,th_ears_cross};

%             BA.ba_forbb{entry} = th_ba_forb;
%             BA.ba_nconvex{entry} = th_ba_ncvx;
%             BA.soaIE(entry) = input.soaIE;
%             BA.nfold(entry) = input.nfold;
%             BA.Omega(entry) = input.Omega;
            
%             save('tmp','BA')
            save('dtb_border_angles','BA')
%             disp('Entry not found. Border angles computed and saved to database.')
            
            if strcmp(outspec,'forbb')
                th_m = th_ba_forb;
            elseif strcmp(outspec,'nconvex')
                th_m = th_ba_ncvx;
            elseif strcmp(outspec,'earscross')
                th_m = th_ears_cross;
            end% switch outspec
            
        end % if entry not found
        
    end % if Omega <=1
    
end % func main


%% border_angles_segm1
% assumes no offset angle
function th_m = border_angles_segm1(input,plotting)
    
    nfold = [input.nfold];
    nfold = nfold(1);
    delta = [input.soaIE];
    delta = delta(1);
    assert(all(([input.nfold]-nfold)==0),'too many inclination dpendent interfaces specified. Check also border_angles>border_angles_segm1.')
    assert(all(([input.soaIE]-delta)==0),'too many inclination dpendent interfaces specified. Check also border_angles>border_angles_segm1.')
    
    syms th tangline(th) dtangline(th)
    tangline = cos(th) / (1+delta*cos(nfold*th) );
    dtangline = diff(tangline);

    g= matlabFunction(tangline);
    dg= matlabFunction(dtangline);
    
    th_m(1,1) = fzero(dg,-pi/nfold); % first border angle negative
    th_m(1,2) = -th_m(1,1); % second to be positive
    
    if plotting 
        ang_all = linspace(-pi/nfold,pi/nfold,300)';
        figure(22)
            plot(ang_all*180/pi,g(ang_all),'.',ang_all*180/pi,dg(ang_all),'o')
            grid on
            hold on
            plot(th_m*180/pi,g(th_m),'dg')
%             plot(ang*180/pi,polymatrix*coeff,'LineWidth',1.5)
%             plot([th_m;th_m]*180/pi,[g(th_m);dg(th_m)],'kd')
%             legend('cos(th)/f(th)','d/d th ( cos(th)/f(th) )','fit poly 9','maxima','Location','south')
            xlim(180/nfold*[-1,1])
            xticks([-180/nfold,0,180/nfold])
            xticklabels({'-\pi/n','0','\pi/n'})
            hold off
            set(gca,'FontSize',14)
    end
    
end % func border_angles_segm1

%%
function th_m_extd = ext_ba_to_other_segs(input,th_m)
    th_m_extd = th_m;
    
    nfold = input.nfold;
%         nfold = nfold(1);
    is_nfold_odd = mod(nfold ,2)~=0;
    segment_width = 2*pi/nfold;
    
    if is_nfold_odd
        segments_in_halfspace = floor(nfold/2);
        for k = 1:segments_in_halfspace
            th_m_extd(k+1,:) = th_m(1,:) + k*segment_width;
            th_m_extd(segments_in_halfspace+k+1,:) = th_m(1,:) - k*segment_width;
        end
        th_m_extd = sort(th_m_extd,2,'ascend');
    else % nfold is even
        segments_in_halfspace = nfold/2 -1;
        for k = 1:segments_in_halfspace
            th_m_extd(k+1,:) = th_m(1,:) + k*segment_width;
            th_m_extd(segments_in_halfspace+k+1,:) = th_m(1,:) - k*segment_width;
        end
        th_m_neg = th_m_extd(1,th_m_extd(1,:)<0);
        th_m_pos = th_m_extd(1,th_m_extd(1,:)>0);
        % to the last rotation I must add length of the ALLOWED interval 
        allowed_segment_width = 2*(pi/nfold-th_m_pos);
        th_m_extd(nfold,1) = th_m_neg - segment_width*segments_in_halfspace - allowed_segment_width;
        th_m_extd(nfold,2) = th_m_pos + segment_width*segments_in_halfspace + allowed_segment_width;
%         th_m(nfold,:) = th_m(1,:) + nfold/2*segment_width;
    end % if is_nfold_odd
    
end% func

%%
function [th_m, entry_found] = lookup_in_dtb(input,outspec,BA)

    % looks for soa and nfold
    is_entry_present = abs(BA.soaIE-input.soaIE)<1e-6 & abs(BA.nfold - input.nfold)<0.1;
    entry_found = any(is_entry_present);
    
    if entry_found

        entry = find(is_entry_present,1); %takes the first matching entry
        
        if strcmp(outspec,'forbb')
            th_m = BA.ba_forbb{entry};

        elseif strcmp(outspec,'nconvex')
            th_m = BA.ba_nconvex{entry};
            
        elseif strcmp(outspec,'earscross')
            th_m = BA.ears_cross{entry};

        end% switch outspec

%         disp('Entry found. Border angles read from database.')

    else % not found, to compute elsewhere
        th_m = [];
        return
        
    end % if entry found

end% func

%%
function th_ears_cross = find_ears_cross(input,ba_forbb,plotting)
    n = input.nfold;
    soaIE = input.soaIE;
    if abs(soaIE-0.62)<1e-3
        soaIE;
    end
    Omg = input.Omega;
    
    if ~any(n==[2,3,4,6])
        warning('Unexpected order of symmetry for crossed ears solution. Only even symmetries optimized for 1st and 2nd order solutions.')
        warning('The sub-function ''find_ears_cross''  exited.')
        th_ears_cross = nan;
        return
    end
    
    Omg_lim_fourfold = 9;
    Omg_lim1_sixfold = 5.7477270848675;
    Omg_lim2_sixfold= 21.7477270848675;

    if n<=3 || ( n==4 && Omg<Omg_lim_fourfold )  || ( n==6 && Omg<Omg_lim1_sixfold ) || n==5 || n>6
        th_ears_cross = nan;
        return

    else
        % __ to bring an earcross on x axis
        rot_earcross(1) = 2*pi/2/n; % 1st order solution, at least for the 4- and 6-fold
        if n ==6 && ( Omg>=Omg_lim2_sixfold) 
            rot_earcross(2) = 0; % 2nd order ears cross solution
        end
        
        th_ears_cross = zeros(n,2,length(rot_earcross));
        
        for ll = 1:length(rot_earcross)
            f = @(th) 1+soaIE*cos(n*(th-rot_earcross(ll)));
            df = @(th) -n*soaIE*sin(n*(th-rot_earcross(ll)));
            wx = @(ang) f(ang).*cos(ang) - df(ang).*sin(ang);
            wy = @(ang) f(ang).*sin(ang) + df(ang).*cos(ang);

            %__ see all to right from isolated Wulff 
            ba_forbb_rot = rotate_to_first_inetrval(ba_forbb,rot_earcross(ll));
            [th_allwd, ba_allwd] = GetAllowedAngles(ba_forbb_rot,360,false);
            ba_allwd = rotate_to_first_inetrval(ba_allwd,0);
            
            % __ select the 2 branches reaching farthest to the +x dir
            wx_max = max(max(wx(ba_allwd)));
            wx_max = wx_max(1); % in case of degenerate value
            % __ the two ears crossing on x axis
            ind_cross_branch = any(abs(wx(ba_allwd)-wx_max)<0.01*wx_max,2); 
            % __ due to symmetry is 1 branch enough
            % __ top cross branch to have both border angles > 0 and end in negative y-value
            if ll == 1 && n == 6  && Omg < Omg_lim2_sixfold
                ind_top_cross_branch = all(ba_allwd>0,2) & ind_cross_branch & ...
                            wy(min(ba_allwd,[],2))<=0 & wy(max(ba_allwd,[],2))>=0;
            elseif ll == 1 && n == 6  && Omg >= Omg_lim2_sixfold
                ind_top_cross_branch = all(ba_allwd>0,2) & ind_cross_branch & ...
                            wy(max(ba_allwd,[],2))<=0 & wy(min(ba_allwd,[],2))>=0;
            else
                ind_top_cross_branch = all(ba_allwd>0,2) & ind_cross_branch;
            end

            th_br = linspace(ba_allwd(ind_top_cross_branch,1),ba_allwd(ind_top_cross_branch,2),100);
%                 plot(wx(th_allwd),wy(th_allwd),'.'), axis equal, grid on
%                 hold on, plot(wx(th_br),wy(th_br),'o')
            
            % __spline fit
            gxmx_spl1 = spline(th_br,wy(th_br));
            earcross = fnzeros(gxmx_spl1);
            % __ the single root
            earcross = earcross(1,1);

            % __ the other roots from symmetry
            th_ears_cross(:,:,ll) = earcross*[-ones(n,1),ones(n,1)] + 2*pi/n*(0:(n-1))';

            th_ears_cross(:,:,ll) = rotate_to_first_inetrval(th_ears_cross(:,:,ll),-rot_earcross(ll));
%             th_ears_cross = rotate_to_first_inetrval(th_ears_cross,0);
        end % for num of orders of crossears solution
        
        if plotting
            rot = 0;
            f = @(th) 1+soaIE*cos(n*(th-rot));
            df = @(th) -n*soaIE*sin(n*(th-rot));
            wx = @(ang) f(ang).*cos(ang) - df(ang).*sin(ang);
            wy = @(ang) f(ang).*sin(ang) + df(ang).*cos(ang);
            [th_allwd, ~] = GetAllowedAngles(ba_forbb,360,false);
            figure(1)
                plot(wx(th_allwd),wy(th_allwd),'.')
                axis equal
                hold on
                for ll = 1:size(th_ears_cross,3)
                    plot(wx(th_ears_cross(:,1,ll)),wy(th_ears_cross(:,1,ll)),'o','linewidth',1.5)    
                end
                hold off
        end
        
    end
end% func