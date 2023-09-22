% COPYRIGHT NOTICE 
 % This file is part of a dataset <Minar, Martin (2023), “Influence of surface energy anisotropy on nucleation and crystallographic texture of polycrystalline deposits”, Mendeley Data, V1, doi: 10.17632/bsdff8shbz.1>, coupled to publication of the same name by Minar, Moelans submitted to Computational Materials Science in September 2023. 
 % Distributed under GPLv3 license.
function [numgr,s_firstrow ] = suppl_MC_initialize(oris,sx,MCspec)
    
    spec = MCspec.init_code;
    numoris = length(oris);
    
    if strcmp(spec,'rand')
        numgr=sx;
        s_firstrow = ceil(50*rand(1,sx));

    elseif strcmp(spec,'3pm2') %mean grain width=3px, mrsd 2
        sumxpx=0;
        numgr = 0;
        while sumxpx<sx 
            pxbin = ceil(3 + 2.*randn);
            if sumxpx+pxbin<=sx && pxbin>0
                colind = sumxpx+[1:pxbin];
                s_firstrow(1,colind) = int8(ceil(rand*numoris));
                sumxpx = sumxpx+pxbin;
                numgr = numgr +1;
            end
        end %while
    
    elseif strcmp(spec,'rand_31pts_center')
        numgr=sx;
        % assuming 50 orientations
        s_firstrow = 10+ceil(30*rand(1,sx));
        
    elseif strcmp(spec,'rand_31pts_edges')
        numgr=sx;
        % assuming 50 orientations
        s_firstrow = 35 + ceil(30*rand(1,sx));
        over_interval = s_firstrow>50;
        s_firstrow(over_interval) = s_firstrow(over_interval)-50;
%         histogram(s_firstrow)
        
        
    elseif strcmp(spec,'2grains') % one with maximal
        maxaniso = unique(max(MCspec.aniso(pi/2,oris)));
        ind_maxaniso = find(MCspec.aniso(pi/2,oris)==maxaniso,1);
        minaniso = unique(min(MCspec.aniso(pi/2,oris)));
        ind_minaniso = find(MCspec.aniso(pi/2,oris)==minaniso,1);
        
        halfind = 1:sx<ceil(sx/2);
        numgr = 2;
        
        if minaniso~=maxaniso    
            s_firstrow(halfind) = ind_maxaniso;
            s_firstrow(~halfind) = ind_minaniso;
        else
            s_firstrow(halfind) = ceil(rand*50);
            s_firstrow(~halfind) = ceil(rand*50);
        end
        
%         disp(['2 grains init.'])
%         disp(['orientations: [' num2str([ind_minaniso ,ind_maxaniso]) ']'])
%         disp(['aniso: [' num2str([minaniso ,maxaniso]) ']'])
    
    elseif strcmp(spec,'2grains1px') % one with maximal
        maxaniso = unique(max(MCspec.aniso(pi/2,oris)));
        ind_maxaniso = find(MCspec.aniso(pi/2,oris)==maxaniso,1);
        minaniso = unique(min(MCspec.aniso(pi/2,oris)));
        ind_minaniso = find(MCspec.aniso(pi/2,oris)==minaniso,1);
        
        s_firstrow(1,1:sx) = ind_maxaniso;
        s_firstrow(1,ceil(sx/2)) = ind_minaniso;
        numgr = 2;
        
    elseif strcmp(spec,'1grain')
        numgr=1;
        if strcmp(MCspec.init_spec,'max')
            maxaniso = unique(max(MCspec.aniso(pi/2,oris)));
            ind_maxaniso = find(MCspec.aniso(pi/2,oris)==maxaniso,1);
            s_firstrow(1,1:sx) = ind_maxaniso;
        else
            s_firstrow = ceil(50*rand) ;% randomly selected orientation;
        end
        
    end % if spec

end% func