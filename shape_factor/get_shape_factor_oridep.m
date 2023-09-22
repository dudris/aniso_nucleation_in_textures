% COPYRIGHT NOTICE 
 % This file is part of a dataset <Minar, Martin (2023), “Influence of surface energy anisotropy on nucleation and crystallographic texture of polycrystalline deposits”, Mendeley Data, V1, doi: 10.17632/bsdff8shbz.1>, coupled to publication of the same name by Minar, Moelans submitted to Computational Materials Science in September 2023. 
 % Distributed under GPLv3 license.
function out = get_shape_factor_oridep(in,outputspec)

%     addpath('../../solver/','../../results_processing/','../')
    
    addtoresults = outputspec{1};
    if addtoresults
        resultsfile = outputspec{2};
        if exist(resultsfile,'file')
            load(resultsfile)
            numentries = length(results);
            newentry = numentries+1;
        end
    end
    
    params.nfold = in.n;
    params.Omega = in.Omega;
    params.soaIE = params.Omega/(params.nfold^2-1);
    params.codeIEaniso = 'IEanisofun_1';

    [TORI,BORI] = meshgrid(in.tori,in.bori);
    
    % __ need to assure periodicity of RS
    % __ https://en.wikipedia.org/wiki/Modulo ; modulo with rounded division convention
    mymod = @(x,y) x - round(x./y).*y ;
    misor = mymod(BORI-TORI,2*pi/in.n);

%     % __ need to assure periodicity of RS
    rsGBE = calc_rsGBE(in.GBE,in.misor_0,misor);

    
    params.offset_ang = 0;

    % sSLE ... substrate Solid-liquid energy
    % instead of rotating the anisofun I rotate here the interface in
    % equivalent way  (checked)
    [fIEijfun,~,~] = AssignAnisotropyFunction(params,false);
    sSLE =fIEijfun(pi/2-BORI,params.soaIE);

    Cshift = sSLE-rsGBE./in.SLE;

    stabcond =  cos(params.nfold*(pi/2-BORI))<(rsGBE./in.SLE+1)/params.Omega;
    
    [A,S,contangs,h,sol_type] = parser_get_shape_factor(TORI,Cshift,params.nfold,params.Omega,false);
    
    stable_points = find(stabcond==1);
    sp = stable_points(1);
    dummy_spec = gen_input_stable_sol(TORI(sp),Cshift(sp),in.n,in.Omega);
    A_hom = dummy_spec.A_hom;
    clear dummy_spec sp stable_points
    
    % notice transpose of the data to have bori on x axis in plotting
    res_data = {S', 'S'; nan, 'bori_check_allwd'; datestr(clock),'timestamp'; sol_type' , 'sol_type'; stabcond', 'stab_cond'
        A_hom, 'A_hom';permute(contangs,[2,1,3]),'contangs';Cshift','Cshift';h','h';params.soaIE,'soaIE'};
    
    out = in;

    for kk = 1:length(res_data)
        out.(res_data{kk,2}) = res_data{kk,1};
    end
    
    if addtoresults
        if exist(resultsfile,'file')
            [results,out] = CreateNonexistentFields(results,out);
            results(newentry) = out;
        else
            results = out;
        end
        save(resultsfile,'results');
        out = results;
    end
end % func

%%
function [results,currententry] = CreateNonexistentFields(results,currententry)
    fieldnamesOLD = fieldnames(results);
    fieldnamesNEW = fieldnames(currententry);
    uniqueNEW = setdiff(fieldnamesNEW,fieldnamesOLD);
    uniqueOLD = setdiff(fieldnamesOLD,fieldnamesNEW);
    
    for k = 1:length(uniqueNEW) % the missing field
        for row = length(results) % the number of rows of results - that mamy times the new field must be added
            results(row).(uniqueNEW{k}) = [];
        end
    end
    
    for k = 1:length(uniqueOLD)
        currententry.(uniqueOLD{k}) = [];
    end

end