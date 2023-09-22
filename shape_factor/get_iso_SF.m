% COPYRIGHT NOTICE 
 % This file is part of a dataset <Minar, Martin (2023), “Influence of surface energy anisotropy on nucleation and crystallographic texture of polycrystalline deposits”, Mendeley Data, V1, doi: 10.17632/bsdff8shbz.1>, coupled to publication of the same name by Minar, Moelans submitted to Computational Materials Science in September 2023. 
 % Distributed under GPLv3 license.
function S_iso = get_iso_SF(in)

    [TORI,BORI] = meshgrid(in.tori,in.bori);
    
    % __ need to assure periodicity of RS
    % __ https://en.wikipedia.org/wiki/Modulo ; modulo with rounded division convention
    mymod = @(x,y) x - round(x./y).*y ;
    misor = mymod(BORI-TORI,2*pi/in.n);

    rsGBE = calc_rsGBE(in.GBE,in.misor_0,misor);
    
    Cshift = ones(size(TORI))-rsGBE./in.SLE;
    alph = 2*acos(Cshift);
    S_iso = (alph - sin(alph))/(2*pi);
end % func


% test 
% clear
% load('Results_shape_factors3.mat')
% imagesc(get_iso_SF(results(62))),set(gca,'YDir','normal'), colorbar

