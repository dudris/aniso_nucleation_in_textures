%% tori = rand_tori(tori_all,tori_prob)
% - numeric sampling from pdf specified in 'tori_prob'
% - principle: inversion and interpolation of CDF (as in
% https://www.av8n.com/physics/arbitrary-probability.htm)
% - visualize the distrubution and test the sampling by setting 'testfun=true' in the script
% 
% INPUT
%   - tori_all ... vector of available top grain oientations
%   - tori_prob ... vector of probability corresponding to points tori_all
%       - does not have to be normalized
% OUTPUT
%   - tori ... top grain orientation sampled from pdf 'tori_prob'
% 
% EXAMPLEs
% x= linspace(-5,5,100);
% y = 1./(1+(3*x).^2);
% y = a*(x-3).^2;
% y = cosh(x);
% y = sqrt(max(x)^2-x.^2);
% y = cos(x/max(x)*pi/2);
% P = peaks(100);
%     ind = 70;
%     y = P(:,ind)-min(P(:,ind));
% 
% tori = rand_tori(x,y);

function tori = rand_tori(tori_all,tori_prob)
% plot(tori_prob,'o')
    
    tori_spacing = tori_all(2)-tori_all(1);
    normfac = 1/trapz(tori_prob)/tori_spacing;
    % normalize to have the pdf
    yn = normfac*tori_prob;
%     plot(yn,'o')
%     trapz(yn)*tori_spacing

    % from definition of cumulative distribution function
    cdf = arrayfun(@(t) trapz(yn(tori_all<=t))*tori_spacing,tori_all);
%     semilogy(distf)
%     plot(distf,tori_all,'.')
    
    [cdf_uq,ind_uq] = unique(cdf);
    tori = interp1(cdf_uq,tori_all(ind_uq),rand);

    testfun=false;
    if testfun
        figure(14)
            pp = rand(10000,1);
            xx = interp1(cdf,tori_all,pp);
            [~,edges] = discretize(xx,100);
            histogram(xx,edges,'Normalization','pdf')
            hold on
            plot(tori_all,yn,'-k')
            hold off
    end % testing
end % func
