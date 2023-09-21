%%
% sometimes wy(th) has extremum on 0 and probably two super close but
% different roots are found
% when found, one will be removed
% removing when 2 roots closer than 1e-4 rad
% happened with n=4, soa=0.4, tori=-30/180*pi, yshift=0.8
% __ orijump is either corner of isolated Wulff or cross-ear of 2 branches
function rootss = remove_doubled(rootss)    

%     % __ differnce between x-coord of wulff points in rootss
%     wx_roots = wx(rootss);
%     distmat = abs(wx(rootss)-wx(rootss'));
%     tooclose_pairwise = distmat<3e-3 - eye(numel(rootss));

    % __ difference between angles - more robust when corner is NEAR 0 and
    % non-convex allowed angle gets confused with convex allowed angle =>
    % identification of truncated Wulff solution may be hampered
    lim_angdist = 5e-3;
    distmat = triu(abs(rootss-rootss')) + tril(ones((numel(rootss)))); 
    tooclose_pairwise = distmat <lim_angdist   ... 
                                 |  ( abs(distmat-2*pi) <lim_angdist) ; %  when angle jump near 0
    
    if any(any(tooclose_pairwise))
        [J,~] = ind2sub(size(tooclose_pairwise),find(tooclose_pairwise));
        rootss(J)=[]; 
    end
    
end% function
