%% MCspec = gen_MCspec_assign_uniform_S(MCspec, varargin)
% - assigns a value to MCspec.nucl_mode_uniform_S
% - variable argument contains the value to be used in the uniform shape
% factor to be used in nucleation, overrides other options when used
% - when omitted, assigns either NaN when nucl_mode is 'aniso' and a mean value
% of anisotorpic S, when nucl_mode is 'uniform'
function MCspec = gen_MCspec_assign_uniform_S(MCspec, varargin)
    if ~isempty(varargin)
        MCspec.nucl_mode_uniform_S = varargin{1}; % shape factor in 'uniform' nucleation
        return
    end

    if strcmp(MCspec.nucl_mode, 'aniso')
        MCspec.nucl_mode_uniform_S = nan; % shape factor in 'uniform' nucleation
    elseif strcmp(MCspec.nucl_mode, 'uniform')
        MCspec.nucl_mode_uniform_S = nan; % shape factor in 'uniform' nucleation
%         % __ mean value of anisotropic SF - deprecated
%         load(MCspec.sf_resultsfile,'results')
%         res = results(MCspec.sf_ind);
%         cond = res.stab_cond & ~isnan(res.S);
%         MCspec.nucl_mode_uniform_S = mean(res.S(cond)); % shape factor in 'uniform' nucleation
    end
end