function [spe, ppd] = cpm_prf_get_ppd(PRF, xy, idx, nsamp, posterior)
% Sample from the prior (or posterior) of a pRF model and run samples 
% through the likelihood, to give the prior (or posterior) predictive 
% density (PPD).
%
% PRF   - estimated PRF
% xy    - coordinates in stimulus space to sample
% idx   - index of the PRF within the PRF structure
% nsamp - the number of samples to draw
%
% Returns:
% spE      - samples from the prior [parameters x samples]
% sEp      - samples from the posterior [parameters x samples]
% prior_pd - average sampled PRF response (prior predictive density)
% prior_pd - average sampled PRF response (posterior predictive density)
%
% ---------------------------------------------------------------------
% Copyright (C) 2016 Peter Zeidman
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.   
% ---------------------------------------------------------------------
% Edited SRSteinkamp 2022

if nargin < 4
    nsamp = 1;
end

if nargin < 5
    posterior = true;
end

% Sample from priors and posteriors
% -------------------------------------------------------------------------

pE_struct = PRF.M.pE{idx};

% Get priors and posteriors
if posterior
    ep = spm_vec(PRF.Ep{idx});
    cp = full(PRF.Cp{idx});
else
    ep = spm_vec(PRF.M.pE{idx});
    cp = spm_vec(PRF.M.pC{idx});
end

if isvector(cp), cp = diag(cp); end

% Sample
spe = spm_normrnd(ep, cp, nsamp);

if nargout < 2
    return;
end

% Integrate model with sampled parameters
% -------------------------------------------------------------------------

% Create progress window

ppd  = 0;
p_bad = 0;

for i = 1 : nsamp    
        
    % Integrate model using prior sample
    spe_struct = spm_unvec(spe(:,i), pE_struct);
    
    g = feval(PRF.M.IS, spe_struct, PRF.M, PRF.U, 'get_response', xy);

    if all(isfinite(g))
        ppd = ppd + g;
    else
        p_bad = p_bad + 1;
    end
    
end

fprintf('Bad samples ppd %d\n', p_bad);

% Average
ppd = ppd .* (1 /  (nsamp - p_bad));