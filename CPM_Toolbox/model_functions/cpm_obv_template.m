function [y,pE,pC] = cpm_obv_template(varargin)
%CPM_OBV_IDENTITY This function defines the observation model and its
% priors
%   The observation models further transforms the weighted pRF signal into the final
%   measurable signal which will be fitted to actual data. 
%
%   Here, the user must define observation_function and observation_priors
%
%   cpm_obv_int uses an extdended Balloon hemodynamic model which is adapted from the
%   Ziedman study. This particular model is powerful but it's also
%   significantly slower than other methods.
%
%   cpm_obv_identity is the nul model, it downsamples the pRF signal down to nscans. 
%   This is used in case you want to precompute the paramaters of the observation model 
%   along with the neuronal parameters.
%   
%   You can easily modify cpm_obv_identity to create a fast canonical HRF model
%   by adding a convolution with the canonical HRF just before the
%   downsampling.
%


%% ignore this
if nargin ==3
    
    P = varargin{1};
    M = varargin{2};
    U = varargin{3};
    
    y = observation_function(P,M,U);
else
    y = [];
    [pE,pC] = observation_priors() ;
    
end

end


function [y] = observation_function(P,M,Z)
% function that transform pRF signal (z) into final prediction (y)
    
    % your observation function here

        % M.f = some user defined dynamics (defaults to Balloon model)
        % M.g = some user defined non linearity (defaults to Balloon model)
        % y=spm_int_sparse(P,M,Z)

        % for BOLD precomputation:
        % n = int32(floor(length(Z.u)/M.ns));
        % y = downsample(Z.u,n);

        % for canonical HRF:
        % Z.u = convolve Z.u with canonical HRF
        % n = int32(floor(length(Z.u)/M.ns));
        % y = downsample(Z.u,n);

    
end

function [pE,pC] = observation_priors()
% function that defines latent priors pE and pC for the observation model
% (transit decay and epsilon for Balloon model)
% these parameters are used inside M.f and M.g (see spm_prf_fx and
% spm_prf_gx as an example)



end



