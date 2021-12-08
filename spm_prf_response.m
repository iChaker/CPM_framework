function varargout = spm_prf_response(P,M,U,varargin)
% Template pRF response function. This example models neuronal
% activity as a scaled version of the inputs:
%
% z(t) = alpha * u(t)
%
% Where alpha is a parameter estimated from the data.
% 
%
% Inputs:
%
% P      - parameter structure
% M      - model structure
% U      - experimental timing
% action - (optional) the action to performm. If not provided, the
%          predicted BOLD and neuronal timeseries are returned.
%
% -------------------------------------------------------------------------
% FORMAT [y,Z] = spm_prf_fcn_template(P,M,U)
% Return the BOLD and neuronal predicted timeseries
%
% P         parameters
% M,U       model, inputs
%
% y         fMRI time series
% Z         neuronal response
% -------------------------------------------------------------------------
% FORMAT P = spm_prf_fcn_template(P,M,U,'get_parameters')
% Return the given parameters corrected for display
%
% P         parameters
% M,U       model, inputs
% -------------------------------------------------------------------------
% FORMAT S = spm_prf_fcn_template(P,M,U,'get_summary')
% Summarises the pRF with simple (Gaussian) parameters x,y,width,beta
%
% S         structure with fields x,y,width,beta
% M,U       model, inputs
% -------------------------------------------------------------------------
% FORMAT tf = spm_prf_fcn_template(P,M,U,'is_above_threshold',Cp,v)
% Return whether the model with parameters P and covariance Cp passes an
% arbitrary threshold for display
%
% P         parameters
% M,U       model, inputs
% Cp        parameter covariance matrix
% v         voxel index
% -------------------------------------------------------------------------
% FORMAT x = spm_prf_fcn_template(P,M,U,'get_response',xy)
% Return the instantaneous response of the PRF at coordinates xy
%
% P         parameters
% M,U       model, inputs
% xy        [2xn] vector of coordinates to evaluate the PRF
% -------------------------------------------------------------------------
% FORMAT [pE,pC] = spm_prf_fcn_template(P,M,U,'get_priors')
% Return the priors for the model. Importantly, this defines which
% parameters are in the model.
%
% pE        structure or vector of prior expectations
% pC        prior covariance maitrx
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

if nargin < 4
    % Integrate the model over time and return neuronal timeseries z and 
    % BOLD timeseries y. As an example, here we have neuronal model 
    % z(t) = alpha, where alpha is an estimated parameter.
    
    % Number of volumes = inputs
    n  = length(U); 

    % Neural timeseries. A vector with one entry per microtime bin. The 
    % U.nbins field is injected automatially by spm_prf_analyse
    nbins = max(U(1).nbins, max(U(end).ind));
    z = zeros(1,nbins);
    
    
    
    
    W = get_response(P,M,U);
    
        
    for t = 1:n    
        % Microtime index for this volume
        ind = U(t).ind;  
        
        resp = U(t).signals1D' .* W;
        z(ind) = z(ind) + sum(resp);
        
    end
            
    % Integrate the BOLD model
    Z.u=z';
    Z.dt=M.dt;
    y=feval(M.cpm.obv_fun,P,M,Z);   
    
    varargout{1} = y;
    varargout{2} = Z;
else
    % This section of the code provides information on the model, primarily
    % for plotting purposes in spm_prf_review()
    
    action = varargin{1};

    switch action         
        case 'get_parameters'
            % Get the parameters with any corrections needed for 
            % display      
%             try
%                 P.mu_tau
%             catch
%                 lP = P;
%                 P = transform_latent_parameters(lP,U);
%                 try
%                 P.transit = lP.transit;
%                 P.decay   = lP.decay;
%                 P.epsilon = lP.epsilon;
%                 catch
%                 end
%             end
%           
            P = cpm_get_true_parameters(P,M,U);
            varargout{1} = P;
        case 'get_summary'
            % Get a summary of the pRF shape under Gaussian assumptions
            varargout{1} = ...
                struct('x',P.x,'y',P.y,'width',P.width,'beta',P.beta);
        case 'is_above_threshold'
            % Return binary vector identifying whether each voxel is
            % above some threshold for display            
            varargout{1} = 1;
        case 'get_response'            
            % Return the prediction of the model at coordinates xy            
            xy = varargin{2};            
            varargout{1} = ones(1,length(xy));            
        case 'get_priors'
            % Return a structure containing the priors
            pE.alpha = 1;
            pC.alpha = 1;
            varargout{1} = pE;
            varargout{2} = pC;            
        case 'glm_initialize'                        
            % (Optional) Return parameters initialized using some
            % rapid initial search on timeseries y
            y = varargin{2};
            
            varargout = P; 
        otherwise
            error('Unknown action');
    end
end
end

% -------------------------------------------------------------------------
function W = get_response(P,M,U)
% Get PRF response to a stimulus 
%
% P  - parameters
% coords - coordinates (tau,eps,eta)

[true_P, tm_names, ts_names] = cpm_get_true_parameters(P,M,U);

p_names = fieldnames(U(1).grid);

true_parameters = struct();
for i=1:length(tm_names)
    moment_vec = [];
    for j=1:length(p_names)
        VL_true = [ tm_names{i} '_' p_names{j} ];
        moment_vec(j)=true_P.(VL_true);
    end
    true_parameters.(tm_names{i}) = moment_vec;
end

for i=1:length(ts_names);
    true_parameters.(ts_names{i})=true_P.(ts_names{i});
end

receptive_field = M.cpm.RF_fun;
get_W = feval(receptive_field);
coords = U(1).gridpoints;

W = get_W(coords',true_parameters);

end
