function PRF = cpm_specify(SPM,options,y,XYZmm,U,population_field,observation_function,outpath,varargin)
% cpm_specify  specifies a PRF model over a precomputed parameter space
% using a  population receptive field model
%   INPUTS:
%       -SPM: SPM structure
%       -options: structure contaning 'name', 'TE'... See spm_prf_analyse
%       for other options
%       -y: fMRI timeseries
%       -XYZmm: voxel locations
%       -U: the population response over which we will estimate the
%       population field
%       -population_field: name of population field model 'cpm_RF_Gaussian'
%       or 'cpm_RF_SoG'. Alternatively, you can implement your own 
%       population field model by following 'cpm_RF_template'  
%       -observation_function: name of observation function (usually
%       hemodynamic function) 'spm_obv_int' is the default, it contains an
%       extended Balloon model as implemented by the BayespRF framework.
%       Alternatively, you can implement your own observation model by 
%       following 'cpm_obv_template'  
%       -outpath: folder in which the PRF will be saved
%       -varargin: optional parameters to be passed into PRF.M.cpm.other
%   
%   OUTPUTS:
%       -PRF: BayespRF PRF structure, the default PRF model is 'spm_prf_response'. The
%       user may implement 'get_summary' function inside spm_prf_response
%       in order to have full access to feautures provided by BayespRF.
%


%SPM
SPM.swd = outpath;

% VOI
xY.y = y;
xY.XYZmm=XYZmm;
VOI.Y=y; % no eigen summary method
VOI.xY=xY;

% priors -------------------------------------------------
if ~isfield(options, 'mu')
    warning("Assuming that maximal parameter values (grid) are equal to probable paramter space (mu).")
    p_names = fieldnames(U(1).grid);
    for ii = 1 : length(p_names)
        options.mu.(p_names{ii}) = [U(1).grid.(p_names{ii})(1), U(1).grid.(p_names{ii})(2)];
    end

elseif isfield(options, 'mu')
    
    if length(fieldnames(options.mu)) ~= length(fieldnames(U(1).grid))
        error("Check for not fullyspecified mu options not yet implemented")
    end
    
end

if isempty(population_field)
    population_field='cpm_RF_Gaussian';
end

[~, ~, ~, fn_lmoment_priors, fn_lscaling_priors] =  feval(population_field);
[mpE,mpC] = fn_lmoment_priors();
[spE,spC] = fn_lscaling_priors();
pE = struct();
pC = struct();
p_names = fieldnames(U(1).grid);
m_names = fieldnames(mpE);
s_names = fieldnames(spE);

for pidx = 1 : length(p_names)
    model_mu = options.mu.(p_names{pidx});
    grid_mu  = U(1).grid.(p_names{pidx});
    
    [sigma_min, sigma_max,  latent_prior_sigma] = define_sigma(grid_mu, model_mu);
   
    options.sigma.(p_names{pidx}) = [sigma_min, sigma_max];
    options.lsigma_pE.(p_names{pidx}) = latent_prior_sigma;
end
     

for j=1:length(m_names)
    for i=1:length(p_names)
        VLparam = [ m_names{j} '_' p_names{i} ];
        if contains(m_names{j}, 'sigma')
            pE.(VLparam) = options.lsigma_pE.(p_names{i});
        else
            pE.(VLparam) = mpE.(m_names{j});
        end
        pC.(VLparam) = mpC.(m_names{j});
    end
end

for i=1:length(s_names)
    pE.(s_names{i}) = spE.(s_names{i});
    pC.(s_names{i}) = spC.(s_names{i});
end

%options
if isempty(observation_function)
    observation_function='cpm_obv_identity';
end

[~,pE_obv,pC_obv] = feval(observation_function);
o_names= fieldnames(pE_obv);
for i=1:length(o_names)
    pE.(o_names{i}) = pE_obv.(o_names{i});
    pC.(o_names{i}) = pC_obv.(o_names{i});
end

options.pE = pE;
options.pC = pC;
try options.model; 
    
catch
   options.model  = 'spm_prf_response';
end
options.voxel_wise = true;
cpm.obv_fun=observation_function;
cpm.RF_fun=population_field;
cpm.mu = options.mu;
cpm.sigma = options.sigma;
try cpm.other=varargin; catch end
options.cpm=cpm;


PRF  = spm_prf_analyse('specify',SPM,VOI,U,options);

end


function [sigma_min, sigma_max, prior_sigma_latent] = define_sigma(grid_mu, model_mu)
    % Grid_mu = Grid over parameters (e.g. U(1).grid.eta = [-1, 1, 10]
    % Model_mu = Reasonable parameter space.
    
    % Find r the minimal radius given the resolution of the grid over P (e.g.
    % covering 95 % of the space around the space.
    
    % Rewriting variables for understanding:
    grid_min = grid_mu(1);
    grid_max = grid_mu(2);
    grid_n = grid_mu(3);
    
    mu_min = model_mu(1);
    mu_max = model_mu(2);
    
    % Check if symmetric:
    grid_dist = grid_max - grid_min;
    mu_dist = mu_max - mu_min;
    
    if abs((grid_max - grid_dist / 2) - (mu_max - mu_dist / 2)) > 1e-10
        warning('Currently only defined under the assumption of symmetric parameter spaces')
    end
    
    r = grid_dist / grid_n;
    sigma_min = r / 2;
    
    % Optimal parameter space:  Sigma max is defined as the distance from mu to
    % P and is large enough to cover the mu space with 95 % probability. If not
    % using an sub-optimal definition of sigma max. 
    % Sigma_max = P_max - Mu_max
    
   if (grid_max - mu_max) >= (mu_dist / 2)
        sigma_max = (grid_max - mu_max)  / 2;
  
   else
       sigma_max = mu_dist / 4;
       warning("Distance between mu and P not large enough, setting sigma_max to span 95 % of mu.")
   end

    if sigma_min > sigma_max
        error('Sigma_min is smaller than sigma_max, try to increase grid resolution.')
    end
 
    prior_sigma = (mu_dist / 4) - eps; % numerical stability
    prior_sigma_latent =  norminv( (prior_sigma - sigma_min) ./ (sigma_max - sigma_min));
    
    if ~isfinite(prior_sigma_latent)
        warning('Prior over sigma is nonfinite, setting latent_prior to 0. ')
        prior_sigma_latent = 0;
    end
end
