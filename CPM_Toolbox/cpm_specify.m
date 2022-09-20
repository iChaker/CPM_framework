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

[~,~,~,fn_lmoment_priors,fn_lscaling_priors] =  feval(population_field);
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
    
    r = (grid_mu(2) - grid_mu(1)) / grid_mu(3);
    sigma_min = r / 2;
    
    if ~ (grid_mu(2) == model_mu(2))
        sigma_max = (grid_mu(2) - model_mu(2)) / 2;
    
    else

        warning("Possible parameter space for mu spans entire grid, setting sigma_max to half the grid space")
        sigma_max = (grid_mu(2) - grid_mu(1))  / 2;        
    end

    if sigma_min >= sigma_max
        warning('Sigma_min is smaller than sigma_max, try to increase grid resolution. Reversing order')
        [sigma_min, sigma_max] = deal(sigma_max, sigma_min);
    end
 
    options.sigma.(p_names{pidx}) = [sigma_min, sigma_max];
    
    prior_sigma = ( ((model_mu(2) - model_mu(1)) / 2) * 0.975) / 2;
    
    prior_sigma_latent =  norminv( (prior_sigma - sigma_min) ./ (sigma_max - sigma_min));
    
     if ~isfinite(prior_sigma_latent)
         warning('Setting latent prior over sigma to 0, the grid has to be twice as large to cover a 95 % prior over mu. Trying prior to cover 95 % of  grid.')
         prior_sigma = ( ((grid_mu(2) - grid_mu(1)) / 2) * 0.975) / 2;
         prior_sigma_latent = norminv( (prior_sigma - sigma_min) ./ (sigma_max - sigma_min));
         
        if ~isfinite(prior_sigma_latent)
            warning('Still not possible find to finite prior over sigma, setting latent_prior to 0. ')
            prior_sigma_latent = 0;
        end
 
     end
     
     options.lsigma_pE.(p_names{pidx}) = prior_sigma_latent;
     
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

