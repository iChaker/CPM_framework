function [ fn_RF_weights, fn_true_moments,fn_true_scaling, fn_latent_moment_priors,fn_scaling_priors,fn_latent_moments,fn_latent_scaling  ] = cpm_RF_Gaussian(varargin)
%cpm_RF_Gaussian One Gaussians as population field
fn_RF_weights=@get_RF_weights;
fn_true_moments=@get_true_moments;
fn_latent_moment_priors=@get_latent_moment_priors;
fn_scaling_priors =@get_latent_scaling_priors;
fn_true_scaling=@get_true_scaling;
fn_latent_scaling=@get_latent_scaling; %optional
fn_latent_moments=@get_latent_moments; %optional
end

%% Receptive field definition
function W = get_RF_weights(coords,true_parameters)

mu = true_parameters.mu;
sigma = true_parameters.sigma;
beta = true_parameters.beta;

sigma = sigma .* sigma;

V = diag(sigma);
x = spm_mvNpdf(coords, mu, V) ; %* sqrt(((2*pi)^3)*det(V));

%xx = 1 / sum(x) * sqrt(((2*pi)^3)*det(V));

x = x ./ sum(x + eps);

if ~isfinite(sum(x))
    warnstring =[];
    for sig = sigma
        warnstring = [warnstring sprintf(' %4.5f', sig)];
    end    
    warning([ 'spm_mvNpdf sums to 0 with sigmas =', warnstring])
end

% Scale
W = beta .* x;
end
%% Priors
function [pE,pC] = get_latent_moment_priors()
% latent prior distributions for receptive field moments 
    pE.lmu = 0;
    pC.lmu = 1.0;
    
    pE.lsigma = -2.5;
    pC.lsigma = 1.0;
    
end

function [pE,pC] = get_latent_scaling_priors()
% latent priors for beta 
    pE.lbeta=0;
    pC.lbeta=5;
end

%% true transforms
function true_moments = get_true_moments(latent_moments, minp, maxp, mins, maxs)
% 

lmu = latent_moments.lmu;
lsigma = latent_moments.lsigma;

% THESE SHOULD BE HYPERPARAMETERS
SIGMA_MIN = mins;
SIGMA_MAX = maxs;

mu = ((maxp-minp) .* spm_Ncdf(lmu,0,1)) + minp ;
sigma= ((SIGMA_MAX-SIGMA_MIN) .* spm_Ncdf(lsigma,0,1)) + SIGMA_MIN ;

true_moments.mu=mu;
true_moments.sigma=sigma;
end

function true_scaling = get_true_scaling(latent_scaling)
    % true_scaling structure must have the same order as pE and pC defined
    % in get_latent_scaling_priors

    lbeta = latent_scaling.lbeta;
    beta = exp(lbeta);
    true_scaling.beta = beta;
end

%% latent transforms

function latent_scaling = get_latent_scaling(true_scaling)

    % not required, but useful
    latent_scaling=[];
    
end

function latent_moments = get_latent_moments(true_moments,grid)
    
    % not required, but useful
    latent_moments=[];
    
end
