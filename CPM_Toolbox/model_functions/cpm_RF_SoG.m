function [ fn_RF_weights, fn_true_moments,fn_true_scaling, fn_latent_moment_priors,fn_scaling_priors,fn_latent_moments,fn_latent_scaling  ] = cpm_RF_SoG(varargin)
%cpm_RF_SoG sum of two Gaussians as population field

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

mu1 = true_parameters.mu1;
sigma1 = true_parameters.sigma1;
mu2 = true_parameters.mu2;
sigma2 = true_parameters.sigma2;
beta1 = true_parameters.beta1;
beta2 = true_parameters.beta2;

V1= diag(sigma1);
x1 = spm_mvNpdf(coords, mu1, V1)*sqrt(((2*pi)^3)*det(V1));

V2= diag(sigma2);
x2 = spm_mvNpdf(coords, mu2, V2)*sqrt(((2*pi)^3)*det(V2));

xx1 = 1/sum(x1);
xx2 = 1/sum(x2);

% Scale
W = xx1*beta1 .* x1 + xx2*beta2 .* x2;
end
%% Priors
function [pE,pC] = get_latent_moment_priors()
% latent priors for receptive field moments 
    pE.lmu1 = 0.1;
    pC.lmu1 = 1;
    
    pE.lsigma1=0;
    pC.lsigma1=1;
    
    pE.lmu2 = -0.1;
    pC.lmu2 = 1;
    
    pE.lsigma2=0;
    pC.lsigma2=1;
    
end

function [pE,pC] = get_latent_scaling_priors()
% latent priors for beta 
    pE.lbeta1=0;
    pC.lbeta1=5;
    
    pE.lbeta2=0;
    pC.lbeta2=5;
end

%% true transforms
function true_moments = get_true_moments(latent_moments, minp,maxp)
% 

lmu1 = latent_moments.lmu1;
lsigma1 = latent_moments.lsigma1;
lmu2 = latent_moments.lmu2;
lsigma2 = latent_moments.lsigma2;

SIGMA_MIN = 0;
SIGMA_MAX = 1.5;

mu1 = ((maxp-minp) .* spm_Ncdf(lmu1,0,1)) + minp ;
sigma1= ((SIGMA_MAX-SIGMA_MIN) .* spm_Ncdf(lsigma1,0,1)) + SIGMA_MIN ;

true_moments.mu1=mu1;
true_moments.sigma1=sigma1;

mu2 = ((maxp-minp) .* spm_Ncdf(lmu2,0,1)) + minp ;
sigma2= ((SIGMA_MAX-SIGMA_MIN) .* spm_Ncdf(lsigma2,0,1)) + SIGMA_MIN ;

true_moments.mu2=mu2;
true_moments.sigma2=sigma2;
end

function true_scaling = get_true_scaling(latent_scaling)
    % true_scaling structure must have the same order as pE and pC defined
    % in get_latent_scaling_priors

    lbeta1 = latent_scaling.lbeta1;
    beta1 = exp(lbeta1);
    true_scaling.beta1 = beta1;
    
    lbeta2 = latent_scaling.lbeta2;
    beta2 = exp(lbeta2);
    true_scaling.beta2 = beta2;
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
