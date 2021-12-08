function [ fn_RF_weights, fn_true_moments,fn_true_scaling, fn_latent_moment_priors,fn_scaling_priors,fn_latent_moments,fn_latent_scaling  ] = cpm_RF_template(varargin)
%CPM_RECEPTIVE_FIELD Summary of this function goes here
%   Detailed explanation goes here
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
% return weight vector W (1 x npoints) for the coords
%
% coords (nparams x npoints) is the list of coordinates in the parameter grid
%
% true_parameters's fieldnames are defined in get_true_moments and
% get_true_scalo,g

mu = true_parameters.mu ; % contains vector [ mu_param1 ... mu_paramN ]
beta = true_parameters.beta ; % contains scalar


W = [];
end
%% Priors
function [pE,pC] = get_latent_moment_priors()
% define priors for the receptive field moments
% chosen fieldnames will indicate how latent parameters are suffixed
    pE.lmu = 0;
    pC.lmu = 1;

end

function [pE,pC] = get_latent_scaling_priors()
% define priors for the scaling priors (such as beta)
% chosen fieldnames will be the names of latent scaling factors 
% you may also use this to define other scalar parameters for the neuronal
% model z
    pE.lbeta=0;
    pC.lbeta=5;
end

%% true transforms
function true_moments = get_true_moments(latent_moments, minp,maxp)
% transforms latent parameters into true parameters
%
% latent_moments's fieldnames are defined in get_latent_moment_priors
% To avoid confusion it's recommended to use different fieldnames for true
% and latent moments
%
%
% minp and maxp are the bounds set by the initial grid specification 

lmu = latent_moments.lmu;

mu = ((maxp-minp) .* spm_Ncdf(lmu,0,1)) + minp ;

true_moments.mu=mu;
end

function true_scaling = get_true_scaling(latent_scaling)
    % latent_scaling's fieldnames are defined in get_latent_scaling_priors

    lbeta = latent_scaling.lbeta;
    beta = exp(lbeta);
    true_scaling.beta = beta;
end

%% latent transforms

function latent_scaling = get_latent_scaling(true_scaling)

    % not required, but useful
    latent_scaling=[];
    
end

function latent_moments = get_latent_moments(true_moments,minp,maxp)
    
    % not required, but useful
    latent_moments=[];
    
end
