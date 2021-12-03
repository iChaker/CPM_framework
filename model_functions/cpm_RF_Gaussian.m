function [ W ,inverse_latent_transform, latent_transform, pE,pC] = cpm_RF_Gaussian(varargin)
%CPM_RECEPTIVE_FIELD Summary of this function goes here
%   Detailed explanation goes here

if nargin ==3 
    
    P = varargin{1};
    M = varargin{2};
    U = varargin{3};
    
    W = get_response(P,M,U);
    latent_transform=@latent;
    inverse_latent_transform=@inverse_latent_transform;
    pE=[];
    pC=[];
else
    p_names=varargin{1};
    W = [];
    latent_transform=[];
    inverse_latent_transform=[];
    [pE,pC]=neural_priors(p_names);
end

end

function x = get_response(P,M,U)


% Get PRF response to a stimulus 
%
% P  - parameters
% coords - coordinates (tau,eps,eta)
P = inverse_latent_transform(P,U);
p_names = fieldnames(U(1).grid);
E_names = add_suffix(p_names,'mu_');
C_names = add_suffix(p_names,'sigma_');

coords = U(1).gridpoints;
%coords = flip(coords,2); % flip??
V = [];
Mu= [];
for i=1:size(p_names,1)
   Mu(i) = P.(E_names{i});
   V(i) = P.(C_names{i}).^2;
    
end

V = diag( V );

% normalized multivariate gaussian
x = spm_mvNpdf(coords', Mu, V)*sqrt(((2*pi)^3)*det(V));

xx = 1/sum(x);

% Scale
x = xx*P.beta .* x;



end


function [pE,pC]= neural_priors(p_names)

latent_mus = add_suffix(p_names,'lmu_');
latent_sigmas = add_suffix(p_names,'lsigma_');
neural_parameter_names = { latent_mus{:} latent_sigmas{:}};

E_values = zeros(length(p_names)*2,1);
E_values = num2cell(E_values,2);

C_values = ones(length(p_names)*2,1);
C_values = num2cell(C_values,2);

pE = cell2struct(E_values,neural_parameter_names,1);
pC = cell2struct(C_values,neural_parameter_names,1);

pE.lbeta    = -2;    pC.lbeta    = 5;


end


function correct_P = inverse_latent_transform(P,U)
% Prepare neuronal parameters for use e.g. exponentiate log parameters

% P - parameters
% M - model

% change this variables to change the "exploration zone" of the parameter
% space
SIGMA_MIN = 0.1;
SIGMA_MAX = 1;

latent_parameter_names =  fieldnames(P);
correct_P = struct();
grid = U(1).grid;
for i=1:length(latent_parameter_names)
    lparam_name = latent_parameter_names{i};
    if startsWith(lparam_name,'lmu_')
        param_name= (lparam_name(2:end));
        theta_name= (lparam_name(5:end));
        maxp = grid.(theta_name)(2);
        minp =grid.(theta_name)(1);
        
        correct_P.(param_name)= ((maxp-minp) .* spm_Ncdf(P.(lparam_name),0,1)) + minp ;
    
    elseif startsWith(lparam_name,'lsigma_')
        param_name= (lparam_name(2:end));
         correct_P.(param_name)= ((SIGMA_MAX-SIGMA_MIN) .* spm_Ncdf(P.(lparam_name),0,1)) + SIGMA_MIN ;
    
    elseif isequal(lparam_name,'lbeta')
        correct_P.beta = exp(P.lbeta);

        
        
    end
    
end


end

function P = latent_transform(P,U)

    %only useful if you want to make simulations, no need to implement
    P;



end

function RF_names = add_suffix(p_names,suffix)

    RF_names = cellfun(@(x) [suffix x],p_names,'UniformOutput',false);
    
end
