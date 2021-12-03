function [ W ,inverse_latent_transform, latent_transform, pE,pC] = cpm_RF_template(varargin)
%cpm_RF_template receptive field shape and latent parameter definitions
%   This function bundles together all user definitnios related to the
%   receptive field. When no RF function is included PRF specification, the RF defaults 
%   to a Gaussian. See cpm_RF_Gaussian for more detail.
%
%   The user should define:
%
%       get_response(P,M,U): This function returns the weights of each
%                            point in the meshgrid according to the parameters P.
%
%       neural_priors(p_names): this function takes in the list of grid
%                               parameter names and outputs the relevant latent 
%                               parameter priors.
%
%       latent_transform(P,U): transforms original parameters back to
%                              their latent form.
%
%       inverse_latent_transform(P,U): transforms latent parameters back to
%                                      their original form.
%                                      NOTE: all P arguments are assumed to
%                                            be latent parameters
%       

%% ignore
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



function [pE,pC]= neural_priors(p_names)
%this function takes in the list of grid parameter names and outputs the 
% relevant latent parameter priors which will be fed into VL
%
%   INPUT:
%       p_names: cell contaning strings of the Theta parameter names
%   
%   OUTPUT:
%       pE: prior expectations of the parameters to be estimated by VL
%       pC: prior certainty 
%
% example for a simple Gaussian prf: p_names= { 'Theta_1' 'Theta_2' ... 'Theta_N'}
% priors would be: pE.lmu_Theta_1 = 0      pC.lmu_Theta_1 = 1       
%                  pE.lmu_Theta_2 = 0      pC.lmu_Theta_2 = 1
%                      ...                     ...           
%                  pE.lmu_Theta_2 = 0      pC.lmu_Theta_N = 1
%                  pE.lbeta = -2           pC.lbeta = 5
% 
% use add_suffix(p_names,'lsomethin_') to add a suffix to all Thetas in one
% line

% define your priors here









pE.lbeta    = -2;    pC.lbeta    = 5;


end



function x = get_response(P,M,U)
% This function returns the weights of each point in the meshgrid according 
% to the latent VL parameters P.

P = inverse_latent_transform(P,U); % important step

p_names = fieldnames(U(1).grid); % useful cell contaning Theta names
coords = U(1).gridpoints; % grid coordonates are automatically saved in U(1).gridpoints


% define your receptive field here







% x = some distribution or function  






% Normalize and Scale with free parameter beta
xx = 1/sum(x);
x = xx*P.beta .* x;



end



function correct_P = inverse_latent_transform(P,U)
% Function to transforms latent VL parameters back to their original form.
% Note: there is no need to transform latent parameters of the observation model.

latent_parameter_names =  fieldnames(P);
correct_P = struct();
grid = U(1).grid; % meshgrid structure contaning min max and stepN for each Theta

for i=1:length(latent_parameter_names)
    lparam_name = latent_parameter_names{i};
    
    % your code here

        %if startsWith(lparam_name,'lmu_')

        %elseif startsWith(lparam_name,'lsigma_')

        %elseif isequal(lparam_name,'lbeta') 

        %end
    
end

end

function P = latent_transform(P,U)

    %only useful if you want to make simulations, no need to implement
    P;



end

function RF_names = add_suffix(p_names,suffix)

    RF_names = cellfun(@(x) [suffix x],p_names,'UniformOutput',false);
    
end
