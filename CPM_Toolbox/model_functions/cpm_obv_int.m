function [y,pE,pC] = cpm_obv_int(varargin)
%cpm_obv_int Extended Balloon model as observation function
% see spm_prf_fx and spm_prf_fx in BayespRF toolbox for more details
if nargin ==3
    
    P = varargin{1};
    M = varargin{2};
    Z = varargin{3};
    
    y = observation_function(P,M,Z);
else
    y = [];
    [pE,pC] = observation_priors() ;
    
end

end


function [y] = observation_function(P,M,Z)
    % These are the default settings for M, uncomment and change to customize hemodynamics:
    %    M = struct()
    %    M.x = [0;0;0;0];
    %    M.f = 'spm_prf_fx';
    %    M.g = 'spm_prf_gx';
    %    M.n=4;
    %    M.u=1;
    %    M.l=1; 
    %    M.m=1;
    %    M.u=0;
    %    M.ns = ...;
    
    y=spm_int_sparse(P,M,Z);
end

function [pE,pC] = observation_priors()


pE.transit= 0 ;
pE.decay= 0 ;
pE.epsilon= 0.4584;

pC.transit= 1/128;
pC.decay= 1/128;
pC.epsilon= 0.24;

end


