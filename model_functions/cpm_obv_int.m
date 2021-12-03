function [y,pE,pC] = cpm_obv_int(varargin)
%CPM_OBV_IDENTITY Summary of this function goes here
%   Detailed explanation goes here
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


function [y] = observation_function(P,M,U)

    y=spm_int_sparse(P,M,U);
end

function [pE,pC] = observation_priors()


pE.transit= 1 ;
pE.decay= 1 ;
pE.epsilon= 0.4584;

pC.transit= 1/128;
pC.decay= 1/128;
pC.epsilon= 0.24;

end


