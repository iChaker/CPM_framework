function [y,pE,pC] = cpm_obv_identity(varargin)
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


function [y] = observation_function(P,M,Z)
    n = int32(floor(length(Z.u)/M.ns));
    y = downsample(Z.u,n);
end

function [pE,pC] = observation_priors()

pE = struct();
pC = struct();

end



