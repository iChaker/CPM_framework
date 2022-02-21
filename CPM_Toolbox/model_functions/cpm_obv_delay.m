function [y,pE,pC] = cpm_obv_delay(varargin)
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
    y_true = M.cpm.other.current_signal;
    yn = y/max(abs(y));
    yn_true = y_true/max(abs(y_true));
    [corrs, lags] = xcorr(yn_true,yn);
    [maxc, argmax]= max(corrs);
    
    global delay; 
    delay = lags(argmax);
    try
    y = circshift(y,delay);
    catch
        
    end
    
end

function [pE,pC] = observation_priors()


pE.transit= 1 ;
pE.decay= 1 ;
pE.epsilon= 0.4584;

pC.transit= 1/128;
pC.decay= 1/128;
pC.epsilon= 0.24;

end


