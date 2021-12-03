function output = cpm_grid_BOLD(freeparams,fixedparams,data)


    % load args

    tau = freeparams.tau;
    try eps = freeparams.eps; catch eps=0; end
    try eta = freeparams.eta; catch eta=0; end
    
    
    gamma = fixedparams.gamma;
    Lambda = fixedparams.Lambda;
    TR= fixedparams.TR;
    nscans = fixedparams.nscans;
    
    C = data.C;
    S = data.S;
    
    % define model
    
    alpha = [ logit_inv(tau - eps), logit_inv(tau + eps)];

    
    
    RPE = cpm_TD_learning_RPE_generation(alpha ,eta, gamma, Lambda, S, C);

    RPEs = zeros( [length(RPE)*3 , 1]);
    RPEs(1:3:end) = RPE(:, 2);
    RPEs(2:3:end) = RPE(:, 3);
    RPEs(3:3:end) = RPE(:, 4);

    % If necessary rescale RPE to -1 and 1 - This can easily be
    % optimizied to be run after the for-loop. Keeping it for
    % readability.
    if max(abs(RPEs(:))) > 1
        abs_max = max(abs(RPEs(:)));
       RPEs(:) = RPEs(:) ./ abs_max;  
    end
    
    
    try P.transit = freeparams.transit; catch P.transit=0; end
    try P.decay = freeparams.decay; catch P.decay=0; end
    try P.epsilon = freeparams.epsilon; catch P.epsilon=0.07; end
    
    M.x = [0;0;0;0];
    M.f = 'spm_prf_fx';
    M.g = 'spm_prf_gx';
    M.n=4;
    M.u=1;
    M.l=1; 
    M.m=1;
    M.u=0;
    
    % Convert input timing from seconds -> microtime bins
    bins_per_second = 1 / data.dt(1);
    bins_per_TR     = bins_per_second * TR;
    nbins = nscans*bins_per_TR;
    z = zeros(1,nbins);
    
            
    for t = 1:length(RPEs)   
        % Microtime index for this volume
        start_bin = ceil( data.ons(t) / data.dt(t)) + 1;
        end_bin   = start_bin + (data.dur(t) * bins_per_second) - 1;
        ind = start_bin : end_bin;  
        
        z(ind) = z(ind) + RPEs(t);
        
    end
    Z.u=z';
    Z.dt=data.dt(1);

    output=spm_int_sparse(P,M,Z);
    
    
end

function out = logit_inv(a)
%LOGIT_INV Summary calculates the inverse logit of a
out = exp(a) ./ (1 + exp(a));
end

