function output = cpm_grid_RPE(freeparams,fixedparams,data)
% Reward prediction error model

    % load args
    
    
    try eps = freeparams.eps; catch eps=0; end
    try eta = freeparams.eta; catch eta=0; end
    
    gamma = fixedparams.gamma;
    Lambda = fixedparams.Lambda;
    
    C = data.C;
    S = data.S;
    
    % define model
    try
                alpha = freeparams.alpha;
    catch
        try
                tau = freeparams.tau;
                alpha = [ logit_inv(tau - eps), logit_inv(tau + eps)];
        catch
            alpha = 0;
         end
    end

    
    
    RPE = cpm_TD_learning_RPE_generation(alpha ,eta, gamma, Lambda, S, C);

    RPEs = zeros( [length(RPE)*3 , 1]);
    RPEs(1:3:end) = RPE(:, 2);
    RPEs(2:3:end) = RPE(:, 3);
    RPEs(3:3:end) = RPE(:, 4);

    % If necessary rescale RPE to -1 and 1 - This can easily be
    % optimizied to be run after the for-loop. Keeping it for
    % readability.
    abs_max = max(abs(RPEs(:)));
    RPEs(:) = RPEs(:) ./ abs_max;  

    
    
    output=RPEs;
    
    
end

function out = logit_inv(a)
%LOGIT_INV Summary calculates the inverse logit of a
out = exp(a) ./ (1 + exp(a));
end



