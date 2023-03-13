function [RPEs, Xs, Ws] = cpm_TD_learning_RPE_generation(alpha, eta, gamma, Lambda, S, C)

% TODO documentation, testing

if length(alpha) == 2
    alpha_neg = alpha(1);
    alpha_pos = alpha(2);
else
    alpha_neg = alpha;
    alpha_pos = alpha;
end

% Recovering run parameters from stimuli.
[trial_n, timestep_n, stimuli_n] = size(S);

% Weightings for each timestep, using CSC
Xs = zeros(trial_n, timestep_n, timestep_n * stimuli_n);

% When there is Payout
primary_reward_timestep = timestep_n;


%% loop over trials
Xs = zeros(trial_n, timestep_n, timestep_n * stimuli_n);

for l =  1 : trial_n
    X_temp = [];
    
    for s = 1 : stimuli_n
       
        stimulus_scr = zeros(timestep_n, timestep_n);
        
        
        
        for t = 1 : timestep_n
            current_stim = S(l, :, s);
            current_stim = current_stim(1 : t);
            current_stim = flip(current_stim);
            stimulus_scr(t, 1 : t) = current_stim;
        end
        
        stimulus_scr(primary_reward_timestep:end, :) = 0;
        
        X_temp = cat(2, X_temp, stimulus_scr);
    end

    X_temp = X_temp ./ max(1, max(sum(X_temp, 2)));
    
    Xs(l, :, :) = X_temp;
    
end

Vs = zeros(trial_n, timestep_n);
W = zeros(timestep_n .* stimuli_n, 1);
RPEs = zeros(trial_n, timestep_n);

Ws = zeros(timestep_n .* stimuli_n, timestep_n, trial_n);

for l = 1 : trial_n
   
    Z = zeros(timestep_n, timestep_n .* stimuli_n);
    
    for t = 1 : timestep_n - 1
       
        Z(t + 1, :) = Lambda .* gamma .* Z(t, :) + squeeze(Xs(l, t + 1, :))';
        c_delta = cpm_isoelastic_utility(C(l, t + 1), eta) - cpm_isoelastic_utility(C(l, t), eta);
        w_delta = gamma .* squeeze(Xs(l, t + 1, :))' * W - squeeze(Xs(l, t, :))' * W;
        
        RPEs(l, t + 1) = c_delta + w_delta;
        
        if RPEs(l, t + 1) >= 0
            W = W + (alpha_pos .* RPEs(l, t + 1) .* Z(t, :))';
        else
            W = W + (alpha_neg .* RPEs(l, t + 1) .* Z(t, :))';
        end
        
        Ws(:, t, l) = W;
    end
    
    Vs(l, :) = squeeze(Xs(l, :, :)) * W;
    
    if l < trial_n - 1
        RPEs(l + 1) = RPEs(l);
    end
    
end
end