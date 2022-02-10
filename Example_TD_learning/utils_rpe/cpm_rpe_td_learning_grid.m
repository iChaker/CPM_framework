function U = cpm_rpe_td_learning_grid(file, options)
% precompute TD learning parameter
%
%% inputs:
% file: tsv file containing event timings 
%
% options.tau_min
% options.tau_max  
% options.tau_stepN
% 
% options.eps_min
% options.eps_max
% options.eps_stepN
%                        
% options.eta_min
% options.eta_max
% options.eta_stepN

% options.TR                 (SPM.xY.RT)
% options.nmicrotime         (SPM.xBF.T)

% options.Lambda;            [optional]  (defaults to .99)
% options.gamma;             [optional]  (defaults to 1) 
%% outputs
% U structure contaning mesh grid for each event:
% U(i).RPE : meshgrid
% U(i).ons
% U(i).dur
% U(i).dt

%%
% Unpacking option values here, so that possible defaults can be added. 
TR = options.TR;
nmicrotime = options.nmicrotime;
dt = TR / nmicrotime;

try
    gamma = options.gamma;
catch
    disp('Setting discount rate (gamma) to 1')
    gamma = 1;
    options.gamma = 1;
end

try 
    Lambda = options.Lambda;
catch
    disp('Setting decay rate (Lambda) to 0.99')
    Lambda = 0.99;
    options.Lambda= 0.99;
end

% display parameters
disp('Calculating RPE with grid parameters: ');
disp(options);

% Tau values
tau_min = options.tau_min;
tau_max = options.tau_max;
tau_stepN = options.tau_stepN;
tau_values = linspace(tau_min, tau_max, tau_stepN);
% eps
try
    eps_min = options.eps_min;
    eps_max = options.eps_max;
    eps_stepN = options.eps_stepN;
    eps_values = linspace(eps_min, eps_max, eps_stepN);
catch
    disp('No eps defined, assuming 0')
    eps_stepN = 1;
    eps_values = [0];
end
% eta

try
    eta_min = options.eta_min;
    eta_max = options.eta_max;
    eta_stepN = options.eta_stepN;
    eta_values = linspace(eta_min, eta_max, eta_stepN);
catch
    disp('No eta defined, assuming 0')
    eta_stepN = 1;
    eta_values = [0];
end    
 

[trial_n, timestep_n, stimuli_n, S, stimuli, C, onsets] = cpm_prepare_experiment_data(file);

% Building up U-struct

durations = ones(3 * trial_n, 1) .* 2 .* dt;
dts = ones(3 * trial_n, 1) .* dt;

RPEs = zeros(3 * trial_n, tau_stepN, eps_stepN, eta_stepN);



for tau_i = 1 : tau_stepN
    for eps_i = 1 : eps_stepN
        for eta_i = 1 : eta_stepN
            
            alpha = [cpm_logit_inv(tau_values(tau_i) - eps_values(eps_i)), ...
                     cpm_logit_inv(tau_values(tau_i) + eps_values(eps_i))];
            
            RPE = cpm_TD_learning_RPE_generation(alpha ,eta_values(eta_i), gamma, Lambda, S, C);
            
            RPEs(1:3:end, tau_i, eps_i, eta_i) = RPE(:, 2);
            RPEs(2:3:end, tau_i, eps_i, eta_i) = RPE(:, 3);
            RPEs(3:3:end, tau_i, eps_i, eta_i) = RPE(:, 4);
            
            % If necessary rescale RPE to -1 and 1 - This can easily be
            % optimizied to be run after the for-loop. Keeping it for
            % readability.
            if max(abs(RPEs(:, tau_i, eps_i, eta_i))) > 1
                abs_max = max(abs(RPEs(:, tau_i, eps_i, eta_i)));
               RPEs(:, tau_i, eps_i, eta_i) = RPEs(:, tau_i, eps_i, eta_i) ./ abs_max;  
            end
        end
    end
end

% Packing everything into U struct

U(3 * trial_n) = struct();

for nt = 1 : 3 * trial_n
    U(nt).RPE = squeeze(RPEs(nt, :, :, :));
    U(nt).ons = onsets(nt);
    U(nt).dur = durations(nt);
    U(nt).dt = dts(nt); 
end



end