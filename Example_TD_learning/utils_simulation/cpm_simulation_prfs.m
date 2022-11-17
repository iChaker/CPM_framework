function [PRFname, PRF] = cpm_simulation_prfs(VOI, cfg_file, REDO)

arguments
    VOI
    cfg_file = 'one_tau_recovery.json'
    REDO = true
end


cfg = spm_jsonread(cfg_file);

tmpdir = cfg.tmpdir;
simulationdir = cfg.simulationdir;
recoveryname = cfg.recoveryname;

mkdir(tmpdir); % Create a directory for temporary files which will be ignored.
mkdir(simulationdir);
%% ======== SIMULATION SETUP =======================


% We can freely chose for which TR we want to simulate data, we use the TR from
% a previous project of ours. 
TR = 0.592;
% The file containing the experimental inputs.
file="tsvfile1.tsv";
% dt is the micro resolution, which is often 8 or 20, but for example when slice-time
% correction has been produced should be set to the number of slices.
dt = TR / 28; % Assuming classical SPM approach 

grid = {};
% These are the options used by BayesPRF:
options = {};
options.TE=0.03;      % Echo time (you may add other options from BayespRF)
options.model='spm_prf_response'; % We use a modified version of this function.
options.name =recoveryname;


for pn = cfg.param_names(:)'
    grid.(pn{1}) = [cfg.(['grid_' pn{1}])(:)', cfg.resolution];

    try
        options.mu.(pn{1}) = cfg.(['mu_' pn{1}])(:)';
    catch
        sprintf("Mu boundaries not defined for recovery. %s\n", pn)
    end
end

%% Select grid model & precompute
model = "cpm_grid_RPE"; % user defined model, see cpm_grid_template for details

fixedparams.gamma=0.97;  % these fieldnames must be the same fieldnames used in the cpm_grid
fixedparams.lambda=1.0;

% Load adata
data = [];
[trial_n, ~, ~, S, ~, C, onsets] = cpm_prepare_experiment_data(file);
data.dur = ones(3 * trial_n, 1) .* 2 .* dt;
data.dt = ones(3 * trial_n, 1) .* dt;
data.ons = onsets;
data.C = C;
data.S = S;

% The BayesPRF requires an SPM struct, but only a few fields from there, which
% we generate here:
SPM = {};
SPM.xY.RT = TR;
SPM.swd = ''; % this field will be overwritten by CPM.

outpath=simulationdir;  % PRF file path
% Now we generate / precompute the recovery grid:

U_recovery = cpm_precompute(model, grid, fixedparams, ...
                                                    data, fullfile(simulationdir, [recoveryname '_U_recovery']), REDO);

% And specify the PRF for recovery:
PRF = cpm_specify(SPM, options, VOI.xY.y, VOI.xY.XYZmm, ...
                   U_recovery, 'cpm_RF_Gaussian', 'cpm_obv_int', outpath);

PRF.M.noprint = 0; % to suppress spm outputs
PRFname = fullfile(simulationdir, [recoveryname '_PRFn.mat']);

save(PRFname, 'PRF'),