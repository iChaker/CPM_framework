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
TR = 2.55;
% The file containing the experimental inputs.
% dt is the micro resolution, which is often 8 or 20, but for example when slice-time
% correction has been produced should be set to the number of slices.
dt = TR / 22; % Assuming classical SPM approach 

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
model = "cpm_grid_SAMSRF"; % user defined model, see cpm_grid_template for details

fixedparams.resolution = cfg.resolution;
% Load example data:
appertures = load('example_data/aps_Bars.mat');
appertures = appertures.ApFrm;

appertures_post = zeros(size(appertures, 3), cfg.resolution, cfg.resolution);

% We will now resize and pad the appertures, so they can be better handled.
% Because we use the default CPM model, we have to draw a inner box, which
% limits the location. We further assume that the appertures a slightly larger
% than the inner field.
for ii = 1  : size(appertures, 3)
        tmp_image = squeeze(appertures(:, :, ii));
        tmp_image = imresize(tmp_image, [cfg.resolution - 16, cfg.resolution - 16]);
        tmp_image = padarray(tmp_image, [8, 8], 0, 'both');
        appertures_post(ii, :, :) = tmp_image;
end
%%
trial_n = size(appertures, 3);
data.dur = ones(3 * trial_n, 1) .* TR;  % Duration of stimuli is assumed to be as long as each TR
data.dt = ones(3 * trial_n, 1) .* dt; % dt is the same across
data.ons = [0 : size(appertures, 3) - 1] .* TR; % Duration is assumed to be as long as number of stimuli and onsets every TR.
data.appertures = appertures_post; 


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