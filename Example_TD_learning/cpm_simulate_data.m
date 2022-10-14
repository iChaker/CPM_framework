function [PRFname, PRF] = cpm_simulate_data(cfg_file, REDO, use_par)

arguments
    cfg_file = 'ideal_grid_cfg.json'
    REDO = true
    use_par =true
end


%%  ===========================================
% =================== SETUP ==================== =
% =============================================
% Random seed for reproducibility of noise etc.
rng(2022)
%% Creating directories to store simulation files 
cfg = spm_jsonread(cfg_file);

tmpdir = cfg.tmpdir;
simulationdir = cfg.simulationdir;
simulationname = cfg.simulationname;

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
% As this is a simulation, nscans are inferred from the data, but you can change
% it here manually:
nscans = nan;
% I see two possibilities in generating the data for our recovery study
% 1. Using a grid and moving the parameters to certain areas within the grid,
% where the grid is similar to the one we use to recover the values.
% 2. Using the point estimates (from U) and generating the BOLD signal directly
% from the RPE trajectory.


nnoise = length(cfg.sim_noise);
nnames = length(cfg.param_names);

p_grid = {};
rec_grid = {};
% These are hte options used by BayesPRF:
options = {};
options.name='simulation'; % PRF file tag
options.TE=0.03;      % Echo time (you may add other options from BayespRF)
options.model='spm_prf_response'; % We use a modified version of this function.
%
rec_options = options;
rec_options.name = 'recovery';

mu_simulations = {};
sigma_simulations ={};
point_simulations = {};

pc = 1;
for pn = cfg.param_names(:)'
    p_grid.(pn{1}) = [cfg.(['p_grid_' pn{1}])(:)', cfg.p_resolution];
    rec_grid.(pn{1}) = [cfg.(['rec_grid_' pn{1}])(:)', cfg.rec_resolution];

    try 
        options.mu.(pn{1}) = cfg.(['p_mu_' pn{1}])(:)';
    catch
        disp("Mu boundaries not defined for simulation.")
    end
    
    try
        rec_options.mu.(pn{1}) = cfg.(['rec_mu_' pn{1}])(:)';
    catch
        disp("Mu boundaries not defined for recovery.")
    end
    mu_simulations{end + 1, 1} = cfg.(['sim_mu_' pn{1}])(:);
    sigma_simulations{end + 1, 1} = cfg.(['sim_sigma_' pn{1}])(:);
    point_simulations{end + 1, 1} = cfg.(['point_mu_' cfg.point_names{pc}])(:);
    pc = pc + 1;
end

pointParams = {};
point_mus = cell(length(point_simulations), 1);
[point_mus{:}] = ndgrid(point_simulations{:});

nvoxels = length(point_mus);

pc = 1;
for  pn = cfg.point_names(:)'
    tmp_mu = point_mus{pc};
    for ii = 1 : length(tmp_mu(:))
        pointParams{ii}.([pn{1}]) = tmp_mu(ii);
    end
    pc = pc + 1;
end

dimParams = {};
for  pn = cfg.param_names(:)'
    [tmp_mu, tmp_sig] = ndgrid(cfg.(['sim_mu_' pn{1}])(:), cfg.(['sim_sigma_' pn{1}])(:));
    nvoxels = nvoxels + length(tmp_mu(:));
    tmp_mu = tmp_mu(:);
    tmp_sig = tmp_sig(:);
    for ii = 1 : length(tmp_mu)
        dimParams.(pn{1}){ii}.(['mu_' pn{1}]) = tmp_mu(ii);
        dimParams.(pn{1}){ii}.(['sigma_' pn{1}]) = tmp_sig(ii);
    end
end

% The resolution of the different grids, generation grid will only be used when
mu_grid = cell(nnames, 1);
sigma_grid  = cell(nnames, 1);
[mu_grid{1 : nnames}, sigma_grid{1:  nnames}] = ndgrid(mu_simulations{:}, sigma_simulations{:});

nvoxels = nvoxels  + length(mu_grid{1}(:));

fullParams = {};
for np = 1 : nnames
    tmp_mu = mu_grid{np}(:);
    tmp_sig = sigma_grid{np}(:);
        for ii = 1 : length(tmp_mu(:))
            fullParams{ii}.(['mu_' cfg.param_names{np}]) = tmp_mu(ii);
            fullParams{ii}.(['sigma_' cfg.param_names{np}]) = tmp_sig(ii);
        end
end

%% =============== SIMULATING DATA ===============
% Preparing behavioural data for the simulation:
% Implementation note here: cpm_precompute requires data to only have the fields
% dur, dt, and ons. Everything else is further passed to the generation
% functions. In this case S and C which contain the complete serial compound
% vector and the wealth trajectory.
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

% As we are using a simulation we can infer the number of scans from the onsets
% and the TR:
if isnan(nscans)
    nscans = ceil(max(onsets) / TR);
end

%% Select grid model & precompute
model = "cpm_grid_RPE"; % user defined model, see cpm_grid_template for details

fixedparams.gamma=0.97;  % these fieldnames must be the same fieldnames used in the cpm_grid
fixedparams.Lambda=1.0;

% The number of voxels is hte number of combinations and the number of noise
% levels:

y = zeros(nscans, nvoxels * nnoise);
xyz = zeros(3, nvoxels * nnoise); 
xyz_def = cell(nvoxels, 1);

nv = 1;
mv = 1;

for pidx = 1 : length(pointParams)
    y(:, nv) =  cpm_generative_point_process(pointParams{pidx}, model, fixedparams, data, TR, nscans, tmpdir, simulationname);
    xyz(1, nv) = mv;
    xyz_def{nv} = pointParams{pidx};
    nv = nv + 1; 
end


for pn = cfg.param_names(:)'
    mv = mv + 1;
    part_p_grid = {};
    part_options = options;
    part_options.mu = {};
    part_options.mu.(pn{1}) = options.mu.(pn{1});
    part_p_grid.(pn{1}) = p_grid.(pn{1});
    U_part = cpm_precompute(model, part_p_grid, fixedparams, data, fullfile(tmpdir, [simulationname '_simU_' pn{1}  '_grid.mat']), true);
    for pidx = 1 : length(dimParams.(pn{1}))
        y(:, nv) = cpm_generative_grid_process(dimParams.(pn{1}){pidx}, SPM, U_part, nscans, part_options, tmpdir);
        xyz(1, nv) = mv;
        xyz_def{nv} = dimParams.(pn{1}){pidx};
        nv = nv + 1;
    end
end

U_full = cpm_precompute(model, p_grid, fixedparams, data, fullfile(tmpdir, [simulationname '_simU_full_grid.mat']), true);

mv = mv + 1;
for pidx = 1 : length(fullParams)
    y(:, nv) = cpm_generative_grid_process(fullParams{pidx}, SPM, U_full, nscans, options, tmpdir);
    xyz(1, nv) = mv;
    xyz_def{nv} =fullParams{pidx};
    nv = nv + 1;
end


for nidx = 1 : nnoise
    
    for vidx =  1 : nvoxels
        new_idx = vidx + nvoxels * (nidx - 1);
    % normrand with 0 variance returns the distributions mean:
        y(:, new_idx) = y(:, vidx) +  normrnd(0, cfg.sim_noise(nidx), [nscans,1]);
        xyz(1,  new_idx) = xyz(1, vidx);
        xyz(2,  new_idx) = nidx;
    end
end

%% 
% Storing and saving the VOI in SPM format:
VOI.xY.y = y;
VOI.xY.XYZmm = xyz;
VOI.xyz_def = xyz_def;
% Mean values for completeness
VOI.y = mean(y, 2); % for completeness
save(fullfile(simulationdir, [simulationname, '_simVOI.mat']), 'VOI')


%% ===================== MODEL INVERSION ======================
outpath=simulationdir;  % PRF file path
% Now we generate / precompute the recovery grid:

U_recovery = cpm_precompute(model, rec_grid, fixedparams, data, fullfile(simulationdir, [simulationname '_U_recovery']), REDO);

% And specify the PRF for recovery:
PRF = cpm_specify(SPM, options, VOI.xY.y, VOI.xY.XYZmm, ...
                   U_recovery, 'cpm_RF_Gaussian', 'cpm_obv_int', outpath);

PRF.M.noprint = 0; % to suppress spm outputs
PRFname = fullfile(simulationdir, [simulationname '_PRFn.mat']);

save(PRFname, 'PRF'),

end