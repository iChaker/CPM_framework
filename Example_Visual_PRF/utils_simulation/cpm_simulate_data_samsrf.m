function [VOI, voiname] = cpm_simulate_data_samsrf(cfg_file)

arguments
    cfg_file = 'ideal_grid_cfg.json'
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
% Similar to 

% We can freely chose for which TR we want to simulate data, we use the TR from
% a previous project of ours. 
TR = 2.55;
% dt is the micro resolution, which is often 8 or 20, but for example when slice-time
% correction has been produced should be set to the number of slices.
dt = TR / 22; % Assuming classical SPM approach 
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
% These are hte options used by BayesPRF:
options = {};
options.name='simulation'; % PRF file tag
options.TE=0.03;      % Echo time (you may add other options from BayespRF)
options.model='spm_prf_response'; % We use a modified version of this function.


mu_simulations = {};
sigma_simulations ={};
point_simulations = {};

pc = 1;
for pn = cfg.param_names(:)'
    p_grid.(pn{1}) = [cfg.(['p_grid_' pn{1}])(:)', cfg.p_resolution];

    try 
        options.mu.(pn{1}) = cfg.(['p_mu_' pn{1}])(:)';
    catch
        disp("Mu boundaries not defined for simulation.")
    end

    mu_simulations{end + 1, 1} = cfg.(['sim_mu_' pn{1}])(:);
    sigma_simulations{end + 1, 1} = cfg.(['sim_sigma_' pn{1}])(:);
    point_simulations{end + 1, 1} = cfg.(['point_mu_' cfg.point_names{pc}])(:);
    pc = pc + 1;
end

pointParams = {};
point_mus = cell(length(point_simulations), 1);
[point_mus{:}] = ndgrid(point_simulations{:});

nvoxels = length(point_mus{1}(:));

pc = 1;
for  pn = cfg.point_names(:)'
    tmp_mu = point_mus{pc};
    for ii = 1 : length(tmp_mu(:))
        pointParams{ii}.([pn{1}]) = tmp_mu(ii);
    end
    pc = pc + 1;
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
% We will use CPM precompute, although we pretty much have U already, however,
% doing somw "awkward" preprocessing here, will help down the line.

% Load example data:
appertures = load('example_data/aps_Bars.mat');
appertures = appertures.ApFrm;

appertures_post = zeros(size(appertures, 3), cfg.p_resolution, cfg.p_resolution);

% We will now resize and pad the appertures, so they can be better handled.
% Because we use the default CPM model, we have to draw a inner box, which
% limits the location. We further assume that the appertures a slightly larger
% than the inner field.
for ii = 1  : size(appertures, 3)
        tmp_image = squeeze(appertures(:, :, ii));
        tmp_image = imresize(tmp_image, [cfg.p_resolution - 16, cfg.p_resolution - 16]);
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

% As we are using a simulation we can infer the number of scans from the onsets
% and the TR:
if isnan(nscans)
    nscans = ceil(max(data.ons) / TR);
end

%% Select grid model & precompute
model = "cpm_grid_SAMSRF"; % user defined model, see cpm_grid_template for details

fixedparams.resolution = cfg.p_resolution;

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
        y(:, new_idx) = y(:, vidx) +  normrnd(0, cfg.sim_noise(nidx), [nscans,1]);
        xyz(1,  new_idx) = xyz(1, vidx);
        xyz(2,  new_idx) = nidx;
        xyz_def{new_idx} = xyz_def{vidx};
    end
end

%% 
% Storing and saving the VOI in SPM format:
VOI.xY.y = y;
VOI.xY.XYZmm = xyz;
VOI.xyz_def = xyz_def;
% Mean values for completeness
VOI.y = mean(y, 2); % for completeness

voiname =  [simulationname, '_simVOI.mat'];
save(fullfile(simulationdir, voiname), 'VOI')

end