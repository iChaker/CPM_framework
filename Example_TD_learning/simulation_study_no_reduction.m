% Clear console etc.
clear;
clc;
close all;

%%  ===========================================
% =================== SETUP ==================== =
% =============================================
% Random seed for reproducibility of noise etc.
rng(2022)
% Adding paths
addpath('utils_rpe')
addpath('utils_simulation')
addpath(genpath('../CPM_Toolbox'))
addpath(genpath('../BayespRF-master'))
% Adding toolboxes, check that you have the right path here!
addpath('/mnt/projects/LogUtil/Neuroimaging_RewardCoding/scripts/cpm.bu/toolboxes/spm12/')
BASEDIR = pwd;
cd('/mnt/projects/LogUtil/Neuroimaging_RewardCoding/scripts/cpm.bu/toolboxes/VBA-toolbox/')
VBA_setup
cd(BASEDIR);
% Creating directories to store simulation files 
mkdir("tmpfiles"); % Create a directory for temporary files which will be ignored.
mkdir("simulationfiles");
mkdir("results");
%% ======== SIMULATION SETUP =======================
REDO = true; % Overwrites everything instead of loading temporary files.
use_par = true; % to use a parallel for loop or not in voxel estimation;

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
DATA_GENERATION = 'grid'; % Type of data generation process.
assert(strcmp(DATA_GENERATION, 'grid') || strcmp(DATA_GENERATION, 'point'), ...
           'Data generation process has to be point or grid')

% Setup of simulation parameters of the VOI:
voi_alpha = [0, 1/3, 2/3, 1]; 
voi_eta = [-1, 0, 1];
voi_noise = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3];
nnoise = length(voi_noise);
nparams = length(voi_alpha) * length(voi_eta);
% The resolution of the different grids, generation grid will only be used when
% DATA_GENERATION = grid.
recovery_grid_resolution = 41;
generation_grid_resolution = 7;

% Support over the parameters in the generation and recovery grid.
alpha_support = [0, 1];
eta_support = [-1.5, 1.5];

% Creating the generation grid, if necessary
if strcmp(DATA_GENERATION, 'grid')
    gen_grid.alpha = [alpha_support, generation_grid_resolution];
    gen_grid.eta = [eta_support, generation_grid_resolution];
end

%% =============== SIMULATING DATA ===============
% Preparing behavioural data for the simulation:
% Implementation note here: cpm_precompute requires data to only have the fields
% dur, dt, and ons. Everything else is further passed to the generation
% functions. In this case S and C which contain the complete serial compound
% vector and the wealth trajectory.
data = [];
[trial_n, ~, ~, S, stimuli, C, onsets] = cpm_prepare_experiment_data(file);
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

% These are hte options used by BayesPRF:
options.name='recovery'; % PRF file tag
options.TE=0.03;      % Echo time (you may add other options from BayespRF)
options.model='spm_prf_response'; % We use a modified version of this function.

% As we are using a simulation we can infer the number of scans from the onsets
% and the TR:
if isnan(nscans)
    nscans = ceil(max(onsets) / TR);
end

%% Select grid model & precompute
model = "cpm_grid_RPE"; % user defined model, see cpm_grid_template for details

fixedparams.gamma=0.97;  % these fieldnames must be the same fieldnames used in the cpm_grid
fixedparams.Lambda=1.0;

% We create all combinations of parameters eta and alpha here for our simulated VOI:
[d1, d2] = ndgrid(voi_alpha, voi_eta);
voi_params = [d1(:), d2(:)];
% The number of voxels is hte number of combinations and the number of noise
% levels:
nvoxels = nparams * nnoise;
% Pre allocating:
y = zeros(nscans, nvoxels);
states = [];
% The voxel location is meaningless in our simulation, 
% but we used it to store noise level and indexes in the parameter space
xyz = zeros(3, nvoxels); 

if ~ exist('simulationfiles/simVOI.mat', 'file') || REDO

    for vidx = 1 : nparams
        % For safety clearing previously set values
        grid = [];
        U_voi = [];
        % BOLD generation block:
        if strcmp(DATA_GENERATION, 'point')
            % In the point process of data generation we create an RPE sequence
            % for a single point (this could be cleaner, but staying like this
            % for now). 
            % Create a grid from the parameter settings:
            grid.alpha = [voi_params(vidx, 1), voi_params(vidx, 1), 1];
            grid.eta = [voi_params(vidx, 2), voi_params(vidx, 2), 1];
            % Generate a grid:
            U_voi = cpm_precompute(model, grid, fixedparams, data, 'tmpfiles/simU.mat', true);
            % Inferring BOLD response from BOLD priors using the RPE sequence
            % directly as neural input.
            [tmpy, stat] = cpm_generate_from_U(U_voi, 1, TR, nscans, 'cpm_obv_int');
            states = [states, stat];
            
        elseif strcmp(DATA_GENERATION, 'grid')
            % in the grid generation we use a grid over a number of points and
            % generat the BOLD response over the grid using the generative
            % process described by the CPM model
            if REDO && vidx == 1 % the Grid needs to be oncly created once, so here we load it or recreat it if necessary
                U_voi = cpm_precompute(model, gen_grid, fixedparams, data, 'tmpfiles/simU_grid.mat', true);
            else
                 U_voi = cpm_precompute(model, gen_grid, fixedparams, data, 'tmpfiles/simU_grid.mat', false);
            end
            % We generate a temporary PRF file, to generate prior values and some additional options:
            tmpPRF = cpm_specify(SPM, options, zeros(nscans, 1), ...
                                                  zeros(3, 1), U_voi, 'cpm_RF_Gaussian', 'cpm_obv_int', 'tmpfiles/');
            % We get the prior values of the latent(!!) parameters
            tmpPE = tmpPRF.M.pE{1};
            % As we define the to be simulated values in native space we have to
            % transform them into the laten space:
            lalpha = norminv( (voi_params(vidx, 1)  - gen_grid.alpha(1)) ./ (gen_grid.alpha(2) - gen_grid.alpha(1)));
            leta = norminv((voi_params(vidx, 2) - gen_grid.eta(1)) ./ (gen_grid.eta(2) - gen_grid.eta(1)));
            % we overwrite values accordingly in the prior structure:
             tmpPE.lmu_alpha = lalpha;
            tmpPE.lmu_eta = leta;
            tmpPE.lsigma_alpha = norminv(exp(-5) ./ sqrt(2)); % Sigma is bounded in our case between 0.001 and sqrt(2) 
            tmpPE.lsigma_eta = norminv(exp(-5) ./ sqrt(2)); % 
            disp(cpm_get_true_parameters(tmpPE, tmpPRF.M, tmpPRF.U))
            [tmpy, stat] =  spm_prf_response(tmpPE, tmpPRF.M, tmpPRF.U);
            % Saving hidden states of the model for later investigation
            states = [states, stat.u];

        end
        % After having generated the BOLD signal we add simulated noise:
        for nidx = 1 : nnoise
            tmp_idx = vidx +  (nidx - 1) * length(voi_params);
            % normrand with 0 variance returns the distributions mean:
            y(:, tmp_idx) = tmpy +  normrnd(0, voi_noise(nidx), [nscans,1]);
            xyz(1, tmp_idx) = vidx;
            xyz(2, tmp_idx) = nidx;
        end
    end
    
    % Storing and saving the VOI in SPM format:
    VOI.xY.y = y;
    VOI.xY.XYZmm = xyz;
    % Mean values for completeness
    VOI.y = mean(y, 2); % for completeness
    save('simulationfiles/simVOI.mat', 'VOI')

else
    VOI = load('simulationfiles/simVOI.mat');
    VOI = VOI.VOI;
end

%% ===================== MODEL INVERSION ======================
outpath='simulationfiles';  % PRF file path
% Now we generate / precompute the recovery grid:
grid = [];
grid.alpha = [alpha_support, recovery_grid_resolution];
grid.eta = [eta_support, recovery_grid_resolution];

U_recovery = cpm_precompute(model, grid, fixedparams, data, 'simulationfiles/U_recovery', REDO);

% And specify the PRF for recovery:
PRF = cpm_specify(SPM, options, VOI.xY.y, VOI.xY.XYZmm, ...
                   U_recovery, 'cpm_RF_Gaussian', 'cpm_obv_int', outpath);

voxels = []; % or select subset of voxels of interest
PRF.M.noprint = 0; % to suppress spm outputs

if ~exist('simulationfiles/PRFn.mat', 'file') || REDO
    PRFn = cpm_estimate(PRF, voxels, use_par);
    save('simulationfiles/PRFn.mat', 'PRFn'),
else
    PRFn = load('simulationfiles/PRFn.mat');
    PRFn = PRFn.PRFn;
end

%% ===================== VISUALIZATION OF RECOVERED PARAMS
nlevel = 0.1;
visualize_idx =  1 + (nparams * (find(voi_noise == nlevel) - 1)) : nparams + (nparams * (find(voi_noise == nlevel) - 1));
fig1 = visualize_recovery(PRFn, visualize_idx, voi_params, 4, true, 40, {'alpha', 'eta'}, 500);
sgtitle(sprintf('Parameter recovery - Gaussian noise with SD %4.2f', nlevel))

cpm_savefig(fig1,  sprintf('results/parameter_recovery_noise_%4.2f_%s.pdf', nlevel, DATA_GENERATION))

%%
%% Perform model reduction
% Models in comparison
% 1 no-learning / no-utility
% 2 learning / no-utility
% 3 no-learning / utility
% 4 full model

% I think this is how the current model comparison is implemented and I have
% quite some issues with it. We shold discuss this!!! 
reducedF = zeros(4, nvoxels);

reduced_sigma = -4; 2.5; %4;
rc_sigma = 1; %0.5;
reduced_alpha = -4;

for vidx = 1 : nvoxels
        pE = PRFn.M.pE{vidx};
        pC = PRFn.M.pC{vidx};
        % Null model 
        % FIXME : I have huge issues with this - this is how we have coded it
        % so far afaik.
        rE0 = spm_vec(pE); 
        rC0 = spm_vec(pC);
        rE0([1]) =  reduced_alpha;  % setting all parameters of interest to small values except for lmu_eta, as 0 = no utility modulation
        rE0([3, 4]) = reduced_sigma; 
        rC0([1, 2]) = 0; % Allowing no covariance.
         rC0([3, 4]) = rc_sigma; % Allowing no covariance.
        reducedF(1, vidx) = cpm_model_reduction(PRFn, vidx, rE0, rC0) + PRFn.F(vidx);
        % learning no utility:
        rE10 = spm_vec(pE); 
        rC10 = spm_vec(pC);
        rC10([2]) = 0; % Allowing no covariance.
        rC10([4]) = rc_sigma; % Allowing no covariance.
        rE10([4]) = reduced_sigma; 

        reducedF(2, vidx) = cpm_model_reduction(PRFn, vidx, rE10, rC10)  + PRFn.F(vidx);
        % no learning utility
        rE01 = spm_vec(pE); 
        rC01 = spm_vec(pC);
        rE01([1]) = reduced_alpha;   
         rE01([3]) =reduced_sigma; 

        rC01([1]) = 0; % Allowing no covariance
        rC01([3]) = rc_sigma; % Allowing no covariance
        reducedF(3, vidx) = cpm_model_reduction(PRFn, vidx, rE01, rC01) + PRFn.F(vidx);
        % Storing full model
        reducedF(4, vidx) =  PRFn.F(vidx);
end
%
% Model Comparison for nullmodel:
[logBF, exceedenceP, comparisons] = deal({}, {}, {});
% Param idx which generated 0 models:
nullidx = find(voi_params(:, 1) == 0 & voi_params(:, 2) == 0);
% small trick to expand initial index to noise 
null_voxels =  repmat(nullidx, 1 , nnoise) +   (0 : nnoise - 1) * nparams;
[logBF{1}, exceedenceP{1}, comparisons{1}]  = cpm_model_comparison(reducedF, null_voxels, 1);

% Model comparisons for no utility model:
% Param idx for no-utility - including no learning, i.e. 0 model is included
noutilityidx =  find(voi_params(:, 2) == 0 & voi_params(:, 1) ~= 0);
nou_voxels =  repmat(noutilityidx, 1 , nnoise) +   (0 : nnoise - 1) * nparams;
[logBF{2}, exceedenceP{2}, comparisons{2}]  = cpm_model_comparison(reducedF, nou_voxels, 2);

% Model comparisons for no learning model:
% Param idx for no-learning - including no utility, i.e. 0 model is included
nolearningidx =  find(voi_params(:, 1) == 0 & voi_params(:, 2) ~= 0 );
nol_voxels =  repmat(nolearningidx, 1 , nnoise) +   (0 : nnoise - 1) * nparams;
[logBF{3}, exceedenceP{3}, comparisons{3}]  = cpm_model_comparison(reducedF, nol_voxels, 3);

% Model comparisons for full model:
% Param idx for full models excluding all 0 models
fullidx = find(voi_params(:, 1) ~= eps & voi_params(:, 2) ~= 0);
full_voxels =  repmat(fullidx, 1 , nnoise) +   (0 : nnoise - 1) * nparams;
[logBF{4}, exceedenceP{4}, comparisons{4}]  = cpm_model_comparison(reducedF, full_voxels, 4);

% Result plotting
fig2 = figure('Color', 'none', 'Units', 'pixels', 'Position', [0, 0, 1600, 1200]);
ep_titles = {'Null-Model', 'No-utility', 'No-learning', 'Full'};

% pre specify colors
colors = [[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; ...
                [0.9290, 0.6940, 0.1250]; [0.4940, 0.1840, 0.5560]];
markers = {'o', 'x', 's', 'p'};

% Plot RFX exceedence probabilities
for ii = 1 : 4
    subplot(2, 4, ii)
    bar(voi_noise, exceedenceP{ii});
    legend({'Null', 'learning', 'utility', 'full'});
    title(ep_titles{ii})
    xlabel('Gaussian noise - SD')
    ylabel('Exceendence Probability')
    ylim([0, 1]);
end
% Plot Bayes Factor for the different comparisons.
for ii = 1  : 4
    subplot(2, 4, ii + 4)
    h = {};
    for jj = 1 : 3
        h{jj} = semilogy(voi_noise, exp(squeeze(logBF{ii}(jj, :, :))), markers{comparisons{ii}(jj)}, ... 
                       'Color', colors(comparisons{ii}(jj), :));
        hold on
    end
    
    larray = [h{1}; h{2}; h{3}];
    
    legend(larray(1 : length(h{1}): end), ep_titles(comparisons{ii}));
    xlabel('Bayes Factor')
    ylabel('Gaussian noise - SD')
    
    title(sprintf('Bayes factor in Favor of %s', ep_titles{ii}));
end
sgtitle( 'Model Comparisons')

cpm_savefig(fig2, sprintf('results/model_recovery_%s.pdf', DATA_GENERATION))

%% Parameter recovery statistics
% it gets a bit messy here

% Extracting and saving the true parameters, which generated the data in the
% simulation.
rec_params = zeros(8, nvoxels);
true_params = zeros(8, nvoxels);

for vidx = 1 : nvoxels
    % Posterior values
    tmpparams = cpm_get_true_parameters(PRFn.Ep{vidx}, PRF.M, PRF.U);
    tmpparams = spm_vec(tmpparams);
    rec_params(:, vidx) = tmpparams;
    
    % Prior values (they've mostly been used in the simulation)
    tmpparams = cpm_get_true_parameters(PRFn.M.pE{vidx}, PRF.M, PRF.U);
    tmpparams = spm_vec(tmpparams);
    true_params(:, vidx) = tmpparams;
    
    true_params(1 : 2, vidx) = voi_params(mod(vidx -1, nparams) + 1, :);
    true_params(3 : 4, vidx) = [exp(-5), exp(-5)];
end

% Reshaping to shape: parameters in model, voxels, noise levels
 true_params = reshape(true_params, 8, [], nnoise);
 rec_params = reshape(rec_params, 8, [], nnoise);

 %  RMSE (Euclidean distance)
rmse = squeeze(sqrt(mean((true_params - rec_params).^2, 2)));
fig3 = figure('Color', 'none', 'Units', 'pixels', 'Position', [0, 0, 1800, 1600]);
param_names =  {'\mu_\alpha', '\mu_\eta', '\sigma_\alpha', '\sigma_\eta', '\beta', 'transit', 'decay', '\epsilon'};
noise_names =strsplit(num2str(voi_noise));

heatmap(param_names, noise_names, rmse');
xlabel('Parameter');
ylabel('Noise level');
title('Euclidean distane to true parameters');

cpm_savefig(fig3, sprintf('results/parameter_distance_noise_%4.2f_%s.pdf', nlevel, DATA_GENERATION))


%% Visualization sugar
% here it gets messy, not sure yet what these parts could be really used for,
% but it allows to plot and estimate predicted hidden states for the model and
% comparison to the true values, stored in "y" and states, states is maximally
% of length nparams, because they are not influence by noise.
% y_pred = [];
% z_pred = [];
% 
% for vidx = 1 : nvoxels
%     % genrating predictions
%     [tmpy, tmpz] = spm_prf_response(PRFn.Ep{vidx}, PRFn.M, PRFn.U);
%     y_pred = [y_pred, tmpy];
%     z_pred = [z_pred, tmpz.u];
% end
%%
% figure;
% 
% for ii = 1% : nparams
% subplot(3, 4, ii)
% plot(z_pred(:, 12) - states(:, ii))
% hold on
% end