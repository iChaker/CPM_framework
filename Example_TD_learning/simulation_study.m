clear;
clc;
close all;
%% SETUP
REDO = true; % Overwrites everything instead of loading temporary files.
DATA_GENERATION = 'grid'; % Type of data generation process.

addpath('utils_rpe')
addpath('utils_simulation')
addpath(genpath('../CPM_Toolbox'))
addpath(genpath('../BayespRF-master'))
%% Simulation study to generate the data necessary for the computational
% parametric mapping paper. 
mkdir("tmpfiles"); % Create a directory for temporary files which will be ignored.
mkdir("simulationfiles");
% add spmpath:
addpath('/mnt/projects/LogUtil/Neuroimaging_RewardCoding/scripts/cpm/toolboxes/spm12/')
BASEDIR = pwd;
cd('/mnt/projects/LogUtil/Neuroimaging_RewardCoding/scripts/cpm/toolboxes/VBA-toolbox/')
VBA_setup
cd(BASEDIR);
% First step - simulate data
% We want data for the full model with values. For this we create an artificial
% VOI, for this grid over the following values for tau and eta:
voi_tau = [0, 1/3, 2/3, 1]; %-36.0437,  -0.6931, 0.6931, 36.0437];
voi_eta = [-1, 0, 1];
voi_noise =  [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3];

if strcmp(DATA_GENERATION, 'grid')
    gen_grid.alpha = [0, 1, 3];
    gen_grid.eta = [-1.5, 1.5, 3];
end

recovery_grid_points = 40;

% We do not yet really need to use the whole PRF structure, but can use our
% functions to simulate a U structure and generate data from it.
TR = 0.592;

file="tsvfile1.tsv";
dt = TR / 8; % Assuming classical SPM approach 
data = [];

[trial_n, ~, ~, S, stimuli, C, onsets] = cpm_prepare_experiment_data(file);
data.dur = ones(3 * trial_n, 1) .* 2 .* dt;
data.dt = ones(3 * trial_n, 1) .* dt;
data.ons = onsets;
data.C = C;
data.S = S;

% We need to create a simple SPM struct to create our PRF structure:
SPM = {};
SPM.xY.RT = TR;
SPM.swd = '';
% And we have to set a few options;
options.name='recovery'; % PRF file tag
options.TE=0.03;      % Echo time (you may add other options from BayespRF)
options.model='spm_prf_response';

% estimate dependent on TR:
nscans = ceil(max(onsets) / TR);
disp(nscans)
% any input data structure must contain the fields: ons dur dt
% Other than that, you may add fields in which ever way you want (in this
% case C and S)

%% Select grid model & precompute
model = "cpm_grid_RPE"; % user defined model, see cpm_grid_template for details

fixedparams.gamma=1;  % these fieldnames must be the same fieldnames used in the cpm_grid
fixedparams.Lambda=0.99;

% To generate real point estimates for the values of the paper (namely learning
% rates of 0, 1/3, 2/3, and 1, we have to resort to a little trick in
% precomputation:
[d1, d2] = ndgrid(voi_tau, voi_eta);
voi_params = [d1(:), d2(:)];


nvoxels = length(voi_params) * length(voi_noise);

y = zeros(nscans, nvoxels);
xyz = zeros(3, nvoxels);


if ~ exist('simulationfiles/simVOI.mat', 'file') || REDO

    for vidx = 1 : length(voi_params)
            grid = [];
            U_voi = [];

            if strcmp(DATA_GENERATION, 'point')
                grid.alpha = [voi_params(vidx, 1), voi_params(vidx, 1), 1];
                grid.eta = [voi_params(vidx, 2), voi_params(vidx, 2), 1];

                U_voi = cpm_precompute(model, grid, fixedparams, data, 'tmpfiles/simU.mat', true);
                tmpy = cpm_generate_from_U(U_voi, 1, TR, nscans, 'cpm_obv_int'); % simulating bold using priors

            elseif strcmp(DATA_GENERATION, 'grid')
 
                U_voi = cpm_precompute(model, gen_grid, fixedparams, data, 'tmpfiles/simU.mat', true);
                tmpPRF = cpm_specify(SPM, options, zeros(nscans, 1), ...
                                                       zeros(3, 1), U_voi, 'cpm_RF_Gaussian', 'cpm_obv_int', 'tmpfiles/');

                tmpPE = spm_vec(tmpPRF.M.pE{1});
                % Inverse of latent transformation
                lalpha = norminv((voi_params(vidx, 1) - gen_grid.alpha(1)) ./ (gen_grid.alpha(2) - gen_grid.alpha(1)));
                leta = norminv((voi_params(vidx, 2) - gen_grid.eta(1)) ./ (gen_grid.eta(2) - gen_grid.eta(1)));
                tmpPE([1, 2]) = [lalpha, leta];
                tmpPE = spm_unvec(tmpPE, tmpPRF.M.pE{1});
                tmpy =  spm_prf_response(tmpPE, tmpPRF.M, tmpPRF.U);

            end
            
            for nidx = 1 : length(voi_noise)
                tmp_idx = vidx +  (nidx - 1) * length(voi_params);
                y(:, tmp_idx) = tmpy +  normrnd(0, voi_noise(nidx),[nscans,1]);
                xyz(1, tmp_idx) = vidx;
                xyz(2, tmp_idx) = nidx;
            end
    end

    VOI.xY.y = y;
    VOI.xY.XYZmm = xyz;
    VOI.y = mean(y, 2); % for completeness

    save('simulationfiles/simVOI.mat', 'VOI')

else
    VOI = load('simulationfiles/simVOI.mat');
    VOI = VOI.VOI;
end

%% Parameter recovery

outpath='simulationfiles';  % PRF file path
% Now we can use precompute for the real grid:
grid = [];
grid.alpha = [0, 1, 40];
grid.eta = [-1.5, 1.5, 40];

U_recovery = cpm_precompute(model, grid, fixedparams, data, 'simulationfiles/U_recovery', REDO);

PRF = cpm_specify(SPM, options, VOI.xY.y, VOI.xY.XYZmm, U_recovery, 'cpm_RF_Gaussian', 'cpm_obv_int', outpath);

%% Estimation
voxels = []; % or select subset of voxels of interest
use_par = true; % to use a parallel for loop or not in voxel estimation;
PRF.M.noprint = 1; % to suppress spm outputs

if ~exist('simulationfiles/PRFn.mat', 'file') || REDO
PRFn = cpm_estimate(PRF, voxels, use_par);
save('simulationfiles/PRFn.mat', 'PRFn'),
else
    PRFn = load('simulationfiles/PRFn.mat');
    PRFn = PRFn.PRFn;
end
%%
% TODO add titles etc.
visualize_recovery(PRFn, [1 : 12] + (12 * (find(voi_noise == 0.0) - 1)), voi_params)

visualize_recovery(PRFn, [1 : 12] + (12 * (find(voi_noise == 0.3) - 1)), voi_params)
%% Perform model reduction

% Models in comparison
% 1 no-learning / no-utility
% 2 learning / no-utility
% 3 no-learning / utility
% 4 full model

if false
reducedF = zeros(3, nvoxels);
% Storing prior estimates (same across voxels)
pE = PRFn.M.pE{1};
pC = PRFn.M.pC{1};

% Using the following order of values in the vectors:
disp(pE)
for vidx = 1 : nvoxels
        % Null model
        rE = spm_vec(pE); 
        rC = spm_vec(pC);
        rE([3, 4]) = -4;  % setting all parameters of interest to small values except for lmu_eta, as 0 = no utility modulation
        rC([1, 2]) = 0; % Allowing no covariance.
        reducedF(1, vidx) = cpm_model_reduction(PRFn, vidx, rE, rC); %+ PRFn.F(vidx);
        % learning no utility:
        rE = spm_vec(pE); 
        rC = spm_vec(pC);
        rE([4]) = -4;  
        rC([2]) = 0; % Allowing no covariance.
        reducedF(2, vidx) = cpm_model_reduction(PRFn, vidx, rE, rC); % + PRFn.F(vidx);
        % no learning utility
        %         rE = spm_vec(pE); 
        %         rC = spm_vec(pC);
        %         rE([1, 3]) = -8.5;  
        %         rC([1, 3]) = 0; % Allowing no covariance
        %         reducedF(3, vidx) = cpm_model_reduction(PRFn, vidx, rE, rC) + PRFn.F(vidx);
        
%         reducedF(3, vidx) =  PRFn.F(vidx);
end

%% Extract Comparisons
% Param idx for null model:
nullidx = find(voi_params(:, 1) == 0 & voi_params(:, 2) == 0);
% Param idx for no-learning
nolearningidx =  find(voi_params(:, 1) == 0);
% Param idx for no-utility
noutilityidx =  find(voi_params(:, 2) == 0);
% Param idx for full models
fullidx = find(voi_params(:, 1) ~= 0 & voi_params(:, 2) ~= 0);

% Model Comparison for nullmodel:
null_voxels =  [nullidx : 12 : nvoxels];

o = {};

cc = 1;

logBF = zeros(2, length(null_voxels));
exceedenceP = zeros(3, length(null_voxels));
vba_options.DisplayWin = 0;

for nvidx = null_voxels
    
    logBF(1, cc)  = reducedF(1, nvidx) - reducedF(2, nvidx);
    logBF(2, cc)  = reducedF(1, nvidx) - reducedF(3, nvidx);
    [~, o] = VBA_groupBMC(reducedF(:, nvidx) + PRFn.F(:, nvidx));
    exceedenceP(:, cc) = o.ep;
    
    cc = cc + 1;
end

%%
figure;
bar(voi_noise, exceedenceP);
%%
x = voi_noise;
figure;
plot(x,logBF(1, :),'*','Color','#D95319');
hold on
plot(x,logBF(2, :),'^','Color','#EDB120');
xlim([0 0.35]);
% ylim([0 inf]);
title('Bayes factors of null model');
legend('reward-learning','utility-learning');
%%
nou_voxels =  repmat(noutilityidx, 1 , length(voi_noise)) +   [0 : length(voi_noise) - 1] * 12;
end