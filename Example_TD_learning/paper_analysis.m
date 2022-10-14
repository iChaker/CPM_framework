% Clear console etc.
clear;
clc;
close all;
% Adding paths
addpath('utils_rpe')
addpath('utils_simulation')
addpath(genpath('../CPM_Toolbox'))
addpath(genpath('../BayespRF-master'))
% Adding toolboxes, check that you have the right path here!
addpath('/mnt/projects/CPM/toolboxes/spm12/')
BASEDIR = pwd;
cd('/mnt/projects/CPM/toolboxes/VBA-toolbox/')
VBA_setup
cd(BASEDIR);

%% Generative Process:
% Point generative process 

[prf_path, ~] = cpm_simulate_data('test_grid_cfg.json');

PRF = load( 'simulationfiles/test_PRFn.mat');
PRF = PRF.PRF;
PRFn = cpm_estimate(PRF, [], true);

save([prf_path(1 : end - 4) '_estimated' prf_path(end-3 : end)], 'PRFn')

%%
VOI = load(fullfile('simulationfiles/test_simVOI.mat'));

simulation_params =VOI.VOI.xyz_def;

mid = VOI.VOI.xY.XYZmm;

noise = 1;
nnoise = 7;
nvoxels = length(mid) / nnoise;


points = find(mid(1, :) == 1 &mid(2, :) == noise);
%%
figure; 
tiles = tiledlayout(3, 4); 

for kk  =points
    nexttile()
    plot_single_voxel(PRFn, kk, {'tau', 'eta'}, { @cpm_logit_inv, []}, { 'alpha', []}, 1000, true)
    
    % params
    x = simulation_params{kk - (noise - 1) * nvoxels}.alpha;
    y = simulation_params{kk - (noise - 1) * nvoxels}.eta;
   
    scatter(x, y,'filled', 'MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerFaceColor',[1 1 1],  'LineWidth',1.5)
end

%%
reducedF = zeros(4, length(mid));

reduced_val = -8.125;

for vidx = 1 : length(mid)
        pE = PRFn.M.pE{vidx};
        pC = PRFn.M.pC{vidx};
        % Null model 
        % FIXME : I have huge issues with this - this is how we have coded it
        % so far afaik.
        rE0 = spm_vec(pE); 
        rC0 = spm_vec(pC);
        rE0([1, 3, 4]) = reduced_val;  % setting all parameters of interest to small values except for lmu_eta, as 0 = no utility modulation
        rC0([1, 2, 3, 4]) = 0; % Allowing no covariance.
        reducedF(1, vidx) = cpm_model_reduction(PRFn, vidx, rE0, rC0) + PRFn.F(vidx);
        % learning no utility:
        rE10 = spm_vec(pE); 
        rC10 = spm_vec(pC);
        rE10([4]) = reduced_val;  
        rC10([2, 4]) = 0; % Allowing no covariance.
        reducedF(2, vidx) = cpm_model_reduction(PRFn, vidx, rE10, rC10)  + PRFn.F(vidx);
        % no learning utility
        rE01 = spm_vec(pE); 
        rC01 = spm_vec(pC);
        rE01([1, 3]) = reduced_val;  
        rC01([1, 3]) = 0; % Allowing no covariance
        reducedF(3, vidx) = cpm_model_reduction(PRFn, vidx, rE01, rC01) + PRFn.F(vidx);
        % Storing full model
        reducedF(4, vidx) =  PRFn.F(vidx);
end

%
nparams = nvoxels;

[logBF, exceedenceP, comparisons] = deal({}, {}, {});
% Param idx which generated 0 models:
nullidx =5;
% small trick to expand initial index to noise 
null_voxels =  repmat(nullidx, 1 , nnoise) +   (0 : nnoise - 1) * nparams;
[logBF{1}, exceedenceP{1}, comparisons{1}]  = cpm_model_comparison(reducedF, null_voxels, 1);

% Model comparisons for no utility model:
% Param idx for no-utility - including no learning, i.e. 0 model is included
noutilityidx = find(mid(1, :) == 2 & mid(2, :) == 1)';
nou_voxels =  repmat(noutilityidx, 1 , nnoise) +   (0 : nnoise - 1) * nparams;
[logBF{2}, exceedenceP{2}, comparisons{2}]  = cpm_model_comparison(reducedF, nou_voxels, 2);

% Model comparisons for no learning model:
% Param idx for no-learning - including no utility, i.e. 0 model is included
nolearningidx = find(mid(1, :) == 3 & mid(2, :) == 1)';
nol_voxels =  repmat(nolearningidx, 1 , nnoise) +   (0 : nnoise - 1) * nparams;
[logBF{3}, exceedenceP{3}, comparisons{3}]  = cpm_model_comparison(reducedF, nol_voxels, 3);

% Model comparisons for full model:
% Param idx for full models excluding all 0 models
fullidx = find(mid(1, : ) == 4 & mid(2, :) == 1)';
full_voxels =  repmat(fullidx, 1 , nnoise) +   (0 : nnoise - 1) * nparams;
[logBF{4}, exceedenceP{4}, comparisons{4}]  = cpm_model_comparison(reducedF, full_voxels, 4);

%% Result plotting
fig2 = figure('Color', 'none', 'Units', 'pixels', 'Position', [0, 0, 1600, 1200]);
ep_titles = {'Null-Model', 'No-utility', 'No-learning', 'Full'};

% pre specify colors
colors = [[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; ...
                [0.9290, 0.6940, 0.1250]; [0.4940, 0.1840, 0.5560]];
markers = {'o', 'x', 's', 'p'};

voi_noise = 1 : 7;
% Plot RFX exceedence probabilities
for ii = 1 : 4
    subplot(1, 4, ii)
    bar(voi_noise, exceedenceP{ii});
    legend({'Null', 'learning', 'utility', 'full'});
    title(ep_titles{ii})
    xlabel('Gaussian noise - SD')
    ylabel('Exceendence Probability')
    ylim([0, 1]);
end
% Plot Bayes Factor for the different comparisons.
% for ii = 1  : 4
%     subplot(2, 4, ii + 4)
%     h = {};
%     for jj = 1 : 3
%         h{jj} = semilogy(voi_noise, exp(squeeze(logBF{ii}(jj, :, :))), markers{comparisons{ii}(jj)}, ... 
%                        'Color', colors(comparisons{ii}(jj), :));
%         hold on
%     end
%     
%     larray = [h{1}; h{2}; h{3}];
%     
%     legend(larray(1 : length(h{1}): end), ep_titles(comparisons{ii}));
%     xlabel('Bayes Factor')
%     ylabel('Gaussian noise - SD')
%     
%     title(sprintf('Bayes factor in Favor of %s', ep_titles{ii}));
% end
sgtitle( 'Model Comparisons')
