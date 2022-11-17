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

REDO = false;

if ~isfile('simulationfiles/tau_PRFn.mat')  || REDO
	[VOI, ~] = cpm_simulate_data('tau_simulation_cfg.json');
else
    VOI = load('simulationfiles/tau_simVOI.mat');
    VOI = VOI.VOI;
end

%%
[prf_onetau] = cpm_simulation_prfs(VOI, 'one_tau_recovery_cfg.json');
[prf_twotau] =  cpm_simulation_prfs(VOI, 'two_tau_recovery_cfg.json');

%%

cc = 1;
PRFn = {};

for prf_path = {prf_onetau, prf_twotau}

    prf_path = prf_path{1};
PRF = load(prf_path);
PRF = PRF.PRF;


    if ~isfile([prf_path(1 : end - 4) '_estimated' prf_path(end-3 : end)]) || REDO
        PRFn{cc} = cpm_estimate(PRF, [], true);
        save([prf_path(1 : end - 4) '_estimated' prf_path(end-3 : end)], 'PRFn')
    else
        PRFn = load([prf_path(1 : end - 4) '_estimated' prf_path(end-3 : end)]);
        PRFn{cc}  = PRFn.PRFn;
    end

    cc = cc + 1;
end

%%
VOI = load(fullfile('simulationfiles/tau_simVOI.mat'));

simulation_params =VOI.VOI.xyz_def;

mid = VOI.VOI.xY.XYZmm;

noise = 1;
nnoise = 2;
nvoxels = length(mid) / nnoise;


points = find(mid(1, :) == 1 &mid(2, :) == noise);
%%


figure; 
tiles = tiledlayout(4, 4); 
inverse_fun = @(x, xmin, xmax) (norminv((x - xmin) ./ (xmax - xmin)));


for kk  =points
    nexttile()
    %plot_single_voxel(PRFn{2}, kk, {'tauneg', 'taupos'}, { @cpm_logit_inv, @cpm_logit_inv}, { 'alphaneg', 'alphapos'}, 1000, true)
    plot_single_voxel(PRFn{2}, kk, {'tauneg', 'taupos'}, { [], []}, {[], []}, 1000, true)
    % params
    x = inverse_fun(simulation_params{kk - (noise - 1) * nvoxels}.alphaneg, -eps, 1 + eps);
    y = inverse_fun(simulation_params{kk - (noise - 1) * nvoxels}.alphapos, -eps, 1 + eps); 
   
    scatter(x, y,'filled', 'MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerFaceColor',[1 1 1],  'LineWidth',1.5)
    xlim(PRFn{2}.U(1).grid.tauneg(1 : 2))
    ylim(PRFn{2}.U(1).grid.tauneg(1 : 2))

end
%%

%%
reducedF = zeros(2, length(mid));
% 
% reduced_val = -8.125;
% 
% reductions_re = {[1, 3, 4    ], [4    ], [1, 3]};
% reductions_rc = {[1, 2, 3, 4], [2, 4], [1, 3]};
% reduction_names = {'null', 'no-utility', 'no-learning', 'full'};
% 
for rd  = 1 :  length(PRFn) % length(reductions_re) + 1
    for vidx = 1 : length(mid)

    %         pE = PRFn.M.pE{vidx};
    %         pC = PRFn.M.pC{vidx};
    %         
    %         if rd <= length(reductions_re)
    %             rE0 = spm_vec(pE); 
    %             rC0 = spm_vec(pC);
    %             rE0(reductions_re{rd}) = reduced_val;  % setting all parameters of interest to small values except for lmu_eta, as 0 = no utility modulation
    %             rC0(reductions_rc{rd}) = 0; % Allowing no covariance.
    %             reducedF(rd, vidx) = cpm_model_reduction(PRFn, vidx, rE0, rC0) + PRFn.F(vidx);
    %         else
    %             reducedF(rd, vidx) =  PRFn.F(vidx);
    %         end
    reducedF(rd, vidx) =  PRFn{rd}.F(vidx);
    end
end

%%
nparams = nvoxels;

[logBF, exceedenceP, comparisons] = deal({}, {}, {});
% Param idx which generated 0 models:
nullidx =5;
% small trick to expand initial index to noise 
null_voxels =  repmat(nullidx, 1 , nnoise) +   (0 : nnoise - 1) * nparams;
[logBF{1}, exceedenceP{1}, comparisons{1}]  = cpm_model_comparison(reducedF, null_voxels, 1);

% Model comparisons for no utility model:
% Param idx for no-utility - including no learning, i.e. 0 model is included
for kk = 2 : 4
    no_idx = find(mid(1, :) == kk & mid(2, :) == 1)';
    no_vox =  repmat(no_idx, 1 , nnoise) +   (0 : nnoise - 1) * nparams;
    [logBF{kk}, exceedenceP{kk}, comparisons{kk}]  = cpm_model_comparison(reducedF, no_vox, 1);
end


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
    subplot(2, 4, ii)
    bar(voi_noise, exceedenceP{ii});
    legend({'Null', 'learning', 'utility', 'full'});
    title(ep_titles{ii})
    xlabel('Gaussian noise - SD')
    ylabel('Exceendence Probability')
    ylim([0, 1]);
end

sgtitle( 'Model Comparisons')

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

cpm_savefig(fig2, sprintf('results/model_recovery_%s.pdf', 'all'))

end
