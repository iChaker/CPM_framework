% Clear console etc.
clear;
clc;
close all;
% Adding paths
addpath('utils_rpe')
addpath('utils_simulation')
addpath(genpath('../CPM_Toolbox'))
addpath(genpath('../BayespRF-master'))
addpath(genpath('../Example_TD_learning/'))
% Adding toolboxes, check that you have the right path here!
addpath('/mnt/projects/CPM/toolboxes/spm12/')
BASEDIR = pwd;
cd('/mnt/projects/CPM/toolboxes/VBA-toolbox/')
VBA_setup
cd(BASEDIR);
%%
% Generative Process:
% Point generative process 
%  Finding the noise variance for the BOLD signal to simulate a given SNR
%  levels. Previously estimates variance of simulated BOLD signal is 0.0083.

% Formula = SNR = 10 * log10(var_signal / var_noise). 
% Target SNRs = [20, 10, 2, -2, -10, -20].

var_sig = 7.1695;
target_snr = [-20, -10, -2, 2, 10, 20];
t_sd = sqrt(var_sig * 10 .^ (-target_snr / 10)); % Solving for standardeviation
% Result =   26.7759    8.4673    3.3709    2.1269    0.8467    0.2678
%%
REDO = false;

if ~isfile('simulationfiles/xy_simVOI.mat')  || REDO
	[VOI1, ~] = cpm_simulate_data_samsrf('samsrf_sim.json');
else
    VOI = load('simulationfiles/xy_simVOI.mat');
    VOI1 = VOI.VOI;
end

%% Combine VOI
VOI1.xY.XYZmm(3, :) = 1; % Model indicator

VOI.xY.y = [VOI1.xY.y];
VOI.xY.XYZmm = [VOI1.xY.XYZmm];
VOI.xyz_def  = [VOI1.xyz_def];
VOI.y = VOI1.y; % Just to be save;

%% Make PRFs
if ~isfile('simulationfiles/xyrec_PRFn.mat')  || REDO
    [prf_onetau] = cpm_simulation_prfs_samsrf(VOI, 'samsrf_rec.json');
else
    prf_onetau = 'simulationfiles/xyrec_PRFn.mat';
end

%% Collate PRFs
cc = 1;
PRFn = {};

for prf_path = {prf_onetau}

prf_path = prf_path{1};
PRF = load(prf_path);
PRF = PRF.PRF;

    if ~isfile([prf_path(1 : end - 4) '_estimated' prf_path(end - 3 : end)]) || REDO
        PRFn{cc} = cpm_estimate(PRF, [], true);
        save([prf_path(1 : end - 4) '_estimated' prf_path(end - 3 : end)], 'PRFn')
    else
        tmp = load([prf_path(1 : end - 4) '_estimated' prf_path(end - 3 : end)]);
        PRFn{cc}  =tmp.PRFn{cc};
    end

    cc = cc + 1;
end

%%  Extract indices
model_idx =VOI.xY.XYZmm(3, :);
noise_idx = VOI.xY.XYZmm(2, :);
noises = unique(noise_idx); % Noise levels
nnoise = length(noises);
num_models = sum(noise_idx == 1); % Number of different models
generating_params = VOI.xyz_def;
%%  Calculate SNR
simY = VOI.xY.y;
base = 1; % Index, with 0 noise
signal_var = var(simY(:, noise_idx == base)); % Variance at each voxel

snrs = zeros(size(simY, 2), 1); % pre allocate snrs 
mean_snr = zeros(length(noises), 1); % pre - allocate 

for nidx = noises
    noise_var = var(simY(:, noise_idx == nidx) - simY(:, noise_idx == base));
    snrs(noise_idx == nidx) = signal_var ./ (noise_var + eps);
    mean_snr(nidx) = mean(10 * log10(snrs(noise_idx == nidx)));
end
snr_label = round(mean_snr, 2);

%%
if true
post_samples = 400;
%
noise_level = 5;
unq_models = find(model_idx == 1 & noise_idx == noise_level);

% make axes
pad = 40;
fig_x = 500; fig_y = 500; 
fig2 = figure('Color', 'white', 'Units', 'pixels', 'Position', ...
                    [0, 0, fig_x + pad, fig_y + 2 * pad]);
axis('off')
hold on;

tiledlayout(4, 4,"TileSpacing","tight", "Padding", "tight");

for ii = unq_models

    nexttile()

        hold on;
        plot_single_voxel(PRFn{1}, ii, {'x', 'y'}, { [], []}, {[], []}, post_samples, 'response')
        xlim(PRFn{1}.U(1).grid.x(1 : 2))
        ylim(PRFn{1}.U(1).grid.y(1 : 2))
        xlabel('x')
        ylabel('y')
        hold on

        x = generating_params{ii}.mu_x;
        x_sig =generating_params{ii}.sigma_x;
        y = generating_params{ii}.mu_y;
        y_sig =generating_params{ii}.sigma_y;
        
        
        tmp_title1 = sprintf('\\mu_{x}= %4.2f, \\sigma_{x}= %4.2f', ...
                                        x, x_sig);
       
        tmp_title2 = sprintf('\\mu_{y}= %4.2f, \\sigma_{y}= %4.2f', ...
                                        y, y_sig);        
        tmpt =  title({tmp_title1; tmp_title2});
        tmpt.FontSize = 6;
        scatter(x, y, 5, 'filled', 'MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerFaceColor',[1 1 1],  'LineWidth',1.0)
        
        tp = -pi:0.01:pi;
        x_post = x + 2 *  x_sig .* cos(tp);
        y_post =y +  2 * y_sig .* sin(tp);
        plot(x_post, y_post)

end
sgtitle({'Parameter Recovery:', ['SNR:', num2str(snr_label(noise_level))]});

cpm_savefig(fig2, 'results/figS1_retino_parameter_recovery.png')
end

%%

[trues, preds] = deal({}, {});
for kk = 1
    mod_idx = find(model_idx == kk);
    trues{kk} = zeros(length(spm_vec(cpm_get_true_parameters(PRFn{kk}, 1))), length(mod_idx));
    preds{kk} = zeros(size(trues{kk}));
    
    jj = 1;
    for cidx = mod_idx
        tmp_true = cpm_get_true_parameters(PRFn{kk}.M.pE{cidx}, PRFn{kk}.M, PRFn{kk}.U);
        fn = fieldnames(generating_params{cidx});
        for fi=1:length(fn)
            tmp_true.(fn{fi}) = generating_params{cidx}.(fn{fi});
        end
        trues{kk}(:, jj) = spm_vec(tmp_true);
        preds{kk}(:, jj) = spm_vec(cpm_get_true_parameters(PRFn{kk}, cidx));
    jj = jj + 1;
    end
end
% Estimate MSE
error_fun =  @(true, pred) squeeze(sqrt(mean((true - pred).^2, 2)));

mses{1} = error_fun(reshape(trues{1}, 8, [], nnoise), reshape(preds{1}, 8, [], nnoise));
%%
fig4 = figure('Position', [0, 0, 1200, 600]);
sub_titles = {'Retinotopic Mapping Simulation'};
for nc = 1
    subplot(1, 1, nc)
    labels = fieldnames(cpm_get_true_parameters(PRFn{nc}, 1));
    labels = strrep(labels, '_', ' ');
    h = heatmap(round(mses{nc}, 4), 'YDisplayLabels', labels, 'XDisplayLabels', snr_label, 'XLabel', 'SNR', 'YLabel', 'Parameter');
    title(sub_titles{nc})
end

sgtitle('Parameter Recovery: RMSE')

cpm_savefig(fig4, 'results/fig_S2_rmse_retino_parameter_recovery.png')