% Clear console etc.
clear;
clc;
close all;
% Adding paths
addpath('../utils_rpe')
addpath('../utils_simulation')
addpath(genpath('../../CPM_Toolbox'))
addpath(genpath('../../BayespRF-master'))
% Adding toolboxes, check that you have the right path here!
addpath('/mnt/projects/CPM/toolboxes/spm12/')
BASEDIR = pwd;
cd('/mnt/projects/CPM/toolboxes/VBA-toolbox/')
VBA_setup
cd(BASEDIR);

% Generative Process:
%% Point generative process 
forward_fun = @(x, xmin, xmax) ((xmax-xmin) .* spm_Ncdf(x,0,1)) + xmin;
inverse_fun = @(x, xmin, xmax) (norminv((x - xmin) ./ (xmax - xmin)));

%%
REDO = false;

if ~isfile('simulationfiles/learn_utility_simVOI.mat')  || REDO
	[VOI1, ~] = cpm_simulate_data('learn_utility_simulation.json');
else
    VOI = load('simulationfiles/learn_utility_simVOI.mat');
    VOI1 = VOI.VOI;
end

if ~isfile('simulationfiles/utility_simVOI.mat')  || REDO
	[VOI2, ~] = cpm_simulate_data('utility_simulation.json');
else
    VOI = load('simulationfiles/utility_simVOI.mat');
    VOI2 = VOI.VOI;
end

if ~isfile('simulationfiles/learn_simVOI.mat')  || REDO
	[VOI3, ~] = cpm_simulate_data('learn_simulation.json');
else
    VOI = load('simulationfiles/learn_simVOI.mat');
    VOI3 = VOI.VOI;
end
%%
% Combine VOI
VOI1.xY.XYZmm(3, :) = 1;
VOI2.xY.XYZmm(3, :) = 2;
VOI3.xY.XYZmm(3, :) = 3;

VOI.xY.y = [VOI1.xY.y, VOI2.xY.y, VOI3.xY.y];
VOI.xY.XYZmm = [VOI1.xY.XYZmm, VOI2.xY.XYZmm, VOI3.xY.XYZmm];
VOI.xyz_def  = [VOI1.xyz_def; VOI2.xyz_def; VOI3.xyz_def];
VOI.y = VOI1.y; % Just to be save;

%%
if ~isfile('simulationfiles/learn_utility_PRFn.mat')  || REDO
    [prf_learn_utility] = cpm_simulation_prfs(VOI, 'learn_utility_recovery_cfg.json');
else
    prf_learn_utility = 'simulationfiles/learn_utility_PRFn.mat';
end

if ~isfile('simulationfiles/utility_PRFn.mat') || REDO
    [prf_utility] =  cpm_simulation_prfs(VOI, 'utility_recovery_cfg.json');
else
    prf_utility = 'simulationfiles/utility_PRFn.mat';
end

if ~isfile('simulationfiles/learn_PRFn.mat') || REDO
    [prf_learn] =  cpm_simulation_prfs(VOI, 'learn_recovery_cfg.json');
else
    prf_learn = 'simulationfiles/learn_PRFn.mat';
end
%%
cc = 1;
PRFn = {};

for prf_path = {prf_learn_utility, prf_utility, prf_learn}

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

%% Plot parameter recovery
gen_idx = VOI.xY.XYZmm;
model_idx = gen_idx(3, :);
noise_idx = gen_idx(2, :);
noises = unique(gen_idx(2, :));
nnoise = length(noises);
generating_params = VOI.xyz_def;

%%
fig5 = figure('Color', 'white', 'Units', 'pixels', 'Position', ...
                    [0, 0, 1000, 1000]);

tiledlayout(3, 3,  'TileSpacing', 'compact', 'Padding', 'compact')
sets{1} = [7:9, 4: 6, 1: 3] + 36; % [1: 9];

for jj = 1 : 9
    nexttile()
plot_single_voxel(PRFn{1}, sets{1}(jj), {'tau', 'eta'}, { @cpm_logit_inv, []}, {'alpha', []}, 1000, true)

hold on
x = generating_params{sets{1}(jj)}.mu_tau;
x_sig =generating_params{sets{1}(jj)}.sigma_tau;
y = generating_params{sets{1}(jj)}.mu_eta;
y_sig = generating_params{sets{1}(jj)}.sigma_eta;

scatter(cpm_logit_inv(x), y,'filled', 'MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerFaceColor',[1 1 1],  'LineWidth',1.0)

tp = -pi:0.01:pi;
x_post = x + 2 *  x_sig .* cos(tp);
y_post =y +  2 * y_sig .* sin(tp);
plot(cpm_logit_inv(x_post), y_post)
xlabel('\alpha')
ylabel('\eta')

xlim(cpm_logit_inv(PRFn{1}.U(1).grid.tau(1 : 2)));
ylim(PRFn{1}.U(1).grid.eta(1 : 2))

title(['\alpha = ' num2str(round(cpm_logit_inv(x), 2))  ', \eta = ' num2str(round(y, 2))]);


end

exportgraphics(fig5,['recovery' '.png'],'Resolution',600)

%%
fig1 = figure('Color', 'white', 'Units', 'pixels', 'Position', ...
                    [0, 0, 1000, 1000]);

tiledlayout(3, 3,  'TileSpacing', 'compact', 'Padding', 'compact')
sets{1} = [7:9, 4: 6, 1: 3] + 36; % [1: 9];

for jj = 1 : 9
    nexttile()
plot_single_voxel(PRFn{1}, sets{1}(jj), {'tau', 'eta'}, {[], []}, {[], []}, 1000, true)

hold on
x = generating_params{sets{1}(jj)}.mu_tau;
x_sig =generating_params{sets{1}(jj)}.sigma_tau;
y = generating_params{sets{1}(jj)}.mu_eta;
y_sig = generating_params{sets{1}(jj)}.sigma_eta;

scatter(x, y,'filled', 'MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerFaceColor',[1 1 1],  'LineWidth',1.0)

tp = -pi:0.01:pi;
x_post = x + 2 *  x_sig .* cos(tp);
y_post =y +  2 * y_sig .* sin(tp);
plot(x_post, y_post)
xlabel('\tau')
ylabel('\eta')

xlim(PRFn{1}.U(1).grid.tau(1 : 2));
ylim(PRFn{1}.U(1).grid.eta(1 : 2))
pbaspect([1 1 1])

title(['\alpha = ' num2str(round(cpm_logit_inv(x), 2))  ', \eta = ' num2str(round(y, 2))]);

end

exportgraphics(fig1,['recovery' '.png'],'Resolution', 72)
%%
num_models = sum(noise_idx == 1);
simY = VOI.xY.y;
% SNRs
base = 1;

signal_var = var(simY(:, noise_idx == base));

snrs = zeros(size(simY, 2), 1);
mean_snr = zeros(length(noises), 1);

for nidx = noises
    noise_var = var(simY(:, noise_idx == nidx) - simY(:, noise_idx == base));
    snrs(noise_idx == nidx) = signal_var ./ (noise_var + eps);
    mean_snr(nidx) = mean(log10(snrs(noise_idx == nidx)));
end

%
voi_sims = {zeros(size(simY)), zeros(size(simY))};
for pidx = 1 : size(gen_idx,2)
    voi_sims{1}(:, pidx) = spm_prf_response(PRFn{1}.Ep{pidx}, PRFn{1}.M, PRFn{1}.U);
    voi_sims{2}(:, pidx) = spm_prf_response(PRFn{2}.Ep{pidx}, PRFn{2}.M, PRFn{2}.U);
    voi_sims{3}(:, pidx) = spm_prf_response(PRFn{3}.Ep{pidx}, PRFn{3}.M, PRFn{3}.U);
end
%
generating_params = VOI.xyz_def;

modelF = zeros(3, size(model_idx, 2));

for vidx = 1 : size(model_idx, 2)
    for midx = 1 : size(PRFn, 2)
        modelF(midx, vidx) = PRFn{midx}.F(vidx);
    end
end

[logBF, exceedenceP, comparisons] = deal({}, {}, {});

for kk = 1 : 3
    no_idx = find(gen_idx(3, :) == kk & gen_idx(2, :) == 1)';
    no_vox =  repmat(no_idx, 1 , nnoise) +   (0 : nnoise - 1) * length(no_idx);
    [logBF{kk}, exceedenceP{kk}, comparisons{kk}]  = cpm_model_comparison(modelF, no_vox, kk);
end

%%
ep_titles = {'Learn Utility', 'Utility', 'Learn'};

% pre specify colors
colors = [[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; ...
                [0.9290, 0.6940, 0.1250]; [0.4940, 0.1840, 0.5560]];
markers = {'o', 'x', 's', 'p'};
%
figure;
tiledlayout(1, 3)
% Plot RFX exceedence probabilities
for ii = 1 : 3
    nexttile()
    bar(noises, exceedenceP{ii}');
    legend({'Learn Utility', 'Utility', 'Learn'});
    title(ep_titles{ii})
    xlabel('Gaussian noise - SD')
    ylabel('Exceedance Probability')
    ylim([0, 1]);
    xticklabels(round(mean_snr, 2))
end

sgtitle('Model Recovery')
%%
t = -8.1259 : 0.1982 : 8.1259;

fig1 = figure('Color', 'white', 'Units', 'pixels', 'Position', ...
                    [0, 0, 500, 800]);

tiledlayout(12, 12, 'TileSpacing', 'none', 'Padding', 'compact')

sets = {[1, 2] + 4, [3, 4] + 4};

for ii = 1 :2

for jj = 1 : 2
    nexttile([2, 1])
    axis('off')
    nexttile([2, 4]);
    plot_single_voxel(PRFn{1}, sets{ii}(jj), {'tau'}, { []}, {[]}, 10, true)
    hold on;
    y = normpdf(t, generating_params{sets{ii}(jj)}.mu_tau, generating_params{sets{ii}(jj)}.sigma_tau);
    y = y ./ sum(y);
    plot(t, y, 'LineWidth', 1.5);
    xlim(PRFn{1}.U(1).grid.tau(1 : 2));
    %ylim([0, 1]);
    tmp_title = sprintf('\\mu_\\tau= %4.2f, \\sigma_\\tau= %4.2f', ...
        generating_params{sets{ii}(jj)}.mu_tau, ...
    generating_params{sets{ii}(jj)}.sigma_tau);
   
title(tmp_title);
    xlabel('\tau')
    nexttile([2, 1])
    axis('off')
end

for jj = 1 : 2
nexttile([4, 6])

plot_single_voxel(PRFn{2}, sets{ii}(jj), {'tauneg', 'taupos'}, { [], []}, {[], []}, 10, true)
          xlim(PRFn{2}.U(1).grid.tauneg(1 : 2))
          ylim(PRFn{2}.U(1).grid.tauneg(1 : 2))
          xlabel('\tau^-')
          ylabel('\tau^+')
end

end

sgtitle('Parameter Recovery: Classic RL')

t = -8.1259 : 0.1982 : 8.1259;
%%
fig2 = figure('Color', 'white', 'Units', 'pixels', 'Position', ...
                    [0, 0, 1000, 1200]);
tiledlayout(24, 16, 'TileSpacing', 'none', 'Padding', 'none')


sets = {[13:16] + 16, [17:20] + 16, [21 : 24] + 16, [25 : 28] + 16};

for ii = 1 :4

for jj = 1 : 4
    nexttile([2, 1])
    axis('off')
    nexttile([2, 2]);
    plot_single_voxel(PRFn{1}, sets{ii}(jj), {'tau'}, { []}, {[]}, 10, true)
    
    xlim(PRFn{1}.U(1).grid.tau(1 : 2));
    
        tmp_title = sprintf('\\mu_{\\tau^+}= %4.2f, \\sigma_{\\tau^+}= %4.2f', ...
        generating_params{sets{ii}(jj)}.mu_taupos, ...
    generating_params{sets{ii}(jj)}.sigma_taupos);
   
        tmp_title = [tmp_title, sprintf('\\mu_{\\tau^-}= %4.2f, \\sigma_{\\tau^-}= %4.2f', ...
        generating_params{sets{ii}(jj)}.mu_tauneg, ...
    generating_params{sets{ii}(jj)}.sigma_tauneg)];
   
title(tmp_title);
    
    nexttile([2, 1])
    axis('off')
end

for jj = 1 : 4
nexttile([4, 4])

plot_single_voxel(PRFn{2}, sets{ii}(jj), {'tauneg', 'taupos'}, { [], []}, {[], []}, 10, true)

hold on
x = generating_params{sets{ii}(jj)}.mu_tauneg;
x_sig =generating_params{sets{ii}(jj)}.sigma_tauneg;
y = generating_params{sets{ii}(jj)}.mu_taupos;
y_sig =generating_params{sets{ii}(jj)}.sigma_taupos;

scatter(x, y,'filled', 'MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerFaceColor',[1 1 1],  'LineWidth',1.0)

tp = -pi:0.01:pi;
x_post = x + 2 *  x_sig .* cos(tp);
y_post =y +  2 * y_sig .* sin(tp);
plot(x_post, y_post)
xlabel('\tau^-')
ylabel('\tau^+')

xlim(PRFn{2}.U(1).grid.tauneg(1 : 2))
ylim(PRFn{2}.U(1).grid.tauneg(1 : 2))

end

end

sgtitle('Parameter Recovery: Distributional RL')

% Important (but later)
% Label plots

%% Estimation error
% Classic
classic_index = find(model_idx == 1);
classic_noise = noise_idx(classic_index);

classic_error = zeros(length(spm_vec(cpm_get_true_parameters(PRFn{1}, 1))), length(classic_index));

for cidx = classic_index
    tmp_gen_vals = spm_vec(cpm_get_true_parameters(PRFn{1}.M.pE{cidx}, PRFn{1}.M, PRFn{1}.U));
    tmp_gen_vals(1 : 2) = [generating_params{cidx}.mu_tau,  generating_params{cidx}.sigma_tau];
    classic_error(:, cidx) = tmp_gen_vals - spm_vec(cpm_get_true_parameters(PRFn{1}, cidx));
end

%% BPA
cl_lambda = zeros(6, 6);

for cidx = classic_index(classic_noise == 1)
    cl_lambda = cl_lambda + PRFn{1}.Cp{cidx};
end

cl_c = pinv(cl_lambda);

%% 
% Classic
classic_index = find(model_idx == 1);
classic_noise = noise_idx(classic_index);

classic_error = zeros(length(spm_vec(cpm_get_true_parameters(PRFn{1}, 1))), length(classic_index));

for cidx = classic_index
    tmp_gen_vals = spm_vec(cpm_get_true_parameters(PRFn{1}.M.pE{cidx}, PRFn{1}.M, PRFn{1}.U));
    tmp_gen_vals(1 : 2) = [generating_params{cidx}.mu_tau,  generating_params{cidx}.sigma_tau];
    classic_error(:, cidx) = tmp_gen_vals - spm_vec(cpm_get_true_parameters(PRFn{1}, cidx));
end

%%
% Distributional
dist_index = find(model_idx == 2);
dist_noise = noise_idx(dist_index);

dist_error = zeros(length(spm_vec(cpm_get_true_parameters(PRFn{2}, 1))), length(dist_index));

for cidx = dist_index
    tmp_gen_vals = spm_vec(cpm_get_true_parameters(PRFn{2}.M.pE{cidx}, PRFn{2}.M, PRFn{2}.U));
    tmp_gen_vals(1 : 2) = [generating_params{cidx}.mu_tauneg,  generating_params{cidx}.mu_taupos];
     tmp_gen_vals(3 : 4) = [generating_params{cidx}.sigma_tauneg,  generating_params{cidx}.sigma_taupos];
    dist_error(:, cidx) = tmp_gen_vals - spm_vec(cpm_get_true_parameters(PRFn{2}, cidx));
end

%% Summarize errors for over noise
dist_mse = zeros(size(dist_error, 1), nnoise);
for nidx =  1 : nnoise
    dist_mse(:, nidx) = mean(dist_error(:, dist_noise == nidx).^2, 2);
end

classic_mse = zeros(size(classic_error, 1), nnoise);
for nidx =  1 : nnoise
    classic_mse(:, nidx) = mean(classic_error(:, classic_noise == nidx).^2, 2);
end
%%
figure;
h = heatmap(sqrt(dist_mse));
labels = fieldnames(cpm_get_true_parameters(PRFn{2}, 1));
h.YDisplayLabels = labels;
%% Classic BPA

included = classic_index(classic_noise == 1);
nincluded = length(included);
GCM = cell(nincluded,1);
i = 1;
for v = included
    GCM{i}.Cp = PRFn{1}.Cp{v};
    GCM{i}.Ep = PRFn{1}.Ep{v};
    GCM{i}.M.pC = PRFn{1}.M.pC{v};
    GCM{i}.M.pE = PRFn{1}.M.pE{v};
    i = i + 1;
end
classic_BPA = spm_dcm_bpa(GCM);

%% Distributional BPA
included = dist_index(classic_noise == 2);
nincluded = length(included);
GCM = cell(nincluded,1);
i = 1;
for v = included
    GCM{i}.Cp = PRFn{2}.Cp{v};
    GCM{i}.Ep = PRFn{2}.Ep{v};
    GCM{i}.M.pC = PRFn{2}.M.pC{v};
    GCM{i}.M.pE = PRFn{2}.M.pE{v};
    i = i + 1;
end
dist_BPA = spm_dcm_bpa(GCM);
