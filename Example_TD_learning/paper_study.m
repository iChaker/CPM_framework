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
%% ======================= Noise levles ========================================
% Finding the noise variance for the BOLD signal to simulate a given SNR
%  levels. Previously estimates variance of simulated BOLD signal is 0.0083.
% var_sig = 0.0083;
% target_snr = [-20, -10, -2, 2, 10, 20];
% t_sd = sqrt(var_sig * 10 .^ (-target_snr / 10)); % Solving for standardeviation
% Result = [0.9110, 0.2881, 0.1147, 0.0724, 0.0288, 0.0091]
%% ======================= Paths for VOIs, PRFs, etc. ==========================
REDO = false;
voi_cfgs = {'one_tau_simulation.json', 'two_tau_simulation.json'};
voi_names = {'one_tau_simVOI.mat', 'two_tau_simVOI.mat'};
prf_names = {'onetau_PRFn.mat', 'twotau_PRFn.mat'};
prf_cfgs = {'one_tau_recovery_cfg.json', 'two_tau_recovery_cfg.json'};

%% ==================== Load VOIs and build them together ======================
VOI.xY.XYZmm = [];
VOI.xyz_def = [];
VOI.xY.y = [];

for cc = 1 : 2
    if ~isfile(fullfile('simulationfiles', voi_names{cc})) || REDO
        [VOIt, ~] = cpm_simulate_data(voi_cfgs{cc});
    else
        VOIt = load(fullfile('simulationfiles', voi_names{cc}));
        VOIt = VOIt.VOI;
    end
    VOIt.xY.XYZmm(3, :) = cc;
    VOI.xY.XYZmm = [VOI.xY.XYZmm, VOIt.xY.XYZmm];
    VOI.xyz_def = [VOI.xyz_def; VOIt.xyz_def];
    VOI.xY.y = [VOI.xY.y, VOIt.xY.y];
end

%% ============================= Make PRFs =====================================
PRFn = {};
for cc = 1 : 2

    if ~isfile(fullfile('simulationfiles', prf_names{cc}))  || REDO
        [prf_tmp] = cpm_simulation_prfs(VOI, prf_cfgs{cc});
    else
        prf_tmp = fullfile('simulationfiles', prf_names{cc});
    end

    PRF = load(prf_tmp);
    PRF = PRF.PRF;

    if ~isfile([prf_tmp(1 : end - 4) '_estimated' prf_tmp(end - 3 : end)]) || REDO
        PRFn{cc} = cpm_estimate(PRF, [], true);
        save([prf_tmp(1 : end - 4) '_estimated' prf_tmp(end - 3 : end)], 'PRFn')
    else
        tmp = load([prf_tmp(1 : end - 4) '_estimated' prf_tmp(end - 3 : end)]);
        PRFn{cc}  = tmp.PRFn{cc};
    end
end
%%======================== Extract data ========================================
model_idx = VOI.xY.XYZmm(3, :);
noise_idx = VOI.xY.XYZmm(2, :);
noises = unique(noise_idx); % Noise levels
nnoise = length(noises); % Number of noise levels
nmodels = sum(noise_idx == 1); % Number of different models
generating_params = VOI.xyz_def; % Parameters used for generation
simY = VOI.xY.y; % Simulated data
%%=================== Calculate actual SNR =====================================
base = 1; % Index, with 0 noise
signal_var = var(simY(:, noise_idx == base)); % Voxel wise, signal variance
snrs = zeros(size(simY, 2), 1); % pre allocate snrs
mean_snr = zeros(length(noises), 1); % pre - allocate 

for nidx = noises
    noise_var = var(simY(:, noise_idx == nidx) - simY(:, noise_idx == base));
    snrs(noise_idx == nidx) = signal_var ./ (noise_var + eps);
    mean_snr(nidx) = mean(10 * log10(snrs(noise_idx == nidx)));
end

snr_label = round(mean_snr, 2);
%% ===================== Extract F values ======================================
modelF = zeros(2, size(model_idx, 2));
for midx = 1 : 2
        modelF(midx, :) = PRFn{midx}.F;
end
%% ======================== Model comparison and diagonal idx ===================
model_idx_ext = model_idx;
for dd = find(model_idx == 2)
    if generating_params{dd}.mu_tauneg == generating_params{dd}.mu_taupos 
        model_idx_ext(dd) = 3;
    end
end

exceedanceP = {};
for kk = 1 : 3
    vox_idx = reshape(find(model_idx_ext == kk), [], nnoise);
    [~, exceedanceP{kk}, ~]  = cpm_model_comparison(modelF, vox_idx, (kk > 1) + 1); 
    % Weird hack, to use k as iterator, but only have two models to compare against.
end
%% ============== PLOT exceedance Probabilities ================================
if true 
ep_titles = {'Classical RL', 'Distributional RL', '\tau^- = \tau^+'};
fig1 = figure('Color', 'white', 'Units', 'pixels', 'Position', [0, 0, 800, 400]);
% Plot RFX exceedence probabilities
for rows = 1 : 3
    subplot(1, 3, rows)
    bar(noises(2:end), exceedanceP{rows}(:, 2:end)');
    legend({'Classic RL Model', 'Distributional RL Model'});
    title(ep_titles{rows});
    xlabel('SNR'); 
    ylabel('Exceedance Probability')
    ylim([0, 1]);
    xticklabels(snr_label(2 : end))
end

sgtitle('Model Recovery')
cpm_savefig(fig1, 'results/fig1_model_recovery.png')
end
%% ==================== Recovery Plots options =================================
ppd_samples = 100;
plot_dimentions = 'response';
plot_noise = 5;
pads = 40;
height_dims = [120, 225];
row_height = sum(height_dims);

onedim_t  = linspace(PRFn{1}.U(1).grid.tau(1), PRFn{1}.U(1).grid.tau(2), PRFn{1}.U(1).grid.tau(3));
%% =================== Plot Recovery classical ==================================
if true
voi_idx = reshape(find(model_idx == 1), [], nnoise); 
sets = reshape(voi_idx(:, plot_noise), 2, 2)';

% make axes
fig_x = 500; fig_y = 700;
fig2 = figure('Color', 'white', 'Units', 'pixels', 'Position',  [0, 0, fig_x + pads, fig_y + 2 * pads]);

for rows = 1 :2
    for cols = 1 : 2
        cl_axes{rows, cols} = axes('Units', 'pixels', 'Position', [0 + (fig_x / 2) * (cols - 1) + pads, ...
                                                                 fig_y - row_height * (rows - 1) - height_dims(1) - pads, ...
                                                                 fig_x / 2 - pads, ...
                                                                 height_dims(1) - 1.5 * pads]);
        hold on;
        gen_mu_tau = generating_params{sets(rows, cols)}.mu_tau;
        gen_sigma_tau = generating_params{sets(rows, cols)}.sigma_tau;

        plot_single_voxel(PRFn{1}, sets(rows, cols), {'tau'}, {[]}, {[]}, ppd_samples, plot_dimentions)
        y = normpdf(onedim_t, gen_mu_tau, gen_sigma_tau);
        y = y ./ sum(y);
        
        plot(onedim_t, y, 'LineWidth', 1.5);
        xlim(PRFn{1}.U(1).grid.tau(1 : 2));
        tmp_title = sprintf('\\mu_\\tau= %4.2f, \\sigma_\\tau= %4.2f', gen_mu_tau, gen_sigma_tau);
        title(tmp_title);
        xticklabels([])
        ylabel('Probability')

        dl_axes{rows, cols} = axes('Units', 'pixels',  'Position', [0 + (fig_x / 2) * (cols - 1) + pads, ...
                                                                 fig_y - rows * row_height - 0.25* pads, ...
                                                                 fig_x / 2 - pads, ...
                                                                 height_dims(2) - pads ]);
        hold on;
        plot_single_voxel(PRFn{2}, sets(rows, cols), {'tauneg', 'taupos'}, ...
                                            { [], []}, {[], []}, ppd_samples, plot_dimentions);

        xlim(PRFn{2}.U(1).grid.tauneg(1 : 2))
        ylim(PRFn{2}.U(1).grid.tauneg(1 : 2))
        xlabel('\tau^-')
        ylabel('\tau^+')
        ci = [gen_mu_tau - 2 * gen_sigma_tau, gen_mu_tau + 2 * gen_sigma_tau] ;
        plot(ci, ci, 'color', 'yellow', 'LineWidth', 1.5);
        plot(gen_mu_tau, gen_mu_tau, 'o', 'color', 'white', 'MarkerFaceColor','white')
    end

end

for rows = 1 :2
    for cols = 1 :2
        cl_axes{rows, cols}.Position(2) = cl_axes{rows, cols}.Position(2) + 2 * pads;
        dl_axes{rows, cols}.Position(2) = dl_axes{rows, cols}.Position(2) + 2 * pads;
    end
end

sgt = sgtitle({'Parameter Recovery: Classic RL', ['SNR:', num2str(snr_label(plot_noise))]});
cpm_savefig(fig2, 'results/fig2_parameter_recovery_classic_rl.png')
end
%% ============================ Plot recovery Distributional ===================
if true 
fig_x = 800 ; fig_y = 4.3 * row_height;
pads = 35;
fig3 = figure('Color', 'white', 'Units', 'pixels', 'Position',  [0, 0, fig_x + pads, fig_y + 6 * pads]);
axis('off')
normalize_vec = [fig_x + pads, fig_y + 6 * pads, fig_x + pads, fig_y + 6 * pads];

voi_idx = reshape(find(model_idx == 2), [], nnoise); 
sets = reshape(voi_idx(:, plot_noise), 4, 4)';

for cols = 1 :4
    for rows = 1 : 4
        hold on
        cl_axes{rows, cols} = axes('Units', 'normalized',  'Position', [0 + (fig_x / 4) * (cols - 1) + pads, ...
                                                                 fig_y - row_height * (rows - 1) - height_dims(1) - 1.5 * pads, ...
                                                                 fig_x / 4 - pads, ...
                                                                 height_dims(1) - 1.75 * pads] ./ normalize_vec);

        plot_single_voxel(PRFn{1}, sets(rows, cols), {'tau'}, { []}, {[]}, ppd_samples, plot_dimentions)
        
        xlim(PRFn{1}.U(1).grid.tau(1 : 2));
        
        gen_mu_taupos = generating_params{sets(rows, cols)}.mu_taupos;
        gen_mu_tauneg = generating_params{sets(rows, cols)}.mu_tauneg;
        gen_sigma_taupos = generating_params{sets(rows, cols)}.sigma_taupos;
        gen_sigma_tauneg = generating_params{sets(rows, cols)}.sigma_tauneg;

        tmp_title1 = sprintf('\\mu_{\\tau^+}= %4.2f, \\sigma_{\\tau^+}= %4.2f', ...
                                        gen_mu_taupos, gen_sigma_taupos);
       
        tmp_title2 = sprintf('\\mu_{\\tau^-}= %4.2f, \\sigma_{\\tau^-}= %4.2f', ...
                                        gen_mu_tauneg, gen_sigma_tauneg);        
        tmpt =  title({tmp_title1; tmp_title2});
        tmpt.FontSize = 6;
        axis('off')
        
        dl_axes{rows, cols} = axes('Units', 'normalized',  'Position', [0 + (fig_x / 4) * (cols - 1) + pads, ...
                                                                 fig_y - rows * row_height - 0.25* pads, ...
                                                                 fig_x / 4 - pads, ...
                                                                 height_dims(2) - 0.25 * pads ] ./ normalize_vec);
        
        plot_single_voxel(PRFn{2}, sets(rows, cols), {'tauneg', 'taupos'}, ...
                                    { [], []}, {[], []}, ppd_samples, plot_dimentions)
        
        scatter(gen_mu_tauneg, gen_mu_taupos, ...
                    'filled', 'MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerFaceColor',[1 1 1],  'LineWidth',1.0)
        
        tp = -pi : 0.01 : pi;
        x_post = gen_mu_tauneg + 2 * gen_sigma_tauneg .* cos(tp);
        y_post = gen_mu_taupos + 2 * gen_sigma_taupos .* sin(tp);
        plot(x_post, y_post)
        xlabel('\tau^-')
        ylabel('\tau^+')
        
        xlim(PRFn{2}.U(1).grid.tauneg(1 : 2))
        ylim(PRFn{2}.U(1).grid.taupos(1 : 2))
    
    end
end

pads_move =  [120; 80; 40; 0] - 15; 

for rows = 1 :4
    for cols = 1 :4
%         cl_axes{rows, cols}.Position(2) = cl_axes{rows, cols}.Position(2) + (pads_move(rows, 1) + 2.6 * pads) / normalize_vec(2);
%         dl_axes{rows, cols}.Position(2) = dl_axes{rows, cols}.Position(2) + (pads_move(rows, 1) + 1.75 * pads) / normalize_vec(2);
        cl_axes{rows, cols}.Position(2) = cl_axes{rows, cols}.Position(2) + (40 + pads_move(rows, 1) ) / normalize_vec(2);
        dl_axes{rows, cols}.Position(2) = dl_axes{rows, cols}.Position(2) + (pads_move(rows, 1) ) / normalize_vec(2);

    end
end

sgtitle({'Parameter Recovery: Distributional RL', ['SNR:', num2str(snr_label(plot_noise))]});

cpm_savefig(fig3, 'results/fig3_parameter_recovery_dist_rl.png')
end
%% ========================== Estimation error ================================
% Classic
[trues, preds] = deal({}, {});
for kk = 1 : 2
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
mses{1} = error_fun(reshape(trues{1}, 6, [], nnoise), reshape(preds{1}, 6, [], nnoise));
mses{2} = error_fun(reshape(trues{2}, 8, [], nnoise), reshape(preds{2}, 8, [], nnoise));
%%
fig4 = figure('Position', [0, 0, 1200, 600]);
sub_titles = {'Classical RL', 'Distributional RL'};
for nc = 1 : 2
    subplot(1, 2, nc)
    labels = fieldnames(cpm_get_true_parameters(PRFn{nc}, 1));
    labels = strrep(labels, '_', ' ');
    h = heatmap(round(mses{nc}, 4), 'YDisplayLabels', labels, 'XDisplayLabels', snr_label, 'XLabel', 'SNR', 'YLabel', 'Parameter');
    title(sub_titles{nc})
end

sgtitle('Parameter Recovery: RMSE')

cpm_savefig(fig4, 'results/fig4_rmse_parameter_recovery.png')

%% Recover Tau* 
alphas_pred = zeros(length(generating_params), 2);
alphas_true = zeros(length(generating_params), 2);

for nn = 1 : length(generating_params)
   
        tmp_pred = spm_vec(cpm_get_true_parameters(PRFn{2}, nn));
        alphas_pred(nn, :) = cpm_logit_inv(tmp_pred(1 : 2));
        try
        tmp_true = [generating_params{nn}.mu_tauneg,  generating_params{nn}.mu_taupos];
        catch
        tmp_true = [generating_params{nn}.mu_tau,  generating_params{nn}.mu_tau];
        end
        alphas_true(nn, :) = cpm_logit_inv(tmp_true(1 : 2));
end
%
laterality_true = alphas_true(:, 1) ./ sum(alphas_true, 2);
laterality_pred = alphas_pred(:, 1) ./ sum(alphas_pred, 2);
laterality_error = laterality_true - laterality_pred;
%
laterality_true = [reshape(laterality_true(model_idx==1, :), [], 7); reshape(laterality_true(model_idx==2, :),  [], 7)];
laterality_pred = [reshape(laterality_pred(model_idx==1, :), [], 7); reshape(laterality_pred(model_idx==2, :),  [], 7)];
laterality_error = [reshape(laterality_error(model_idx==1, :), [], 7); reshape(laterality_error(model_idx==2, :),  [], 7)];
%
noise_mat = ones(size(laterality_error)) .* [1 : 7];
model_index_res = [reshape(model_idx(model_idx ==1), [], 7); reshape(model_idx(model_idx ==2), [], 7)];
%%
fig5 = figure('Position', [0, 0, 1200, 400]); 
subplot(1, 2, 1)
for bl = unique(laterality_true)'
    scatter(reshape(noise_mat(laterality_true(:, 1)==bl, 2:end), [], 1), reshape(abs(laterality_error(laterality_true(:, 1)==bl, 2:end )), [], 1), ...
        30 * reshape(model_index_res(laterality_true(:, 1)==bl, 2:end), [], 1), 'filled')
    hold on
end

xlabel('SNR')
ylabel('Error')

xticklabels(mean_snr(2 : end))
legend({'\tau^*=0.25', '\tau^*=0.50', '\tau^*=0.75'})

title('Absolute error of learning rate asymmetry')
subplot(1, 2, 2)
plot(laterality_true(:, plot_noise), laterality_true(:, plot_noise),  'Color', [0, 0, 0] + 0.5)
hold on

markers = {'o', 'square'};
for nc = 1 : 2
    scatter(laterality_true(model_index_res(:, plot_noise) == nc, plot_noise), ...
                laterality_pred(model_index_res(:, plot_noise) == nc, plot_noise), [], 'filled', 'Marker', markers{nc})
end
xlabel('True values')
ylabel('Estimated values')
xlim([0.2, 0.8]);
ylim([0.2, 0.8]);
title(sprintf('Estimation at SNR %4.2f', snr_label(plot_noise)))
legend({'', 'Classic RL', 'Distributional RL'}, 'Location','northwest')

cpm_savefig(fig5, 'results/fig5_learning_assymetry_tau.png')

%% Classic BPA
fig5 = figure('Position', [0, 0, 1900, 700]); 
prf_names = {'Classic RL', 'Distributional RL'};
tiledlayout(2, 6, 'TileSpacing', 'tight', 'Padding', 'compact')

for nc = 1 : 2
    for nn = 2 : 7
        nexttile()
        included =  find(model_idx == nc & noise_idx == nn);
        nincluded = length(included);
        GCM = cell(nincluded,1);
        i = 1;
        for v = included
            GCM{i}.Cp = PRFn{nc}.Cp{v};
            GCM{i}.Ep = PRFn{nc}.Ep{v};
            GCM{i}.M.pC = PRFn{nc}.M.pC{v};
            GCM{i}.M.pE = PRFn{nc}.M.pE{v};
            i = i + 1;
        end
        classic_BPA = spm_dcm_bpa(GCM);
        
        cl_labels = fieldnames(PRFn{nc}.M.pE{1});
        cl_labels = strrep(cl_labels, '_', ' ');
        tmp_mat = VBA_cov2corr(classic_BPA.Cp);
        idx = tril(tmp_mat);
        tmp_mat(~idx) = nan;
        t = heatmap(tmp_mat, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', ...
                       " ", 'ColorbarVisible', 'off', 'Colormap', colormap('parula'), 'XDisplayLabels', cl_labels, ...
                        'YDisplayLabels', cl_labels, 'FontSize', 10 - 2 *  nc, 'CellLabelFormat', '%0.2f');
        t.InnerPosition = [0, 0, 1, 1];

        if nc == 1
            title(sprintf('SNR %4.2f', snr_label(nn)));
        end
        if nn == 2
            ylabel(prf_names{nc});
        end

    end
end

sgtitle('Posterior Correlation after BPA')

cpm_savefig(fig5, 'results/fig5_posterior_correlation_bpa.png')

%% Predicted Y
predYs = {};
for nc =  1 : 2
    predYs{nc} = zeros(size(PRFn{nc}.Y.y));
    for vx = 1 : nnoise * nmodels
        predYs{nc}(:, vx) = spm_prf_response(PRFn{nc}.Ep{vx}, PRFn{nc}.M, PRFn{nc}.U);
    end
end

%% Simulated HRF accuracy RMSE / R2

fits_r2 = zeros(2, nmodels * nnoise);
fits_mse = zeros(2, nmodels * nnoise);

for nc = 1 : 2
    for vx = 1 : nnoise * nmodels
        fits_r2(nc, vx) = VBA_r2(predYs{nc}(:, vx), simY(:, vx));
        fits_mse(nc, vx) = sqrt(mean((simY(:, vx) - predYs{nc}(:, vx)).^2));
    end
end

%% Plot classical RL

fig6 = figure('Position', [0, 0, 1200, 500]);

for nc = 1 : 2
    subplot(2, 2, nc)
    tmp_r2 = (fits_r2(:, model_idx==nc));
    tmp_r2 =  [tmp_r2(1, :), tmp_r2(2, :)];
    model_g =[zeros(1, sum(model_idx==nc)), zeros(1, sum(model_idx==nc)) + 1];
    noise_g = [noise_idx(model_idx==nc), noise_idx(model_idx==nc)];

boxchart(noise_g(noise_g > 1) - 2, tmp_r2(noise_g > 1), 'GroupByColor', model_g(noise_g > 1))
legend({'Classic RL', 'Distrubtional RL'}, 'Location', 'SouthEast')
xticks(0 : 5)
xticklabels(snr_label(2 : 7))
xlabel('SNR')
ylabel('R^2')
title('Generative Process:', prf_names{nc})
end

subplot(2, 2, 3)

% Fits of interest: At Model SNR largest difference
tmp_idx =  find(model_idx == 2 & noise_idx == plot_noise);
[~, max_diff] = max(diff(fits_r2(:, tmp_idx)));
max_diff = 6;

hold on
t = (0 : size(simY, 1) - 1) * 0.592;
lh  = plot(t, simY(:, tmp_idx(max_diff)), 'Color', 'black');
lh.Color(4) = 0.1;
lh2 = plot(t, predYs{1}(:, tmp_idx(max_diff)), 'LineStyle',':', 'LineWidth', 1.5);
lh3 = plot(t, predYs{2}(:, tmp_idx(max_diff)), 'Color', 'black');

lh2.Color  = [ 0, 0.4470, 0.7410, 1];
lh3.Color = [0.8500    0.3250    0.0980, 1];

xlabel('time')
ylabel('signal, a.u.')
title('Simulated and Recovered Signals', sprintf('SNR %4.2f', snr_label(plot_noise)));
legend({'Simulated', 'Classic RL', 'Distributional RL'}, 'Location', 'SouthEast')

subplot(2, 2, 4)
hold on

sc1 = scatter(simY(:, tmp_idx(max_diff)), predYs{1}(:, tmp_idx(max_diff)), 0.8, 'filled');
sc2 = scatter(simY(:, tmp_idx(max_diff)), predYs{2}(:, tmp_idx(max_diff)), 0.8, 'filled');
l2  = lsline;

l2(1).Color = sc1.CData;
l2(2).Color = sc2.CData;

xlim([min(simY(:, tmp_idx(max_diff))), max(simY(:, tmp_idx(max_diff)))])
ylim([ min(predYs{2}(:, tmp_idx(max_diff))), max(predYs{2}(:, tmp_idx(max_diff)))])

ylabel('Simulated')
xlabel('Recovered')
rf = refline(0.5, 0);
rf.Color = [0.2, 0.2, 0.2, 0.4];

title_r2 = round(fits_r2(:, tmp_idx(max_diff)), 2);

legend({['Classic RL (R^2 = ', num2str(round(title_r2(1), 2)) ')'], ...
             ['Distributional RL (R^2 = ', num2str(round(title_r2(2), 2)) ')']}, 'Location', 'SouthEast')

disp(generating_params{tmp_idx(max_diff)})

title('Simulated vs Recovered')

sgtitle('Classical Model fit');

cpm_savefig(fig6, 'results/fig6_signal_recovery_r2.png')
