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


%%
if true
post_samples = 400;
%
noise_level = 2;
unq_models = find(model_idx == 1 & noise_idx == noise_level);

% make axes
pad = 40;
fig_x = 500; fig_y = 500; 
fig2 = figure('Color', 'white', 'Units', 'pixels', 'Position', ...
                    [0, 0, fig_x + pad, fig_y + 2 * pad]);
axis('off')
hold on;


tiledlayout(4, 4,"TileSpacing","tight", "Padding", "tight");


for ii = 1 : unq_models

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
        
  
        scatter(x, y,'filled', 'MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerFaceColor',[1 1 1],  'LineWidth',1.0)
        
        tp = -pi:0.01:pi;
        x_post = x + 2 *  x_sig .* cos(tp);
        y_post =y +  2 * y_sig .* sin(tp);
        plot(x_post, y_post)

end

end

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
dist_error_alpha = zeros(length(spm_vec(cpm_get_true_parameters(PRFn{2}, 1))) + 2, length(dist_index));

for cidx = dist_index
    tmp_gen_vals = spm_vec(cpm_get_true_parameters(PRFn{2}.M.pE{cidx}, PRFn{2}.M, PRFn{2}.U));
    tmp_gen_vals(1 : 2) = [generating_params{cidx}.mu_tauneg,  generating_params{cidx}.mu_taupos];
    tmp_gen_vals(3 : 4) = [generating_params{cidx}.sigma_tauneg,  generating_params{cidx}.sigma_taupos];
    dist_error(:, cidx) = spm_vec(cpm_get_true_parameters(PRFn{2}, cidx)) - tmp_gen_vals;
    
    tmp_gen_alpha = zeros(10, 1);
    tmp_gen_alpha(1 : 2) = cpm_logit_inv([generating_params{cidx}.mu_tauneg,  generating_params{cidx}.mu_taupos]);
    tmp_gen_alpha(3 : 4) = cpm_logit_inv([generating_params{cidx}.mu_tauneg - generating_params{cidx}.sigma_tauneg, ...
                                                                  generating_params{cidx}.mu_tauneg + generating_params{cidx}.sigma_tauneg]);
    tmp_gen_alpha(5 : 6) = cpm_logit_inv([generating_params{cidx}.mu_taupos - generating_params{cidx}.sigma_taupos, ...
                                                                  generating_params{cidx}.mu_taupos + generating_params{cidx}.sigma_taupos]);
    tmp_gen_alpha(7 : 10) = tmp_gen_vals(5 : end);
        
    tmp_inv_cpm = spm_vec(cpm_get_true_parameters(PRFn{2}, cidx));
    tmp_inv_alpha = zeros(10, 1);
    tmp_inv_alpha(1 : 2) = cpm_logit_inv(tmp_inv_cpm(1 : 2));
    tmp_inv_alpha(3 : 4) = cpm_logit_inv([tmp_inv_cpm(1) -  tmp_inv_cpm(3), tmp_inv_cpm(1) +  tmp_inv_cpm(3)]);
    tmp_inv_alpha(5 : 6) = cpm_logit_inv([tmp_inv_cpm(2) -  tmp_inv_cpm(4), tmp_inv_cpm(2) +  tmp_inv_cpm(4)]);
    tmp_inv_alpha(7 : 10) = tmp_inv_cpm(5 : end);
    dist_error_alpha(:, cidx)  = tmp_inv_alpha - tmp_gen_alpha;
end

%% Create some plot for errors?!


%% Summarize errors for over noise
dist_mse = zeros(size(dist_error, 1), nnoise);
dist_alpha_mse = zeros(size(dist_error_alpha, 1), nnoise);

for nidx =  1 : nnoise
    dist_mse(:, nidx) = mean(dist_error(:, dist_noise == nidx).^2, 2);
    dist_alpha_mse(:, nidx) = mean(dist_error_alpha(:, dist_noise == nidx).^2, 2);
end

classic_mse = zeros(size(classic_error, 1), nnoise);
for nidx =  1 : nnoise
    classic_mse(:, nidx) = mean(classic_error(:, classic_noise == nidx).^2, 2);
end
%%
figure;
h = heatmap(sqrt(dist_alpha_mse));
labels = fieldnames(cpm_get_true_parameters(PRFn{2}, 1));
%h.YDisplayLabels = labels;


%% Recover Tau* 
alphas_pred = zeros(num_models, 2, nnoise);
alphas_true = zeros(num_models, 2, nnoise);
color_idx = zeros(num_models, nnoise);
for nn = 1 : nnoise
    tmp_idx = find(noise_idx == nn);
    cc = 1;
    for jj = tmp_idx
        tmp_pred = spm_vec(cpm_get_true_parameters(PRFn{2}, jj));
        alphas_pred(cc, :, nn) = cpm_logit_inv(tmp_pred(1 : 2));
        try
        tmp_true = [generating_params{jj}.mu_tauneg,  generating_params{jj}.mu_taupos];
        color_idx(cc, nn) = 1;
        catch
        tmp_true = [generating_params{jj}.mu_tau,  generating_params{jj}.mu_tau];
        color_idx(cc, nn) = 2;
        end
        alphas_true(cc, :, nn) = cpm_logit_inv(tmp_true(1 : 2));
        cc = cc + 1;
    end
end

laterality_true = squeeze(alphas_true(:, 1, :)) ./ squeeze(sum(alphas_true, 2));
laterality_pred = squeeze(alphas_pred(:, 1, :)) ./ squeeze(sum(alphas_pred, 2));
laterality_error = laterality_true - laterality_pred;
%% make table
noise_mat = ones(size(laterality_error)) .* [1 : 7];
%%

figure; 
subplot(2, 1, 1)
for bl = unique(laterality_true)'
    scatter(reshape(noise_mat(laterality_true(:, 1)==bl, 2:end), [], 1), reshape(sqrt(laterality_error(laterality_true(:, 1)==bl, 2:end ).^2), [], 1), 30 * reshape(color_idx(laterality_true(:, 1)==bl, 2:end), [], 1), 'filled')
    hold on
end
xticklabels(mean_snr(2 : end))
legend({'\tau^*=0.25', '\tau^*=0.50', '\tau^*=0.75'})

subplot(2, 1, 2)
plot(laterality_true(:, 5), laterality_true(:, 5),  'Color', [0, 0, 0] + 0.5)
hold on

scatter(laterality_true(:,5), laterality_pred(:, 5), [], color_idx(:,5) ./ 2 , 'filled')
%% Classic BPA
figure('Position', [0, 0, 1800, 600]); 

tiledlayout(2, 6, 'TileSpacing', 'tight', 'Padding', 'tight')
for nn = 2 : 7
nexttile()
included = classic_index(classic_noise == nn);
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

cl_labels = fieldnames(PRFn{1}.M.pE{1});
tmp_mat = VBA_cov2corr(classic_BPA.Cp);
idx = tril(tmp_mat);
tmp_mat(~idx) = nan;
t = heatmap(tmp_mat, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', ...
               " ", 'ColorbarVisible', 'off', 'Colormap', colormap('parula'), 'XDisplayLabels', cl_labels, ...
                'YDisplayLabels', cl_labels, 'FontSize', 8, 'CellLabelFormat', '%0.2f');

t.InnerPosition = [0, 0, 1, 1];
%% Distributional BPA
end

for nn = 2 : 7
nexttile()
included = dist_index(classic_noise == nn);
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

cl_labels = fieldnames(PRFn{2}.M.pE{1});

dist_BPA = spm_dcm_bpa(GCM);
tmp_mat = VBA_cov2corr(dist_BPA.Cp);
idx = tril(tmp_mat);
tmp_mat(~idx) = nan;

t = heatmap(tmp_mat, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', ...
               " ", 'ColorbarVisible', 'off', 'Colormap', colormap('parula'), 'XDisplayLabels', cl_labels, ...
                'YDisplayLabels', cl_labels, 'FontSize', 6, 'CellLabelFormat', '%0.2f');

t.InnerPosition = [0, 0, 1, 1];
end