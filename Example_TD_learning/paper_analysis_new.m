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

if ~isfile('simulationfiles/tau_simVOI.mat')  || REDO
	[VOI, ~] = cpm_simulate_data('tau_simulation_cfg.json');
else
    VOI = load('simulationfiles/tau_simVOI.mat');
    VOI = VOI.VOI;
end

%
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
        tmp = load([prf_path(1 : end - 4) '_estimated' prf_path(end-3 : end)]);
        PRFn{cc}  =tmp.PRFn{cc};
    end

    cc = cc + 1;
end

%%

% TODOs:
% Create table on differences between recovered values and simulation values for
% non-point and point estimates.

%%
VOI = load(fullfile('simulationfiles/tau_simVOI.mat'));
VOI = VOI.VOI;

gen_idx = VOI.xY.XYZmm;
model_idx = gen_idx(1, :);
noise_idx = gen_idx(2, :);
noises = unique(gen_idx(2, :));

num_models = sum(noise_idx == 1);

point_idx = find(gen_idx == 1);
grid_idx = find(gen_idx ~= 1);

simY = VOI.xY.y;
%% SNRs
base = 1;

signal_var = var(simY(:, noise_idx == base));

snrs = zeros(size(simY, 2), 1);
mean_snr = zeros(length(noises), 1);

for nidx = noises
    noise_var = var(simY(:, noise_idx == nidx) - simY(:, noise_idx == base));
    snrs(noise_idx == nidx) = signal_var ./ (noise_var + eps);
    mean_snr(nidx) = mean(log10(snrs(noise_idx == nidx)));
end

%%
voi_sims = {zeros(size(simY)), zeros(size(simY))};
for pidx = 1 : size(gen_idx,2)
    voi_sims{1}(:, pidx) = spm_prf_response(PRFn{1}.Ep{pidx}, PRFn{1}.M, PRFn{1}.U);
    voi_sims{2}(:, pidx) = spm_prf_response(PRFn{2}.Ep{pidx}, PRFn{2}.M, PRFn{2}.U);
end
%%
generating_params = repmat(VOI.xyz_def, length(noises), 1);
sigma_params = nan(2, size(gen_idx, 2));
mu_params = zeros(2, size(gen_idx, 2));

for gidx = 1 : size(gen_idx, 2)
    try generating_params{gidx}.mu_tauneg;
        sigma_params(1, gidx) = generating_params{gidx}.sigma_tauneg;
        sigma_params(2, gidx) = generating_params{gidx}.sigma_taupos;
        mu_params(1, gidx) = generating_params{gidx}.mu_tauneg;
        mu_params(2, gidx) = generating_params{gidx}.mu_taupos;
    catch
        mu_params(1, gidx) = generating_params{gidx}.alphaneg;
        mu_params(2, gidx) = generating_params{gidx}.alphapos;
    end
end
%%
sigma_combs = zeros(1, size(sigma_params, 2));
unq_sigmas = unique(sigma_params(:));
unq_sigmas = unq_sigmas(~isnan(unq_sigmas)) ;

for kk = 1 : length(unq_sigmas)
    sigma_combs(sigma_params(1, :) == unq_sigmas(kk)) = sigma_combs(sigma_params(1, :) == unq_sigmas(kk)) +  kk *  10;
    sigma_combs(sigma_params(2, :) == unq_sigmas(kk)) = sigma_combs(sigma_params(2, :) == unq_sigmas(kk)) +  kk;
end

% Important (but later)
% Label plots

%% Massive Plotting Recovery

unq_sig_combs = unique(sigma_combs);
unq_sig_combs = unq_sig_combs(unq_sig_combs ~= 0);
forward_fun = @(x, xmin, xmax) ((xmax-xmin) .* spm_Ncdf(x,0,1)) + xmin;

for midx = 1 : 2
    for nn =  noises 
       for sig = unq_sig_combs
            
           figure; 
           tiles = tiledlayout(3, 3); 
        
           points = find(noise_idx == nn & sigma_combs == sig);
        for kk  =points
                nexttile()
              
                if midx == 2
                    plot_single_voxel(PRFn{midx}, kk, {'tauneg', 'taupos'}, { [], []}, {[], []}, 1000, true)
                else
                    plot_single_voxel(PRFn{midx}, kk, {'tau', 'tau'}, { [], []}, {[], []}, 1000, true)
                end
              % params
              x = generating_params{kk}.mu_tauneg;
              x_sig = forward_fun(generating_params{kk}.sigma_tauneg, PRFn{2}.M.cpm.sigma.tauneg(1), PRFn{2}.M.cpm.sigma.tauneg(2));
              y = generating_params{kk}.mu_taupos;
              y_sig =forward_fun(generating_params{kk}.sigma_taupos, PRFn{2}.M.cpm.sigma.taupos(1), PRFn{2}.M.cpm.sigma.taupos(2));

              scatter(x, y,'filled', 'MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerFaceColor',[1 1 1],  'LineWidth',1.5)
              
              if midx == 2
                  xlim(PRFn{2}.U(1).grid.tauneg(1 : 2))
                  ylim(PRFn{2}.U(1).grid.tauneg(1 : 2))
              else
                  xlim(PRFn{1}.U(1).grid.tau(1 : 2))
                  ylim(PRFn{1}.U(1).grid.tau(1 : 2))
              end

                t = -pi:0.01:pi;
                x_post = x + 2 *  x_sig .* cos(t);
                y_post =y +  2 * y_sig .* sin(t);
                plot(x_post, y_post)
                xlabel('\tau^-')
                ylabel('\tau^+')
                
                title(sprintf('Generating: \\mu_{\\tau^-} = %4.2f, \\mu_{\\tau^+} = %4.2f ', x, y))
        end
        sgtitle(sprintf('Recovery, with SNR: %4.2f\n \\sigma_{\\tau^-} = %4.2f, \\sigma_{\\tau^+} = %4.2f', mean(log10(snrs(points))), x_sig, y_sig))
           
        end
    end
end

%% Plotting BOLD
for midx = 2 %1 : 2
       for sig = 11 %unq_sig_combs
            figure;
            points = find(sigma_combs == sig & noise_idx == 1);
            for kk  =points
                
                nexttile()
                %plot(VOI.xY.y(:, kk))
                hold on;
                for nn = 2 %noises
                    tmp_idx = kk + (nn - 1) * num_models;
                    plot(VOI.xY.y(:, kk) - voi_sims{midx}(:, tmp_idx));
                end
                hold off
            end
    end
end

%% Parameters
n_params = length(spm_vec(PRFn{2}.Ep{1}));

original_parameters = nan(n_params, num_models);
recovered_parameters = nan(n_params, num_models);

static_params  = spm_vec(cpm_get_true_parameters(PRFn{2}.M.pE{1}, ...
                                            PRFn{2}.M, PRFn{2}.U));
static_params = static_params(end-3 : end);

for pidx = 1 : size(gen_idx, 2)
    
    if isfield(generating_params{pidx}, 'mu_tauneg')
        sigma_neg = forward_fun(generating_params{pidx}.sigma_tauneg, ...
                                                  PRFn{2}.M.cpm.sigma.tauneg(1), PRFn{2}.M.cpm.sigma.tauneg(2));
        sigma_pos = forward_fun(generating_params{pidx}.sigma_taupos, ...
                                                  PRFn{2}.M.cpm.sigma.taupos(1), PRFn{2}.M.cpm.sigma.taupos(2));
         mu_neg = generating_params{pidx}.mu_tauneg;
         mu_pos = generating_params{pidx}.mu_taupos;
         
         original_parameters(:, pidx) = [mu_neg; mu_pos; sigma_neg; sigma_pos; static_params];
         
         recovered_parameters(:, pidx) = spm_vec(cpm_get_true_parameters(PRFn{2}, pidx));
    else
        sigma_neg = forward_fun(-8.1259, ...
                                          PRFn{2}.M.cpm.sigma.tauneg(1), PRFn{2}.M.cpm.sigma.tauneg(2));
        sigma_pos = forward_fun(-8.1259, ...
                                                  PRFn{2}.M.cpm.sigma.taupos(1), PRFn{2}.M.cpm.sigma.taupos(2));
                                              
        alpha_neg = forward_fun(cpm_logit(generating_params{pidx}.alphaneg),-8.1259, 8.1259);
        alpha_pos = forward_fun(cpm_logit(generating_params{pidx}.alphapos),-8.1259, 8.1259);
        
         original_parameters(:, pidx) = [alpha_neg; alpha_pos; sigma_neg; sigma_pos; static_params];
         recovered_parameters(:, pidx) = spm_vec(cpm_get_true_parameters(PRFn{2}, pidx));
    end
end
%%
original_parameters = reshape(original_parameters, 8, num_models, []);
recovered_parameters = reshape(recovered_parameters, 8, num_models, []);
%%
parameter_errors = sqrt((original_parameters - recovered_parameters).^2);

figure;
tiledlayout(2, 4)

for ii = 1 : 8
nexttile()
hold on
plot(squeeze(parameter_errors(ii, :, :)))
    
end
%%
noise = 1;
nnoise = 5;
nvoxels = length(gen_idx) / nnoise;



%%

figure; 
tiles = tiledlayout(3, 3); 
inverse_fun = @(x, xmin, xmax) (norminv((x - xmin) ./ (xmax - xmin)));

points = find(gen_idx(1, :) == 1 & gen_idx(2, :) == noise);

for kk  =points
    nexttile()
    %plot_single_voxel(PRFn{2}, kk, {'tauneg', 'taupos'}, { @cpm_logit_inv, @cpm_logit_inv}, { 'alphaneg', 'alphapos'}, 1000, true)
    plot_single_voxel(PRFn{2}, kk, {'tauneg', 'taupos'}, { [], []}, {[], []}, 1000, true)
    % params
    x = inverse_fun(generating_params{kk - (noise - 1) * num_models}.alphaneg, -eps, 1 + eps);
    y = inverse_fun(generating_params{kk - (noise - 1) * num_models}.alphapos, -eps, 1 + eps); 
   
    scatter(x, y,'filled', 'MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerFaceColor',[1 1 1],  'LineWidth',1.5)
    xlim(PRFn{2}.U(1).grid.tauneg(1 : 2))
    ylim(PRFn{2}.U(1).grid.tauneg(1 : 2))

end
%%

reducedF = zeros(2, length(gen_idx));
% 
% reduced_val = -8.125;
% 
% reductions_re = {[1, 3, 4    ], [4    ], [1, 3]};
% reductions_rc = {[1, 2, 3, 4], [2, 4], [1, 3]};
% reduction_names = {'null', 'no-utility', 'no-learning', 'full'};
% 
for rd  = 1 :  length(PRFn) % length(reductions_re) + 1
    for vidx = 1 : length(gen_idx)

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
nparams = num_models;

[logBF, exceedenceP, comparisons] = deal({}, {}, {});
% Param idx which generated 0 models:
nullidx =1;
% small trick to expand initial index to noise 
null_voxels =  repmat(nullidx, 1 , nnoise) +   (0 : nnoise - 1) * nparams;
[logBF{1}, exceedenceP{1}, comparisons{1}]  = cpm_model_comparison(reducedF, null_voxels, 1);

% Model comparisons for no utility model:
% Param idx for no-utility - including no learning, i.e. 0 model is included
for kk = 2 : 3
    if kk == 2
        idx = 2;
    elseif kk == 3
        idx = 1;
    end
    no_idx = find(gen_idx(1, :) == kk & gen_idx(2, :) == 1)';
    no_vox =  repmat(no_idx, 1 , nnoise) +   (0 : nnoise - 1) * nparams;
    [logBF{kk}, exceedenceP{kk}, comparisons{kk}]  = cpm_model_comparison(reducedF, no_vox, idx);
end



%%

 points = gen_idx(1, gen_idx(2,:) == 1);
 points = find(points ~=1);
 params = VOI.xyz_def(points);

models_grid = cell(3, 3);

for ii = 1 : length(params)
    
    if params{ii}.sigma_tauneg == params{ii}.sigma_taupos
    if params{ii}.mu_tauneg == -4 
        rw_idx = 1;
    elseif params{ii}.mu_tauneg == 0
        rw_idx = 2;
    elseif    params{ii}.mu_tauneg == 4 
        rw_idx = 3;
    end  
    
        if params{ii}.mu_taupos == -4 
        cl_idx = 1;
    elseif params{ii}.mu_taupos == 0
        cl_idx = 2;
    elseif    params{ii}.mu_taupos == 4 
        cl_idx = 3;
        end  
    
        models_grid{rw_idx, cl_idx} = [models_grid{rw_idx, cl_idx}, points(ii)];
    end
end
%%
figure; 
tiles = tiledlayout(3, 3); 

vba_options.DisplayWin = 0;

for c = 1 : 3
    for r = 1 : 3
                 nexttile()

         [~, o] = VBA_groupBMC(reducedF(:, models_grid{r, c}), vba_options);
         bar(o.ep);
    end
end

%%
% noutilityidx = find(mid(1, :) == 2 & mid(2, :) == 1)';
% nou_voxels =  repmat(noutilityidx, 1 , nnoise) +   (0 : nnoise - 1) * nparams;
% [logBF{2}, exceedenceP{2}, comparisons{2}]  = cpm_model_comparison(reducedF, nou_voxels, 2);
% % Model comparisons for no learning model:
% % Param idx for no-learning - including no utility, i.e. 0 model is included
% nolearningidx = find(mid(1, :) == 3 & mid(2, :) == 1)';
% nol_voxels =  repmat(nolearningidx, 1 , nnoise) +   (0 : nnoise - 1) * nparams;
% [logBF{3}, exceedenceP{3}, comparisons{3}]  = cpm_model_comparison(reducedF, nol_voxels, 3);
% 
% % Model comparisons for full model:
% % Param idx for full models excluding all 0 models
% fullidx = find(mid(1, : ) == 4 & mid(2, :) == 1)';
% full_voxels =  repmat(fullidx, 1 , nnoise) +   (0 : nnoise - 1) * nparams;
% [logBF{4}, exceedenceP{4}, comparisons{4}]  = cpm_model_comparison(reducedF, full_voxels, 4);

%% Result plotting
%fig2 = figure('Color', 'none', 'Units', 'pixels', 'Position', [0, 0, 1600, 1200]);
fig2 = figure;
ep_titles = {'Null-Model', 'two taus', 'one tau'};

% pre specify colors
colors = [[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; ...
                [0.9290, 0.6940, 0.1250]; [0.4940, 0.1840, 0.5560]];
markers = {'o', 'x', 's', 'p'};

voi_noise = 1 : 5;
% Plot RFX exceedence probabilities
for ii = 1 : 3
    subplot(2, 3, ii)
    bar(voi_noise, exceedenceP{ii}');
    legend({'One Tau', 'Two Tau'});
    title(ep_titles{ii})
    xlabel('Gaussian noise - SD')
    ylabel('Exceendence Probability')
    ylim([0, 1]);
    xticklabels(round(mean_snr, 2))
end

sgtitle( 'Model Comparisons')


% Plot Bayes Factor for the different comparisons.
for ii = 1  : 3
    subplot(2, 3, ii + 3)
    h = {};
    for jj = 1 
        if length(size(logBF{ii})) == 2
        h{jj} = semilogy(voi_noise, (squeeze(logBF{ii}(jj, :))), markers{comparisons{ii}(jj)}, ... 
                       'Color', colors(comparisons{ii}(jj), :));
        hold on
        else
                h{jj} = semilogy(voi_noise, (squeeze(logBF{ii}(jj, :, :))), markers{comparisons{ii}(jj)}, ... 
                       'Color', colors(comparisons{ii}(jj), :));
        end
    end
    
    larray = [h{1}];
    
    legend(larray(1 : length(h{1}): end), ep_titles(comparisons{ii}));
    ylabel('Bayes Factor')
    xlabel('Gaussian noise - SD')
    
    title(sprintf('Bayes factor in Favor of %s', ep_titles{ii}));
end

% cpm_savefig(fig2, sprintf('results/model_recovery_%s.pdf', 'all'))

