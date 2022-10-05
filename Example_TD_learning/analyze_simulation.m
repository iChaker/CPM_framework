
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
