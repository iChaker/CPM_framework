function [y] = cpm_generative_grid_process(params, SPM, U, nscans, options, tmpdir, verbose)

arguments
    params
    SPM
    U
    nscans
    options
    tmpdir
    verbose = true
end

inverse_fun = @(x, xmin, xmax) (norminv((x - xmin) ./ (xmax - xmin)));

tmpPRF = cpm_specify(SPM, options, zeros(nscans, 1), ...
                                              zeros(3, 1), U, 'cpm_RF_Gaussian', 'cpm_obv_int', [tmpdir filesep]);
% We get the prior values of the latent(!!) parameters
tmpPE = tmpPRF.M.pE{1};
% As we define the to be simulated values in native space we have to
% transform them into the laten space:
pnames = fieldnames(params);
param_names = {};
for ii = 1 : length(pnames)
    tmp = strsplit(pnames{ii}, '_');
    param_names{ii} = tmp{2};
end

param_names = unique(param_names);

pc = 1;
for pn = param_names(:)'
    tmpPE.(['lmu_' pn{1}]) = inverse_fun(params.(['mu_' pn{1}]),  tmpPRF.options.cpm.mu.(pn{1})(1), ...
                                                                                   tmpPRF.options.cpm.mu.(pn{1})(2));
    tmpPE.(['lsigma_' pn{1}]) = params.(['sigma_' pn{1}]); 
    pc = pc + 1;
end


% we overwrite values accordingly in the prior structure:
if verbose
    disp(cpm_get_true_parameters(tmpPE, tmpPRF.M, tmpPRF.U))
end

[y, ~] =  spm_prf_response(tmpPE, tmpPRF.M, tmpPRF.U);
