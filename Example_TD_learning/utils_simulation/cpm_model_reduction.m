function [dF, sE, sC] = cpm_model_reduction(PRF, voxel, rE, rC)
    % Simple wrapper for model reduction, uses vecors of covariance and
    % parameters to perform the reduction, extracts priors and posterior
    % estimates from PRF

    if isempty(rE)
       rE = spm_vec(PRF.M.pE{voxel});
    end
    
    if isempty(rC)
       rC = spm_vec(PRF.M.pC{voxel}); 
    end
    % posterior expectation & covariance of the full model
    Cp = PRF.Cp{voxel};
    Ep = PRF.Ep{voxel};

    % prior expectation & covariance of the full model
    pC = PRF.M.pC{voxel};
    pE = PRF.M.pE{voxel};
   
    [dF, sE, sC] = spm_log_evidence(Ep, Cp, pE, diag(spm_vec(pC)), rE, diag(rC));

end