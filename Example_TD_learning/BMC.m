%% 1 - comparion between the full model and 3 reduced models.


% number of voxels in the pRF model
function T= BMC(PRF)
    n_voxels = size(PRF.xY.y,2);
    %n_voxels = length(PRF.xY.params)

    % performing the model reduction
    for i=1:n_voxels

    % posterior expectation & covariance of the full model
    Cp = PRF.Cp{i};
    Ep = PRF.Ep{i};

    % prior expectation & covariance of the full model
    pC = PRF.M.pC{i};
    pE = PRF.M.pE{i};

    % prior expectation & covariance of the reduced models

    rE_1 = spm_vec(pE); % nul model
    rE_2 = spm_vec(pE); % learning model
    rE_3 = spm_vec(pE); % learning utility model
    
    % to put mean of latent param prior of sigmas to -10 or -100 is enough, so
    % it litterally is equal to sigma_min.
    rE_1(3) = -5;
    rE_1(4) = -5;
    rE_2(4) = -5;

    rC_1 = spm_vec(pC); % nul model
    rC_2 = spm_vec(pC); % learning model
    rC_3 = spm_vec(pC); % learning utility model
    
    rC_1(1)=0;
    rC_1(2)=0;
    rC_2(2)=0;

    % F - reduced log-evidence: ln p(y|reduced model) - ln p(y|full model)
    F_1(i) = spm_log_evidence(Ep,Cp,pE,diag(spm_vec(pC)),rE_1,diag(rC_1));
    F_2(i) = spm_log_evidence(Ep,Cp,pE,diag(spm_vec(pC)),rE_2,diag(rC_2));
    F_3(i) = spm_log_evidence(Ep,Cp,pE,diag(spm_vec(pC)),rE_3,diag(rC_3));

    end

    % we create a nmodelxnvoxels array of log-evidence
    F = [F_1+PRF.F;F_2+PRF.F;F_3+PRF.F];

    % we call VBA toolbox to invert the RFX generative model, corresponding to
    % the full model and its reduced models as concurrent models and voxels as
    % subjects.

    VBA_groupBMC(F) ;


    % create table to display reduced log-evidences of models  in matlab.
    F2 = [F_1;F_2;F_3];

    val = cell(size(PRF.Ep,2),1);
    for i=1:size(PRF.Ep,2)
        val{i} = num2str(i);
    end

    T = array2table(F2,...
                    'RowNames',{'model 1 (tau) ','model 2 (tau eps)','model 3 (tau eta)'},...
                    'VariableNames',val);

    disp("reduced models evidence for each voxel : (tau eps eta)");           
    disp(T);
end





