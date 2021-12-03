function onevoxel = cpm_simulate(PRF,nscan,noise)
% This is a placeholder function that simulates one voxel , make sure to manually change
% parameter fields (see PRF.M.pE) and values.

    try PRF.M.pE{1,1}.lmu_transit
        
           %demoBOLD
        P.lmu_tau = 0.674489750196125;
        P.lsigma_tau =  -0.139710298879709;

        P.lmu_transit = 0.3;
        P.lmu_decay   = 0.2;
        P.lmu_epsilon = -0.5;
        P.lsigma_transit = -0.33;
        P.lsigma_decay   = -0.739710298879709;
        P.lsigma_epsilon = -0.139710298879709;
        P.lbeta = 0.4;
    catch
        
            %demoRPE
        P.lmu_tau = 0.674489750196125;
        P.lmu_eps = -0.474489750196125;
        P.lmu_eta = -0.33;

        P.lsigma_tau =  -0.139710298879709;
        P.lsigma_eps =  -0.739710298879709;
        P.lsigma_eta =  -0.33;

        P.transit = PRF.M.pE{1}.transit+0.3;
        P.decay   = PRF.M.pE{1}.decay+0.2;
        P.epsilon = PRF.M.pE{1}.epsilon-0.5;
        P.lbeta = 0.4;

    end
    
    BOLD_signal = feval(PRF.M.IS, P, PRF.M, PRF.U);
    
    % adding gaussian noise to the BOLD signal
    onevoxel = BOLD_signal + normrnd(0,noise,[nscan,1]) ;

 
end