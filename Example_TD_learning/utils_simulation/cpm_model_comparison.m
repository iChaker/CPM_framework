function [logBF, exceedenceP, comparisons] =  cpm_model_comparison(reducedF, voxel_idx, base_idx)
    % Helper function to do model comparison in a more efficient way than
    % previously. 
    
    % Assume that F has size reduced model x voxels
    % so compare all models ~= base_idx vs base_idx
    comparisons = 1 : size(reducedF, 1);
    comparisons = comparisons(comparisons ~= base_idx);
    
    logBF  = nan([length(comparisons), size(voxel_idx, 2), size(voxel_idx, 1)]);
    
    exceedenceP = nan(size(reducedF, 1), size(voxel_idx, 2));
   
    % Options to supress screens
    vba_options.DisplayWin = 0;
    
    cc = 1;
    
   
    for kk = voxel_idx
        
        for jj = 1 : length(comparisons)
            logBF(jj, cc, :) =  reducedF(base_idx, kk) - reducedF(comparisons(jj), kk);
        end
        
        % We also use the VBA to compute RFX exceedence probabilties for all
        % models included in reduced F
        [~, o] = VBA_groupBMC(reducedF(:, kk), vba_options);
        exceedenceP(:, cc) = o.ep;
        cc = cc + 1;
    end
    
end