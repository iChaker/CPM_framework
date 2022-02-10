function PRF = cpm_estimate(PRF,voxels)
%CPM_ESTIMATE Summary of this function goes here
%   Wrapper function to call spm_prf_analyse
%   Here you may consider splitting and combining PRFs for parallel computation 
%   using BayespRF's feautures 



% Here you would use prf's splitting and parallelize execution 


if isempty(voxels)
    options  = struct('init','NONE',...
                      'nograph',true);
else
    options  = struct('init','NONE',...
                      'nograph',false,...
                      'voxels', voxels);
end

PRF = spm_prf_analyse('estimate',PRF,options);


end

