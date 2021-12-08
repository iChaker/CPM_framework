function PRF = cpm_estimate(PRF,voxels)
%CPM_ESTIMATE Summary of this function goes here
%   Detailed explanation goes here

if isempty(voxels)
    options  = struct('init','NONE',...
                      'nograph',true);
else
    options  = struct('init','NONE',...
                      'nograph',true,...
                      'voxels', voxels);
end

PRF = spm_prf_analyse('estimate',PRF,options);


end

