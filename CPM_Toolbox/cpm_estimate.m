function PRF = cpm_estimate(PRF,voxels, nograph, use_parfor)
%CPM_ESTIMATE Summary of this function goes here
%   Wrapper function to call spm_prf_analyse
%   Here you may consider splitting and combining PRFs for parallel computation 
%   using BayespRF's feautures 


% Here you would use prf's splitting and parallelize execution 

if nargin < 2
    nograph = true;
    use_parfor = false;
end

if nargin < 3
    use_parfor = false;
end


if isempty(voxels)
    options  = struct('init', 'NONE',...
                      'nograph', nograph, ...
                       'use_parfor', use_parfor);
else
    options  = struct('init','NONE',...
                      'nograph', nograph,...
                      'voxels', voxels, ...
                       'use_parfor', use_parfor);
end

PRF = spm_prf_analyse('estimate',PRF,options);


end

