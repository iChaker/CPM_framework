function [offset,PRF] = cpm_calculate_offset(PRF,voxels)
%CPM_CALCULATE_OFFSET Calculates offset between measured BOLD signal and
%experimental timings
%   
%   INPUT: -PRF: structure specified using cpm_specify
%          -voxels: the list of voxels to be used in the calculation
%
%   OUTPUT: -offset: the measured BOLD signal offset in milliseconds
%           -PRF: structure with estimated voxels + per voxel offsets in
%           PRF.xY.offsets
PRF.M.cpm.obv_fun = 'cpm_obv_delay';
PRF.M.cpm.obv_fun = 'cpm_RF_SoG';
Y=PRF.xY.y;
nvoxels = size(PRF.xY.y,2);
try voxels; catch voxels=''; end
if isempty(voxels)
   voxels = 1:nvoxels; 
end
global delay;
delays= zeros(1,nvoxels);

for i=1:length(voxels)
    delay = 0;
    y = Y(:,voxels(i));
    PRF.M.cpm.other.current_signal=y;
    PRF = cpm_estimate(PRF,[voxels(i)]);
    delays(voxels(i)) = delay;
end

offset = mean(delays(delays>0));
PRF.xY.offsets = delays;
PRF.xY.estimated_offset = offset; 

end

