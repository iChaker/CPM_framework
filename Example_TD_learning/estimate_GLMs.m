% This is written for my puny 4 core CPU, you'll have to make your own
% script to run on SCRUM, pRF framework allows for the splitting of the PRF
% so that should be useful.


fnames = ls('GLMs');
name_tag = 'RPE_20';

voxels = []; % choose voxels to estimate

% iterate over noise levels
for i=1:size(fnames,1);
    
   if startsWith(fnames(i,:),'PRF_RPE_20_')
       
       % load PRF structure for current noise level
       load( [ './GLMs/' fnames(i,:) ] );
       
       try PRF.Ep % check if already estimated
          
           disp( [ 'aleady estimated ' fnames(i,:)] ); 
       catch
           %estimate
           PRF = cpm_estimate(PRF,voxels);
           
           %save
           save([ './GLMs/' fnames(i,:) ],'PRF');
           disp(fnames(i,:));
       end
       
   end
    
    
end



