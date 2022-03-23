% This is written for my puny 4 core CPU, you'll have to make your own
% script to run on SCRUM, pRF framework allows for the splitting of the PRF
% so that should be useful.


fnames = ls('GLMs');

for i=1:size(fnames,1);
    
   if startsWith(fnames(i,:),'PRF_RPE_alpha_0.m')
       load( [ './GLMs/' fnames(i,:) ] );
       try PRF.Ep
          
           disp( [ 'aleady estimated ' fnames(i,:)] ); 
       catch
           voxels = [];
           PRF = cpm_estimate(PRF,voxels);
           save([ './GLMs/' fnames(i,:) ],'PRF');
           disp(fnames(i,:));
       end
       
   end
    
    
end



