clear()

% this is how to specify a cpm model using experimental data and BOLD
% signals
% we will use voxels simulated in simulate_fields script
% the first two section are exactly similar to simulate_fields script


%% specify a cpm models using simulated data

% cpm uses the same PRF structure as BayespRF
% It contains
%   U: input data, the precomputed grid, name of computational model
%   M: observation model, population field model and scanning parameters
%   xY, Y: the measured fMRI timeseries (in out case simulated, see simulate_fields script) 

%iterate over simulation files (noise levels)
fnames = ls('sim');
name_tag = 'RPE_20'; % files name tag 


for i=1:size(fnames,1)
    if startsWith(fnames(i,:),['sim_' name_tag]) 
        load( [ 'sim/' fnames(i,:)]); % load simuation file
        
        % 1.SPM structure (can be outside the loop)
        load('SPM.mat');
        
        % 2.BayespRF / scanning parameters
        options.name=[ name_tag '_' num2str(sim.noise) ]; % PRF file name
        options.TE=0.03;      % Echo time (you may add other options from BayespRF)
        outpath='./GLMs';  % PRF file path
        
        % 3.measured fMRI data (BOLD signals), here they are simulated by
        % simulate_fields script
        y=xY.y;               % timeseries
        XYZmm=xY.XYZmm;       % locations
        
        % 4.precomputed parameter space U (can be outside the loop)
        load(['U/U_' name_tag  '.mat'])
        
        % 5,6.cpm model functions (can be outside the loop)
        obfun= 'cpm_obv_int'; % observation model (Balloon model), if empty, defaults to a simple downsampling
        RFfun = ''; % population field model, defaults to Gaussian

        % specify a PRF for current noise level 
        % cpm_specify(1,2,3,4,5,6)
        PRF=cpm_specify(SPM,options,y,XYZmm,U,RFfun,obfun,outpath);
    end
end