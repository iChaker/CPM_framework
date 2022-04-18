% make sure to run script when inside Example_TD_learning folder
% make sure to include SPM and CPM folders and subfolders to MATLAB path
clear()

%% load dummy PRF

% this PRF structure is prededfined outside the given scripts
% This is a necessary step because:
% 1. we simulate BOLD signals using the PRF structure
% 2. the PRF structure specification requires a BOLD signal

% To solve this chicken or egg problem, we provide this PRF structure which was
% specified using an arbitrary BOLD signal having the same dimentions as the measured fMRI
load('GLMs\PRF_dummy_GLM.mat'); 


%% choose population response:

% prefix to add to all relevant files
name_tag = 'RPE_20'; % we seperate filenames using population response grid,  noise levels will be added automatically to filenames

% load population response
load(['U/U_' name_tag  '.mat'])

% fill in some missing values
for i=1:length(U)
    U(i).nbins=PRF.U(i).nbins;
    U(i).ind=PRF.U(i).ind;
end
PRF.U=U;
clear U;


%% choose observation model and population field model:
obv_fun= 'cpm_obv_int'; % observation model (Balloon model), if empty, defaults to a simple downsampling
RF_fun = ''; % population field model, defaults to Gaussian

%% choose voxels to simulate
% we choose values to be combinations of 0,1/3,2/3,1 and -1,0,1
% These are latent values, because we will save them directly into the PRF structure

lalphas = [-5,-0.4308,0.4306,5];
letas = [-0.9674,0.9674];

% population field model values
lsigma = -5; % for all points
lbeta = 0;   % for all points

% observation (hemodynamic) model values
% cpm doesn't require these parameter to be latent, it depends on the
% implementation
% In our case, these are latent parameters that will be exponentiated within the Balloon model
% implementation
transit=0.3;
decay=0.2;
epsilon=-0.5;


%% run simulation
nvoxels = length(lalphas) * length(letas);
nscans = size(PRF.xY.y,1);

% initialize location matrix
location = zeros(3,nvoxels);
location(1,:) = 1:nvoxels;

% all possible computational parameter combinations
P_voxels = combvec(lalphas,letas)';

% noise values
noise_variance = [0 0.3 0.25 0.2 0.15 0.1 0.05];

% simulating voxels for each noise value
for i =1:size(noise_variance,2)
    
    disp(['generating timeseries with noise level '  num2str(noise_variance(i))]);
    
    % timeseries
    y = zeros(nscans,nvoxels);
    Params = {};

     for v =1:nvoxels
            
            % P: all PRF parameters for voxel v
            P.lmu_alpha = P_voxels(v,1); 
            P.lmu_eta = P_voxels(v,2);

            P.lsigma_alpha = lsigma;
            P.lsigma_eta = lsigma;
            
            P.lbeta=lbeta;
            
            P.transit=transit;
            P.decay=decay;
            P.epsilon=epsilon;
            BOLD_signal = feval(PRF.M.IS, P, PRF.M, PRF.U);
            

            % adding gaussian noise to the BOLD signal

            y(:,v) = BOLD_signal + normrnd(0,noise_variance(i),[nscans,1]);
            Params{v} = P;
            
            % sanity check: display true parameters for current voxel v
            cpm_get_true_parameters(P,PRF.M,PRF.U)
            
            % sanity check: draw population field in real parameter space
            % points are too tiny to notice, raise lsigma value for more clarity
            % cpm_draw_voxel(PRF,P,'alpha','eta','',100); 
            
            clear P;
     end    

    % saving xY structure:
    
    % fMRI timeseries
    Y = mean(y,2);
    xY.y = y; 
    
    % simulation parameter values
    xY.XYZmm =  location;
    sim.Params= Params; % all PRF parameters (latent) 
    sim.params = P_voxels; % NRF parameters
  
    sim.PRF=PRF; % another sanity check, this isn't required
    
    % current noise value
    sim.noise = noise_variance(i);
    
    % saving to file with name tag and noise tag
    save(['./sim/sim_' name_tag '_' num2str(noise_variance(i)) '.mat'],'xY','Y','sim');
            
end


    
