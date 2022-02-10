% make sure to run script when inside Example_TD_learning folder
clear()
%% Load data
file="tsvfile1.tsv";
dt = 0.074;
[trial_n, timestep_n, stimuli_n, S, stimuli, C, onsets] = cpm_prepare_experiment_data(file);
durations = ones(3 * trial_n, 1) .* 2 .* dt;
dts = ones(3 * trial_n, 1) .* dt;

data.ons = onsets;
data.dur = durations;
data.dt = dts;
data.C = C;
data.S = S;

% any input data structure must contain the fields: ons dur dt
% Other than that, you may add fields in which ever way you want (in this
% case C and S)



%% Select grid model & precompute
model = "cpm_grid_RPE"; % user defined model, see cpm_grid_template for details



grid.tau = [-36.0437 36.0437 4]; % grid.Theta_i = [min max Nsteps]
%grid.eps = [-4 4 20]; % these fieldnames must be the same fieldnames used in the cpm_grid
grid.eta = [-1.5 1.5 4]; % you can use different stepsizes for each parameter

fixedparams.gamma=1;  % these fieldnames must be the same fieldnames used in the cpm_grid
fixedparams.Lambda=0.99;
% make sure to change output file name if you change grid structure
output_file = './U/U_restr.mat'; % cpm_precompute does not overwrite an existing file

U = cpm_precompute(model,grid,fixedparams,data,output_file);
%% specify a PRF 

load('SPM.mat');
load('./sim_VOI_point.mat');
y=xY.y;               % timeseries
XYZmm=xY.XYZmm;       % locations
options.name='pRPE'; % PRF file tag
options.TE=0.03;      % Echo time (you may add other options from BayespRF)
options.model = 'spm_prf_response';

outpath='GLMs';  % PRF file path

obfun= 'cpm_obv_int'; % used defined observation model, if empty, defaults to a simple downsampling
RFfun = ''; % user defined receptive field model, defaults to Gaussian

PRF=cpm_specify(SPM,options,y,XYZmm,U,RFfun,obfun,outpath);

%% specify locations

% values are combinations of 0.5 0 and -0.5
% These are latent values
taus = [0,-1,-0.024115,0.024115,1];

etas = [-0.9674,0,0.9674];
%% specify other values
lsigma = -5;
lbeta = 0;  

transit=0.3;
decay=0.2;
epsilon=-0.5;
            
%%
nvoxels = length(taus) * length(etas);

% we set the location of our voxels
location = zeros(3,nvoxels);
location(1,:) = 1:nvoxels;

% array to store parameter set for each voxel
P_voxels = combvec(taus,etas)';

nscans = size(PRF.xY.y,1);

noise_variance = [0 0.3 0.25 0.2 0.15 0.1 0.05];

% computing the timeseries for each voxels,for each noise level.
for i =1:size(noise_variance,2)
    
    disp(['generating timeseries with noise level '  num2str(noise_variance(i))]);
    % timeseries
    y = zeros(nscans,nvoxels);
    Params = {};

     for s =1:nvoxels
            
            P.lmu_tau = P_voxels(s,1); 
            P.lmu_eta = P_voxels(s,2);

            P.lsigma_tau = lsigma;
            P.lsigma_eta = lsigma;
            
            P.lbeta=lbeta;
            
            P.transit=transit;
            P.decay=decay;
            P.epsilon=epsilon;
            BOLD_signal = feval(PRF.M.IS, P, PRF.M, PRF.U);
            

            % adding gaussian noise to the BOLD signal

            y(:,s) = BOLD_signal + normrnd(0,noise_variance(i),[nscans,1]);
            Params{s} = P;
            cpm_get_true_parameters(P,PRF.M,PRF.U)
            %cpm_draw_voxel(PRF,P,'tau','eta','',100); 
            clear P;
     end    

    % timeserie of each voxels
    xY.y = y;
    xY.PRF=PRF;
    % location of the voxels
    xY.XYZmm =  location;
    xY.Params= Params;
    % parameter of the voxels nrf
    xY.params = P_voxels;
    % mean value of the voxels timeseries
    Y = mean(y,2);
    
    xY.noise = noise_variance(i);
    save(['sim/sim_BMR_alpha_' num2str(noise_variance(i)) '.mat'],'xY','Y');
            
end


    
