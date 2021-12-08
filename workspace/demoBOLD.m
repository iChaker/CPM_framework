

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
model = "cpm_grid_BOLD"; % user defined model, see cpm_grid_template for details



grid.tau = [-2 2 3]; % grid.Theta_i = [min max Nsteps]
%grid.eps = [-2 2 5]; % these fieldnames must be the same fieldnames used in the cpm_grid
%grid.eta = [-2 2 5]; % you can use different stepsizes for each parameter

grid.transit = [0.01 2 3];
grid.decay = [0.01 2 3];
grid.epsilon = [0.01 2 3];


fixedparams.gamma=1;  % these fieldnames must be the same fieldnames used in the cpm_grid
fixedparams.Lambda=0.99;
fixedparams.TR=SPM.xY.RT;
fixedparams.nscans=SPM.nscan;

% make sure to change output file name if you change grid structure
output_file = 'workspace/U_BOLD_0.mat'; % cpm_precompute does not overwrite an existing file

U = cpm_precompute(model,grid,fixedparams,data,output_file);
%% specify a PRF 

load('SPM.mat');
load('simulation_VOI_0.mat');
y=xY.y;               % timeseries
XYZmm=xY.XYZmm;       % locations
options.name='BOLD'; % PRF file tag
options.TE=0.03;      % Echo time (you may add other options from BayespRF)
outpath='workspace';  % PRF file path

obfun= 'cpm_obv_identity'; % used defined observation model, if empty, defaults to a simple downsampling
RFfun = []; % user defined receptive field model, defaults to Gaussian

PRF=cpm_specify(SPM,options,y,XYZmm,U,RFfun,obfun,outpath);

%% simulating one voxel using CPM
% noise = 0;
% %IMPORTANT: change the variables inside cpm_simulate to conform to the ones in
% %the PRF structure
% onevoxel = cpm_simulate(PRF,2707,noise); %this function is hardcoded for 1 voxel, I use latent parameters as much as possible to avoid headaches
% 
% y(:,1)=onevoxel;
% PRF=cpm_specify(SPM,options,y,XYZmm,U,RFfun,obfun,outpath); % updating timeseries

%% estimating one voxel

voxels = [1];
PRFn = cpm_estimate(PRF,voxels);
RE= spm_prf_response(PRFn.Ep{1,1},PRFn.M,PRFn.U,'get_parameters') % display estimated latent parameters for voxel 1









