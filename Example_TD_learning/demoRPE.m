

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



grid.tau = [-36.0437 36.0437 3]; % grid.Theta_i = [min max Nsteps]
%grid.eps = [-4 4 20]; % these fieldnames must be the same fieldnames used in the cpm_grid
grid.eta = [-1.5 1.5 3]; % you can use different stepsizes for each parameter

fixedparams.gamma=1;  % these fieldnames must be the same fieldnames used in the cpm_grid
fixedparams.Lambda=0.99;
% make sure to change output file name if you change grid structure
output_file = './U/U_RPEttq.mat'; % cpm_precompute does not overwrite an existing file

U = cpm_precompute(model,grid,fixedparams,data,output_file);
%% specify a PRF 


options.name='tmp'; % PRF file tag
options.TE=0.03;      % Echo time (you may add other options from BayespRF)
options.model='spm_prf_response';
outpath='GLMs';  % PRF file path

obfun= 'cpm_obv_int'; % used defined observation model, if empty, defaults to a simple downsampling
RFfun = ''; % user defined receptive field model, defaults to Gaussian

%% BMR
voxels = [];
load('SPM.mat');

load('./sim/sim_BMR_null_0.mat');
y=xY.y;               % timeseries
XYZmm=xY.XYZmm;       % locations
PRF=cpm_specify(SPM,options,y,XYZmm,U,RFfun,obfun,outpath);
PRF_nul = cpm_estimate(PRF,voxels);
save('GLMs/PRF_nul','PRF_nul');
clear PRF;

load('./sim/sim_BMR_tau_0.mat');
y=xY.y;               % timeseries
XYZmm=xY.XYZmm;       % locations
PRF=cpm_specify(SPM,options,y,XYZmm,U,RFfun,obfun,outpath);
PRF_tau = cpm_estimate(PRF,voxels);
save('GLMs/PRF_tau','PRF_tau');
clear PRF;

load('./sim/sim_BMR_taueta_0.mat');
y=xY.y;               % timeseries
XYZmm=xY.XYZmm;       % locations
PRF=cpm_specify(SPM,options,y,XYZmm,U,RFfun,obfun,outpath);
PRF_taueta = cpm_estimate(PRF,voxels);
save('GLMs/PRF_taueta','PRF_taueta');
clear PRF;


%% estimating one voxel

% voxels = [];
% PRF = cpm_estimate(PRF,voxels);
% 
% try
%     
%    r1=cpm_draw_voxel(PRF,voxels(1),'tau','eta','',100);
%    
%    PRFgen = xY.PRF;
%    Params = xY.Params;
%    
%    r2=cpm_draw_voxel(PRFgen,Params{voxels(1)},'tau','eta','',100);
%    
%    KL=sum(r1(:) .* (log2(r1(:))-log2(r2(:))))
%     
% catch
%     
% end








