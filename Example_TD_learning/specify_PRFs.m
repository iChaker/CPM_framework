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



grid.tau = [-36.0437 36.0437 80]; % grid.Theta_i = [min max Nsteps]
%grid.eps = [-4 4 20]; % these fieldnames must be the same fieldnames used in the cpm_grid
grid.eta = [-1.5 1.5 15]; % you can use different stepsizes for each parameter

fixedparams.gamma=1;  % these fieldnames must be the same fieldnames used in the cpm_grid
fixedparams.Lambda=0.99;
% make sure to change output file name if you change grid structure
output_file = './U/U_restr.mat'; % cpm_precompute does not overwrite an existing file

U = cpm_precompute(model,grid,fixedparams,data,output_file);
%% specify a PRF 

load('SPM.mat');
fnames = ls('sim');
for i=1:size(fnames,1)
    if startsWith(fnames(i,:),'sim_BMR')
        load( [ 'sim/' fnames(i,:)]);
        load('PRF_INIT_RPE');
        y=xY.y;               % timeseries
        XYZmm=xY.XYZmm;       % locations
        noise = xY.noise;     % noise
        options.name=[ 'RPE_alpha_' num2str(noise) ]; % PRF file tag
        options.TE=0.03;      % Echo time (you may add other options from BayespRF)
        outpath='./GLMs';  % PRF file path

        obfun= 'cpm_obv_int'; % user defined observation model, if empty, defaults to a simple downsampling
        RFfun = []; % user defined probability function, defaults to Gaussian

        PRF=cpm_specify(SPM,options,y,XYZmm,U,RFfun,obfun,outpath);
    end
end