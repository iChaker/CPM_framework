% make sure to run script when inside Example_TD_learning folder
% make sure to include SPM and CPM folders and subfolders to MATLAB path
clear()

%% Load experimental inputs
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

%% defined grid & precompute

% 1.name of computational model
model = "cpm_grid_RPE"; % see cpm_grid_template for details

% 2.define precomputation grid
grid.alpha = [0 1 20]; % grid.parameter = [min max Nsteps]
                       % these fieldnames must be the same fieldnames used in the cpm_grid
grid.eta = [-1.5 1.5 20]; % you can use different stepsizes for each parameter

% 3.computational model hyperparameters
fixedparams.gamma=1;  
fixedparams.Lambda=0.99;

% 4. experimental input: 
% already loaded in the previous section

% 5. output file name
% prefix to add to all relevant files
name_tag = 'RPE_20'; % we seperate filenames using population response grid,  noise levels will be added automatically to filenames

output_file = ['./U/U_' name_tag '.mat']; 

% run precomputation
% cpm_precompute does not overwrite an existing file
% cpm_precompute(1,2,3,4,5,6)
U = cpm_precompute(model,grid,fixedparams,data,output_file);